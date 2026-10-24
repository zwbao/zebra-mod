"""The heavy sequence-to-function backends, run by zebra itself.

Each function scores one hg38 SNV and returns a *summary* (the same shape for every
model: status, model, inputs, the headline score's name and value, its units and
definition, the claim level, preflight checks). `zebra.s2f` writes that summary to a
run receipt in the case's evidence folder and reads it back into a model row.

  gpn_msa      the authors' published genome-wide hg38 table (songlab/gpn-msa-hg38-scores),
               read over HTTP with a pure-Python tabix reader: byte ranges of a BGZF file
               located through its .tbi index. No key, no download of the table, no tabix.
  evo2         NVIDIA's hosted Evo 2 `forward` endpoint (NVCF_RUN_KEY or EVO2_API_KEY):
               REF and ALT windows from the Ensembl reference, the unembed logits decoded
               from the returned .npy, the causal log-softmax summed over the window.
  alphagenome  Google DeepMind's hosted AlphaGenome API through its official Python
               client (`pip install alphagenome`, gRPC), run in a child interpreter
               ($ZEBRA_ALPHAGENOME_PYTHON, else this one) by zebra/runners/alphagenome_variant.py.

Why GPN-MSA is a lookup and not a forward pass: songlab/gpn-msa-sapiens is conditioned on
an 89-vertebrate alignment (445 auxiliary channels); called on sequence alone those
channels are zero, a state never seen in training, and the score comes out wrong (wrong
sign on KRAS G12D). The published table is the model run with its alignment.
"""

from __future__ import annotations

import ast
import base64
import io
import json
import math
import os
import struct
import subprocess
import sys
import time
import urllib.error
import urllib.request
import zipfile
import zlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from zebra.http import _opener, cache_dir, deadline_seconds

# ------------------------------------------------------------------ GPN-MSA

GPN_SCORES_URL = os.environ.get(
    "ZEBRA_GPN_SCORES_URL",
    "https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/main/scores.tsv.bgz",
)
GPN_CUTOFF = -7.0
# Published values the reader checks itself against when one of them is asked for.
GPN_KNOWN = {
    ("12", 25245350, "C", "T"): (-12.06, "KRAS G12D"),
    ("7", 140753336, "A", "T"): (-12.63, "BRAF V600E"),
    ("17", 7675088, "C", "T"): (-6.50, "TP53 R175H"),
}
BGZF_MAX_BLOCK = 65536


class RangeError(Exception):
    """A byte-range read failed (transport, not 'absent')."""


def http_range(url: str, start: int, end: int, timeout: float) -> bytes:
    """Bytes [start, end] (inclusive) of `url`; refuses a server that ignores Range."""
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}", "User-Agent": "zebra-mod"})
    last = None
    for attempt in range(3):
        try:
            with _opener().open(req, timeout=timeout) as resp:
                status = getattr(resp, "status", resp.getcode())
                if status != 206:
                    raise RangeError(f"the server answered {status} to a byte-range request (expected 206)")
                return resp.read()
        except RangeError:
            raise
        except (urllib.error.URLError, OSError, TimeoutError) as err:
            last = err
            if attempt < 2:
                time.sleep(1.5 * (attempt + 1))
    raise RangeError(f"byte range {start}-{end} of {url} failed: {last}")


def http_get_bytes(url: str, timeout: float) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "zebra-mod"})
    with _opener().open(req, timeout=timeout) as resp:
        return resp.read()


def bgzf_decompress(data: bytes) -> Tuple[bytes, List[Tuple[int, int]]]:
    """Inflate consecutive BGZF blocks. Returns the bytes and, per block,
    (compressed offset within `data`, uncompressed offset within the result).
    A truncated last block is dropped, not half-read."""
    out = bytearray()
    blocks: List[Tuple[int, int]] = []
    p = 0
    while p + 18 <= len(data):
        if data[p:p + 4] != b"\x1f\x8b\x08\x04":
            raise ValueError(f"not a BGZF block at byte {p}")
        xlen = struct.unpack_from("<H", data, p + 10)[0]
        bsize = None
        q = p + 12
        while q < p + 12 + xlen:
            si1, si2, slen = data[q], data[q + 1], struct.unpack_from("<H", data, q + 2)[0]
            if si1 == 66 and si2 == 67 and slen == 2:
                bsize = struct.unpack_from("<H", data, q + 4)[0]
            q += 4 + slen
        if bsize is None:
            raise ValueError("BGZF block without its BC size field")
        if p + bsize + 1 > len(data):
            break
        blocks.append((p, len(out)))
        out += zlib.decompress(data[p:p + bsize + 1], 31)
        p += bsize + 1
    return bytes(out), blocks


def parse_tbi(raw: bytes) -> Dict[str, Any]:
    """A tabix index (already inflated): names, column layout, bins and linear index."""
    if raw[:4] != b"TBI\x01":
        raise ValueError("not a tabix index")
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack_from("<8i", raw, 4)
    p = 36
    names = raw[p:p + l_nm].split(b"\x00")[:n_ref]
    p += l_nm
    refs = []
    for _ in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", raw, p)
        p += 4
        bins: Dict[int, List[Tuple[int, int]]] = {}
        for _ in range(n_bin):
            b, n_chunk = struct.unpack_from("<Ii", raw, p)
            p += 8
            chunks = list(struct.iter_unpack("<QQ", raw[p:p + 16 * n_chunk]))
            p += 16 * n_chunk
            bins[b] = chunks
        (n_intv,) = struct.unpack_from("<i", raw, p)
        p += 4
        ioff = list(struct.unpack_from(f"<{n_intv}Q", raw, p))
        p += 8 * n_intv
        refs.append({"bins": bins, "ioff": ioff})
    return {"names": [n.decode() for n in names], "format": fmt, "col_seq": col_seq, "col_beg": col_beg,
            "col_end": col_end, "meta": chr(meta), "skip": skip, "refs": refs}


def reg2bins(beg: int, end: int) -> List[int]:
    """Bins overlapping the 0-based half-open [beg, end) (UCSC binning, 14-bit leaves)."""
    end -= 1
    bins = [0]
    for shift, offset in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        bins.extend(range(offset + (beg >> shift), offset + (end >> shift) + 1))
    return bins


def _tbi(url: str, timeout: float) -> Dict[str, Any]:
    """The .tbi of `url`, cached on disk (it is fixed for a published table)."""
    import hashlib

    path = cache_dir() / "s2f" / f"{hashlib.sha256(url.encode()).hexdigest()[:16]}.tbi"
    if not path.is_file():
        path.parent.mkdir(parents=True, exist_ok=True)
        data = http_get_bytes(url + ".tbi", timeout)
        tmp = path.with_suffix(".part")
        tmp.write_bytes(data)
        tmp.replace(path)
    raw, _ = bgzf_decompress(path.read_bytes())
    return parse_tbi(raw)


def tabix_rows(url: str, chrom: str, pos: int, timeout: float = 120.0) -> List[List[str]]:
    """Rows of a remote bgzipped, tabix-indexed table at one 1-based position."""
    idx = _tbi(url, timeout)
    names = idx["names"]
    key = chrom if chrom in names else (f"chr{chrom}" if f"chr{chrom}" in names else None)
    if key is None:
        return []
    ref = idx["refs"][names.index(key)]
    beg, end = pos - 1, pos
    min_off = ref["ioff"][beg >> 14] if (beg >> 14) < len(ref["ioff"]) else 0
    chunks = sorted(c for b in reg2bins(beg, end) for c in ref["bins"].get(b, []) if c[1] > min_off)
    zero_based = bool(idx["format"] & 0x10000)
    rows: List[List[str]] = []
    seen = set()
    for vbeg, vend in chunks:
        vbeg = max(vbeg, min_off)
        cstart, ustart = vbeg >> 16, vbeg & 0xFFFF
        cend, uend = vend >> 16, vend & 0xFFFF
        data = http_range(url, cstart, cend + BGZF_MAX_BLOCK - 1, timeout)
        text, blocks = bgzf_decompress(data)
        stop = next((u for c, u in blocks if c == cend - cstart), len(text)) + uend
        for line in text[ustart:stop].decode("utf-8", "replace").split("\n"):
            if not line or line.startswith(idx["meta"]):
                continue
            f = line.split("\t")
            try:
                start = int(f[idx["col_beg"] - 1]) + (1 if zero_based else 0)
            except (IndexError, ValueError):
                continue
            if f[idx["col_seq"] - 1] == key and start == pos and line not in seen:
                seen.add(line)
                rows.append(f)
            elif start > pos:
                break
    return rows


def gpn_msa_score(chrom: str, pos: int, ref: str, alt: str, timeout: float = 120.0) -> Dict[str, Any]:
    ref, alt = ref.upper(), alt.upper()
    inputs = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "assembly": "hg38"}
    summary: Dict[str, Any] = {
        "model": {"id": "songlab/gpn-msa-sapiens", "checkpoint": "published genome-wide score table",
                  "version": "hg38 scores.tsv.bgz"},
        "source_url": GPN_SCORES_URL,
        "inputs": inputs,
        "headline_score": "score",
        "units": "nats (log-likelihood ratio)",
        "score_definitions": {
            "score": "log P(ALT) - log P(REF) under GPN-MSA, both strands averaged; negative = alternate less "
                     "likely than reference. Effective range about [-13.7, +11].",
        },
        "readout_class": "function_measured_proxy",
        "claim_level": "molecular",
        "higher_is_better": False,
    }
    try:
        rows = tabix_rows(GPN_SCORES_URL, chrom, pos, timeout=timeout)
    except Exception as err:  # transport, not absence
        summary.update(status="error", lookup_status="error",
                       error=f"reading the published table failed ({type(err).__name__}: {err}); "
                             "this is not evidence that the variant is absent")
        return summary
    score = None
    for f in rows:
        if len(f) >= 5 and f[2].upper() == ref and f[3].upper() == alt:
            try:
                score = float(f[4])
            except ValueError:
                score = None
            break
    lookup_status = "scored" if score is not None else ("present_but_null" if rows else "absent")
    check = None
    known = GPN_KNOWN.get((chrom, pos, ref, alt))
    if known and score is not None:
        check = {"variant": known[1], "expected": known[0], "observed": score, "agrees": abs(score - known[0]) <= 0.05}
    summary.update(
        status="ok", lookup_status=lookup_status, score=score, rows_at_position=len(rows),
        below_suggested_cutoff=(score < GPN_CUTOFF) if score is not None else None,
        self_check=check,
        preflight_steps=[
            {"step": "assembly", "ok": True, "detail": "hg38 SNV, as the published table requires"},
            {"step": "lookup_status", "ok": lookup_status == "scored", "detail": lookup_status},
            {"step": "published_value_self_check", "ok": bool(check and check["agrees"]),
             "detail": check or "no published calibration value for this variant"},
        ],
    )
    return summary


# -------------------------------------------------------------------- Evo 2

EVO2_BASE = os.environ.get("ZEBRA_EVO2_BASE", "https://health.api.nvidia.com/v1/biology/arc")
EVO2_MODELS = ("evo2-7b", "evo2-40b")
NVCF_STATUS = "https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/"


def evo2_key() -> Optional[str]:
    return os.environ.get("NVCF_RUN_KEY") or os.environ.get("EVO2_API_KEY") or None


def read_npy(raw: bytes) -> Tuple[Tuple[int, ...], List[float]]:
    """A .npy array as (shape, flat C-order floats). float32/float16/float64 and bfloat16
    (stored as 2-byte void, the upper half of a float32)."""
    if raw[:6] != b"\x93NUMPY":
        raise ValueError("not a .npy file")
    major = raw[6]
    if major == 1:
        (hlen,) = struct.unpack_from("<H", raw, 8)
        start = 10
    else:
        (hlen,) = struct.unpack_from("<I", raw, 8)
        start = 12
    header = ast.literal_eval(raw[start:start + hlen].decode("latin1"))
    if header.get("fortran_order"):
        raise ValueError("Fortran-ordered arrays are not supported")
    shape = tuple(header["shape"])
    n = 1
    for d in shape:
        n *= d
    descr = header["descr"]
    body = raw[start + hlen:]
    order = ">" if descr.startswith(">") else "<"
    kind = descr.lstrip("<>|=")
    if kind in ("f4", "f2", "f8"):
        code = {"f4": "f", "f2": "e", "f8": "d"}[kind]
        vals = list(struct.unpack_from(f"{order}{n}{code}", body))
    elif kind in ("V2", "u2", "bfloat16"):
        halves = struct.unpack_from(f"{order}{n}H", body)
        vals = list(struct.unpack(f"<{n}f", struct.pack(f"<{n}I", *[h << 16 for h in halves])))
    else:
        raise ValueError(f"unsupported .npy dtype {descr!r}")
    return shape, vals


def token_logprobs(shape: Sequence[int], flat: Sequence[float], sequence: str) -> List[float]:
    """log P(base_i | prefix) for every base of `sequence` (position 0 unconditioned = 0.0).
    Evo 2 is causal and byte-level: logits at i predict the token at i+1, a base's token id
    is its ASCII code."""
    vocab = shape[-1]
    length = shape[-2]
    out = [0.0] * len(sequence)
    for i in range(min(len(sequence) - 1, length - 1)):
        row = flat[i * vocab:(i + 1) * vocab]
        m = max(row)
        lse = m + math.log(sum(math.exp(x - m) for x in row))
        out[i + 1] = row[ord(sequence[i + 1])] - lse
    return out


def _post_json(url: str, payload: Dict[str, Any], key: str, timeout: float) -> Tuple[int, Dict[str, str], bytes]:
    req = urllib.request.Request(url, data=json.dumps(payload).encode(), method="POST",
                                 headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json",
                                          "Accept": "application/json", "NVCF-POLL-SECONDS": "60"})
    try:
        with _opener().open(req, timeout=timeout) as resp:
            return getattr(resp, "status", resp.getcode()), dict(resp.headers), resp.read()
    except urllib.error.HTTPError as err:
        return err.code, dict(err.headers or {}), err.read()


def evo2_forward(sequence: str, model: str, key: str, timeout: float) -> Tuple[Tuple[int, ...], List[float]]:
    """The unembed logits for `sequence` from one hosted model; raises on failure."""
    status, headers, body = _post_json(f"{EVO2_BASE}/{model}/forward",
                                       {"sequence": sequence, "output_layers": ["unembed"]}, key, timeout)
    deadline = time.monotonic() + timeout
    while status == 202:  # NVCF queued the request: poll its status URL
        reqid = headers.get("NVCF-REQID") or headers.get("nvcf-reqid")
        if not reqid or time.monotonic() > deadline:
            raise RuntimeError(f"{model}: request queued and not finished in {timeout:g}s")
        req = urllib.request.Request(NVCF_STATUS + reqid, headers={"Authorization": f"Bearer {key}",
                                                                   "NVCF-POLL-SECONDS": "60"})
        try:
            with _opener().open(req, timeout=timeout) as resp:
                status, headers, body = getattr(resp, "status", resp.getcode()), dict(resp.headers), resp.read()
        except urllib.error.HTTPError as err:
            status, headers, body = err.code, dict(err.headers or {}), err.read()
    if status != 200:
        raise RuntimeError(f"{model}: HTTP {status}: {' '.join(body.decode('utf-8', 'replace').split())[:200]}")
    data = json.loads(body)
    raw = base64.b64decode(data["data"])
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        name = next((n for n in zf.namelist() if n.startswith("unembed") and n.endswith(".npy")), None)
        if name is None:
            raise RuntimeError(f"{model}: no unembed layer in the response ({zf.namelist()})")
        return read_npy(zf.read(name))


def evo2_score(chrom: str, pos: int, ref: str, alt: str, window: int = 2048, timeout: float = 240.0,
               reference=None) -> Dict[str, Any]:
    """Evo 2's zero-shot variant score: ALT minus REF summed log-likelihood over a window
    centred on the variant. `reference(chrom, start, end)` returns the hg38 sequence."""
    ref, alt = ref.upper(), alt.upper()
    if window % 2 or not 64 <= window <= 8192:
        raise ValueError("the Evo 2 window must be an even number of bases between 64 and 8192")
    half = window // 2
    start = pos - half + 1
    end = start + window - 1
    inputs = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "assembly": "hg38",
              "window_start": start, "window_end": end, "window_len": window}
    summary: Dict[str, Any] = {
        "inputs": inputs,
        "headline_score": "delta_loglik",
        "units": "nats (summed log-probability difference over the window)",
        "score_definitions": {
            "delta_loglik": f"sum over the {window:,} bp window of log P(base_i | prefix) for the ALT sequence "
                            "minus the REF sequence, from the causal log-softmax of the unembed logits. Negative = "
                            "ALT window less likely under the model; a constraint proxy, no tissue, no direction "
                            "of effect.",
            "delta_logprob_at_variant": "the same difference at the variant position only (one base).",
        },
        "readout_class": "function_measured_proxy",
        "claim_level": "molecular",
        "higher_is_better": False,
    }
    key = evo2_key()
    if not key:
        summary.update(status="error", error="no NVCF_RUN_KEY or EVO2_API_KEY")
        return summary
    if reference is None:
        from zebra.sources import ensembl

        def reference(c, s, e):
            return ensembl.sequence(c, s, e, assembly="GRCh38").result
    ref_window = str(reference(chrom, start, end)).upper()
    idx = pos - start
    if len(ref_window) != window:
        summary.update(status="error", error=f"reference window has {len(ref_window)} bases, expected {window}")
        return summary
    if ref_window[idx] != ref:
        summary.update(status="error", error=f"REF mismatch at chr{chrom}:{pos} (hg38): given {ref}, genome has {ref_window[idx]}")
        return summary
    alt_window = ref_window[:idx] + alt + ref_window[idx + 1:]
    errors = []
    preferred = os.environ.get("ZEBRA_EVO2_MODEL")
    models = [preferred] if preferred else list(EVO2_MODELS)
    for model in models:
        try:
            began = time.monotonic()
            rs, rv = evo2_forward(ref_window, model, key, timeout)
            left = max(5.0, timeout - (time.monotonic() - began))
            as_, av = evo2_forward(alt_window, model, key, left)
        except Exception as err:
            errors.append(f"{type(err).__name__}: {err}")
            continue
        r_lp = token_logprobs(rs, rv, ref_window)
        a_lp = token_logprobs(as_, av, alt_window)
        delta = sum(a_lp) - sum(r_lp)
        summary.update(
            status="ok",
            # the model that actually answered, never the one first asked for
            model={"id": model, "checkpoint": None, "version": "NVIDIA hosted (forward)"},
            delta_loglik=delta,
            delta_logprob_at_variant=a_lp[idx] - r_lp[idx],
            ref_loglik=sum(r_lp), alt_loglik=sum(a_lp),
            fallback_from=models[0] if model != models[0] else None,
            preflight_steps=[
                {"step": "reference_base", "ok": True, "detail": f"REF {ref} matches Ensembl GRCh38"},
                {"step": "window", "ok": True, "detail": f"{window} bp centred on the variant, causal log-softmax"},
                {"step": "model", "ok": model == models[0],
                 "detail": model if model == models[0] else f"{models[0]} failed ({errors[-1]}); answered by {model}"},
            ],
        )
        return summary
    summary.update(status="error", model={"id": models[0], "checkpoint": None, "version": "NVIDIA hosted"},
                   error="; ".join(errors))
    return summary


# -------------------------------------------------------------- AlphaGenome

ALPHAGENOME_RUNNER = Path(__file__).resolve().parent / "runners" / "alphagenome_variant.py"


def alphagenome_python() -> str:
    return os.path.expanduser(os.environ.get("ZEBRA_ALPHAGENOME_PYTHON") or sys.executable)


def alphagenome_importable(python: Optional[str] = None) -> bool:
    try:
        r = subprocess.run([python or alphagenome_python(), "-c", "import alphagenome"],
                           capture_output=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return r.returncode == 0


def alphagenome_variant(chrom: str, pos: int, ref: str, alt: str, ontology: str, outputs: str = "RNA_SEQ",
                        timeout: float = 240.0) -> Dict[str, Any]:
    """Run the AlphaGenome child; a run cut by `timeout` comes back with timed_out=True."""
    import signal

    args = {"chrom": f"chr{chrom}", "pos": pos, "ref": ref.upper(), "alt": alt.upper(),
            "ontology": [t.strip() for t in ontology.split(",") if t.strip()],
            "outputs": [o.strip().upper() for o in outputs.split(",") if o.strip()],
            "timeout": min(120.0, timeout)}
    base: Dict[str, Any] = {
        "inputs": dict(args, chrom=chrom, assembly="hg38"),
        "headline_score": "log2fc",
        "units": "log2((ALT+1)/(REF+1)), pseudocount 1 in model-native track units",
        "score_definitions": {
            "log2fc": "mean over ALL positions of the 16,384 bp window and ALL tracks returned for the ontology "
                      "term(s) of log2((ALT+1)/(REF+1)). A full-window mean dilutes a local effect; a value near "
                      "0 means no window-wide shift, not no effect.",
            "local_max_abs_log2fc": "the largest |log2((ALT+1)/(REF+1))| within +/-500 bp of the variant, over "
                                    "all tracks (signed; which track and where are reported beside it).",
            "mean_diff": "mean over the window and tracks of (ALT - REF) in model-native units; positive = ALT higher.",
        },
        "readout_class": "cellular_annotation",
        "claim_level": "cellular",
        "higher_is_better": None,
    }
    proc = subprocess.Popen([alphagenome_python(), str(ALPHAGENOME_RUNNER), json.dumps(args)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, errors="replace",
                            start_new_session=True)
    try:
        stdout, stderr = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError, OSError):
            proc.kill()
        proc.communicate()
        return dict(base, status="error", timed_out=True, error=f"timed out after {timeout:g}s")
    line = next((ln for ln in reversed((stdout or "").splitlines()) if ln.startswith("{")), None)
    if line is None:
        tail = " ".join((stderr or stdout or f"exit {proc.returncode}").split())[-300:]
        return dict(base, status="error", error=f"the AlphaGenome runner returned nothing readable: {tail}")
    got = json.loads(line)
    got.pop("inputs", None)
    out = dict(base, **got)
    if out.get("status") == "ok":
        out["preflight_steps"] = [
            {"step": "ontology_terms", "ok": True, "detail": f"tracks restricted to {args['ontology']}"},
            {"step": "aggregation", "ok": True,
             "detail": "log2fc is a full-window mean; local_max_abs_log2fc is the +/-500 bp readout"},
        ]
    return out
