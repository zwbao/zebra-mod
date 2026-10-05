"""Sequence-to-function models for one variant: splicing, regulation, likelihood, constraint.

Each model answers on its own axis and is reported on its own, with the units and
the definition of its headline number and the highest claim it can support. Nothing
here averages, votes or builds a "consensus score": a model that did not run is
`not_run` with a reason, never "no effect".

Backends
  spliceai, pangolin   the Broad SpliceAI-lookup REST service (Google Cloud Run).
                       Interactive use only, a few requests per minute; answers are
                       cached here for 30 days and live calls are paced. The SpliceAI
                       weights are CC BY-NC 4.0 (Illumina), the hosted lookup is for
                       research use.
  alphagenome, evo2,
  gpn_msa              the `s2f` CLI from s2f-penguin (receipt per run). Found through
                       $S2F_BIN, else on PATH. alphagenome needs ALPHAGENOME_API_KEY
                       (non-commercial, explicitly not for clinical decision-making),
                       evo2 needs NVCF_RUN_KEY or EVO2_API_KEY, gpn_msa needs `tabix`
                       and answers hg38 SNVs only. A missing key is `not_run`; this
                       module never calls a hosted model without its key.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import signal
import subprocess
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from zebra.core import Outcome, UsageError
from zebra.http import get_json, source_record
from zebra.sources import ensembl

MODELS = ("spliceai", "pangolin", "alphagenome", "evo2", "gpn_msa")
SPLICE_MODELS = ("spliceai", "pangolin")
S2F_MODELS = ("alphagenome", "evo2", "gpn_msa")

# One Cloud Run service per tool per genome (the lookup page's own table; verified live
# 2026-10-05). spliceailookup-api.broadinstitute.org no longer answers from here.
LOOKUP_BASE = {
    ("spliceai", "GRCh38"): "https://spliceai-38-xwkwwwxdwq-uc.a.run.app",
    ("spliceai", "GRCh37"): "https://spliceai-37-xwkwwwxdwq-uc.a.run.app",
    ("pangolin", "GRCh38"): "https://pangolin-38-xwkwwwxdwq-uc.a.run.app",
    ("pangolin", "GRCh37"): "https://pangolin-37-xwkwwwxdwq-uc.a.run.app",
}
LOOKUP_TERMS = (
    "Broad SpliceAI-lookup (TGG): interactive research use, a few requests per minute; "
    "SpliceAI weights are CC BY-NC 4.0 (Illumina, GPLv3 code). Batch work needs a local instance."
)
# Minimum seconds between two live calls to the lookup service from this process.
LOOKUP_INTERVAL = float(os.environ.get("ZEBRA_SPLICEAI_INTERVAL", "2.0"))
_last_live_lookup = 0.0

CHROM_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT|M)$", re.I)
PRIORITY_RANK = {"MS": 3, "MP": 2, "C": 1, "N": 0}
PRIORITY_LABEL = {
    "MS": "MANE Select",
    "MP": "MANE Plus Clinical",
    "C": "Ensembl canonical",
    "N": "other transcript",
}
INSTALL_HINT = 'uv tool install "git+https://github.com/zwbao/s2f-penguin"'
DEFAULT_TIMEOUT = 240.0

CLAIM_CEILINGS = {
    "spliceai": "molecular (a predicted change in splice-site strength; RNA from the patient outranks it)",
    "pangolin": "molecular (a predicted change in splice-site probability; RNA from the patient outranks it)",
    "alphagenome": "cellular (a predicted change in a track in one tissue/cell ontology term)",
    "evo2": "molecular (sequence likelihood under a genome language model; a constraint proxy)",
    "gpn_msa": "molecular (cross-vertebrate conservation under GPN-MSA; a constraint proxy)",
}
COMBINE_RULE = (
    "Read each model on its own axis: never average, vote or rank raw scores across models. "
    "status=not_run means the axis was not measured, not that there is no effect."
)
NO_CLINICAL = "No model here reaches a clinical claim; the highest any of them speaks to is a cellular change."


# ---------------------------------------------------------------- variant input


def _std_chrom(token: str) -> Optional[str]:
    m = CHROM_RE.match((token or "").strip())
    if not m:
        return None
    c = m.group(1).upper()
    return "MT" if c == "M" else c


def _vcf_from_recoder(entry: Dict[str, Any]) -> Optional[Tuple[str, int, str, str]]:
    """Pick the primary-assembly VCF spelling out of variant_recoder's vcf_string list."""
    for s in entry.get("vcf_string") or []:
        parsed = ensembl.parse_vcf_like(str(s))
        if parsed and _std_chrom(parsed[0]):  # skips LRG_663-179178-C-T and friends
            return parsed
    return None


def normalise_variant(variant: str, assembly: str = "GRCh38") -> Outcome:
    """Resolve any accepted spelling to VCF coordinates on `assembly`.

    chrom-pos-ref-alt is taken as given. HGVS and rsIDs go through Ensembl
    variant_recoder and are read from its `vcf_string`, which is on the forward
    genomic strand -- VEP's `allele_string` is on the transcript strand and is wrong
    for a minus-strand gene (SCN1A NM_001165963.4:c.2134C>T is 2-166042334-G-A).
    """
    text = (variant or "").strip()
    if not text:
        raise UsageError("give a variant: chrom-pos-ref-alt, transcript HGVS (NM_...:c....) or an rsID")
    kind = ensembl.classify_input(text)
    if kind == "vcf":
        chrom, pos, ref, alt = ensembl.parse_vcf_like(text)  # type: ignore[misc]
        return Outcome(
            {"input": text, "kind": "vcf", "assembly": assembly, "chrom": chrom, "pos": pos,
             "ref": ref, "alt": alt, "resolved_from": "given as chrom-pos-ref-alt", "ids": [], "hgvs_g": None}
        )
    if kind == "unknown":
        raise UsageError(
            f"cannot read {text!r}: use chrom-pos-ref-alt, transcript HGVS (NM_...:c.../ENST...:c...) or an rsID"
        )
    out = ensembl.recode(text, assembly=assembly)
    warnings: List[str] = list(out.warnings)
    alleles: List[Tuple[str, Dict[str, Any]]] = []
    for item in out.result if isinstance(out.result, list) else []:
        if not isinstance(item, dict):
            continue
        for allele, entry in item.items():
            if isinstance(entry, dict) and entry.get("vcf_string"):
                alleles.append((allele, entry))
    if not alleles:
        raise UsageError(f"Ensembl variant_recoder returned no VCF spelling for {text!r}")
    if len(alleles) > 1:
        warnings.append(
            f"{text} recodes to {len(alleles)} alternate alleles ({', '.join(a for a, _ in alleles)}); "
            f"using {alleles[0][0]} -- give chrom-pos-ref-alt to pick one"
        )
    allele, entry = alleles[0]
    parsed = _vcf_from_recoder(entry)
    if not parsed:
        raise UsageError(f"variant_recoder gave no primary-assembly VCF spelling for {text!r}")
    chrom, pos, ref, alt = parsed
    hgvs_g = next((g for g in (entry.get("hgvsg") or []) if str(g).startswith("NC_")), None)
    ids = [i for i in (entry.get("id") or []) if str(i).startswith("rs")][:3]
    return Outcome(
        {"input": text, "kind": kind, "assembly": assembly, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
         "resolved_from": f"Ensembl variant_recoder ({assembly}) vcf_string", "ids": ids, "hgvs_g": hgvs_g},
        sources=out.sources,
        warnings=warnings,
    )


def vcf_key(v: Dict[str, Any], with_chr: bool = True) -> str:
    chrom = v["chrom"]
    prefix = "chr" if with_chr else ""
    return f"{prefix}{chrom}-{v['pos']}-{v['ref']}-{v['alt']}"


def is_snv(v: Dict[str, Any]) -> bool:
    return len(v["ref"]) == 1 and len(v["alt"]) == 1 and v["ref"] != v["alt"]


# ------------------------------------------------------- SpliceAI / Pangolin


def _pace_lookup() -> None:
    global _last_live_lookup
    wait = _last_live_lookup + LOOKUP_INTERVAL - time.monotonic()
    if wait > 0:
        time.sleep(wait)


def lookup(tool: str, v: Dict[str, Any], assembly: str = "GRCh38", distance: int = 500,
           mask: int = 0, gene_set: str = "basic") -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """One call to the SpliceAI-lookup service; returns (payload, source record)."""
    base = LOOKUP_BASE.get((tool, assembly))
    if base is None:
        raise UsageError(f"{tool} is served for GRCh38 and GRCh37 only (asked for {assembly})")
    params = {
        "hg": "38" if assembly == "GRCh38" else "37",
        "bc": gene_set,
        "distance": int(distance),
        "mask": int(mask),
        "variant": vcf_key(v),
    }

    def call(ttl):
        global _last_live_lookup
        _pace_lookup()
        got = get_json(
            f"{base}/{tool}/",
            source=f"{tool} (SpliceAI-lookup)",
            params=params,
            timeout=180,
            retries=1,
            cache_ttl=ttl,
            # the service answers input errors and rate limits with a JSON body; read it
            # instead of retrying, so a rate limit is reported rather than hammered
            ok_statuses=(200, 400, 429),
        )
        if not got.cached:
            _last_live_lookup = time.monotonic()
        return got

    resp = call(30 * 86400)
    data = resp.json()
    # A rate limit is about this minute, not about this variant: never serve one from the
    # cache. An input error (a REF the genome does not have) is deterministic, so it stays.
    transient = isinstance(data, dict) and (resp.status == 429 or "rate limit" in str(data.get("error", "")).lower())
    if transient and resp.cached:
        resp = call(0)
        data = resp.json()
    rec = source_record(
        f"{tool} via Broad SpliceAI-lookup", f"{params['variant']} {assembly} distance={distance} mask={mask}",
        resp, note=LOOKUP_TERMS,
    )
    if not isinstance(data, dict):
        raise ValueError(f"{tool}: unexpected response shape")
    data["_http_status"] = resp.status
    return data, rec


def _num(x: Any) -> Optional[float]:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _sort_transcripts(rows: Sequence[Dict[str, Any]], magnitude) -> List[Dict[str, Any]]:
    return sorted(rows, key=lambda r: (-PRIORITY_RANK.get(str(r.get("t_priority") or "N"), 0), -magnitude(r)))


SPLICEAI_DEFS = {
    "DS_AG": "acceptor gain: increase in the probability that a position becomes a splice acceptor",
    "DS_AL": "acceptor loss: decrease in the probability of an existing acceptor",
    "DS_DG": "donor gain: increase in the probability that a position becomes a splice donor",
    "DS_DL": "donor loss: decrease in the probability of an existing donor",
}
SPLICEAI_READING = (
    "delta >= 0.2 supports a splice effect (PP3 at ClinGen/Walker 2023 calibration), >= 0.5 is high confidence, "
    "<= 0.1 supports no effect (BP4/BP7). Read which score moved and its position, then predict the transcript: "
    "exon skipping, intron retention, cryptic exon, shifted site; in-frame or frameshift -> NMD?"
)
PANGOLIN_DEFS = {
    "DS_SG": "splice gain: increase in P(splice) at the reported position (the lookup service's score, not usage)",
    "DS_SL": "splice loss: change in P(splice) at an existing site; negative = the site weakens",
}


def _spliceai_transcript(row: Dict[str, Any], pos: int) -> Dict[str, Any]:
    deltas = []
    for key in ("DS_AG", "DS_AL", "DS_DG", "DS_DL"):
        val = _num(row.get(key))
        if val is None:
            continue
        offset = row.get("DP_" + key.split("_")[1])
        deltas.append({
            "score": key,
            "label": SPLICEAI_DEFS[key].split(":")[0],
            "delta": val,
            "offset": offset,
            "position": (pos + int(offset)) if isinstance(offset, int) else None,
            "ref_prob": _num(row.get(key + "_REF")),
            "alt_prob": _num(row.get(key + "_ALT")),
        })
    top = max(deltas, key=lambda d: d["delta"]) if deltas else None
    return {
        "gene": row.get("g_name"),
        "transcript": row.get("t_id"),
        "refseq": (row.get("t_refseq_ids") or [None])[0],
        "priority": row.get("t_priority"),
        "priority_label": PRIORITY_LABEL.get(str(row.get("t_priority") or "N"), "other transcript"),
        "strand": row.get("t_strand"),
        "biotype": row.get("t_type"),
        "deltas": deltas,
        "max_delta": top,
    }


def _aberrations(payload: Dict[str, Any], limit: int = 3) -> List[Dict[str, Any]]:
    """SAI-10k's transcript-level reading that the lookup service adds (not SpliceAI itself)."""
    block = payload.get("sai10kPredictions") or {}
    out = []
    for ab in (block.get("aberrations") or [])[:limit]:
        desc = ab.get("description") or {}
        region = ab.get("affected_region") or {}
        out.append({
            "type": ab.get("aberration_type"),
            "label": desc.get("label"),
            "size_bp": desc.get("size_bp"),
            "affects_coding": ab.get("affects_coding"),
            "frameshift": ab.get("frameshift"),
            "introduces_stop": desc.get("introduces_stop_codon"),
            "delta_type": ab.get("delta_type"),
            "region": {k: region.get(k) for k in ("region_type", "region_number", "distance_to_boundary", "nearest_boundary")},
        })
    return out


def run_spliceai(v: Dict[str, Any], assembly: str, distance: int, mask: int = 0) -> Dict[str, Any]:
    payload, rec = lookup("spliceai", v, assembly=assembly, distance=distance, mask=mask)
    out: Dict[str, Any] = {
        "model": "spliceai",
        "version": "SpliceAI via Broad SpliceAI-lookup (GENCODE basic)",
        "claim_ceiling": CLAIM_CEILINGS["spliceai"],
        "_source": rec,
        "url": rec.get("url"),
        "terms": LOOKUP_TERMS,
    }
    if payload.get("error"):
        out.update(status="error", reason=str(payload["error"]),
                   input_error=bool(payload.get("inputError")) or payload.get("_http_status") == 400)
        return out
    rows = [r for r in payload.get("scores") or [] if isinstance(r, dict)]
    if not rows:
        out.update(status="error", reason="the service returned no per-transcript scores")
        return out
    pos = int(payload.get("pos") or v["pos"])
    tx = [_spliceai_transcript(r, pos) for r in _sort_transcripts(
        rows, lambda r: max([_num(r.get(k)) or 0.0 for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")] or [0.0]))]
    chosen = tx[0]
    top = chosen.get("max_delta") or {}
    out.update(
        status="ran",
        scored_as={"chrom": payload.get("chrom"), "pos": pos, "ref": payload.get("ref"), "alt": payload.get("alt")},
        distance=payload.get("distance"),
        mask=payload.get("mask"),
        gene_set=payload.get("bc"),
        source_note=payload.get("source"),
        transcript=chosen,
        transcripts=tx[:6],
        transcripts_total=len(rows),
        headline={"name": "max_delta", "score": top.get("score"), "value": top.get("delta"),
                  "position": top.get("position"), "transcript": chosen.get("transcript"),
                  "refseq": chosen.get("refseq"), "priority": chosen.get("priority_label")},
        units="delta score, 0-1 (change in the probability that a position is a splice site)",
        definitions=dict(SPLICEAI_DEFS, DP_x="distance in bases from the variant to the position that score refers to"),
        reading=SPLICEAI_READING,
        aberrations=_aberrations(payload),
        aberrations_note="SAI-10k reading added by the lookup service (transcript-level consequence of the deltas), not a SpliceAI output",
    )
    return out


def run_pangolin(v: Dict[str, Any], assembly: str, distance: int, mask: int = 0) -> Dict[str, Any]:
    payload, rec = lookup("pangolin", v, assembly=assembly, distance=distance, mask=mask)
    out: Dict[str, Any] = {
        "model": "pangolin",
        "version": "Pangolin via Broad SpliceAI-lookup (GENCODE basic)",
        "claim_ceiling": CLAIM_CEILINGS["pangolin"],
        "_source": rec,
        "url": rec.get("url"),
        "terms": LOOKUP_TERMS,
    }
    if payload.get("error"):
        out.update(status="error", reason=str(payload["error"]),
                   input_error=bool(payload.get("inputError")) or payload.get("_http_status") == 400)
        return out
    rows = [r for r in payload.get("scores") or [] if isinstance(r, dict)]
    if not rows:
        out.update(status="error", reason="the service returned no per-transcript scores")
        return out
    pos = int(payload.get("pos") or v["pos"])
    tx = []
    for r in _sort_transcripts(rows, lambda r: max(abs(_num(r.get("DS_SG")) or 0.0), abs(_num(r.get("DS_SL")) or 0.0))):
        deltas = []
        for key, dp in (("DS_SG", "DP_SG"), ("DS_SL", "DP_SL")):
            val = _num(r.get(key))
            if val is None:
                continue
            offset = r.get(dp)
            prefix = key.split("_")[1]
            deltas.append({
                "score": key,
                "label": PANGOLIN_DEFS[key].split(":")[0],
                "delta": val,
                "offset": offset,
                "position": (pos + int(offset)) if isinstance(offset, int) else None,
                "ref_prob": _num(r.get(prefix + "_REF")),
                "alt_prob": _num(r.get(prefix + "_ALT")),
            })
        top = max(deltas, key=lambda d: abs(d["delta"])) if deltas else None
        tx.append({
            "gene": r.get("g_name"),
            "transcript": r.get("t_id"),
            "refseq": (r.get("t_refseq_ids") or [None])[0],
            "priority": r.get("t_priority"),
            "priority_label": PRIORITY_LABEL.get(str(r.get("t_priority") or "N"), "other transcript"),
            "strand": r.get("t_strand"),
            "biotype": r.get("t_type"),
            "deltas": deltas,
            "max_delta": top,
        })
    chosen = tx[0]
    top = chosen.get("max_delta") or {}
    out.update(
        status="ran",
        scored_as={"chrom": payload.get("chrom"), "pos": pos, "ref": payload.get("ref"), "alt": payload.get("alt")},
        distance=payload.get("distance"),
        mask=payload.get("mask"),
        gene_set=payload.get("bc"),
        source_note=payload.get("source"),
        transcript=chosen,
        transcripts=tx[:6],
        transcripts_total=len(rows),
        headline={"name": "largest_abs_delta", "score": top.get("score"), "value": top.get("delta"),
                  "position": top.get("position"), "transcript": chosen.get("transcript"),
                  "refseq": chosen.get("refseq"), "priority": chosen.get("priority_label")},
        units="change in P(splice) at the reported position, -1..1 (sign kept: negative DS_SL = site weakened)",
        definitions=dict(PANGOLIN_DEFS, DP_x="distance in bases from the variant to the position that score refers to"),
        reading="Pangolin is a second splicing model, not a replica of SpliceAI: agreement between the two is "
                "evidence, a numeric difference between them is not. These are changes in P(splice), not in usage.",
    )
    return out


# ------------------------------------------------------------------ s2f CLI


def s2f_bin() -> Optional[str]:
    """$S2F_BIN if it is runnable, else `s2f` on PATH."""
    env = os.environ.get("S2F_BIN")
    if env:
        p = Path(os.path.expanduser(env))
        if p.is_file() and os.access(str(p), os.X_OK):
            return str(p)
    return shutil.which("s2f")


def workspace_dir(case: Optional[str] = None, workspace: Optional[str] = None) -> Path:
    """Where s2f run directories and receipts go."""
    if workspace:
        return Path(os.path.expanduser(workspace))
    case = case or os.environ.get("ZEBRA_CASE")
    if case:
        return Path(os.path.expanduser(case)).resolve() / "evidence" / "s2f"
    return Path(os.path.expanduser("~")) / ".cache" / "zebra-mod" / "s2f-runs"


def _s2f_call(binary: str, args: Sequence[str], ws: Path, timeout: float) -> Dict[str, Any]:
    """Run `s2f run --json --workspace <ws> <args>`; cwd=ws so stray index files land there."""
    ws.mkdir(parents=True, exist_ok=True)
    cmd = [binary, "run", "--json", "--workspace", str(ws)] + list(args)
    started = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            errors="replace", cwd=str(ws), start_new_session=True)
    try:
        stdout, stderr = proc.communicate(timeout=timeout)
        timed_out = False
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError, OSError):
            proc.kill()
        stdout, stderr = proc.communicate()
        timed_out = True
    wall = round(time.time() - started, 2)
    manifest: Optional[Dict[str, Any]] = None
    text = stdout or ""
    brace = text.find("{")
    if brace >= 0:
        try:
            manifest = json.loads(text[brace:])
        except ValueError:
            manifest = None
    return {"manifest": manifest, "exit_code": proc.returncode, "wall_s": wall, "timed_out": timed_out,
            "stdout_tail": " ".join(text.split())[-400:], "stderr_tail": " ".join((stderr or "").split())[-400:],
            "command": " ".join(cmd)}


def _headline_value(summary: Dict[str, Any], name: Optional[str]) -> Any:
    """The headline scalar: top level, else inside head_stats (alphagenome puts it per head)."""
    if not name:
        return None
    if name in summary:
        return summary[name]
    heads = summary.get("head_stats")
    if isinstance(heads, dict):
        for head, stats in heads.items():
            if isinstance(stats, dict) and name in stats:
                return {"head": head, "value": stats[name]}
    return None


def _from_receipt(model: str, call: Dict[str, Any], ws: Path, timeout: float) -> Dict[str, Any]:
    out: Dict[str, Any] = {"model": model, "claim_ceiling": CLAIM_CEILINGS[model],
                           "wall_s": call["wall_s"], "s2f_command": call["command"]}
    if call["timed_out"]:
        out.update(status="error", reason=f"timed out after {timeout:g}s (raise --timeout or run it in the background)")
        return out
    manifest = call["manifest"]
    if not isinstance(manifest, dict):
        detail = call["stderr_tail"] or call["stdout_tail"] or f"exit {call['exit_code']}"
        out.update(status="error", reason=f"s2f run produced no receipt: {detail}")
        return out
    summary = manifest.get("summary") if isinstance(manifest.get("summary"), dict) else {}
    receipt = Path(manifest.get("run_dir") or (ws / "output" / "s2f" / str(manifest.get("run_id")))) / "manifest.json"
    name = summary.get("headline_score")
    out.update(
        run_id=manifest.get("run_id"),
        receipt=str(receipt),
        url=f"file://{receipt}",
        model_id=(summary.get("model") or {}).get("id"),
        model_version=(summary.get("model") or {}).get("version"),
        stack=manifest.get("stack"),
        exit_code=manifest.get("exit_code", call["exit_code"]),
        units=summary.get("units"),
        definitions=summary.get("score_definitions"),
        readout_class=summary.get("readout_class"),
        declared_claim_level=summary.get("claim_level"),
        higher_is_better=summary.get("higher_is_better"),
        inputs=summary.get("inputs"),
        preflight=summary.get("preflight_steps"),
    )
    if summary.get("lookup_status"):
        out["lookup_status"] = summary["lookup_status"]
        if summary["lookup_status"] != "scored":
            out["lookup_status_note"] = "not a score; never read as 'no effect'"
    if summary.get("status") != "ok":
        out.update(status="error", reason=str(summary.get("error") or f"s2f reported status {summary.get('status')!r}"))
        return out
    value = _headline_value(summary, name)
    out.update(status="ran", headline={"name": name, "value": value})
    if value is None:
        out["reason"] = f"the receipt carries no value for its declared headline score {name!r}"
    return out


def run_gpn_msa(v: Dict[str, Any], assembly: str, binary: Optional[str], ws: Path, timeout: float) -> Dict[str, Any]:
    base = {"model": "gpn_msa", "claim_ceiling": CLAIM_CEILINGS["gpn_msa"]}
    if not binary:
        return dict(base, status="not_run", reason=f"the s2f CLI was not found ($S2F_BIN or PATH); install: {INSTALL_HINT}")
    if assembly != "GRCh38":
        return dict(base, status="not_run", reason="GPN-MSA is published for hg38 only; lift the variant over first")
    if not is_snv(v):
        return dict(base, status="not_run", reason="GPN-MSA's published table holds SNVs only (this is an indel)")
    if not shutil.which("tabix"):
        return dict(base, status="not_run", reason="gpn_msa reads the published table with `tabix`, which is not on PATH (brew install htslib)")
    args = ["gpn_msa", "variant", "--chrom", f"chr{v['chrom']}", "--pos", str(v["pos"]), "--ref", v["ref"], "--alt", v["alt"]]
    out = _from_receipt("gpn_msa", _s2f_call(binary, args, ws, timeout), ws, timeout)
    out["note"] = ("log-likelihood ratio from the authors' published hg38 table (not a model call). "
                   "The suggested -7 cutoff is not a pathogenicity test.")
    return out


def run_evo2(v: Dict[str, Any], assembly: str, binary: Optional[str], ws: Path, timeout: float,
             window: int = 2048) -> Dict[str, Any]:
    base = {"model": "evo2", "claim_ceiling": CLAIM_CEILINGS["evo2"]}
    if not binary:
        return dict(base, status="not_run", reason=f"the s2f CLI was not found ($S2F_BIN or PATH); install: {INSTALL_HINT}")
    if not (os.environ.get("NVCF_RUN_KEY") or os.environ.get("EVO2_API_KEY")):
        return dict(base, status="not_run", reason="Evo 2 runs on NVIDIA's hosted endpoint: set NVCF_RUN_KEY (or EVO2_API_KEY) in the environment")
    if assembly != "GRCh38":
        return dict(base, status="not_run", reason="the Evo 2 runner takes hg38 coordinates only")
    if not is_snv(v):
        return dict(base, status="not_run", reason="the Evo 2 variant runner scores single-base substitutions only")
    args = ["evo2", "score", "--chrom", f"chr{v['chrom']}", "--pos", str(v["pos"]), "--ref", v["ref"],
            "--alt", v["alt"], "--variant-window-len", str(int(window))]
    out = _from_receipt("evo2", _s2f_call(binary, args, ws, timeout), ws, timeout)
    out["note"] = ("delta_loglik is ALT minus REF summed over the window, in nats: negative = the ALT window is less "
                   "likely under the model. A constraint hypothesis, with no tissue and no direction of effect. "
                   "Known upstream defect: the runner can fall back between evo2-7b and 40b while reporting 7b.")
    return out


def run_alphagenome(v: Dict[str, Any], assembly: str, binary: Optional[str], ws: Path, timeout: float,
                    ontology: Optional[str] = None, outputs: str = "RNA_SEQ") -> Dict[str, Any]:
    base = {"model": "alphagenome", "claim_ceiling": CLAIM_CEILINGS["alphagenome"]}
    if not binary:
        return dict(base, status="not_run", reason=f"the s2f CLI was not found ($S2F_BIN or PATH); install: {INSTALL_HINT}")
    if not os.environ.get("ALPHAGENOME_API_KEY"):
        return dict(base, status="not_run", reason="AlphaGenome is a hosted API: set ALPHAGENOME_API_KEY in the environment")
    if not ontology:
        return dict(base, status="not_run", reason=(
            "AlphaGenome needs the tissue or cell type to answer for: pass --ontology with a UBERON/CL CURIE "
            "(brain UBERON:0000955, liver UBERON:0002107, skeletal muscle UBERON:0001134, heart UBERON:0000948). "
            "Left out, the runner silently defaults to UBERON:0001157 (transverse colon) and would answer for the wrong tissue."))
    if assembly != "GRCh38":
        return dict(base, status="not_run", reason="the AlphaGenome API refuses hg19; lift the variant over to hg38 first")
    if not is_snv(v):
        return dict(base, status="not_run", reason="this AlphaGenome runner scores single-base substitutions only")
    args = ["alphagenome", "variant", "--assembly", "hg38", "--chrom", f"chr{v['chrom']}", "--position", str(v["pos"]),
            "--ref", v["ref"], "--alt", v["alt"], "--ontology", ontology, "--outputs", outputs]
    out = _from_receipt("alphagenome", _s2f_call(binary, args, ws, timeout), ws, timeout)
    out["ontology"] = ontology
    out["requested_outputs"] = outputs
    out["terms"] = ("AlphaGenome's hosted API is for non-commercial use and, in Google DeepMind's own terms, "
                    "must not be used for clinical decision-making.")
    out["note"] = ("log2fc is a mean over the whole 16,384 bp window and all tracks for the ontology term: a value "
                   "near 0 means no window-wide shift, not no effect. Look at the local change in the tracks.")
    return out


# ------------------------------------------------------------------- predict


def _select(models: Optional[Sequence[str]], v: Dict[str, Any], assembly: str, ontology: Optional[str],
            binary: Optional[str]) -> Tuple[List[str], List[Dict[str, str]]]:
    """Which models to run. Default: the two splice models plus every s2f model ready here."""
    if models:
        unknown = [m for m in models if m not in MODELS]
        if unknown:
            raise UsageError(f"unknown model(s): {', '.join(unknown)}; choose from {', '.join(MODELS)}")
        return list(dict.fromkeys(models)), []
    chosen = ["spliceai", "pangolin"]
    skipped: List[Dict[str, str]] = []
    snv38 = assembly == "GRCh38" and is_snv(v)
    checks = (
        ("gpn_msa", bool(binary) and snv38 and bool(shutil.which("tabix")),
         "needs the s2f CLI, `tabix`, and an hg38 SNV"),
        ("evo2", bool(binary) and snv38 and bool(os.environ.get("NVCF_RUN_KEY") or os.environ.get("EVO2_API_KEY")),
         "needs the s2f CLI, NVCF_RUN_KEY or EVO2_API_KEY, and an hg38 SNV"),
        ("alphagenome", bool(binary) and snv38 and bool(os.environ.get("ALPHAGENOME_API_KEY")) and bool(ontology),
         "needs the s2f CLI, ALPHAGENOME_API_KEY, --ontology (the disease tissue), and an hg38 SNV"),
    )
    for name, ready, why in checks:
        if ready:
            chosen.append(name)
        else:
            skipped.append({"model": name, "reason": f"not in the default set here: {why}"})
    return chosen, skipped


def predict(variant: str, assembly: str = "GRCh38", models: Optional[Sequence[str]] = None, distance: int = 500,
            workspace: Optional[str] = None, ontology: Optional[str] = None, timeout: float = DEFAULT_TIMEOUT,
            case: Optional[str] = None, mask: int = 0, evo2_window: int = 2048) -> Outcome:
    """Run the chosen sequence-to-function models on one variant, each reported on its own."""
    if assembly not in ("GRCh38", "GRCh37"):
        raise UsageError("assembly must be GRCh38 or GRCh37")
    if not 1 <= int(distance) <= 10000:
        raise UsageError("--distance must be between 1 and 10000 (the lookup service's limit)")
    norm = normalise_variant(variant, assembly=assembly)
    v = norm.result
    warnings: List[str] = list(norm.warnings)
    sources: List[Dict[str, Any]] = list(norm.sources)

    binary = s2f_bin()
    chosen, skipped = _select(models, v, assembly, ontology, binary)
    ws = workspace_dir(case=case, workspace=workspace)

    rows: List[Dict[str, Any]] = []
    for name in chosen:
        try:
            if name == "spliceai":
                row = run_spliceai(v, assembly, int(distance), mask)
            elif name == "pangolin":
                row = run_pangolin(v, assembly, int(distance), mask)
            elif name == "gpn_msa":
                row = run_gpn_msa(v, assembly, binary, ws, timeout)
            elif name == "evo2":
                row = run_evo2(v, assembly, binary, ws, timeout, window=evo2_window)
            else:
                row = run_alphagenome(v, assembly, binary, ws, timeout, ontology=ontology)
        except UsageError:
            raise
        except Exception as err:  # one model failing must not take the others down
            row = {"model": name, "status": "error", "claim_ceiling": CLAIM_CEILINGS[name],
                   "reason": f"{type(err).__name__}: {err}"}
        http_source = row.pop("_source", None)
        if isinstance(http_source, dict):
            sources.append(http_source)
        receipt = row.get("receipt")
        if isinstance(receipt, str):
            sources.append(source_record(f"s2f-penguin {name}", row.get("run_id"), url=f"file://{receipt}",
                                         note=f"run receipt; model {row.get('model_id')}"))
        if row.get("status") in ("error", "not_run"):
            warnings.append(f"{name} {row['status']}: {row.get('reason')}")
        rows.append(row)

    result = {
        "variant": v,
        "assembly": assembly,
        "models": rows,
        "models_requested": chosen,
        "models_not_selected": skipped,
        "workspace": str(ws),
        "s2f_cli": binary or f"not found; install with: {INSTALL_HINT}",
        "claim_ceilings": {r["model"]: r.get("claim_ceiling") for r in rows},
        "how_to_combine": COMBINE_RULE,
        "ceiling_note": NO_CLINICAL,
        "confirm_in_the_lab": ("RNA from a tissue that expresses the gene (blood if expressed, else fibroblasts) by "
                               "RT-PCR or RNA-seq, a minigene assay, or allele-specific expression. RNA evidence "
                               "outranks every model here."),
    }
    return Outcome(result, sources=sources, warnings=warnings,
                   query={"variant": variant, "assembly": assembly, "models": chosen, "distance": distance,
                          "ontology": ontology, "timeout": timeout})


# --------------------------------------------------------------- rendering


def _headline_text(row: Dict[str, Any]) -> str:
    h = row.get("headline") or {}
    if row["model"] in SPLICE_MODELS:
        val, score = h.get("value"), h.get("score")
        if val is None:
            return "-"
        where = f" at {row['variant_chrom']}:{h['position']}" if h.get("position") and row.get("variant_chrom") else (
            f" at pos {h['position']}" if h.get("position") else "")
        tx = h.get("refseq") or h.get("transcript")
        return f"{score} {val:+.3f}{where} on {tx} ({h.get('priority')})"
    value = h.get("value")
    if isinstance(value, dict):
        value = f"{value.get('value')} [{value.get('head')}]"
    if isinstance(value, float):
        value = f"{value:.4g}"
    return f"{h.get('name')} = {value}" if h.get("name") else "-"


def render(result: Dict[str, Any]) -> str:
    v = result["variant"]
    lines = [f"{v['input']} -> {v['chrom']}-{v['pos']}-{v['ref']}-{v['alt']} ({result['assembly']}; {v['resolved_from']})"]
    if v.get("ids"):
        lines[0] += "  " + ", ".join(v["ids"])
    for row in result["models"]:
        row = dict(row, variant_chrom=v["chrom"])
        status = row.get("status")
        head = _headline_text(row) if status == "ran" else (row.get("reason") or "")
        lines.append(f"{row['model']:<12} {status:<8} {head}")
        if status == "ran" and row["model"] == "spliceai":
            for ab in row.get("aberrations") or []:
                lines.append(f"{'':<21} SAI-10k: {ab.get('label')} ({ab.get('type')}, "
                             f"coding={ab.get('affects_coding')}, frameshift={ab.get('frameshift')})")
        if status == "ran" and row.get("lookup_status_note"):
            lines.append(f"{'':<21} lookup_status={row.get('lookup_status')}: {row['lookup_status_note']}")
    lines.append("")
    lines.append("claim ceilings:")
    for model, ceiling in result["claim_ceilings"].items():
        lines.append(f"  {model:<12} {ceiling}")
    for item in result.get("models_not_selected") or []:
        lines.append(f"  {item['model']:<12} {item['reason']}")
    terms = []
    for row in result["models"]:
        if row.get("status") == "ran" and row.get("terms") and row["terms"] not in terms:
            terms.append(row["terms"])
    if terms:
        lines.append("")
        lines.append("terms of the services used:")
        for t in terms:
            lines.append(f"  - {t}")
    lines.append("")
    lines.append(result["how_to_combine"])
    lines.append(result["ceiling_note"])
    lines.append("confirm: " + result["confirm_in_the_lab"])
    return "\n".join(lines)
