"""Variant card: one variant, normalised and annotated from primary sources.

Ensembl VEP (consequence on an overlapping transcript, MANE Select preferred;
REVEL, AlphaMissense, CADD, SpliceAI on GRCh38) -> VCF-style left-normalised
form checked against the Ensembl reference -> gnomAD (frequencies, coverage),
ClinVar (exact-allele match) and LitVar (mention count) -> `acmg_inputs` for
`zebra acmg suggest`.

Coordinates: VEP reports start/end on the forward strand but, for HGVS input
on a minus-strand transcript, the alleles on the transcript strand
(strand = -1). Alleles are reverse-complemented before the VCF form is built,
then anchored and left-aligned against the reference sequence, which is also
used to check REF.

REF is checked, not assumed (E7). For `chrom-pos-ref-alt` input the given REF
is compared with the Ensembl reference of the stated build *before* VEP is
called, and a mismatch is refused with both bases named; the other build is
checked too and named when it matches. A REF that VEP's own coordinates
contradict is refused in the same way. Without this, a variant given in the
wrong build was annotated on a fabricated allele and then reported as "absent
from gnomAD at a covered site", which is exactly the input `acmg suggest`
turns into PM2.

Transcript choice goes by overlap first (P1e). VEP reports a consequence for
every transcript near the variant and marks the ones the variant does *not*
touch with a `distance`; preferring MANE Select over the whole list therefore
annotated m.3243A>G (MELAS) on MT-ND1 64 bp away instead of the MT-TL1 tRNA it
sits in. `overlapping_transcripts` keeps only transcripts the variant is
inside, and `pick_transcript` prefers, among those, the one carrying VEP's own
`most_severe_consequence`, then MANE Select, then Ensembl canonical
protein-coding. When several genes overlap, the others are listed in
`transcript.also_overlapping` and a warning names them.
"""

from __future__ import annotations

import importlib
import re
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Dict, List, Optional, Tuple

from zebra.core import Outcome, UsageError
from zebra.http import SourceError, request
from zebra.sources import attempt, clinvar, ensembl, gnomad, record, validated_text

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
LEFT_PAD = 100
RIGHT_PAD = 20
# gnomAD v4 grpmax group sets, as gnomAD defines them per dataset; VEP's
# colocated frequencies are keyed `gnomade_*` (v4 exomes) and `gnomadg_*`
# (v4 genomes), and gnomAD excludes the Middle Eastern group from the *genome*
# grpmax only (see zebra.sources.gnomad for the quoted definition).
CONTINENTAL_EXOME = ("afr", "amr", "eas", "mid", "nfe", "sas")
CONTINENTAL_GENOME = ("afr", "amr", "eas", "nfe", "sas")
CONTINENTAL = CONTINENTAL_EXOME  # kept for callers that do not distinguish the two
# VEP marks a transcript the variant lies outside of with `distance`; these are
# the consequence terms that go with it.
NON_OVERLAP_TERMS = frozenset(("upstream_gene_variant", "downstream_gene_variant", "intergenic_variant"))


class OutOfWindow(Exception):
    """A base outside the fetched reference window was needed."""


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


class Window:
    """Forward-strand reference sequence chrom:start..end (1-based, inclusive)."""

    def __init__(self, chrom: str, start: int, seq: str):
        self.chrom, self.start, self.seq = chrom, start, seq.upper()
        self.end = start + len(seq) - 1

    def base(self, p: int) -> str:
        if p < self.start or p > self.end:
            raise OutOfWindow(p)
        return self.seq[p - self.start]

    def slice(self, a: int, b: int) -> str:
        if a < self.start or b > self.end:
            raise OutOfWindow((a, b))
        return self.seq[a - self.start:b - self.start + 1]


def left_normalize(pos: int, ref: str, alt: str, base_at: Callable[[int], str]) -> Tuple[int, str, str]:
    """Parsimonious, left-aligned VCF representation (Tan et al. 2015)."""
    ref, alt = ref.upper(), alt.upper()
    if ref == alt:
        raise ValueError("REF equals ALT")
    for _ in range(100000):
        changed = False
        if not ref or not alt:
            b = base_at(pos - 1)
            ref, alt, pos = b + ref, b + alt, pos - 1
            changed = True
        if ref and alt and ref[-1] == alt[-1] and (len(ref) > 1 or len(alt) > 1):
            ref, alt = ref[:-1], alt[:-1]
            changed = True
        if not changed:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def vep_alleles(rec: Dict[str, Any]) -> Tuple[str, List[str]]:
    parts = str(rec.get("allele_string") or "").split("/")
    if len(parts) < 2:
        raise ValueError(f"VEP allele_string {rec.get('allele_string')!r} has no alternate allele")
    return parts[0], parts[1:]


def raw_vcf(rec: Dict[str, Any], alt_raw: str, base_at: Callable[[int], str]) -> Tuple[int, str, str, Optional[str]]:
    """VEP record + one alternate allele -> (pos, ref, alt) on the forward strand, before normalisation.

    The fourth value is a REF-check message when the forward-strand REF that VEP
    implies differs from the reference sequence at start..end, else None.
    """
    ref_raw, _ = vep_alleles(rec)
    strand = rec.get("strand", 1)
    ref = "" if ref_raw == "-" else ref_raw.upper()
    alt = "" if alt_raw == "-" else alt_raw.upper()
    if strand == -1:
        ref, alt = revcomp(ref), revcomp(alt)
    start, end = int(rec["start"]), int(rec["end"])
    if not ref:  # insertion between end and start (start = end + 1)
        b = base_at(end)
        return end, b, b + alt, None
    found = "".join(base_at(p) for p in range(start, end + 1))
    mismatch = None if found == ref else f"VEP implies REF {ref} at {start}-{end}, reference has {found}"
    if not alt:  # deletion of start..end
        b = base_at(start - 1)
        return start - 1, b + ref, b, mismatch
    return start, ref, alt, mismatch


def parse_spdi(spdi: str) -> Tuple[str, int, str, str]:
    acc, pos0, dele, ins = spdi.split(":")
    return acc, int(pos0), dele, ins


def spdi_to_vcf(pos0: int, dele: str, ins: str, base_at: Callable[[int], str]) -> Tuple[int, str, str]:
    if dele.isdigit():
        n = int(dele)
        dele = "".join(base_at(p) for p in range(pos0 + 1, pos0 + n + 1))
    if dele and ins:
        return pos0 + 1, dele.upper(), ins.upper()
    b = base_at(pos0)  # the base before the event (1-based pos0)
    return pos0, (b + dele).upper(), (b + ins).upper()


def _num(x: Any) -> Optional[float]:
    if isinstance(x, (int, float)):
        return float(x)
    if isinstance(x, str):
        for part in x.split(","):
            try:
                return float(part)
            except ValueError:
                continue
    return None


def overlaps(tc: Dict[str, Any]) -> bool:
    """True when the variant lies inside this transcript.

    VEP sets `distance` only for a transcript the variant is outside of, and
    gives it an `upstream_gene_variant`/`downstream_gene_variant` term. Both
    signals are checked: the field, because it is what VEP documents, and the
    terms, because a record without `distance` must still not be read as
    overlapping if its only consequence is being near the gene.
    """
    if tc.get("distance") is not None:
        return False
    terms = set(tc.get("consequence_terms") or [])
    return bool(terms) and not terms.issubset(NON_OVERLAP_TERMS)


def overlapping_transcripts(tcs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return [t for t in tcs or [] if overlaps(t)]


def pick_transcript(tcs: List[Dict[str, Any]], gene: Optional[str] = None,
                    most_severe: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """The transcript to annotate on: overlapping first, then MANE Select.

    Order, applied to the transcripts the variant actually overlaps (all of
    them only when none overlaps):
      1. MANE Select carrying VEP's `most_severe_consequence`
      2. Ensembl canonical carrying `most_severe_consequence`
      3. MANE Select
      4. Ensembl canonical and protein-coding
      5. Ensembl canonical
      6. protein-coding
      7. the first one
    Steps 1-2 are what separates MT-TL1 (`non_coding_transcript_exon_variant`,
    the most severe term for m.3243A>G) from the MANE-less MT-ND1 64 bp away,
    and RNU4ATAC from CLASP1's MANE transcript, whose intron the same base sits
    in. `most_severe` is optional so a caller with a hand-built transcript list
    keeps the old behaviour.
    """
    pool = tcs or []
    if gene:
        same = [t for t in pool if str(t.get("gene_symbol", "")).upper() == gene.upper()]
        if same:
            pool = same
    over = overlapping_transcripts(pool)
    if over:
        pool = over

    def carries(t: Dict[str, Any]) -> bool:
        return bool(most_severe) and most_severe in (t.get("consequence_terms") or [])

    preds = [
        lambda t: t.get("mane_select") and carries(t),
        lambda t: t.get("canonical") and carries(t),
        lambda t: t.get("mane_select"),
        lambda t: t.get("canonical") and t.get("biotype") == "protein_coding",
        lambda t: t.get("canonical"),
        lambda t: t.get("biotype") == "protein_coding",
        lambda t: True,
    ]
    for pred in preds:
        for t in pool:
            if pred(t):
                return t
    return None


def nearest_transcript(tcs: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """The closest non-overlapping transcript, for a variant inside no transcript."""
    cands = [t for t in tcs or [] if isinstance(t.get("distance"), (int, float))]
    if not cands:
        return None
    return min(cands, key=lambda t: t["distance"])


def predictors(tc: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    tc = tc or {}
    am = tc.get("alphamissense") if isinstance(tc.get("alphamissense"), dict) else {}
    sp = tc.get("spliceai") if isinstance(tc.get("spliceai"), dict) else None
    splice = None
    if sp:
        splice = {k: sp.get(k) for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")}
        splice["max"] = ensembl.spliceai_max(tc)
        splice["gene"] = sp.get("SYMBOL")
    out = {
        "revel": _num(tc.get("revel", tc.get("revel_score"))),
        "alphamissense": {"score": _num(am.get("am_pathogenicity")), "class": am.get("am_class")} if am else None,
        "cadd_phred": _num(tc.get("cadd_phred")),
        "spliceai": splice,
        "sift": {"prediction": tc.get("sift_prediction"), "score": _num(tc.get("sift_score"))} if tc.get("sift_prediction") else None,
        "polyphen": {"prediction": tc.get("polyphen_prediction"), "score": _num(tc.get("polyphen_score"))} if tc.get("polyphen_prediction") else None,
        "source": "Ensembl VEP plugins on the selected transcript",
    }
    return out


def vep_frequencies(rec: Dict[str, Any], alt_raw: str) -> Optional[Dict[str, Any]]:
    """gnomAD frequencies VEP attaches to the colocated variant (fallback when gnomAD is down)."""
    alt_key = alt_raw if alt_raw else "-"
    for cv in rec.get("colocated_variants") or []:
        freqs = cv.get("frequencies") or {}
        groups = freqs.get(alt_key) or (next(iter(freqs.values())) if len(freqs) == 1 else None)
        if not groups:
            continue
        cont = {}
        for k, v in groups.items():
            m = re.match(r"^gnomad([eg])_(\w+)$", k)
            if not m or not isinstance(v, (int, float)):
                continue
            allowed = CONTINENTAL_EXOME if m.group(1) == "e" else CONTINENTAL_GENOME
            if m.group(2) in allowed:
                cont[f"{m.group(1)}:{m.group(2)}"] = v
        best = max(cont.items(), key=lambda kv: kv[1]) if cont else None
        return {"rsid": cv.get("id"), "exome_af": groups.get("gnomade"), "genome_af": groups.get("gnomadg"),
                "grpmax_af": best[1] if best else None, "grpmax_group": best[0] if best else None,
                "grpmax_basis": "highest group AF over gnomAD v4 exome groups "
                                f"{'/'.join(CONTINENTAL_EXOME)} ('e:') and genome groups "
                                f"{'/'.join(CONTINENTAL_GENOME)} ('g:'); exome and genome are not pooled",
                "note": "from VEP colocated_variants (gnomAD exome 'e:'/genome 'g:' groups; no allele numbers)"}
    return None


def _clinvar_candidates(rec: Dict[str, Any]) -> Tuple[List[str], Optional[str]]:
    vcvs, rsid = [], None
    for cv in rec.get("colocated_variants") or []:
        cid = str(cv.get("id") or "")
        if cid.startswith("rs") and rsid is None:
            rsid = cid
        for s in (cv.get("var_synonyms") or {}).get("ClinVar", []) or []:
            if str(s).startswith("VCV"):
                vcvs.append(str(s))
    return list(dict.fromkeys(vcvs))[:10], rsid


def match_clinvar(records: List[Dict[str, Any]], vcf: Optional[Tuple[str, int, str, str]], assembly: str,
                  window: Optional[Window], caid: Optional[str]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Split ClinVar records into those for exactly this allele and the others at the locus."""
    exact, others = [], []
    for r in records:
        if r.get("compound"):
            continue
        hit = False
        if caid and r.get("caid") and r["caid"] == caid:
            hit = True
        elif vcf and r.get("canonical_spdi"):
            try:
                _, pos0, dele, ins = parse_spdi(r["canonical_spdi"])
                if assembly == "GRCh37":
                    l38, l37 = r["locations"].get("GRCh38"), r["locations"].get("GRCh37")
                    if not (l38 and l37 and l38.get("start") and l37.get("start")):
                        raise ValueError("no GRCh37 location")
                    pos0 += l37["start"] - l38["start"]
                if window is not None:
                    got = left_normalize(*spdi_to_vcf(pos0, dele, ins, window.base), window.base)
                elif len(dele) == 1 and len(ins) == 1:
                    got = (pos0 + 1, dele.upper(), ins.upper())
                else:
                    raise ValueError("no reference window")
                hit = got == (vcf[1], vcf[2], vcf[3])
            except (ValueError, KeyError, OutOfWindow):
                hit = False
        (exact if hit else others).append(r)
    exact.sort(key=lambda r: -(r.get("stars") or 0))
    return exact, others


def _litvar(text: str, gene: Optional[str]) -> Outcome:
    from zebra.sources import litvar  # owned by another work package; may be absent

    return litvar.lookup(text, gene=gene)


def _litvar_compact(res: Any, rsids: List[str], hgvs_p: Optional[str] = None,
                    hgvs_c: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """LitVar records for this variant: LitVar's first suggestion per spelling, and any with our rsID.

    CP1-9: a record whose own spelling names ANOTHER allele at the same residue or base (LitVar's
    rsID-level record merges all alleles at an rsID) is excluded and reported as excluded, so a
    VUS does not inherit the papers of a common pathogenic allele at the same position.
    """
    if not isinstance(res, dict):
        return None
    matches = res.get("matches")
    if not isinstance(matches, list):
        return None
    from zebra.sources import litvar  # owned by this work package; imported lazily like `_litvar`

    ours = [m for m in matches if isinstance(m, dict) and (m.get("top") or (m.get("rsid") and m.get("rsid") in rsids))]
    recs: List[Dict[str, Any]] = []
    excluded: List[Dict[str, Any]] = []
    for m in ours:
        row = {"litvar_id": m.get("litvar_id"), "rsid": m.get("rsid"), "name": m.get("name"), "hgvs": m.get("hgvs"),
               "pmid_count": m.get("pmid_count")}
        verdict = litvar.allele_match(m, hgvs_p, hgvs_c)
        if verdict == "different":
            row["reason"] = (f"spelled {m.get('hgvs') or m.get('name')}: another allele at the same position"
                             + (f" (LitVar merges every allele at {m.get('rsid')})" if m.get("rsid") else ""))
            excluded.append(row)
            continue
        row["allele"] = "this allele" if verdict == "same" else "allele not verified from LitVar's spelling"
        recs.append(row)
    recs = recs[:4]
    counts = [r["pmid_count"] for r in recs if isinstance(r.get("pmid_count"), int)]
    out = {"records": recs, "pmid_count_max": max(counts) if counts else None,
           "excluded": excluded[:4],
           "pmids_excluded_max": max([e["pmid_count"] for e in excluded if isinstance(e.get("pmid_count"), int)] or [0]) or None,
           "note": "LitVar keeps unlinked spellings as separate records; counts overlap, do not add them"}
    if excluded:
        out["note"] += (f". {len(excluded)} record(s) excluded because they are about another allele at the same "
                        "position; their PMIDs are not this variant's literature")
    return out


OTHER_ASSEMBLY = {"GRCh38": "GRCh37", "GRCh37": "GRCh38"}


class RefMismatch(UsageError, ValueError):
    """The REF implied by the input does not match the reference at that position.

    A `UsageError` so every command that calls `card()` reports it as bad input
    (exit 2, message shown as is) without having to catch it: `zebra acmg
    suggest` does not wrap `card()`, and a plain ValueError reached the CLI as
    an `InternalError` with a traceback. Also a `ValueError`, so the callers
    that do catch ValueError keep working.
    """


def reference_bases(chrom: str, pos: int, length: int, assembly: str) -> Optional[str]:
    """The `length` reference bases at chrom:pos in `assembly`, or None if Ensembl did not answer."""
    if length < 1:
        return None
    try:
        return ensembl.sequence(chrom, pos, pos + length - 1, assembly).result
    except (SourceError, ValueError, KeyError, TypeError, AttributeError):
        return None


def ref_mismatch_message(chrom: str, pos: int, given: str, found: str, assembly: str,
                         other: Optional[str]) -> str:
    """The refusal text for a REF that does not match, naming both bases and the other build."""
    other_name = OTHER_ASSEMBLY[assembly]
    msg = (f"REF {given} does not match the {assembly} reference at {chrom}:{pos}, "
           f"which has {found}.")
    if other is None:
        msg += f" {other_name} could not be checked (Ensembl sequence unavailable)."
    elif other.upper() == given.upper():
        msg += f" {other_name} has {other} there, which matches: rerun with --assembly {other_name}."
    else:
        msg += (f" {other_name} has {other} there, which does not match either: check the position, "
                "the reference allele and the build.")
    return msg + (" Refusing to annotate: a REF that is not the reference describes an allele that does "
                  "not exist, and its absence from gnomAD or ClinVar is not evidence.")


def verify_input_ref(chrom: str, pos: int, ref: str, assembly: str,
                     sources: List[Dict[str, Any]], warnings: List[str]) -> None:
    """Refuse `chrom-pos-ref-alt` input whose REF is not the reference at that position (E7).

    Runs before VEP, because VEP's `/vep/human/region` endpoint accepts any REF
    and echoes the coordinates back, so by the time the card is built the
    fabricated allele is indistinguishable from a real one.
    """
    got = attempt("Ensembl reference sequence (REF check)",
                  lambda: ensembl.sequence(chrom, pos, pos + len(ref) - 1, assembly), warnings)
    if got is None:
        warnings.append(f"REF {ref} not checked against the {assembly} reference "
                        "(Ensembl sequence unavailable): the build is unverified")
        return
    sources += got.sources
    found = got.result
    if found.upper() == ref.upper():
        return
    other = reference_bases(chrom, pos, len(ref), OTHER_ASSEMBLY[assembly])
    raise RefMismatch(ref_mismatch_message(chrom, pos, ref.upper(), found.upper(), assembly, other))


# ------------------------------------------------------------ mtDNA: heteroplasmy, MITOMAP (CP1-4)

HETEROPLASMY_RE = re.compile(
    r"^(?P<variant>.+?)[\s,;]+(?:heteroplasmy[\s:=]*|het[\s:=]*|异质性[\s:=：]*)?(?P<level>\d+(?:\.\d+)?)\s*%$", re.I)


def split_heteroplasmy(text: str) -> Tuple[str, Optional[float]]:
    """`m.3243A>G 35%` -> (`m.3243A>G`, 0.35). Only an mtDNA variant takes a heteroplasmy level."""
    m = HETEROPLASMY_RE.match(text.strip())
    if not m:
        return text.strip(), None
    var = m.group("variant").strip()
    parsed = ensembl.parse_vcf_like(var)
    if parsed and parsed[0] != "MT":
        # a nuclear coordinate is refused here; an rsID or an NC_012920.1 HGVS is checked once VEP has placed it
        raise UsageError(f"{text!r}: a percentage is read as mtDNA heteroplasmy, and {var!r} is not a mitochondrial "
                         "variant. For a nuclear variant give the variant alone (record a mosaic allele fraction "
                         "in the case instead)")
    level = float(m.group("level")) / 100.0
    if not 0 < level <= 1:
        raise UsageError(f"heteroplasmy {m.group('level')}% is outside 0-100%")
    return var, level


def heteroplasmy_note(level: float, gnomad_mt: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        "level": level, "percent": round(level * 100, 2),
        "class": "homoplasmic (>= 95%)" if level >= 0.95 else ("heteroplasmic" if level >= 0.10 else
                                                              "below gnomAD's 10% calling floor"),
        "reading": ("Heteroplasmy is the fraction of mtDNA copies carrying the variant in the tissue tested. It "
                    "differs between tissues and can change over time, so the tested tissue matters when this "
                    "number is compared with thresholds; on its own it does not predict severity."),
    }
    if gnomad_mt and gnomad_mt.get("found") and gnomad_mt.get("max_heteroplasmy") is not None:
        out["gnomad_max_heteroplasmy"] = gnomad_mt["max_heteroplasmy"]
        out["gnomad_carriers"] = {"homoplasmic": gnomad_mt.get("ac_hom"), "heteroplasmic": gnomad_mt.get("ac_het"),
                                  "samples": gnomad_mt.get("an")}
    return out


MITOMASTER = "https://mitomap.org/mitomaster/websrvc.cgi"
MITOMAP_PAGE = "https://www.mitomap.org/foswiki/bin/view/MITOMAP/WebHome"
_MM_BOUNDARY = "zebra-mod-mitomaster-boundary"


def _multipart(fields: Dict[str, str], file_name: str, file_text: str) -> bytes:
    parts = []
    for k, v in fields.items():
        parts.append(f"--{_MM_BOUNDARY}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n")
    parts.append(f"--{_MM_BOUNDARY}\r\nContent-Disposition: form-data; name=\"file\"; filename=\"{file_name}\"\r\n"
                 f"Content-Type: text/plain\r\n\r\n{file_text}\r\n--{_MM_BOUNDARY}--\r\n")
    return "".join(parts).encode("utf-8")


def mitomap(pos: int, ref: str, alt: str) -> Outcome:
    """MITOMAP's disease and GenBank-frequency annotation for one mtDNA SNV, via the MITOMASTER web service.

    MITOMAP's own pages sit behind a browser challenge (HTTP 403 to scripts, checked 2026-10-06); the
    MITOMASTER service (`websrvc.cgi`, an SNV list in, a tab-separated table out) answers scripts.
    It does not return MITOMAP's curation status (Reported / Cfrm), so the result says so.
    """
    if len(ref) != 1 or len(alt) != 1:
        raise ValueError("MITOMASTER's SNV list takes single-base changes only")
    body = _multipart({"fileType": "snvlist", "output": "detail"}, "zebra.txt", f"q\t{int(pos)}{alt.upper()}\n")
    resp = request(MITOMASTER, source="MITOMAP (MITOMASTER)", method="POST", body=body,
                   headers={"Content-Type": f"multipart/form-data; boundary={_MM_BOUNDARY}"},
                   accept="text/plain,*/*", cache_ttl=30 * 86400, timeout=60,
                   # an HTML challenge page or an error text is never cached as the answer
                   validate=lambda t: None if "patientphenotype" in (t or "") else "not the MITOMASTER result table")
    text = validated_text(resp, "MITOMAP (MITOMASTER)", must_contain="patientphenotype")
    rows = [line.split("\t") for line in text.splitlines() if line.strip()]
    header = rows[0]
    data = [dict(zip(header, r)) for r in rows[1:]]
    hit = next((d for d in data if d.get("tpos") == str(int(pos)) and d.get("qnt", "").upper() == alt.upper()), None)
    rec = record("MITOMAP (MITOMASTER)", f"m.{int(pos)}{ref.upper()}>{alt.upper()}", resp, url=MITOMASTER,
                 note=f"MITOMASTER web service (POST, SNV list); browse MITOMAP at {MITOMAP_PAGE}")
    if hit is not None and hit.get("tnt") and hit["tnt"].upper() != ref.upper():
        raise ValueError(f"MITOMASTER's reference base at m.{int(pos)} is {hit['tnt']}, not {ref.upper()}")
    if hit is None:
        return Outcome({"found": False, "note": "MITOMASTER returned no row for this change"}, sources=[rec])

    def clean(s: Optional[str]) -> Optional[str]:
        s = re.sub(r"<[^>]+>", " ", s or "").strip()
        return " ".join(s.split()) or None

    def flag(s: Optional[str]) -> Optional[bool]:
        return {"true": True, "false": False}.get(str(s).strip().lower()) if s is not None else None

    disease = clean(hit.get("patientphenotype"))
    gb = _num(hit.get("gb_cnt"))
    out = {
        "found": True,
        # MITOMASTER writes "-<br>L(UUA/G)" for a tRNA (no gene symbol, then the tRNA)
        "locus": (clean(hit.get("calc_locus")) or "").lstrip("- ").strip() or None,
        # calc_aachange: the amino-acid change for a protein gene, MITOMASTER's MitoTIP note for a tRNA
        "change_note": clean(hit.get("calc_aachange")),
        "disease_reported": disease,
        "listed_as_disease_mutation": (None if flag(hit.get("is_mmut")) is None and flag(hit.get("is_rtmut")) is None
                                       else bool(flag(hit.get("is_mmut")) or flag(hit.get("is_rtmut")))),
        "listed_as_polymorphism": flag(hit.get("is_polymorphism")),
        "genbank_sequences_with_variant": int(gb) if gb is not None else None,
        "genbank_percent": _num(hit.get("gb_perc")),
        "conservation": clean(hit.get("conservation")),
        "status_note": ("MITOMASTER does not return MITOMAP's curation status (Reported vs Confirmed, 'Cfrm'); "
                        "read it on the MITOMAP page before relying on the disease association"),
    }
    return Outcome(out, sources=[rec])


# ------------------------------------------------------------ GRCh37 input: the GRCh38 view (CP1-2)

def grch38_view(vcf37: Tuple[str, int, str, str], gene: Optional[str]) -> Outcome:
    """Map a GRCh37 VCF allele to GRCh38 (Ensembl assembly map) and annotate it there.

    GRCh37 VEP serves no AlphaMissense or SpliceAI and gnomAD has only v2 on GRCh37, although the
    GRCh38 coordinates are computable. The REF is re-checked against the GRCh38 reference: a base
    that changed between the builds is reported and the view is not used.
    """
    chrom, pos, ref, alt = vcf37
    sources: List[Dict[str, Any]] = []
    warnings: List[str] = []
    m = ensembl.map_assembly(chrom, pos, pos + len(ref) - 1, "GRCh37", "GRCh38")
    sources += m.sources
    c38, p38 = m.result["chrom"], m.result["start"]
    seq = ensembl.sequence(c38, p38, p38 + len(ref) - 1, "GRCh38")
    sources += seq.sources
    if seq.result.upper() != ref.upper():
        raise ValueError(f"GRCh37 REF {ref} maps to {c38}:{p38}, where GRCh38 has {seq.result}: the reference "
                         "changed between builds, so the allele is not carried over")
    vid = f"{c38}-{p38}-{ref.upper()}-{alt.upper()}"
    v = ensembl.vep(vid, "GRCh38")
    sources += v.sources
    warnings += v.warnings
    rec = v.result
    tcs = [t for t in rec.get("transcript_consequences") or [] if t.get("variant_allele") in (vep_alleles(rec)[1][0], None)] \
        or rec.get("transcript_consequences") or []
    tc = pick_transcript(tcs, gene, rec.get("most_severe_consequence"))
    preds = predictors(tc)
    out: Dict[str, Any] = {
        "assembly": "GRCh38", "vcf": {"chrom": c38, "pos": p38, "ref": ref.upper(), "alt": alt.upper(), "id": vid},
        "mapped_by": "Ensembl assembly map (GRCh37 -> GRCh38), REF re-checked on GRCh38",
        "transcript": {"ensembl": (tc or {}).get("transcript_id"), "refseq": (tc or {}).get("mane_select"),
                       "gene": (tc or {}).get("gene_symbol")} if tc else None,
        "consequence_terms": (tc or {}).get("consequence_terms") or [],
        "intron": (tc or {}).get("intron"), "exon": (tc or {}).get("exon"),
        "hgvsc": (tc or {}).get("hgvsc"), "hgvsp": (tc or {}).get("hgvsp"),
        "predictors": preds,
    }
    w: List[str] = []
    gn = attempt("gnomAD v4 (GRCh38 view)", lambda: gnomad.variant(c38, p38, ref, alt, "GRCh38"), w)
    warnings += w
    if gn is not None:
        sources += gn.sources
        warnings += gn.warnings
        out["population"] = gn.result
    return Outcome(out, sources=sources, warnings=warnings)


# ------------------------------------------------------------ contracts with other work packages

def _optional_source(module: str, label: str, call: Callable[[Any], Outcome],
                     warnings: List[str]) -> Optional[Outcome]:
    """Call `zebra.sources.<module>` if this build has it; a missing module is a named gap, not silence."""
    try:
        mod = importlib.import_module(f"zebra.sources.{module}")
    except ImportError:
        warnings.append(f"{label}: not checked — zebra.sources.{module} is not available in this build")
        return None
    except Exception as err:  # noqa: BLE001 - another package's module failing at import must not sink the card
        warnings.append(f"{label}: not checked — zebra.sources.{module} failed to load ({type(err).__name__}: {err})")
        return None
    try:
        got = attempt(label, lambda: call(mod), warnings)
    except Exception as err:  # noqa: BLE001 - incl. a UsageError the module raises for an input it cannot take
        warnings.append(f"{label}: not checked — {type(err).__name__}: {str(err)[:300]}")
        return None
    if got is not None and not isinstance(got, Outcome):
        warnings.append(f"{label}: not used — zebra.sources.{module} returned {type(got).__name__}, not an Outcome")
        return None
    return got


def china_frequencies(vcf: Optional[Tuple[str, int, str, str]], warnings: List[str],
                      assembly: str = "GRCh38") -> Optional[Outcome]:
    """W4's `zebra.sources.china_freq.lookup(chrom, pos, ref, alt, assembly)` for this allele."""
    if not vcf:
        warnings.append("Chinese population frequencies: not checked — no VCF coordinates for this variant")
        return None
    c, p, r, a = vcf
    return _optional_source("china_freq", "Chinese population frequencies",
                            lambda mod: mod.lookup(c, p, r, a, assembly=assembly), warnings)


def mavedb_scores(gene: Optional[str], hgvs_p: Optional[str], hgvs_c: Optional[str],
                  warnings: List[str]) -> Optional[Outcome]:
    if not gene or not (hgvs_p or hgvs_c):
        return None
    return _optional_source("mavedb", "MaveDB functional scores",
                            lambda mod: mod.lookup(gene, hgvs_p=hgvs_p, hgvs_c=hgvs_c), warnings)


def card(variant: str, assembly: str = "GRCh38", gene: Optional[str] = None,
         heteroplasmy: Optional[float] = None, contracts: bool = True) -> Outcome:
    """The variant card. `heteroplasmy` (0-1) is an mtDNA level; `m.3243A>G 35%` carries it in the text.

    `contracts=False` skips the Chinese-cohort frequencies and MaveDB scores (other packages' modules,
    each a further 15-70 s cold): `acmg suggest` uses neither for a code and leaves them out.
    """
    text, het_text = split_heteroplasmy(variant)
    if heteroplasmy is not None and het_text is not None and abs(heteroplasmy - het_text) > 1e-9:
        raise UsageError(f"two heteroplasmy levels given ({het_text:.0%} in the variant text, {heteroplasmy:.0%} as a flag)")
    het_level = heteroplasmy if heteroplasmy is not None else het_text
    if het_level is not None and not 0 < het_level <= 1:
        raise UsageError("heteroplasmy must be a fraction in (0, 1] (or a percentage after the variant: m.3243A>G 35%)")
    if assembly not in ("GRCh38", "GRCh37"):
        raise ValueError("assembly must be GRCh38 or GRCh37")
    # A report often gives "SCN1A c.2134C>T" as a gene and a bare c. change: read it as GENE:c.…, which
    # VEP places on the gene's canonical transcript, and the gene_hgvs warning below says which one
    if gene and re.match(r"^\s*[cn]\.\S", text) and re.fullmatch(r"[A-Za-z0-9-]{2,20}", gene.strip()):
        text = f"{gene.strip()}:{text.strip()}"
    kind = ensembl.classify_input(text)
    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []

    if kind == "vcf":
        parsed = ensembl.parse_vcf_like(text)
        if parsed:
            c_in, p_in, r_in, _a_in = parsed
            verify_input_ref(c_in, p_in, r_in, assembly, sources, warnings)

    v = ensembl.vep(text, assembly)  # the anchor of the card: failure propagates
    rec = v.result
    sources += v.sources
    warnings += v.warnings
    ref_raw, alts_raw = vep_alleles(rec)
    alt_raw = alts_raw[0]
    if len(alts_raw) > 1:
        warnings.append(f"{text} has {len(alts_raw)} alternate alleles ({'/'.join(alts_raw)}); annotated {alt_raw} — "
                        "give HGVS or chrom-pos-ref-alt to choose another")
    chrom = str(rec.get("seq_region_name"))

    # reference window: REF check, anchoring and left-alignment
    window: Optional[Window] = None
    vcf: Optional[Tuple[str, int, str, str]] = None
    start, end = int(rec["start"]), int(rec["end"])
    for pad in (LEFT_PAD, 2000):
        lo = max(1, min(start, end) - pad)
        hi = max(start, end) + RIGHT_PAD
        seq = attempt("Ensembl reference sequence", lambda: ensembl.sequence(chrom, lo, hi, assembly), warnings)
        if seq is None:
            break
        sources += seq.sources
        window = Window(chrom, lo, seq.result)
        try:
            pos0, ref0, alt0, mismatch = raw_vcf(rec, alt_raw, window.base)
            if mismatch:
                # E7: never continue on a fabricated allele. `mismatch` means the
                # forward-strand REF that VEP's coordinates imply is not the
                # reference there, which is what a wrong genome build looks like.
                m = re.search(r"REF ([ACGTN]+) at (\d+)-\d+, reference has ([ACGTN]+)", mismatch)
                if m:
                    raise RefMismatch(ref_mismatch_message(
                        chrom, int(m.group(2)), m.group(1), m.group(3), assembly,
                        reference_bases(chrom, int(m.group(2)), len(m.group(1)), OTHER_ASSEMBLY[assembly])))
                raise RefMismatch(f"reference mismatch on chr{chrom} ({assembly}): {mismatch} — "
                                  "refusing to annotate; check the assembly and the HGVS")
            npos, nref, nalt = left_normalize(pos0, ref0, alt0, window.base)
            vcf = (chrom, npos, nref, nalt)
            break
        except OutOfWindow:
            vcf = None
            continue
    if vcf is None:
        if len(ref_raw) == 1 and len(alt_raw) == 1 and "-" not in (ref_raw, alt_raw):
            r, a = (revcomp(ref_raw), revcomp(alt_raw)) if rec.get("strand") == -1 else (ref_raw, alt_raw)
            vcf = (chrom, start, r.upper(), a.upper())
            warnings.append("REF not checked against the reference sequence (Ensembl sequence unavailable)")
        else:
            warnings.append("could not build the VCF form of this indel (reference sequence unavailable): "
                            "gnomAD lookup skipped; ClinVar allele match unverified")

    is_mt = gnomad.is_mito(chrom)
    if is_mt and assembly == "GRCh37" and re.match(r"^chrM\b|^chrM[-:_]", text, re.I):
        warnings.append("chrM on GRCh37: UCSC hg19's chrM is the Yoruba sequence (NC_001807), not the rCRS "
                        "(NC_012920) that Ensembl GRCh37, gnomAD and MITOMAP use; positions can differ. This was "
                        "read as rCRS — give the m. position from the report to be sure")
    if het_level is not None and not is_mt:
        raise UsageError(f"a heteroplasmy level applies to mtDNA variants only; {text} is on chromosome {chrom}")

    # transcript: overlap first, then MANE Select (P1e)
    most_severe = rec.get("most_severe_consequence")
    tcs_all = rec.get("transcript_consequences") or []
    tcs = [t for t in tcs_all if t.get("variant_allele") in (alt_raw, None)] or tcs_all
    overlapping = overlapping_transcripts(tcs)
    tc = pick_transcript(tcs, gene, most_severe)
    gene_symbol = (tc or {}).get("gene_symbol")
    overlap_genes = sorted({str(t.get("gene_symbol")) for t in overlapping if t.get("gene_symbol")})
    if not tcs:
        where = "intergenic" if most_severe == "intergenic_variant" else f"reported only as {most_severe or 'no consequence'}"
        warnings.append(f"this position lies in no transcript ({where}): there is no gene, transcript or protein "
                        "consequence to report, and no gene-level evidence applies")
    elif not overlapping:
        near = nearest_transcript(tcs) or {}
        warnings.append(f"no transcript overlaps this variant: the nearest is {near.get('gene_symbol') or '?'} "
                        f"{near.get('transcript_id') or '?'} at {near.get('distance')} bp "
                        f"({', '.join(near.get('consequence_terms') or []) or '-'}); annotated on it, but the "
                        "variant is outside every transcript — treat the gene as a neighbour, not the gene hit")
    elif len(overlap_genes) > 1:
        other_genes = [g for g in overlap_genes if g != str(gene_symbol)]
        detail = "; ".join(
            f"{g}: " + ", ".join(sorted({c for t in overlapping if str(t.get('gene_symbol')) == g
                                         for c in (t.get('consequence_terms') or [])}))
            for g in other_genes[:4])
        warnings.append(f"{len(overlap_genes)} genes overlap this position; annotated on {gene_symbol} "
                        f"({', '.join(tc.get('consequence_terms') or []) if tc else '-'}). Also overlapping — "
                        f"{detail} — use --gene to annotate on another")
    if gene and tc is not None and not overlaps(tc):
        warnings.append(f"--gene {gene} names a gene the variant does not lie in: annotated on "
                        f"{tc.get('transcript_id')} ({', '.join(tc.get('consequence_terms') or []) or '-'}"
                        + (f", {tc['distance']} bp away" if tc.get("distance") is not None else "") + "). "
                        + (f"The variant is inside {', '.join(overlap_genes)}" if overlap_genes
                           else "It lies inside no transcript") + " — drop --gene to annotate where it sits")
    if gene and gene_symbol and gene.upper() != str(gene_symbol).upper():
        genes_hit = sorted({t.get("gene_symbol") for t in tcs if t.get("gene_symbol")})
        warnings.append(f"gene mismatch: you gave {gene}; VEP places this variant in {', '.join(genes_hit) or '-'} "
                        f"(annotated on {gene_symbol}) — check the transcript and assembly")
    elif gene and not tc:
        warnings.append(f"gene {gene} given, but VEP reports no transcript consequence")

    hgvs: Dict[str, Any] = {}
    transcript = None
    consequence: Dict[str, Any] = {"most_severe": rec.get("most_severe_consequence")}
    if tc:
        mane = tc.get("mane_select")
        hgvsc = tc.get("hgvsc")
        cpart = hgvsc.split(":", 1)[1] if hgvsc and ":" in hgvsc else None
        ppart = tc["hgvsp"].split(":", 1)[1] if tc.get("hgvsp") and ":" in tc["hgvsp"] else None
        picked_overlaps = overlaps(tc)
        basis = "MANE Select" if mane else ("Ensembl canonical" if tc.get("canonical") else "first protein-coding")
        if not picked_overlaps:
            basis += ", OUTSIDE the variant"
        transcript = {"ensembl": hgvsc.split(":")[0] if hgvsc else tc.get("transcript_id"), "refseq": mane,
                      "mane_select": bool(mane or "MANE_Select" in (tc.get("mane") or [])),
                      "canonical": bool(tc.get("canonical")), "gene": gene_symbol, "hgnc_id": tc.get("hgnc_id"),
                      "biotype": tc.get("biotype"), "overlaps_variant": picked_overlaps,
                      "distance_bp": tc.get("distance"),
                      "also_overlapping": [{"gene": g, "consequences": sorted(
                          {c for t in overlapping if str(t.get("gene_symbol")) == g
                           for c in (t.get("consequence_terms") or [])})}
                          for g in overlap_genes if g != str(gene_symbol)][:6],
                      "basis": basis}
        if transcript["hgnc_id"] is not None and not str(transcript["hgnc_id"]).startswith("HGNC:"):
            transcript["hgnc_id"] = f"HGNC:{transcript['hgnc_id']}"
        hgvs = {"c": f"{mane}:{cpart}" if mane and cpart else hgvsc, "c_ensembl": hgvsc,
                "p": ppart, "p_ensembl": tc.get("hgvsp")}
        if kind == "gene_hgvs":
            # F37: `GENE:c.…` carries no transcript. VEP resolves it on the gene's
            # canonical transcript, and c. numbering is transcript-specific, so the
            # same c. position on the transcript a report used can be a different base.
            warnings.append(f"you gave only a gene name ({text}): annotated on "
                            f"{(transcript or {}).get('refseq') or (transcript or {}).get('ensembl') or '?'} "
                            f"({(transcript or {}).get('basis')}), VEP's normalised form is {hgvs.get('c')}. "
                            "c. numbering is transcript-specific — if the report used another transcript, the same "
                            "c. position is a different base. Give the transcript (NM_…:c.…) to be sure")
        if kind in ("hgvs", "gene_hgvs") and cpart:
            acc_in, _, change_in = text.partition(":")
            if mane and acc_in.split(".")[0].upper() == mane.split(".")[0].upper():
                notes = []
                if acc_in.upper() != mane.upper():
                    notes.append(f"input transcript version {acc_in}, MANE Select is {mane}")
                if change_in != cpart:
                    notes.append(f"VEP's normalised form is {mane}:{cpart} (input {text})")
                if notes:
                    hgvs["note"] = "; ".join(notes)
            elif mane and acc_in.upper().startswith(("NM_", "ENST")) and acc_in.split(".")[0] != (hgvsc or "").split(".")[0]:
                hgvs["note"] = f"input is on {acc_in}; MANE Select is {mane} ({hgvs['c']})"
        consequence.update({"terms": tc.get("consequence_terms") or [], "impact": tc.get("impact"),
                            "exon": tc.get("exon"), "intron": tc.get("intron"),
                            "protein_position": tc.get("protein_start"), "amino_acids": tc.get("amino_acids"),
                            "codons": tc.get("codons"), "biotype": tc.get("biotype")})
        if not mane:
            warnings.append(f"no MANE Select transcript in VEP output ({assembly}); annotated {transcript['basis']} "
                            f"{transcript['ensembl']}" + (" — GRCh37 VEP does not report MANE" if assembly == "GRCh37" else ""))
    else:
        consequence["terms"] = [rec.get("most_severe_consequence")] if rec.get("most_severe_consequence") else []
    preds = predictors(tc)
    if assembly == "GRCh37":
        pass  # filled from the GRCh38 view below, or warned about there
    elif is_mt:
        pass  # SpliceAI is not computed for the mitochondrial genome
    elif tc and preds["spliceai"] is None:
        warnings.append("SpliceAI: no precomputed score from VEP for this variant (VEP serves SNVs and short indels); see `zebra s2f`")
    if preds["spliceai"] and preds["spliceai"].get("gene") and gene_symbol and preds["spliceai"]["gene"] != gene_symbol:
        warnings.append(f"SpliceAI scores are for {preds['spliceai']['gene']}, not {gene_symbol}")

    vcvs, rsid = _clinvar_candidates(rec)
    rsids = [cv.get("id") for cv in rec.get("colocated_variants") or [] if str(cv.get("id", "")).startswith("rs")]

    # gnomAD, ClinVar, LitVar (and the GRCh38 view, MITOMAP, MaveDB) in parallel (different hosts)
    w_gn: List[str] = []
    w_cv: List[str] = []
    w_lv: List[str] = []
    w_38: List[str] = []
    w_mm: List[str] = []
    w_mave: List[str] = []
    w_cn: List[str] = []

    def do_gnomad():
        if not vcf:
            return None
        return attempt("gnomAD", lambda: gnomad.variant(vcf[0], vcf[1], vcf[2], vcf[3], assembly), w_gn)

    def do_clinvar():
        if not (vcvs or rsid or vcf or kind == "hgvs"):
            w_cv.append("ClinVar not searched: no VCV id, rsID, HGVS or position to search with")
            return None

        def run():
            pos_q = None
            if vcf:
                p = vcf[1] + 1 if (len(vcf[2]) > len(vcf[3]) and vcf[2][0] == vcf[3][0]) else vcf[1]
                pos_q = p
            term = clinvar.term_for(vcv=vcvs, rsid=rsid, hgvs=text if kind == "hgvs" else None,
                                    chrom=chrom if pos_q else None, pos=pos_q, assembly=assembly)
            found = clinvar.search(term)
            out = clinvar.summaries(found.result["ids"])
            out.sources = found.sources + out.sources
            return out
        return attempt("ClinVar", run, w_cv)

    lit_text = hgvs.get("p") if (gene_symbol and hgvs.get("p")) else (rsid or hgvs.get("c") or text)

    def do_litvar():
        try:
            return _litvar(lit_text, gene_symbol)
        except Exception as err:  # noqa: BLE001 - sibling module may be missing or change shape
            w_lv.append(f"LitVar unavailable: {type(err).__name__}: {err}")
            return None

    def do_grch38():
        if assembly != "GRCh37" or not vcf or is_mt:
            return None
        try:
            return grch38_view(vcf, gene_symbol)
        except SourceError as err:
            w_38.append(f"GRCh38 view unavailable: {err.message} (HTTP {err.status or '-'})")
        except ValueError as err:  # a refusal (REF changed between builds, split mapping): a decision
            w_38.append(f"GRCh38 view not built: {err}")
        except (KeyError, TypeError, IndexError, AttributeError) as err:
            w_38.append(f"GRCh38 view: unexpected response shape ({type(err).__name__}: {err})")
        return None

    def do_mitomap():
        if not (is_mt and vcf):
            return None
        if len(vcf[2]) != 1 or len(vcf[3]) != 1:
            w_mm.append("MITOMAP: not checked — the MITOMASTER service takes single-base changes only")
            return None
        return attempt("MITOMAP (MITOMASTER)", lambda: mitomap(vcf[1], vcf[2], vcf[3]), w_mm)

    def do_mavedb():
        if not contracts:
            return None
        if is_mt:
            w_mave.append("MaveDB functional scores: not checked for an mtDNA variant")
            return None
        if not (hgvs.get("p") or hgvs.get("c")):
            return None
        return mavedb_scores(gene_symbol, hgvs.get("p"), hgvs.get("c"), w_mave)

    def do_china():
        if not contracts:
            return None
        if is_mt:
            w_cn.append("Chinese population frequencies: not checked for an mtDNA variant (the cohorts queried "
                        "report nuclear genotypes)")
            return None
        return china_frequencies(vcf, w_cn, assembly)

    with ThreadPoolExecutor(max_workers=7) as pool:
        f_gn, f_cv, f_lv = pool.submit(do_gnomad), pool.submit(do_clinvar), pool.submit(do_litvar)
        f_38, f_mm, f_mave = pool.submit(do_grch38), pool.submit(do_mitomap), pool.submit(do_mavedb)
        f_cn = pool.submit(do_china)
        gn, cv, lv = f_gn.result(), f_cv.result(), f_lv.result()
        v38, mm, mave, cn = f_38.result(), f_mm.result(), f_mave.result(), f_cn.result()
    warnings += w_gn + w_cv + w_lv + w_38 + w_mm + w_mave + w_cn

    # population
    population: Optional[Dict[str, Any]] = None
    if gn is not None:
        sources += gn.sources
        warnings += gn.warnings
        population = gn.result
    fallback = None
    if population is None and is_mt:
        warnings.append("no mtDNA population frequency available (gnomAD's mtDNA callset did not answer); "
                        "an unanswered lookup is not absence")
    elif population is None:
        fallback = vep_frequencies(rec, alt_raw)
        if fallback:
            warnings.append("gnomAD API unavailable: frequencies below come from VEP's colocated-variant data "
                            "(no allele numbers, no coverage; BA1/BS1 need the allele number — recheck in gnomAD)")
        else:
            warnings.append("no population frequency available (gnomAD unavailable and VEP reports none)")

    # ClinVar
    clin = None
    others: List[Dict[str, Any]] = []
    if cv is not None:
        sources += cv.sources
        warnings += cv.warnings
        caid = (population or {}).get("caid")
        exact, others = match_clinvar(cv.result, vcf, assembly, window, caid)
        if exact:
            clin = clinvar.compact(exact[0])
            if len(exact) > 1:
                warnings.append(f"ClinVar: {len(exact)} records for this allele ({', '.join(r['vcv'] for r in exact)}); showing the best reviewed")
        elif vcf is None and len([r for r in cv.result if not r.get('compound')]) == 1:
            only = [r for r in cv.result if not r.get("compound")][0]
            clin = clinvar.compact(only)
            warnings.append(f"ClinVar {only['vcv']}: allele match not verified (no VCF form)")

    literature = None
    if lv is not None:
        try:
            sources += list(getattr(lv, "sources", []) or [])
            warnings += list(getattr(lv, "warnings", []) or [])
            literature = {"query": f"{gene_symbol} {lit_text}" if gene_symbol and lit_text == hgvs.get("p") else lit_text,
                          "litvar": _litvar_compact(getattr(lv, "result", lv), rsids, hgvs.get("p"), hgvs.get("c"))}
        except Exception as err:  # noqa: BLE001
            warnings.append(f"LitVar result unreadable: {err}")

    # GRCh38 view of a GRCh37 input (CP1-2): predictors VEP does not serve on GRCh37, gnomAD v4
    view38: Optional[Dict[str, Any]] = None
    if v38 is not None:
        sources += v38.sources
        warnings += v38.warnings
        view38 = v38.result
        p38 = view38.get("predictors") or {}
        filled = []
        for key in ("alphamissense", "spliceai", "revel", "cadd_phred"):
            if preds.get(key) is None and p38.get(key) is not None:
                preds[key] = p38[key]
                filled.append(key)
        if filled:
            preds["filled_from_grch38"] = filled
            preds["source"] += f"; {', '.join(filled)} from GRCh38 VEP on the mapped coordinates ({view38['vcf']['id']})"
    elif assembly == "GRCh37" and not is_mt:
        warnings.append("GRCh37 VEP serves no AlphaMissense or SpliceAI scores and the GRCh38 view could not be built "
                        "(see the warning above); for splicing use `zebra s2f`")

    # ACMG inputs
    terms = consequence.get("terms") or []
    revel = preds["revel"]
    sp_max = (preds["spliceai"] or {}).get("max")
    insil = (population or {}).get("in_silico") or {}
    if revel is None and insil.get("revel_max") not in (None, ""):
        revel = _num(insil.get("revel_max"))
        if revel is not None:
            preds["revel_gnomad"] = revel
    if sp_max is None and insil.get("spliceai_ds_max") not in (None, ""):
        sp_max = _num(insil.get("spliceai_ds_max"))
        if sp_max is not None:
            preds["spliceai_gnomad_max"] = sp_max
    # Which population answer feeds the ACMG inputs: gnomAD v4 on the mapped GRCh38 coordinates
    # when a GRCh37 input could be mapped (joint exome+genome FAF, 1.6 M alleles), else the input build's.
    freq_pop = population
    freq_note = None
    pop38 = (view38 or {}).get("population")
    if pop38 is not None and pop38.get("dataset") == "gnomad_r4":
        freq_pop = pop38
        freq_note = f"gnomAD v4 on the GRCh38 coordinates {view38['vcf']['id']} (mapped from GRCh37)"
    grp = (freq_pop or {}).get("grpmax") or {}
    site_covered: Optional[bool] = None
    coverage_text: Optional[str] = None
    mt_pop = population if (population or {}).get("mitochondrial") else None
    if mt_pop is not None:
        gac = None  # mtDNA: no nuclear allele count; homoplasmic/heteroplasmic counts below
    elif freq_pop is not None:
        site_covered = freq_pop.get("covered")
        fr = ((freq_pop.get("coverage_detail") or {}).get("fraction_over_20") or {})
        coverage_text = ", ".join(f"{k} {f:.0%} of samples at >=20x" for k, f in fr.items()) or None
        if freq_pop.get("found"):
            gac = (freq_pop.get("total") or {}).get("ac")
        else:
            gac = 0 if freq_pop.get("covered") else None
            if freq_pop.get("covered") is not True:
                warnings.append("absent from gnomAD but site coverage is low or unknown: PM2 not supported by absence alone")
    else:
        gac = None
    acmg_inputs = {
        "revel": revel,
        "alphamissense": (preds.get("alphamissense") or {}).get("score"),
        "spliceai_max": sp_max,
        "spliceai_source": ("SpliceAI precomputed scores served by Ensembl VEP"
                            + (" (GRCh38 coordinates mapped from GRCh37)" if "spliceai" in (preds.get("filled_from_grch38") or []) else "")
                            if sp_max is not None and preds.get("spliceai") else
                            ("SpliceAI DS max from gnomAD's in-silico table" if sp_max is not None else None)),
        "grpmax_af": grp.get("af") if freq_pop and not mt_pop else (fallback or {}).get("grpmax_af"),
        "grpmax_an": grp.get("an") if freq_pop and not mt_pop else None,
        "gnomad_ac": gac,
        "faf95": None if mt_pop else ((freq_pop or {}).get("faf95") or {}).get("value"),
        "site_covered": site_covered,
        "coverage_text": coverage_text,
        "consequence": list(terms),
        "is_missense": any("missense" in t for t in terms),
        "hgvs_c": hgvs.get("c"),
        "frequency_source": freq_note or (freq_pop or {}).get("dataset") or ("VEP colocated gnomAD" if fallback else None),
    }
    if is_mt:
        # set from the chromosome, never from whether gnomAD answered: with gnomAD down, the nuclear
        # frequency rules and PVS1 inputs must still stay away from an mtDNA variant
        acmg_inputs.update({"mitochondrial": True, "mt_af_hom": (mt_pop or {}).get("af_hom"),
                            "mt_af_het": (mt_pop or {}).get("af_het"),
                            "mt_max_heteroplasmy": (mt_pop or {}).get("max_heteroplasmy"), "heteroplasmy": het_level,
                            "gnomad_ac": None, "grpmax_af": None, "grpmax_an": None, "faf95": None})

    # Chinese population frequencies (W4's contract), looked up in parallel above
    china = None
    if cn is not None:
        sources += cn.sources
        warnings += cn.warnings
        china = cn.result
    functional = None
    if mave is not None:
        sources += mave.sources
        warnings += mave.warnings
        functional = {"mavedb": mave.result,
                      "note": ("MaveDB scores are raw assay readouts; their strength as PS3/BS3 evidence depends on "
                               "the assay's calibration (Brnich et al. 2019), which is not done here")}
    mito_part = None
    if is_mt:
        mito_part = {"heteroplasmy": heteroplasmy_note(het_level, mt_pop) if het_level is not None else None,
                     "mitomap": mm.result if mm is not None else None}
        if mm is not None:
            sources += mm.sources
            warnings += mm.warnings

    other_assembly = (population or {}).get("liftover")
    label_parts = [hgvs.get("c") or text]
    if hgvs.get("p"):
        label_parts.append(f"({hgvs['p']})")
    if gene_symbol:
        label_parts.append(gene_symbol)
    if vcf:
        label_parts.append(f"| {assembly} {vcf[0]}-{vcf[1]}-{vcf[2]}-{vcf[3]}")
    result = {
        "variant": " ".join(label_parts),
        "input": text,
        "kind": kind,
        "assembly": assembly,
        "vcf": {"assembly": assembly, "chrom": vcf[0], "pos": vcf[1], "ref": vcf[2], "alt": vcf[3],
                "id": f"{vcf[0]}-{vcf[1]}-{vcf[2]}-{vcf[3]}"} if vcf else None,
        "other_assembly": other_assembly,
        "rsids": rsids,
        "caid": (population or {}).get("caid"),
        "gene": gene_symbol,
        "transcript": transcript,
        "hgvs": hgvs,
        "consequence": consequence,
        "predictors": preds,
        "population": population if population is not None else ({"fallback": fallback} if fallback else None),
        "clinvar": clin,
        "clinvar_others_at_locus": [{"vcv": r["vcv"], "title": r["title"], "classification": r["classification"],
                                     "stars": r["stars"], "url": r["url"]} for r in others[:5]],
        "literature": literature,
        "acmg_inputs": acmg_inputs,
    }
    if view38 is not None:
        result["grch38"] = view38
        result["builds"] = {"GRCh37": result["vcf"], "GRCh38": dict(view38["vcf"], mapped_by=view38["mapped_by"])}
    if china is not None:
        result["population_china"] = china
    if functional is not None:
        result["functional_scores"] = functional
    if mito_part is not None:
        result["mitochondrial"] = mito_part
    if cv is not None and clin is None:
        result["clinvar_note"] = "no ClinVar record for this exact allele" + (
            f" ({len(others)} other record(s) at the locus)" if others else "")
    return Outcome(result, sources=sources, warnings=warnings)


# ------------------------------------------------------------ PVS1 / PS1 / PM5 inputs (CP1-3)

def _p_part(hgvsp: Optional[str]) -> Optional[str]:
    return hgvsp.split(":", 1)[1] if hgvsp and ":" in hgvsp else hgvsp


def judgement_inputs(r: Dict[str, Any]) -> Outcome:
    """What the PVS1 decision tree and PS1/PM5 ask about, computed for one variant card.

    Structure, never codes: for a null variant, where the premature stop falls on the exon structure
    (NMD prediction), how much of the protein follows it, the canonical splice site's exon and its
    frame, and the loss-of-function mechanism evidence (ClinGen HI score, gnomAD LOEUF/pLI); for a
    missense variant, the ClinVar Pathogenic/Likely pathogenic records at the same codon (same change
    -> PS1 input, different change -> PM5 input) with their review stars.
    """
    from zebra import acmg
    from zebra.sources import clingen  # gene-level ClinGen dosage (HI score)

    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []
    terms = list((r.get("consequence") or {}).get("terms") or [])
    gene = r.get("gene")
    asm = r.get("assembly") or "GRCh38"
    vcf = r.get("vcf") or {}
    tx_id = (r.get("transcript") or {}).get("ensembl")
    hgvsp = _p_part((r.get("hgvs") or {}).get("p"))
    v38 = r.get("grch38")
    intron = (r.get("consequence") or {}).get("intron")
    if v38 and (v38.get("transcript") or {}).get("ensembl"):
        # the GRCh38 (MANE) transcript is used as a whole: its exon structure, consequence terms,
        # intron number and protein change; never one transcript's structure with another's numbering
        asm, vcf, tx_id = "GRCh38", v38["vcf"], v38["transcript"]["ensembl"]
        hgvsp = _p_part(v38.get("hgvsp"))
        terms = list(v38.get("consequence_terms") or terms)
        intron = v38.get("intron")
    null = any(t in acmg.NULL_CONSEQUENCES for t in terms)
    missense = any("missense" in t for t in terms)
    out: Dict[str, Any] = {}
    if (r.get("acmg_inputs") or {}).get("mitochondrial"):
        return Outcome({"note": "PVS1/PS1/PM5 inputs are not computed for mtDNA variants (McCormick et al. 2020 "
                                "specify their own rules)"})
    if not (null or missense):
        return Outcome({"note": f"no PVS1 or PS1/PM5 inputs for this variant class ({', '.join(terms) or '-'})"})

    tasks: Dict[str, Tuple[str, Callable[[], Outcome]]] = {}
    if null and tx_id:
        tasks["structure"] = ("Ensembl transcript structure", lambda: ensembl.transcript(tx_id, asm))
    if null and gene:
        tasks["dosage"] = ("ClinGen dosage (HI score)", lambda: clingen.dosage(gene, "GRCh38"))
        tasks["constraint"] = ("gnomAD constraint", lambda: gnomad.gene_constraint(gene, "GRCh38"))
    if missense and vcf.get("pos"):
        lo, hi = int(vcf["pos"]) - 2, int(vcf["pos"]) + 2
        tasks["codon"] = ("ClinVar records at the codon",
                          lambda: clinvar.at_positions(str(vcf["chrom"]), lo, hi, asm))
    got: Dict[str, Any] = {}

    def go(key: str):
        w: List[str] = []
        return key, attempt(tasks[key][0], tasks[key][1], w), w

    with ThreadPoolExecutor(max_workers=max(1, len(tasks))) as pool:
        for key, res, w in pool.map(go, list(tasks)):
            warnings += w
            if res is not None:
                sources += res.sources
                warnings += res.warnings
                got[key] = res.result

    if null:
        pv: Dict[str, Any] = {"variant_class": [t for t in terms if t in acmg.NULL_CONSEQUENCES],
                              "transcript": tx_id, "assembly": asm,
                              "decision_tree": ("Abou Tayoun et al. 2018 (Hum Mutat 39:1517) PVS1 decision tree. "
                                                "These are its inputs, not a PVS1 call: the strength (VeryStrong, "
                                                "Strong, Moderate, Supporting or not applicable) is the skill's "
                                                "judgement.")}
        st = got.get("structure")
        if st and st.get("translation") and st["translation"].get("start"):
            strand = int(st["strand"])
            exons = acmg.transcript_exons(st["exons"], strand)
            tl = st["translation"]
            cds_start_g = tl["start"] if strand == 1 else tl["end"]
            cds_end_g = tl["end"] if strand == 1 else tl["start"]
            cds = (acmg.cdna_position(exons, cds_start_g, strand), acmg.cdna_position(exons, cds_end_g, strand))
            pv["exon_count"] = len(exons)
            pv["protein_length"] = tl.get("length")
            first, ptc, how = acmg.ptc_codon_from_hgvsp(hgvsp, terms)
            pv["premature_stop"] = {"read_as": how, "first_altered_residue": first, "stop_codon": ptc}
            if first and tl.get("length"):
                pv["premature_stop"]["fraction_of_protein_from_first_altered_residue"] = round(
                    max(0, tl["length"] - first + 1) / tl["length"], 3)
            if ptc:
                try:
                    pv["nmd"] = acmg.nmd_inputs(exons, strand, cds_start_g, ptc, tl.get("length"))
                except ValueError as err:
                    warnings.append(f"NMD prediction not computed: {err}")
            elif any(t in ("stop_gained", "frameshift_variant") for t in terms):
                pv["nmd"] = {"nmd_predicted": None, "position": f"not computed: {how}", "rule": acmg.NMD_RULE}
            splice = [t for t in terms if t in ("splice_donor_variant", "splice_acceptor_variant")]
            if splice and intron and None not in cds:
                try:
                    n = int(str(intron).split("/")[0])
                    k = n if splice[0] == "splice_donor_variant" else n + 1
                    ex = exons[k - 1]
                    coding = acmg.coding_overlap(ex, cds)  # type: ignore[arg-type]
                    pv["canonical_splice_site"] = {
                        "site": "donor" if splice[0] == "splice_donor_variant" else "acceptor",
                        "intron": intron, "exon_affected": f"{k}/{len(exons)}", "exon_length": ex["length"],
                        "exon_coding_length": coding,
                        "skipping_keeps_frame": (coding % 3 == 0) if coding else None,
                        "note": ("if the exon is skipped: an in-frame skip removes coding_length/3 residues; an "
                                 "out-of-frame skip shifts the frame (NMD depends on where the new stop falls). "
                                 "Cryptic sites nearby can rescue — see the SpliceAI/Pangolin result"),
                    }
                except (ValueError, IndexError) as err:
                    warnings.append(f"canonical splice site inputs not computed: {err}")
            if "start_lost" in terms:
                pv["start_lost"] = ("PVS1 start-loss branch: an in-frame downstream methionine and pathogenic "
                                    "variants upstream of it decide the strength (Abou Tayoun 2018); not computed here")
        elif null:
            pv["structure_note"] = "transcript exon structure unavailable: NMD and splice-site inputs not computed"
        dos = got.get("dosage")
        con = got.get("constraint")
        hi = (dos or {}).get("haploinsufficiency") or {}
        pv["lof_mechanism"] = {
            "clingen_hi_score": hi.get("score") if dos else None,
            "clingen_hi_description": hi.get("description") if dos else None,
            "clingen_url": (dos or {}).get("url"),
            "gnomad_loeuf": (con or {}).get("loeuf"), "gnomad_pli": (con or {}).get("pLI"),
            "gnomad_version": (con or {}).get("version"),
            "note": acmg.PVS1_MECHANISM_NOTE,
        }
        if "dosage" in tasks and dos is None and "dosage" in got:
            pv["lof_mechanism"]["clingen_note"] = "no ClinGen dosage curation for this gene"
        out["pvs1_inputs"] = pv
    if missense:
        recs = got.get("codon")
        if recs is not None and hgvsp:
            try:
                own_spdi = None
                if len(str(vcf.get("ref") or "")) == 1 and len(str(vcf.get("alt") or "")) == 1:
                    own_spdi = f":{int(vcf['pos']) - 1}:{vcf['ref']}:{vcf['alt']}"  # SPDI is 0-based
                cr = acmg.codon_records(recs, hgvsp, exclude_vcv=(r.get("clinvar") or {}).get("vcv"), gene=gene,
                                        exclude_c=(r.get("hgvs") or {}).get("c"), exclude_spdi=own_spdi)
                cr["searched"] = (f"ClinVar records within {vcf['chrom']}:{int(vcf['pos']) - 2}-{int(vcf['pos']) + 2} "
                                  f"({asm}); a codon split by an intron is only partly covered")
                out["ps1_pm5_inputs"] = cr
            except ValueError as err:
                warnings.append(f"PS1/PM5 inputs not computed: {err}")
    return Outcome(out, sources=sources, warnings=warnings)
