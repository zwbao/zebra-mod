"""Ensembl REST: VEP (with AlphaMissense, REVEL, CADD, SpliceAI), sequence, ID recoding, lookup.

GRCh38: https://rest.ensembl.org   GRCh37: https://grch37.rest.ensembl.org
No key. ~15 requests/s; POST endpoints take up to 200 variants.
"""

from __future__ import annotations

import re
import urllib.parse
from typing import Any, Dict, List, Optional, Sequence, Tuple

from zebra.core import Outcome
from zebra.http import SourceError, get_json, post_json, request
from zebra.sources import record as source_record

HOSTS = {"GRCh38": "https://rest.ensembl.org", "GRCh37": "https://grch37.rest.ensembl.org"}
# plugin flags as the REST service names them; GRCh37 serves fewer plugins
VEP_FLAGS_38 = {
    "canonical": 1, "mane": 1, "hgvs": 1, "numbers": 1, "protein": 1, "domains": 0, "variant_class": 1,
    "AlphaMissense": 1, "REVEL": 1, "CADD": 1, "SpliceAI": 1, "pick_order": "mane_select,canonical",
}
VEP_FLAGS_37 = {"canonical": 1, "hgvs": 1, "numbers": 1, "protein": 1, "variant_class": 1, "CADD": 1, "REVEL": 1}
VCF_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT|M)[-:_\s](\d+)[-:_\s]([ACGTNacgtn]+)[-:_>\s]([ACGTNacgtn]+)$")
RSID_RE = re.compile(r"^rs\d+$", re.I)
# Anchored at both ends (B-P2-7): an unanchored pattern accepted
# "NM_000492.4:c.1521_1523del --assembly GRCh37" as one argv string and silently
# ignored the trailing text, annotating on the default build.
ACCESSION = r"(N[MRCGP]_\d+(?:\.\d+)?|ENS[TGP]\d+(?:\.\d+)?|LRG_\d+(?:t\d+)?)"
HGVS_RE = re.compile(rf"^{ACCESSION}(?:\([A-Za-z0-9-]+\))?:[cgnmrp]\.\S+$", re.I)
GENE_HGVS_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9.@-]*:[cp]\.\S+$")
# Ensembl's message when a RefSeq accession carries a version the current
# release does not hold (verified live: NM_001165963.1 fails, NM_001165963 works).
NO_TRANSCRIPT_RE = re.compile(r"Could not get a Transcript object|Could not fetch a Transcript", re.I)


def host(assembly: str) -> str:
    if assembly not in HOSTS:
        raise ValueError("assembly must be GRCh38 or GRCh37")
    return HOSTS[assembly]


# Mitochondrial results are reported in m. notation (m.3243A>G, m.8344A>G), the
# form every laboratory letter and patient group uses. It is unambiguous — the
# mitochondrial genome has one coordinate system — so it is read as a position on MT
# rather than sent to VEP, which has no m. reference sequence for it.
MITO_RE = re.compile(r"^m\.(\d+)\s*([ACGTacgt]+)\s*(?:>|->|/)\s*([ACGTacgt]+)$", re.I)
MITO_DEL_RE = re.compile(r"^m\.(\d+)(?:_(\d+))?\s*del([ACGTacgt]*)$", re.I)


def parse_vcf_like(text: str) -> Optional[Tuple[str, int, str, str]]:
    t = text.strip()
    m = MITO_RE.match(t)
    if m:
        return "MT", int(m.group(1)), m.group(2).upper(), m.group(3).upper()
    m = VCF_RE.match(t)
    if not m:
        return None
    chrom = m.group(1).upper()
    chrom = "MT" if chrom == "M" else chrom
    return chrom, int(m.group(2)), m.group(3).upper(), m.group(4).upper()


def classify_input(text: str) -> str:
    t = text.strip()
    if parse_vcf_like(t):
        return "vcf"
    if MITO_DEL_RE.match(t):
        return "mito_indel"
    if RSID_RE.match(t):
        return "rsid"
    if HGVS_RE.match(t):
        return "hgvs"
    if GENE_HGVS_RE.match(t):
        return "gene_hgvs"  # e.g. SCN1A:c.2134C>T — VEP resolves via the gene's canonical transcript
    return "unknown"


def strip_version(accession: str) -> str:
    """`NM_001165963.1` -> `NM_001165963`; an accession with no version is returned unchanged."""
    return accession.split(".")[0]


def _vep_params(assembly: str) -> Dict[str, Any]:
    return dict(VEP_FLAGS_38 if assembly == "GRCh38" else VEP_FLAGS_37)


def _vep_hgvs(base: str, text: str, params: Dict[str, Any]):
    return get_json(f"{base}/vep/human/hgvs/{urllib.parse.quote(text, safe=':>()')}",
                    source="Ensembl VEP", params=params, cache_ttl=14 * 86400, timeout=60)


def vep(variant: str, assembly: str = "GRCh38") -> Outcome:
    """VEP for one variant given as HGVS, rsID, or chrom-pos-ref-alt.

    A legacy transcript version is resolved, not refused (F37). Ensembl answers
    an accession whose version the current release does not hold with
    "Could not get a Transcript object" (verified live for NM_001165963.1,
    while NM_001165963 without the version is accepted); the call is retried
    without the version and the answer carries a warning naming the version VEP
    actually used, because c. numbering can differ between versions.
    """
    kind = classify_input(variant)
    base = host(assembly)
    params = _vep_params(assembly)
    warnings: List[str] = []
    if kind == "vcf":
        chrom, pos, ref, alt = parse_vcf_like(variant)  # type: ignore[misc]
        resp = post_json(f"{base}/vep/human/region", {"variants": [vcf_line(chrom, pos, ref, alt)], **params},
                         source="Ensembl VEP", cache_ttl=14 * 86400, timeout=60)
    elif kind == "rsid":
        resp = get_json(f"{base}/vep/human/id/{urllib.parse.quote(variant.strip(), safe='')}", source="Ensembl VEP",
                        params=params, cache_ttl=14 * 86400, timeout=60)
    elif kind in ("hgvs", "gene_hgvs"):
        text = variant.strip()
        try:
            resp = _vep_hgvs(base, text, params)
        except SourceError as err:
            acc, sep, change = text.partition(":")
            bare = strip_version(acc)
            if not (sep and bare != acc and NO_TRANSCRIPT_RE.search(err.message or "")):
                raise
            resp = _vep_hgvs(base, f"{bare}{sep}{change}", params)
            warnings.append(f"{acc} is not in this Ensembl release; resolved {bare}{sep}{change} instead. "
                            "c. numbering can differ between transcript versions — check the reference base and "
                            "the position against the version your report used "
                            "(VariantValidator or Mutalyzer map between versions)")
    else:
        raise ValueError(f"cannot read {variant!r}: use transcript HGVS (NM_...:c.), an rsID, or chrom-pos-ref-alt "
                         "(one value only — flags such as --assembly must be separate arguments)")
    data = resp.json()
    if not isinstance(data, list) or not data:
        raise ValueError(f"VEP returned nothing for {variant!r}")
    if len(data) > 1:
        alleles = [d.get("allele_string") for d in data]
        warnings.append(f"{variant} maps to {len(data)} alleles ({', '.join(map(str, alleles))}); using {alleles[0]} — give HGVS or chrom-pos-ref-alt to pick one")
    return Outcome(data[0], warnings=warnings,
                   sources=[source_record("Ensembl VEP", variant, resp, note=f"{assembly}; plugins AlphaMissense/REVEL/CADD/SpliceAI where served")])


def vcf_line(chrom: str, pos: int, ref: str, alt: str) -> str:
    """VEP's region input: VCF columns, with the shared leading base kept as VCF does."""
    return f"{chrom} {pos} . {ref} {alt} . . ."


def vep_batch(variants: Sequence[Tuple[str, int, str, str]], assembly: str = "GRCh38") -> Outcome:
    """VEP for up to thousands of VCF-style variants, 200 per request."""
    base = host(assembly)
    params = _vep_params(assembly)
    out: List[Dict[str, Any]] = []
    sources = []
    for i in range(0, len(variants), 200):
        chunk = variants[i:i + 200]
        resp = post_json(f"{base}/vep/human/region", {"variants": [vcf_line(*v) for v in chunk], **params},
                         source="Ensembl VEP", cache_ttl=14 * 86400, timeout=180)
        out.extend(resp.json())
        sources.append(source_record("Ensembl VEP", f"batch {i // 200 + 1} ({len(chunk)} variants)", resp, note=assembly))
    return Outcome(out, sources=sources)


def pick_transcript(rec: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """The transcript to annotate on: one the variant overlaps, MANE Select preferred.

    Delegates to `zebra.sources.variant.pick_transcript`, so the editing and
    triage paths choose the same transcript as the variant card. Before this,
    ordering MANE -> canonical -> protein_coding over the whole consequence list
    could return a transcript the variant lies *outside* (VEP reports one for
    every nearby gene), which is how m.3243A>G was annotated on MT-ND1.
    """
    from zebra.sources.variant import pick_transcript as _pick

    return _pick(rec.get("transcript_consequences") or [], None, rec.get("most_severe_consequence"))


def spliceai_max(tc: Optional[Dict[str, Any]]) -> Optional[float]:
    if not tc or not isinstance(tc.get("spliceai"), dict):
        return None
    s = tc["spliceai"]
    vals = [s.get(k) for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL") if isinstance(s.get(k), (int, float))]
    return max(vals) if vals else None


def vep_allele(ref: str, alt: str) -> str:
    """The allele as VEP keys colocated frequencies: shared leading bases trimmed, '-' when nothing is left."""
    i = 0
    while i < min(len(ref), len(alt)) and ref[i] == alt[i]:
        i += 1
    rest = alt[i:]
    return rest if rest else "-"


def frequencies(rec: Dict[str, Any], alt: Optional[str] = None, ref: Optional[str] = None) -> Dict[str, Any]:
    """gnomAD frequencies VEP reports for the colocated known variant (exomes 'gnomade_*', genomes 'gnomadg_*').

    With `ref` and `alt` only that exact allele is accepted (a multi-allelic site's other alleles are ignored);
    with `alt` alone it must equal the VEP allele key; with neither, the best-populated entry is taken.
    """
    want = vep_allele(ref, alt) if (ref and alt) else alt
    best: Dict[str, Any] = {}
    for cv in rec.get("colocated_variants") or []:
        freqs = cv.get("frequencies") or {}
        for allele, groups in freqs.items():
            if want is not None and allele != want:
                continue
            if groups and len(groups) > len(best.get("groups", {})):
                best = {"id": cv.get("id"), "allele": allele, "groups": groups,
                        "clin_sig": cv.get("clin_sig"), "pubmed": cv.get("pubmed")}
    return best


def sequence(chrom: str, start: int, end: int, assembly: str = "GRCh38", strand: int = 1) -> Outcome:
    """Reference sequence, 1-based inclusive."""
    if end < start or end - start > 1_000_000:
        raise ValueError("bad region")
    region = f"{chrom}:{start}..{end}:{strand}"
    resp = request(f"{host(assembly)}/sequence/region/human/{region}", source="Ensembl sequence",
                   accept="text/plain", cache_ttl=90 * 86400)
    return Outcome(resp.text.strip().upper(), sources=[source_record("Ensembl sequence", f"{assembly} {region}", resp)])


def recode(variant: str, assembly: str = "GRCh38") -> Outcome:
    """variant_recoder: every equivalent notation (HGVS g./c./p., SPDI, VCF strings, rsIDs)."""
    resp = get_json(f"{host(assembly)}/variant_recoder/human/{urllib.parse.quote(variant.strip(), safe=':>()')}",
                    source="Ensembl variant_recoder", params={"vcf_string": 1}, cache_ttl=30 * 86400, timeout=60)
    return Outcome(resp.json(), sources=[source_record("Ensembl variant_recoder", variant, resp)])


def lookup_symbol(symbol: str, assembly: str = "GRCh38", expand: bool = False) -> Outcome:
    resp = get_json(f"{host(assembly)}/lookup/symbol/homo_sapiens/{urllib.parse.quote(symbol, safe='')}", source="Ensembl lookup",
                    params={"expand": 1 if expand else 0}, cache_ttl=30 * 86400)
    return Outcome(resp.json(), sources=[source_record("Ensembl lookup", symbol, resp)])
