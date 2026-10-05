"""Ensembl REST: VEP (with AlphaMissense, REVEL, CADD, SpliceAI), sequence, ID recoding, lookup.

GRCh38: https://rest.ensembl.org   GRCh37: https://grch37.rest.ensembl.org
No key. ~15 requests/s; POST endpoints take up to 200 variants.
"""

from __future__ import annotations

import re
import urllib.parse
from typing import Any, Dict, List, Optional, Sequence, Tuple

from zebra.core import Outcome
from zebra.http import get_json, post_json, request, source_record

HOSTS = {"GRCh38": "https://rest.ensembl.org", "GRCh37": "https://grch37.rest.ensembl.org"}
# plugin flags as the REST service names them; GRCh37 serves fewer plugins
VEP_FLAGS_38 = {
    "canonical": 1, "mane": 1, "hgvs": 1, "numbers": 1, "protein": 1, "domains": 0, "variant_class": 1,
    "AlphaMissense": 1, "REVEL": 1, "CADD": 1, "SpliceAI": 1, "pick_order": "mane_select,canonical",
}
VEP_FLAGS_37 = {"canonical": 1, "hgvs": 1, "numbers": 1, "protein": 1, "variant_class": 1, "CADD": 1, "REVEL": 1}
VCF_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT|M)[-:_\s](\d+)[-:_\s]([ACGTNacgtn]+)[-:_>\s]([ACGTNacgtn]+)$")
RSID_RE = re.compile(r"^rs\d+$", re.I)
HGVS_RE = re.compile(r"^(N[MRCGP]_\d+(?:\.\d+)?|ENS[TGP]\d+(?:\.\d+)?|LRG_\d+(?:t\d+)?)(?:\([A-Za-z0-9-]+\))?:[cgnmrp]\.", re.I)


def host(assembly: str) -> str:
    if assembly not in HOSTS:
        raise ValueError("assembly must be GRCh38 or GRCh37")
    return HOSTS[assembly]


def parse_vcf_like(text: str) -> Optional[Tuple[str, int, str, str]]:
    m = VCF_RE.match(text.strip())
    if not m:
        return None
    chrom = m.group(1).upper()
    chrom = "MT" if chrom == "M" else chrom
    return chrom, int(m.group(2)), m.group(3).upper(), m.group(4).upper()


def classify_input(text: str) -> str:
    t = text.strip()
    if parse_vcf_like(t):
        return "vcf"
    if RSID_RE.match(t):
        return "rsid"
    if HGVS_RE.match(t):
        return "hgvs"
    if re.match(r"^[A-Za-z0-9-]+:[cp]\.", t):
        return "gene_hgvs"  # e.g. SCN1A:c.2134C>T — VEP resolves via the gene's canonical transcript
    return "unknown"


def _vep_params(assembly: str) -> Dict[str, Any]:
    return dict(VEP_FLAGS_38 if assembly == "GRCh38" else VEP_FLAGS_37)


def vep(variant: str, assembly: str = "GRCh38") -> Outcome:
    """VEP for one variant given as HGVS, rsID, or chrom-pos-ref-alt."""
    kind = classify_input(variant)
    base = host(assembly)
    params = _vep_params(assembly)
    if kind == "vcf":
        chrom, pos, ref, alt = parse_vcf_like(variant)  # type: ignore[misc]
        resp = post_json(f"{base}/vep/human/region", {"variants": [vcf_line(chrom, pos, ref, alt)], **params},
                         source="Ensembl VEP", cache_ttl=14 * 86400, timeout=60)
    elif kind == "rsid":
        resp = get_json(f"{base}/vep/human/id/{urllib.parse.quote(variant.strip())}", source="Ensembl VEP",
                        params=params, cache_ttl=14 * 86400, timeout=60)
    elif kind in ("hgvs", "gene_hgvs"):
        resp = get_json(f"{base}/vep/human/hgvs/{urllib.parse.quote(variant.strip(), safe=':>()')}",
                        source="Ensembl VEP", params=params, cache_ttl=14 * 86400, timeout=60)
    else:
        raise ValueError(f"cannot read {variant!r}: use transcript HGVS (NM_...:c.), an rsID, or chrom-pos-ref-alt")
    data = resp.json()
    if not isinstance(data, list) or not data:
        raise ValueError(f"VEP returned nothing for {variant!r}")
    warnings = []
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
    """MANE Select, else canonical, else the first protein-coding consequence."""
    tcs = rec.get("transcript_consequences") or []
    for pred in (lambda t: t.get("mane_select"), lambda t: t.get("canonical"),
                 lambda t: t.get("biotype") == "protein_coding", lambda t: True):
        for t in tcs:
            if pred(t):
                return t
    return None


def spliceai_max(tc: Optional[Dict[str, Any]]) -> Optional[float]:
    if not tc or not isinstance(tc.get("spliceai"), dict):
        return None
    s = tc["spliceai"]
    vals = [s.get(k) for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL") if isinstance(s.get(k), (int, float))]
    return max(vals) if vals else None


def frequencies(rec: Dict[str, Any], alt: Optional[str] = None) -> Dict[str, Any]:
    """gnomAD frequencies VEP reports for the colocated known variant (exomes 'gnomade_*', genomes 'gnomadg_*')."""
    best: Dict[str, Any] = {}
    for cv in rec.get("colocated_variants") or []:
        freqs = cv.get("frequencies") or {}
        for allele, groups in freqs.items():
            if alt and allele not in (alt, "-") and len(alt) == 1 and allele != alt:
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
    resp = get_json(f"{host(assembly)}/lookup/symbol/homo_sapiens/{urllib.parse.quote(symbol)}", source="Ensembl lookup",
                    params={"expand": 1 if expand else 0}, cache_ttl=30 * 86400)
    return Outcome(resp.json(), sources=[source_record("Ensembl lookup", symbol, resp)])
