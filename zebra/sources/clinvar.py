"""ClinVar through NCBI E-utilities (esearch + esummary).

No key needed (3 requests/s); NCBI_API_KEY in the environment raises it to 10/s.
Records are found by VariationID (VCV), rsID, HGVS name, canonical SPDI or
genomic position; HGVS text search is fuzzy, so callers that know the exact
allele should match the returned records on SPDI (see zebra.sources.variant).
"""

from __future__ import annotations

import os
import re
from typing import Any, Dict, Iterable, List, Optional

from zebra.core import Outcome
from zebra.http import get_json, source_record

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CACHE_TTL = 7 * 86400  # ClinVar releases weekly
VCV_URL = "https://www.ncbi.nlm.nih.gov/clinvar/variation/{}/"

# ClinVar review status -> gold stars (https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, multiple submitters": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "no classifications from unflagged records": 0,
    "no classification for the individual variant": 0,
}
COMPOUND_TYPES = {"haplotype", "compoundheterozygote", "diplotype", "distinct chromosomes", "phase unknown"}
VCV_RE = re.compile(r"^(?:VCV)?0*(\d+)(?:\.\d+)?$", re.I)


def _params(extra: Dict[str, Any]) -> Dict[str, Any]:
    p = {"db": "clinvar", "retmode": "json", "tool": "zebra-mod"}
    key = os.environ.get("NCBI_API_KEY")
    if key:
        p["api_key"] = key
    p.update(extra)
    return p


def stars(review_status: Optional[str]) -> Optional[int]:
    if not review_status:
        return None
    return STARS.get(review_status.strip().lower())


def vcv_uid(vcv: str) -> str:
    m = VCV_RE.match(vcv.strip())
    if not m:
        raise ValueError(f"{vcv!r} is not a ClinVar VariationID / VCV accession")
    return m.group(1)


def search(term: str, retmax: int = 40) -> Outcome:
    """esearch: ClinVar VariationIDs (uids) for an Entrez query."""
    params = _params({"term": term, "retmax": retmax})
    resp = get_json(f"{EUTILS}/esearch.fcgi", source="ClinVar", params=params, cache_ttl=CACHE_TTL)
    res = resp.json().get("esearchresult") or {}
    if "ERROR" in res and resp.cached:  # do not keep serving a cached server-side error
        resp = get_json(f"{EUTILS}/esearch.fcgi", source="ClinVar", params=params, cache_ttl=0)
        res = resp.json().get("esearchresult") or {}
    if "ERROR" in res:
        raise ValueError(f"ClinVar esearch error: {res['ERROR']}")
    ids = list(res.get("idlist") or [])
    return Outcome({"term": term, "count": int(res.get("count") or 0), "ids": ids},
                   sources=[source_record("ClinVar esearch", term, resp)])


def term_for(*, vcv: Optional[Iterable[str]] = None, rsid: Optional[str] = None, hgvs: Optional[str] = None,
             spdi: Optional[str] = None, chrom: Optional[str] = None, pos: Optional[int] = None,
             assembly: str = "GRCh38") -> str:
    """One OR-ed Entrez query for every handle we hold (one request instead of four)."""
    parts: List[str] = []
    for v in vcv or []:
        parts.append(f"{vcv_uid(v)}[VID]")
    if rsid:
        parts.append(f"{rsid.lower()}[VRID]")
    if hgvs:
        parts.append(f'"{hgvs}"[VRNM]')
    if spdi:
        parts.append(f'"{spdi}"[CSPDI]')
    if chrom and pos:
        field = "CPOS" if assembly == "GRCh38" else "C37"
        parts.append(f"({chrom}[CHR] AND {int(pos)}[{field}])")
    if not parts:
        raise ValueError("nothing to search ClinVar with")
    return " OR ".join(parts)


def summaries(uids: List[str]) -> Outcome:
    """esummary for up to 40 VariationIDs, parsed."""
    uids = [u for u in dict.fromkeys(str(u) for u in uids)][:40]
    if not uids:
        return Outcome([])
    resp = get_json(f"{EUTILS}/esummary.fcgi", source="ClinVar", params=_params({"id": ",".join(uids)}),
                    cache_ttl=CACHE_TTL, timeout=60)
    res = resp.json().get("result") or {}
    out = [parse_summary(res[u]) for u in res.get("uids", []) if u in res and not res[u].get("error")]
    return Outcome(out, sources=[source_record("ClinVar", ",".join(f"VCV{int(u):09d}" for u in uids), resp)])


def lookup(*, vcv: Optional[str] = None, rsid: Optional[str] = None, hgvs: Optional[str] = None) -> Outcome:
    """Records for a VCV id, an rsID or an HGVS name (several when ambiguous; the caller picks)."""
    if vcv:
        out = summaries([vcv_uid(vcv)])
        return out
    term = term_for(rsid=rsid, hgvs=hgvs)
    found = search(term)
    out = summaries(found.result["ids"])
    out.sources = found.sources + out.sources
    if found.result["count"] > len(found.result["ids"]):
        out.warnings.append(f"ClinVar: {found.result['count']} records match {term}; first {len(found.result['ids'])} returned")
    return out


def _date(s: Optional[str]) -> Optional[str]:
    if not s or s.startswith("1/01/01"):
        return None
    d = s.split(" ")[0].replace("/", "-")
    return d or None


def parse_summary(rec: Dict[str, Any]) -> Dict[str, Any]:
    uid = str(rec.get("uid"))
    germ = rec.get("germline_classification") or {}
    vset = rec.get("variation_set") or []
    first = vset[0] if vset else {}
    xrefs = first.get("variation_xrefs") or []
    rsid = next((f"rs{x['db_id']}" for x in xrefs if x.get("db_source") == "dbSNP"), None)
    caid = next((x["db_id"] for x in xrefs if x.get("db_source") == "ClinGen"), None)
    locs = {}
    for loc in first.get("variation_loc") or []:
        if loc.get("assembly_name") in ("GRCh38", "GRCh37") and loc.get("status") in ("current", "previous", None):
            locs[loc["assembly_name"]] = {"chr": loc.get("chr"), "start": _int(loc.get("start")), "stop": _int(loc.get("stop"))}
    conditions = []
    for t in (germ.get("trait_set") or [])[:12]:
        ids = [f"{x['db_source']}:{x['db_id']}" if not str(x.get("db_id", "")).startswith(x.get("db_source", "") + ":")
               else x["db_id"] for x in t.get("trait_xrefs") or [] if x.get("db_source") in ("MONDO", "OMIM", "Orphanet", "MedGen", "HP")]
        conditions.append({"name": t.get("trait_name"), "ids": ids})
    review = germ.get("review_status") or None
    subs = rec.get("supporting_submissions") or {}
    return {
        "variation_id": uid,
        "vcv": rec.get("accession") or f"VCV{int(uid):09d}",
        "version": rec.get("accession_version"),
        "url": VCV_URL.format(uid),
        "title": rec.get("title"),
        "type": rec.get("obj_type"),
        "compound": (str(rec.get("obj_type") or "").replace(" ", "").lower() in COMPOUND_TYPES) or len(vset) > 1,
        "classification": germ.get("description") or None,
        "review_status": review,
        "stars": stars(review),
        "last_evaluated": _date(germ.get("last_evaluated")),
        "conditions": conditions,
        "submissions_scv": len(subs.get("scv") or []),
        "rcv_count": len(subs.get("rcv") or []),
        "genes": [g.get("symbol") for g in rec.get("genes") or [] if g.get("symbol")],
        "canonical_spdi": first.get("canonical_spdi") or None,
        "cdna_change": first.get("cdna_change") or None,
        "protein_change": rec.get("protein_change") or None,
        "locations": locs,
        "rsid": rsid,
        "caid": caid,
        "molecular_consequences": rec.get("molecular_consequence_list") or [],
    }


def _int(x: Any) -> Optional[int]:
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def compact(rec: Dict[str, Any]) -> Dict[str, Any]:
    """The fields a variant card shows."""
    keep = ("vcv", "url", "title", "classification", "review_status", "stars", "last_evaluated", "conditions",
            "submissions_scv", "rcv_count", "canonical_spdi", "rsid", "caid")
    out = {k: rec.get(k) for k in keep}
    out["conditions"] = (rec.get("conditions") or [])[:6]
    return out
