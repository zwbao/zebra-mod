"""Open Targets Platform GraphQL: diseases and targets → drugs and clinical candidates, tractability.

API: https://api.platform.opentargets.org/api/v4/graphql (no key). Schema as
verified live (API 26.9, data 26.09): `knownDrugs` is gone; drugs come from
`drugAndClinicalCandidates { rows { maxClinicalStage drug clinicalReports } }`
on Disease and on Target. Clinical reports carry their own sources (FDA label,
EMA EPAR, PMDA, ClinicalTrials.gov, TTD, USAN, DailyMed). Disease ids are
MONDO_/EFO_/Orphanet_/HP_ (underscore); OMIM and ORPHA numbers are only
cross-references here, not ids you can query by.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Sequence

from zebra.core import Outcome
from zebra.http import SourceError, post_json, source_record

API = "https://api.platform.opentargets.org/api/v4/graphql"
WEB = "https://platform.opentargets.org"
ENSG_RE = re.compile(r"^ENSG\d{11}$", re.I)
DISEASE_ID_RE = re.compile(r"^(MONDO|EFO|Orphanet|HP|DOID|OTAR|GO|NCIT)[_:]\d+$", re.I)
_PREFIX_CASE = {"mondo": "MONDO", "efo": "EFO", "orphanet": "Orphanet", "hp": "HP", "doid": "DOID", "otar": "OTAR",
                "go": "GO", "ncit": "NCIT"}
# highest first; anything else ranks last
STAGE_ORDER = ["APPROVAL", "PREAPPROVAL", "PHASE_4", "PHASE_3", "PHASE_2_3", "PHASE_2", "PHASE_1_2", "PHASE_1",
               "EARLY_PHASE_1", "PHASE_0", "PRECLINICAL", "UNKNOWN"]
MODALITIES = {"SM": "small_molecule", "AB": "antibody", "PR": "protac", "OC": "other_modalities"}

_DRUG_FIELDS = """
  maxClinicalStage
  drug { id name drugType maximumClinicalStage
         mechanismsOfAction { rows { mechanismOfAction actionType targets { id approvedSymbol } references { source ids urls } } } }
  clinicalReports { id source clinicalStage trialOverallStatus url origin }
"""

DISEASE_Q = """
query($id: String!, $nt: Int!) {
  meta { dataVersion { year month } }
  disease(efoId: $id) {
    id name description dbXRefs
    synonyms { relation terms }
    drugAndClinicalCandidates { count rows { %s } }
    associatedTargets(page: {index: 0, size: $nt}) { count rows { score target { id approvedSymbol approvedName } } }
  }
}""" % _DRUG_FIELDS

TARGET_Q = """
query($id: String!, $nd: Int!) {
  meta { dataVersion { year month } }
  target(ensemblId: $id) {
    id approvedSymbol approvedName biotype
    tractability { label modality value }
    drugAndClinicalCandidates { count rows { %s diseases { diseaseFromSource disease { id name } } } }
    associatedDiseases(page: {index: 0, size: $nd}) { count rows { score disease { id name } } }
  }
}""" % _DRUG_FIELDS

SEARCH_Q = """
query($q: String!, $entities: [String!], $size: Int!) {
  search(queryString: $q, entityNames: $entities, page: {index: 0, size: $size}) {
    total hits { id name entity description }
  }
}"""


def gql(query: str, variables: Dict[str, Any], record: str):
    """POST a GraphQL query; GraphQL errors (HTTP 200 with `errors`) raise SourceError."""
    resp = post_json(API, {"query": query, "variables": variables}, source="Open Targets", cache_ttl=7 * 86400, timeout=90)
    data = resp.json()
    if data.get("errors"):
        msg = "; ".join(str(e.get("message")) for e in data["errors"])[:300]
        raise SourceError("Open Targets", API, resp.status, f"GraphQL error for {record}: {msg}")
    return resp, data.get("data") or {}


def normalise_disease_id(text: str) -> Optional[str]:
    t = text.strip()
    if not DISEASE_ID_RE.match(t):
        return None
    prefix, num = re.split(r"[_:]", t, maxsplit=1)
    return f"{_PREFIX_CASE[prefix.lower()]}_{num}"


def search(text: str, entities: Sequence[str] = ("disease", "target"), limit: int = 10) -> Outcome:
    resp, data = gql(SEARCH_Q, {"q": text, "entities": list(entities), "size": max(1, min(limit, 50))}, text)
    s = data.get("search") or {}
    hits = [{"id": h.get("id"), "name": h.get("name"), "entity": h.get("entity"),
             "description": (h.get("description") or "")[:200] or None} for h in s.get("hits") or []]
    return Outcome({"query": text, "total": s.get("total"), "hits": hits},
                   sources=[source_record("Open Targets search", text, resp)])


def _stage_rank(stage: Optional[str]) -> int:
    return STAGE_ORDER.index(stage) if stage in STAGE_ORDER else len(STAGE_ORDER)


def _report_rank(r: Dict[str, Any]) -> tuple:
    return (r.get("origin") != "REGULATORY_AGENCY", _stage_rank(r.get("clinicalStage")), r.get("source") or "")


def _evidence(r: Dict[str, Any]) -> Dict[str, Any]:
    ev = {"source": r.get("source"), "stage": r.get("clinicalStage"), "status": r.get("trialOverallStatus"), "url": r.get("url")}
    if not ev["url"]:
        ev["id"] = r.get("id")
    return ev


def parse_drug_row(row: Dict[str, Any], with_diseases: bool = False) -> Dict[str, Any]:
    drug = row.get("drug") or {}
    moas = []
    for m in ((drug.get("mechanismsOfAction") or {}).get("rows") or [])[:3]:
        syms = [t.get("approvedSymbol") for t in (m.get("targets") or []) if t.get("approvedSymbol")]
        refs = []
        for ref in m.get("references") or []:
            refs.extend(ref.get("urls") or [])
        moas.append({"mechanism": m.get("mechanismOfAction"), "action": m.get("actionType"),
                     "targets": syms[:5] + ([f"+{len(syms) - 5} more"] if len(syms) > 5 else []),
                     "references": list(dict.fromkeys(refs))[:2]})
    reports = row.get("clinicalReports") or []
    ordered = sorted(reports, key=_report_rank)
    out: Dict[str, Any] = {
        "drug": drug.get("name"),
        "chembl_id": drug.get("id"),
        "type": drug.get("drugType"),
        "stage": row.get("maxClinicalStage"),
        "drug_max_stage": drug.get("maximumClinicalStage"),
        "mechanisms": moas,
        "reports": len(reports),
        "trials": sum(1 for r in reports if r.get("source") == "ClinicalTrials.gov"),
        "report_sources": sorted({r.get("source") for r in reports if r.get("source")}),
        "evidence": [_evidence(r) for r in ordered[:4]],
        "url": f"{WEB}/drug/{drug.get('id')}" if drug.get("id") else None,
    }
    if with_diseases:
        names = []
        for d in row.get("diseases") or []:
            n = (d.get("disease") or {}).get("name") or d.get("diseaseFromSource")
            if n and n not in names:
                names.append(n)
        out["indications"] = names[:6] + ([f"+{len(names) - 6} more"] if len(names) > 6 else [])
    return out


def _drugs(block: Dict[str, Any], limit: int, with_diseases: bool) -> Dict[str, Any]:
    rows = [r for r in (block or {}).get("rows") or [] if r.get("drug")]
    parsed = [parse_drug_row(r, with_diseases) for r in rows]
    parsed.sort(key=lambda d: (_stage_rank(d["stage"]), d["drug"] or ""))
    return {"count": (block or {}).get("count", len(rows)), "shown": min(limit, len(parsed)),
            "order": "by clinical stage for this disease/target, then name", "rows": parsed[:limit]}


def _version(data: Dict[str, Any]) -> Optional[str]:
    v = ((data.get("meta") or {}).get("dataVersion")) or {}
    return f"{v.get('year')}.{v.get('month')}" if v else None


def disease(efo_id: str, drug_limit: int = 40, n_targets: int = 10) -> Outcome:
    """A disease's drugs and clinical candidates (with their sources) and its top associated targets."""
    did = normalise_disease_id(efo_id) or efo_id
    resp, data = gql(DISEASE_Q, {"id": did, "nt": n_targets}, did)
    d = data.get("disease")
    if not d:
        return Outcome(None, sources=[source_record("Open Targets", did, resp)],
                       warnings=[f"Open Targets has no disease {did}"])
    exact = []
    for s in d.get("synonyms") or []:
        if s.get("relation") == "hasExactSynonym":
            exact.extend(s.get("terms") or [])
    result = {
        "id": d.get("id"), "name": d.get("name"),
        "description": (d.get("description") or "")[:400] or None,
        "xrefs": [x for x in d.get("dbXRefs") or [] if x.split(":")[0] in ("OMIM", "Orphanet", "MESH", "GARD", "UMLS",
                                                                           "ICD10CM", "NCIT", "DOID", "MEDGEN", "NORD")][:20],
        "synonyms": exact[:12],
        "drugs": _drugs(d.get("drugAndClinicalCandidates"), drug_limit, False),
        "top_targets": [{"symbol": (r.get("target") or {}).get("approvedSymbol"), "ensembl": (r.get("target") or {}).get("id"),
                         "score": round(r.get("score") or 0, 3)}
                        for r in ((d.get("associatedTargets") or {}).get("rows") or [])],
        "targets_associated": (d.get("associatedTargets") or {}).get("count"),
        "url": f"{WEB}/disease/{d.get('id')}",
        "data_version": _version(data),
    }
    return Outcome(result, sources=[source_record("Open Targets", did, resp, url=result["url"],
                                                  note=f"GraphQL API v4, data {result['data_version']}")])


def target(ensembl_id: str, drug_limit: int = 40, n_diseases: int = 5) -> Outcome:
    """A target's tractability by modality, drugs and clinical candidates acting on it, and top diseases."""
    resp, data = gql(TARGET_Q, {"id": ensembl_id, "nd": n_diseases}, ensembl_id)
    t = data.get("target")
    if not t:
        return Outcome(None, sources=[source_record("Open Targets", ensembl_id, resp)],
                       warnings=[f"Open Targets has no target {ensembl_id}"])
    tract: Dict[str, List[str]] = {v: [] for v in MODALITIES.values()}
    for row in t.get("tractability") or []:
        if row.get("value"):
            tract.setdefault(MODALITIES.get(row.get("modality"), row.get("modality") or "other"), []).append(row.get("label"))
    result = {
        "id": t.get("id"), "symbol": t.get("approvedSymbol"), "name": t.get("approvedName"), "biotype": t.get("biotype"),
        "tractability": tract,
        "drugs": _drugs(t.get("drugAndClinicalCandidates"), drug_limit, True),
        "top_diseases": [{"name": (r.get("disease") or {}).get("name"), "id": (r.get("disease") or {}).get("id"),
                          "score": round(r.get("score") or 0, 3)}
                         for r in ((t.get("associatedDiseases") or {}).get("rows") or [])],
        "url": f"{WEB}/target/{t.get('id')}",
        "data_version": _version(data),
    }
    return Outcome(result, sources=[source_record("Open Targets", ensembl_id, resp, url=result["url"],
                                                  note=f"GraphQL API v4, data {result['data_version']}")])
