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
from zebra.http import SourceError, post_json
from zebra.sources import record as source_record
from zebra.sources import validated_json

API = "https://api.platform.opentargets.org/api/v4/graphql"
WEB = "https://platform.opentargets.org"
ENSG_RE = re.compile(r"^ENSG\d{11}$", re.I)
DISEASE_ID_RE = re.compile(r"^(MONDO|EFO|Orphanet|HP|DOID|OTAR|GO|NCIT)[_:]\d+$", re.I)
_PREFIX_CASE = {"mondo": "MONDO", "efo": "EFO", "orphanet": "Orphanet", "hp": "HP", "doid": "DOID", "otar": "OTAR",
                "go": "GO", "ncit": "NCIT"}
# highest first; anything else ranks last
STAGE_ORDER = ["APPROVAL", "PREAPPROVAL", "PHASE_4", "PHASE_3", "PHASE_2_3", "PHASE_2", "PHASE_1_2", "PHASE_1",
               "EARLY_PHASE_1", "PHASE_0", "PRECLINICAL", "UNKNOWN"]
# P1g: "APPROVAL" in `maxClinicalStage` is the highest stage ever reached, with
# no jurisdiction and no current status, so a withdrawn authorisation showed as
# an approval. Ataluren for Duchenne muscular dystrophy is the case: its only
# regulatory report is `EMA Human Drugs / WITHDRAWAL` (verified live on
# 2026-10-06 via MONDO_0010679), while the row's `maxClinicalStage` is
# APPROVAL. The regulatory reports are therefore read per agency.
# Report `origin` and `source` values observed live across four rare-disease
# queries (MONDO_0010679, MONDO_0100135, MONDO_0019079, MONDO_0018150):
# REGULATORY_AGENCY/{FDA, EMA Human Drugs, PMDA}, DRUG_LABEL/DailyMed,
# CLINICAL_TRIAL/ClinicalTrials.gov, CURATED_RESOURCE/{TTD, USAN, INN}.
JURISDICTIONS = {
    "FDA": "FDA (United States)",
    "EMA Human Drugs": "EMA (European Union)",
    "PMDA": "PMDA (Japan)",
    "DailyMed": "United States (DailyMed label archive)",
}
# clinicalStage values seen on REGULATORY_AGENCY reports, and how to read them
REGULATORY_STAGE_READING = {
    "APPROVAL": "marketing authorisation on record",
    "PREAPPROVAL": "under review, not authorised",
    "WITHDRAWAL": "authorisation withdrawn — NOT currently approved",
    "SUSPENSION": "authorisation suspended — NOT currently usable",
    "REFUSAL": "authorisation refused",
    "UNKNOWN": "the agency holds a record, but Open Targets resolved no stage from it",
}
NOT_APPROVED_STAGES = frozenset(("WITHDRAWAL", "WITHDRAWN", "SUSPENSION", "SUSPENDED", "REFUSAL", "REFUSED",
                                 "REVOCATION"))
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


def _post(payload: Dict[str, Any], *, refresh: bool = False):
    # E-6: `refresh=True` skips the cache read and still writes the good answer back;
    # `cache_ttl=0` (the old recovery) disabled the write, so the bad entry stayed and
    # every later call went to the network for the rest of its 7 days
    return post_json(API, payload, source="Open Targets", cache_ttl=7 * 86400, refresh=refresh, timeout=90)


def gql(query: str, variables: Dict[str, Any], record: str):
    """POST a GraphQL query; GraphQL errors (HTTP 200 with `errors`) raise SourceError.

    A 200 carrying `errors` is an error answer that `zebra.http` has already
    written to its 7-day cache (F10), so it was raised again on every later call
    with no network request at all — a transient upstream timeout became a
    week-long outage. `validated_json` re-requests once when the bad body came
    from the cache, and only then gives up.
    """
    payload = {"query": query, "variables": variables}
    resp = _post(payload)
    try:
        data = validated_json(resp, "Open Targets", refetch=lambda: _post(payload, refresh=True))
    except SourceError as err:
        raise SourceError("Open Targets", API, err.status, f"GraphQL error for {record}: {err.message}") from None
    if data.get("errors"):  # defensive: validated_json already rejects a top-level `errors`
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
    ev = {"source": r.get("source"), "stage": r.get("clinicalStage"), "status": r.get("trialOverallStatus"),
          "origin": r.get("origin"), "url": r.get("url")}
    if not ev["url"]:
        ev["id"] = r.get("id")
    return ev


def regulatory_status(reports: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Per-jurisdiction reading of the regulatory reports on one drug row (P1g).

    Returns `{"by_jurisdiction": [...], "approved_in": [...],
    "withdrawn_or_suspended_in": [...], "unresolved_in": [...], "label_only_in":
    [...], "headline": str}`. "Approved" is only ever said of a jurisdiction
    whose own agency report says so; a withdrawal is reported as a withdrawal,
    not folded into a highest-stage-ever.
    """
    per: Dict[str, Dict[str, Any]] = {}
    labels: Dict[str, Dict[str, Any]] = {}
    for r in reports or []:
        src = r.get("source") or "?"
        who = JURISDICTIONS.get(src, src)
        stage = (r.get("clinicalStage") or "UNKNOWN").upper()
        row = {"jurisdiction": who, "source": src, "stage": stage,
               "reading": REGULATORY_STAGE_READING.get(stage, f"stage {stage} as the agency record reports it"),
               "url": r.get("url")}
        if r.get("origin") == "REGULATORY_AGENCY":
            cur = per.get(who)
            # a withdrawal outranks an approval for the same agency: it is the later fact
            if cur is None or (stage in NOT_APPROVED_STAGES and cur["stage"] not in NOT_APPROVED_STAGES):
                per[who] = row
        elif r.get("origin") == "DRUG_LABEL":
            labels.setdefault(who, row)
    approved = sorted(w for w, r in per.items() if r["stage"] == "APPROVAL")
    gone = sorted(w for w, r in per.items() if r["stage"] in NOT_APPROVED_STAGES)
    review = sorted(w for w, r in per.items() if r["stage"] == "PREAPPROVAL")
    unresolved = sorted(w for w, r in per.items() if r["stage"] not in NOT_APPROVED_STAGES
                        and r["stage"] not in ("APPROVAL", "PREAPPROVAL"))
    parts = []
    if approved:
        parts.append("authorisation on record: " + ", ".join(approved))
    if gone:
        parts.append("WITHDRAWN/SUSPENDED: " + ", ".join(f"{w} ({per[w]['stage']})" for w in gone))
    if review:
        parts.append("under review, not authorised: " + ", ".join(review))
    if unresolved:
        parts.append("agency record with no stage: " + ", ".join(unresolved))
    if not per and labels:
        parts.append("no agency report; a drug label exists for " + ", ".join(sorted(labels)))
    if not parts:
        parts.append("no regulatory report in Open Targets for this drug and disease")
    return {
        "by_jurisdiction": [per[w] for w in sorted(per)],
        "approved_in": approved,
        "withdrawn_or_suspended_in": gone,
        "unresolved_in": unresolved,
        "under_review_in": review,
        "label_only_in": sorted(w for w in labels if w not in per),
        "headline": "; ".join(parts),
    }


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
        "stage_meaning": "highest stage this drug reached FOR THIS DISEASE in any jurisdiction, ever — "
                         "not its current approval status",
        "drug_max_stage": drug.get("maximumClinicalStage"),
        "drug_max_stage_meaning": "highest stage this drug reached for ANY indication, ever",
        "regulatory": regulatory_status(reports),
        "mechanisms": moas,
        "reports": len(reports),
        "trials": sum(1 for r in reports if r.get("source") == "ClinicalTrials.gov"),
        "report_sources": sorted({r.get("source") for r in reports if r.get("source")}),
        "evidence": [_evidence(r) for r in ordered[:3]],
        "url": f"{WEB}/drug/{drug.get('id')}" if drug.get("id") else None,
    }
    reg = out["regulatory"]
    if out["stage"] == "APPROVAL" and not reg["approved_in"]:
        # E-4: "do not read this as approved" is said only when an agency record says
        # withdrawn / suspended / refused. No record at all (cannabidiol and fenfluramine
        # for Dravet syndrome, both FDA- and EMA-approved) means "unknown here".
        if reg["withdrawn_or_suspended_in"]:
            out["stage_warning"] = (
                "Open Targets gives this row maxClinicalStage APPROVAL, but the agency record says: "
                f"{reg['headline']} — do not read this as an approved therapy")
        elif reg["by_jurisdiction"]:
            out["stage_warning"] = (
                "Open Targets gives this row maxClinicalStage APPROVAL; the agency records it holds give no stage "
                f"({reg['headline']}) — approval status unknown here: check Drugs@FDA, the EMA register and NMPA "
                "before calling it approved or not approved")
        else:
            out["stage_warning"] = (
                "Open Targets gives this row maxClinicalStage APPROVAL but holds no agency record for this drug and "
                "disease — approval status unknown here: check Drugs@FDA, the EMA register and NMPA before calling "
                "it approved or not approved")
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
    # the two "meaning" strings are the same on every row: said once for the block, so the
    # 60,000-character budget goes on rows rather than on 40 copies of the same sentence
    meanings = {}
    for p in parsed:
        for k in ("stage_meaning", "drug_max_stage_meaning"):
            if k in p:
                meanings[k] = p.pop(k)
    parsed.sort(key=lambda d: (_stage_rank(d["stage"]), d["drug"] or ""))
    withdrawn = [d["drug"] for d in parsed[:limit] if d["regulatory"]["withdrawn_or_suspended_in"]]
    return {"count": (block or {}).get("count", len(rows)), "shown": min(limit, len(parsed)),
            "order": "by clinical stage for this disease/target, then name",
            "what_this_is": "what Open Targets holds for this disease or target — drugs and clinical candidates "
                            "it has linked, at the highest stage each reached. It is NOT a list of approved "
                            "therapies and it is not complete: read `regulatory` per row for approval status",
            "withdrawn_or_suspended": withdrawn,
            **meanings,
            "rows": parsed[:limit]}


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
