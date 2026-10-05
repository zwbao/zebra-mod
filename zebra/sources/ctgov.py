"""ClinicalTrials.gov API v2: trials by condition, keyword, country and status.

API: https://clinicaltrials.gov/api/v2/studies (no key). The country filter is
`filter.advanced=AREA[LocationCountry]"<country>"`, which matches the
location's country field exactly; `query.locn` is a text search over
facility names too and lets Taiwan's "China Medical University Hospital"
through a China filter. ClinicalTrials.gov lists Hong Kong and Taiwan as their
own countries. Results keep the service's order (relevance).
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from zebra.core import Outcome
from zebra.http import get_json, source_record

BASE = "https://clinicaltrials.gov/api/v2/studies"
STATUSES = ("RECRUITING", "NOT_YET_RECRUITING", "ACTIVE_NOT_RECRUITING", "COMPLETED", "ENROLLING_BY_INVITATION",
            "SUSPENDED", "TERMINATED", "WITHDRAWN", "AVAILABLE", "UNKNOWN", "ANY")
FIELDS = ("NCTId,BriefTitle,Acronym,Phase,StudyType,OverallStatus,Condition,InterventionName,InterventionType,"
          "EnrollmentCount,EnrollmentType,StartDate,PrimaryCompletionDate,LocationCountry,LocationFacility,"
          "LocationCity,LocationStatus,MinimumAge,MaximumAge,StdAge,LeadSponsorName,LastUpdatePostDate,HasResults")
# common spellings → the country names ClinicalTrials.gov uses
COUNTRY_ALIASES = {
    "usa": "United States", "us": "United States", "u.s.": "United States", "u.s.a.": "United States",
    "united states of america": "United States", "america": "United States",
    "uk": "United Kingdom", "great britain": "United Kingdom", "britain": "United Kingdom", "england": "United Kingdom",
    "south korea": "Korea, Republic of", "korea": "Korea, Republic of", "republic of korea": "Korea, Republic of",
    "prc": "China", "mainland china": "China", "people's republic of china": "China", "中国": "China",
    "russia": "Russian Federation", "iran": "Iran, Islamic Republic of", "vietnam": "Vietnam",
    "czech republic": "Czechia", "turkey": "Turkey", "türkiye": "Turkey",
}


def country_name(country: str) -> str:
    c = country.strip()
    return COUNTRY_ALIASES.get(c.lower(), c)


def _params(condition: str, term: Optional[str], country: Optional[str], status: str, page_size: int,
            fields: str = FIELDS) -> Dict[str, Any]:
    p: Dict[str, Any] = {"query.cond": condition, "pageSize": page_size, "countTotal": "true", "format": "json",
                         "fields": fields}
    if term:
        p["query.term"] = term
    if country:
        p["filter.advanced"] = f'AREA[LocationCountry]"{country}"'
    if status != "ANY":
        p["filter.overallStatus"] = status
    return p


def parse_study(st: Dict[str, Any], country: Optional[str] = None) -> Dict[str, Any]:
    ps = st.get("protocolSection") or {}
    ident = ps.get("identificationModule") or {}
    status = ps.get("statusModule") or {}
    design = ps.get("designModule") or {}
    arms = ps.get("armsInterventionsModule") or {}
    elig = ps.get("eligibilityModule") or {}
    locs = (ps.get("contactsLocationsModule") or {}).get("locations") or []
    nct = ident.get("nctId")
    countries = sorted({l.get("country") for l in locs if l.get("country")})
    enrol = design.get("enrollmentInfo") or {}
    row: Dict[str, Any] = {
        "nct_id": nct,
        "title": ident.get("briefTitle"),
        "acronym": ident.get("acronym"),
        "phases": design.get("phases") or [],
        "study_type": design.get("studyType"),
        "status": status.get("overallStatus"),
        "conditions": ((ps.get("conditionsModule") or {}).get("conditions") or [])[:6],
        "interventions": [{"name": i.get("name"), "type": i.get("type")} for i in (arms.get("interventions") or [])][:8],
        "enrollment": enrol.get("count"),
        "enrollment_type": enrol.get("type"),
        "start_date": (status.get("startDateStruct") or {}).get("date"),
        "primary_completion_date": (status.get("primaryCompletionDateStruct") or {}).get("date"),
        "last_update": (status.get("lastUpdatePostDateStruct") or {}).get("date"),
        "countries": countries,
        "n_sites": len(locs),
        "min_age": elig.get("minimumAge"),
        "max_age": elig.get("maximumAge"),
        "age_groups": elig.get("stdAges") or [],
        "sponsor": ((ps.get("sponsorCollaboratorsModule") or {}).get("leadSponsor") or {}).get("name"),
        "has_results": st.get("hasResults"),
        "url": f"https://clinicaltrials.gov/study/{nct}" if nct else None,
    }
    if country:
        here = [l for l in locs if (l.get("country") or "").lower() == country.lower()]
        row["sites_in_country"] = len(here)
        row["sites"] = [{"facility": l.get("facility"), "city": l.get("city"), "status": l.get("status")} for l in here[:5]]
    return row


def search(condition: str, term: Optional[str] = None, country: Optional[str] = None, status: str = "RECRUITING",
           limit: int = 20) -> Outcome:
    """Trials for a condition; `status` is an overall status or ANY; `country` as ClinicalTrials.gov names it."""
    condition = (condition or "").strip()
    if not condition:
        raise ValueError("empty condition")
    status = (status or "RECRUITING").strip().upper()
    if status not in STATUSES:
        raise ValueError(f"status must be one of {', '.join(STATUSES)}")
    cname = country_name(country) if country else None
    limit = max(1, min(int(limit), 100))
    resp = get_json(BASE, source="ClinicalTrials.gov", params=_params(condition, term, cname, status, limit), cache_ttl=86400)
    data = resp.json()
    studies = [parse_study(s, cname) for s in data.get("studies") or []][:limit]
    result: Dict[str, Any] = {"condition": condition, "term": term, "country": cname, "status": status,
                              "total": data.get("totalCount"), "returned": len(studies), "studies": studies}
    sources = [source_record("ClinicalTrials.gov", f"{condition}" + (f" / {term}" if term else ""), resp,
                             note=f"status={status}" + (f"; country={cname}" if cname else ""))]
    warnings: List[str] = []
    if not studies and status != "ANY":
        any_resp = get_json(BASE, source="ClinicalTrials.gov",
                            params=_params(condition, term, cname, "ANY", 1, fields="NCTId"), cache_ttl=86400)
        result["total_any_status"] = any_resp.json().get("totalCount")
        sources.append(source_record("ClinicalTrials.gov", f"{condition} (any status count)", any_resp))
    if cname and not studies and not result.get("total_any_status"):
        warnings.append(f"no trials with a site in '{cname}': the country must be spelled as ClinicalTrials.gov does "
                        "(e.g. United States, Korea, Republic of)")
    return Outcome(result, sources=sources, warnings=warnings)
