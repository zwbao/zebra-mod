"""ClinicalTrials.gov API v2: trials by condition, keyword, country and status.

API: https://clinicaltrials.gov/api/v2/studies (no key). The country filter is
`filter.advanced=AREA[LocationCountry]"<country>"`, which matches the
location's country field exactly; `query.locn` is a text search over
facility names too and lets Taiwan's "China Medical University Hospital"
through a China filter. ClinicalTrials.gov lists Hong Kong and Taiwan as their
own countries. Results keep the service's order (relevance).

Eligibility, who may enrol and who to write to (F38). `FIELDS` asks for
`EligibilityCriteria`, `HealthyVolunteers`, `Sex`, the five central-contact
fields, the four location-contact fields and `StatusVerifiedDate`; every one of
those names was checked live against the v2 API on 2026-10-06 (an unknown field
name makes the service answer 400, so the check is not a formality). They land
in `eligibility`, `healthy_volunteers`, `sex`, `contacts` (`MAX_CONTACTS`
central contacts at most), each site's `contacts` (`MAX_SITE_CONTACTS` at most)
and `status_verified`. `StatusVerifiedDate` is served, as
`statusModule.statusVerifiedDate`, at year-month granularity ("2012-11").
`LocationContactRole` is asked for as well; the captured responses carry
CONTACT, PRINCIPAL_INVESTIGATOR and SUB_INVESTIGATOR, and `role` is None when
the registration leaves it out.

Criteria are long: across the 31 study records captured in
tests/fixtures/ctgov/f38_*.json the median criteria text is 1,638 characters and
the longest 4,543. Parsed with their criteria whole, the 20 trials of
f38_dmd_china.json serialise to 70,831 characters, past the 60,000 a tool result
is cut at; at the budget below the same 20 come to 40,398. So `eligibility` is
whitespace-collapsed and cut to `ELIGIBILITY_BUDGET` (600) characters, with
`eligibility_truncated` true and `eligibility_chars` carrying the full length.
`parse_study(..., full_eligibility=True)` keeps the whole text verbatim instead;
that is what `zebra trials --full-eligibility` passes, and it is the caller's job
to keep `--limit` small enough.

`status_flags` records what makes a record untrustworthy: a `lastUpdatePostDate`
more than `STALE_DAYS` (730) before today, an overall status of UNKNOWN, and a
site still advertised as recruiting under an overall status that says the trial
is not. "Today" is `datetime.now(timezone.utc)`, never a baked-in date.
"""

from __future__ import annotations

from datetime import date, datetime, timezone
from typing import Any, Dict, List, Optional, Sequence

from zebra.core import Outcome
from zebra.http import get_json
from zebra.sources import record as source_record

BASE = "https://clinicaltrials.gov/api/v2/studies"
STATUSES = ("RECRUITING", "NOT_YET_RECRUITING", "ACTIVE_NOT_RECRUITING", "COMPLETED", "ENROLLING_BY_INVITATION",
            "SUSPENDED", "TERMINATED", "WITHDRAWN", "AVAILABLE", "UNKNOWN", "ANY")
FIELDS = ("NCTId,BriefTitle,Acronym,Phase,StudyType,OverallStatus,Condition,InterventionName,InterventionType,"
          "EnrollmentCount,EnrollmentType,StartDate,PrimaryCompletionDate,LocationCountry,LocationFacility,"
          "LocationCity,LocationStatus,MinimumAge,MaximumAge,StdAge,LeadSponsorName,LastUpdatePostDate,HasResults,"
          "EligibilityCriteria,HealthyVolunteers,Sex,CentralContactName,CentralContactRole,CentralContactPhone,"
          "CentralContactEMail,LocationContactName,LocationContactRole,LocationContactEMail,LocationContactPhone,"
          "StatusVerifiedDate")
ELIGIBILITY_BUDGET = 600  # characters of criteria kept per trial unless full_eligibility
STALE_DAYS = 730  # ~2 years without an update posted
MAX_CONTACTS = 5  # central contacts kept per trial
MAX_SITE_CONTACTS = 3  # contacts kept per site
# an overall status that says the trial is not taking patients, against a site that says it is
SITE_OPEN_STATUSES = ("RECRUITING", "NOT_YET_RECRUITING")
OVERALL_SHUT_STATUSES = ("UNKNOWN", "COMPLETED", "TERMINATED", "WITHDRAWN", "SUSPENDED")
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
        p["filter.advanced"] = f'AREA[LocationCountry]"{country.replace(chr(34), "")}"'
    if status != "ANY":
        p["filter.overallStatus"] = status
    return p


def _as_date(text: Optional[str]) -> Optional[date]:
    """A ClinicalTrials.gov date string → date. The service serves YYYY-MM-DD and YYYY-MM; unparseable → None."""
    for fmt in ("%Y-%m-%d", "%Y-%m", "%Y"):
        try:
            return datetime.strptime(str(text), fmt).date()
        except (TypeError, ValueError):
            continue
    return None


def _contact(c: Dict[str, Any]) -> Dict[str, Any]:
    return {"name": c.get("name"), "role": c.get("role"), "phone": c.get("phone"), "email": c.get("email")}


def _contacts(raw: Any, limit: int) -> List[Dict[str, Any]]:
    rows = [_contact(c) for c in (raw or []) if isinstance(c, dict)]
    return [c for c in rows if any(c.values())][:limit]


def _eligibility(text: Optional[str], full: bool) -> Dict[str, Any]:
    """`eligibility` plus its two bookkeeping keys; see the module docstring for the budget."""
    raw = (text or "").strip()
    out: Dict[str, Any] = {"eligibility": None, "eligibility_truncated": False, "eligibility_chars": len(raw)}
    if not raw:
        return out
    if full:
        out["eligibility"] = raw
        return out
    flat = " ".join(raw.split())
    if len(flat) > ELIGIBILITY_BUDGET:
        out["eligibility"] = flat[:ELIGIBILITY_BUDGET].rstrip() + "…"
        out["eligibility_truncated"] = True
    else:
        out["eligibility"] = flat
    return out


def status_flags(overall: Optional[str], last_update: Optional[str], site_statuses: Sequence[Optional[str]],
                 today: Optional[date] = None) -> List[str]:
    """Machine-readable reasons to distrust a record. `today` defaults to the UTC date now."""
    if today is None:
        today = datetime.now(timezone.utc).date()
    overall = (overall or "").strip().upper()
    flags: List[str] = []
    posted = _as_date(last_update)
    if posted is not None and (today - posted).days > STALE_DAYS:
        flags.append(f"stale: last update {last_update}, more than 2 years ago")
    if overall == "UNKNOWN":
        flags.append("status UNKNOWN: the sponsor has not verified this record")
    if overall in OVERALL_SHUT_STATUSES and any((s or "").strip().upper() in SITE_OPEN_STATUSES
                                                for s in site_statuses):
        flags.append("site status conflicts with overall status")
    return flags


def parse_study(st: Dict[str, Any], country: Optional[str] = None, full_eligibility: bool = False,
                today: Optional[date] = None) -> Dict[str, Any]:
    ps = st.get("protocolSection") or {}
    ident = ps.get("identificationModule") or {}
    status = ps.get("statusModule") or {}
    design = ps.get("designModule") or {}
    arms = ps.get("armsInterventionsModule") or {}
    elig = ps.get("eligibilityModule") or {}
    cl = ps.get("contactsLocationsModule") or {}
    locs = cl.get("locations") or []
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
        "sex": elig.get("sex"),
        "healthy_volunteers": elig.get("healthyVolunteers"),
        "sponsor": ((ps.get("sponsorCollaboratorsModule") or {}).get("leadSponsor") or {}).get("name"),
        "has_results": st.get("hasResults"),
        "status_verified": status.get("statusVerifiedDate"),
        "contacts": _contacts(cl.get("centralContacts"), MAX_CONTACTS),
        "url": f"https://clinicaltrials.gov/study/{nct}" if nct else None,
    }
    row.update(_eligibility(elig.get("eligibilityCriteria"), full_eligibility))
    row["status_flags"] = status_flags(row["status"], row["last_update"],
                                       [l.get("status") for l in locs], today=today)
    if country:
        here = [l for l in locs if (l.get("country") or "").lower() == country.lower()]
        row["sites_in_country"] = len(here)
        row["sites"] = [{"facility": l.get("facility"), "city": l.get("city"), "status": l.get("status"),
                         "contacts": _contacts(l.get("contacts"), MAX_SITE_CONTACTS)} for l in here[:5]]
    return row


def _stale_warning(studies: List[Dict[str, Any]]) -> Optional[str]:
    """One warning naming the trials that carry a status flag, so the Outcome says so before the text does."""
    flagged = [s for s in studies if s.get("status_flags")]
    if not flagged:
        return None
    ids = [s.get("nct_id") or "?" for s in flagged]
    shown = ", ".join(ids[:10]) + (f" +{len(ids) - 10} more" if len(ids) > 10 else "")
    return (f"{len(flagged)} of {len(studies)} trials may be out of date (status UNKNOWN or last updated over "
            f"2 years ago): {shown}")


def search(condition: str, term: Optional[str] = None, country: Optional[str] = None, status: str = "RECRUITING",
           limit: int = 20, full_eligibility: bool = False, today: Optional[date] = None) -> Outcome:
    """Trials for a condition; `status` is an overall status or ANY; `country` as ClinicalTrials.gov names it.

    Eligibility criteria are cut to ELIGIBILITY_BUDGET characters unless `full_eligibility`.
    `today` is the day staleness is counted from; it defaults to the UTC date now and
    exists so a test against a frozen response can ask the question on a frozen day.
    """
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
    if today is None:
        today = datetime.now(timezone.utc).date()
    studies = [parse_study(s, cname, full_eligibility=full_eligibility, today=today)
               for s in data.get("studies") or []][:limit]
    result: Dict[str, Any] = {"condition": condition, "term": term, "country": cname, "status": status,
                              "total": data.get("totalCount"), "returned": len(studies), "studies": studies,
                              "eligibility_budget": None if full_eligibility else ELIGIBILITY_BUDGET,
                              "checked_on": today.isoformat()}
    sources = [source_record("ClinicalTrials.gov", f"{condition}" + (f" / {term}" if term else ""), resp,
                             note=f"status={status}" + (f"; country={cname}" if cname else ""))]
    warnings: List[str] = []
    stale = _stale_warning(studies)
    if stale:
        warnings.append(stale)
    if not studies and status != "ANY":
        any_resp = get_json(BASE, source="ClinicalTrials.gov",
                            params=_params(condition, term, cname, "ANY", 1, fields="NCTId"), cache_ttl=86400)
        result["total_any_status"] = any_resp.json().get("totalCount")
        sources.append(source_record("ClinicalTrials.gov", f"{condition} (any status count)", any_resp))
    if cname and not studies and not result.get("total_any_status"):
        warnings.append(f"no trials with a site in '{cname}': the country must be spelled as ClinicalTrials.gov does "
                        "(e.g. United States, Korea, Republic of)")
    return Outcome(result, sources=sources, warnings=warnings)
