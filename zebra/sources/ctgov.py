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

import re
import unicodedata
from datetime import date, datetime, timezone
from typing import Any, Dict, List, Optional, Sequence, Tuple

from zebra.core import Outcome
from zebra.http import get_json
from zebra.sources import record as source_record
from zebra.sources import validated_json

BASE = "https://clinicaltrials.gov/api/v2/studies"
STATUSES = ("RECRUITING", "NOT_YET_RECRUITING", "ACTIVE_NOT_RECRUITING", "COMPLETED", "ENROLLING_BY_INVITATION",
            "SUSPENDED", "TERMINATED", "WITHDRAWN", "AVAILABLE", "UNKNOWN", "ANY")
FIELDS = ("NCTId,BriefTitle,Acronym,Phase,StudyType,OverallStatus,Condition,InterventionName,InterventionType,"
          "EnrollmentCount,EnrollmentType,StartDate,PrimaryCompletionDate,LocationCountry,LocationFacility,"
          "LocationCity,LocationStatus,MinimumAge,MaximumAge,StdAge,LeadSponsorName,LastUpdatePostDate,HasResults,"
          "EligibilityCriteria,HealthyVolunteers,Sex,CentralContactName,CentralContactRole,CentralContactPhone,"
          "CentralContactEMail,LocationContactName,LocationContactRole,LocationContactEMail,LocationContactPhone,"
          "StatusVerifiedDate,OfficialTitle,Keyword,ConditionMeshTerm")
ELIGIBILITY_BUDGET = 600  # characters of criteria kept per trial unless full_eligibility
STALE_DAYS = 730  # ~2 years without an update posted
MAX_CONTACTS = 5  # central contacts kept per trial
MAX_SITE_CONTACTS = 3  # contacts kept per site
# an overall status that says the trial is not taking patients, against a site that says it is
SITE_OPEN_STATUSES = ("RECRUITING", "NOT_YET_RECRUITING")
OVERALL_SHUT_STATUSES = ("UNKNOWN", "COMPLETED", "TERMINATED", "WITHDRAWN", "SUSPENDED")
# E-3: one canonical key per country, for every spelling ClinicalTrials.gov has used.
# The server's country filter accepts the old ISO-style names ("Korea, Republic of")
# and matches the trials, but the records now spell the location "South Korea",
# "Iran" and "Turkey (Türkiye)", so comparing the strings dropped every site.
# Checked live 2026-10-06 (DMD/Korea 25 trials, thalassemia/Iran 9, FMF/Turkey 54).
COUNTRY_SYNONYMS = (
    ("south korea", "korea, republic of", "republic of korea", "korea", "korea (south)", "korea, south", "한국"),
    ("north korea", "korea, democratic people's republic of"),
    ("iran", "iran, islamic republic of", "islamic republic of iran"),
    ("turkey", "türkiye", "turkiye", "turkey (türkiye)", "turkey (turkiye)"),
    ("russia", "russian federation"),
    ("czechia", "czech republic"),
    ("vietnam", "viet nam"),
    ("syria", "syrian arab republic"),
    ("moldova", "moldova, republic of", "republic of moldova"),
    ("tanzania", "tanzania, united republic of"),
    ("venezuela", "venezuela, bolivarian republic of"),
    ("bolivia", "bolivia, plurinational state of"),
    ("laos", "lao people's democratic republic"),
    ("north macedonia", "macedonia, the former yugoslav republic of", "macedonia"),
    ("democratic republic of the congo", "congo, the democratic republic of the", "dr congo", "drc",
     "congo-kinshasa"),
    ("cote d'ivoire", "côte d'ivoire", "ivory coast"),
    ("palestine", "palestinian territory, occupied", "palestinian territories", "state of palestine"),
    ("republic of the congo", "congo", "congo, republic of", "congo-brazzaville"),
    ("macao", "macau", "macao sar", "macau sar"),
    ("taiwan", "taiwan, province of china"),
    ("hong kong", "hong kong sar", "hong kong, china"),
    ("united states", "usa", "us", "u.s.", "u.s.a.", "united states of america", "america"),
    ("united kingdom", "uk", "great britain", "britain", "england"),
    ("china", "prc", "mainland china", "people's republic of china", "中国"),
    ("netherlands", "the netherlands", "holland"),
)
_COUNTRY_KEY = {name: group[0] for group in COUNTRY_SYNONYMS for name in group}


def country_key(text: Optional[str]) -> str:
    """One key per country however it is spelled: 'Korea, Republic of' == 'South Korea' == 'korea'."""
    t = " ".join(str(text or "").replace("’", "'").split()).lower()
    if t in _COUNTRY_KEY:
        return _COUNTRY_KEY[t]
    base = t.split(" (")[0].strip()  # "Turkey (Türkiye)" -> "turkey"
    return _COUNTRY_KEY.get(base, base)


# common spellings → the country names ClinicalTrials.gov's filter accepts
COUNTRY_ALIASES = {
    "usa": "United States", "us": "United States", "u.s.": "United States", "u.s.a.": "United States",
    "united states of america": "United States", "america": "United States",
    "uk": "United Kingdom", "great britain": "United Kingdom", "britain": "United Kingdom", "england": "United Kingdom",
    "south korea": "South Korea", "korea": "South Korea", "republic of korea": "South Korea",
    "korea, republic of": "South Korea",
    "prc": "China", "mainland china": "China", "people's republic of china": "China", "中国": "China",
    "russia": "Russian Federation", "iran": "Iran", "iran, islamic republic of": "Iran", "vietnam": "Vietnam",
    "czech republic": "Czechia", "turkey": "Turkey", "türkiye": "Turkey", "turkiye": "Turkey",
}


def country_name(country: str) -> str:
    c = country.strip()
    return COUNTRY_ALIASES.get(c.lower(), c)


def _plain(text: Optional[str]) -> Optional[str]:
    """Words only: quotes, backslashes and the search operators AND/OR/NOT make the v2 API answer 400."""
    if text is None:
        return None
    t = re.sub(r'["\\()\[\]{}]+', " ", text)
    t = re.sub(r"\b(AND|OR|NOT)\b", " ", t)
    return " ".join(t.split())


def _params(condition: str, term: Optional[str], country: Optional[str], status: str, page_size: int,
            fields: str = FIELDS) -> Dict[str, Any]:
    p: Dict[str, Any] = {"query.cond": _plain(condition), "pageSize": page_size, "countTotal": "true",
                         "format": "json", "fields": fields}
    if term:
        p["query.term"] = _plain(term)
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
    # the relevance check reads every condition, keyword and MeSH term (PEARL lists Pompe disease as its
    # 8th condition); `search()` removes this private block before the row is returned
    row["_relevance"] = {
        "conditions": (ps.get("conditionsModule") or {}).get("conditions") or [],
        "official_title": ident.get("officialTitle"),
        "keywords": (ps.get("conditionsModule") or {}).get("keywords") or [],
        "mesh": [m.get("term") for m in ((st.get("derivedSection") or {}).get("conditionBrowseModule") or {}).get("meshes")
                 or [] if m.get("term")]}
    if country:
        want = country_key(country)
        here = [l for l in locs if country_key(l.get("country")) == want]
        row["sites_in_country"] = len(here)
        row["sites"] = [{"facility": l.get("facility"), "city": l.get("city"), "status": l.get("status"),
                         "contacts": _contacts(l.get("contacts"), MAX_SITE_CONTACTS)} for l in here[:5]]
    return row


# words that do not make a condition name specific (CP1-8 relevance check)
_STOP = frozenset(("of", "the", "and", "with", "in", "or", "an", "to", "for", "type", "disease", "diseases",
                   "syndrome", "syndromes", "disorder", "disorders", "condition"))
# British and American spellings compare equal (haemophilia ~ hemophilia, tumour ~ tumor)
_SPELLING = ((re.compile(r"ae"), "e"), (re.compile(r"oe"), "e"), (re.compile(r"our\b"), "or"))


def _fold(text: Optional[str]) -> str:
    t = unicodedata.normalize("NFKD", str(text or ""))
    t = "".join(c for c in t if not unicodedata.combining(c)).lower().replace("’", "'")
    t = re.sub(r"'s\b", "", t)
    for pat, rep in _SPELLING:
        t = pat.sub(rep, t)
    return t


def _words(text: Optional[str], keep_short: bool = True) -> set:
    """Specific words, folded: plural endings dropped; one-character words kept (haemophilia A vs B, type 1)."""
    out = set()
    for w in re.split(r"[^0-9a-z]+", _fold(text)):
        if not w or w in _STOP or (len(w) < 2 and not keep_short):
            continue
        if w.endswith("ies") and len(w) > 4:
            w = w[:-3] + "y"
        elif w.endswith("s") and len(w) > 3 and not w.endswith("ss"):
            w = w[:-1]
        out.add(w)
    return out


def _has_acronym(condition: str) -> bool:
    """MPS II, SMA type 1, CDKL5: a capitalised short token the registry may spell out."""
    return any(len(w) >= 2 and sum(c.isupper() for c in w) >= 2 for w in re.split(r"[^0-9A-Za-z]+", condition or ""))


def _relevance_fields(study: Dict[str, Any]) -> List[Tuple[str, str]]:
    full = study.get("_relevance") or {}
    return [("condition", c) for c in full.get("conditions") or study.get("conditions") or []] + \
        [("title", study.get("title")), ("official title", full.get("official_title") or study.get("official_title")),
         ("acronym", study.get("acronym"))] + \
        [("keyword", k) for k in full.get("keywords") or study.get("keywords") or []] + \
        [("MeSH term", m) for m in full.get("mesh") or study.get("condition_mesh") or []]


def relevance(study: Dict[str, Any], condition: str) -> Optional[str]:
    """Where a trial names the condition asked for, or None when it does not (CP1-8).

    ClinicalTrials.gov expands a condition through MeSH: "spinal muscular atrophy"
    brings back androgen-receptor prostate and breast cancer trials, because their
    keyword "androgen receptor" maps to "Bulbo-Spinal Atrophy, X-Linked". A trial
    names the condition when every specific word of the query appears in one of its
    conditions, its title, official title, acronym, keywords or MeSH condition terms
    (all of them, not the shortened lists shown), after folding accents, British
    spellings, possessives and plurals; single letters and numbers count (haemophilia
    A is not haemophilia B). `search()` adds a second rule: a trial sharing a MeSH
    term with most of the trials that do name the condition is kept too (Hunter
    syndrome ~ MPS II), and nothing is set aside when most trials word it otherwise.
    """
    want = _words(condition)
    if not want:
        return "query has no specific words"
    # the initialism of a multi-word name (Duchenne muscular dystrophy -> DMD), 3+ letters
    words = [w for w in re.split(r"[^0-9A-Za-z]+", condition) if w]
    initials = "".join(w[0] for w in words).upper() if len(words) >= 3 else ""
    for label, text in _relevance_fields(study):
        if not text:
            continue
        have = _words(text)
        if all(any(_same_word(w, h) for h in have) for w in want):
            return f"{label}: {text}"
        if initials and re.search(rf"(?<![A-Za-z0-9]){initials}(?![A-Za-z])", str(text)):
            return f"{label}: {text} (initialism {initials})"
    return None


def _same_word(a: str, b: str) -> bool:
    """Equal, or one spelling slip apart on a long stem (dystrophy / dystrophin)."""
    if a == b:
        return True
    if min(len(a), len(b)) < 6:
        return False
    common = 0
    for x, y in zip(a, b):
        if x != y:
            break
        common += 1
    return common >= max(6, min(len(a), len(b)) - 2)


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
           limit: int = 20, full_eligibility: bool = False, today: Optional[date] = None,
           filter_relevance: bool = True) -> Outcome:
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
    asked_limit = int(limit)
    limit = max(1, min(asked_limit, 100))
    # CP1-8: ask for more than `limit` so trials set aside as unrelated still leave a full page
    page = min(100, limit * 2) if filter_relevance else limit
    resp = get_json(BASE, source="ClinicalTrials.gov", params=_params(condition, term, cname, status, page), cache_ttl=86400)
    data = validated_json(resp, "ClinicalTrials.gov")
    if today is None:
        today = datetime.now(timezone.utc).date()
    parsed = [parse_study(s, cname, full_eligibility=full_eligibility, today=today)
              for s in data.get("studies") or []]
    warnings: List[str] = []
    if asked_limit > 100:
        warnings.append(f"--limit {asked_limit} is above the 100 trials one page can hold: 100 were asked for")
    filtered: List[Dict[str, Any]] = []
    check = filter_relevance and not _CJK.search(condition) and not _has_acronym(condition)
    skipped_reason = None
    if check:
        named = {st["nct_id"]: relevance(st, condition) for st in parsed}
        passing = [st for st in parsed if named[st["nct_id"]]]
        if parsed and len(passing) * 2 < len(parsed):
            # most records word the condition otherwise (Lou Gehrig disease ~ ALS): a word test cannot
            # tell relevant from unrelated here, so nothing is set aside
            check, skipped_reason = False, (f"not applied: only {len(passing)} of {len(parsed)} trials name "
                                            f"'{condition}' in those words, so the registry words it otherwise")
        else:
            # MeSH terms carried by at least half of the trials that name the condition
            counts: Dict[str, int] = {}
            for st in passing:
                for m in set((st.get("_relevance") or {}).get("mesh") or []):
                    counts[m] = counts.get(m, 0) + 1
            core_mesh = {m for m, n in counts.items() if n * 2 >= len(passing) and len(passing) >= 2}
            kept = []
            for st in parsed:
                why = named[st["nct_id"]]
                if not why:
                    shared = sorted(core_mesh & set((st.get("_relevance") or {}).get("mesh") or []))
                    if shared:
                        why = f"MeSH term shared with most trials that name the condition: {shared[0]}"
                if why:
                    st["relevance"] = why
                    kept.append(st)
                else:
                    filtered.append({"nct_id": st["nct_id"], "title": st["title"], "conditions": st["conditions"],
                                     "url": st["url"]})
            parsed = kept
    for st in parsed:
        st.pop("_relevance", None)
    studies = parsed[:limit]
    result: Dict[str, Any] = {"condition": condition, "term": term, "country": cname, "status": status,
                              "total": data.get("totalCount"), "returned": len(studies), "studies": studies,
                              "eligibility_budget": None if full_eligibility else ELIGIBILITY_BUDGET,
                              "checked_on": today.isoformat()}
    if check:
        result["relevance_check"] = ("a trial is listed when every specific word of the condition appears in its "
                                     "conditions, title, keywords or MeSH condition terms, or when it shares a MeSH "
                                     "term with most of the trials that do")
        result["filtered"] = filtered
        if filtered:
            warnings.append(
                f"{len(filtered)} trial(s) that ClinicalTrials.gov returned for '{condition}' do not name it in their "
                "conditions, title, keywords or MeSH terms and were set aside as possibly unrelated "
                "(result.filtered; they may still be relevant — --keep-unrelated lists them in full): "
                + "; ".join(f"{f['nct_id']} ({', '.join(f['conditions'][:2]) or 'no condition listed'})"
                            for f in filtered[:6])
                + (f"; +{len(filtered) - 6} more" if len(filtered) > 6 else ""))
    elif filter_relevance:
        result["relevance_check"] = skipped_reason or (
            "not applied: the condition is or contains an acronym, or is not in English, so a trial that spells "
            "it out cannot be told from an unrelated one")
    sources = [source_record("ClinicalTrials.gov", f"{condition}" + (f" / {term}" if term else ""), resp,
                             note=f"status={status}" + (f"; country={cname}" if cname else ""))]
    stale = _stale_warning(studies)
    if stale:
        warnings.append(stale)
    if not studies and status != "ANY":
        any_resp = get_json(BASE, source="ClinicalTrials.gov",
                            params=_params(condition, term, cname, "ANY", 1, fields="NCTId"), cache_ttl=86400)
        result["total_any_status"] = any_resp.json().get("totalCount")
        sources.append(source_record("ClinicalTrials.gov", f"{condition} (any status count)", any_resp))
    if cname and not studies and not result.get("total_any_status") and not filtered:
        warnings.append(f"no trials with a site in '{cname}': the country must be spelled as ClinicalTrials.gov does "
                        "(e.g. United States, South Korea)")
    if cname:
        lost = [s["nct_id"] for s in studies if s.get("sites_in_country") == 0]
        if lost:
            seen = sorted({c for s in studies if s["nct_id"] in lost for c in s.get("countries") or []})
            warnings.append(f"{len(lost)} trial(s) matched the '{cname}' filter on ClinicalTrials.gov but none of their "
                            f"sites reads as that country here (sites list: {', '.join(seen[:8])}): "
                            + ", ".join(lost[:8]) + " — open the trial page for its sites")
    return Outcome(result, sources=sources, warnings=warnings)


_CJK = re.compile(r"[㐀-鿿]")
