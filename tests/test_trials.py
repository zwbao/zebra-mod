"""ClinicalTrials.gov v2 source and `zebra trials` (offline on captured responses, plus live)."""

from __future__ import annotations

import datetime as dt
import json
import os

import pytest

from zebra import cli
from zebra.http import Response
from zebra.sources import ctgov

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "ctgov")


def resp(name, url="https://clinicaltrials.gov/api/v2/studies"):
    if isinstance(name, str):
        with open(os.path.join(FIX, name), encoding="utf-8") as fh:
            text = fh.read()
    else:
        text = json.dumps(name)
    return Response(url, 200, text, "2026-10-05T00:00:00+00:00", False)


def test_params_use_exact_country_filter_and_status():
    p = ctgov._params("Duchenne muscular dystrophy", "exon 51", "Korea, Republic of", "RECRUITING", 20)
    assert p["filter.advanced"] == 'AREA[LocationCountry]"Korea, Republic of"'
    assert "query.locn" not in p  # text search over facility names leaks Taiwan into "China"
    assert p["filter.overallStatus"] == "RECRUITING" and p["query.term"] == "exon 51" and p["query.cond"].startswith("Duchenne")
    p = ctgov._params("Dravet syndrome", None, None, "ANY", 5)
    assert "filter.overallStatus" not in p and "filter.advanced" not in p and "query.term" not in p


def test_country_aliases():
    assert ctgov.country_name("USA") == "United States"
    assert ctgov.country_name("south korea") == "South Korea"  # E-3: the spelling records now carry
    assert ctgov.country_name("中国") == "China"
    assert ctgov.country_name("Taiwan") == "Taiwan"


def test_parse_study_fields():
    data = json.loads(open(os.path.join(FIX, "dravet_recruiting.json"), encoding="utf-8").read())
    rows = [ctgov.parse_study(s) for s in data["studies"]]
    fen = next(r for r in rows if r["nct_id"] == "NCT06598449")
    assert fen["phases"] == ["PHASE4"] and fen["status"] == "RECRUITING"
    assert fen["interventions"] == [{"name": "fenfluramine", "type": "DRUG"}]
    assert fen["min_age"] == "12 Months" and fen["max_age"] == "24 Months"
    assert fen["countries"] == ["United States"] and fen["n_sites"] == 1
    assert fen["url"] == "https://clinicaltrials.gov/study/NCT06598449" and fen["sponsor"]
    assert fen["start_date"] and fen["primary_completion_date"] and "sites" not in fen


def test_parse_study_country_sites():
    data = json.loads(open(os.path.join(FIX, "dmd_china.json"), encoding="utf-8").read())
    rows = [ctgov.parse_study(s, "China") for s in data["studies"]]
    assert all("China" in r["countries"] and r["sites_in_country"] >= 1 for r in rows)
    multi = next(r for r in rows if len(r["countries"]) > 1)
    assert multi["sites_in_country"] < multi["n_sites"]
    assert all(s["facility"] for s in multi["sites"])


def test_search_reports_any_status_count_when_empty(monkeypatch):
    calls = []

    def fake(url, source, params=None, **kw):
        calls.append(dict(params))
        if "filter.overallStatus" in params:
            return resp({"totalCount": 0, "studies": []})
        return resp({"totalCount": 7, "studies": [{"protocolSection": {"identificationModule": {"nctId": "NCT0"}}}]})

    monkeypatch.setattr(ctgov, "get_json", fake)
    out = ctgov.search("some rare disease", status="RECRUITING")
    assert out.result["studies"] == [] and out.result["total_any_status"] == 7
    assert len(calls) == 2 and calls[1]["pageSize"] == 1 and len(out.sources) == 2
    with pytest.raises(ValueError):
        ctgov.search("x", status="OPEN")


def test_search_country_zero_warns(monkeypatch):
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp({"totalCount": 0, "studies": []}))
    out = ctgov.search("Dravet syndrome", country="Narnia", status="ANY")
    assert any("spelled as ClinicalTrials.gov" in w for w in out.warnings)


def test_cli_trials(monkeypatch, capsys):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    seen = {}

    def fake(url, source, params=None, **kw):
        seen.update(params)
        return resp("dmd_china.json")

    monkeypatch.setattr(ctgov, "get_json", fake)
    assert cli.main(["trials", "Duchenne", "muscular", "dystrophy", "--country", "China", "--status", "any", "--limit", "3", "--json"]) == 0
    env = json.loads(capsys.readouterr().out)
    assert env["ok"] and env["result"]["country"] == "China" and env["result"]["status"] == "ANY"
    assert env["query"] == {"condition": "Duchenne muscular dystrophy", "term": None, "country": "China",
                            "status": "ANY", "limit": 3, "full_eligibility": False, "keep_unrelated": False}
    assert seen["filter.advanced"] == 'AREA[LocationCountry]"China"' and "filter.overallStatus" not in seen
    assert env["sources"][0]["db"] == "ClinicalTrials.gov" and env["sources"][0]["url"].startswith("https://clinicaltrials.gov/")


def test_cli_trials_bad_status(monkeypatch, capsys):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    assert cli.main(["trials", "Dravet syndrome", "--status", "OPEN", "--json"]) == 2
    assert json.loads(capsys.readouterr().out)["error"]["type"] == "UsageError"


# --------------------------------------------------------------- F38: eligibility, contacts, staleness
# Fixtures f38_*.json were fetched live on 2026-10-06 with the field list `_params` sends.

NEW_F38_FIELDS = ("EligibilityCriteria", "HealthyVolunteers", "Sex", "CentralContactName", "CentralContactRole",
                  "CentralContactPhone", "CentralContactEMail", "LocationContactName", "LocationContactRole",
                  "LocationContactEMail", "LocationContactPhone", "StatusVerifiedDate")
TODAY = dt.date(2026, 10, 6)  # the day the fixtures were captured


def load(name):
    with open(os.path.join(FIX, name), encoding="utf-8") as fh:
        return json.load(fh)


def by_id(name, nct, country=None, **kw):
    data = load(name)
    st = next(s for s in data["studies"] if s["protocolSection"]["identificationModule"]["nctId"] == nct)
    return ctgov.parse_study(st, country, today=TODAY, **kw)


def test_F38_statuses_include_available():
    """--status AVAILABLE is accepted by the mod's enum (expanded-access records)."""
    assert "AVAILABLE" in ctgov.STATUSES


def test_F38_fields_are_requested():
    p = ctgov._params("Dravet syndrome", None, None, "ANY", 5)
    for f in NEW_F38_FIELDS:
        assert f in ctgov.FIELDS, f
        assert f in p["fields"], f


def test_F38_eligibility_truncated():
    long = by_id("f38_dravet_recruiting.json", "NCT06585605")
    assert long["eligibility_chars"] == 855 and long["eligibility_truncated"] is True
    assert len(long["eligibility"]) == ctgov.ELIGIBILITY_BUDGET + 1 and long["eligibility"].endswith("…")
    assert "\n" not in long["eligibility"]  # whitespace collapsed so the budget is spent on words
    assert long["eligibility"].startswith("Inclusion Criteria")

    short = by_id("f38_dravet_recruiting.json", "NCT05126914")
    assert short["eligibility_chars"] == 264 and short["eligibility_truncated"] is False
    assert short["eligibility"] and "…" not in short["eligibility"]

    full = by_id("f38_dravet_recruiting.json", "NCT06585605", full_eligibility=True)
    assert full["eligibility_truncated"] is False and len(full["eligibility"]) == 855
    assert "\n" in full["eligibility"]  # the API's own line breaks survive


def test_F38_eligibility_absent_keeps_keys():
    """The pre-F38 fixture has no eligibilityModule: the keys are still there, with no invented text."""
    row = by_id("dmd_china.json", "NCT01610440")
    assert row["eligibility"] is None and row["eligibility_chars"] == 0 and row["eligibility_truncated"] is False
    assert row["sex"] is None and row["healthy_volunteers"] is None and row["contacts"] == []
    bare = ctgov.parse_study({"protocolSection": {"identificationModule": {"nctId": "NCT0"}}})
    assert bare["eligibility"] is None and bare["status_flags"] == [] and bare["contacts"] == []


def test_F38_healthy_volunteers_and_sex():
    row = by_id("f38_dravet_recruiting.json", "NCT06585605")
    assert row["sex"] == "ALL" and row["healthy_volunteers"] is False
    assert row["status_verified"] == "2026-03"


def test_F38_stale_status_flagged():
    row = by_id("f38_nct01610440.json", "NCT01610440")
    assert row["last_update"] == "2012-11-30"
    assert "stale: last update 2012-11-30, more than 2 years ago" in row["status_flags"]
    assert "status UNKNOWN: the sponsor has not verified this record" in row["status_flags"]
    fresh = by_id("f38_dravet_recruiting.json", "NCT06585605")
    assert fresh["status_flags"] == []  # updated 2026-03-18, RECRUITING, no conflict


def test_F38_stale_uses_today_not_a_fixed_date():
    """730 days is counted from the day the question is asked, so the boundary moves."""
    assert ctgov.status_flags("RECRUITING", "2024-10-06", [], today=dt.date(2026, 10, 5)) == []
    assert ctgov.status_flags("RECRUITING", "2024-10-06", [], today=dt.date(2026, 10, 7))[0].startswith("stale:")
    # ClinicalTrials.gov also serves year-month dates; neither form may raise
    assert ctgov.status_flags("RECRUITING", "2012-11", [], today=dt.date(2026, 10, 6))[0].startswith("stale:")
    assert ctgov.status_flags("RECRUITING", None, [], today=dt.date(2026, 10, 6)) == []
    assert ctgov.status_flags("RECRUITING", "not a date", [], today=dt.date(2026, 10, 6)) == []
    # with no `today` the threshold is the UTC date now: a 1900 update is stale whenever this runs
    assert ctgov.status_flags("RECRUITING", "1900-01-01", [])[0].startswith("stale:")


def test_F38_site_status_conflict_flagged():
    row = by_id("f38_nct01610440.json", "NCT01610440")
    assert row["status"] == "UNKNOWN" and row["status_flags"][-1] == "site status conflicts with overall status"
    for shut in ctgov.OVERALL_SHUT_STATUSES:
        assert "site status conflicts with overall status" in ctgov.status_flags(shut, "2026-10-01", ["RECRUITING"],
                                                                                today=TODAY)
    # a recruiting trial with closed sites is not a conflict, and a closed trial with closed sites is not either
    assert ctgov.status_flags("RECRUITING", "2026-10-01", ["ACTIVE_NOT_RECRUITING", "WITHDRAWN"], today=TODAY) == []
    assert ctgov.status_flags("COMPLETED", "2026-10-01", ["COMPLETED", None], today=TODAY) == []


def test_F38_contacts_parsed():
    row = by_id("f38_dravet_recruiting.json", "NCT06585605")
    assert row["contacts"] == [
        {"name": "Darius Ebrahimi-Fakhari, MD, PhD.", "role": "CONTACT", "phone": "617-355-0097",
         "email": "movementdisorders@childrens.harvard.edu"},
        {"name": "Vicente Quiroz, MD", "role": "CONTACT", "phone": None,
         "email": "movementdisorders@childrens.harvard.edu"},
    ]
    site = by_id("f38_nct01610440.json", "NCT01610440", country="China")
    assert site["sites_in_country"] == 1
    assert site["sites"][0]["facility"] == "The Second Affiliated Hospital of Kunming Medical College"
    assert site["sites"][0]["contacts"] == [
        {"name": "Liqing Yao", "role": "CONTACT", "phone": None, "email": "yaoliqing98731@yahoo.com.cn"},
        {"name": "Liqing Yao", "role": "PRINCIPAL_INVESTIGATOR", "phone": None, "email": None},
    ]
    # a site the API gives no contacts for gets an empty list, never a made-up one
    quiet = by_id("f38_dmd_china.json", "NCT06114056", country="China")
    assert all(s["contacts"] == [] for s in quiet["sites"])


def test_F38_search_warns_about_flagged_trials(monkeypatch):
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp("f38_dmd_china.json"))
    # today is pinned to the day the fixture was captured: a frozen response has to be asked about a frozen day,
    # or a trial updated 729 days before capture crosses the threshold and the count below drifts.
    out = ctgov.search("Duchenne muscular dystrophy", country="China", status="ANY", limit=20, today=TODAY)
    warn = next(w for w in out.warnings if "may be out of date" in w)
    assert warn.startswith("9 of 20 trials may be out of date (status UNKNOWN or last updated over 2 years ago): ")
    assert "NCT01610440" in warn
    flagged = [s["nct_id"] for s in out.result["studies"] if s["status_flags"]]
    assert len(flagged) == 9 and out.result["eligibility_budget"] == ctgov.ELIGIBILITY_BUDGET
    assert out.result["checked_on"] == "2026-10-06"
    assert ctgov.search("Dravet syndrome", status="ANY", limit=1).result["checked_on"] == \
        dt.datetime.now(dt.timezone.utc).date().isoformat()  # default is today, not a baked-in date


def test_F38_one_stale_trial_in_a_recruiting_page(monkeypatch):
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp("f38_dravet_recruiting.json"))
    out = ctgov.search("Dravet syndrome", status="RECRUITING", limit=10, today=TODAY)
    # NCT05651204 was last updated 2022-12-14, so it is stale from 2024-12-14 onward
    assert any("may be out of date" in w for w in out.warnings)
    assert [s["nct_id"] for s in out.result["studies"] if s["status_flags"]] == ["NCT05651204"]
    assert ctgov._stale_warning([{"nct_id": "NCT1", "status_flags": []}]) is None


def test_F38_stale_warning_caps_the_list():
    rows = [{"nct_id": f"NCT{i:08d}", "status_flags": ["status UNKNOWN: the sponsor has not verified this record"]}
            for i in range(14)]
    warn = ctgov._stale_warning(rows)
    assert warn.startswith("14 of 14 trials") and warn.endswith("+4 more") and warn.count("NCT") == 10


def test_F38_cli_renders_flags_eligibility_and_contacts(monkeypatch, capsys):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp("f38_nct01610440.json"))
    assert cli.main(["trials", "Duchenne muscular dystrophy", "--country", "China", "--status", "ANY"]) == 0
    cap = capsys.readouterr()
    out = cap.out
    assert "! stale: last update 2012-11-30, more than 2 years ago" in out
    assert "! status UNKNOWN: the sponsor has not verified this record" in out
    assert "! site status conflicts with overall status: a site still reads recruiting, but the overall status " \
           "(UNKNOWN) is the one to trust" in out
    assert "eligibility: Inclusion Criteria: * Aged 5-12 years" in out
    assert "[cut from 855 characters; --full-eligibility for all of it]" in out
    assert "contact: Liqing Yao (contact) · yaoliqing98731@yahoo.com.cn" in out
    assert "sex all · patients only" in out
    assert "Kunming Medical College (Kunming, RECRUITING) [Liqing Yao (contact)" in out
    assert "may be out of date" in cap.err  # the Outcome warning reaches the reader too


def test_F38_cli_full_eligibility(monkeypatch, capsys):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    seen = {}

    def fake(url, source, params=None, **kw):
        seen.update(params)
        return resp("f38_dmd_china.json")

    monkeypatch.setattr(ctgov, "get_json", fake)
    assert cli.main(["trials", "Duchenne muscular dystrophy", "--country", "China", "--status", "ANY",
                     "--limit", "20", "--full-eligibility", "--json"]) == 0
    env = json.loads(capsys.readouterr().out)
    r = env["result"]
    assert env["query"]["full_eligibility"] is True and r["eligibility_budget"] is None
    assert all(s["eligibility_truncated"] is False for s in r["studies"])
    full = next(s for s in r["studies"] if s["nct_id"] == "NCT01610440")
    assert len(full["eligibility"]) == full["eligibility_chars"] == 855
    over = next(w for w in env["warnings"] if w.startswith("--full-eligibility produced "))
    assert "over the 50000-character guide" in over and "lower --limit" in over


def test_F38_full_eligibility_quiet_when_small(monkeypatch, capsys):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp("f38_nct01610440.json"))
    assert cli.main(["trials", "Duchenne muscular dystrophy", "--full-eligibility", "--status", "ANY", "--json"]) == 0
    env = json.loads(capsys.readouterr().out)
    assert not any(w.startswith("--full-eligibility produced ") for w in env["warnings"])


def test_F38_chictr_note_on_every_query(monkeypatch, capsys):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp("f38_dravet_recruiting.json"))
    assert cli.main(["trials", "Dravet syndrome"]) == 0
    cap = capsys.readouterr()
    assert "ChiCTR (chictr.org.cn), which is not searched here" in cap.out  # no --country at all
    assert "ChiCTR" not in cap.err  # only a note, not a warning, outside China/Hong Kong/Taiwan


@pytest.mark.parametrize("country,warns", [("China", True), ("中国", True), ("Hong Kong", True), ("Taiwan", True),
                                           ("United States", False), (None, False)])
def test_F38_chictr_warning_for_chinese_registries(monkeypatch, capsys, country, warns):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp("f38_dmd_china.json"))
    argv = ["trials", "Duchenne muscular dystrophy", "--status", "ANY"]
    if country:
        argv += ["--country", country]
    assert cli.main(argv) == 0
    cap = capsys.readouterr()
    assert "ChiCTR (chictr.org.cn), which is not searched here" in cap.out  # the note, on every query
    assert ("ChiCTR" in cap.err) is warns  # raised to a warning only for the Chinese registries


@pytest.mark.live
def test_F38_live_fields_are_accepted_by_the_api():
    """An unknown field name makes v2 answer 400, so a 200 here is the field list being valid."""
    out = ctgov.search("Dravet syndrome", status="ANY", limit=1)
    assert out.result["returned"] == 1
    assert "StatusVerifiedDate" in ctgov.FIELDS
    # CP1-8: the relevance check reads these; a 200 above is the API accepting them
    assert ctgov.FIELDS.endswith("OfficialTitle,Keyword,ConditionMeshTerm")


@pytest.mark.live
def test_F38_live_eligibility_and_contacts_present():
    r = ctgov.search("Dravet syndrome", status="RECRUITING", limit=10).result
    assert r["eligibility_budget"] == ctgov.ELIGIBILITY_BUDGET
    with_elig = [s for s in r["studies"] if s["eligibility"]]
    assert len(with_elig) >= 5
    assert all(len(s["eligibility"]) <= ctgov.ELIGIBILITY_BUDGET + 1 for s in r["studies"] if s["eligibility"])
    cut = [s for s in r["studies"] if s["eligibility_truncated"]]
    assert cut and all(s["eligibility_chars"] > ctgov.ELIGIBILITY_BUDGET for s in cut)
    assert any(s["contacts"] and s["contacts"][0]["name"] for s in r["studies"])
    assert any(s["sex"] for s in r["studies"])


@pytest.mark.live
def test_F38_live_full_eligibility_is_whole():
    r = ctgov.search("Dravet syndrome", status="RECRUITING", limit=10, full_eligibility=True).result
    assert r["eligibility_budget"] is None
    assert all(s["eligibility_truncated"] is False for s in r["studies"])
    assert any(len(s["eligibility"]) > ctgov.ELIGIBILITY_BUDGET for s in r["studies"] if s["eligibility"])
    assert all(len(s["eligibility"]) == s["eligibility_chars"] for s in r["studies"] if s["eligibility"])


@pytest.mark.live
def test_F38_live_unknown_stem_cell_trial_is_flagged():
    """NCT01610440: overall UNKNOWN since 2012-11-30, one Chinese site still reading RECRUITING."""
    out = ctgov.search("Duchenne muscular dystrophy", country="China", status="ANY", limit=50)
    s = next(x for x in out.result["studies"] if x["nct_id"] == "NCT01610440")
    assert s["status"] == "UNKNOWN"
    assert any(f.startswith("stale: last update 2012-11-30") for f in s["status_flags"])
    assert "status UNKNOWN: the sponsor has not verified this record" in s["status_flags"]
    assert "site status conflicts with overall status" in s["status_flags"]
    assert s["sites"][0]["status"] == "RECRUITING"
    assert any("may be out of date" in w and "NCT01610440" in w for w in out.warnings)


@pytest.mark.live
def test_live_recruiting_dravet():
    out = ctgov.search("Dravet syndrome", status="RECRUITING", limit=20)
    r = out.result
    assert r["total"] >= 3 and r["studies"]
    assert all(s["status"] == "RECRUITING" for s in r["studies"])
    assert all(s["nct_id"].startswith("NCT") for s in r["studies"])


@pytest.mark.live
def test_live_country_filter_is_exact():
    r = ctgov.search("Duchenne muscular dystrophy", country="China", status="ANY", limit=50).result
    assert r["total"] >= 10
    assert all("China" in s["countries"] and s["sites_in_country"] >= 1 for s in r["studies"])


# ---------------------------------------------------------------- E-3: country spellings

@pytest.mark.parametrize("a,b", [("Korea, Republic of", "South Korea"), ("korea", "South Korea"),
                                 ("Iran, Islamic Republic of", "Iran"), ("Turkey", "Turkey (Türkiye)"),
                                 ("Türkiye", "Turkey (Türkiye)"), ("Russian Federation", "Russia"),
                                 ("Czech Republic", "Czechia"), ("Viet Nam", "Vietnam")])
def test_e_3_country_spellings_compare_equal(a, b):
    assert ctgov.country_key(a) == ctgov.country_key(b)
    assert ctgov.country_key("China") != ctgov.country_key("Taiwan")


@pytest.mark.parametrize("fixture,country", [("e_3_dmd_korea.json", "Korea"), ("e_3_dmd_korea.json", "South Korea"),
                                             ("e_3_thal_iran.json", "Iran"), ("e_3_fmf_turkey.json", "Turkey"),
                                             ("e_3_fmf_turkey.json", "Türkiye")])
def test_e_3_sites_are_found_whatever_the_country_spelling(fixture, country, monkeypatch):
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp(fixture))
    cond = {"e_3_dmd_korea.json": "Duchenne muscular dystrophy", "e_3_thal_iran.json": "thalassemia",
            "e_3_fmf_turkey.json": "familial Mediterranean fever"}[fixture]
    out = ctgov.search(cond, country=country, status="ANY", limit=4)
    studies = out.result["studies"]
    assert studies and all(s["sites_in_country"] >= 1 and s["sites"] for s in studies), \
        [(s["nct_id"], s["countries"], s["sites_in_country"]) for s in studies]
    assert not any("none of their sites reads as that country" in w for w in out.warnings)


def test_e_3_a_site_mismatch_is_warned_not_silent(monkeypatch):
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp("e_3_thal_iran.json"))
    out = ctgov.search("thalassemia", country="Mongolia", status="ANY", limit=4)
    assert any("none of their sites reads as that country" in w for w in out.warnings)


# ---------------------------------------------------------------- CP1-8: relevance

def test_cp1_8_trials_that_do_not_name_the_condition_are_set_aside(monkeypatch, capsys):
    seen = {}

    def fake(url, source, params=None, **kw):
        seen.update(params)
        return resp("cp1_8_sma_china.json")

    monkeypatch.setattr(ctgov, "get_json", fake)
    out = ctgov.search("spinal muscular atrophy", country="China", status="ANY", limit=20)
    kept = {s["nct_id"] for s in out.result["studies"]}
    gone = {f["nct_id"] for f in out.result["filtered"]}
    assert "NCT07190300" in gone and "NCT02348281" in gone and "NCT02294461" in gone and "NCT06628583" in gone
    assert {"NCT04089566", "NCT05614531", "NCT06971094"} <= kept
    assert "NCT05631418" in kept  # condition 'Recruitment', but its title and keywords name SMA
    assert all(s.get("relevance") for s in out.result["studies"])
    assert any("NCT07190300" in w and "set aside" in w for w in out.warnings)
    assert seen["pageSize"] == 40  # asks for more so the page still fills
    code = cli.main(["--json", "trials", "spinal muscular atrophy", "--country", "China", "--status", "ANY",
                     "--keep-unrelated"])
    env = json.loads(capsys.readouterr().out)
    assert "NCT07190300" in {s["nct_id"] for s in env["result"]["studies"]} and "filtered" not in env["result"]


def test_cp1_8_relevance_rules():
    s = {"conditions": ["Duchenne Muscular Dystrophin (DMD)"], "title": "CRISPR in DMD Patients", "keywords": []}
    assert ctgov.relevance(s, "Duchenne muscular dystrophy")  # spelling slip and initialism
    s = {"conditions": ["Muscular Atrophy, Spinal"], "title": "x"}
    assert ctgov.relevance(s, "spinal muscular atrophy")  # word order
    s = {"conditions": ["Spinal Muscular Atrophies of Childhood"], "title": "x"}
    assert ctgov.relevance(s, "spinal muscular atrophy")  # plural
    s = {"conditions": ["Breast Cancer"], "keywords": ["androgen receptor"],
         "condition_mesh": ["Bulbo-Spinal Atrophy, X-Linked"], "title": "Bicalutamide"}
    assert ctgov.relevance(s, "spinal muscular atrophy") is None
    s = {"conditions": ["Dravet Syndrome"], "title": "x"}
    assert ctgov.relevance(s, "Dravet syndrome")


def test_cp1_8_an_acronym_query_is_not_filtered(monkeypatch):
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp("cp1_8_sma_china.json"))
    out = ctgov.search("SMA", country="China", status="ANY", limit=20)
    assert "filtered" not in out.result and out.result["relevance_check"].startswith("not applied")


def test_chictr_note_says_what_blocks_it(monkeypatch, capsys):
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp("dmd_china.json"))
    cli.main(["trials", "Duchenne muscular dystrophy", "--country", "China"])
    out = capsys.readouterr().out
    assert "HTTP 405" in out and "WAF" in out and "https://www.chictr.org.cn/searchproj.html" in out


@pytest.mark.live
def test_live_e_3_korea_and_turkey_sites_are_listed():
    for cond, country in (("Duchenne muscular dystrophy", "Korea"), ("familial Mediterranean fever", "Turkey")):
        out = ctgov.search(cond, country=country, status="ANY", limit=5)
        assert out.result["studies"] and all(s["sites_in_country"] >= 1 for s in out.result["studies"]), country


@pytest.mark.live
def test_live_cp1_8_sma_china_has_no_prostate_cancer_trial():
    out = ctgov.search("spinal muscular atrophy", country="China", status="ANY", limit=40)
    assert "NCT07190300" not in {s["nct_id"] for s in out.result["studies"]}


# ---------------------------------------------------------------- adversarial review (W4 round), relevance

def test_rev_p0_4_spelling_accents_and_subtypes():
    s = {"conditions": ["Hemophilia A"], "title": "x"}
    assert ctgov.relevance(s, "haemophilia A")
    assert ctgov.relevance({"conditions": ["Haemophilia B"], "title": "x"}, "haemophilia A") is None
    assert ctgov.relevance({"conditions": ["Sjogren's Syndrome"], "title": "x"}, "Sjögren syndrome")
    assert ctgov.relevance({"conditions": ["Tumor"], "title": "x"}, "tumour")


def test_rev_p0_4_full_lists_are_read_and_mesh_rescues(monkeypatch):
    def study(nct, conds, mesh):
        return {"protocolSection": {"identificationModule": {"nctId": nct, "briefTitle": "t"},
                                    "conditionsModule": {"conditions": conds}},
                "derivedSection": {"conditionBrowseModule": {"meshes": [{"term": m} for m in mesh]}}}

    many = [f"C{i}" for i in range(7)]
    data = {"totalCount": 4, "studies": [
        study("NCT1", ["Hunter Syndrome"], ["Mucopolysaccharidosis II"]),
        study("NCT2", ["Hunter syndrome"], ["Mucopolysaccharidosis II"]),
        study("NCT3", ["MPS II"], ["Mucopolysaccharidosis II"]),             # rescued by the shared MeSH term
        study("NCT4", many + ["Hunter syndrome"], []),                        # its 8th condition names it
        study("NCT5", ["Breast Cancer"], ["Breast Neoplasms"])]}
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp(data))
    out = ctgov.search("Hunter syndrome", status="ANY", limit=10)
    kept = {s["nct_id"] for s in out.result["studies"]}
    assert kept == {"NCT1", "NCT2", "NCT3", "NCT4"} and [f["nct_id"] for f in out.result["filtered"]] == ["NCT5"]
    assert all("_relevance" not in s and "keywords" not in s for s in out.result["studies"])
    assert any("possibly unrelated" in w and "may still be relevant" in w for w in out.warnings)


def test_rev_p0_4_nothing_is_set_aside_when_most_records_word_it_otherwise(monkeypatch):
    def study(nct, conds):
        return {"protocolSection": {"identificationModule": {"nctId": nct, "briefTitle": "t"},
                                    "conditionsModule": {"conditions": conds}}}

    data = {"totalCount": 3, "studies": [study("NCT1", ["Amyotrophic Lateral Sclerosis"]),
                                         study("NCT2", ["ALS"]), study("NCT3", ["Lou Gehrig's disease"])]}
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp(data))
    out = ctgov.search("Lou Gehrig disease", status="ANY", limit=10)
    assert out.result["returned"] == 3 and "filtered" not in out.result
    assert out.result["relevance_check"].startswith("not applied: only 1 of 3")


@pytest.mark.parametrize("cond", ["MPS II", "SMA type 1", "CDKL5 deficiency"])
def test_rev_p0_4_a_query_with_an_acronym_is_not_filtered(cond, monkeypatch):
    monkeypatch.setattr(ctgov, "get_json", lambda url, source, params=None, **kw: resp("cp1_8_sma_china.json"))
    out = ctgov.search(cond, country="China", status="ANY", limit=20)
    assert "filtered" not in out.result


def test_rev_p2_query_syntax_is_neutralised_and_congo_macao_spellings():
    p = ctgov._params('Duchenne "muscular" dystrophy AND (x)', 'a\\\\b', None, "ANY", 5)
    assert '"' not in p["query.cond"] and " AND " not in p["query.cond"] and "\\\\" not in p["query.term"]
    assert ctgov.country_key("Congo") == ctgov.country_key("Republic of the Congo")
    assert ctgov.country_key("Congo") != ctgov.country_key("Congo, The Democratic Republic of the")
    assert ctgov.country_key("Macao") == ctgov.country_key("Macau")
