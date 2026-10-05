"""ClinicalTrials.gov v2 source and `zebra trials` (offline on captured responses, plus live)."""

from __future__ import annotations

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
    assert ctgov.country_name("south korea") == "Korea, Republic of"
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
    assert env["query"] == {"condition": "Duchenne muscular dystrophy", "term": None, "country": "China", "status": "ANY", "limit": 3}
    assert seen["filter.advanced"] == 'AREA[LocationCountry]"China"' and "filter.overallStatus" not in seen
    assert env["sources"][0]["db"] == "ClinicalTrials.gov" and env["sources"][0]["url"].startswith("https://clinicaltrials.gov/")


def test_cli_trials_bad_status(monkeypatch, capsys):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    assert cli.main(["trials", "Dravet syndrome", "--status", "OPEN", "--json"]) == 2
    assert json.loads(capsys.readouterr().out)["error"]["type"] == "UsageError"


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
