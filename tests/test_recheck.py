"""zebra case recheck: what changed since the last check (CP1-12). Offline: the sources are stubbed."""

from __future__ import annotations

import json

from zebra import case as case_mod
from zebra import recheck
from zebra.core import Outcome


def _case(tmp_path, identifiers=()):
    d = str(tmp_path / "c")
    case_mod.init(d, title="t")
    with case_mod.editing(d) as data:
        data["variants"].append({"id": "v1", "kind": "small", "gene": "SCN1A", "hgvs_c": "NM_001165963.4:c.2134C>T"})
        data["variants"].append({"id": "v2", "kind": "cnv", "region": "chr22:18900000-21500000", "cnv_type": "loss"})
        data["hypotheses"].append({"id": "h1", "disease": "Dravet syndrome", "status": "leading"})
        data["hypotheses"].append({"id": "h2", "disease": "Excluded thing", "status": "excluded"})
        data["privacy"]["identifiers"] = list(identifiers)
    return d


def _stub(monkeypatch, clin="Pathogenic", stars=2, trials=("NCT1", "NCT2"), papers=5):
    monkeypatch.setattr(recheck, "_ask_variant", lambda item: (
        {"clinvar": {"classification": clin, "stars": stars, "url": "u"}, "grpmax_af": None, "gnomad_ac": 0},
        Outcome({}, sources=[{"db": "ClinVar", "record": item["query"]}])))
    monkeypatch.setattr(recheck, "_ask_gene", lambda item: (
        {"validity": [{"disease": "Dravet", "mondo": "MONDO:1", "classification": "Definitive", "moi": "AD"}]},
        Outcome([], sources=[{"db": "ClinGen", "record": item["symbol"]}])))
    monkeypatch.setattr(recheck, "_ask_trials", lambda item: (
        {"recruiting": sorted(trials), "titles": {t: f"title {t}" for t in trials}},
        Outcome({}, sources=[{"db": "ClinicalTrials.gov", "record": item["disease"]}])))
    monkeypatch.setattr(recheck, "_ask_papers", lambda item, since, until: (
        {"since": since.isoformat(), "count": papers, "latest": []},
        Outcome({}, sources=[{"db": "Europe PMC", "record": item["disease"]}])))


def test_cp1_12_plan_asks_only_about_sequence_variants_and_open_hypotheses(tmp_path):
    d = _case(tmp_path)
    todo = recheck.plan(case_mod.load(d))
    assert [v["id"] for v in todo["variants"]] == ["v1"]  # the CNV is not a card query
    assert [g["symbol"] for g in todo["genes"]] == ["SCN1A"]
    assert [h["id"] for h in todo["hypotheses"]] == ["h1"]  # excluded hypotheses are not rechecked


def test_cp1_12_first_run_is_a_baseline_then_changes_are_reported(tmp_path, monkeypatch):
    d = _case(tmp_path)
    _stub(monkeypatch, clin="Uncertain significance", stars=1, trials=("NCT1",))
    out, current, previous = recheck.run(d, case_mod.load(d))
    assert out.result["first_check"] and out.result["changes"] == [] and previous is None
    recheck.save_snapshot(d, current, previous)
    # the world moves
    _stub(monkeypatch, clin="Pathogenic", stars=2, trials=("NCT1", "NCT2"), papers=3)
    out, current, previous = recheck.run(d, case_mod.load(d))
    whats = {c["what"]: c for c in out.result["changes"]}
    assert whats["ClinVar"]["before"].startswith("Uncertain significance") and whats["ClinVar"]["after"].startswith("Pathogenic")
    assert [t["nct_id"] for t in whats["new recruiting trials"]["trials"]] == ["NCT2"]
    assert whats["new papers"]["count"] == 3
    assert "ClinGen gene-disease validity" not in whats  # unchanged
    assert len(out.sources) == 4  # one per question asked


def test_cp1_12_a_failed_question_costs_only_its_line(tmp_path, monkeypatch):
    d = _case(tmp_path)
    _stub(monkeypatch)

    def broken(item):
        raise RuntimeError("ClinGen down")

    monkeypatch.setattr(recheck, "_ask_gene", broken)
    out, _current, _prev = recheck.run(d, case_mod.load(d))
    assert "gene:SCN1A" in out.result["not_checked"]
    assert any("ClinGen down" in w for w in out.warnings)
    assert "variant:v1" not in out.result["not_checked"]


def test_cp1_12_a_query_holding_a_protected_identifier_is_never_sent(tmp_path, monkeypatch):
    d = _case(tmp_path, identifiers=["Dravet"])  # contrived: the identifier appears in the hypothesis text
    sent = []
    _stub(monkeypatch)
    monkeypatch.setattr(recheck, "_ask_trials", lambda item: sent.append(item) or ({}, Outcome({})))
    out, _c, _p = recheck.run(d, case_mod.load(d))
    assert sent == []
    assert "disease:h1" in out.result["not_checked"]


def test_cp1_12_cli_records_baseline_then_timeline(tmp_path, monkeypatch, capsys):
    from zebra.cli import main

    d = _case(tmp_path)
    _stub(monkeypatch)
    assert main(["--json", "case", "recheck", d]) == 0
    body = json.loads(capsys.readouterr().out)
    assert body["ok"] and body["result"]["first_check"]
    assert (tmp_path / "c" / "evidence" / "recheck.json").exists()
    assert case_mod.load(d)["timeline"][-1]["event"] == "recheck: baseline recorded"
    assert len(case_mod.read_ledger(d)) == 4  # the answers' sources went into the ledger
