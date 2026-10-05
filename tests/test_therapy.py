"""Open Targets source, EMA orphan designations and `zebra therapy` (offline on captured responses, plus live)."""

from __future__ import annotations

import copy
import json
import os

import pytest

from zebra import cli
from zebra.commands import therapy as therapy_cmd
from zebra.core import Outcome, UsageError
from zebra.http import Response, SourceError
from zebra.sources import opentargets as ot
from zebra.sources import orphan

FIX = os.path.join(os.path.dirname(__file__), "fixtures")


def load(rel):
    with open(os.path.join(FIX, rel), encoding="utf-8") as fh:
        return json.load(fh)


def resp(obj, url=ot.API):
    return Response(url, 200, json.dumps(obj), "2026-10-05T00:00:00+00:00", False)


def ot_router(overrides=None):
    """Replay captured GraphQL answers keyed by the query kind and variables."""
    overrides = overrides or {}

    def fake(url, payload, source, **kw):
        q, v = payload["query"], payload["variables"]
        if "search(" in q:
            name = {"dravet syndrome": "search_dravet", "cftr": "search_cftr"}[v["q"].lower()]
        elif "disease(efoId" in q:
            name = {"MONDO_0100135": "disease_dravet"}.get(v["id"], "disease_null")
        else:
            name = "target_cftr"
        data = overrides.get(name) or load(f"opentargets/{name}.json")
        return resp(data, url)
    return fake


def ema_fake(url, source, **kw):
    return Response(url, 200, json.dumps(load("ema/orphan_sample.json")), "2026-10-05T00:00:00+00:00", False)


# --------------------------------------------------------------------------- Open Targets

def test_normalise_disease_id():
    assert ot.normalise_disease_id("MONDO:0100135") == "MONDO_0100135"
    assert ot.normalise_disease_id("mondo_0100135") == "MONDO_0100135"
    assert ot.normalise_disease_id("Orphanet_98896") == "Orphanet_98896"
    assert ot.normalise_disease_id("OMIM:607208") is None
    assert ot.normalise_disease_id("Dravet syndrome") is None


def test_disease_drugs_ordered_by_stage_with_regulatory_evidence(monkeypatch):
    monkeypatch.setattr(ot, "post_json", ot_router())
    out = ot.disease("MONDO:0100135")
    d = out.result
    assert d["id"] == "MONDO_0100135" and d["name"] == "Dravet syndrome" and d["data_version"]
    names = [x["drug"] for x in d["drugs"]["rows"]]
    assert names[:2] == ["FENFLURAMINE", "STIRIPENTOL"]  # APPROVAL first, then by name
    assert [x["stage"] for x in d["drugs"]["rows"]][-1] == "PHASE_2"
    sti = next(x for x in d["drugs"]["rows"] if x["drug"] == "STIRIPENTOL")
    assert sti["evidence"][0]["source"] == "PMDA" and sti["evidence"][0]["stage"] == "APPROVAL"
    assert sti["mechanisms"] == [] and sti["chembl_id"] == "CHEMBL1983350"
    ata = next(x for x in d["drugs"]["rows"] if x["drug"] == "ATALUREN")
    assert len(ata["mechanisms"][0]["targets"]) == 6 and ata["mechanisms"][0]["targets"][-1].startswith("+")
    cloba = next(x for x in d["drugs"]["rows"] if x["drug"] == "CLOBAZAM")
    assert cloba["stage"] == "PHASE_3" and cloba["drug_max_stage"] == "APPROVAL"
    assert d["top_targets"][0]["symbol"] == "SCN1A" and 0 < d["top_targets"][0]["score"] <= 1
    assert any(x.startswith("GARD:") for x in d["xrefs"]) and "Dravet" in d["synonyms"]
    src = out.sources[0]
    assert src["db"] == "Open Targets" and src["record"] == "MONDO_0100135"
    assert src["url"] == "https://platform.opentargets.org/disease/MONDO_0100135"


def test_null_drug_and_mechanism_rows_are_guarded(monkeypatch):
    data = copy.deepcopy(load("opentargets/disease_dravet.json"))
    rows = data["data"]["disease"]["drugAndClinicalCandidates"]["rows"]
    rows[0]["drug"] = None
    rows[1]["drug"]["mechanismsOfAction"] = None
    monkeypatch.setattr(ot, "post_json", ot_router({"disease_dravet": data}))
    d = ot.disease("MONDO_0100135").result
    assert len(d["drugs"]["rows"]) == len(rows) - 1
    assert all(isinstance(x["mechanisms"], list) for x in d["drugs"]["rows"])


def test_target_tractability_and_drugs(monkeypatch):
    monkeypatch.setattr(ot, "post_json", ot_router())
    t = ot.target("ENSG00000001626").result
    assert t["symbol"] == "CFTR" and t["biotype"] == "protein_coding"
    assert "Approved Drug" in t["tractability"]["small_molecule"]
    assert set(t["tractability"]) == {"small_molecule", "antibody", "protac", "other_modalities"}
    iva = next(x for x in t["drugs"]["rows"] if x["drug"] == "IVACAFTOR")
    assert iva["drug_max_stage"] == "APPROVAL" and "cystic fibrosis" in iva["indications"]
    assert iva["evidence"][0]["source"] == "EMA Human Drugs"  # regulatory reports first
    assert t["drugs"]["rows"][0]["drug"] == "DEUTIVACAFTOR"
    assert t["top_diseases"][0]["name"] == "cystic fibrosis"


def test_graphql_error_raises_source_error(monkeypatch):
    monkeypatch.setattr(ot, "post_json", lambda url, payload, source, **kw: resp(load("opentargets/graphql_error.json"), url))
    with pytest.raises(SourceError) as err:
        ot.disease("MONDO_0010679")
    assert "knownDrugs" in err.value.message


def test_unknown_disease_is_a_warning(monkeypatch):
    monkeypatch.setattr(ot, "post_json", ot_router())
    out = ot.disease("MONDO_0100079")
    assert out.result is None and "no disease" in out.warnings[0]


# --------------------------------------------------------------------------- query resolution

def _resolve(query, monkeypatch):
    monkeypatch.setattr(ot, "post_json", ot_router())
    out = Outcome({})
    return therapy_cmd._resolve(query, out), out


def test_resolve_gene_symbol_as_target(monkeypatch):
    matched, out = _resolve("CFTR", monkeypatch)
    assert matched == [{"as": "target", "id": "ENSG00000001626", "name": "CFTR", "how": "exact gene symbol"}]
    assert out.result["search_hits"][0]["entity"] == "target"


def test_resolve_disease_name(monkeypatch):
    matched, _ = _resolve("dravet syndrome", monkeypatch)
    assert matched == [{"as": "disease", "id": "MONDO_0100135", "name": "Dravet syndrome", "how": "exact disease name"}]


def test_resolve_ids_and_refuses_xref_numbers(monkeypatch):
    assert _resolve("MONDO:0100135", monkeypatch)[0][0]["id"] == "MONDO_0100135"
    assert _resolve("ENSG00000001626", monkeypatch)[0][0]["as"] == "target"
    with pytest.raises(UsageError):
        _resolve("OMIM:607208", monkeypatch)
    with pytest.raises(UsageError):
        _resolve("ORPHA:33069", monkeypatch)


# --------------------------------------------------------------------------- EMA orphan designations

def test_match_terms_skips_abbreviations():
    assert orphan.match_terms("Dravet syndrome", ["DS", "Dravet", "SME", "SMEB", "dravet syndrome"]) == ["Dravet syndrome", "Dravet"]


def test_ema_designations_word_bounded(monkeypatch):
    monkeypatch.setattr(orphan, "get_json", ema_fake)
    o = orphan.ema_designations(["spinal muscular atrophy"]).result
    assert o["count"] >= 5
    assert not any("bulbar" in d["intended_use"] for d in o["designations"])  # SBMA is another disease
    statuses = [d["status"] for d in o["designations"]]
    assert statuses == sorted(statuses, key=lambda s: s != "Positive")
    d = orphan.ema_designations(["Dravet syndrome"]).result
    assert d["count"] == 6 and all(x["url"].startswith("https://www.ema.europa.eu/") for x in d["designations"])
    assert d["dataset_timestamp"] and d["note"]


def test_date_key_sorts_newest_first():
    keys = sorted(["01/02/2013", "22/05/2023", "15/10/2014"], key=orphan._date_key)
    assert keys == ["22/05/2023", "15/10/2014", "01/02/2013"]


# --------------------------------------------------------------------------- zebra therapy (offline end to end)

def test_cli_therapy_disease(monkeypatch, capsys):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    monkeypatch.setattr(ot, "post_json", ot_router())
    monkeypatch.setattr(orphan, "get_json", ema_fake)
    assert cli.main(["therapy", "Dravet", "syndrome", "--json"]) == 0
    env = json.loads(capsys.readouterr().out)
    r = env["result"]
    assert r["matched"][0]["as"] == "disease" and r["disease"]["name"] == "Dravet syndrome" and r["target"] is None
    assert r["orphan_designations"]["count"] == 6
    assert any("FDA" in n for n in r["notes"])
    assert {"Open Targets search", "Open Targets", "EMA orphan designations"} <= {s["db"] for s in env["sources"]}


def test_cli_therapy_target(monkeypatch, capsys):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    monkeypatch.setattr(ot, "post_json", ot_router())
    monkeypatch.setattr(orphan, "get_json", ema_fake)
    assert cli.main(["therapy", "CFTR"]) == 0
    text = capsys.readouterr().out
    assert "== Target CFTR ENSG00000001626" in text and "tractability, small molecule: Approved Drug" in text
    assert "IVACAFTOR" in text


def test_cli_therapy_ema_down_is_a_warning(monkeypatch, capsys):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    monkeypatch.setattr(ot, "post_json", ot_router())

    def down(url, source, **kw):
        raise SourceError("EMA orphan designations", url, None, "network error: timed out")

    monkeypatch.setattr(orphan, "get_json", down)
    assert cli.main(["therapy", "Dravet syndrome", "--json"]) == 0
    env = json.loads(capsys.readouterr().out)
    assert env["result"]["disease"] and env["result"]["orphan_designations"] is None
    assert any(w.startswith("EMA orphan designations unavailable") for w in env["warnings"])


# --------------------------------------------------------------------------- live

def _drug_names(block):
    return {x["drug"] for x in block["drugs"]["rows"]}


@pytest.mark.live
def test_live_therapy_dravet(capsys, monkeypatch):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    out = ot.disease(ot.search("Dravet syndrome", ("disease",), 3).result["hits"][0]["id"])
    names = _drug_names(out.result)
    assert {"STIRIPENTOL", "CANNABIDIOL"} <= names and any(n.startswith("FENFLURAMINE") for n in names)


@pytest.mark.live
def test_live_therapy_cftr_target():
    t = ot.target("ENSG00000001626").result
    assert "IVACAFTOR" in _drug_names(t) and "Approved Drug" in t["tractability"]["small_molecule"]


@pytest.mark.live
def test_live_therapy_sma():
    d = ot.disease("MONDO_0001516").result
    names = _drug_names(d)
    assert {"RISDIPLAM", "ONASEMNOGENE ABEPARVOVEC"} <= names and any(n.startswith("NUSINERSEN") for n in names)


@pytest.mark.live
def test_live_ema_orphan_register():
    o = orphan.ema_designations(["Dravet syndrome"]).result
    assert o["records_scanned"] > 3000 and o["count"] >= 3
