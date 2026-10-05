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
    assert matched == [{"as": "target", "id": "ENSG00000001626", "name": "CFTR", "how": "exact gene symbol",
                        "exact": True}]
    assert out.result["search_hits"][0]["entity"] == "target"


def test_resolve_disease_name(monkeypatch):
    matched, _ = _resolve("dravet syndrome", monkeypatch)
    assert matched == [{"as": "disease", "id": "MONDO_0100135", "name": "Dravet syndrome",
                        "how": "exact disease name", "exact": True}]


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


# -------------------------------------------- P1g: approval status is per jurisdiction and current

def _dmd_rows():
    data = load("opentargets/p1g_dmd_drugs.json")["data"]
    return (data["disease"]["drugAndClinicalCandidates"]["rows"], data)


def test_P1g_jurisdiction_mapping_and_stage_reading():
    reports = [
        {"origin": "REGULATORY_AGENCY", "source": "FDA", "clinicalStage": "APPROVAL", "url": "u1"},
        {"origin": "REGULATORY_AGENCY", "source": "EMA Human Drugs", "clinicalStage": "WITHDRAWAL", "url": "u2"},
        {"origin": "REGULATORY_AGENCY", "source": "PMDA", "clinicalStage": "APPROVAL", "url": "u3"},
        {"origin": "DRUG_LABEL", "source": "DailyMed", "clinicalStage": "APPROVAL", "url": "u4"},
        {"origin": "CLINICAL_TRIAL", "source": "ClinicalTrials.gov", "clinicalStage": "PHASE_3", "url": "u5"},
    ]
    reg = ot.regulatory_status(reports)
    assert reg["approved_in"] == ["FDA (United States)", "PMDA (Japan)"]
    assert reg["withdrawn_or_suspended_in"] == ["EMA (European Union)"]
    assert reg["label_only_in"] == ["United States (DailyMed label archive)"]
    assert {j["jurisdiction"] for j in reg["by_jurisdiction"]} == {
        "FDA (United States)", "EMA (European Union)", "PMDA (Japan)"}
    withdrawn = next(j for j in reg["by_jurisdiction"] if j["stage"] == "WITHDRAWAL")
    assert withdrawn["reading"] == "authorisation withdrawn — NOT currently approved"
    assert "WITHDRAWN/SUSPENDED: EMA (European Union)" in reg["headline"]


def test_P1g_a_withdrawal_outranks_an_approval_for_the_same_agency():
    reg = ot.regulatory_status([
        {"origin": "REGULATORY_AGENCY", "source": "EMA Human Drugs", "clinicalStage": "APPROVAL", "url": "a"},
        {"origin": "REGULATORY_AGENCY", "source": "EMA Human Drugs", "clinicalStage": "WITHDRAWAL", "url": "b"},
    ])
    assert reg["approved_in"] == [] and reg["withdrawn_or_suspended_in"] == ["EMA (European Union)"]


def test_P1g_ataluren_is_not_presented_as_approved():
    """The row's maxClinicalStage is APPROVAL while its only agency report is an EMA withdrawal."""
    rows, _ = _dmd_rows()
    row = next(r for r in rows if (r["drug"] or {}).get("name") == "ATALUREN")
    parsed = ot.parse_drug_row(row)
    assert parsed["stage"] == "APPROVAL"  # the raw field is unchanged ...
    assert parsed["regulatory"]["approved_in"] == []   # ... and contradicted per jurisdiction
    assert parsed["regulatory"]["withdrawn_or_suspended_in"] == ["EMA (European Union)"]
    assert "do not read this as an approved therapy" in parsed["stage_warning"]
    assert "highest stage this drug reached FOR THIS DISEASE" in parsed["stage_meaning"]


def test_P1g_a_genuinely_approved_drug_says_where():
    rows, _ = _dmd_rows()
    row = next(r for r in rows if (r["drug"] or {}).get("name") == "CASIMERSEN")
    parsed = ot.parse_drug_row(row)
    assert parsed["regulatory"]["approved_in"] == ["FDA (United States)"]
    assert "stage_warning" not in parsed
    vam = ot.parse_drug_row(next(r for r in rows if (r["drug"] or {}).get("name") == "VAMOROLONE"))
    assert vam["regulatory"]["approved_in"] == ["EMA (European Union)", "FDA (United States)"]


def test_P1g_the_list_is_labelled_as_what_the_source_holds():
    _, data = _dmd_rows()
    drugs = ot._drugs(data["disease"]["drugAndClinicalCandidates"], 40, False)
    assert "NOT a list of approved therapies" in drugs["what_this_is"]
    assert "ATALUREN" in drugs["withdrawn_or_suspended"]


def test_P1g_orphan_designations_missing_from_the_drug_list_are_named():
    """Elevidys has an EU orphan designation for DMD and no row in Open Targets' 62 drugs."""
    rows, data = _dmd_rows()
    parsed = ot._drugs(data["disease"]["drugAndClinicalCandidates"], 40, False)["rows"]
    designations = [
        {"substance": "adeno-associated virus serotype rh74 containing the human micro-dystrophin gene",
         "medicine": "Elevidys", "status": "Positive", "eu_number": "EU/3/18/x", "url": "u"},
        {"substance": "ATALUREN", "medicine": None, "status": "Positive", "eu_number": "EU/3/05/y", "url": "v"},
    ]
    missing = therapy_cmd._designations_not_in_drug_list(designations, parsed)
    names = [m["medicine"] or m["substance"] for m in missing]
    assert "Elevidys" in names
    assert "ATALUREN" not in [m["substance"] for m in missing]  # it IS in the list, so not a gap
    assert therapy_cmd._gap_label(missing[0]) == "Elevidys"


def test_P1g_gap_label_shortens_a_chemical_description():
    label = therapy_cmd._gap_label({"substance": "palmitoyl-conjugated tricyclo-DNA antisense oligonucleotide "
                                                 "5'-Palm-C6-*GGA GAT GgC AGT TTC-3", "medicine": None})
    assert len(label) <= 61 and label.endswith("…")


def test_P1g_rendering_shows_the_jurisdictions_and_the_warning():
    rows, _ = _dmd_rows()
    row = ot.parse_drug_row(next(r for r in rows if (r["drug"] or {}).get("name") == "ATALUREN"))
    line = therapy_cmd._drug_line(row, False)
    assert "status by jurisdiction:" in line
    assert "EMA (European Union): WITHDRAWAL — authorisation withdrawn — NOT currently approved" in line
    assert "! Open Targets gives this row maxClinicalStage APPROVAL" in line
    assert "any indication:" in line or row["drug_max_stage"] == row["stage"]


# ---------------------------------------------------------------- F13: a fuzzy disease match warns

def test_F13_an_exact_disease_name_raises_no_warning(monkeypatch):
    monkeypatch.setattr(ot, "post_json", ot_router())
    out = Outcome({})
    matched = therapy_cmd._resolve("Dravet Syndrome", out)
    assert matched[0]["exact"] is True and out.warnings == []


def test_F13_first_hit_resolution_is_a_warning_not_a_footnote(monkeypatch):
    hits = {"data": {"search": {"total": 4, "hits": [
        {"id": "MONDO_0019079", "name": "proximal spinal muscular atrophy", "entity": "disease", "description": ""},
        {"id": "EFO_0010970", "name": "severe malarial anemia", "entity": "disease", "description": ""},
        {"id": "MONDO_0001516", "name": "spinal muscular atrophy", "entity": "disease", "description": ""},
    ]}}}

    def fake(url, payload, source, **kw):
        return resp(hits, url)

    monkeypatch.setattr(ot, "post_json", fake)
    out = Outcome({})
    matched = therapy_cmd._resolve("SMA", out)
    assert matched[0]["id"] == "MONDO_0019079" and matched[0]["exact"] is False
    assert len(out.warnings) == 1
    w = out.warnings[0]
    assert "'SMA' is not an exact Open Targets disease name" in w
    assert "proximal spinal muscular atrophy" in w
    assert "spinal muscular atrophy (MONDO_0001516)" in w
    assert "`zebra disease` resolves a name to an id without guessing" in w


# ---------------------------------------------------------------- live

@pytest.mark.live
def test_live_P1g_ataluren_for_dmd_is_withdrawn_not_approved():
    d = ot.disease("MONDO_0010679").result
    row = next(r for r in d["drugs"]["rows"] if r["drug"] == "ATALUREN")
    assert row["stage"] == "APPROVAL"
    assert row["regulatory"]["withdrawn_or_suspended_in"] == ["EMA (European Union)"]
    assert row["regulatory"]["approved_in"] == []
    assert "ATALUREN" in d["drugs"]["withdrawn_or_suspended"]


@pytest.mark.live
def test_live_F13_sma_warns_that_the_name_was_not_exact(capsys):
    code = cli.main(["--json", "therapy", "SMA"])
    env = json.loads(capsys.readouterr().out)
    assert code == 0
    assert any("is not an exact Open Targets disease name" in w for w in env["warnings"])
    assert env["result"]["matched"][0]["exact"] is False


# ---------------------------------------------------------------- F10: never serve a cached error

def test_F10_a_cached_graphql_error_is_refetched_once(monkeypatch):
    """B-P1-1's first repro: a 200 carrying `errors` was cached for 7 days and raised on every call."""
    ttls = []

    def fake_post(url, payload, source, **kw):
        ttls.append(kw.get("cache_ttl"))
        cached = len(ttls) == 1
        body = ('{"errors":[{"message":"upstream timeout"}]}' if cached
                else '{"data":{"disease":{"id":"MONDO_1","name":"X"}}}')
        return Response(url, 200, body, "2026-10-06T00:00:00+00:00", cached)

    monkeypatch.setattr(ot, "post_json", fake_post)
    resp, data = ot.gql(ot.DISEASE_Q, {"id": "MONDO_1", "nt": 1}, "MONDO_1")
    assert ttls == [7 * 86400, 0]  # the second call bypasses the cache
    assert data["disease"]["name"] == "X"


def test_F10_a_fresh_graphql_error_is_raised_once_and_names_the_record(monkeypatch):
    calls = []

    def fake_post(url, payload, source, **kw):
        calls.append(1)
        return Response(url, 200, '{"errors":[{"message":"upstream timeout"}]}',
                        "2026-10-06T00:00:00+00:00", False)

    monkeypatch.setattr(ot, "post_json", fake_post)
    with pytest.raises(SourceError) as err:
        ot.gql(ot.DISEASE_Q, {"id": "MONDO_1", "nt": 1}, "MONDO_1")
    assert len(calls) == 1  # a fresh error is not retried
    assert "GraphQL error for MONDO_1" in err.value.message and "upstream timeout" in err.value.message
