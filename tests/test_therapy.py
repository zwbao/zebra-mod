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


@pytest.fixture(autouse=True)
def _no_agency_network(request, monkeypatch):
    """Offline tests do not reach openFDA or the EMA export; tests of the agency check replay their own answers."""
    if "live" in request.keywords or "agency_network" in request.fixturenames:
        return
    from zebra.sources import regulators

    def stub(drug, terms, with_fda_date=True):
        return Outcome({"drug": drug, "terms": list(terms), "FDA": {"status": "unavailable"},
                        "EMA": {"status": "unavailable"}, "approved_in": [], "refused_or_withdrawn_in": [],
                        "headline": "not checked (offline test)"})

    monkeypatch.setattr(regulators, "check", stub)


@pytest.fixture
def agency_network():
    """Marker fixture: this test drives zebra.sources.regulators itself."""
    return True


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
        ttls.append((kw.get("cache_ttl"), kw.get("refresh", False)))
        cached = len(ttls) == 1
        body = ('{"errors":[{"message":"upstream timeout"}]}' if cached
                else '{"data":{"disease":{"id":"MONDO_1","name":"X"}}}')
        return Response(url, 200, body, "2026-10-06T00:00:00+00:00", cached)

    monkeypatch.setattr(ot, "post_json", fake_post)
    resp, data = ot.gql(ot.DISEASE_Q, {"id": "MONDO_1", "nt": 1}, "MONDO_1")
    # the second call skips the cache read; E-6: it keeps a TTL, so the good answer is written back
    assert ttls == [(7 * 86400, False), (7 * 86400, True)]
    assert data["disease"]["name"] == "X"


def test_e_6_recovery_from_a_bad_cached_body_writes_the_good_answer_back(monkeypatch, tmp_path):
    """E-6, end to end through zebra.http: seed a bad cached 200, recover once, then the cache serves the good body."""
    import zebra.http as zh

    calls = []

    class FakeResp:
        status = 200

        def __init__(self, body):
            self.body = body

        def read(self, *a):
            b, self.body = self.body, b""
            return b

        read1 = read

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    class FakeOpener:
        def open(self, req, timeout=None):
            calls.append(req.full_url)
            return FakeResp(b'{"data":{"disease":{"id":"MONDO_1","name":"X"}}}')

    monkeypatch.setattr(zh, "_opener", lambda: FakeOpener())
    payload = {"query": ot.DISEASE_Q, "variables": {"id": "MONDO_1", "nt": 1}}
    key = f"POST {ot.API} {json.dumps(payload)}"
    path = zh._cache_path(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"stored": __import__("time").time(), "status": 200,
                                "text": '{"errors":[{"message":"upstream timeout"}]}',
                                "retrieved_at": "2026-10-01T00:00:00+00:00"}), "utf-8")
    for _ in range(3):
        _resp, data = ot.gql(ot.DISEASE_Q, {"id": "MONDO_1", "nt": 1}, "MONDO_1")
        assert data["disease"]["name"] == "X"
    assert len(calls) == 1  # recovered once; the next two calls were served from the rewritten cache
    assert "upstream timeout" not in path.read_text("utf-8")


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


# ---------------------------------------------------------------- E-4: agency records before "approved"

def _reg_fixture(name):
    with open(os.path.join(FIX, "regulators", name), encoding="utf-8") as fh:
        return fh.read()


def _regulators_router(calls=None):
    from zebra.sources import regulators

    def fake(url, source, params=None, **kw):
        if calls is not None:
            calls.append((url, dict(params or {})))
        q = (params or {}).get("search", "")
        if url == regulators.EMA_MEDICINES:
            return Response(url, 200, _reg_fixture("ema_medicines_subset.json"), "2026-10-06T00:00:00+00:00", False)
        if url == regulators.FDA_DRUGSFDA and "NDA210365" in q:
            return Response(url, 200, _reg_fixture("drugsfda_nda210365.json"), "2026-10-06T00:00:00+00:00", False)
        if url == regulators.FDA_LABEL and "cannabidiol" in q and "Dravet" in q:
            return Response(url, 200, _reg_fixture("label_cannabidiol_dravet.json"), "2026-10-06T00:00:00+00:00", False)
        return Response(url, 404, _reg_fixture("label_not_found.json"), "2026-10-06T00:00:00+00:00", False)

    return fake


def test_e_4_no_agency_record_reads_unknown_not_unapproved():
    row = {"maxClinicalStage": "APPROVAL", "drug": {"id": "CHEMBL190461", "name": "CANNABIDIOL"}, "clinicalReports": [
        {"id": "x", "source": "ClinicalTrials.gov", "clinicalStage": "PHASE_3", "origin": "CLINICAL_TRIAL"}]}
    parsed = ot.parse_drug_row(row)
    assert "approval status unknown here" in parsed["stage_warning"]
    assert "do not read this as an approved therapy" not in parsed["stage_warning"]


def test_e_4_cannabidiol_for_dravet_is_approved_by_the_fda_and_ema_records(agency_network, monkeypatch):
    from zebra.sources import regulators

    calls = []
    monkeypatch.setattr(regulators, "get_json", _regulators_router(calls))
    out = regulators.check("CANNABIDIOL", ["Dravet syndrome"])
    r = out.result
    assert r["FDA"]["status"] == "approved_for_disease" and r["FDA"]["application"] == "NDA210365"
    assert r["FDA"]["approved_on"] == "2018-06-25" and "Epidiolex" in (r["FDA"]["brand"] or "")
    assert r["EMA"]["status"] == "approved_for_disease" and r["EMA"]["medicine"] == "Epidyolex"
    assert r["EMA"]["date"] == "2019-09-19"
    assert r["approved_in"] == ["FDA (United States)", "EMA (European Union)"]
    # the disease is part of the FDA query (Epidiolex is not in a name-only first page)
    assert any("indications_and_usage" in p.get("search", "") for _u, p in calls)
    assert {s["db"] for s in out.sources} >= {"openFDA drug label", "openFDA Drugs@FDA",
                                              "EMA medicines (EPAR register export)"}


def test_e_4_a_refusal_is_reported_and_nothing_found_is_unknown(agency_network, monkeypatch):
    from zebra.sources import regulators

    monkeypatch.setattr(regulators, "get_json", _regulators_router())
    ete = regulators.check("ETEPLIRSEN", ["Duchenne muscular dystrophy"]).result
    assert ete["EMA"]["status"] == "refused" and ete["EMA"]["date"] == "2018-12-06"
    assert ete["FDA"]["status"] == "no_record"
    assert "EMA (European Union) (refused)" in ete["refused_or_withdrawn_in"]
    none = regulators.check("NEWDRUGAMAB", ["Dravet syndrome"]).result
    assert none["approved_in"] == [] and "approval status unknown here" in none["headline"]
    assert "not 'not approved'" in none["headline"]


def test_e_4_therapy_row_drops_the_warning_when_an_agency_approved_it(agency_network, monkeypatch, capsys):
    from zebra.sources import regulators

    monkeypatch.setattr(regulators, "get_json", _regulators_router())
    rows = [ot.parse_drug_row({"maxClinicalStage": "APPROVAL", "drug": {"id": "CHEMBL190461", "name": "CANNABIDIOL"},
                               "clinicalReports": []})]
    assert "approval status unknown here" in rows[0]["stage_warning"]
    out = Outcome({})
    therapy_cmd._agency_checks(rows, ["Dravet syndrome"], out)
    assert "stage_warning" not in rows[0]
    assert rows[0]["agency_check"]["FDA"]["status"] == "approved_for_disease"
    assert rows[0]["agency_check"]["EMA"]["status"] == "approved_for_disease"
    line = therapy_cmd._drug_line(dict(rows[0], evidence=[], mechanisms=[], reports=0, trials=0), False)
    assert "FDA check: approved for disease" in line and "EMA check: approved for disease" in line
    assert "dailymed" in line and "epidyolex" in line


def test_e_4_moiety_strips_salt_words():
    from zebra.sources import regulators

    assert regulators.moiety("FENFLURAMINE HYDROCHLORIDE") == "fenfluramine"
    assert regulators.moiety("NUSINERSEN SODIUM") == "nusinersen"
    assert regulators.moiety("SODIUM PHENYLBUTYRATE") == "phenylbutyrate"
    assert regulators.moiety("sodium") == "sodium"  # nothing left: keep the word


# ---------------------------------------------------------------- E-8: the coverage-gap check

def test_e_8_salt_and_modality_words_do_not_hide_a_gap():
    rows = [{"drug": "VALPROATE SODIUM"}, {"drug": "ONASEMNOGENE ABEPARVOVEC"}]
    assert therapy_cmd._designations_not_in_drug_list([{"substance": "Newdrugamab sodium"}], rows)
    assert therapy_cmd._designations_not_in_drug_list([{"substance": "Otherdrug hydrochloride"}], rows)
    aav = [{"substance": "adeno-associated viral vector serotype 9 containing the human XYZ1 gene"}]
    assert therapy_cmd._designations_not_in_drug_list(aav, rows)
    # a shared distinctive word still counts as present
    assert not therapy_cmd._designations_not_in_drug_list([{"substance": "valproate"}], rows)


def test_e_8_the_gap_warning_says_possibly_missing(monkeypatch, capsys):
    import zebra.commands.therapy as T

    out = Outcome({"query": "x", "matched": [], "disease": None, "target": None, "orphan_designations": None,
                   "notes": []})
    missing = T._designations_not_in_drug_list([{"substance": "Newdrugamab sodium", "status": "Positive"}],
                                               [{"drug": "VALPROATE SODIUM"}])
    assert missing and missing[0]["substance"] == "Newdrugamab sodium"
    src = open(T.__file__, encoding="utf-8").read()
    assert "possibly missing from the" in src and "are NOT in the Open Targets" not in src


@pytest.mark.live
def test_live_e_4_dravet_approvals_come_from_the_agencies(capsys):
    code = cli.main(["--json", "therapy", "Dravet syndrome"])
    env = json.loads(capsys.readouterr().out)
    assert code == 0
    rows = {r["drug"]: r for r in env["result"]["disease"]["drugs"]["rows"]}
    for drug in ("CANNABIDIOL", "FENFLURAMINE"):
        chk = rows[drug]["agency_check"]
        assert chk["FDA"]["status"] == "approved_for_disease" and chk["EMA"]["status"] == "approved_for_disease"
        assert "stage_warning" not in rows[drug]


def test_e_4_a_label_that_names_the_disease_only_to_exclude_it_is_not_an_approval():
    from zebra.sources import regulators

    assert regulators.names_disease("EPIDIOLEX is indicated for the treatment of seizures associated with "
                                    "Lennox-Gastaut syndrome (LGS), Dravet syndrome (DS).", ["Dravet syndrome"])
    assert regulators.names_disease("X is indicated for focal seizures. Limitations of Use: X is not indicated "
                                    "for Dravet syndrome.", ["Dravet syndrome"]) is None
    assert regulators.names_disease("Safety and effectiveness in Dravet syndrome have not been established.",
                                    ["Dravet syndrome"]) is None
    assert regulators.names_disease("indicated for Duchenne muscular dystrophy (DMD).", ["Duchenne muscular dystrophy"])


# ---------------------------------------------------------------- adversarial review (W4 round), agency checks

@pytest.mark.parametrize("text", [
    "Azithromycin should not be used in patients with pneumonia who are judged to be inappropriate for oral therapy "
    "because of moderate to severe illness or risk factors such as any of the following: • patients with cystic "
    "fibrosis, • patients with nosocomially acquired infections.",
    "X is contraindicated in patients with cystic fibrosis.",
    "X was not effective in cystic fibrosis.",
    "X is indicated for bronchiectasis but not for cystic fibrosis.",
    "Uses temporary relief from symptoms associated with cystic fibrosis.",
    "OJJAARA is indicated for myelofibrosis, including secondary MF [post‑cystic fibrosis (CF)], in adults.",
])
def test_rev_p0_1_a_label_that_does_not_indicate_the_disease_is_not_an_approval(text):
    from zebra.sources import regulators

    assert regulators.names_disease(text, ["cystic fibrosis"]) is None


def test_rev_p0_1_real_indications_still_count():
    from zebra.sources import regulators

    assert regulators.names_disease("X is indicated as adjunctive therapy in patients whose seizures have not been "
                                    "controlled, including Dravet syndrome.", ["Dravet syndrome"])
    assert regulators.names_disease("Jakavi is indicated for adults with polycythaemia vera who are resistant to "
                                    "hydroxyurea.", ["polycythemia vera"])
    assert regulators.names_disease("indicated for late-onset Pompe’s disease.", ["Pompe disease"])


def test_rev_p0_2_p0_3_a_label_only_counts_for_its_own_active_ingredient():
    from zebra.sources import regulators

    cayston = {"openfda": {"generic_name": ["AZTREONAM"], "substance_name": ["AZTREONAM", "SODIUM CHLORIDE"]}}
    assert not regulators._active_here("sodium chloride", cayston)
    evrysdi = {"spl_product_data_elements": ["EVRYSDI Risdiplam RISDIPLAM MANNITOL ISOMALT"]}
    assert regulators._active_here("risdiplam", evrysdi) and not regulators._active_here("mannitol", evrysdi)


def test_rev_p1_ema_matches_the_same_moiety_not_a_word_inside_another():
    from zebra.sources import regulators

    assert "phenylbutyrate" not in regulators._components("glycerol phenylbutyrate")
    assert regulators.moiety("LITHIUM CARBONATE") == "lithium carbonate"
    assert regulators.moiety("LITHIUM CARBONATE") not in regulators._components("sevelamer carbonate")
    assert regulators.moiety("SODIUM BENZOATE") == "sodium benzoate"
    assert regulators._components("elexacaftor / tezacaftor / ivacaftor") == ["elexacaftor", "tezacaftor", "ivacaftor"]


def test_rev_p1_ema_date_is_the_date_of_the_status(agency_network, monkeypatch):
    from zebra.sources import regulators

    def fake(url, source, params=None, **kw):
        body = {"meta": {"total_records": 1, "timestamp": "x"}, "data": [{
            "category": "Human", "name_of_medicine": "Translarna", "international_non_proprietary_name_common_name": "ataluren",
            "active_substance": "ataluren", "medicine_status": "Expired", "therapeutic_indication": "Duchenne muscular dystrophy",
            "marketing_authorisation_date": "31/07/2014",
            "withdrawal_expiry_revocation_lapse_of_marketing_authorisation_date": "28/03/2025",
            "medicine_url": "u"}]}
        if url == regulators.EMA_MEDICINES:
            return Response(url, 200, json.dumps(body), "t", False)
        return Response(url, 404, '{"error":{"code":"NOT_FOUND"}}', "t", False)

    monkeypatch.setattr(regulators, "get_json", fake)
    r = regulators.check("ATALUREN", ["Duchenne muscular dystrophy"]).result
    assert r["EMA"]["status"] == "withdrawn" and r["EMA"]["date"] == "2025-03-28"


def test_rev_p2_an_unexpected_openfda_404_is_an_error_not_no_label(agency_network, monkeypatch):
    from zebra.sources import regulators

    monkeypatch.setattr(regulators, "get_json", lambda url, source, params=None, **kw:
                        Response(url, 404, "<html>not found</html>", "t", False))
    out = regulators.check("CANNABIDIOL", ["Dravet syndrome"])
    assert out.result["FDA"]["status"] == "unavailable"
    assert "could not be checked" in out.result["headline"]


def test_rev_p1_one_agency_check_per_moiety(monkeypatch):
    from zebra.sources import regulators

    calls = []
    monkeypatch.setattr(regulators, "check", lambda d, t, with_fda_date=True: calls.append(d) or Outcome({
        "FDA": {"status": "approved_for_disease"}, "EMA": {"status": "no_record"}, "approved_in": ["FDA"],
        "refused_or_withdrawn_in": [], "headline": "h"}))
    rows = [{"drug": "GIVINOSTAT", "stage": "APPROVAL", "regulatory": {"approved_in": []}},
            {"drug": "GIVINOSTAT HYDROCHLORIDE", "stage": "APPROVAL", "regulatory": {"approved_in": []}}]
    therapy_cmd._agency_checks(rows, ["Duchenne muscular dystrophy"], Outcome({}))
    assert len(calls) == 1 and all("agency_check" in r for r in rows)


def test_rev_p1_alfa_does_not_hide_a_gap():
    rows = [{"drug": "AGALSIDASE ALFA"}]
    assert therapy_cmd._designations_not_in_drug_list([{"substance": "cipaglucosidase alfa"}], rows)
    generic = therapy_cmd._designations_not_in_drug_list([{"substance": "adeno-associated viral vector"}], rows)
    assert generic and "generically" in generic[0]["why"]


def test_rev_p1_pompe_designations_match_the_possessive_wording():
    assert orphan._fold("glycogen storage disease type II (Pompe's disease)").find("pompe disease") >= 0
    assert "glycogen storage disease ii" in orphan._fold("Glycogen storage disease type II")


def test_rev_p2_preapproval_is_under_review_not_no_stage():
    reg = ot.regulatory_status([{"origin": "REGULATORY_AGENCY", "source": "EMA Human Drugs", "clinicalStage": "PREAPPROVAL"}])
    assert reg["under_review_in"] == ["EMA (European Union)"] and "under review, not authorised" in reg["headline"]
    assert "no stage" not in reg["headline"]


def test_rev_advisor_ema_bare_noun_phrase_indications_still_count():
    """The EMA export often writes "Treatment of X." with no "indicated": the FDA-only anchor must not apply."""
    from zebra.sources import regulators

    assert regulators.names_disease("Treatment of Duchenne muscular dystrophy.", ["Duchenne muscular dystrophy"],
                                    anchor=False)
    assert regulators.names_disease("Treatment of Duchenne muscular dystrophy.", ["Duchenne muscular dystrophy"]) is None
    assert regulators.names_disease("Treatment of seizures associated with Dravet syndrome as an add-on therapy for "
                                    "patients 2 years of age and older.Treatment of seizures associated with "
                                    "Lennox-Gastaut syndrome.", ["Lennox-Gastaut syndrome"], anchor=False)
    # negation still applies without the anchor
    assert regulators.names_disease("Not for the treatment of Dravet syndrome.", ["Dravet syndrome"], anchor=False) is None


def test_rev_advisor_ema_record_with_a_bare_indication_is_approved(agency_network, monkeypatch):
    from zebra.sources import regulators

    def fake(url, source, params=None, **kw):
        body = {"meta": {"total_records": 1, "timestamp": "x"}, "data": [{
            "category": "Human", "name_of_medicine": "Fintepla", "international_non_proprietary_name_common_name": "fenfluramine",
            "active_substance": "fenfluramine hydrochloride", "medicine_status": "Authorised",
            "therapeutic_indication": "Treatment of seizures associated with Dravet syndrome as an add-on therapy.",
            "marketing_authorisation_date": "18/12/2020", "medicine_url": "u"}]}
        if url == regulators.EMA_MEDICINES:
            return Response(url, 200, json.dumps(body), "t", False)
        return Response(url, 404, '{"error":{"code":"NOT_FOUND"}}', "t", False)

    monkeypatch.setattr(regulators, "get_json", fake)
    r = regulators.check("FENFLURAMINE HYDROCHLORIDE", ["Dravet syndrome"]).result
    assert r["EMA"]["status"] == "approved_for_disease" and r["EMA"]["date"] == "2020-12-18"
