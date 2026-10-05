import json
from pathlib import Path

import pytest

from zebra import cli
from zebra.commands import disease as D
from zebra.core import Outcome
from zebra.sources import genereviews, monarch, ols, orphanet

FIX = Path(__file__).parent / "fixtures"


def _json(rel):
    return json.loads((FIX / rel).read_text("utf-8"))


def _text(rel):
    return (FIX / rel).read_text("utf-8")


# ---------------------------------------------------------------- Orphanet parsers


def test_orphanet_disorder_refs_and_relations():
    row = orphanet.parse_disorder(_json("orphanet/cross_33069_en.json")["data"]["results"])
    assert row["id"] == "ORPHA:33069" and row["name"] == "Dravet syndrome"
    assert row["definition"].startswith("A rare, genetic, developmental and epileptic encephalopathy")
    assert orphanet.exact_refs(row, "OMIM") == ["607208"]  # 612164 / 615744 are BTNT, not exact
    rel = {(r["source"], r["id"]): r["relation"] for r in row["refs"]}
    assert rel[("OMIM", "612164")] == "BTNT" and rel[("ICD-10", "G40.4")] == "NTBT" and rel[("ICD-11", "8A61.11")] == "E"
    assert orphanet.exact_refs(row, "MONDO") == ["0011794"]
    refs = orphanet.curie_refs(row)
    assert {"id": "ICD10:G40.4", "relation": "NTBT"} in refs["ICD-10"]


def test_orphanet_chinese_name():
    zh = orphanet.parse_disorder(_json("orphanet/cross_33069_zh.json")["data"]["results"])
    assert zh["name"] == "Dravet综合征" and "婴儿严重肌阵挛型癫痫" in zh["synonyms"]
    assert zh["dataset_date"].startswith("2020-06-01")


def test_orphanet_epi_natural_history_genes():
    prev = orphanet.parse_epidemiology(_json("orphanet/epidemiology_33069.json"))
    assert prev[0] == {"type": "Prevalence at birth", "class": "1-9 / 100 000", "value": "3.3",
                       "qualification": "Value and class", "region": "Europe", "validation": "Validated",
                       "source": "22719002[PMID]_25772213[PMID]_ 31302675[PMID]"}
    nh = orphanet.parse_natural_history(_json("orphanet/natural_history_33069.json"))
    assert nh == {"onset": ["Infancy", "Neonatal"], "inheritance": ["Autosomal dominant"]}
    genes = orphanet.parse_genes(_json("orphanet/genes_33069.json"))
    assert genes[0]["symbol"] == "SCN1A" and genes[0]["hgnc"] == "HGNC:10585"
    assert genes[0]["association"] == "Disease-causing germline mutation(s) in" and genes[0]["status"] == "Assessed"
    assert genes[0]["pmids"][0] == "20301494"


def test_orphanet_name_endpoint_is_checked():
    # the name endpoint returns its single closest match: "Gaucher's Disease" comes back as Alexander disease
    row = orphanet.parse_disorder(_json("orphanet/name_gauchers_disease.json")["data"]["results"])
    assert row["name"] == "Alexander disease"
    assert orphanet.name_matches(row, "Gaucher’s Disease") is False
    dravet = orphanet.parse_disorder(_json("orphanet/cross_33069_en.json")["data"]["results"])
    assert orphanet.name_matches(dravet, "dravet syndrome")
    assert orphanet.name_matches(dravet, "Severe myoclonic epilepsy of infancy")  # synonym


def test_orphanet_codes():
    assert orphanet.orpha_code("ORPHA:33069") == "33069" == orphanet.orpha_code("Orphanet_33069") == orphanet.orpha_code(33069)
    with pytest.raises(ValueError):
        orphanet.orpha_code("OMIM:607208")


# ---------------------------------------------------------------- Monarch / OLS / GeneReviews parsers


def test_monarch_entity_and_search():
    e = monarch.parse_entity(_json("monarch/entity_MONDO_0100135.json"))
    assert e["name"] == "Dravet syndrome" and e["omim"] == [] and e["orpha"] == []  # Mondo dropped them
    assert "UMLS:C0751122" in e["xref"]
    hits = monarch.parse_search(_json("monarch/search_dravet.json"), 3)
    assert hits[0]["id"] == "MONDO:0100135"


def test_monarch_mappings_and_associations():
    rows = monarch.parse_mappings(_json("monarch/mappings_dravet.json"))
    pairs = {(r["object"], r["subject"]) for r in rows}
    assert ("OMIM:607208", "MONDO:0100079") in pairs and ("ORPHA:33069", "MONDO:0011794") in pairs
    assocs = monarch.parse_associations(_json("monarch/assoc_HGNC_10585.json"), "object")
    assert assocs[0]["id"].startswith("MONDO:") and assocs[0]["source"] in ("omim", "orphanet", "clingen")
    assert monarch.normalize_curie("Orphanet:33069") == "ORPHA:33069" and monarch.to_monarch_curie("ORPHA:1") == "Orphanet:1"


def test_ols_obsolete_term_and_oxo():
    t = ols.parse_term(_json("ols/term_MONDO_0011794.json"))
    assert t["obsolete"] is True and t["replaced_by"] == "MONDO:0100135"
    cur = ols.parse_term(_json("ols/term_MONDO_0100135.json"))
    assert cur["obsolete"] is False and cur["label"] == "Dravet syndrome"
    assert "DOID:0080422" in cur["exact_xrefs"]
    ox = ols.parse_oxo(_json("ols/oxo_Orphanet_33069.json"), "Orphanet:33069")
    assert ox == [{"id": "MONDO:0011794", "label": "obsolete Dravet syndrome", "scope": "EXACT", "mapping_set": "mondo.sssom.tsv"}]


def test_genereviews_tables_and_esummary():
    omap = genereviews.parse_omim_map(_text("genereviews/NBKid_shortname_OMIM.txt"))
    rows = genereviews.parse_gene_map(_text("genereviews/GRshortname_NBKid_genesymbol_dzname.txt"))
    hits = genereviews.match_chapters(omap, rows, omim_ids=["OMIM:607208"])
    assert list(hits) == ["NBK1318"] and hits["NBK1318"]["title"] == "SCN1A Seizure Disorders"
    by_gene = genereviews.match_chapters(omap, rows, genes=["SCN1A"])
    assert set(by_gene) == {"NBK1318", "NBK1388"}
    by_name = genereviews.match_chapters(omap, rows, name="Gaucher disease")
    assert list(by_name) == ["NBK1269"]
    assert genereviews.match_chapters(omap, rows, name="disease") == {}  # one word must equal the title
    meta = genereviews.parse_esummary(_json("genereviews/esummary_NBK1318.json"))
    assert list(meta) == ["NBK1318"]  # the table record is dropped
    assert "Updated 2022" in meta["NBK1318"]["dates"] and meta["NBK1318"]["authors"].startswith("Miller IO")


# ---------------------------------------------------------------- command pieces


@pytest.mark.parametrize("text,expected", [
    ("ORPHA:33069", ("ORPHA", "ORPHA:33069")), ("orphanet:33069", ("ORPHA", "ORPHA:33069")),
    ("OMIM:607208", ("OMIM", "OMIM:607208")), ("MIM 607208", ("OMIM", "OMIM:607208")),
    ("MONDO:0100135", ("MONDO", "MONDO:0100135")), ("Dravet syndrome", ("name", "Dravet syndrome")),
    ("戈谢病", ("name", "戈谢病")),
])
def test_parse_query(text, expected):
    assert D.parse_query(text) == expected


def test_frequency_weights():
    assert D._freq_weight("Very frequent") > D._freq_weight("Frequent") > D._freq_weight("Occasional")
    assert D._freq_weight("7/7") == 1.0 and D._freq_weight("25%") == 0.25


def test_resolver_from_orpha_offline(monkeypatch):
    row = orphanet.parse_disorder(_json("orphanet/cross_33069_en.json")["data"]["results"])
    monkeypatch.setattr(orphanet, "disorder", lambda code, lang="en": Outcome(row, sources=[{"db": "Orphanet", "record": "33069"}]))
    monkeypatch.setattr(monarch, "mappings", lambda object_ids=(), subject_ids=(), predicate=None: Outcome(
        monarch.parse_mappings(_json("monarch/mappings_dravet.json")), sources=[{"db": "Monarch mappings"}]))
    t = ols.parse_term(_json("ols/term_MONDO_0011794.json"))
    monkeypatch.setattr(ols, "resolve_obsolete", lambda c: Outcome({x: t["replaced_by"] for x in c}, sources=[{"db": "OLS"}]))
    r = D.Resolver()
    assert r.from_orpha("ORPHA:33069")
    assert r.orpha == "ORPHA:33069"
    assert r.omim["OMIM:607208"].startswith("exact") and "BTNT" in r.omim["OMIM:615744"]
    assert list(r.mondo)[0] == "MONDO:0100135"  # the ORPHA route first
    assert "obsolete" in r.mondo["MONDO:0100135"] and r.mondo["MONDO:0100079"].startswith("exactMatch of OMIM:607208")
    assert "MONDO:0011794" not in r.mondo
    assert len(r.out.sources) == 3


# ---------------------------------------------------------------- live


def _card(argv, capsys):
    code = cli.main(["--json", "disease", *argv])
    env = json.loads(capsys.readouterr().out)
    assert code == 0 and env["ok"]
    return env


@pytest.mark.live
@pytest.mark.parametrize("query", ["Dravet syndrome", "ORPHA:33069", "OMIM:607208", "MONDO:0100135"])
def test_live_dravet_cards_are_coherent(query, capsys):
    env = _card([query], capsys)
    c = env["result"]
    assert c["status"] == "resolved"
    assert c["ids"]["ORPHA"] == "ORPHA:33069"
    assert any(o["id"] == "OMIM:607208" for o in c["ids"]["OMIM"])
    assert any(m["id"] == "MONDO:0100135" for m in c["ids"]["MONDO"])
    assert c["name"] == "Dravet syndrome" and c["name_zh"] == "Dravet综合征"
    assert "Autosomal dominant" in c["inheritance"]
    assert any(g["symbol"] == "SCN1A" for g in c["genes"])
    assert c["genereviews"] and c["genereviews"][0]["nbk"] == "NBK1318"
    assert c["china_rare_list"]["on_list"] and c["china_rare_list"]["matches"][0]["no"] == 105
    assert c["hpo_annotations"]["top"]
    dbs = {s["db"] for s in env["sources"]}
    assert {"Orphanet cross-referencing", "Orphanet epidemiology", "HPO annotations"} <= dbs


@pytest.mark.live
def test_live_chinese_name_and_gene_card(capsys):
    c = _card(["戈谢病"], capsys)["result"]
    assert c["ids"]["ORPHA"] == "ORPHA:355" and c["china_rare_list"]["matches"][0]["no"] == 31
    assert any(g["symbol"] == "GBA1" for g in c["genes"])


@pytest.mark.live
def test_live_misspelt_name_lists_candidates(capsys):
    env = _card(["Dravet sindrome"], capsys)
    assert env["result"]["status"] in ("ambiguous", "not_found")
    if env["result"]["status"] == "ambiguous":
        assert any(x["id"] for x in env["result"]["candidates"])


@pytest.mark.live
def test_live_monarch_gene_diseases_scn1a():
    out = monarch.gene_diseases("SCN1A")
    assert out.result["gene"] == "HGNC:10585"
    assert any(d["id"] == "MONDO:0100135" for d in out.result["diseases"])


@pytest.mark.live
def test_live_genereviews_ngly1():
    out = genereviews.chapters(omim_ids=["OMIM:615273"])
    assert out.result["chapters"][0]["nbk"] == "NBK481554"


def test_sources_down_reports_unavailable_not_absent(monkeypatch, capsys):
    from zebra.http import SourceError

    def boom(*a, **k):
        raise SourceError("Orphadata", "https://api.orphadata.com/x", 503, "service unavailable")

    monkeypatch.setattr(orphanet, "disorder", boom)
    monkeypatch.setattr(monarch, "mappings", boom)
    code = cli.main(["--json", "disease", "ORPHA:33069"])
    env = json.loads(capsys.readouterr().out)
    assert code == 0
    assert env["result"]["status"] == "unavailable"
    assert any(w.startswith("Orphanet unavailable") for w in env["warnings"])
