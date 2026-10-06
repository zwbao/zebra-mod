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


# ------------------------------------- F39: GeneReviews text, plain language, support pointers

def test_F39_genereviews_sections_are_parsed_from_the_published_summary():
    body = _json("genereviews/f39_epmc_book_text.json")
    by_id = {r["bookid"]: r for r in body["resultList"]["result"]}
    secs = genereviews.parse_sections(by_id["NBK1119"]["abstractText"])
    assert list(secs) == ["Clinical characteristics", "Diagnosis/testing", "Management", "Genetic counseling"]
    assert "Duchenne muscular dystrophy" in secs["Clinical characteristics"]
    assert "Corticosteroid therapy" in secs["Management"]
    assert "<" not in secs["Management"]  # tags stripped
    assert by_id["NBK1119"]["bookOrReportDetails"]["comprisingTitle"].startswith("GeneReviews")


def test_F39_chapter_text_keyed_by_nbk(monkeypatch):
    from zebra.http import Response

    body = _json("genereviews/f39_epmc_book_text.json")
    seen = {}

    def fake_get(url, source, **kw):
        seen["query"] = kw["params"]["query"]
        return Response(url, 200, json.dumps(body), "2026-10-06T00:00:00+00:00", True)

    monkeypatch.setattr(genereviews, "get_json", fake_get)
    out = genereviews.chapter_text(["NBK1119", "NBK1318"])
    assert seen["query"] == "BOOK_ID:NBK1119 OR BOOK_ID:NBK1318"
    assert set(out.result) == {"NBK1119", "NBK1318"}
    assert out.result["NBK1318"]["pmid"] == "20301494"
    assert "Clinical characteristics" in out.result["NBK1318"]["sections"]
    assert "the full chapter is at https://www.ncbi.nlm.nih.gov/books/NBK1318/" in out.result["NBK1318"]["source"]


def test_F39_a_chapter_with_no_book_record_is_warned_about(monkeypatch):
    from zebra.http import Response

    body = _json("genereviews/f39_epmc_book_text.json")

    def fake_get(url, source, **kw):
        return Response(url, 200, json.dumps(body), "2026-10-06T00:00:00+00:00", True)

    monkeypatch.setattr(genereviews, "get_json", fake_get)
    out = genereviews.chapter_text(["NBK1119", "NBK99999"])
    assert "NBK99999" not in out.result
    assert any("NBK99999" in w and "chapter text not retrieved" in w for w in out.warnings)


def test_F39_medlineplus_parses_a_plain_language_summary():
    from zebra.sources import medlineplus

    data = _json("medlineplus/f39_osteogenesis_imperfecta.json")
    r = medlineplus.parse_condition(data, "osteogenesis-imperfecta", ["OMIM:166200"])
    assert r["name"] == "Osteogenesis imperfecta"
    assert r["url"].startswith("https://medlineplus.gov/genetics/condition/")
    body = next(b for b in r["text"] if b["role"] == "description")
    assert "bones that break (fracture) easily" in body["text"] and "<" not in body["text"]
    assert "COL1A1" in r["genes"]
    assert r["register"].startswith("MedlinePlus Genetics")


def test_F39_medlineplus_says_when_a_match_rests_on_the_name_alone():
    from zebra.sources import medlineplus

    data = _json("medlineplus/f39_osteogenesis_imperfecta.json")
    confirmed = medlineplus.parse_condition(data, "osteogenesis-imperfecta", data_omim(data))
    assert confirmed["id_confirmed"] is True and confirmed["omim_overlap"]
    name_only = medlineplus.parse_condition(data, "osteogenesis-imperfecta", ["OMIM:999999"])
    assert name_only["id_confirmed"] is False and name_only["omim_overlap"] == []


def data_omim(data):
    return ["OMIM:" + k["db-key"]["key"] for k in data["db-key-list"] if k["db-key"]["db"] == "OMIM"][:1]


def test_F39_medlineplus_slug_building_and_candidate_order():
    from zebra.sources import medlineplus

    assert medlineplus.slug("Duchenne and Becker muscular dystrophy") == "duchenne-and-becker-muscular-dystrophy"
    assert medlineplus.slug("Gaucher's disease") == "gauchers-disease"
    cands = medlineplus.slug_candidates(["Osteogenesis imperfecta", "OI", "", "Osteogenesis imperfecta",
                                         "Brittle bone disease"])
    assert cands == ["osteogenesis-imperfecta", "brittle-bone-disease"]  # "OI" is too short, duplicates dropped


def test_F39_medlineplus_miss_is_a_404_not_an_error(monkeypatch):
    from zebra.http import Response
    from zebra.sources import medlineplus

    def fake_get(url, source, **kw):
        assert kw["ok_statuses"] == (200, 404)
        return Response(url, 404, "<!DOCTYPE html><html>not found", "2026-10-06T00:00:00+00:00", False)

    monkeypatch.setattr(medlineplus, "get_json", fake_get)
    out = medlineplus.condition(["Dravet syndrome"])
    assert out.result is None
    assert any("no MedlinePlus Genetics page found" in w and "dravet-syndrome" in w for w in out.warnings)


def test_F39_support_pointers_name_their_source_and_say_they_are_not_retrieved():
    sup = D.support_pointers("ORPHA:98896", {"url": "https://www.orpha.net/en/disease/detail/98896"},
                             {"GARD": [{"id": "GARD:6291", "relation": "E"}]},
                             {"url": "https://medlineplus.gov/genetics/condition/x"})
    whats = {p["what"]: p for p in sup["pointers"]}
    assert "patient organisations and expert centres" in whats
    assert whats["patient organisations and expert centres"]["retrieved"] is False
    gard = whats["patient organisations and advocacy groups"]
    assert gard["url"] == "https://rarediseases.info.nih.gov/diseases/6291/index"
    assert "serves no patient-organisation or expert-centre dataset" in sup["note"]
    assert "rd-cross-referencing" in sup["note"]


# ---------------------------------------------------------------- F12: a deterministic 404 is a miss

def test_F12_orphanet_name_search_does_not_use_the_bare_json_accept(monkeypatch):
    from zebra.http import Response

    seen = {}

    def fake_get(url, source, **kw):
        seen.update(kw)
        return Response(url, 404, "<html>404</html>", "2026-10-06T00:00:00+00:00", False)

    monkeypatch.setattr(orphanet, "get_json", fake_get)
    orphanet._BY_NAME_CACHE.clear()
    out = orphanet.by_name("Hyperphenylalaninemia/PKU")
    # the HTML-where-JSON-was-expected retry only fires for accept == "application/json"
    assert seen["accept"] == "application/json, */*"
    assert seen["ok_statuses"] == (200, 404)
    assert out.result is None


def test_F12_slash_is_replaced_before_quoting_and_said_so(monkeypatch):
    from zebra.http import Response

    urls = []

    def fake_get(url, source, **kw):
        urls.append(url)
        return Response(url, 404, "<html>404</html>", "2026-10-06T00:00:00+00:00", False)

    monkeypatch.setattr(orphanet, "get_json", fake_get)
    orphanet._BY_NAME_CACHE.clear()
    out = orphanet.by_name("Hyperphenylalaninemia/PKU")
    assert "%2F" not in urls[0] and "Hyperphenylalaninemia%20PKU" in urls[0]
    assert any("rejects '/' in a name" in w for w in out.warnings)


def test_F12_by_name_is_memoized_so_one_run_makes_one_request(monkeypatch):
    from zebra.http import Response

    calls = []

    def fake_get(url, source, **kw):
        calls.append(url)
        return Response(url, 404, "<html>404</html>", "2026-10-06T00:00:00+00:00", False)

    monkeypatch.setattr(orphanet, "get_json", fake_get)
    orphanet._BY_NAME_CACHE.clear()
    first = orphanet.by_name("Some Name")
    second = orphanet.by_name("Some Name")
    assert len(calls) == 1
    assert second.result == first.result and second.warnings == []  # not repeated per caller


@pytest.mark.live
def test_live_F39_genereviews_text_and_support_on_the_dmd_card(capsys):
    code = cli.main(["--json", "disease", "ORPHA:98896"])
    env = json.loads(capsys.readouterr().out)
    assert code == 0
    r = env["result"]
    ch = r["genereviews"][0]
    assert ch["nbk"] == "NBK1119" and ch["pmid"] == "20301298"
    assert "Clinical characteristics" in ch["sections"] and "Management" in ch["sections"]
    assert "Corticosteroid" in ch["sections"]["Management"]
    pointers = {p["source"] for p in r["support"]["pointers"]}
    assert any("Orphanet disease page" in s for s in pointers)
    assert any("GARD" in s for s in pointers)
    assert "serves no patient-organisation or expert-centre dataset" in r["support"]["note"]


@pytest.mark.live
def test_live_F39_medlineplus_plain_language_on_osteogenesis_imperfecta(capsys):
    code = cli.main(["--json", "disease", "ORPHA:666"])
    r = json.loads(capsys.readouterr().out)["result"]
    assert code == 0
    mp = r["plain_language"]
    assert mp["name"] == "Osteogenesis imperfecta"
    assert mp["url"].startswith("https://medlineplus.gov/genetics/condition/osteogenesis-imperfecta")
    assert any("bones" in b["text"] for b in mp["text"])
    assert mp["register"].startswith("MedlinePlus Genetics")


@pytest.mark.live
def test_live_F12_a_disease_name_with_a_slash_is_answered_fast_and_not_adopted():
    """Before: an HTML 404 retried 4 times with 10.5 s of backoff, twice per run (62 s total).

    With `/` replaced by a space the endpoint answers normally, so there is no
    retry storm at all — and the answer is Orphadata's *closest* name, which
    `name_matches` must still reject so the resolver does not adopt it.
    """
    import time

    from zebra.sources import orphanet as O

    O._BY_NAME_CACHE.clear()
    t0 = time.monotonic()
    out = O.by_name("Hyperphenylalaninemia/PKU")
    elapsed = time.monotonic() - t0
    assert elapsed < 20  # was 4 attempts plus 10.5 s of backoff, twice per run
    assert any("rejects '/' in a name" in w for w in out.warnings)
    assert O.name_matches(out.result, "Hyperphenylalaninemia/PKU") is False
    # and the memoized second call costs nothing
    t1 = time.monotonic()
    again = O.by_name("Hyperphenylalaninemia/PKU")
    assert time.monotonic() - t1 < 0.05 and again.result == out.result


@pytest.mark.live
def test_live_F12_the_command_calls_it_a_miss_not_an_unavailable_source(capsys):
    code = cli.main(["--json", "disease", "Hyperphenylalaninemia/PKU"])
    env = json.loads(capsys.readouterr().out)
    assert code == 0
    assert env["result"]["status"] in ("ambiguous", "not_found")  # never "unavailable"
    assert not any("Orphanet name search unavailable" in w for w in env["warnings"])


# ---------------------------------------------------------------- E-2: a generic Chinese word never resolves silently

@pytest.mark.parametrize("query,near", [("糖尿病", "ORPHA:511"), ("白内障", "ORPHA:98989"),
                                        ("智力障碍", "ORPHA:87277"), ("心肌病", "ORPHA:167848")])
def test_e_2_a_generic_chinese_word_is_ambiguous_not_a_card(query, near, monkeypatch, capsys):
    called = []
    monkeypatch.setattr(orphanet, "disorder", lambda *a, **k: called.append(a) or Outcome(None))
    code = cli.main(["--json", "disease", query])
    env = json.loads(capsys.readouterr().out)
    r = env["result"]
    assert code == 0 and r["status"] == "ambiguous", r
    assert near in {c["id"] for c in r["candidates"]}
    assert called == []  # no card was fetched for a near miss
    assert any("not an exact Chinese disease name" in w for w in env["warnings"])


def test_a_near_exact_chinese_name_still_resolves(monkeypatch):
    seen = []

    def from_orpha(self, orpha_id, how="given"):
        seen.append((orpha_id, how))
        self.orpha = orpha_id
        return True

    monkeypatch.setattr(D.Resolver, "from_orpha", from_orpha)
    r, status, _ = D.resolve("杜氏肌营养不良")
    assert status == "resolved" and seen[0][0] == "ORPHA:98896"  # bigram 0.923 >= 0.85


# ---------------------------------------------------------------- CP1-7: the Chinese entry is as strong as the English

def test_cp1_7_pompe_in_chinese_resolves_to_pompe_disease_not_the_gsd_group(monkeypatch):
    seen = []

    def from_orpha(self, orpha_id, how="given"):
        seen.append((orpha_id, how))
        self.orpha = orpha_id
        return True

    monkeypatch.setattr(D.Resolver, "from_orpha", from_orpha)
    r, status, _ = D.resolve("庞贝病")
    assert status == "resolved" and seen[0][0] == "ORPHA:365"
    assert "names ORPHA:365 itself" in seen[0][1]


def test_cp1_7_a_list_name_with_no_orphacode_follows_the_english_name(monkeypatch):
    """脊髓性肌萎缩症 (list #110, no ORPHAcode) used to stop at ORPHA:70 by bigram 0.769 and lost GeneReviews."""
    tried = []

    def by_name(name, lang="en"):
        tried.append(name)
        return Outcome({"id": "ORPHA:83330", "name": "x"})

    monkeypatch.setattr(orphanet, "by_name", by_name)
    monkeypatch.setattr(orphanet, "name_matches", lambda row, name: False)
    monkeypatch.setattr(monarch, "search", lambda *a, **k: Outcome({"hits": [{"id": "MONDO:0001516",
                                                                             "name": "spinal muscular atrophy"}]}))
    monkeypatch.setattr(D.Resolver, "from_mondo", lambda self, mid: setattr(self, "name", "spinal muscular atrophy") or True)
    r, status, _ = D.resolve("脊髓性肌萎缩症")
    assert status == "resolved" and "Spinal Muscular Atrophy" in tried


def test_cp1_7_genereviews_is_found_by_disease_gene_when_name_and_omim_miss(monkeypatch):
    calls = []

    def chapters(omim_ids=(), genes=(), name=None, limit=5, with_text=True):
        calls.append({"omim": list(omim_ids), "genes": list(genes), "name": name})
        if genes:
            return Outcome({"chapters": [{"nbk": "NBK1352", "title": "Spinal Muscular Atrophy", "url": "u",
                                          "matched_by": ["gene SMN1"]}]})
        return Outcome({"chapters": []})

    from zebra.sources import hpo, medlineplus

    monkeypatch.setattr(genereviews, "chapters", chapters)
    monkeypatch.setattr(orphanet, "disorder", lambda c, lang="en": Outcome(None))
    monkeypatch.setattr(orphanet, "epidemiology", lambda c: Outcome([]))
    monkeypatch.setattr(orphanet, "natural_history", lambda c: Outcome({}))
    monkeypatch.setattr(orphanet, "genes", lambda c: Outcome([{"symbol": "SMN1", "hgnc": "HGNC:11117",
                                                               "association": "Disease-causing germline mutation(s) in",
                                                               "status": "Assessed", "pmids": []}]))
    monkeypatch.setattr(hpo, "disease_annotations", lambda i: Outcome({"phenotypes": []}))
    monkeypatch.setattr(medlineplus, "condition", lambda names, omim=(): Outcome(None))
    r = D.Resolver()
    r.orpha, r.orpha_row, r.name = "ORPHA:70", {"id": "ORPHA:70", "name": "Proximal spinal muscular atrophy"}, \
        "Proximal spinal muscular atrophy"
    card = D.build_card(r)
    assert [c["nbk"] for c in card["genereviews"]] == ["NBK1352"]
    assert calls[-1]["genes"] == ["SMN1"]


# ---------------------------------------------------------------- E-1 on the disease card

def test_e_1_the_card_does_not_match_the_list_through_an_acronym_synonym(monkeypatch):
    """catastrophic antiphospholipid syndrome has the synonym CAPS, which is also cryopyrin-associated periodic syndrome."""
    from zebra.sources import hpo, medlineplus

    monkeypatch.setattr(genereviews, "chapters", lambda **k: Outcome({"chapters": []}))
    monkeypatch.setattr(orphanet, "disorder", lambda c, lang="en": Outcome(None))
    monkeypatch.setattr(orphanet, "epidemiology", lambda c: Outcome([]))
    monkeypatch.setattr(orphanet, "natural_history", lambda c: Outcome({}))
    monkeypatch.setattr(orphanet, "genes", lambda c: Outcome([]))
    monkeypatch.setattr(hpo, "disease_annotations", lambda i: Outcome({"phenotypes": []}))
    monkeypatch.setattr(medlineplus, "condition", lambda names, omim=(): Outcome(None))
    r = D.Resolver()
    r.orpha = "ORPHA:464343"
    r.orpha_row = {"id": "ORPHA:464343", "name": "Catastrophic antiphospholipid syndrome",
                   "synonyms": ["CAPS", "Catastrophic APS", "caps"]}
    r.name = "Catastrophic antiphospholipid syndrome"
    card = D.build_card(r)
    assert card["china_rare_list"]["on_list"] is False
    assert all(m.get("matched_on") not in ("CAPS", "caps") for m in card["china_rare_list"]["matches"])


def test_e_1_the_card_matches_by_orphacode_and_says_qualified(monkeypatch):
    from zebra.sources import hpo, medlineplus

    monkeypatch.setattr(genereviews, "chapters", lambda **k: Outcome({"chapters": []}))
    monkeypatch.setattr(orphanet, "disorder", lambda c, lang="en": Outcome(None))
    monkeypatch.setattr(orphanet, "epidemiology", lambda c: Outcome([]))
    monkeypatch.setattr(orphanet, "natural_history", lambda c: Outcome({}))
    monkeypatch.setattr(orphanet, "genes", lambda c: Outcome([]))
    monkeypatch.setattr(hpo, "disease_annotations", lambda i: Outcome({"phenotypes": []}))
    monkeypatch.setattr(medlineplus, "condition", lambda names, omim=(): Outcome(None))
    r = D.Resolver()
    r.orpha, r.orpha_row, r.name = "ORPHA:79201", {"id": "ORPHA:79201", "name": "Glycogen storage disease"}, \
        "Glycogen storage disease"
    card = D.build_card(r)
    cn = card["china_rare_list"]
    assert cn["status"] == "qualified" and cn["on_list"] is False and cn["matches"][0]["no"] == 35
    assert any("covers glycogen storage disease types I and II only" in w for w in r.warnings)
    assert "only a subtype is listed" in D.render(dict(card, status="resolved"))


# ---------------------------------------------------------------- E-7: a 404 slug is not evidence

def test_e_7_medlineplus_records_only_the_page_that_answered(monkeypatch):
    from zebra.http import Response
    from zebra.sources import medlineplus

    good = json.dumps({"name": "Spinal muscular atrophy", "ghr_page": "https://medlineplus.gov/genetics/condition/spinal-muscular-atrophy",
                       "text-list": [], "db-key-list": []})

    def fake_get(url, source, **kw):
        if url.endswith("/spinal-muscular-atrophy.json"):
            return Response(url, 200, good, "2026-10-06T00:00:00+00:00", False)
        return Response(url, 404, "<!DOCTYPE html><html>not found", "2026-10-06T00:00:00+00:00", False)

    monkeypatch.setattr(medlineplus, "get_json", fake_get)
    out = medlineplus.condition(["Proximal spinal muscular atrophy", "SMA type 1", "Spinal muscular atrophy"])
    assert out.result["slug"] == "spinal-muscular-atrophy"
    assert [s["record"] for s in out.sources] == ["spinal-muscular-atrophy"]
    assert out.result["tried"] == ["proximal-spinal-muscular-atrophy", "sma-type-1", "spinal-muscular-atrophy"]
    miss = medlineplus.condition(["Dravet syndrome", "SMEI"])
    assert miss.result is None and miss.sources == []
    assert any("dravet-syndrome" in w for w in miss.warnings)


@pytest.mark.live
def test_live_e_2_cp1_7_chinese_entries(capsys):
    code = cli.main(["--json", "disease", "糖尿病"])
    assert json.loads(capsys.readouterr().out)["result"]["status"] == "ambiguous"
    code = cli.main(["--json", "disease", "脊髓性肌萎缩症"])
    r = json.loads(capsys.readouterr().out)["result"]
    assert r["status"] == "resolved" and "NBK1352" in {c["nbk"] for c in r["genereviews"]}
    assert r["plain_language"] and r["china_rare_list"]["status"] == "on_list"
    code = cli.main(["--json", "disease", "庞贝病"])
    r = json.loads(capsys.readouterr().out)["result"]
    assert r["ids"]["ORPHA"] == "ORPHA:365" and r["china_rare_list"]["matches"][0]["no"] == 35


@pytest.mark.parametrize("raw", ["0.0", "0", 0, 0.0, "0.00", None])
def test_e_13_orphanet_mean_value_zero_is_no_value(raw):
    data = {"data": {"results": {"Prevalence": [{"PrevalenceType": "Point prevalence", "ValMoy": raw}]}}}
    rows = orphanet.parse_epidemiology(data)
    assert rows and rows[0]["value"] is None
    data["data"]["results"]["Prevalence"][0]["ValMoy"] = "3.3"
    assert orphanet.parse_epidemiology(data)[0]["value"] == "3.3"


# ---------------------------------------------------------------- adversarial review (W4 round), Chinese resolution

@pytest.mark.parametrize("query,wrong", [("甲基丙二酸血症", "ORPHA:289504"), ("神经纤维瘤病", "ORPHA:252183"),
                                         ("地中海贫血", "ORPHA:846"), ("酪氨酸血症Ⅲ型", "ORPHA:28378"),
                                         ("脊髓小脑性共济失调3型", "ORPHA:276183"), ("黏多糖贮积症Ⅲ型", "ORPHA:217085"),
                                         ("神经元蜡样脂褐质沉积症", "ORPHA:228329")])
def test_rev_p0_3_a_near_chinese_name_never_resolves_to_a_different_disease(query, wrong, monkeypatch):
    seen = []

    def from_orpha(self, orpha_id, how="given"):
        seen.append(orpha_id)
        self.orpha = orpha_id
        return True

    monkeypatch.setattr(D.Resolver, "from_orpha", from_orpha)
    r = D.Resolver()
    D.resolve_chinese(r, query)
    assert wrong not in seen, (query, seen)


def test_rev_p0_3_zh_near_exact_rule():
    assert D.zh_near_exact("杜氏肌营养不良", "杜氏肌营养不良症")
    assert D.zh_near_exact("成骨不全症", "成骨不全")
    assert not D.zh_near_exact("神经纤维瘤病", "神经纤维瘤")
    assert not D.zh_near_exact("酪氨酸血症Ⅲ型", "酪氨酸血症Ⅱ型")


def test_rev_p0_3_a_misspelt_list_name_is_followed_by_its_corrected_spelling(monkeypatch):
    tried = []
    monkeypatch.setattr(orphanet, "by_name", lambda name, lang="en": tried.append(name) or Outcome(None))
    monkeypatch.setattr(orphanet, "name_matches", lambda row, name: False)
    monkeypatch.setattr(monarch, "search", lambda *a, **k: Outcome({"hits": []}))
    D.resolve("甲基丙二酸血症")
    assert tried and tried[0] == "Methylmalonic Acidemia"


def test_rev_p1_genereviews_outage_is_not_written_as_no_chapter(monkeypatch):
    from zebra.http import SourceError
    from zebra.sources import hpo, medlineplus

    def boom(**k):
        raise SourceError("GeneReviews (NCBI FTP)", "u", 503, "down")

    monkeypatch.setattr(genereviews, "chapters", boom)
    monkeypatch.setattr(orphanet, "disorder", lambda c, lang="en": Outcome(None))
    monkeypatch.setattr(orphanet, "epidemiology", lambda c: Outcome([]))
    monkeypatch.setattr(orphanet, "natural_history", lambda c: Outcome({}))
    monkeypatch.setattr(orphanet, "genes", lambda c: Outcome([]))
    monkeypatch.setattr(hpo, "disease_annotations", lambda i: Outcome({"phenotypes": []}))
    monkeypatch.setattr(medlineplus, "condition", lambda names, omim=(): Outcome(None))
    r = D.Resolver()
    r.orpha, r.orpha_row, r.name = "ORPHA:33069", {"id": "ORPHA:33069", "name": "Dravet syndrome"}, "Dravet syndrome"
    card = D.build_card(r)
    text = D.render(dict(card, status="resolved"))
    assert "GeneReviews: not checked — the source was unavailable" in text
    assert "no chapter found" not in text
