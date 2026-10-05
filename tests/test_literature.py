"""Literature sources (Europe PMC, LitVar2, PubTator3) and `zebra lit`.

Offline tests replay real responses captured from the services (trimmed) in
tests/fixtures/{europepmc,litvar,pubtator}; live tests call the services.
"""

from __future__ import annotations

import json
import os

import pytest

from zebra import cli
from zebra.http import Response, SourceError
from zebra.sources import europepmc, litvar, pubtator

FIX = os.path.join(os.path.dirname(__file__), "fixtures")


def fixture(rel: str) -> str:
    with open(os.path.join(FIX, rel), encoding="utf-8") as fh:
        return fh.read()


def resp(src, url: str = "https://example.test/") -> Response:
    """A Response replaying a fixture file (path under fixtures/) or a literal JSON value."""
    text = fixture(src) if isinstance(src, str) else json.dumps(src)
    return Response(url, 200, text, "2026-10-05T00:00:00+00:00", False)


# --------------------------------------------------------------------------- Europe PMC

def test_parse_hit_fields_and_open_access_links():
    data = json.loads(fixture("europepmc/search_dravet_core.json"))
    hits = [europepmc.parse_hit(r) for r in data["resultList"]["result"]]
    first = hits[0]
    assert first["pmid"] == "42814106" and first["source"] == "MED"
    assert first["authors"].endswith("et al.")
    assert first["journal"] == "Epilepsia"
    assert first["pubType"] and isinstance(first["pubType"], list)
    assert first["fullTextUrl"] is None  # subscription DOI link only
    assert first["url"] == "https://europepmc.org/article/MED/42814106"
    preprint = next(h for h in hits if h["source"] == "PPR")
    assert preprint["pmid"] is None and preprint["url"].startswith("https://europepmc.org/article/PPR/PPR")
    pmc = next(h for h in hits if h["source"] == "PMC")
    assert pmc["isOpenAccess"] is True
    assert pmc["fullTextUrl"] == "https://europepmc.org/articles/PMC13579307"  # Europe PMC html preferred
    oa = next(h for h in hits if h["pmid"] == "42632955")
    assert oa["fullTextUrl"] == "https://europepmc.org/articles/PMC13498680"


def test_parse_hit_lite_derives_full_text_from_pmcid():
    rec = {"id": "1", "source": "MED", "pmid": "1", "pmcid": "PMC9", "inEPMC": "Y", "isOpenAccess": "Y",
           "authorString": "Smith A, Jones B.", "journalTitle": "J X", "pubType": "review; journal article", "title": "<i>T</i>"}
    hit = europepmc.parse_hit(rec)
    assert hit["fullTextUrl"] == "https://europepmc.org/articles/PMC9"
    assert hit["authors"] == "Smith A et al." and hit["pubType"] == ["review", "journal article"] and hit["title"] == "T"


def test_search_passes_query_through_and_maps_sort(monkeypatch):
    seen = {}

    def fake(url, source, params=None, **kw):
        seen.update(params)
        return resp("europepmc/search_dravet_core.json", url)

    monkeypatch.setattr(europepmc, "get_json", fake)
    q = '"Dravet syndrome" AND (case report) AND PUB_YEAR:[2020 TO 2026]'
    out = europepmc.search(q, limit=4, sort="date")
    assert seen["query"] == q and seen["sort"] == "P_PDATE_D desc" and seen["resultType"] == "core"
    assert len(out.result["hits"]) == 4 and out.result["hitCount"] > 1000
    assert out.sources[0]["db"] == "Europe PMC" and out.sources[0]["record"] == q
    europepmc.search("x", limit=2, sort="relevance")
    assert seen.get("sort") is None


def test_by_pmids_keeps_input_order_and_reports_missing(monkeypatch):
    monkeypatch.setattr(europepmc, "get_json", lambda url, source, params=None, **kw: resp("europepmc/by_pmids_core.json", url))
    out = europepmc.by_pmids(["34055682", "38785537", "27397505", "99999999", "34055682", "not-a-pmid"])
    assert [p["pmid"] for p in out.result] == ["34055682", "38785537", "27397505"]
    assert any("99999999" in w for w in out.warnings)


def test_abstract_strips_markup(monkeypatch):
    monkeypatch.setattr(europepmc, "get_json", lambda url, source, params=None, **kw: resp("europepmc/abstract_35490361.json", url))
    out = europepmc.abstract("35490361")
    a = out.result
    assert a["found"] and a["pmid"] == "35490361" and "Dravet" in a["title"]
    assert a["abstract"] and "<" not in a["abstract"]
    assert a["fullTextUrl"]
    with pytest.raises(ValueError):
        europepmc.abstract("Dravet")


def test_abstract_not_found(monkeypatch):
    monkeypatch.setattr(europepmc, "get_json", lambda url, source, params=None, **kw: resp("europepmc/empty.json", url))
    out = europepmc.abstract("1")
    assert out.result == {"id": "1", "found": False} and out.warnings


# --------------------------------------------------------------------------- LitVar2

def test_spellings_cover_stop_codon_forms():
    assert litvar.spellings("p.Arg712*") == ["p.Arg712*", "p.Arg712Ter", "p.R712X"]
    assert litvar.spellings("R712X") == ["R712X", "p.Arg712Ter", "p.R712*"]
    assert litvar.spellings("p.Arg712Ter") == ["p.Arg712Ter", "p.R712X", "p.R712*"]
    assert litvar.spellings("p.Phe508del") == ["p.Phe508del"]
    assert litvar.spellings("rs794726730") == ["rs794726730"]


def _litvar_router(mapping):
    def fake(url, source, params=None, **kw):
        q = (params or {}).get("query")
        if q in mapping:
            return resp(mapping[q], url)
        return resp([], url)
    return fake


SCN1A_MAP = {"SCN1A p.Arg712*": "litvar/auto_scn1a_star.json", "SCN1A p.Arg712Ter": "litvar/auto_scn1a_ter.json",
             "SCN1A p.R712X": "litvar/auto_scn1a_x.json"}


def test_lookup_merges_records_found_under_different_spellings(monkeypatch):
    monkeypatch.setattr(litvar, "get_json", _litvar_router(SCN1A_MAP))
    out = litvar.lookup("p.Arg712*", gene="SCN1A")
    r = out.result
    assert r["queries"] == ["SCN1A p.Arg712*", "SCN1A p.Arg712Ter", "SCN1A p.R712X"]
    ids = [m["litvar_id"] for m in r["matches"]]
    assert ids == ["litvar@#6323#p.R712*", "litvar@rs794726730##"]
    rs = r["matches"][1]
    assert rs["rsid"] == "rs794726730" and rs["gene"] == "SCN1A" and rs["pmid_count"] >= 20 and rs["top"]
    assert rs["spellings"] == ["SCN1A p.Arg712Ter", "SCN1A p.R712X"]
    assert {"litvar_id", "rsid", "gene", "hgvs", "pmid_count"} <= set(rs)
    assert len(out.sources) == 3 and all(s["db"] == "LitVar2 autocomplete" for s in out.sources)


def test_lookup_marks_suggestions_that_are_not_the_first_hit(monkeypatch):
    monkeypatch.setattr(litvar, "get_json", _litvar_router({"CFTR p.Phe508del": "litvar/auto_cftr_f508del.json"}))
    r = litvar.lookup("p.Phe508del", gene="CFTR").result
    assert r["matches"][0]["rsid"] == "rs113993960" and r["matches"][0]["top"]
    assert not any(m["top"] for m in r["matches"][1:])  # other CFTR variants LitVar suggests


def test_lookup_gene_filter_drops_other_genes(monkeypatch):
    # the bare-change response (CNKSR2, TRAPPC9 records) replayed for a gene-prefixed query
    monkeypatch.setattr(litvar, "get_json", lambda url, source, params=None, **kw: resp("litvar/auto_bare_r712.json", url))
    out = litvar.lookup("p.Arg712*", gene="TRAPPC9")
    assert [m["gene"] for m in out.result["matches"]] == ["TRAPPC9"]
    assert any("other genes" in w for w in out.warnings)
    unfiltered = litvar.lookup("p.Arg712*").result["matches"]
    assert {m["gene"] for m in unfiltered} >= {"CNKSR2", "TRAPPC9"}


def test_lookup_does_not_prefix_gene_to_rsid_or_transcript_hgvs():
    assert litvar._with_gene("rs794726730", "SCN1A") == "rs794726730"
    assert litvar._with_gene("NM_001165963.4:c.2134C>T", "SCN1A") == "NM_001165963.4:c.2134C>T"
    assert litvar._with_gene("SCN1A p.R712*", "SCN1A") == "SCN1A p.R712*"
    assert litvar._with_gene("p.R712*", "SCN1A") == "SCN1A p.R712*"


def test_lookup_partial_and_total_failure(monkeypatch):
    def flaky(url, source, params=None, **kw):
        if params["query"] == "SCN1A p.Arg712*":
            raise SourceError("LitVar2", url, None, "network error: reset")
        return _litvar_router(SCN1A_MAP)(url, source, params)

    monkeypatch.setattr(litvar, "get_json", flaky)
    out = litvar.lookup("p.Arg712*", gene="SCN1A")
    assert [m["litvar_id"] for m in out.result["matches"]] == ["litvar@rs794726730##"]
    assert any("unavailable" in w for w in out.warnings)

    def down(url, source, params=None, **kw):
        raise SourceError("LitVar2", url, 503, "down")

    monkeypatch.setattr(litvar, "get_json", down)
    with pytest.raises(SourceError):
        litvar.lookup("p.Arg712*", gene="SCN1A")


def test_pmids_and_get_quote_the_id(monkeypatch):
    urls = []

    def fake(url, source, params=None, **kw):
        urls.append(url)
        return resp("litvar/pubs_rs794726730.json" if url.endswith("/publications") else "litvar/get_rs794726730.json", url)

    monkeypatch.setattr(litvar, "get_json", fake)
    out = litvar.pmids("litvar@rs794726730##", limit=5)
    assert out.result["pmid_count"] == 21 and len(out.result["pmids"]) == 5 and all(p.isdigit() for p in out.result["pmids"])
    assert urls[0].endswith("/variant/get/litvar%40rs794726730%23%23/publications")
    info = litvar.get("litvar@rs794726730##").result
    assert info["rsid"] == "rs794726730" and info["clingen_ids"] == ["CA274966"]
    assert not urls[1].endswith("/")  # the record endpoint 404s with a trailing slash


# --------------------------------------------------------------------------- PubTator3

def test_resolve_gene_entity_prefers_exact_ncbi_gene(monkeypatch):
    monkeypatch.setattr(pubtator, "get_json", lambda url, source, params=None, **kw: resp("pubtator/auto_gene_scn1a.json", url))
    out = pubtator.resolve_entity("scn1a", "GENE")
    assert out.result["entity"]["id"] == "@GENE_SCN1A" and out.result["entity"]["db_id"] == "6323"
    assert "@GENE_SCN1A.S" in out.result["alternatives"]


def test_disease_autocomplete_shows_mesh_mapping(monkeypatch):
    monkeypatch.setattr(pubtator, "get_json", lambda url, source, params=None, **kw: resp("pubtator/auto_disease_dravet.json", url))
    ents = pubtator.autocomplete("Dravet syndrome", concept="disease").result["entities"]
    assert ents[0]["db"] == "ncbi_mesh" and "<m>" not in (ents[0]["match"] or "")
    with pytest.raises(ValueError):
        pubtator.autocomplete("x", concept="protein")


def test_search_pages_until_limit(monkeypatch):
    pages = []

    def fake(url, source, params=None, **kw):
        pages.append(params["page"])
        return resp(f"pubtator/search_variant_p{params['page']}.json", url)

    monkeypatch.setattr(pubtator, "get_json", fake)
    out = pubtator.search("@VARIANT_p.R712X_SCN1A_human", limit=15)
    assert pages == [1, 2] and out.result["count"] == 12 and len(out.result["hits"]) == 12
    h = out.result["hits"][0]
    assert h["pmid"].isdigit() and h["url"] == f"https://pubmed.ncbi.nlm.nih.gov/{h['pmid']}/" and h["year"]
    pages.clear()
    out = pubtator.search("@VARIANT_p.R712X_SCN1A_human", limit=5)
    assert pages == [1] and len(out.result["hits"]) == 5


# --------------------------------------------------------------------------- zebra lit (offline end to end)

def _all_routes(monkeypatch):
    def epmc(url, source, params=None, **kw):
        q = params["query"]
        return resp("europepmc/by_pmids_core.json" if q.startswith("EXT_ID:(") else "europepmc/search_dravet_core.json", url)

    def lv(url, source, params=None, **kw):
        if url.endswith("/publications"):
            return resp("litvar/pubs_rs794726730.json", url)
        return _litvar_router(SCN1A_MAP)(url, source, params)

    def pt(url, source, params=None, **kw):
        if "autocomplete" in url:
            return resp("pubtator/auto_variant_rs794726730.json" if params.get("concept") == "VARIANT"
                        else "pubtator/auto_gene_scn1a.json", url)
        return resp(f"pubtator/search_variant_p{params['page']}.json", url)

    monkeypatch.setattr(europepmc, "get_json", epmc)
    monkeypatch.setattr(litvar, "get_json", lv)
    monkeypatch.setattr(pubtator, "get_json", pt)
    monkeypatch.delenv("ZEBRA_CASE", raising=False)


def test_cli_lit_sections(monkeypatch, capsys):
    _all_routes(monkeypatch)
    assert cli.main(["lit", "Dravet syndrome", "--gene", "SCN1A", "--variant", "p.Arg712*", "--limit", "5", "--json"]) == 0
    env = json.loads(capsys.readouterr().out)
    assert env["ok"] and env["command"] == "lit"
    r = env["result"]
    assert set(r) == {"europepmc", "litvar", "pubtator"}
    assert len(r["europepmc"]["hits"]) == 5
    assert r["litvar"]["papers"] and all("litvar_ids" in p for p in r["litvar"]["papers"])
    assert r["pubtator"]["entity"]["id"] == "@VARIANT_p.R712X_SCN1A_human"
    assert env["query"] == {"query": "Dravet syndrome", "gene": "SCN1A", "variant": "p.Arg712*", "limit": 5,
                            "sort": "relevance", "abstract": None}
    dbs = {s["db"] for s in env["sources"]}
    assert {"Europe PMC", "LitVar2 autocomplete", "LitVar2 publications", "PubTator3 autocomplete", "PubTator3 search"} <= dbs


def test_cli_lit_gene_only_uses_pubtator_gene_entity(monkeypatch, capsys):
    _all_routes(monkeypatch)
    seen = []
    orig = pubtator.search

    def spy(text, limit=10):
        seen.append(text)
        return orig(text, limit)

    monkeypatch.setattr(pubtator, "search", spy)
    assert cli.main(["lit", "--gene", "SCN1A"]) == 0
    assert seen == ["@GENE_SCN1A"]
    assert "PubTator3: @GENE_SCN1A" in capsys.readouterr().out


def test_cli_lit_needs_input(monkeypatch, capsys):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    assert cli.main(["lit", "--json"]) == 2
    assert json.loads(capsys.readouterr().out)["error"]["type"] == "UsageError"


def test_cli_lit_source_failure_is_a_warning(monkeypatch, capsys):
    _all_routes(monkeypatch)

    def down(url, source, params=None, **kw):
        raise SourceError("Europe PMC", url, 503, "Service Unavailable")

    monkeypatch.setattr(europepmc, "get_json", down)
    assert cli.main(["lit", "Dravet syndrome", "--gene", "SCN1A", "--json"]) == 0
    env = json.loads(capsys.readouterr().out)
    assert "europepmc" not in env["result"] and "pubtator" in env["result"]
    assert any(w.startswith("Europe PMC unavailable") for w in env["warnings"])


# --------------------------------------------------------------------------- live

@pytest.mark.live
def test_live_europepmc_search_and_abstract():
    out = europepmc.search('"Dravet syndrome"', limit=3)
    assert out.result["hitCount"] > 1000 and len(out.result["hits"]) == 3
    a = europepmc.abstract("35490361").result
    assert "Dravet" in a["title"] and a["abstract"]


@pytest.mark.live
def test_live_litvar_scn1a_stop_gain():
    r = litvar.lookup("p.Arg712*", gene="SCN1A").result
    rs = [m for m in r["matches"] if m["rsid"] == "rs794726730"]
    assert rs and rs[0]["top"] and rs[0]["pmid_count"] >= 10
    assert len(litvar.pmids("litvar@rs794726730##", 5).result["pmids"]) == 5


@pytest.mark.live
def test_live_litvar_cftr_f508del_and_rsid():
    r = litvar.lookup("p.Phe508del", gene="CFTR").result
    assert r["matches"][0]["rsid"] == "rs113993960" and r["matches"][0]["pmid_count"] > 1000
    r = litvar.lookup("NM_001165963.4:c.2134C>T").result
    assert r["matches"] and r["matches"][0]["rsid"] == "rs794726730"


@pytest.mark.live
def test_live_pubtator_gene_entity_and_search():
    ent = pubtator.resolve_entity("NGLY1", "GENE").result["entity"]
    assert ent["id"] == "@GENE_NGLY1"
    res = pubtator.search(ent["id"], limit=3).result
    assert res["count"] > 50 and len(res["hits"]) == 3


@pytest.mark.live
def test_live_cli_lit(capsys, monkeypatch):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    assert cli.main(["lit", "Dravet syndrome", "--limit", "3", "--json"]) == 0
    env = json.loads(capsys.readouterr().out)
    assert len(env["result"]["europepmc"]["hits"]) == 3 and not env["warnings"]
