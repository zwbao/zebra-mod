"""Literature sources (Europe PMC, LitVar2, PubTator3) and `zebra lit`.

Offline tests replay real responses captured from the services (trimmed) in
tests/fixtures/{europepmc,litvar,pubtator}; live tests call the services.
"""

from __future__ import annotations

import json
import math
import os

import pytest

from zebra import cli
from zebra.commands import lit as lit_cmd
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
        # the 15 PMIDs LitVar links to rs794726730, with the abstracts resultType=core returns
        return resp("europepmc/by_pmids_scn1a_r712_core.json" if q.startswith("EXT_ID:(")
                    else "europepmc/search_dravet_core.json", url)

    def lv(url, source, params=None, **kw):
        if url.endswith("/publications"):
            # each LitVar record lists its own papers: the p.R712* record one, the rsID record 21
            return resp("litvar/pubs_6323_p_r712star.json" if "6323" in url else "litvar/pubs_rs794726730.json", url)
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


# --------------------------------------------------------------------------- F40: is the paper about the thing asked for?
# `zebra lit --gene SCN1A --variant "p.Arg712*"` lists PMID 27397505 (Iorio 2016, Cell, cancer cell-line
# pharmacogenomics) because LitVar normalises p.R712* onto rs794726730, whose PubTator counterpart is
# numbered p.R701X on another transcript. Every hit must now carry the handle that produced it and what
# the paper's own title and abstract say, so a reader can judge without the tool's word for it.

SCN1A_FORMS = ["Arg712*", "Arg712Ter", "Arg712X", "R712*", "R712Ter", "R712X", "rs794726730"]


def _lit_json(capsys, monkeypatch, *argv):
    _all_routes(monkeypatch)
    assert cli.main(["lit", *argv, "--json"]) == 0
    return json.loads(capsys.readouterr().out)


def test_F40_variant_forms_cover_both_letter_codes_and_the_rsid():
    assert lit_cmd._variant_forms("p.Arg712*", ["rs794726730"]) == SCN1A_FORMS
    assert lit_cmd._variant_forms("R712X") == ["R712*", "R712Ter", "R712X", "Arg712*", "Arg712Ter", "Arg712X"]
    assert lit_cmd._variant_forms("p.Arg712Cys") == ["Arg712Cys", "R712C"]
    assert lit_cmd._variant_forms("NM_001165963.4:c.2134C>T") == ["NM_001165963.4:c.2134C>T"]  # no protein codes to twin
    assert lit_cmd._variant_forms(None, ["rs794726730"]) == ["rs794726730"]


def test_F40_matched_entity_attached(capsys, monkeypatch):
    env = _lit_json(capsys, monkeypatch, "--gene", "SCN1A", "--variant", "p.Arg712*")
    papers = env["result"]["litvar"]["papers"]
    assert papers and len(papers) == 15
    for p in papers:
        e = p["matched_entity"]
        assert e["source"] == "LitVar2"
        assert e["id"] in [m["litvar_id"] for m in env["result"]["litvar"]["matches"]]
        assert e["id"] == p["litvar_ids"][0]
        assert e["via"] in env["result"]["litvar"]["queries"]
    iorio = next(p for p in papers if p["pmid"] == "27397505")
    assert iorio["matched_entity"] == {"source": "LitVar2", "id": "litvar@rs794726730##", "name": "c.2134C>T",
                                       "rsid": "rs794726730", "via": "SCN1A p.Arg712Ter"}
    hits = env["result"]["pubtator"]["hits"]
    assert hits and all(h["matched_entity"]["source"] == "PubTator3" for h in hits)
    e = hits[0]["matched_entity"]
    # the id says p.R712X, the entity's own name says p.R701X: the transcript renumbering, on the record
    assert e == {"source": "PubTator3", "id": "@VARIANT_p.R712X_SCN1A_human", "name": "p.R701X",
                 "rsid": "rs794726730", "via": "rs794726730"}


def test_F40_mentions_detects_absent_gene(capsys, monkeypatch):
    env = _lit_json(capsys, monkeypatch, "--gene", "SCN1A", "--variant", "p.Arg712*")
    lit = env["result"]["litvar"]
    assert lit["variant_forms_checked"] == SCN1A_FORMS
    by_pmid = {p["pmid"]: p for p in lit["papers"]}
    iorio = by_pmid["27397505"]
    assert iorio["mentions"] == {"gene": False, "variant": False, "variant_forms_found": []}
    assert iorio["abstract_excerpt"].startswith("Systematic studies of cancer genomes")
    assert "SCN1A" not in iorio["title"]
    # a paper that does name the gene is told apart from one that does not
    assert by_pmid["38785537"]["mentions"]["gene"] is True
    assert by_pmid["31253177"]["mentions"]["gene"] is False  # colorectal neoantigen paper
    assert lit["mention_counts"] == {"papers_checked": 15, "no_abstract": 0, "gene_mentioned": 6,
                                     "gene_not_mentioned": 9, "gene_unknown": 0, "variant_mentioned": 0,
                                     "variant_not_mentioned": 15, "variant_unknown": 0}
    # no paper's abstract names the variant under any spelling — a fact, not a classification
    assert all(p["mentions"]["variant_forms_found"] == [] for p in lit["papers"])
    assert "title and abstract only" in lit["note"] and "supplementary tables" in lit["note"]


def test_F40_rendering_shows_tag_and_entity(capsys, monkeypatch):
    _all_routes(monkeypatch)
    assert cli.main(["lit", "--gene", "SCN1A", "--variant", "p.Arg712*"]) == 0
    text = capsys.readouterr().out
    assert "PMID 27397505" in text
    iorio_line = next(ln for ln in text.splitlines() if "PMID 27397505" in ln)
    assert "[gene✗, variant✗]" in iorio_line and "via litvar@rs794726730##" in iorio_line
    assert "[gene✓, variant✗] via litvar@" in text  # a paper that does name SCN1A
    assert "6 mention SCN1A, 9 do not" in text and "title and abstract only" in text
    assert lit_cmd._mention_tag({}) == "" and lit_cmd._mention_tag({"title": "x"}) == ""  # EPMC search hits untagged


def test_F40_no_abstract_is_null_not_false(monkeypatch):
    monkeypatch.setattr(europepmc, "get_json",
                        lambda url, source, params=None, **kw: resp("europepmc/by_pmids_no_abstract_core.json", url))
    papers = europepmc.by_pmids(["38327537", "27397505"]).result
    assert [p["pmid"] for p in papers] == ["38327537", "27397505"]
    assert papers[0]["abstract"] is None  # a Comment; Europe PMC holds no abstract for it
    assert papers[1]["abstract"] and "SCN1A" not in papers[1]["abstract"]
    lit_cmd._attach_mentions(papers, "SCN1A", SCN1A_FORMS)
    assert papers[0]["mentions"] == {"gene": None, "variant": None, "variant_forms_found": []}
    assert papers[0]["abstract_excerpt"] is None
    assert lit_cmd._mention_tag(papers[0]) == "[no abstract]"
    # the same answer for the same gene, when there was an abstract to look in, is false — not null
    assert papers[1]["mentions"]["gene"] is False and papers[1]["mentions"]["variant"] is False
    assert lit_cmd._mention_tag(papers[1]) == "[gene✗, variant✗]"
    assert lit_cmd._mention_counts(papers) == {"papers_checked": 2, "no_abstract": 1, "gene_mentioned": 0,
                                               "gene_not_mentioned": 1, "gene_unknown": 1, "variant_mentioned": 0,
                                               "variant_not_mentioned": 1, "variant_unknown": 1}
    # a no-abstract paper whose title names the gene is a `true`, so gene_unknown alone would understate it
    titled_no_abstract = [{"title": "SCN1A review", "abstract": None}]
    lit_cmd._attach_mentions(titled_no_abstract, "SCN1A", SCN1A_FORMS)
    counted = lit_cmd._mention_counts(titled_no_abstract)
    assert counted["no_abstract"] == 1 and counted["gene_unknown"] == 0 and counted["gene_mentioned"] == 1
    # a title Europe PMC does supply still counts: no abstract does not blind the title
    titled = [{"title": "SCN1A-related epilepsy with a p.Arg712Ter allele", "abstract": None}]
    lit_cmd._attach_mentions(titled, "SCN1A", SCN1A_FORMS)
    assert titled[0]["mentions"] == {"gene": True, "variant": True, "variant_forms_found": ["Arg712Ter"]}
    assert lit_cmd._mention_tag(titled[0]) == "[gene✓, variant✓, no abstract]"
    # whole-word matching: SCN1AB is a different gene, SCN1A-related is the same one
    assert lit_cmd._mentions("SCN1A", [], "SCN1AB variants", "a b")["gene"] is False
    assert lit_cmd._mentions("SCN1A", [], "SCN1A-related epilepsy", "a b")["gene"] is True
    # nothing asked about is null too, and distinguishable by the absence of forms
    assert lit_cmd._mentions(None, [], "t", "a") == {"gene": None, "variant": None, "variant_forms_found": []}


def test_F40_warning_when_most_hits_unmatched(capsys, monkeypatch):
    env = _lit_json(capsys, monkeypatch, "--gene", "SCN1A", "--variant", "p.Arg712*")
    warning = next(w for w in env["warnings"] if w.startswith("LitVar2: 9 of 15"))
    assert "do not mention SCN1A in their title or abstract" in warning
    assert "normalised variant records that can merge different transcript numbering" in warning
    assert "check each PMID before citing it" in warning
    # the threshold is a third of the papers that could be checked, with a floor of two
    assert lit_cmd.UNMATCHED_WARN_FRACTION == 1.0 / 3 and lit_cmd.UNMATCHED_WARN_FLOOR == 2
    assert env["result"]["litvar"]["mention_counts"]["gene_not_mentioned"] >= math.ceil(15 / 3)


def test_F40_warning_silent_when_the_list_is_clean(monkeypatch, capsys):
    """One stray paper out of fifteen is shown per line, not shouted about."""
    _all_routes(monkeypatch)
    clean = [{"title": "SCN1A paper", "abstract": "about SCN1A"} for _ in range(14)]
    clean.append({"title": "something else", "abstract": "no gene here"})
    lit_cmd._attach_mentions(clean, "SCN1A", SCN1A_FORMS)
    counts = lit_cmd._mention_counts(clean)
    assert counts["gene_not_mentioned"] == 1
    assert counts["gene_not_mentioned"] < max(lit_cmd.UNMATCHED_WARN_FLOOR,
                                              math.ceil(counts["papers_checked"] * lit_cmd.UNMATCHED_WARN_FRACTION))


def test_F40_gene_only_query_has_no_variant_half(capsys, monkeypatch):
    env = _lit_json(capsys, monkeypatch, "--gene", "SCN1A")
    hits = env["result"]["pubtator"]["hits"]
    assert env["result"]["pubtator"]["variant_forms_checked"] == []
    assert all(h["mentions"]["variant"] is None for h in hits)
    assert all(h["mentions"]["variant_forms_found"] == [] for h in hits)
    assert any(h["mentions"]["gene"] is True for h in hits)
    assert lit_cmd._mention_tag(hits[0]).startswith("[gene")
    assert "variant" not in lit_cmd._mention_tag(hits[0])


def test_F40_abstracts_are_trimmed_out_of_the_output(capsys, monkeypatch):
    env = _lit_json(capsys, monkeypatch, "--gene", "SCN1A", "--variant", "p.Arg712*")
    for p in env["result"]["litvar"]["papers"] + env["result"]["pubtator"]["hits"]:
        assert "abstract" not in p  # the full text is read, then dropped
        excerpt = p["abstract_excerpt"]
        assert excerpt is None or len(excerpt) <= lit_cmd.ABSTRACT_EXCERPT + 3


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
def test_live_F40_by_pmids_returns_abstracts_only_with_result_type_core():
    """The claim in europepmc.by_pmids' docstring, checked against the service."""
    core = europepmc.by_pmids(["27397505", "38785537"]).result
    assert [p["pmid"] for p in core] == ["27397505", "38785537"]
    assert all(p["abstract"] and len(p["abstract"]) > 500 for p in core)
    lite = europepmc.by_pmids(["27397505"], result_type="lite").result
    assert lite[0]["title"] and lite[0]["abstract"] is None


@pytest.mark.live
def test_live_F40_scn1a_stop_gain_labels_the_irrelevant_pmid(capsys, monkeypatch):
    """The original F40 reproduction: the Iorio 2016 cancer paper is now labelled, not just listed."""
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    assert cli.main(["lit", "--gene", "SCN1A", "--variant", "p.Arg712*", "--json"]) == 0
    env = json.loads(capsys.readouterr().out)
    lit = env["result"]["litvar"]
    iorio = next(p for p in lit["papers"] if p["pmid"] == "27397505")
    assert "Pharmacogenomic" in iorio["title"]
    assert iorio["mentions"]["gene"] is False  # not falsy: false, because there was an abstract to look in
    assert iorio["mentions"]["variant"] is False
    assert iorio["mentions"]["variant_forms_found"] == []
    assert iorio["matched_entity"]["source"] == "LitVar2"
    assert iorio["matched_entity"]["id"] == "litvar@rs794726730##"
    assert iorio["matched_entity"]["rsid"] == "rs794726730"
    assert "712" in iorio["matched_entity"]["via"]
    counts = lit["mention_counts"]
    assert counts["gene_not_mentioned"] >= max(lit_cmd.UNMATCHED_WARN_FLOOR,
                                               math.ceil(counts["papers_checked"] * lit_cmd.UNMATCHED_WARN_FRACTION))
    warning = next(w for w in env["warnings"] if w.startswith("LitVar2:") and "do not mention" in w)
    assert "SCN1A" in warning and "transcript numbering" in warning
    # at least one genuinely SCN1A paper is in the same list and is labelled differently
    assert any(p["mentions"]["gene"] is True for p in lit["papers"])


@pytest.mark.live
def test_live_F40_pubtator_entity_name_shows_the_renumbering(capsys, monkeypatch):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    assert cli.main(["lit", "--gene", "SCN1A", "--variant", "p.Arg712*", "--json"]) == 0
    hits = json.loads(capsys.readouterr().out)["result"]["pubtator"]["hits"]
    assert hits
    e = hits[0]["matched_entity"]
    assert e["source"] == "PubTator3" and e["id"] == "@VARIANT_p.R712X_SCN1A_human"
    assert e["name"] == "p.R701X"  # the same change numbered on another transcript
    assert all(isinstance(h["mentions"]["variant_forms_found"], list) for h in hits)
    assert any(h["mentions"]["gene"] is True for h in hits)


@pytest.mark.live
def test_live_cli_lit(capsys, monkeypatch):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    assert cli.main(["lit", "Dravet syndrome", "--limit", "3", "--json"]) == 0
    env = json.loads(capsys.readouterr().out)
    assert len(env["result"]["europepmc"]["hits"]) == 3 and not env["warnings"]
