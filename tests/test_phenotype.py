import json
import os
from pathlib import Path

import pytest

from zebra import cli
from zebra.commands import phenotype as P
from zebra.core import Outcome
from zebra.sources import hpo as hpo_source
from zebra.sources import monarch, ols, pubcasefinder

FIX = Path(__file__).parent / "fixtures"
QUERY = ["HP:0002373", "HP:0007359", "HP:0002133", "HP:0001263", "HP:0007207"]
REAL_HPO = Path(os.path.expanduser("~/.cache/zebra-mod/hpo"))


def _json(rel):
    return json.loads((FIX / rel).read_text("utf-8"))


def _text(rel):
    return (FIX / rel).read_text("utf-8")


# ---------------------------------------------------------------- parsers


def test_parse_semsim_dravet_rank_and_matched_terms():
    hits = monarch.parse_semsim(_json("monarch/semsim_dravet_query.json"), 10)
    assert [h["rank"] for h in hits] == [1, 2, 3, 4]
    dravet = hits[2]
    assert dravet["id"] == "MONDO:0100135" and dravet["name"] == "Dravet syndrome"
    assert dravet["score"] == pytest.approx(13.884, abs=1e-3)
    by_q = {m["query"]: m for m in dravet["matched"]}
    assert set(by_q) == set(QUERY)
    assert by_q["HP:0002373"]["exact"] is True
    assert by_q["HP:0007359"]["via"] == "HP:0001250" and by_q["HP:0007359"]["exact"] is False


def test_parse_pubcasefinder_tsv_keeps_served_ties():
    omim = pubcasefinder.parse_tsv(_text("pubcasefinder/rank_dravet_query_omim.tsv"), "omim", 15)
    assert len(omim) == 15
    assert [r["rank"] for r in omim[:4]] == [1, 1, 3, 3]
    assert omim[2]["id"] == "OMIM:615744" and omim[2]["genes"] == ["GABRA1"]
    assert omim[0]["matched"][0].startswith("HP:")
    orpha = pubcasefinder.parse_tsv(_text("pubcasefinder/rank_dravet_query_orphanet.tsv"), "orphanet", 5)
    assert orpha[0]["id"] == "ORPHA:307" and orpha[0]["name"] == "juvenile myoclonic epilepsy"
    genes = pubcasefinder.parse_tsv(_text("pubcasefinder/rank_dravet_query_gene.tsv"), "gene", 3)
    assert genes[0]["gene_id"].startswith("GENEID:") and genes[0]["symbol"]


def test_parse_tsv_rejects_non_tsv():
    with pytest.raises(ValueError):
        pubcasefinder.parse_tsv("<html>error</html>", "omim", 5)


# ---------------------------------------------------------------- consensus / id reconciliation


def _items():
    return [
        {"source": "local", "list": "local", "id": "ORPHA:33069", "name": "Dravet syndrome", "rank": 1},
        {"source": "local", "list": "local", "id": "OMIM:615744", "name": "DEE 19", "rank": 2},
        {"source": "local", "list": "local", "id": "OMIM:607208", "name": "EIEE6 (Dravet syndrome)", "rank": 8},
        {"source": "monarch", "list": "monarch", "id": "MONDO:0100135", "name": "Dravet syndrome", "rank": 3, "xref": []},
        {"source": "pubcasefinder", "list": "pubcasefinder omim", "id": "OMIM:615744", "name": "DEE 19", "rank": 3},
    ]


def test_consensus_joins_only_through_exact_links():
    rows = monarch.parse_mappings(_json("monarch/mappings_dravet.json"))
    keys = P.resolve_keys(_items(), rows, {"MONDO:0011794": "MONDO:0100135"})
    # ORPHA:33069 -> obsolete MONDO:0011794 -> replaced by MONDO:0100135 (Monarch's hit)
    assert keys["ORPHA:33069"][0] == "MONDO:0100135"
    assert "obsolete" in keys["ORPHA:33069"][1] and "OLS" in keys["ORPHA:33069"][1]
    # OMIM:607208 is MONDO:0100079 (DEE6A), which Mondo keeps distinct from Dravet syndrome
    assert keys["OMIM:607208"][0] == "MONDO:0100079"
    cons = P.build_consensus(_items(), keys)
    by_key = {c["key"]: c for c in cons}
    assert by_key["MONDO:0100135"]["best_rank"] == {"local": 1, "monarch": 3}
    # DEE19 (Orphanet maps it only as BTNT to ORPHA:33069) stays its own row
    assert by_key["MONDO:0014328"]["best_rank"] == {"local": 2, "pubcasefinder": 3}
    assert "MONDO:0100079" not in by_key  # only one source
    assert cons[0]["key"] == "MONDO:0100135"  # ordered by #sources then worst rank: (1,3) before (2,3)


def test_consensus_without_mappings_uses_same_ids_only():
    keys = P.resolve_keys(_items(), None, None)
    cons = P.build_consensus(_items(), keys)
    assert [c["key"] for c in cons] == ["OMIM:615744"]
    assert cons[0]["linked_by"] == ["same id"]


def test_obsolete_mapping_without_replacement_is_not_linked():
    rows = monarch.parse_mappings(_json("monarch/mappings_dravet.json"))
    keys = P.resolve_keys(_items(), rows, {})
    assert keys["ORPHA:33069"][0] == "ORPHA:33069"


def test_gene_consensus():
    per = {"local": {"genes": [{"rank": 1, "symbol": "SCN1A"}, {"rank": 2, "symbol": "PCDH19"}]},
           "monarch": {"genes": [{"rank": 3, "symbol": "PCDH19"}]},
           "pubcasefinder": {"genes": [{"rank": 7, "symbol": "scn1a"}, {"rank": 1, "symbol": "PCDH19"}]}}
    rows = P.gene_consensus(per)
    assert rows[0]["symbol"] == "PCDH19" and rows[0]["n_sources"] == 3
    assert rows[1] == {"symbol": "SCN1A", "n_sources": 2, "rank": {"local": 1, "pubcasefinder": 7}}


# ---------------------------------------------------------------- inputs


def test_parse_sources_and_bad_input():
    assert P.parse_sources(None) == ["local", "monarch", "pubcasefinder"]
    assert P.parse_sources("monarch, local") == ["monarch", "local"]
    with pytest.raises(P.UsageError):
        P.parse_sources("local,omim")


def test_cli_usage_errors(capsys):
    assert cli.main(["phenotype", "rank", "--present", "HP:12"]) == 2
    assert cli.main(["phenotype", "rank"]) == 2
    assert cli.main(["--json", "phenotype", "rank", "--from-case"]) == 2
    err = json.loads(capsys.readouterr().out.strip().splitlines()[-1])
    assert err["ok"] is False and "active case" in err["error"]["message"]


# ---------------------------------------------------------------- whole command, offline (fixtures)


@pytest.fixture
def offline_sources(monkeypatch):
    semsim = _json("monarch/semsim_dravet_query.json")

    def fake_semsim(present, limit=15, group="Human Diseases", metric="ancestor_information_content"):
        hits = monarch.parse_semsim(semsim, limit)
        if group == "Human Genes":
            hits = [dict(h, symbol="GENE%d" % h["rank"]) for h in hits]
        return Outcome({"group": group, "metric": metric, "hits": hits},
                       sources=[{"db": "Monarch semsim search", "record": group, "url": "fixture"}])

    def fake_pcf(hpo_ids, target="omim", limit=15):
        rows = pubcasefinder.parse_tsv(_text(f"pubcasefinder/rank_dravet_query_{target}.tsv"), target, limit)
        return Outcome({"target": target, "hits": rows}, sources=[{"db": "PubCaseFinder", "record": target, "url": "fixture"}])

    def fake_mappings(object_ids=(), subject_ids=(), predicate="skos:exactMatch"):
        return Outcome(monarch.parse_mappings(_json("monarch/mappings_dravet.json")),
                       sources=[{"db": "Monarch mappings", "record": "fixture", "url": "fixture"}])

    def fake_obsolete(curies):
        t = ols.parse_term(_json("ols/term_MONDO_0011794.json"))
        return Outcome({c: (t["replaced_by"] if c == t["id"] else c) for c in curies},
                       sources=[{"db": "OLS", "record": "fixture", "url": "fixture"}])

    def fake_term(hpo_id, prefer_local=True):
        # F14 added an HPO API check for ids that no local release can resolve;
        # an offline test must answer it from here rather than over the network.
        return Outcome({"id": hpo_id, "name": f"term {hpo_id}", "obsolete": False, "replaced_by": None})

    monkeypatch.setattr(monarch, "semsim_rank", fake_semsim)
    monkeypatch.setattr(pubcasefinder, "ranked", fake_pcf)
    monkeypatch.setattr(monarch, "mappings", fake_mappings)
    monkeypatch.setattr(hpo_source, "term", fake_term)
    monkeypatch.setattr(ols, "resolve_obsolete", fake_obsolete)


def test_rank_offline_without_local_hpo_warns_and_keeps_sources_separate(offline_sources, monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(tmp_path / "no-hpo"))
    code = cli.main(["--json", "phenotype", "rank", "--present", *QUERY, "--exclude", "HP:0001252", "--top", "10"])
    env = json.loads(capsys.readouterr().out)
    assert code == 0 and env["ok"]
    r = env["result"]
    assert r["sources"] == ["monarch", "pubcasefinder"]
    assert any("zebra hpo fetch" in w for w in env["warnings"])
    assert any("does not take excluded" in n for n in r["notes"])
    assert set(r["per_source"]) == {"monarch", "pubcasefinder"}
    assert r["per_source"]["monarch"]["excluded_supported"] is False
    assert len(r["per_source"]["pubcasefinder"]["diseases_omim"]) == 10
    assert "score" not in json.dumps(r["consensus"])  # no combined score
    assert {s["db"] for s in env["sources"]} >= {"Monarch semsim search", "PubCaseFinder"}


@pytest.mark.skipif(not (REAL_HPO / "phenotype.hpoa").exists(), reason="local HPO release not fetched")
def test_rank_offline_with_local_hpo_puts_dravet_in_consensus(offline_sources, monkeypatch, capsys):
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(REAL_HPO))
    code = cli.main(["--json", "phenotype", "rank", *QUERY, "--top", "15"])
    r = json.loads(capsys.readouterr().out)["result"]
    assert code == 0
    assert r["per_source"]["local"]["diseases"][0]["id"] == "ORPHA:33069"
    dravet = next(c for c in r["consensus"] if c["key"] == "MONDO:0100135")
    assert dravet["best_rank"]["local"] == 1 and dravet["best_rank"]["monarch"] == 3


def test_from_case_reads_present_and_excluded(offline_sources, monkeypatch, tmp_path, capsys):
    from zebra import case as case_mod

    monkeypatch.setenv("ZEBRA_HPO_DIR", str(tmp_path / "no-hpo"))
    cdir = tmp_path / "c1"
    case_mod.init(str(cdir))
    for t in QUERY[:3]:
        case_mod.add_phenotype(str(cdir), t, "x")
    case_mod.add_phenotype(str(cdir), "HP:0001252", "Hypotonia", status="excluded")
    code = cli.main(["--json", "phenotype", "rank", "--from-case", "--case", str(cdir), "--sources", "monarch"])
    env = json.loads(capsys.readouterr().out)
    assert code == 0
    assert env["query"]["present"] == QUERY[:3] and env["query"]["excluded"] == ["HP:0001252"]
    assert env["ledger"]  # sources appended to the case's evidence ledger


# ---------------------------------------------------------------- live


@pytest.mark.live
def test_live_monarch_semsim_dravet_top5():
    out = monarch.semsim_rank(QUERY, limit=10)
    ids = [h["id"] for h in out.result["hits"][:5]]
    assert "MONDO:0100135" in ids
    assert out.sources[0]["url"].startswith("https://api-v3.monarchinitiative.org/")


@pytest.mark.live
def test_live_pubcasefinder_lists():
    out = pubcasefinder.ranked(["HP:0001250", "HP:0001263"], "orphanet", 5)
    assert len(out.result["hits"]) == 5 and all(h["id"].startswith("ORPHA:") for h in out.result["hits"])


@pytest.mark.live
@pytest.mark.skipif(not (REAL_HPO / "phenotype.hpoa").exists(), reason="local HPO release not fetched")
def test_live_acceptance_dravet_top5_in_two_sources(monkeypatch, capsys):
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(REAL_HPO))
    code = cli.main(["--json", "phenotype", "rank", "--present", *QUERY])
    r = json.loads(capsys.readouterr().out)["result"]
    assert code == 0
    in_top5 = []
    if any(d["id"] == "ORPHA:33069" for d in r["per_source"]["local"]["diseases"][:5]):
        in_top5.append("local")
    if any(d["id"] == "MONDO:0100135" for d in r["per_source"]["monarch"]["diseases"][:5]):
        in_top5.append("monarch")
    pcf = r["per_source"].get("pubcasefinder") or {}
    if any(d["id"] in ("OMIM:607208", "ORPHA:33069") for d in pcf.get("diseases_omim", [])[:5] + pcf.get("diseases_orphanet", [])[:5]):
        in_top5.append("pubcasefinder")
    assert len(in_top5) >= 2, in_top5
    assert any(c["key"] == "MONDO:0100135" for c in r["consensus"])


# ---------------------------------------------------------------- E9: HPO search without local files

def test_E9_ols_ranking_puts_the_exact_label_first():
    """The JAX search endpoint never surfaces HP:0001250 for "seizure"; OLS4, re-ranked, does."""
    from zebra.sources import hpo as hpo_src

    ols_raw = _json("hpo/e9_ols_seizure.json")
    docs = ols_raw["response"]["docs"]
    hits = [{"id": d.get("obo_id"), "label": d.get("label"), "synonyms": (d.get("synonym") or [])[:5]}
            for d in docs if str(d.get("obo_id") or "").startswith("HP:") and not d.get("is_obsolete")]
    ranked = hpo_src.rank_hits("seizure", hits, 5)
    assert ranked[0]["id"] == "HP:0001250" and ranked[0]["label"] == "Seizure"
    assert ranked[0]["matched_on"] == "label"


def test_E9_the_jax_endpoint_does_not_serve_the_exact_term_at_all():
    """Why the online backend is OLS4: HP:0001250 is on no JAX page for "seizure"."""
    from zebra.sources import hpo as hpo_src

    jax = _json("hpo/e9_jax_seizure.json")
    hits = [{"id": t.get("id"), "label": t.get("name"), "synonyms": (t.get("synonyms") or [])[:5]}
            for t in jax["terms"]]
    assert len(hits) == 50 and not any(h["id"] == "HP:0001250" for h in hits)
    assert hits[0]["label"] == "Maternal seizure"
    # re-ranking cannot invent a term the service never returned, so it must only
    # ever be the fallback; every hit it does return says what matched
    ranked = hpo_src.rank_hits("seizure", hits, 3)
    assert all(h["matched"] and h["matched_on"] in ("label", "synonym") for h in ranked)


def test_E9_rank_hits_order_is_exact_then_synonym_then_prefix_then_word():
    from zebra.sources import hpo as hpo_src

    hits = [
        {"id": "HP:5", "label": "Hypocalcemic fits", "synonyms": ["Fits with low calcium and seizure risk"]},
        {"id": "HP:4", "label": "Maternal seizure", "synonyms": []},
        {"id": "HP:3", "label": "Seizure cluster", "synonyms": []},
        {"id": "HP:2", "label": "Epilepsy", "synonyms": ["Seizure"]},
        {"id": "HP:1", "label": "Seizure", "synonyms": []},
    ]
    order = [h["id"] for h in hpo_src.rank_hits("seizure", hits, 5)]
    assert order[:2] == ["HP:1", "HP:2"]          # exact label, then exact synonym
    assert order.index("HP:3") < order.index("HP:4")  # prefix before mid-label word
    assert order[-1] == "HP:5"                    # substring inside a much longer wording, last


def test_E9_search_falls_back_to_jax_with_a_warning(monkeypatch):
    from zebra.http import Response, SourceError
    from zebra.sources import hpo as hpo_src

    def dead_ols(text, rows=100):
        raise SourceError("OLS (HPO search)", "u", 503, "down")

    def jax(text, rows=50):
        jx = _json("hpo/e9_jax_seizure.json")
        hits = [{"id": t["id"], "label": t["name"], "synonyms": []} for t in jx["terms"]]
        return Response("https://ontology.jax.org/api/hp/search", 200, "{}", "2026-10-06T00:00:00+00:00", True), hits

    monkeypatch.setattr(hpo_src, "_ols_search", dead_ols)
    monkeypatch.setattr(hpo_src, "_jax_search", jax)
    out = hpo_src.search("seizure", 3)
    assert out.result["hits"]
    assert any("OLS4 HPO search unavailable" in w for w in out.warnings)


def test_E9_non_ascii_query_is_not_sent_to_the_jax_endpoint(monkeypatch):
    from zebra.http import Response
    from zebra.sources import hpo as hpo_src

    def empty_ols(text, rows=100):
        return Response("https://www.ebi.ac.uk/ols4/api/search", 200, "{}", "2026-10-06T00:00:00+00:00", True), []

    def boom(text, rows=50):
        raise AssertionError("the JAX endpoint rejects non-ASCII and must not be called")

    monkeypatch.setattr(hpo_src, "_ols_search", empty_ols)
    monkeypatch.setattr(hpo_src, "_jax_search", boom)
    out = hpo_src.search("发育迟缓", 3)
    assert out.result["hits"] == []
    assert any("is not ASCII" in w for w in out.warnings)


# -------------------------------------- F14: `phenotype rank` with no usable terms fails usefully

def test_F14_no_usable_term_raises_instead_of_an_empty_monarch_query(monkeypatch):
    from zebra.core import UsageError
    from zebra.sources import hpo as hpo_src

    def not_found(hpo_id, prefer_local=True):
        return Outcome({"id": hpo_id, "name": None, "obsolete": None, "note": "not found in HPO"})

    def boom(*a, **kw):
        raise AssertionError("a ranker must not be called with an empty term set")

    monkeypatch.setattr(hpo_src, "term", not_found)
    monkeypatch.setattr(P, "run_monarch", boom)
    monkeypatch.setattr(P, "run_pubcasefinder", boom)
    args = _rank_args(["HP:9999999"], "local,monarch")
    with pytest.raises(UsageError) as err:
        P._rank(args)
    text = str(err.value)
    assert "HP:9999999" in text and "cannot be queried with an empty set" in text
    assert "zebra hpo term" in text


def test_F14_a_resolvable_term_still_runs(monkeypatch):
    from zebra.sources import hpo as hpo_src

    def found(hpo_id, prefer_local=True):
        return Outcome({"id": hpo_id, "name": "Seizure", "obsolete": False, "replaced_by": None})

    called = {}

    def fake_monarch(terms, top, warn):
        called["terms"] = list(terms)
        return Outcome({"diseases": [], "genes": []})

    monkeypatch.setattr(hpo_src, "term", found)
    monkeypatch.setattr(P, "run_monarch", fake_monarch)
    out = P._rank(_rank_args(["HP:0001250"], "monarch"))
    assert called["terms"] == ["HP:0001250"] and out.result["sources"] == ["monarch"]


def _rank_args(present, sources):
    import argparse

    return argparse.Namespace(hpo=[], present=list(present), exclude=[], from_case=False,
                              sources=sources, top=5, case=None)


@pytest.mark.live
def test_live_E9_seizure_without_local_files(tmp_path, monkeypatch):
    """The acceptance for E9: HP:0001250 first, from the web, on a fresh install."""
    from zebra.sources import hpo as hpo_src

    monkeypatch.setenv("ZEBRA_HPO_DIR", str(tmp_path / "no-hpo"))
    out = hpo_src.search("seizure", 3)
    assert out.result["hits"][0]["id"] == "HP:0001250"
    assert out.result["hits"][0]["label"] == "Seizure"
    assert "OLS4" in out.result["backend"]
    assert out.sources[0]["db"] == "HPO search (EBI OLS4)"


@pytest.mark.live
@pytest.mark.skipif(not (REAL_HPO / "phenotype.hpoa").exists(), reason="local HPO release not fetched")
def test_live_E9_seizure_with_local_files(monkeypatch, capsys):
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(REAL_HPO))
    code = cli.main(["--json", "hpo", "search", "seizure", "--limit", "3"])
    env = json.loads(capsys.readouterr().out)
    assert code == 0 and env["result"]["hits"][0]["id"] == "HP:0001250"


def test_F10_an_empty_pubcasefinder_answer_is_not_no_match(monkeypatch):
    """An empty 200 was cached for 7 days and read as "PubCaseFinder found nothing"."""
    from zebra.http import Response, SourceError

    ttls = []

    def fake_request(url, source, **kw):
        ttls.append(kw.get("cache_ttl"))
        return Response(url, 200, "", "2026-10-06T00:00:00+00:00", len(ttls) == 1)

    monkeypatch.setattr(pubcasefinder, "request", fake_request)
    with pytest.raises(SourceError) as err:
        pubcasefinder.ranked(["HP:0001250"], "omim", 5)
    assert ttls == [7 * 86400, 0]  # the cached empty body was refetched before giving up
    assert "empty body" in err.value.message


def test_F10_an_html_pubcasefinder_answer_is_rejected(monkeypatch):
    from zebra.http import Response, SourceError

    def fake_request(url, source, **kw):
        return Response(url, 200, "<html>maintenance</html>", "2026-10-06T00:00:00+00:00", False)

    monkeypatch.setattr(pubcasefinder, "request", fake_request)
    with pytest.raises(SourceError) as err:
        pubcasefinder.ranked(["HP:0001250"], "omim", 5)
    assert "HTML page" in err.value.message


def test_F10_a_row_with_no_score_is_dropped_not_emitted_as_nan():
    rows = pubcasefinder.parse_tsv(
        "Rank\tScore\tOMIM_ID\tDisease_Name\tMatched_Phenotype\tCausative_Gene\n"
        "1\t0.9\tOMIM:607208\tDravet\tHP:0001250\tSCN1A\n"
        "2\t\tOMIM:111111\tNo score\t\t\n", "omim", 5)
    assert [r["id"] for r in rows] == ["OMIM:607208"]
    json.dumps(rows, allow_nan=False)  # would raise if a NaN had got through
