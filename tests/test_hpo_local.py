import os
import shutil

import pytest

from zebra import hpo_local

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "hpo")


@pytest.fixture
def idx(tmp_path, monkeypatch):
    for f in hpo_local.FILES:
        shutil.copy(os.path.join(FIX, f), tmp_path / f)
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(tmp_path))
    return hpo_local.load(rebuild=True)


def test_parse(idx):
    assert idx.version == "fixture-1"
    assert idx.names["HP:0001250"] == "Seizure"
    assert idx.obsolete["HP:0009999"] == "HP:0001250"
    assert "HP:0000118" in idx.ancestors("HP:0002373")
    assert len(idx.disease_list) == 4
    assert idx.disease_genes["OMIM:1"] == ["SCN1A"]


def test_not_qualifier_is_ignored(idx):
    assert "HP:0001250" not in idx.disease_terms["OMIM:3"]


def test_ic_is_monotone(idx):
    assert idx.ic["HP:0002373"] >= idx.ic["HP:0001250"] >= idx.ic.get("HP:0000707", 0)
    assert idx.ic.get("HP:0000118", 0) == pytest.approx(0.0)


def test_rank_puts_best_match_first(idx):
    res = hpo_local.rank(idx, ["HP:0002373", "HP:0007359", "HP:0001263"], top=4)
    assert res["diseases"][0]["disease"] == "OMIM:1"
    assert res["diseases"][0]["relative"] == pytest.approx(1.0)
    assert all(m["exact"] for m in res["diseases"][0]["matches"])
    assert res["genes"][0]["gene"] == "SCN1A"


def test_excluded_terms_penalise(idx):
    plain = hpo_local.rank(idx, ["HP:0001250", "HP:0001263"], top=4)
    excl = hpo_local.rank(idx, ["HP:0001250", "HP:0001263"], excluded=["HP:0001252"], top=4)
    score = {d["disease"]: d["score"] for d in plain["diseases"]}
    score_x = {d["disease"]: d["score"] for d in excl["diseases"]}
    assert score_x["OMIM:2"] < score["OMIM:2"]
    assert score_x["OMIM:1"] == score["OMIM:1"]


def test_obsolete_and_non_phenotype_terms(idx):
    res = hpo_local.rank(idx, ["HP:0009999", "HP:0003593"], top=2)
    assert [q["id"] for q in res["query"]] == ["HP:0001250"]
    assert any("obsolete" in n for n in res["notes"]) and any("not a phenotypic abnormality" in n for n in res["notes"])


def test_search(idx):
    hits = hpo_local.search(idx, "seizures")
    assert hits[0]["id"] == "HP:0001250"
    assert hpo_local.search(idx, "febrile")[0]["id"] == "HP:0002373"
