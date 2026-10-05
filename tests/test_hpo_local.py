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


@pytest.fixture
def real_release(monkeypatch):
    """The HPO release installed on this machine, for the live-marked F36 checks.

    conftest's autouse fixture redirects ZEBRA_CACHE_DIR to a tmp dir; these tests need
    the real files, so the override is removed here rather than in conftest.py.
    """
    monkeypatch.delenv("ZEBRA_CACHE_DIR", raising=False)
    monkeypatch.delenv("ZEBRA_HPO_DIR", raising=False)
    try:
        return hpo_local.load()
    except hpo_local.HpoDataMissing as err:
        pytest.skip(f"{err} (run: zebra hpo fetch)")


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


# --------------------------------------------------------------------- regressions
# One test per backlog id in docs/ROADMAP.md.


# ---- F29: the excluded-term penalty was the max, not the sum; NOT annotations were dropped

def test_f29_excluded_penalty_is_the_sum_over_distinct_terms(idx):
    """Two excluded terms must cost more than either on its own.

    The docstring always said "each excluded term ... costs its IC times the annotation
    frequency"; the code kept only the largest. OMIM:2 is annotated with both HP:0001252
    and HP:0001263, so excluding both must cost the sum.
    """
    one = hpo_local.rank(idx, ["HP:0001250"], excluded=["HP:0001252"], top=5)
    two = hpo_local.rank(idx, ["HP:0001250"], excluded=["HP:0001252", "HP:0001263"], top=5)
    p1 = {d["disease"]: d["penalty"] for d in one["diseases"]}["OMIM:2"]
    p2 = {d["disease"]: d["penalty"] for d in two["diseases"]}["OMIM:2"]
    assert p1 > 0
    assert p2 > p1, "a second contradicted feature must add to the penalty, not be discarded"
    single_other = hpo_local.rank(idx, ["HP:0001250"], excluded=["HP:0001263"], top=5)
    p_other = {d["disease"]: d["penalty"] for d in single_other["diseases"]}["OMIM:2"]
    # the reported penalties are rounded to 4 dp, hence the tolerance
    assert p2 == pytest.approx(p1 + p_other, abs=2e-4)


def test_f29_penalty_is_bounded_by_the_excluded_terms_own_ic(idx):
    """The sum needs no artificial ceiling: it cannot exceed sum(IC) / len(query)."""
    excl = ["HP:0001252", "HP:0001263"]
    res = hpo_local.rank(idx, ["HP:0001250"], excluded=excl, top=5)
    bound = sum(idx.ic.get(e, 0.0) for e in excl) / 1
    for d in res["diseases"]:
        assert d["penalty"] <= bound + 1e-4  # penalties are reported rounded to 4 dp


def test_f29_method_string_says_sum_not_max(idx):
    m = hpo_local.rank(idx, ["HP:0001250"], top=2)["method"]
    assert "SUM" in m and "NOT" in m


def test_f29_hpoa_not_annotations_are_kept(idx):
    """phenotype.hpoa row: OMIM:3 NOT HP:0001250. It must be recorded, not discarded."""
    assert idx.disease_excluded["OMIM:3"]["HP:0001250"] == "HPOA NOT qualifier"
    # and still absent from the positive annotations
    assert "HP:0001250" not in idx.disease_terms["OMIM:3"]


def test_f29_not_annotation_penalises_a_contradicted_candidate(idx):
    """A disease curated NOT to have the patient's present feature is pushed down."""
    res = hpo_local.rank(idx, ["HP:0001250", "HP:0000518"], top=5)
    rows = {d["disease"]: d for d in res["diseases"]}
    omim3 = rows["OMIM:3"]
    assert omim3["penalty"] > 0
    assert omim3["contradicted_by_curation"] == [
        {"term": "HP:0001250", "label": "Seizure", "source": "HPOA NOT qualifier"}
    ]
    # the contradiction costs OMIM:3 exactly IC(Seizure) / len(query)
    assert omim3["penalty"] == pytest.approx(idx.ic["HP:0001250"] / 2, abs=1e-4)
    # and it ranks below where it would have without the NOT row
    no_not = dict(idx.disease_excluded)
    no_not.pop("OMIM:3")
    idx.disease_excluded = no_not
    unpenalised = {d["disease"]: d for d in
                   hpo_local.rank(idx, ["HP:0001250", "HP:0000518"], top=5)["diseases"]}["OMIM:3"]
    assert unpenalised["penalty"] == 0.0 and unpenalised["score"] > omim3["score"]


def test_f29_excluded_frequency_class_is_kept_too(tmp_path, monkeypatch):
    """Frequency class HP:0040285 ("Excluded", 0%) is the other NOT channel."""
    import shutil

    for f in hpo_local.FILES:
        shutil.copy(os.path.join(FIX, f), tmp_path / f)
    hpoa = tmp_path / "phenotype.hpoa"
    hpoa.write_text(hpoa.read_text(encoding="utf-8")
                    + "OMIM:2\tHypotonia syndrome\t\tHP:0002373\tPMID:1\tPCS\t\tHP:0040285\t\t\tP\tHPO:x\n",
                    encoding="utf-8")
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(tmp_path))
    i2 = hpo_local.load(rebuild=True)
    assert "0%" in i2.disease_excluded["OMIM:2"]["HP:0002373"]
    assert "HP:0002373" not in i2.disease_terms["OMIM:2"]


# ---- F35: large ties broken by id order, with no way for the caller to know

def test_f35_ties_are_reported(idx):
    res = hpo_local.rank(idx, ["HP:0001250", "HP:0001263"], top=4)
    assert "ties" in res and isinstance(res["ties"], int) and res["ties"] >= 1
    assert "exact query-term matches" in res["tie_break"]
    for d in res["diseases"]:
        assert d["tied_at_this_score"] >= 1


def test_f35_tie_is_broken_by_exact_matches_then_specificity(idx):
    """OMIM:2 and ORPHA:4 score the same here; the tie-break must be principled."""
    res = hpo_local.rank(idx, ["HP:0001250", "HP:0001263"], top=4)
    rows = [d for d in res["diseases"] if d["disease"] in ("OMIM:2", "ORPHA:4")]
    assert len(rows) == 2
    assert rows[0]["score"] == rows[1]["score"], "this fixture pair is a genuine tie"
    assert res["ties"] >= 1
    first, second = rows
    key_first = (first["exact_matches"], first["annotation_specificity"])
    key_second = (second["exact_matches"], second["annotation_specificity"])
    assert key_first >= key_second, "the better-matched/more specific disease must come first"


def test_f35_ordering_is_deterministic(idx):
    a = [d["disease"] for d in hpo_local.rank(idx, ["HP:0001250", "HP:0001263"], top=4)["diseases"]]
    b = [d["disease"] for d in hpo_local.rank(idx, ["HP:0001250", "HP:0001263"], top=4)["diseases"]]
    assert a == b


def test_f35_exact_match_count_is_reported(idx):
    res = hpo_local.rank(idx, ["HP:0002373", "HP:0007359", "HP:0001263"], top=4)
    assert res["diseases"][0]["disease"] == "OMIM:1"
    assert res["diseases"][0]["exact_matches"] == 3


# ---- F36: HPO search was literal, so lay phrases found nothing or the wrong term

def test_f36_lay_terms_point_at_real_ids_in_the_fixture(idx):
    """Every lay mapping that the fixture ontology contains must resolve to a real term."""
    seen = 0
    for phrase, ids in hpo_local.LAY_TERMS.items():
        for tid in ids:
            if tid in idx.names:
                seen += 1
                assert idx.primary(tid)[0] == tid, f"{phrase} -> {tid} is not a primary id"
    assert seen >= 3, "the fixture should cover at least a few of the mapped terms"


def test_f36_lay_phrase_beats_a_literal_match(idx):
    hits = hpo_local.search(idx, "fits")
    assert hits[0]["id"] == "HP:0001250" and hits[0]["matched_on"] == "lay phrase"


def test_f36_matched_on_is_reported(idx):
    by_label = hpo_local.search(idx, "seizure")
    assert by_label[0]["id"] == "HP:0001250" and by_label[0]["matched_on"] == "label"
    assert all("matched_on" in h for h in by_label)
    syn = [h for h in hpo_local.search(idx, "seizure", limit=10) if h["matched_on"] == "synonym"]
    assert all(h["id"] != "HP:0001250" for h in syn)


def test_f36_lay_hits_are_deduplicated_and_ordered(idx):
    hits = hpo_local.search(idx, "convulsions")
    ids = [h["id"] for h in hits]
    assert ids[0] == "HP:0001250"
    assert len(ids) == len(set(ids))


def test_f36_a_lay_phrase_with_no_term_in_this_release_degrades_quietly(idx):
    # the fixture has no Short stature; the mapping must simply yield nothing, not raise
    assert all(h["id"] != "HP:0004322" for h in hpo_local.search(idx, "不长个"))


def test_f36_english_ranking_prefers_the_shorter_label(idx):
    """A query matched inside a much longer label is usually the wrong term."""
    hits = hpo_local.search(idx, "seizure", limit=10)
    first = hits[0]
    assert first["id"] == "HP:0001250" and first["label"] == "Seizure"
    lengths = [len(h["label"]) for h in hits]
    assert lengths[0] == min(lengths)


@pytest.mark.live
def test_f36_lay_terms_point_at_real_ids_in_the_release(real_release):
    """Every id in LAY_TERMS must exist in the installed HPO release, as a primary id."""
    live = real_release
    bad = []
    for phrase, ids in hpo_local.LAY_TERMS.items():
        assert ids, phrase
        for tid in ids:
            pid, note = live.primary(tid)
            if pid != tid:
                bad.append((phrase, tid, note))
    assert not bad, f"lay mappings no longer resolve: {bad}"


@pytest.mark.live
@pytest.mark.parametrize("phrase,expect", [
    ("走路晚", "HP:0031936"),
    ("抽风", "HP:0001250"),
    ("抽搐", "HP:0001250"),
    ("不会说话", "HP:0001344"),
    ("智力低下", "HP:0001249"),
    ("发热惊厥", "HP:0002373"),
    ("发烧抽筋", "HP:0002373"),
    ("翻身晚", "HP:0032989"),
    ("不长个", "HP:0004322"),
    ("floppy baby", "HP:0008947"),
    ("late walker", "HP:0031936"),
    ("not talking", "HP:0001344"),
    ("small head", "HP:0000252"),
])
def test_f36_live_lay_phrases_rank_the_right_term_first(real_release, phrase, expect):
    live = real_release
    hits = hpo_local.search(live, phrase, limit=3)
    assert hits, f"{phrase!r} found nothing"
    assert hits[0]["id"] == expect, (phrase, [(h["id"], h["label"]) for h in hits])


@pytest.mark.live
def test_f36_live_ranking_no_longer_promotes_the_wrong_zh_term(real_release):
    """智力低下 used to rank 面部肌张力低下 first; 发热惊厥 ranked 假性惊厥 third."""
    live = real_release
    top3 = [h["id"] for h in hpo_local.search(live, "智力低下", limit=3)]
    assert top3[0] == "HP:0001249"
    febrile = [h["id"] for h in hpo_local.search(live, "发热惊厥", limit=3)]
    assert febrile[0] == "HP:0002373"
    assert "HP:0033053" not in febrile, "Pseudoseizure must not be in the top 3"
    # 抽搐 used to rank Self hugging (a 10-character label) above Tetany (4 characters)
    chou = [h["id"] for h in hpo_local.search(live, "抽搐", limit=4)]
    assert chou[0] == "HP:0001250"
    assert chou.index("HP:0001281") < (chou.index("HP:0032521") if "HP:0032521" in chou else 99)
