import json
import math
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
    # D-P1-4 / CP1-7 / D-P2-7 (v0.2): each must rank its intended term FIRST in the real release
    ("头围小", "HP:0040195"),
    ("头围减小", "HP:0040195"),
    ("听力下降", "HP:0000365"),
    ("耳聋", "HP:0000365"),
    ("肌张力低", "HP:0001252"),
    ("肌张力低下", "HP:0001252"),
    ("孩子走路晚", "HP:0031936"),
    ("走路晚了", "HP:0031936"),
    ("发烧抽风", "HP:0002373"),
    ("不会坐", "HP:0025336"),
    ("说话不清楚", "HP:0001350"),
    ("退行", "HP:0002376"),
    ("走路摇晃", "HP:0002317"),
    ("孩子软", "HP:0001252"),
    ("尿有怪味", "HP:0012088"),
    ("CK高", "HP:0003236"),
    ("不会走路", "HP:0002540"),
    ("slow learner", "HP:0001328"),
    ("my son walked late", "HP:0031936"),
    ("fit", "HP:0001250"),
])
def test_f36_live_lay_phrases_rank_the_right_term_first(real_release, phrase, expect):
    live = real_release
    hits = hpo_local.search(live, phrase, limit=3)
    assert hits, f"{phrase!r} found nothing"
    assert hits[0]["id"] == expect, (phrase, [(h["id"], h["label"]) for h in hits])


@pytest.mark.live
def test_d_p1_4_live_generic_chinese_ranking_respects_direction(real_release, monkeypatch):
    """Without the lay table, 头围小 must still not put Increased head circumference first."""
    monkeypatch.setattr(hpo_local, "LAY_TERMS", {})
    hits = [h["id"] for h in hpo_local.search(real_release, "头围小", limit=3)]
    assert hits[0] != "HP:0040194", hits
    assert "HP:0040195" in hits


@pytest.mark.live
def test_d_p1_3_live_aarskog_scott_is_not_penalised_for_cryptorchidism(real_release):
    """OMIM:305400 lists cryptorchidism with counts 6/11, 1/1 and 0/2; the 0/2 made it a curated NOT."""
    live = real_release
    assert "HP:0000028" not in (live.disease_excluded.get("OMIM:305400") or {})
    for method in ("resnik", "lr"):
        res = hpo_local.rank(live, ["HP:0000049", "HP:0000028", "HP:0000316"], top=50, method=method)
        row = next(d for d in res["diseases"] if d["disease"] == "OMIM:305400")
        assert row["contradicted_by_curation"] == []
    both = [(d, t) for d, ex in live.disease_excluded.items() for t in ex if t in live.disease_terms.get(d, {})]
    assert both == [], f"terms both positive and curated NOT in one disease: {both[:5]}"


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


# ---------------------------------------------------------------- v0.2 regressions (review ids)


def _with_rows(tmp_path, monkeypatch, extra_rows):
    """The HPO fixture with extra phenotype.hpoa rows appended."""
    for f in hpo_local.FILES:
        shutil.copy(os.path.join(FIX, f), tmp_path / f)
    hpoa = tmp_path / "phenotype.hpoa"
    hpoa.write_text(hpoa.read_text(encoding="utf-8") + "".join(
        f"{d}\t{name}\t{qual}\t{term}\tPMID:9\tPCS\t\t{freq}\t\t\tP\tHPO:x\n" for d, name, qual, term, freq in extra_rows),
        encoding="utf-8")
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(tmp_path))
    return hpo_local.load(rebuild=True)


def test_d_p1_3_a_zero_case_count_is_not_a_curated_not(tmp_path, monkeypatch):
    """HPOA "0/2" means two reported patients lacked the feature; it is not a curator's NOT."""
    i2 = _with_rows(tmp_path, monkeypatch, [("OMIM:2", "Hypotonia syndrome", "", "HP:0001250", "0/2")])
    assert "HP:0001250" not in (i2.disease_excluded.get("OMIM:2") or {})
    assert "HP:0001250" not in i2.disease_terms["OMIM:2"]
    for method in ("resnik", "lr"):
        rows = {d["disease"]: d for d in hpo_local.rank(i2, ["HP:0001250", "HP:0001252"], top=5, method=method)["diseases"]}
        assert rows["OMIM:2"]["contradicted_by_curation"] == []


def test_d_p1_3_zero_counts_are_pooled_with_the_positive_ones(tmp_path, monkeypatch):
    """Aarskog-Scott and cryptorchidism: 6/11, 1/1 and 0/2 for one pair pool to 7/14, never a NOT."""
    i2 = _with_rows(tmp_path, monkeypatch, [("OMIM:1", "Febrile epilepsy", "", "HP:0001263", "0/2"),
                                            ("OMIM:1", "Febrile epilepsy", "", "HP:0001263", "3/4")])
    # the fixture already has 1/2 for this pair: (1 + 0 + 3) / (2 + 2 + 4)
    assert i2.disease_terms["OMIM:1"]["HP:0001263"] == pytest.approx(0.5)
    assert "OMIM:1" not in i2.disease_excluded
    res = hpo_local.rank(i2, ["HP:0002373", "HP:0001263"], top=5)
    assert res["diseases"][0]["disease"] == "OMIM:1" and res["diseases"][0]["contradicted_by_curation"] == []


def test_d_p1_3_a_not_never_penalises_a_feature_the_disease_also_has(tmp_path, monkeypatch):
    """Two curations disagree (NOT and a positive row for the same pair): the disease keeps its feature."""
    i2 = _with_rows(tmp_path, monkeypatch, [("OMIM:3", "Cataract disorder", "", "HP:0001250", "1/1")])
    assert i2.disease_excluded["OMIM:3"]["HP:0001250"] == "HPOA NOT qualifier"
    for method in ("resnik", "lr"):
        rows = {d["disease"]: d for d in hpo_local.rank(i2, ["HP:0001250", "HP:0000518"], top=5, method=method)["diseases"]}
        assert rows["OMIM:3"]["contradicted_by_curation"] == [] and rows["OMIM:3"]["penalty"] == 0


def test_d_p1_3_the_index_version_was_bumped():
    """A pickle built by 0.1.0 (version 5) holds the 0/n rows as NOT and must be rebuilt."""
    assert hpo_local.INDEX_VERSION >= 6


def test_d_p2_8_not_on_an_ancestor_of_a_present_term_contradicts_it(idx):
    """OMIM:3 is curated NOT Seizure; a patient with febrile seizure has a seizure."""
    rows = {d["disease"]: d for d in hpo_local.rank(idx, ["HP:0002373", "HP:0000518"], top=5)["diseases"]}
    assert rows["OMIM:3"]["contradicted_by_curation"][0]["term"] == "HP:0001250"
    assert rows["OMIM:3"]["penalty"] > 0


def test_d_p2_8_the_excluded_penalty_scales_with_annotation_frequency(idx):
    """ORPHA:4 has Hypotonia at 'very frequent' (0.9); OMIM:2 states no frequency (1.0)."""
    res = hpo_local.rank(idx, ["HP:0001250"], excluded=["HP:0001252"], top=5)
    pen = {d["disease"]: d["penalty"] for d in res["diseases"]}
    assert pen["ORPHA:4"] == pytest.approx(0.9 * pen["OMIM:2"], abs=2e-4)


def test_d_p2_8_exact_matches_break_ties_before_specificity(idx):
    """The first tie-break key: same score, more exact matches first, whatever the specificity."""
    rows = [(1.0, 1.0, 0.0, 0, [], 0, {}), (1.0, 1.0, 0.0, 1, [], 2, {})]
    spec = {0: 9.0, 1: 1.0}
    out = hpo_local._sort_lazy(rows, lambda di: spec[di], idx)
    assert [r[3] for r in out] == [1, 0]


def test_d_p2_8_chou_feng_is_seizure_offline(idx):
    hits = hpo_local.search(idx, "抽风")
    assert hits[0]["id"] == "HP:0001250" and hits[0]["matched_on"] == "lay phrase"


def test_d_p2_1_an_excluded_term_that_contradicts_a_present_one_is_dropped_and_reported(idx):
    res = hpo_local.rank(idx, ["HP:0002373", "HP:0001263"], excluded=["HP:0001250", "HP:0001263", "HP:0001252"])
    assert [e["id"] for e in res["excluded"]] == ["HP:0001252"]
    got = {c["excluded"]: c for c in res["contradictions"]}
    assert set(got) == {"HP:0001250", "HP:0001263"}
    assert got["HP:0001250"]["present"] == "HP:0002373" and "ancestor" in got["HP:0001250"]["why"]
    assert "both present and excluded" in got["HP:0001263"]["why"]


def test_d_p2_7_lay_overcalls_are_corrected():
    assert hpo_local.LAY_TERMS["slow learner"] == ["HP:0001328"]  # Specific learning disability, not ID
    assert hpo_local.LAY_TERMS["不会走路"][0] == "HP:0002540"      # Inability to walk first


def test_d_p1_4_lay_phrases_match_inside_a_longer_query_unless_negated():
    assert hpo_local.lay_lookup("孩子走路晚")[1] == "走路晚"
    assert hpo_local.lay_lookup("走路晚了")[1] == "走路晚"
    assert hpo_local.lay_lookup("我家宝宝头围小")[1] == "头围小"
    assert hpo_local.lay_lookup("my son walked late")[1] == "walked late"
    assert hpo_local.lay_lookup("my 2 year old son walked late")[1] == "walked late"
    for negated in ("没有抽搐", "否认抽搐", "不抽搐", "no febrile convulsions", "never walked late"):
        assert hpo_local.lay_lookup(negated) == ([], None), negated


@pytest.mark.parametrize("phrase", [
    # R-A-P0-2: negated, or part of another finding: a lay key inside them must not fire
    "无明显抽搐", "未出现抽搐", "不伴抽搐", "没有出现过抽搐", "抽搐：无", "抽搐（-）", "听力下降不明显",
    "hearing loss ruled out", "no speech delay", "no speech or language delay",
    "头大小正常", "舌头大", "拳头大小的肿块", "头围大小正常", "额头小", "孩子软骨发育不全", "身体软组织肿块",
    "脊柱退行性变", "神经退行性疾病", "耳聋家族史", "眼睑抽搐", "孩子发育不好走路晚",
])
def test_r_a_p0_2_a_lay_phrase_inside_another_finding_does_not_fire(phrase):
    assert hpo_local.lay_lookup(phrase) == ([], None)


def test_r_a_p0_2_an_official_label_that_is_the_query_beats_a_lay_key_inside_it(tmp_path, monkeypatch):
    i2 = _zh_index(tmp_path, monkeypatch, {
        "HP:0001263": ("Global developmental delay", "全面发育迟缓"),
        "HP:0000750": ("Delayed speech and language development", "语言发育迟缓"),
    })
    assert hpo_local.search(i2, "语言发育迟缓", 2)[0]["id"] == "HP:0000750"
    # even when the lay key inside the query would fire (here forced by widening the filler list),
    # the official label that IS the query comes first
    monkeypatch.setattr(hpo_local, "_ZH_FILLER", hpo_local._ZH_FILLER + ["语言"])
    assert hpo_local.lay_lookup("语言发育迟缓")[1] == "发育迟缓"
    hits = hpo_local.search(i2, "语言发育迟缓", 2)
    assert [h["id"] for h in hits] == ["HP:0000750", "HP:0001263"]


def test_r_a_p1_4_a_negated_direction_word_is_not_a_direction():
    assert hpo_local._polarity("肌张力不高") == 0
    assert hpo_local._polarity("头围不大") == 0
    assert hpo_local._polarity("肌张力高") == 1
    assert hpo_local.lay_lookup("fit and well") == ([], None)  # a single lay word inside a sentence is not used
    assert hpo_local.lay_lookup("ck高")[0] == ["HP:0003236"]


def _zh_index(tmp_path, monkeypatch, labels):
    """A tiny ontology whose terms carry the given Chinese labels: {id: (english, chinese)}."""
    nodes = [{"id": "http://purl.obolibrary.org/obo/HP_0000001", "lbl": "All"},
             {"id": "http://purl.obolibrary.org/obo/HP_0000118", "lbl": "Phenotypic abnormality"}]
    edges = [{"sub": "http://purl.obolibrary.org/obo/HP_0000118", "pred": "is_a",
              "obj": "http://purl.obolibrary.org/obo/HP_0000001"}]
    for tid, (en, _zh) in labels.items():
        iri = "http://purl.obolibrary.org/obo/HP_" + tid.split(":")[1]
        nodes.append({"id": iri, "lbl": en})
        edges.append({"sub": iri, "pred": "is_a", "obj": "http://purl.obolibrary.org/obo/HP_0000118"})
    (tmp_path / "hp.json").write_text(json.dumps({"graphs": [{"nodes": nodes, "edges": edges, "meta": {}}]}), "utf-8")
    (tmp_path / "phenotype.hpoa").write_text(
        "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration\n"
        + "".join(f"OMIM:1\tD\t\t{t}\tPMID:1\tPCS\t\t\t\t\tP\tHPO:x\n" for t in labels), "utf-8")
    (tmp_path / "genes_to_phenotype.txt").write_text("ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n",
                                                     "utf-8")
    (tmp_path / hpo_local.ZH_FILE).write_text(
        "subject_id\tpredicate_id\ttranslation_value\n"
        + "".join(f"{t}\trdfs:label\t{zh}\n" for t, (_en, zh) in labels.items()), "utf-8")
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(tmp_path))
    return hpo_local.load(rebuild=True)


def test_d_p1_4_small_head_circumference_is_not_answered_with_its_opposite(tmp_path, monkeypatch):
    i2 = _zh_index(tmp_path, monkeypatch, {
        "HP:0040194": ("Increased head circumference", "头围增加"),
        "HP:0040195": ("Decreased head circumference", "头围减小"),
        "HP:6000339": ("Decreased fetal abdominal circumference", "胎儿腹围小"),
        "HP:0000252": ("Microcephaly", "小头畸形"),
    })
    assert [h["id"] for h in hpo_local.search(i2, "头围小", 2)] == ["HP:0040195", "HP:0000252"]
    # without the lay table, the generic ranking must still not put the opposite first
    monkeypatch.setattr(hpo_local, "LAY_TERMS", {})
    assert hpo_local.search(i2, "头围小", 3)[0]["id"] == "HP:0040195"


def test_d_p1_4_direction_breaks_an_otherwise_exact_tie(tmp_path, monkeypatch):
    """Same bigram and character overlap, opposite direction: the one that agrees comes first."""
    monkeypatch.setattr(hpo_local, "LAY_TERMS", {})
    i2 = _zh_index(tmp_path, monkeypatch, {
        "HP:0000002": ("Up label", "头围加大"),  # the smaller id: it would win a plain tie
        "HP:0000003": ("Down label", "头围减低"),
    })
    hits = [h["id"] for h in hpo_local.search(i2, "头围小", 3)]
    assert hits.index("HP:0000003") < hits.index("HP:0000002")
    assert hpo_local._polarity("头围增加") == 1 and hpo_local._polarity("头围小") == -1
    assert hpo_local._polarity("小脑萎缩") == 0  # 小 inside 小脑 (cerebellum) is not a direction


def test_d_p1_4_hearing_and_tone_phrases_hit_the_general_term(tmp_path, monkeypatch):
    i2 = _zh_index(tmp_path, monkeypatch, {
        "HP:0007663": ("Reduced visual acuity", "视力下降"),
        "HP:0000365": ("Hearing impairment", "听力受损"),
        "HP:0009900": ("Unilateral deafness", "单侧耳聋"),
        "HP:0001252": ("Hypotonia", "肌张力减退"),
        "HP:0000297": ("Facial hypotonia", "面部肌张力低下"),
        "HP:0002376": ("Developmental regression", "发育倒退"),
        "HP:0005237": ("Degenerative liver disease", "退行性肝病"),
    })
    for phrase, want in (("听力下降", "HP:0000365"), ("耳聋", "HP:0000365"), ("肌张力低下", "HP:0001252"),
                         ("肌张力低", "HP:0001252"), ("退行", "HP:0002376"), ("孩子软", "HP:0001252")):
        assert hpo_local.search(i2, phrase, 3)[0]["id"] == want, phrase


# ---------------------------------------------------------------- the likelihood-ratio method


def test_lr_ranks_the_fully_matching_disease_first(idx):
    res = hpo_local.rank(idx, ["HP:0002373", "HP:0007359", "HP:0001263"], top=4, method="lr")
    assert res["method_id"] == "lr" and res["diseases"][0]["disease"] == "OMIM:1"
    assert res["diseases"][0]["matches"][0]["exact"] is True and "ln_lr" in res["diseases"][0]["matches"][0]


def test_lr_uses_the_annotation_frequency(tmp_path, monkeypatch):
    """Same feature, very frequent in one disease and very rare in the other."""
    i2 = _with_rows(tmp_path, monkeypatch, [("OMIM:5", "Rare-feature disease", "", "HP:0001252", "HP:0040284"),
                                            ("OMIM:6", "Common-feature disease", "", "HP:0001252", "HP:0040281")])
    rows = [d["disease"] for d in hpo_local.rank(i2, ["HP:0001252"], top=10, method="lr")["diseases"]]
    assert rows.index("OMIM:6") < rows.index("OMIM:5")
    # Resnik ignores frequency: the two tie
    res = {d["disease"]: d["score"] for d in hpo_local.rank(i2, ["HP:0001252"], top=10, method="resnik")["diseases"]}
    assert res["OMIM:5"] == res["OMIM:6"]


def test_lr_excluded_term_cost_is_capped(idx):
    res = hpo_local.rank(idx, ["HP:0001250"], excluded=["HP:0001252"], top=5, method="lr")
    pen = {d["disease"]: d["penalty"] for d in res["diseases"]}
    cap = hpo_local.LR_PARAMS["excluded_cap"]
    assert pen["OMIM:2"] == pytest.approx(-math.log(1 - cap), abs=1e-3)  # frequency 1.0 is capped
    assert pen["ORPHA:4"] == pytest.approx(-math.log(1 - min(cap, 0.9)), abs=1e-3)


def test_rank_refuses_an_unknown_method_or_parameter(idx):
    with pytest.raises(ValueError, match="method must be"):
        hpo_local.rank(idx, ["HP:0001250"], method="exomiser")
    with pytest.raises(ValueError, match="unknown ranking parameter"):
        hpo_local.rank(idx, ["HP:0001250"], method="lr", params={"bogus": 1})
