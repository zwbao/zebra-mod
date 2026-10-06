import json

import pytest

from zebra import cli
from zebra.commands import china


def _run(argv, capsys):
    code = cli.main(argv)
    out = capsys.readouterr().out
    return code, out


def test_bundled_lists_are_complete_and_sourced():
    data = china.load()
    prov = data["provenance"]
    assert prov["counts"] == {"1": 121, "2": 86, "total": 207}
    for li, n in ((1, 121), (2, 86)):
        rows = [d for d in data["diseases"] if d["list"] == li]
        assert [d["no"] for d in rows] == list(range(1, n + 1))
        assert all(d["name_zh"] and d["name_en"] for d in rows)
        assert prov["lists"][str(li)]["source_url"].startswith("https://www.gov.cn/")
    assert prov["retrieved_at"]
    # spot checks against the published lists
    by = {(d["list"], d["no"]): d for d in data["diseases"]}
    assert by[(1, 31)]["name_zh"] == "戈谢病"
    assert by[(1, 105)]["name_en"] == "Severe Myoclonic Epilepsy in Infancy (Dravet Syndrome)"
    assert by[(2, 64)]["name_en"] == "Primary growth hormone deficiency"  # gov.cn copy (MIIT copy has a wrong English name)
    assert by[(2, 86)]["name_zh"] == "West综合征/婴儿痉挛综合征"


@pytest.mark.parametrize("query,list_no,kind", [
    ("戈谢病", (1, 31), "exact"),
    ("Gaucher disease", (1, 31), "exact"),
    ("Gaucher’s Disease", (1, 31), "exact"),
    ("GAUCHER'S DISEASE", (1, 31), "exact"),
    ("Dravet syndrome", (1, 105), "part"),
    ("Dravet综合征", (1, 105), "part"),
    ("白塞病", (2, 9), "part"),
    ("贝赫切特综合征", (2, 9), "part"),
    ("Wilson disease", (1, 37), "part"),
    ("CDKL5-deficiency disorder", (2, 11), "exact"),
    ("脊髓性肌萎缩症", (1, 110), "exact"),
])
def test_match_names(query, list_no, kind):
    hits = china.match(query)
    assert hits, query
    assert (hits[0]["list"], hits[0]["no"]) == list_no
    assert hits[0]["match"] == kind


def test_e_9_group_entries_claim_members_only_through_orphanet_classification():
    """E-9: one rule for group entries. The published name alone only overlaps (offered, not claimed);
    the alias layer's Orphanet members (classification under the entry's own ORPHAcode) are claimed, and say so."""
    hits = china.match("Duchenne muscular dystrophy")
    assert hits and hits[0]["no"] == 98 and hits[0]["match"] == "overlap"
    res = china.lookup(["Duchenne muscular dystrophy"])
    assert res["on_list"] is True and res["status"] == "on_list"
    m = res["matches"][0]
    assert m["no"] == 98 and m["match"] == "group_member" and m["alias_orpha"] == "ORPHA:98896"
    assert m["member_of"] == "ORPHA:206644" and "Orphanet classification" in m["alias_source"]
    # Becker is a member too; a partial name is only offered
    assert china.lookup(["Becker muscular dystrophy"])["status"] == "on_list"
    part = china.lookup(["Duchenne"])
    assert part["on_list"] is False and part["status"] == "possible"


def test_not_on_list():
    assert china.lookup(["NGLY1 deficiency"])["on_list"] is False


def test_norm_folds_width_case_and_possessive():
    assert china.norm("Gaucher’s Disease") == china.norm("gaucher disease") == "gaucherdisease"
    assert china.norm("（Ｄｒａｖｅｔ）") == "dravet"


def test_cli_text_and_json(capsys):
    code, out = _run(["china", "戈谢病"], capsys)
    assert code == 0 and "#31" in out and "Gaucher" in out
    code, out = _run(["--json", "china", "Gaucher", "disease"], capsys)
    env = json.loads(out)
    assert env["ok"] and env["command"] == "china"
    assert env["result"]["status"] == "on_list"
    assert env["result"]["matches"][0]["no"] == 31
    assert env["sources"][0]["url"].startswith("https://www.gov.cn/")
    code, out = _run(["--json", "china", "NGLY1 deficiency"], capsys)
    assert json.loads(out)["result"]["status"] == "not_found"


def test_missing_data_file_is_reported(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(china, "DATA_FILE", str(tmp_path / "nope.json"))
    monkeypatch.setattr(china, "_cache", {})
    code, out = _run(["--json", "china", "戈谢病"], capsys)
    env = json.loads(out)
    assert env["result"]["status"] == "unavailable"
    assert any("unavailable" in w for w in env["warnings"])


# ---------------------------------------------------------------- P1f: Chinese and common names

def test_P1f_alias_layer_is_bundled_with_provenance():
    d = china.load_aliases()
    assert d["entries"] and d["schema"] == "zebra.china_aliases/2"
    p = d["provenance"]
    assert p["counts"]["list_entries"] == 207
    # every alias carries a kind and a source; a folk name names the official entry
    kinds = set()
    for e in d["entries"]:
        for a in e["aliases"]:
            assert a.get("text") and a.get("kind") and a.get("source")
            kinds.add(a["kind"])
            if a["kind"] == "folk_name":
                assert a["maps_to"] == e["name_zh"]
                assert "not from a retrieved dataset" in a["source"]
            else:
                assert "China list" in a["source"] or "Orphanet" in a["source"]
    assert {"official_zh", "official_en", "orphanet_zh", "orphanet_synonym_zh", "folk_name"} <= kinds
    # the Orphanet endpoints and the dataset date are recorded
    assert "orphadata.com" in json.dumps(p["sources"])
    # the entries that could not be linked are listed, not guessed at
    assert isinstance(p["unmatched_list_entries"], list) and p["unmatched_list_entries"]


def test_P1f_orphanet_zh_name_index_is_bundled_with_provenance():
    d = china.load_zh_names()
    p = d["provenance"]
    assert d["names"]["98896"] == "杜氏肌营养不良症"
    assert d["names"]["70"] == "近端脊髓性肌萎缩"
    assert p["source_url"].startswith("https://api.orphadata.com/")
    assert p["dataset_date"] == "2020-06-01" and p["licence"].startswith("CC-BY")
    assert p["names_kept"] == len(d["names"]) == 10563


@pytest.mark.parametrize("query,no,kind", [
    ("小胖威利", 93, "folk_name"),
    ("瓷娃娃", 86, "folk_name"),
    ("渐冻症", 4, "folk_name"),
    ("快乐木偶综合征", 5, "folk_name"),
    ("德拉韦综合征", 105, "folk_name"),
    ("普拉德-威利综合征", 93, "folk_name"),
    ("瑞特综合征", 72, "folk_name"),
    ("杜氏肌营养不良", 98, "folk_name"),
    ("脊髓性肌萎缩症", 110, "official_zh"),
    ("成骨不全症", 86, "official_part_zh"),
])
def test_P1f_names_families_use_are_found_on_the_list(query, no, kind):
    hits = china.alias_hits(query, limit=3)
    assert hits, query
    assert hits[0]["no"] == no and hits[0]["match"] in china.ON_LIST_KINDS
    assert hits[0]["alias_kind"] == kind
    assert hits[0]["matched_on"] and hits[0]["alias_source"]
    assert china.lookup([query])["on_list"] is True


def test_P1f_bigram_similarity_is_dice_on_character_bigrams():
    assert china.dice(china.bigrams("成骨不全症"), china.bigrams("成骨不全")) == pytest.approx(0.857, abs=0.001)
    assert china.dice(china.bigrams("杜氏肌营养不良"), china.bigrams("杜氏肌营养不良症")) == pytest.approx(0.923, abs=0.001)
    # the family of one-letter syndrome names must stay below the threshold
    assert china.dice(china.bigrams("德拉韦综合征"), china.bigrams("W综合征")) < china.BIGRAM_MIN
    assert china.dice(china.bigrams(""), china.bigrams("x")) == 0.0


def test_P1f_bigram_matching_is_only_for_chinese_input():
    """On Latin names a shared word is not a near-miss: NGLY1 deficiency ~ T2 deficiency."""
    assert china.dice(china.bigrams("NGLY1 deficiency"), china.bigrams("T2 deficiency")) > china.BIGRAM_MIN
    assert not [h for h in china.alias_hits("NGLY1 deficiency", limit=5) if h["match"] == "alias_bigram"]
    assert china.lookup(["NGLY1 deficiency"])["on_list"] is False


def test_P1f_every_answer_says_what_it_matched(capsys):
    code, out = _run(["china", "瓷娃娃"], capsys)
    assert code == 0
    assert "matched '瓷娃娃'" in out and "folk_name" in out
    assert "maps to the official name 成骨不全症（脆骨病）" in out
    code, out = _run(["--json", "china", "杜氏肌营养不良"], capsys)
    m = json.loads(out)["result"]["matches"][0]
    assert m["matched_on"] == "杜氏肌营养不良" and m["alias_kind"] == "folk_name"


def test_P1f_an_unmatched_name_is_not_called_absent_from_the_list(capsys):
    code, out = _run(["china", "NGLY1 deficiency"], capsys)
    assert "not the same as the disease being absent from the lists" in out
    assert "zebra disease" in out
    assert "not on 第一批" not in out  # the old wording asserted absence


def test_P1f_orphanet_zh_candidates_are_offered_for_a_chinese_query(capsys):
    code, out = _run(["--json", "china", "脊髓性肌萎缩症"], capsys)
    cands = json.loads(out)["result"]["orphanet_zh_candidates"]
    assert cands[0]["orpha"] == "ORPHA:70" and cands[0]["similarity"] >= 0.7
    assert "Orphanet Chinese preferred term" in cands[0]["source"]


def test_P1f_alias_sources_are_recorded_in_the_ledger_rows(capsys):
    code, out = _run(["--json", "china", "瓷娃娃"], capsys)
    dbs = {s["db"] for s in json.loads(out)["sources"]}
    assert "China rare disease list alias layer" in dbs
    assert "China national rare disease list" in dbs


def test_P1f_a_missing_alias_file_falls_back_to_the_published_names(monkeypatch, capsys):
    monkeypatch.setattr(china, "ALIAS_FILE", "/nonexistent/aliases.json")
    monkeypatch.setattr(china, "_cache", {})
    code, out = _run(["--json", "china", "戈谢病"], capsys)
    r = json.loads(out)["result"]
    assert code == 0 and r["status"] == "on_list" and r["matches"][0]["no"] == 31


@pytest.mark.live
def test_live_P1f_chinese_names_resolve_to_the_right_disease_with_the_right_genes(capsys):
    """The acceptance for P1f: 杜氏肌营养不良 -> DMD, 脊髓性肌萎缩症 -> SMN1.

    E-2/CP1-7: 脊髓性肌萎缩症 is the exact name of list entry #110, which has no ORPHAcode; it used to stop at
    ORPHA:70 by a 0.769 bigram match (and lost GeneReviews); it now follows the entry's English name, like the
    English query does."""
    for query, name, orpha, gene in [
        ("杜氏肌营养不良", "Duchenne muscular dystrophy", "ORPHA:98896", "DMD"),
        ("脊髓性肌萎缩症", "spinal muscular atrophy", None, "SMN1"),
    ]:
        code = cli.main(["--json", "disease", query])
        env = json.loads(capsys.readouterr().out)
        assert code == 0, query
        r = env["result"]
        assert r["status"] == "resolved" and r["name"].lower() == name.lower()
        if orpha:
            assert r["ids"]["ORPHA"] == orpha
            # the fuzzy match is always disclosed
            assert any("character-bigram similarity" in n for n in r["notes"]) or \
                any("character-bigram similarity" in w for w in env["warnings"])
        assert gene in {g["symbol"] for g in r["genes"]}, (query, r["genes"])


@pytest.mark.live
def test_live_P1f_folk_names_resolve_to_a_clinical_entity(capsys):
    for query, orpha, gene in [("瓷娃娃", "ORPHA:666", "COL1A1"), ("渐冻症", "ORPHA:803", "SOD1"),
                               ("德拉韦综合征", "ORPHA:33069", "SCN1A"), ("瑞特综合征", "ORPHA:778", "MECP2")]:
        code = cli.main(["--json", "disease", query])
        r = json.loads(capsys.readouterr().out)["result"]
        assert code == 0 and r["ids"]["ORPHA"] == orpha, (query, r.get("ids"))
        assert gene in {g["symbol"] for g in r["genes"]}, (query, sorted({g["symbol"] for g in r["genes"]}))


# ---------------------------------------------------------------- E-1: qualifiers, acronyms, subtypes (P0)

@pytest.mark.parametrize("query,no", [
    ("帕金森病", (1, 87)), ("Parkinson's disease", (1, 87)), ("Parkinson disease", (1, 87)),
    ("地中海贫血", (2, 78)), ("glycogen storage disease", (1, 35)), ("糖原累积病", (1, 35)), ("糖原贮积症", (1, 35)),
])
def test_e_1_a_wider_name_than_a_qualified_entry_is_not_on_the_list(query, no):
    """#87 is young/early-onset PD only, #78 thalassaemia major only, #35 GSD types I and II only."""
    res = china.lookup([query])
    assert res["on_list"] is False and res["status"] == "qualified", (query, res)
    m = res["matches"][0]
    assert (m["list"], m["no"]) == no and m["match"] == "wider"
    assert res["qualifier"][0]["covers"]


def test_e_1_qualified_cli_says_which_subtype_and_warns(capsys):
    code, out = _run(["--json", "china", "地中海贫血"], capsys)
    env = json.loads(out)
    assert env["result"]["status"] == "qualified"
    assert "thalassaemia major only" in env["result"]["qualifier"][0]["covers"]
    assert any("wider than" in w and "do not say the disease is on the list" in w for w in env["warnings"])
    code, out = _run(["china", "帕金森病"], capsys)
    assert "only a subtype is on the national list" in out and "young-onset / early-onset" in out
    assert "on the national list\n" not in out


@pytest.mark.parametrize("query", ["重型", "冷吡啉", "冷炎素", "I型、Ⅱ型", "青年型、早发型"])
def test_e_1_a_bracket_qualifier_or_gloss_is_not_a_disease_name(query):
    res = china.lookup([query])
    assert res["on_list"] is False and res["status"] != "on_list"
    assert all(m["match"] not in ("exact", "part") for m in res["matches"])


@pytest.mark.parametrize("query", ["CAD", "HSP", "PV", "MFS", "GSD", "AE", "ALD", "FD", "MEN", "DMD", "SMA", "CAPS"])
def test_e_1_an_acronym_is_never_on_the_list(query, capsys):
    res = china.lookup([query])
    assert res["on_list"] is False and res["status"] == "possible", (query, res["status"])
    assert res["matches"][0]["match"] == "acronym"
    code, out = _run(["--json", "china", query], capsys)
    env = json.loads(out)
    assert env["result"]["status"] == "possible"
    assert any("acronym" in w for w in env["warnings"])


@pytest.mark.parametrize("query,no,orpha", [
    ("Pompe disease", (1, 35), "ORPHA:365"), ("庞贝病", (1, 35), "ORPHA:365"),
    ("Glycogen storage disease type II", (1, 35), "ORPHA:365"), ("Von Gierke disease", (1, 35), "ORPHA:364"),
    ("β-地中海贫血重型", (2, 78), "ORPHA:231214"), ("Beta-thalassemia major", (2, 78), "ORPHA:231214"),
    ("Young-onset Parkinson disease", (1, 87), "ORPHA:2828"), ("青年发病型帕金森病", (1, 87), "ORPHA:2828"),
])
def test_e_1_cp1_7_a_subtype_inside_the_qualifier_is_on_the_list_with_its_own_orpha(query, no, orpha):
    res = china.lookup([query])
    assert res["on_list"] is True and res["status"] == "on_list", query
    m = res["matches"][0]
    assert (m["list"], m["no"]) == no and m["match"] == "subtype" and m["alias_orpha"] == orpha


def test_e_1_danon_disease_is_not_taken_as_gsd_type_ii():
    """'GSD IIb' is a historical label of Danon disease (LAMP2), not Pompe disease."""
    res = china.lookup(["Danon disease"])
    assert res["status"] != "on_list"


def test_e_1_lookup_by_orphacode_first():
    assert china.lookup([], orpha_ids=["ORPHA:365"])["matches"][0]["match"] == "subtype"
    gsd = china.lookup([], orpha_ids=["ORPHA:79201"])
    assert gsd["status"] == "qualified" and gsd["matches"][0]["no"] == 35
    dmd = china.lookup([], orpha_ids=["ORPHA:98896"])
    assert dmd["status"] == "on_list" and dmd["matches"][0]["match"] == "group_member"
    assert china.lookup([], orpha_ids=["ORPHA:33069"])["matches"][0]["match"] == "orpha_exact"
    assert china.lookup([], orpha_ids=["ORPHA:999999999"])["status"] == "not_found"


@pytest.mark.parametrize("query", ["癫痫", "贫血", "糖尿病", "白内障", "心肌病"])
def test_e_13_a_generic_term_stays_possible(query, monkeypatch):
    """E-13: kills the mutant that adds `contains` / `alias_contains` to the on-list kinds."""
    res = china.lookup([query])
    assert res["on_list"] is False and res["status"] in ("possible", "not_found"), (query, res["status"])


def test_e_13_mutant_contains_on_list_is_caught(monkeypatch):
    monkeypatch.setattr(china, "ON_LIST_KINDS", china.ON_LIST_KINDS + ("contains", "alias_contains"))
    assert china.lookup(["癫痫"])["status"] == "on_list"  # what the guard above protects against


def test_session_probe_12_a_containment_is_not_worded_as_a_published_name_match(capsys):
    code, out = _run(["china", "Dravet"], capsys)
    assert "no entry matched this name exactly" in out
    assert "matched the published name" not in out


# ---------------------------------------------------------------- E-10: the alias layer is rebuildable

def test_e_10_the_alias_builder_is_in_the_repo_and_named_by_the_data():
    import os

    prov = china.load_aliases()["provenance"]
    assert prov["builder"].startswith("tools/china_list/build_china_aliases.py")
    assert "/private/tmp" not in json.dumps(prov) and "scratchpad" not in json.dumps(prov)
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    assert os.path.exists(os.path.join(root, "tools", "china_list", "build_china_aliases.py"))
    # the rebuild reproduces the lost builder's links and aliases; the only differences are the
    # deliberate ones: the obsolete ORPHA:73274 link dropped, Orphanet's "NON RARE IN EUROPE:" label
    # removed from an alias, the ambiguous folk name 玻璃人 dropped, and the Pompe folk names added
    # ... and one new exact link: British spellings now fold, so "Hyperornithinaemia-…" equals Orphanet's
    # "Hyperornithinemia-…" (ORPHA:415)
    chk = prov["rebuild_check"]
    assert chk["links_differing"] == ["1/48: None -> ORPHA:415", "2/2: ORPHA:73274 -> None"]
    assert set(chk["aliases_removed"]) == {"1/36 玻璃人 [folk_name]",
                                          "1/76 NON RARE IN EUROPE: Multiple sclerosis [orphanet_en]",
                                          "2/2 OBSOLETE: Acquired hemophilia [orphanet_en]", "2/2 获得性血友病 [orphanet_zh]"}
    added = set(chk["aliases_added"])
    assert {"1/35 庞贝病 [folk_name]", "1/35 庞贝氏症 [folk_name]", "1/76 Multiple sclerosis [orphanet_en]"} <= added
    assert all(a.startswith(("1/35 ", "1/76 ", "1/48 ")) for a in added)
    assert prov["counts"]["linked_to_orphanet"] == 149


def test_e_10_builder_folding_matches_the_command():
    import importlib.util
    import os

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    import sys

    spec = importlib.util.spec_from_file_location("bca", os.path.join(root, "tools", "china_list",
                                                                      "build_china_aliases.py"))
    bca = importlib.util.module_from_spec(spec)
    old, sys.dont_write_bytecode = sys.dont_write_bytecode, True  # no __pycache__ inside tools/
    try:
        spec.loader.exec_module(bca)
    finally:
        sys.dont_write_bytecode = old
    for t in ("Gaucher’s Disease", "糖原累积病（I型、Ⅱ型）", "West syndrome/Infantile spasms syndrome", "ＤＲＡＶＥＴ"):
        assert bca.norm(t) == china.norm(t)
    for t in ("CAD", "aHUS", "X-ALD", "GEP-NEN", "SMA"):
        assert bca.is_acronym(t) and china.is_acronym(t)
    for t in ("Dravet", "Pompe disease", "戈谢病"):
        assert not bca.is_acronym(t) and not china.is_acronym(t)
    assert set(bca.QUALIFIED) == {(1, 35), (1, 87), (2, 78)}


def test_e_1_every_bracketed_part_of_a_published_name_has_a_reviewed_role():
    """A new list version with a new bracket must be reviewed: a part is a kept alias or a dropped fragment."""
    import re
    import unicodedata

    al = china.load_aliases()
    dropped = " ".join(al["provenance"]["notes"])
    for e in china.load()["diseases"]:
        for field in ("name_zh", "name_en"):
            name = unicodedata.normalize("NFKC", e[field])
            for p in [x.strip() for x in re.split(r"[/()\[\]（）]", name) if x.strip()]:
                if p == name.strip():
                    continue
                kept = (e["list"], e["no"], china.norm(p)) in al["_parts"]
                assert kept or f"{e['list']}/{e['no']} '{p}'" in dropped, (e["list"], e["no"], p)


# ---------------------------------------------------------------- CP0-4: 协作网 hospitals

def test_cp0_4_network_hospitals_are_bundled_with_provenance():
    d = china.load_hospitals()
    p = d["provenance"]
    assert p["counts"]["total"] == 419 and p["counts"]["national_lead"] == 1
    assert p["counts"]["provincial_lead"] == 32 and p["counts"]["member"] == 386
    assert p["counts"]["list_2019"]["total"] == 324
    for sid in ("nhc2019_157", "nhc2024_64"):
        s = p["sources"][sid]
        assert s["url"] and s["retrieved_at"] and len(s["sha256"]) == 64 and s["rows"] in (324, 419)
    assert "国卫办医政函〔2024〕64号" == p["sources"]["nhc2024_64"]["document_no"]
    assert p["builder"].startswith("tools/china_access/build_hospitals.py")
    cur = [h for h in d["hospitals"] if h["source_id"] == "nhc2024_64"]
    assert len(cur) == 419
    assert next(h for h in cur if h["role"] == "national_lead")["name"] == "中国医学科学院北京协和医院"


def test_cp0_4_china_hospitals_by_province(capsys):
    code, out = _run(["--json", "china", "hospitals", "--province", "浙江省"], capsys)
    env = json.loads(out)
    r = env["result"]
    assert code == 0 and r["status"] == "ok" and r["count"] == 15
    assert [h["role"] for h in r["hospitals"]].count("provincial_lead") == 1
    assert all(h["province"] == "浙江" for h in r["hospitals"])
    assert env["sources"] and env["sources"][0]["url"].startswith("https://")
    code, out = _run(["--json", "china", "hospitals", "浙江"], capsys)
    assert json.loads(out)["result"]["count"] == 15
    # without a province: the national and the 32 provincial lead hospitals
    code, out = _run(["--json", "china", "hospitals"], capsys)
    r = json.loads(out)["result"]
    assert r["count"] == 33 and r["not_in_current_list"] == 28
    # a 2019 name absent from the 2024 list is never offered as a network hospital
    assert all(h.get("status") != "removed_2024" for h in r["hospitals"])


def test_cp0_4_an_unknown_province_is_named_not_answered_empty(capsys):
    code, out = _run(["--json", "china", "hospitals", "--province", "火星"], capsys)
    env = json.loads(out)
    assert env["result"]["status"] == "province_not_found" and "浙江" in env["result"]["provinces"]
    assert any("no 协作网 hospital is listed under '火星'" in w for w in env["warnings"])


def test_cp0_4_a_missing_hospital_file_is_not_bundled_not_absent(monkeypatch, capsys):
    monkeypatch.setattr(china, "HOSPITALS_FILE", "/nonexistent/h.json")
    monkeypatch.setattr(china, "_cache", {})
    code, out = _run(["--json", "china", "hospitals"], capsys)
    env = json.loads(out)
    assert env["result"]["status"] == "unavailable" and env["result"]["hospitals"] == []
    assert any("not evidence that a hospital is or is not a member" in w for w in env["warnings"])


def test_province_flag_goes_with_hospitals_only(capsys):
    code, out = _run(["--json", "china", "戈谢病", "--province", "浙江"], capsys)
    assert code == 2 and json.loads(out)["error"]["type"] == "UsageError"


# ---------------------------------------------------------------- adversarial review (W4 round), China list

@pytest.mark.parametrize("query", ["glycogen storage disease type III", "Glycogen Storage Disease Type III",
                                   "糖原累积病III型", "Cori disease"])
def test_rev_p0_1_gsd_type_iii_is_not_folded_into_types_i_and_ii(query):
    assert china.norm("Glycogen Storage Disease (Type I、II）") != china.norm("glycogen storage disease type III")
    assert china.lookup([query])["on_list"] is False


@pytest.mark.parametrize("query,entry", [("Gastroschisis", (2, 73)), ("Whipple disease", (2, 66)),
                                         ("Herpes simplex virus encephalitis", (2, 66)),
                                         ("22q11.2 deletion syndrome", (2, 34)), ("Dent disease", (1, 51)),
                                         ("Acute disseminated encephalomyelitis", (1, 9)),
                                         ("Congenital fibrosis of extraocular muscles", (1, 98))])
def test_rev_p0_2_a_disease_orphanet_only_files_under_a_group_is_not_on_the_list(query, entry):
    res = china.lookup([query])
    assert res["on_list"] is False, (query, res["matches"][:1])
    hit = next((m for m in res["matches"] if (m["list"], m["no"]) == entry), None)
    assert hit is None or hit["match"] == "group_filed"


def test_rev_p0_2_named_members_and_orpha_lookup_still_count():
    assert china.lookup(["Severe combined immunodeficiency"])["status"] == "on_list"
    assert china.lookup(["Hemophilia A"])["status"] == "on_list"
    gas = china.lookup([], orpha_ids=["ORPHA:2368"])  # Gastroschisis
    assert gas["on_list"] is False and gas["matches"][0]["match"] == "group_filed"


@pytest.mark.parametrize("query,status", [("MAD deficiency", "possible"), ("G6P deficiency", "possible"),
                                          ("thalassaemia major", "on_list"), ("parkinsons disease", "qualified"),
                                          ("Pompe病", "on_list"), ("Hyperornithinemia-hyperammonemia-homocitrullinuria syndrome",
                                                                  "on_list"),
                                          ("Methylmalonic acidemia", "on_list"),
                                          ("Mitochondrial encephalomyopathy", "on_list"),
                                          ("糖原累积病I型", "on_list"), ("重型地中海贫血", "on_list"),
                                          ("玻璃人", "not_found")])
def test_rev_p1_acronym_phrases_spelling_and_derived_subtypes(query, status):
    assert china.lookup([query])["status"] == status


def test_rev_p2_obsolete_and_labels_are_not_aliases():
    al = china.load_aliases()["_by_key"]
    assert al[(2, 2)]["orpha"] is None  # ORPHA:73274 is "OBSOLETE: Acquired hemophilia"
    texts = [a["text"] for e in china.load_aliases()["entries"] for a in e["aliases"]]
    assert not any(t.upper().startswith(("NON RARE IN EUROPE", "OBSOLETE")) for t in texts)
    import re as _re
    assert not [a for e in china.load_aliases()["entries"] for a in e["aliases"]
                if a["kind"] == "orphanet_member_zh" and not _re.search(r"[㐀-鿿]", a["text"])]


def test_rev_p2_limit_and_bingtuan(capsys):
    code, out = _run(["--json", "china", "Pompe disease", "--limit", "0"], capsys)
    assert code == 2 and json.loads(out)["error"]["type"] == "UsageError"
    assert china.province_key("兵团") == china.province_key("新疆生产建设兵团")
    assert china.hospitals("兵团").result["count"] == 11
