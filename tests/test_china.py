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


def test_group_entries_are_offered_not_claimed():
    hits = china.match("Duchenne muscular dystrophy")
    assert hits and hits[0]["no"] == 98 and hits[0]["match"] == "overlap"
    res = china.lookup(["Duchenne muscular dystrophy"])
    assert res["on_list"] is False


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
    assert d["entries"] and d["schema"] == "zebra.china_aliases/1"
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
    ("SMA", 110, "folk_name"),
    ("DMD", 98, "folk_name"),
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
    """The acceptance for P1f: 杜氏肌营养不良 -> DMD, 脊髓性肌萎缩症 -> SMN1."""
    for query, name, orpha, gene in [
        ("杜氏肌营养不良", "Duchenne muscular dystrophy", "ORPHA:98896", "DMD"),
        ("脊髓性肌萎缩症", "Proximal spinal muscular atrophy", "ORPHA:70", "SMN1"),
    ]:
        code = cli.main(["--json", "disease", query])
        env = json.loads(capsys.readouterr().out)
        assert code == 0, query
        r = env["result"]
        assert r["status"] == "resolved" and r["name"] == name
        assert r["ids"]["ORPHA"] == orpha
        assert gene in {g["symbol"] for g in r["genes"]}, (query, r["genes"])
        # the fuzzy match is always disclosed
        assert any("character-bigram similarity" in n for n in r["notes"]) or \
               any("character-bigram similarity" in w for w in env["warnings"])


@pytest.mark.live
def test_live_P1f_folk_names_resolve_to_a_clinical_entity(capsys):
    for query, orpha, gene in [("瓷娃娃", "ORPHA:666", "COL1A1"), ("渐冻症", "ORPHA:803", "SOD1"),
                               ("德拉韦综合征", "ORPHA:33069", "SCN1A"), ("瑞特综合征", "ORPHA:778", "MECP2")]:
        code = cli.main(["--json", "disease", query])
        r = json.loads(capsys.readouterr().out)["result"]
        assert code == 0 and r["ids"]["ORPHA"] == orpha, (query, r.get("ids"))
        assert gene in {g["symbol"] for g in r["genes"]}, (query, sorted({g["symbol"] for g in r["genes"]}))
