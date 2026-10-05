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
