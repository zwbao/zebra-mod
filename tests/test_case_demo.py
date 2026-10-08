"""zebra case demo: the bundled synthetic case a person can try zebra-mod on."""

import json
import os
import subprocess
import sys

from zebra import case as C

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def zebra(*args, home=None):
    env = dict(os.environ)
    if home:
        env["HOME"] = str(home)
    out = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "zebra"), "--json", *args], capture_output=True, text=True, env=env)
    return json.loads(out.stdout), out.returncode


def test_demo_creates_the_case_with_records_and_identifiers(tmp_path):
    d = tmp_path / "demo"
    got, code = zebra("case", "demo", str(d), "--lang", "zh")
    assert code == 0 and got["ok"] is True
    r = got["result"]
    assert r["path"] == str(d.resolve()) and r["reused"] is False
    assert sorted(os.listdir(d / "records")) == sorted(r["records"]) and len(r["records"]) == 2
    for name in r["records"]:
        p = d / "records" / name
        assert not p.is_symlink()  # copied: the privacy gate's case-folder check resolves real paths
        assert "合成数据" in p.read_text("utf-8")
    data = C.load(str(d))
    assert data["language"] == "zh" and "合成数据" in data["title"]
    assert r["identifiers"] == len(data["privacy"]["identifiers"]) == 6
    assert "王小雨" in data["privacy"]["identifiers"]
    assert "王小雨" not in json.dumps(r, ensure_ascii=False)  # never printed back
    assert "records" in r["first_prompt"]


def test_demo_is_reused_and_puts_back_deleted_records(tmp_path):
    d = tmp_path / "demo"
    zebra("case", "demo", str(d), "--lang", "en")
    first = sorted(os.listdir(d / "records"))
    os.remove(d / "records" / first[0])
    got, code = zebra("case", "demo", str(d), "--lang", "en")
    assert code == 0 and got["result"]["reused"] is True
    assert sorted(os.listdir(d / "records")) == first
    assert len(C.load(str(d))["privacy"]["identifiers"]) == 5


def test_demo_refuses_a_folder_holding_another_case(tmp_path):
    d = tmp_path / "real"
    C.init(str(d), title="a real case")
    got, code = zebra("case", "demo", str(d))
    assert code != 0 and got["ok"] is False
    assert "not the demo" in got["error"]["message"]
    assert not (d / "records" / "门诊病历.txt").exists()


def test_demo_defaults_to_the_home_folder(tmp_path):
    got, code = zebra("case", "demo", "--lang", "en", home=tmp_path)
    assert code == 0
    assert got["result"]["path"] == str((tmp_path / "zebra-cases" / "demo-lily").resolve())
