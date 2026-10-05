import json
import os
import shutil
import subprocess
import sys

import pytest

from zebra import case as C
from zebra import hpo_local

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIX = os.path.join(os.path.dirname(__file__), "fixtures", "hpo")


@pytest.fixture
def case_dir(tmp_path, monkeypatch):
    hdir = tmp_path / "hpo"
    hdir.mkdir()
    for f in hpo_local.FILES:
        shutil.copy(os.path.join(FIX, f), hdir / f)
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(hdir))
    d = tmp_path / "case1"
    C.init(str(d), title="demo", role="family")
    return str(d)


def zebra(*args, env_extra=None):
    env = dict(os.environ, **(env_extra or {}))
    out = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "zebra"), "--json", *args], capture_output=True, text=True, env=env)
    return json.loads(out.stdout), out.returncode


def test_init_layout(case_dir):
    for sub in ("evidence", "records", "reports"):
        assert os.path.isdir(os.path.join(case_dir, sub))
    assert open(os.path.join(case_dir, ".gitignore")).read().strip().endswith("*")
    with pytest.raises(C.CaseError):
        C.init(case_dir)


def test_phenotype_validation(case_dir):
    with pytest.raises(C.CaseError):
        C.add_phenotype(case_dir, "HP:12", "x")
    C.add_phenotype(case_dir, "HP:0001250", "Seizure")
    C.add_phenotype(case_dir, "HP:0001250", "Seizure", status="excluded")
    data = C.load(case_dir)
    assert len(data["phenotypes"]) == 1 and data["phenotypes"][0]["status"] == "excluded"


def test_hypothesis_merge_and_acmg(case_dir):
    C.add_hypothesis(case_dir, "Dravet syndrome", ids={"ORPHA": "33069"}, support=["E1"])
    h = C.add_hypothesis(case_dir, "dravet syndrome", status="leading", support=["E2"])
    assert h["status"] == "leading" and h["support"] == ["E1", "E2"]
    v = C.add_variant(case_dir, gene="SCN1A", hgvs_c="NM_001165963.4:c.2134C>T")
    C.set_acmg(case_dir, v["id"], "Pathogenic", 13, ["PVS1", "PS2", "PM2_Supporting"])
    s = C.summary(case_dir)
    assert s["variants"][0]["classification"] == "Pathogenic"


def test_ledger_ids_are_sequential(case_dir):
    a = C.append_ledger(case_dir, "x", {}, [{"db": "A", "record": "1"}, {"db": "B", "record": "2"}])
    b = C.append_ledger(case_dir, "y", {}, [{"db": "C", "record": "3"}])
    assert a == ["E1", "E2"] and b == ["E3"]
    assert [r["eid"] for r in C.read_ledger(case_dir)] == ["E1", "E2", "E3"]


def test_identifiers_not_echoed(case_dir):
    env, code = zebra("case", "identifiers", "--case", case_dir, "--add", "张三", "2019-03-02")
    assert code == 0 and "张三" not in json.dumps(env, ensure_ascii=False)
    assert C.load(case_dir)["privacy"]["identifiers"] == ["2019-03-02", "张三"]


def test_apply_verifies_hpo_offline(case_dir):
    ops = {"phenotypes": [{"id": "HP:0002373", "source": "records/a.pdf p1"}, {"id": "HP:0000000"}],
           "variants": [{"gene": "SCN1A", "hgvs_c": "NM_001165963.4:c.2134C>T", "zygosity": "het"}],
           "hypotheses": [{"disease": "Dravet syndrome", "ids": ["ORPHA:33069"], "status": "leading"}],
           "acmg": [{"variant_id": "v1", "codes": ["PVS1", "PM2_Supporting"]}],
           "questions": ["Were the parents tested?"]}
    env, code = zebra("case", "apply", case_dir, "--ops", json.dumps(ops), env_extra={"ZEBRA_HPO_DIR": os.environ["ZEBRA_HPO_DIR"]})
    assert code == 0, env
    res = env["result"]
    assert res["counts"] == {"phenotypes": 1, "variants": 1, "hypotheses": 1, "leads": 0}
    assert any("HP:0000000" in e for e in res["errors"])
    data = C.load(case_dir)
    assert data["phenotypes"][0]["label"].startswith("Febrile seizure")
    assert data["variants"][0]["zebra_acmg"]["classification"] == "Likely pathogenic"
    assert env["ledger"] == []  # case bookkeeping is not evidence


def test_cli_envelope_and_errors(case_dir):
    env, code = zebra("acmg", "classify", "PVS1", "PS2")
    assert code == 0 and env["ok"] and env["result"]["classification"] == "Pathogenic"
    env, code = zebra("acmg", "classify", "PX9")
    assert code == 2 and env["ok"] is False and env["error"]["type"] == "UsageError"
    env, code = zebra("case", "summary", "/nonexistent/zebra-case")
    assert code == 2 and "no case" in env["error"]["message"]


def test_ledger_written_for_evidence_commands(case_dir):
    env, code = zebra("hpo", "rank", "HP:0002373", "HP:0001263", "--case", case_dir, env_extra={"ZEBRA_HPO_DIR": os.environ["ZEBRA_HPO_DIR"]})
    assert code == 0 and env["ledger"] == ["E1"]
    assert C.read_ledger(case_dir)[0]["db"].startswith("HPO annotations")
