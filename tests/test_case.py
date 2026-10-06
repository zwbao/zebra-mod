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
    # the HPO verification is a source lookup, so it gets an evidence id: a
    # hypothesis can then cite the phenotype it rests on
    assert env["ledger"] == ["E1"]


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


# ===================================================================
# Regression tests for the defects found in review (docs/ROADMAP.md).
# Each is named after its backlog id.
# ===================================================================

# ------------------------------------------------- E2: hypothesis status

def test_e2_hypothesis_status_is_kept_when_none_is_given(case_dir):
    """Attaching evidence must not demote a leading or confirmed hypothesis."""
    C.add_hypothesis(case_dir, "Dravet syndrome", status="leading")
    h = C.add_hypothesis(case_dir, "dravet syndrome", note="new evidence")
    assert h["status"] == "leading"
    env, code = zebra("case", "apply", case_dir,
                      "--ops", json.dumps({"hypotheses": [{"disease": "Dravet syndrome", "note": "more"}]}))
    assert code == 0 and env["result"]["applied"]["hypotheses"][0]["status"] == "leading"
    # an explicit status still changes it, and a new hypothesis starts as considered
    assert C.add_hypothesis(case_dir, "Dravet syndrome", status="confirmed")["status"] == "confirmed"
    assert C.add_hypothesis(case_dir, "GLUT1 deficiency")["status"] == "considered"
    env, code = zebra("case", "add-hypothesis", "--case", case_dir, "Dravet syndrome", "--note", "cli")
    assert code == 0 and env["result"]["status"] == "confirmed"


# ------------------------------------------- E3 / F8: validate before writing

def test_e3_wrong_types_are_refused_and_nothing_is_written(case_dir):
    """A string where a list belongs used to be stored one character per item."""
    before = json.dumps(C.load(case_dir), sort_keys=True)
    for ops, needle in (
        ({"hypotheses": [{"disease": "Dravet", "support": "E12"}]}, "must be a list of strings"),
        ({"hypotheses": [{"disease": "Dravet", "ids": "ORPHA:33069"}]}, "must be a list of strings"),
        ({"hypotheses": [{"disease": 123}]}, "disease must be a non-empty string"),
        ({"questions": "why"}, "must be a list of strings"),
        ({"phenotypes": ["HP:0001250"]}, "every item must be an object"),
        ({"remove": ["v1"]}, "every item must be an object"),
        ({"hypotheses": {"disease": "Dravet"}}, "wrap the single object in a list"),
        ({"leads": [{"name": "ketogenic diet"}]}, "kind must be one of"),
        ({"variants": [{"gene": "SCN1A"}]}, "needs at least one of"),
        ({"variants": [{"hgvs_c": "exon 2 deleted"}]}, "is not HGVS"),
        ({"profile": {"nickname": "x"}}, "unknown profile field"),
    ):
        env, code = zebra("case", "apply", case_dir, "--ops", json.dumps(ops))
        assert code == 2, (ops, env)
        assert env["error"]["type"] == "UsageError" and needle in env["error"]["message"], (ops, env)
        assert "nothing was written" in env["error"]["message"] or "accepted" in env["error"]["message"]
    assert json.dumps(C.load(case_dir), sort_keys=True) == before


def test_e3_a_non_string_disease_can_no_longer_poison_later_writes(case_dir):
    env, code = zebra("case", "apply", case_dir, "--ops", json.dumps({"hypotheses": [{"disease": 123}]}))
    assert code == 2
    env, code = zebra("case", "apply", case_dir, "--ops", json.dumps({"hypotheses": [{"disease": "Dravet"}]}))
    assert code == 0 and env["result"]["counts"]["hypotheses"] == 1


def test_f8_unknown_keys_are_refused_with_the_accepted_list(case_dir):
    ops = {"phenotype": [{"id": "HP:0001250"}], "therapy_leads": [{"name": "x", "kind": "trial"}]}
    env, code = zebra("case", "apply", case_dir, "--ops", json.dumps(ops))
    assert code == 2
    message = env["error"]["message"]
    assert "'phenotype' (did you mean 'phenotypes'?)" in message
    assert "'therapy_leads' (did you mean 'leads'?)" in message
    for accepted in ("profile", "phenotypes", "variants", "hypotheses", "leads", "acmg", "questions", "remove"):
        assert accepted in message


def test_f8_a_batch_is_applied_in_one_write(case_dir, monkeypatch):
    """Every op lands in a single locked write: a crash cannot leave half of it behind."""
    writes = []
    original = C._write
    monkeypatch.setattr(C, "_write", lambda root, data: (writes.append(1), original(root, data))[1])
    ops = {"profile": {"title": "demo 2"},
           "phenotypes": [{"id": "HP:0002373"}],
           "variants": [{"gene": "SCN1A", "hgvs_c": "NM_001165963.4:c.2134C>T"}],
           "hypotheses": [{"disease": "Dravet syndrome", "status": "leading"}],
           "leads": [{"name": "ketogenic diet", "kind": "supportive"}],
           "acmg": [{"variant_id": "v1", "codes": ["PVS1", "PM2_Supporting"]}],
           "questions": ["Were the parents tested?"],
           "remove": [{"kind": "lead", "id": "t99"}]}
    from zebra.commands import case as case_cmd

    args = type("A", (), {"ops": json.dumps(ops), "dir": case_dir, "case": None})()
    out = case_cmd._apply(args)
    assert len(writes) == 1, f"{len(writes)} writes instead of one"
    assert out.result["counts"] == {"phenotypes": 1, "variants": 1, "hypotheses": 1, "leads": 1}
    assert out.result["applied"]["removed"] == [{"kind": "lead", "id": "t99", "ok": False}]
    data = C.load(case_dir)
    assert data["title"] == "demo 2" and data["variants"][0]["zebra_acmg"]["classification"] == "Likely pathogenic"


# -------------------------------------------- F34: making the doctrine checkable

def test_f34_identifier_shapes_are_validated(case_dir):
    assert C.check_id("orpha:33069") == "ORPHA:33069"
    assert C.check_id("OMIM:607208") == "OMIM:607208"
    assert C.check_id("MONDO:0100135") == "MONDO:0100135"
    assert C.check_id("NCT:04006210") == "NCT:04006210"
    assert C.check_id("PMID:28919360") == "PMID:28919360"
    for bad, needle in (("OMIM:99", "well-formed OMIM"),
                        ("HP:1250", "well-formed HP"),
                        ("MONDO:135", "well-formed MONDO"),
                        ("NCT:4006210", "well-formed NCT"),
                        ("WEIRD:1", "unknown identifier prefix"),
                        ("33069", "not an identifier")):
        with pytest.raises(C.CaseError, match=needle):
            C.check_id(bad)
    assert C.check_hgvs("NM_000492.4:c.1521_1523del").startswith("NM_000492.4")
    assert C.check_hgvs("NM_004006.3:c.6439-?_7309+?del")  # the exon-boundary form
    for bad in ("SCN1A p.Arg712Ter", "NM_000492.4", "exon 45-50 deleted", "1521del"):
        with pytest.raises(C.CaseError, match="is not HGVS"):
            C.check_hgvs(bad)
    env, code = zebra("case", "apply", case_dir,
                      "--ops", json.dumps({"hypotheses": [{"disease": "Dravet", "ids": ["OMIM:99"]}]}))
    assert code == 2 and "well-formed OMIM" in env["error"]["message"]


def test_f34_evidence_ids_must_exist_in_the_ledger(case_dir):
    """`support: ["E7"]` with an empty ledger used to be accepted in silence."""
    env, code = zebra("case", "apply", case_dir,
                      "--ops", json.dumps({"hypotheses": [{"disease": "Dravet", "support": ["E7"]}]}))
    assert code == 2 and env["error"]["type"] == "UsageError"
    assert "E7 is not in this case's evidence ledger" in env["error"]["message"]
    assert "the ledger is empty" in env["error"]["message"]
    assert C.load(case_dir)["hypotheses"] == []

    C.append_ledger(case_dir, "hpo term", {}, [{"db": "HPO", "record": "HP:0001250"}])
    env, code = zebra("case", "apply", case_dir,
                      "--ops", json.dumps({"hypotheses": [{"disease": "Dravet", "support": ["E1"]}]}))
    assert code == 0 and env["result"]["applied"]["hypotheses"][0]["support"] == ["E1"]
    env, code = zebra("case", "apply", case_dir,
                      "--ops", json.dumps({"hypotheses": [{"disease": "Dravet", "against": ["E2"]}]}))
    assert code == 2 and "E2 is not in this case" in env["error"]["message"]
    env, code = zebra("case", "apply", case_dir,
                      "--ops", json.dumps({"hypotheses": [{"disease": "Dravet", "support": ["seven"]}]}))
    assert code == 2 and "is not an evidence id" in env["error"]["message"]
    # the same check on the single-op commands
    env, code = zebra("case", "add-hypothesis", "--case", case_dir, "Dravet", "--support", "E9")
    assert code == 2 and "E9 is not in this case" in env["error"]["message"]
    env, code = zebra("case", "add-lead", "--case", case_dir, "cannabidiol", "--kind", "approved",
                      "--evidence", "E9")
    assert code == 2 and "E9 is not in this case" in env["error"]["message"]
    env, code = zebra("case", "add-lead", "--case", case_dir, "cannabidiol", "--kind", "approved",
                      "--evidence", "E1", "--id", "NCT:04006210")
    assert code == 0 and env["result"]["evidence"] == ["E1"] and env["result"]["ids"] == {"NCT": "04006210"}


# ------------------------------------------------- A-P2-8: ledger numbering

def test_a_p2_8_ledger_ids_come_from_the_maximum_id(case_dir):
    """Ids must come from the highest id ever issued, and a partial last line must be repaired."""
    C.append_ledger(case_dir, "a", {}, [{"db": "A"}, {"db": "B"}, {"db": "C"}])
    path = os.path.join(case_dir, "evidence", "ledger.jsonl")
    rows = open(path).read().splitlines()
    open(path, "w").write(rows[0] + "\n" + rows[2] + "\n")  # E2 deleted by hand
    assert C.append_ledger(case_dir, "b", {}, [{"db": "D"}]) == ["E4"]  # not E3 again

    with open(path, "a") as fh:
        fh.write('{"eid": "E5", "command": "crash", "db": "E", "url": "htt')  # no newline: a crash mid-write
    assert C.append_ledger(case_dir, "c", {}, [{"db": "F"}]) == ["E6"]
    text = open(path).read()
    assert text.endswith("\n") and '"eid": "E6"' in text
    eids = [r["eid"] for r in C.read_ledger(case_dir)]
    assert eids == ["E1", "E3", "E4", "E6"]  # the truncated E5 row is unreadable, but its id was not re-used
    assert C.ledger_ids(case_dir) == {"E1", "E3", "E4", "E6"}


# ------------------------------------------- A-P2-10: hand-edited case.json

def test_a_p2_10_hand_edited_case_json_is_validated(case_dir):
    path = os.path.join(case_dir, "case.json")
    data = json.loads(open(path).read())
    data.pop("privacy")
    data["phenotypes"] = None
    data.pop("questions")
    json.dump(data, open(path, "w"))
    s = C.summary(case_dir)  # used to raise KeyError / TypeError
    assert s["phenotypes"] == [] and s["identifiers"] == 0 and s["questions"] == []

    data["phenotypes"] = {"id": "HP:0001250"}
    json.dump(data, open(path, "w"))
    with pytest.raises(C.CaseError, match="'phenotypes' must be a list"):
        C.summary(case_dir)
    data["phenotypes"] = ["HP:0001250"]
    json.dump(data, open(path, "w"))
    with pytest.raises(C.CaseError, match="every item in 'phenotypes' must be an object"):
        C.summary(case_dir)
    data["phenotypes"] = []
    data["privacy"] = {"identifiers": "zhang"}
    json.dump(data, open(path, "w"))
    with pytest.raises(C.CaseError, match="privacy.identifiers must be a list"):
        C.summary(case_dir)
    env, code = zebra("case", "summary", case_dir)
    assert code == 2 and env["error"]["type"] == "UsageError"


# ===================================================================
# The envelope contract: F5 (trimming), F6 (deadline), F7 (always JSON,
# no flag smuggling), F17 (no NaN/Infinity).
# ===================================================================

def test_f7_argparse_errors_and_help_are_still_json(case_dir):
    env, code = zebra("case", "add-lead", "--case", case_dir, "x")  # a missing required option
    assert code == 2 and env["ok"] is False and env["error"]["type"] == "UsageError"
    assert "usage: zebra case add-lead" in env["error"]["message"] and "--kind" in env["error"]["message"]
    env, code = zebra("acmg", "--bogus")
    assert code == 2 and env["error"]["type"] == "UsageError"
    env, code = zebra("--help")
    assert code == 0 and env["ok"] and "usage: zebra" in env["result"]["stdout"]
    env, code = zebra("--version")
    assert code == 0 and env["result"]["stdout"].startswith("zebra ")
    env, code = zebra("vcf")  # a named command with no action: that command's own usage
    assert code == 2 and env["ok"] is False and env["error"]["type"] == "UsageError"
    assert env["command"] == "vcf" and env["error"]["message"].startswith("usage: zebra vcf")
    # an abbreviated --json must give JSON too (argparse accepts any unambiguous prefix)
    env, code = zebra("--js", "acmg", "classify", "PVS1")
    assert code == 0 and env["ok"] and env["result"]["classification"]
    env, code = zebra("--jso", "nosuchcommand")
    assert code == 2 and env["ok"] is False


def test_f7_positional_values_cannot_smuggle_a_flag(case_dir):
    env, code = zebra("variant", "--", "--help")
    assert code == 2 and env["error"]["type"] == "UsageError"
    assert "looks like a flag" in env["error"]["message"]
    env, code = zebra("case", "add-question", "--case", case_dir, "--", "-rf /")
    assert code == 2 and "looks like a flag" in env["error"]["message"]
    # a legitimate value that merely contains a dash still works
    env, code = zebra("case", "add-question", "--case", case_dir, "Was the EEG re-read?")
    assert code == 0 and env["result"] == ["Was the EEG re-read?"]


def test_f7_an_unexpected_exception_is_an_internal_error_envelope(monkeypatch, capsys):
    from zebra import acmg as acmg_mod
    from zebra import cli

    def boom(*a, **kw):
        raise RuntimeError("the floor gave way")

    monkeypatch.setattr(acmg_mod, "classify", boom)
    code = cli.main(["--json", "acmg", "classify", "PVS1"])
    env = json.loads(capsys.readouterr().out)
    assert code == 70 and env["ok"] is False
    assert env["error"]["type"] == "InternalError" and env["error"]["exception"] == "RuntimeError"
    assert env["error"]["message"] == "the floor gave way"
    assert env["error"]["traceback_tail"].rstrip().endswith("RuntimeError: the floor gave way")


def test_f17_non_finite_numbers_become_null_with_a_warning(capsys):
    from zebra import cli
    from zebra.core import Outcome

    env = cli._envelope("x", Outcome({"a": float("nan"), "b": [1.0, float("inf")]},
                                     query={"q": float("-inf")}), [])
    text = json.dumps(env)  # would raise ValueError if a non-finite number survived
    assert json.loads(text)["result"] == {"a": None, "b": [1.0, None]}
    assert env["query"] == {"q": None}
    assert any("not a finite number" in w and "result.a" in w for w in env["warnings"])
    with pytest.raises(ValueError):
        json.dumps({"x": float("nan")}, allow_nan=False)
    # argparse refuses a non-finite number for a float option, with a JSON envelope
    env, code = zebra("stats", "maxaf", "--prevalence", "0.0002", "--allelic", "0.1", "--faf95", "nan")
    assert code == 2 and env["ok"] is False


def test_f5_an_oversized_result_is_trimmed_and_stays_valid_json(case_dir, monkeypatch, capsys):
    from zebra import cli
    from zebra.core import Outcome

    rows = [{"n": i, "pad": "x" * 60} for i in range(400)]
    env = cli._envelope("x", Outcome({"hits": rows, "note": "kept"}, warnings=["a real warning"],
                                     sources=[{"db": "Europe PMC"}]), ["E1", "E2"])
    text = cli._trim(env, 4000)
    assert len(text) <= 4000
    assert len(text) > 3600  # the budget is packed, not emptied: as many items as fit are kept
    got = json.loads(text)  # still valid JSON
    assert got["warnings"][0] == "a real warning"
    assert got["sources"] == [{"db": "Europe PMC"}] and got["ledger"] == ["E1", "E2"]
    assert got["result"]["note"] == "kept" and 0 < len(got["result"]["hits"]) < 400
    dropped = [w for w in got["warnings"] if "trimmed" in w]
    assert dropped and "items dropped from result.hits" in dropped[0]
    # warnings, ledger and sources are serialised before the result, so a cut cannot reach them first
    assert list(env) == ["ok", "zebra", "command", "query", "warnings", "ledger", "sources", "result"]
    # and the budget is honoured end to end
    monkeypatch.setenv("ZEBRA_MAX_BYTES", "900")
    assert cli.main(["--json", "case", "summary", case_dir]) == 0
    out = capsys.readouterr().out
    assert len(out.strip()) <= 900 and json.loads(out)["ok"]


def test_f6_a_deadline_clamps_each_attempt_and_stops_retries(monkeypatch):
    import time

    from zebra import http
    from zebra.http import SourceError

    monkeypatch.setenv("ZEBRA_DEADLINE_MS", "60000")
    monkeypatch.setattr(http, "_PROCESS_START", time.monotonic())
    assert 59 < http.deadline_seconds() <= 60
    monkeypatch.setattr(http, "_PROCESS_START", time.monotonic() - 59.5)
    assert http.deadline_seconds() < 1

    seen = {"timeouts": [], "sleeps": []}

    class FakeOpener:
        def open(self, req, timeout=None):
            seen["timeouts"].append(timeout)
            raise TimeoutError("timed out")

    def advancing_sleep(seconds):
        """Spend the budget instead of real time: moving the start back leaves less of it."""
        seen["sleeps"].append(seconds)
        http._PROCESS_START -= seconds

    monkeypatch.setattr(http, "_opener", lambda: FakeOpener())
    monkeypatch.setattr(http.time, "sleep", advancing_sleep)

    monkeypatch.setattr(http, "_PROCESS_START", time.monotonic() - 50)  # 10 s left of 60
    with pytest.raises(SourceError, match="network error"):
        http.request("https://example.invalid/x", source="Example", timeout=180, retries=3, cache_ttl=0)
    assert max(seen["timeouts"]) <= 10  # clamped to the time left, not the 180 s asked for
    assert seen["timeouts"] == sorted(seen["timeouts"], reverse=True)  # each attempt gets less
    assert len(seen["timeouts"]) == 3 and sum(seen["sleeps"]) <= 10  # the 4th retry is past the deadline

    seen["timeouts"].clear()
    monkeypatch.setattr(http, "_PROCESS_START", time.monotonic() - 61)  # the budget is spent
    with pytest.raises(SourceError, match="deadline reached"):
        http.request("https://example.invalid/x", source="Example", timeout=180, retries=3, cache_ttl=0)
    assert seen["timeouts"] == []  # not even one attempt is started

    # a small but positive budget still gets its attempt: a 1.5 s deadline must
    # not turn every source off (the attempt timeout is clamped, not refused)
    seen["timeouts"].clear()
    monkeypatch.setenv("ZEBRA_DEADLINE_MS", "1500")
    monkeypatch.setattr(http, "_PROCESS_START", time.monotonic())
    with pytest.raises(SourceError, match="network error"):
        http.request("https://example.invalid/x", source="Example", timeout=180, retries=0, cache_ttl=0)
    assert len(seen["timeouts"]) == 1 and 1.0 < seen["timeouts"][0] <= 1.5

    # a non-finite or malformed budget is no budget at all, not an infinite one
    for bad in ("nan", "inf", "-nan", "1e9999", "500ms", "", "0", "-1000"):
        monkeypatch.setenv("ZEBRA_DEADLINE_MS", bad)
        assert http.deadline_seconds() is None, bad
    monkeypatch.setenv("ZEBRA_DEADLINE_MS", "60000")

    monkeypatch.delenv("ZEBRA_DEADLINE_MS")
    assert http.deadline_seconds() is None
    seen["timeouts"].clear()
    with pytest.raises(SourceError, match="network error"):
        http.request("https://example.invalid/x", source="Example", timeout=7, retries=1, cache_ttl=0)
    assert seen["timeouts"] == [7, 7]  # no deadline: the given timeout is used as before


def test_f6_a_timed_out_source_becomes_a_warning_not_a_dead_call():
    """`core.attempt` turns the deadline SourceError into a named warning, so the rest of the card stands."""
    from zebra.core import attempt
    from zebra.http import SourceError

    warnings = []

    def hung():
        raise SourceError("ClinicalTrials.gov", "https://clinicaltrials.gov/api", None,
                          "deadline reached (no time left): gave up after 1 attempt(s) without an answer")

    assert attempt("trials", hung, warnings) is None
    assert len(warnings) == 1 and "deadline reached" in warnings[0] and "trials unavailable" in warnings[0]
    assert attempt("shape", lambda: {}["missing"], warnings) is None
    assert "unexpected response shape" in warnings[1]


def test_f34_hgvs_without_a_reference_is_recorded_with_a_warning(case_dir):
    """`p.Phe508del` is how reports write it: record it and say the transcript is missing, do not refuse it."""
    assert C.check_hgvs("p.Phe508del") == "p.Phe508del"
    assert C.hgvs_reference("p.Phe508del") is None
    assert C.hgvs_reference("NM_000492.4:c.1521_1523del") == "NM_000492.4"
    assert C.hgvs_reference("SCN1A:c.2134C>T") == "SCN1A"
    notes = C.variant_notes({"hgvs_p": "p.Phe508del", "hgvs_c": "NM_000492.4:c.1521_1523del"})
    assert len(notes) == 1 and "names no reference sequence" in notes[0]
    assert C.variant_notes({"hgvs_c": "NM_000492.4:c.1521_1523del"}) == []

    env, code = zebra("case", "apply", case_dir,
                      "--ops", json.dumps({"variants": [{"gene": "CFTR", "hgvs_p": "p.Phe508del"}]}))
    assert code == 0 and env["result"]["applied"]["variants"][0]["hgvs_p"] == "p.Phe508del"
    assert any("names no reference sequence" in w for w in env["warnings"])
    env, code = zebra("case", "add-variant", "--case", case_dir, "--gene", "SCN1A", "--hgvs-c", "c.2134C>T")
    assert code == 0 and any("names no reference sequence" in w for w in env["warnings"])
    # a string that is not HGVS at all is still refused
    env, code = zebra("case", "apply", case_dir,
                      "--ops", json.dumps({"variants": [{"gene": "CFTR", "hgvs_c": "exon 2 missing"}]}))
    assert code == 2 and "is not HGVS" in env["error"]["message"]


def test_f7_a_stdout_that_cannot_encode_the_answer_still_gets_json(case_dir):
    """A non-UTF-8 stdout used to lose the envelope entirely and print a traceback."""
    for env_extra in ({"PYTHONIOENCODING": "ascii"},):
        env, code = zebra("china", "白化病", env_extra=env_extra)
        assert code == 0 and env["ok"], env_extra
        assert env["query"]["query"] == "白化病"  # round-trips through \uXXXX escapes
        env, code = zebra("variant", "白化病", env_extra=env_extra)
        assert code in (2, 3) and env["ok"] is False
        assert "白化病" in json.dumps(env, ensure_ascii=False) or env["error"]["message"]


def test_f5_the_budget_covers_the_whole_envelope(case_dir):
    """An oversized query or source list used to sail past the budget untouched."""
    big = "A" * 20_000
    env, code = zebra("china", big, env_extra={"ZEBRA_MAX_BYTES": "2000"})
    assert code == 0 and env["ok"]
    assert len(json.dumps(env)) <= 2000
    assert len(env["query"]["query"]) < 1000
    assert any("was 20000 characters and was cut" in w for w in env["warnings"])
    # an error envelope is bounded as well: it must not echo a whole bad input
    env, code = zebra("variant", big)
    assert code == 2 and len(json.dumps(env)) < 10_000
    assert "cut to stay inside the size budget" in env["error"]["message"]

    from zebra import cli
    from zebra.core import Outcome

    envelope = cli._envelope("x", Outcome({"hits": [1, 2, 3]},
                                          sources=[{"db": f"src{i}", "url": "u" * 200} for i in range(300)]), [])
    text = cli._trim(envelope, 4000)
    assert len(text) <= 4000
    got = json.loads(text)
    assert 0 < len(got["sources"]) < 300 and got["result"]["hits"] == [1, 2, 3]
    assert any("provenance record(s) dropped" in w for w in got["warnings"])


def test_f5_trimming_cuts_the_innermost_list_not_whole_records(case_dir):
    """Ranking candidate lists by bytes could only ever pick the outermost one."""
    from zebra import cli
    from zebra.core import Outcome

    variants = [{"id": f"v{i}", "gene": "SCN1A", "transcripts": [{"t": "x" * 40} for _ in range(40)]}
                for i in range(2)]
    envelope = cli._envelope("x", Outcome({"variants": variants}), [])
    got = json.loads(cli._trim(envelope, 3000))
    assert len(got["result"]["variants"]) == 2  # both variants survive
    assert any(len(v["transcripts"]) < 40 for v in got["result"]["variants"])
    assert any("transcripts" in w for w in got["warnings"])
    # a one-element wrapper around a long list is trimmed at the long list
    wrapped = cli._envelope("x", Outcome([{"hits": list(range(500))}]), [])
    got = json.loads(cli._trim(wrapped, 1200))
    assert len(got["result"]) == 1 and 0 < len(got["result"][0]["hits"]) < 500
    # a result with no list at all says so instead of blaming lists
    blob = cli._envelope("x", Outcome({"text": "z" * 5000}), [])
    got = json.loads(cli._trim(blob, 1500))
    assert len(json.dumps(got)) <= 1500
    assert any("was 5000 characters and was cut" in w for w in got["warnings"])


def test_f17_the_ledger_never_receives_a_non_finite_number(case_dir, monkeypatch):
    """The envelope was sanitised while the durable artifact kept a bare NaN."""
    from zebra import case as case_mod
    from zebra import cli
    from zebra.core import Outcome

    def handler(args):
        return Outcome({"ok": True}, sources=[{"db": "gnomAD", "record": "x", "af": float("nan")}],
                       query={"q": float("inf")})

    parser_build = cli.build_parser

    def patched():
        p = parser_build()
        cli._subparser(cli._subparser(p, "acmg"), "classify").set_defaults(func=handler, no_ledger=False)
        return p

    monkeypatch.setattr(cli, "build_parser", patched)
    assert cli.main(["--json", "acmg", "classify", "PVS1", "--case", case_dir]) == 0
    rows = case_mod.read_ledger(case_dir)
    assert rows and rows[0]["af"] is None and rows[0]["query"]["q"] is None
    raw = open(os.path.join(case_dir, "evidence", "ledger.jsonl")).read()
    assert "NaN" not in raw and "Infinity" not in raw
    json.loads(raw.splitlines()[0], parse_constant=lambda c: (_ for _ in ()).throw(ValueError(c)))


def test_f7_an_interrupted_or_exiting_handler_still_answers(monkeypatch, capsys):
    from zebra import cli

    def make(raiser):
        parser_build = cli.build_parser

        def patched():
            p = parser_build()
            cli._subparser(cli._subparser(p, "acmg"), "classify").set_defaults(func=raiser)
            return p

        return patched

    monkeypatch.setattr(cli, "build_parser", make(lambda args: (_ for _ in ()).throw(SystemExit(7))))
    assert cli.main(["--json", "acmg", "classify", "PVS1"]) == 70
    env = json.loads(capsys.readouterr().out)
    assert env["ok"] is False and env["error"]["type"] == "InternalError"
    assert "exited with status 7" in env["error"]["message"]

    monkeypatch.setattr(cli, "build_parser", make(lambda args: (_ for _ in ()).throw(KeyboardInterrupt())))
    assert cli.main(["--json", "acmg", "classify", "PVS1"]) == 130
    env = json.loads(capsys.readouterr().out)
    assert env["ok"] is False and env["error"]["type"] == "Interrupted"


# -------------- second-round review findings (same ids, deeper defects)

def test_e2_a_phenotype_keeps_its_status_and_its_note(case_dir):
    """`excluded` with a reason must not become `present` when the term is recorded again."""
    ops = {"phenotypes": [{"id": "HP:0002373", "status": "excluded", "note": "ruled out by EEG"}]}
    env, code = zebra("case", "apply", case_dir, "--ops", json.dumps(ops))
    assert code == 0 and env["result"]["applied"]["phenotypes"][0]["status"] == "excluded"
    env, code = zebra("case", "apply", case_dir, "--ops", json.dumps({"phenotypes": [{"id": "HP:0002373"}]}))
    kept = env["result"]["applied"]["phenotypes"][0]
    assert code == 0 and kept["status"] == "excluded" and kept["note"] == "ruled out by EEG"
    # an explicit status still changes it
    env, _ = zebra("case", "apply", case_dir,
                   "--ops", json.dumps({"phenotypes": [{"id": "HP:0002373", "status": "present"}]}))
    assert env["result"]["applied"]["phenotypes"][0]["status"] == "present"
    assert C.apply_phenotype({"phenotypes": []}, "HP:0001250", "Seizure")["status"] == "present"


def test_f34_case_item_ids_are_never_re_issued(case_dir):
    """`v1` quoted from an earlier summary must not come to mean a different variant."""
    first = C.add_variant(case_dir, gene="SCN1A", hgvs_c="NM_001165963.4:c.2134C>T")
    assert first["id"] == "v1"
    assert C.remove(case_dir, "variant", "v1") is True
    second = C.add_variant(case_dir, gene="CFTR", hgvs_c="NM_000492.4:c.1521_1523del")
    assert second["id"] == "v2"  # not v1 again
    assert C.add_hypothesis(case_dir, "A")["id"] == "h1"
    C.remove(case_dir, "hypothesis", "h1")
    assert C.add_hypothesis(case_dir, "B")["id"] == "h2"
    # and the count survives a reload, because the water mark is stored
    assert C.load(case_dir)["_issued_ids"]["v"] == 2


def test_a_p2_8_deleting_the_last_ledger_row_does_not_re_issue_its_id(case_dir):
    C.append_ledger(case_dir, "x", {}, [{"db": "A"}, {"db": "B"}, {"db": "C"}])
    path = os.path.join(case_dir, "evidence", "ledger.jsonl")
    rows = open(path).read().splitlines()
    open(path, "w").write("\n".join(rows[:2]) + "\n")  # E3 deleted by hand
    assert C.append_ledger(case_dir, "y", {}, [{"db": "D"}]) == ["E4"]


def test_a_p2_10_a_ledger_row_that_is_not_an_object_does_not_break_the_case(case_dir):
    path = os.path.join(case_dir, "evidence", "ledger.jsonl")
    with open(path, "w") as fh:
        fh.write('{"eid": "E1", "db": "X"}\nnull\n[1, 2]\n3\n{"db": "no eid"}\n')
    assert [r["eid"] for r in C.read_ledger(case_dir) if r.get("eid")] == ["E1"]
    assert C.ledger_ids(case_dir) == {"E1"}
    env, code = zebra("case", "apply", case_dir, "--ops", json.dumps({"questions": ["still works"]}))
    assert code == 0 and env["ok"]
    env, code = zebra("case", "ledger", "--case", case_dir)
    assert code == 0 and [r.get("eid") for r in env["result"]] == ["E1", None]


def test_a_p2_10_nested_objects_are_validated_too(case_dir):
    path = os.path.join(case_dir, "case.json")
    data = json.loads(open(path).read())
    data["variants"] = [{"id": "v1", "zebra_acmg": "pathogenic"}]
    json.dump(data, open(path, "w"))
    with pytest.raises(C.CaseError, match="zebra_acmg must be an object"):
        C.summary(case_dir)
    data["variants"] = []
    data["hypotheses"] = [{"id": "h1", "disease": "X", "support": "E12"}]
    json.dump(data, open(path, "w"))
    with pytest.raises(C.CaseError, match="must be a list of evidence ids"):
        C.summary(case_dir)


def test_f8_an_unknown_field_inside_an_item_is_refused(case_dir):
    for ops, needle in (
        ({"hypotheses": [{"disease": "Angelman", "statuss": "leading"}]}, "unknown field(s) statuss"),
        ({"leads": [{"name": "X", "kind": "trial", "evidenc": ["E1"]}]}, "unknown field(s) evidenc"),
        ({"phenotypes": [{"id": "HP:0001250", "onsett": "6m"}]}, "unknown field(s) onsett"),
        ({"acmg": [{"variant_id": "v1", "codes": ["PVS1"], "notes": "x"}]}, "unknown field(s) notes"),
        ({"remove": [{"kind": "variant", "ids": "v1"}]}, "unknown field(s) ids"),
        ({"profile": {"role": "doctor"}}, "role must be one of"),
        ({"profile": {"age": {"a": 1}}}, "profile.age must be a string"),
    ):
        env, code = zebra("case", "apply", case_dir, "--ops", json.dumps(ops))
        assert code == 2, (ops, env)
        assert needle in env["error"]["message"], (ops, env)
        assert "nothing was written" in env["error"]["message"]


def test_f8_an_invalid_acmg_code_is_caught_before_anything_is_written(case_dir):
    C.add_variant(case_dir, gene="SCN1A", hgvs_c="NM_001165963.4:c.2134C>T")
    before = open(os.path.join(case_dir, "case.json")).read()
    ops = {"questions": ["is this committed?"], "acmg": [{"variant_id": "v1", "codes": ["PVS9"]}]}
    env, code = zebra("case", "apply", case_dir, "--ops", json.dumps(ops))
    assert code == 2 and "PVS9" in env["error"]["message"]
    assert open(os.path.join(case_dir, "case.json")).read() == before


def test_f34_registry_ids_are_accepted_in_the_form_registries_print(case_dir):
    assert C.check_id("NCT04006210") == "NCT:04006210"
    assert C.check_id("rs113993960") == "DBSNP:rs113993960"
    assert C.check_id("PMC5760072") == "PMC:5760072"
    env, code = zebra("case", "add-lead", "--case", case_dir, "nusinersen", "--kind", "approved",
                      "--id", "NCT04006210")
    assert code == 0 and env["result"]["ids"] == {"NCT": "04006210"}
    with pytest.raises(C.CaseError, match="not a well-formed NCT"):
        C.check_id("NCT123")
    # and the shape check is still a shape check
    for bad in ("ZZZ:c.rubbish", "A:p.X"):
        with pytest.raises(C.CaseError, match="is not HGVS"):
            C.check_hgvs(bad)


def test_f34_a_case_folder_that_cannot_be_written_gives_a_named_error(case_dir, tmp_path):
    readonly = tmp_path / "ro"
    C.init(str(readonly))
    os.chmod(readonly, 0o500)
    try:
        with pytest.raises(C.CaseError, match="could not take the lock|could not write"):
            C.add_question(str(readonly), "hi")
    finally:
        os.chmod(readonly, 0o700)
    # and a mistyped path is not created on the way to failing
    missing = tmp_path / "typo"
    with pytest.raises(C.CaseError, match="no case at"):
        C.add_question(str(missing), "hi")
    assert not missing.exists()


def test_f34_case_apply_ledgers_the_hpo_terms_it_verified(case_dir):
    """A phenotype recorded through `apply` must get an evidence id, or nothing can cite it."""
    env, code = zebra("case", "apply", case_dir, "--case", case_dir,
                      "--ops", json.dumps({"phenotypes": [{"id": "HP:0002373"}]}))
    assert code == 0 and env["ledger"] == ["E1"]
    assert C.read_ledger(case_dir)[0]["db"].startswith("HPO")
    env, code = zebra("case", "apply", case_dir,
                      "--ops", json.dumps({"hypotheses": [{"disease": "Dravet syndrome", "support": ["E1"]}]}))
    assert code == 0 and env["result"]["applied"]["hypotheses"][0]["support"] == ["E1"]


def test_f8_item_values_must_have_the_right_type(case_dir):
    """`region`, `status`, `note` and friends used to accept any JSON at all."""
    for ops, needle in (
        ({"variants": [{"kind": "cnv", "region": "15:28500000-23000000"}]}, "ends before it starts"),
        ({"variants": [{"kind": "cnv", "region": "the big one"}]}, "must look like 15:23000000-28500000"),
        ({"variants": [{"kind": "exon_cnv", "gene": "DMD", "exons": {"a": 1}}]}, "exons must be a string"),
        ({"leads": [{"name": "X", "kind": "trial", "status": {"a": 1}}]}, "status must be a string"),
        ({"hypotheses": [{"disease": "X", "note": [1]}]}, "note must be a string"),
        ({"phenotypes": [{"id": "HP:0001250", "onset": {"a": 1}}]}, "onset must be a string"),
    ):
        env, code = zebra("case", "apply", case_dir, "--ops", json.dumps(ops))
        assert code == 2 and needle in env["error"]["message"], (ops, env)
    # the dict form `case show` prints round-trips back in
    env, code = zebra("case", "apply", case_dir,
                      "--ops", json.dumps({"hypotheses": [{"disease": "Dravet", "ids": {"ORPHA": "33069"}}]}))
    assert code == 0 and env["result"]["applied"]["hypotheses"][0]["ids"] == {"ORPHA": "33069"}


def test_f6_a_source_that_answers_slowly_is_bounded_too(monkeypatch):
    """A peer that sends one byte at a time never trips a socket timeout; the deadline must still bite."""
    import time

    from zebra import http
    from zebra.http import SourceError

    monkeypatch.setenv("ZEBRA_DEADLINE_MS", "3000")
    monkeypatch.setattr(http, "_PROCESS_START", time.monotonic())

    class Trickle:
        """read1 returns at most one raw packet, which is why http.py must use it."""

        def __init__(self):
            self.sent = 0

        def read1(self, n):
            self.sent += 1
            http._PROCESS_START -= 1.0  # a second of the budget goes by per byte
            return b"x"

        def read(self, n=None):  # read() would block until the whole body arrived
            raise AssertionError("_read_body must not call read()")

    with pytest.raises(SourceError, match="deadline reached after"):
        http._read_body(Trickle(), "probe", "http://example.invalid/slow")

    # a body that arrives inside the budget is returned whole
    class Quick:
        def __init__(self):
            self.chunks = [b"abc", b"def", b""]

        def read1(self, n):
            return self.chunks.pop(0)

    monkeypatch.setattr(http, "_PROCESS_START", time.monotonic())
    assert http._read_body(Quick(), "probe", "u") == b"abcdef"
    # with no deadline the body is read in one call, as before
    monkeypatch.delenv("ZEBRA_DEADLINE_MS")

    class Plain:
        def read(self, n=None):
            return b"whole body"

    assert http._read_body(Plain(), "probe", "u") == b"whole body"


def test_f34_apply_ledgers_without_an_explicit_case_flag(case_dir):
    """The mod calls `case apply <dir> --ops …` with no --case: the sources must still be recorded."""
    env, code = zebra("case", "apply", case_dir, "--ops", json.dumps({"phenotypes": [{"id": "HP:0002373"}]}))
    assert code == 0 and env["ledger"] == ["E1"]
    assert C.ledger_ids(case_dir) == {"E1"}


def test_cp1_11_family_tests_timeline_identifiers_ops(tmp_path):
    """CP1-11/SP1: relatives, tests done, a timeline and identifiers are first-class case ops."""
    import json as _json

    from zebra import case as case_mod
    from zebra.commands import case as case_cmd

    d = str(tmp_path / "c")
    case_mod.init(d, title="t")
    ops = {
        "family": [{"relation": "Mother", "affected": False, "genotype": "het for v1"}, {"relation": "brother", "affected": True}],
        "tests": [{"type": "CMA", "date": "2025-03", "result": "normal"}],
        "timeline": [{"date": "2025-01", "event": "regression"}, {"date": "2024-06", "event": "first seizure"}],
        "identifiers": ["王小雨", "MZ0012345"],
    }
    out = case_cmd._apply(_ns(d, _json.dumps(ops, ensure_ascii=False)))
    applied = out.result["applied"]
    assert applied["identifiers"] == {"protected": 2}
    assert "王小雨" not in _json.dumps(out.result, ensure_ascii=False)  # never echoed back
    data = case_mod.load(d)
    assert [m["relation"] for m in data["family"]["members"]] == ["mother", "brother"]
    assert [m["id"] for m in data["family"]["members"]] == ["f1", "f2"]
    assert data["tests"][0]["id"] == "x1" and data["tests"][0]["result"] == "normal"
    assert [e["event"] for e in data["timeline"]] == ["first seizure", "regression"]  # sorted by date
    assert set(data["privacy"]["identifiers"]) == {"王小雨", "MZ0012345"}
    s = case_mod.summary(d)
    assert s["tests"][0]["type"] == "CMA" and s["family"][1]["affected"] is True and s["identifiers"] == 2

    # removal by the new kinds, and an id is never re-used after removal
    case_cmd._apply(_ns(d, _json.dumps({"remove": [{"kind": "relative", "id": "f2"}, {"kind": "test", "id": "x1"}]})))
    case_cmd._apply(_ns(d, _json.dumps({"family": [{"relation": "sister"}], "tests": [{"type": "exome"}]})))
    data = case_mod.load(d)
    assert [m["id"] for m in data["family"]["members"]] == ["f1", "f3"]
    assert [t["id"] for t in data["tests"]] == ["x2"]


def test_cp1_11_bad_family_op_writes_nothing(tmp_path):
    import json as _json

    import pytest

    from zebra import case as case_mod
    from zebra.commands import case as case_cmd
    from zebra.core import UsageError

    d = str(tmp_path / "c")
    case_mod.init(d, title="t")
    with pytest.raises(UsageError):
        case_cmd._apply(_ns(d, _json.dumps({"tests": [{"type": "CMA"}], "family": [{"relation": "neighbour"}]})))
    with pytest.raises(UsageError):
        case_cmd._apply(_ns(d, _json.dumps({"tests": [{"result": "normal"}]})))
    assert case_mod.load(d)["tests"] == []


def test_c_p1_1_bad_ledger_lines_never_break_the_case(tmp_path):
    """C-P1-1: a non-object row or bytes cut inside a character must not stop reading or numbering."""
    from zebra import case as case_mod

    d = str(tmp_path / "c")
    case_mod.init(d, title="t")
    ledger = tmp_path / "c" / "evidence" / "ledger.jsonl"
    ledger.parent.mkdir(parents=True, exist_ok=True)
    with open(ledger, "wb") as fh:
        fh.write(b'{"eid": "E4", "db": "HPO"}\n')
        fh.write(b'null\n[1, 2]\n42\n')
        fh.write('{"eid": "E5", "db": "中文'.encode("utf-8")[:-1] + b"\n")  # cut inside a character
    rows = case_mod.read_ledger(d)
    assert [r["eid"] for r in rows] == ["E4"]
    assert case_mod.append_ledger(d, "x", {}, [{"db": "HPO"}]) == ["E6"]
    assert case_mod.summary(d)["evidence_count"] == 2


def _ns(case_dir, ops_json):
    return type("A", (), {"ops": ops_json, "dir": case_dir, "case": None})()


def test_c_p1_6_readding_an_excluded_term_keeps_it_excluded(tmp_path):
    from zebra import case as case_mod
    from zebra.cli import main

    d = str(tmp_path / "c")
    case_mod.init(d, title="t")
    assert main(["--json", "case", "add-hpo", "HP:0001250", "--label", "Seizure", "--status", "excluded",
                 "--note", "EEG normal", "--case", d]) == 0
    assert main(["--json", "case", "add-hpo", "HP:0001250", "--label", "Seizure", "--source", "records/neuro.pdf",
                 "--case", d]) == 0
    p = case_mod.load(d)["phenotypes"][0]
    assert p["status"] == "excluded" and p["note"] == "EEG normal" and p["source"] == "records/neuro.pdf"
    assert main(["--json", "case", "add-hpo", "HP:0001263", "--label", "Global developmental delay", "--case", d]) == 0
    assert case_mod.load(d)["phenotypes"][1]["status"] == "present"


def test_c_p1_5_ledger_tail_keeps_the_newest_rows_and_the_true_count(tmp_path, capsys):
    import json as _json

    from zebra import case as case_mod
    from zebra.cli import main

    d = str(tmp_path / "c")
    case_mod.init(d, title="t")
    case_mod.append_ledger(d, "x", {}, [{"db": "HPO", "record": f"r{i}", "url": "https://x.org/" + "a" * 200}
                                        for i in range(1000)])
    capsys.readouterr()
    assert main(["--json", "case", "ledger", "--case", d, "--tail", "25"]) == 0
    out = _json.loads(capsys.readouterr().out)
    assert out["result"]["total"] == 1000
    assert [r["eid"] for r in out["result"]["rows"]][-1] == "E1000"
    assert len(out["result"]["rows"]) == 25


def test_c_p1_4_a_bad_cached_body_is_evicted_and_not_found_lives_a_day(tmp_path, monkeypatch):
    import json as _json
    import time as _time

    from zebra import http
    from zebra.sources import validated_json

    monkeypatch.setenv("ZEBRA_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.delenv("ZEBRA_NO_CACHE", raising=False)
    url = "https://rest.genenames.org/fetch/symbol/SCN1A"
    key = f"GET {url} "
    path = http._cache_path(key)
    path.parent.mkdir(parents=True, exist_ok=True)

    def store(status, text, age):
        path.write_text(_json.dumps({"stored": _time.time() - age, "status": status, "text": text,
                                     "retrieved_at": "2026-10-01T00:00:00+00:00"}), "utf-8")

    # an error document cached with HTTP 200: read once, rejected, and gone
    store(200, '{"error": "Service temporarily unavailable"}', 5 * 86400)
    resp = http.request(url, source="HGNC", cache_ttl=30 * 86400)
    assert resp.cached
    try:
        validated_json(resp, "HGNC", require="response")
        raise AssertionError("an error body was accepted")
    except http.SourceError:
        pass
    assert not path.exists()

    # an HTML 200 body in the cache is never served
    def no_network(*a, **k):
        raise http.SourceError("HGNC", url, None, "network error: offline")

    monkeypatch.setattr(http, "_opener", lambda: type("O", (), {"open": staticmethod(no_network)})())
    monkeypatch.setattr(http, "_wait_for_retry", lambda s: False)
    store(200, "<html>maintenance</html>", 60)
    try:
        http.request(url, source="HGNC", cache_ttl=30 * 86400, retries=0)
        raise AssertionError("an HTML cached body was served")
    except http.SourceError:
        pass

    # an accepted 404 is served for a day, not for the record's 30
    store(404, '{"detail": "not found"}', 2 * 86400)
    try:
        http.request(url, source="HGNC", cache_ttl=30 * 86400, ok_statuses=(200, 404), retries=0)
        raise AssertionError("a two-day-old not-found was served from a 30-day cache")
    except http.SourceError:
        pass
    store(404, '{"detail": "not found"}', 3600)
    assert http.request(url, source="HGNC", cache_ttl=30 * 86400, ok_statuses=(200, 404)).status == 404
