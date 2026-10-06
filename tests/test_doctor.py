"""zebra doctor: result shape, failure reporting, deadline, key secrecy; one live run."""

from __future__ import annotations

import json
import time

import pytest

from zebra.commands import doctor as D
from zebra.http import Response, SourceError

BODIES = {
    # trimmed real answers (2026-10-05) for the validators
    "ensembl": '{"ping":1}',
    "spliceai_ok": '{"variant": "chr8-140300616-T-G", "scores": [{"DS_AG": "0.045", "DS_AL": "0.827"}]}',
    "spliceai_bad": '{"variant": "chr8-140300616-A-G", "inputError": true, "error": "ref allele mismatch"}',
    "gnomad": '{"data":{"meta":{"clinvar_release_date":"2026-09-28"}}}',
    "litvar": '[{"_id":"litvar@rs6746030##","rsid":"rs6746030","gene":["SCN9A","SCN1A-AS1"]}]',
}


def _fake_request(status_by_name):
    def request(url, *, source, **kw):
        what = status_by_name.get(source, "ok")
        if what == "ok":
            body = BODIES["ensembl"] if "ensembl" in url else BODIES["gnomad"]
            return Response(url, 200, body, "2026-10-05T00:00:00+00:00", False)
        if what == "slow":
            time.sleep(5)
            return Response(url, 200, "{}", "2026-10-05T00:00:00+00:00", False)
        raise SourceError(source, url, what if isinstance(what, int) else None, "No server is available")
    return request


def test_validators_on_real_bodies():
    assert D._ensembl_ping(BODIES["ensembl"]) == ""
    assert "Cloud Run" in D._spliceai(BODIES["spliceai_ok"])
    with pytest.raises(ValueError):
        D._spliceai(BODIES["spliceai_bad"])
    D._json_has("data", "meta")(BODIES["gnomad"])
    D._json_has("0", "_id")(BODIES["litvar"])
    with pytest.raises(KeyError):
        D._json_has("data", "variant")(BODIES["gnomad"])
    assert "homepage" in D._html("<!DOCTYPE html><html lang='en'></html>")


def test_shape_and_failures(monkeypatch, tmp_path):
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(tmp_path / "nohpo"))
    monkeypatch.setattr(D, "request", _fake_request({"Ensembl REST GRCh37": 503, "PanelApp": "down"}))
    res = D.run_checks(case_dir=None, deadline=10)
    names = [c["name"] for c in res["checks"]]
    assert names[:5] == ["Python ≥ 3.9", "zebra", "cache dir writable", "HPO local files", "case"]
    assert all(set(c) - {"optional"} == {"name", "ok", "detail"} for c in res["checks"])
    by = {c["name"]: c for c in res["checks"]}
    assert len([n for n in names if n.startswith("network: ")]) == len(D.NETWORK) == 18
    assert by["network: Ensembl REST GRCh37"]["ok"] is False and "HTTP 503" in by["network: Ensembl REST GRCh37"]["detail"]
    assert by["network: PanelApp"]["ok"] is False and "unreachable" in by["network: PanelApp"]["detail"]
    assert by["network: Ensembl REST GRCh38"]["ok"] is True
    assert by["HPO local files"]["ok"] is False and "zebra hpo fetch" in by["HPO local files"]["detail"]
    assert by["case"]["ok"] is True and "no active case" in by["case"]["detail"]
    s = res["summary"]
    assert s["ok"] + s["failed"] + s["optional_missing"] == len(res["checks"])


def test_keys_are_never_printed(monkeypatch):
    secret = "sk-zebra-SECRET-1234567890"
    monkeypatch.setenv("ALPHAGENOME_API_KEY", secret)
    monkeypatch.setenv("EVO2_API_KEY", secret + "evo")
    monkeypatch.delenv("NVCF_RUN_KEY", raising=False)
    monkeypatch.delenv("OMIM_API_KEY", raising=False)
    res = D.run_checks(deadline=5, network=[])
    dump = json.dumps(res)
    assert secret not in dump and "SECRET" not in dump
    by = {c["name"]: c for c in res["checks"]}
    assert by["AlphaGenome key"]["ok"] and by["AlphaGenome key"]["detail"] == "ALPHAGENOME_API_KEY is set"
    assert by["Evo 2 key"]["ok"] and "EVO2_API_KEY is set" in by["Evo 2 key"]["detail"]
    assert by["OMIM API key"]["ok"] is False and "optional" in by["OMIM API key"]["detail"]


def test_deadline_bounds_a_hung_host(monkeypatch):
    monkeypatch.setattr(D, "request", _fake_request({"slowhost": "slow"}))
    t0 = time.monotonic()
    res = D.run_checks(deadline=1.0, network=[("slowhost", "https://example.org/", "GET", None, "application/json",
                                                D._json_has())])
    assert time.monotonic() - t0 < 3.0
    slow = [c for c in res["checks"] if c["name"] == "network: slowhost"][0]
    assert slow["ok"] is False and "timed out" in slow["detail"]


def test_case_and_hpo_present(monkeypatch, tmp_path):
    from zebra import case as case_mod

    hpo = tmp_path / "hpo"
    hpo.mkdir()
    (hpo / "hp.json").write_text("x" * 2000)
    (hpo / "genes_to_phenotype.txt").write_text("x" * 2000)
    (hpo / "phenotype.hpoa").write_text("#description: test\n#version: 2026-09-02\n"
                                        "#hpo-version: http://purl.obolibrary.org/obo/hp/releases/2026-09-01/hp.json\n"
                                        + "x" * 2000)
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(hpo))
    case_mod.init(str(tmp_path / "c1"), title="Test family")
    res = D.run_checks(case_dir=str(tmp_path / "c1"), deadline=5, network=[])
    by = {c["name"]: c for c in res["checks"]}
    assert by["HPO local files"]["ok"] and "annotations 2026-09-02, ontology 2026-09-01" in by["HPO local files"]["detail"]
    assert by["case"]["ok"] and "'Test family'" in by["case"]["detail"]
    bad = D.run_checks(case_dir=str(tmp_path / "missing"), deadline=5, network=[])
    assert {c["name"]: c for c in bad["checks"]}["case"]["ok"] is False


def test_cli_json_envelope(monkeypatch, capsys):
    from zebra.cli import main

    monkeypatch.setattr(D, "NETWORK", [])
    assert main(["--json", "doctor"]) == 0
    env = json.loads(capsys.readouterr().out)
    assert env["ok"] and env["command"] == "doctor" and env["sources"] == []
    assert isinstance(env["result"]["checks"], list) and env["result"]["checks"][0]["name"] == "Python ≥ 3.9"


@pytest.mark.live
def test_live_doctor_under_30s():
    t0 = time.monotonic()
    res = D.run_checks()
    took = time.monotonic() - t0
    assert took < 30
    net = [c for c in res["checks"] if c["name"].startswith("network: ")]
    assert len(net) == 18
    assert sum(c["ok"] for c in net) >= 12, [c for c in net if not c["ok"]]


def test_install_optional_items_are_not_failures(monkeypatch):
    """A family installing zebra-mod must not read 'not ok' for a research key they do not need."""
    from zebra.commands import doctor

    for _name, envs, _what in doctor.KEYS:
        for e in envs:
            monkeypatch.delenv(e, raising=False)
    res = doctor.run_checks(network=[], deadline=5)
    keys = [c for c in res["checks"] if c["name"] in {k[0] for k in doctor.KEYS}]
    assert keys and all(c.get("optional") and not c["ok"] for c in keys)
    assert res["summary"]["optional_missing"] >= len(keys)
