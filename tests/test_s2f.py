"""zebra.s2f: parsers on real captured responses, model gating, and live calls.

Fixtures under tests/fixtures/spliceai/ and tests/fixtures/s2f/ are real responses
and real `s2f run --json` receipts, trimmed to the fields the parsers read.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from zebra import s2f
from zebra.core import UsageError
from zebra.sources import ensembl

FIXTURES = Path(__file__).parent / "fixtures"


def load(name):
    return json.loads((FIXTURES / name).read_text("utf-8"))


def fake_lookup(payload):
    def _lookup(tool, v, **kw):
        return dict(payload), {"db": tool, "record": "fixture", "url": "https://example.invalid/fixture"}

    return _lookup


# ------------------------------------------------------------ normalisation


def test_vcf_input_is_taken_as_given():
    out = s2f.normalise_variant("7-117559590-ATCT-A")
    assert (out.result["chrom"], out.result["pos"], out.result["ref"], out.result["alt"]) == ("7", 117559590, "ATCT", "A")
    assert out.result["kind"] == "vcf" and out.sources == []


@pytest.mark.parametrize("spelling", ["chr7-117559590-ATCT-A", "7:117559590:ATCT:A", "7 117559590 ATCT A"])
def test_vcf_spellings(spelling):
    out = s2f.normalise_variant(spelling)
    assert out.result["chrom"] == "7" and out.result["pos"] == 117559590


def test_hgvs_on_minus_strand_gene_uses_vcf_string_not_transcript_alleles(monkeypatch):
    """SCN1A NM_001165963.4:c.2134C>T is 2-166042334-G-A on the forward strand.

    VEP would report allele_string "C/T" with strand -1 for this variant; reading that
    as the VCF alleles gives a REF the reference genome does not have. Regression test
    on the real variant_recoder response.
    """
    captured = {}

    def fake_recode(variant, assembly="GRCh38"):
        captured["variant"] = variant
        from zebra.core import Outcome

        return Outcome(load("s2f/recoder_scn1a_minus_strand.json"),
                       sources=[{"db": "Ensembl variant_recoder", "record": variant, "url": "https://rest.ensembl.org/"}])

    monkeypatch.setattr(ensembl, "recode", fake_recode)
    out = s2f.normalise_variant("NM_001165963.4:c.2134C>T")
    v = out.result
    assert captured["variant"] == "NM_001165963.4:c.2134C>T"
    assert (v["chrom"], v["pos"], v["ref"], v["alt"]) == ("2", 166042334, "G", "A")
    assert v["ref"] != "C" and v["alt"] != "T", "transcript-strand alleles must not reach the VCF form"
    assert v["hgvs_g"] == "NC_000002.12:g.166042334G>A"
    assert "rs794726730" in v["ids"]
    assert out.sources, "the recoder call must be recorded as a source"


def test_lrg_spelling_is_never_chosen(monkeypatch):
    from zebra.core import Outcome

    payload = [{"T": {"vcf_string": ["LRG_8-36306-C-T", "2-166042334-G-A"], "id": [], "hgvsg": []}}]
    monkeypatch.setattr(ensembl, "recode", lambda *a, **k: Outcome(payload))
    v = s2f.normalise_variant("NM_001165963.4:c.2134C>T").result
    assert (v["chrom"], v["pos"], v["ref"], v["alt"]) == ("2", 166042334, "G", "A")


def test_several_alternate_alleles_warn_and_take_the_first(monkeypatch):
    from zebra.core import Outcome

    payload = [{"T": {"vcf_string": ["7-100-A-T"]}, "G": {"vcf_string": ["7-100-A-G"]}}]
    monkeypatch.setattr(ensembl, "recode", lambda *a, **k: Outcome(payload))
    out = s2f.normalise_variant("rs1")
    assert out.result["alt"] in ("T", "G")
    assert any("alternate alleles" in w for w in out.warnings)


def test_unreadable_input_is_a_usage_error():
    with pytest.raises(UsageError):
        s2f.normalise_variant("the gene is broken")
    with pytest.raises(UsageError):
        s2f.normalise_variant("")


def test_is_snv_and_vcf_key():
    v = s2f.normalise_variant("2-166042334-G-A").result
    assert s2f.is_snv(v) and s2f.vcf_key(v) == "chr2-166042334-G-A"
    assert s2f.vcf_key(v, with_chr=False) == "2-166042334-G-A"
    assert not s2f.is_snv(s2f.normalise_variant("7-117559590-ATCT-A").result)


# --------------------------------------------------------------- SpliceAI


def test_spliceai_canonical_acceptor_loss(monkeypatch):
    monkeypatch.setattr(s2f, "lookup", fake_lookup(load("spliceai/spliceai_tp53_acceptor_loss.json")))
    v = s2f.normalise_variant("17-7674291-C-T").result
    row = s2f.run_spliceai(v, "GRCh38", 500)
    assert row["status"] == "ran"
    assert row["transcript"]["refseq"] == "NM_000546.6"
    assert row["transcript"]["priority"] == "MS" and row["transcript"]["priority_label"] == "MANE Select"
    head = row["headline"]
    assert head["score"] == "DS_AL" and head["value"] == pytest.approx(0.998)
    # DP_AL = -1 on this variant: the acceptor that is lost sits one base 5' of it
    assert head["position"] == 7674290
    assert row["units"].startswith("delta score")
    assert set(row["definitions"]) >= {"DS_AG", "DS_AL", "DS_DG", "DS_DL"}
    assert row["claim_ceiling"].startswith("molecular")
    labels = [ab["label"] for ab in row["aberrations"]]
    assert any("exon 7 skipping" in (l or "").lower() for l in labels)
    assert row["aberrations"][0]["frameshift"] is True


def test_spliceai_deep_intronic_donor_gain(monkeypatch):
    monkeypatch.setattr(s2f, "lookup", fake_lookup(load("spliceai/spliceai_cftr_deep_intronic.json")))
    v = s2f.normalise_variant("7-117639961-C-T").result
    row = s2f.run_spliceai(v, "GRCh38", 500)
    assert row["headline"]["score"] == "DS_DG" and row["headline"]["value"] == pytest.approx(0.162)
    assert row["headline"]["position"] == 117639959  # DP_DG = -2
    assert row["transcript"]["gene"] == "CFTR"
    # the lncRNA CFTR-AS2 transcript is in the response but must not outrank MANE Select
    assert row["transcripts"][0]["priority"] == "MS"
    assert any(t["biotype"] == "lncRNA" for t in row["transcripts"])


def test_spliceai_indel_keeps_the_spelling_that_was_scored(monkeypatch):
    monkeypatch.setattr(s2f, "lookup", fake_lookup(load("spliceai/spliceai_cftr_f508del.json")))
    v = s2f.normalise_variant("7-117559590-ATCT-A").result
    row = s2f.run_spliceai(v, "GRCh38", 500)
    assert row["scored_as"] == {"chrom": "7", "pos": 117559590, "ref": "ATCT", "alt": "A"}
    assert row["headline"]["value"] == pytest.approx(0.007)


def test_spliceai_input_error_becomes_status_error(monkeypatch):
    monkeypatch.setattr(s2f, "lookup", fake_lookup(load("spliceai/spliceai_ref_mismatch_error.json")))
    v = s2f.normalise_variant("7-117559590-G-A").result
    row = s2f.run_spliceai(v, "GRCh38", 500)
    assert row["status"] == "error" and row["input_error"] is True
    assert "unexpected reference allele" in row["reason"]
    assert "headline" not in row


def test_spliceai_no_scores_is_an_error_not_a_zero(monkeypatch):
    monkeypatch.setattr(s2f, "lookup", fake_lookup(load("spliceai/spliceai_no_scores_error.json")))
    v = s2f.normalise_variant("1-1000-A-G").result
    row = s2f.run_spliceai(v, "GRCh38", 500)
    assert row["status"] == "error"
    assert "did not return any scores" in row["reason"]


# --------------------------------------------------------------- Pangolin


def test_pangolin_keeps_the_sign_of_a_splice_loss(monkeypatch):
    monkeypatch.setattr(s2f, "lookup", fake_lookup(load("spliceai/pangolin_tp53_acceptor_loss.json")))
    v = s2f.normalise_variant("17-7674291-C-T").result
    row = s2f.run_pangolin(v, "GRCh38", 500)
    assert row["status"] == "ran"
    head = row["headline"]
    assert head["score"] == "DS_SL" and head["value"] == pytest.approx(-0.893)
    assert head["position"] == 7674290
    assert "P(splice)" in row["units"]
    assert row["transcript"]["refseq"] == "NM_000546.6"


def test_pangolin_splice_gain_on_the_deep_intronic_variant(monkeypatch):
    monkeypatch.setattr(s2f, "lookup", fake_lookup(load("spliceai/pangolin_cftr_deep_intronic.json")))
    v = s2f.normalise_variant("7-117639961-C-T").result
    row = s2f.run_pangolin(v, "GRCh38", 500)
    assert row["headline"]["score"] == "DS_SG" and row["headline"]["value"] == pytest.approx(0.333)
    assert row["transcript"]["gene"] == "CFTR" and row["transcripts"][0]["priority"] == "MS"


def test_pangolin_and_spliceai_are_never_merged(monkeypatch):
    calls = []

    def _lookup(tool, v, **kw):
        calls.append(tool)
        name = "spliceai/spliceai_tp53_acceptor_loss.json" if tool == "spliceai" else "spliceai/pangolin_tp53_acceptor_loss.json"
        return load(name), {"db": tool, "record": "fixture", "url": "https://example.invalid/f"}

    monkeypatch.setattr(s2f, "lookup", _lookup)
    monkeypatch.setattr(s2f, "s2f_bin", lambda: None)
    out = s2f.predict("17-7674291-C-T", models=["spliceai", "pangolin"])
    assert calls == ["spliceai", "pangolin"]
    rows = {r["model"]: r for r in out.result["models"]}
    assert rows["spliceai"]["headline"]["value"] != rows["pangolin"]["headline"]["value"]
    text = s2f.render(out.result)
    assert "never average" in out.result["how_to_combine"]
    assert "spliceai" in text and "pangolin" in text
    assert "consensus" not in json.dumps(out.result).lower()


# -------------------------------------------------------------- s2f CLI


def test_headline_value_top_level_and_per_head():
    assert s2f._headline_value({"score": -6.77}, "score") == -6.77
    got = s2f._headline_value({"head_stats": {"RNA_SEQ": {"log2fc": -4e-05}}}, "log2fc")
    assert got == {"head": "RNA_SEQ", "value": -4e-05}
    assert s2f._headline_value({}, None) is None
    assert s2f._headline_value({"a": 1}, "missing") is None


@pytest.mark.parametrize(
    "fixture,model,name,value",
    [("s2f/gpn_msa_variant_manifest.json", "gpn_msa", "score", -6.77),
     ("s2f/evo2_score_manifest.json", "evo2", "delta_loglik", -40.10391114001686)],
)
def test_receipt_parsing(fixture, model, name, value):
    call = {"manifest": load(fixture), "exit_code": 0, "wall_s": 1.0, "timed_out": False,
            "stdout_tail": "", "stderr_tail": "", "command": "s2f run ..."}
    row = s2f._from_receipt(model, call, Path("/tmp/ws"), 240.0)
    assert row["status"] == "ran"
    assert row["headline"] == {"name": name, "value": value}
    assert row["receipt"].endswith("manifest.json") and row["url"].startswith("file://")
    assert row["units"] and row["definitions"] and name in row["definitions"]
    assert row["declared_claim_level"] in ("molecular", "cellular")
    assert row["model_id"]


def test_alphagenome_receipt_reports_the_head_it_answered_for():
    call = {"manifest": load("s2f/alphagenome_variant_manifest.json"), "exit_code": 0, "wall_s": 8.0,
            "timed_out": False, "stdout_tail": "", "stderr_tail": "", "command": "s2f run ..."}
    row = s2f._from_receipt("alphagenome", call, Path("/tmp/ws"), 240.0)
    assert row["status"] == "ran"
    assert row["headline"]["name"] == "log2fc"
    assert row["headline"]["value"]["head"] == "RNA_SEQ"
    assert row["declared_claim_level"] == "cellular"
    assert row["inputs"]["ontology"] == ["UBERON:0000955"]


def test_receipt_lookup_status_other_than_scored_is_flagged():
    manifest = load("s2f/gpn_msa_variant_manifest.json")
    manifest["summary"]["lookup_status"] = "not_in_table"
    manifest["summary"].pop("score", None)
    call = {"manifest": manifest, "exit_code": 0, "wall_s": 1.0, "timed_out": False,
            "stdout_tail": "", "stderr_tail": "", "command": "c"}
    row = s2f._from_receipt("gpn_msa", call, Path("/tmp/ws"), 240.0)
    assert row["lookup_status"] == "not_in_table"
    assert "never read as 'no effect'" in row["lookup_status_note"]


def test_timeout_is_reported_as_an_error():
    call = {"manifest": None, "exit_code": -9, "wall_s": 240.0, "timed_out": True,
            "stdout_tail": "", "stderr_tail": "", "command": "c"}
    row = s2f._from_receipt("evo2", call, Path("/tmp/ws"), 240.0)
    assert row["status"] == "error" and "timed out after 240s" in row["reason"]


def test_missing_receipt_is_an_error_with_the_tail():
    call = {"manifest": None, "exit_code": 2, "wall_s": 0.4, "timed_out": False,
            "stdout_tail": "stack 's2f-core' is not installed", "stderr_tail": "", "command": "c"}
    row = s2f._from_receipt("gpn_msa", call, Path("/tmp/ws"), 240.0)
    assert row["status"] == "error" and "not installed" in row["reason"]


def test_script_level_failure_in_the_receipt_is_an_error():
    manifest = load("s2f/evo2_score_manifest.json")
    manifest["summary"]["status"] = "error"
    manifest["summary"]["error"] = "RuntimeError: endpoint refused the key"
    call = {"manifest": manifest, "exit_code": 1, "wall_s": 3.0, "timed_out": False,
            "stdout_tail": "", "stderr_tail": "", "command": "c"}
    row = s2f._from_receipt("evo2", call, Path("/tmp/ws"), 240.0)
    assert row["status"] == "error" and "endpoint refused" in row["reason"]


# -------------------------------------------------------- gating / defaults


def test_no_s2f_cli_is_not_run_with_the_install_hint(monkeypatch):
    monkeypatch.setattr(s2f, "s2f_bin", lambda: None)
    for fn in (s2f.run_gpn_msa, s2f.run_evo2):
        row = fn(s2f.normalise_variant("2-166042334-G-A").result, "GRCh38", None, Path("/tmp/ws"), 10.0)
        assert row["status"] == "not_run" and "s2f-penguin" in row["reason"]


def test_alphagenome_without_a_key_never_runs(monkeypatch, tmp_path):
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    row = s2f.run_alphagenome(s2f.normalise_variant("2-166042334-G-A").result, "GRCh38", "/bin/false",
                              tmp_path, 10.0, ontology="UBERON:0000955")
    assert row["status"] == "not_run" and "ALPHAGENOME_API_KEY" in row["reason"]


def test_alphagenome_without_an_ontology_never_runs(monkeypatch, tmp_path):
    """The family script's default ontology is transverse colon; silence would be wrong science."""
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "x")
    row = s2f.run_alphagenome(s2f.normalise_variant("2-166042334-G-A").result, "GRCh38", "/bin/false",
                              tmp_path, 10.0, ontology=None)
    assert row["status"] == "not_run"
    assert "--ontology" in row["reason"] and "UBERON:0001157" in row["reason"]


def test_evo2_without_a_key_never_runs(monkeypatch, tmp_path):
    monkeypatch.delenv("NVCF_RUN_KEY", raising=False)
    monkeypatch.delenv("EVO2_API_KEY", raising=False)
    row = s2f.run_evo2(s2f.normalise_variant("2-166042334-G-A").result, "GRCh38", "/bin/false", tmp_path, 10.0)
    assert row["status"] == "not_run" and "NVCF_RUN_KEY" in row["reason"]


def test_gpn_msa_refuses_grch37_and_indels(tmp_path):
    snv = s2f.normalise_variant("2-166042334-G-A").result
    assert "hg38" in s2f.run_gpn_msa(snv, "GRCh37", "/bin/false", tmp_path, 10.0)["reason"]
    indel = s2f.normalise_variant("7-117559590-ATCT-A").result
    row = s2f.run_gpn_msa(indel, "GRCh38", "/bin/false", tmp_path, 10.0)
    assert row["status"] == "not_run" and "SNVs only" in row["reason"]


def test_default_models_without_the_cli_are_the_two_splice_models(monkeypatch):
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    monkeypatch.delenv("NVCF_RUN_KEY", raising=False)
    monkeypatch.delenv("EVO2_API_KEY", raising=False)
    v = s2f.normalise_variant("2-166042334-G-A").result
    chosen, skipped = s2f._select(None, v, "GRCh38", None, None)
    assert chosen == ["spliceai", "pangolin"]
    assert {s["model"] for s in skipped} == {"gpn_msa", "evo2", "alphagenome"}
    assert all(s["reason"] for s in skipped)


def test_explicit_models_are_honoured_and_validated():
    v = s2f.normalise_variant("2-166042334-G-A").result
    chosen, skipped = s2f._select(["evo2", "evo2", "spliceai"], v, "GRCh38", None, None)
    assert chosen == ["evo2", "spliceai"] and skipped == []
    with pytest.raises(UsageError):
        s2f._select(["alphafold"], v, "GRCh38", None, None)


def test_distance_and_assembly_are_validated():
    with pytest.raises(UsageError):
        s2f.predict("2-166042334-G-A", assembly="hg38")
    with pytest.raises(UsageError):
        s2f.predict("2-166042334-G-A", distance=0)
    with pytest.raises(UsageError):
        s2f.predict("2-166042334-G-A", distance=10001)


def test_an_unexpected_error_in_one_model_does_not_sink_the_others(monkeypatch):
    def boom(*a, **k):
        raise RuntimeError("socket closed")

    monkeypatch.setattr(s2f, "run_spliceai", boom)
    monkeypatch.setattr(s2f, "run_pangolin", lambda *a, **k: {"model": "pangolin", "status": "ran",
                                                              "headline": {"name": "x", "value": 1.0},
                                                              "claim_ceiling": "molecular"})
    monkeypatch.setattr(s2f, "s2f_bin", lambda: None)
    out = s2f.predict("17-7674291-C-T", models=["spliceai", "pangolin"])
    rows = {r["model"]: r for r in out.result["models"]}
    assert rows["spliceai"]["status"] == "error" and "socket closed" in rows["spliceai"]["reason"]
    assert rows["pangolin"]["status"] == "ran"
    assert any("spliceai error" in w for w in out.warnings)


def test_workspace_goes_into_the_case_when_there_is_one(tmp_path, monkeypatch):
    monkeypatch.delenv("ZEBRA_CASE", raising=False)
    assert s2f.workspace_dir() == Path(os.path.expanduser("~")) / ".cache" / "zebra-mod" / "s2f-runs"
    assert s2f.workspace_dir(case=str(tmp_path)) == tmp_path.resolve() / "evidence" / "s2f"
    monkeypatch.setenv("ZEBRA_CASE", str(tmp_path))
    assert s2f.workspace_dir() == tmp_path.resolve() / "evidence" / "s2f"
    assert s2f.workspace_dir(workspace=str(tmp_path / "ws")) == tmp_path / "ws"


def test_s2f_bin_prefers_the_env_variable(monkeypatch, tmp_path):
    fake = tmp_path / "s2f"
    fake.write_text("#!/bin/sh\n")
    fake.chmod(0o755)
    monkeypatch.setenv("S2F_BIN", str(fake))
    assert s2f.s2f_bin() == str(fake)
    monkeypatch.setenv("S2F_BIN", str(tmp_path / "missing"))
    assert s2f.s2f_bin() != str(tmp_path / "missing")


def test_render_reports_not_run_rows_and_ceilings(monkeypatch):
    monkeypatch.setattr(s2f, "lookup", fake_lookup(load("spliceai/spliceai_tp53_acceptor_loss.json")))
    monkeypatch.setattr(s2f, "s2f_bin", lambda: None)
    out = s2f.predict("17-7674291-C-T", models=["spliceai", "evo2"])
    text = s2f.render(out.result)
    assert "evo2" in text and "not_run" in text
    assert "claim ceilings:" in text and "clinical" in text
    assert "RNA from a tissue" in text


# --------------------------------------------------------------- live


@pytest.mark.live
def test_live_spliceai_and_pangolin_on_a_deep_intronic_cftr_variant():
    out = s2f.predict("NM_000492.4:c.3718-2477C>T", models=["spliceai", "pangolin"])
    v = out.result["variant"]
    assert (v["chrom"], v["pos"], v["ref"], v["alt"]) == ("7", 117639961, "C", "T")
    rows = {r["model"]: r for r in out.result["models"]}
    assert rows["spliceai"]["status"] == "ran" and rows["pangolin"]["status"] == "ran"
    assert rows["spliceai"]["transcript"]["gene"] == "CFTR"
    assert rows["spliceai"]["headline"]["score"] == "DS_DG"
    assert rows["pangolin"]["headline"]["score"] == "DS_SG"
    assert all(src["url"] for src in out.sources)


@pytest.mark.live
def test_live_hgvs_on_scn1a_resolves_to_the_forward_strand():
    out = s2f.normalise_variant("NM_001165963.4:c.2134C>T")
    v = out.result
    assert (v["chrom"], v["pos"], v["ref"], v["alt"]) == ("2", 166042334, "G", "A")
    seq = ensembl.sequence(v["chrom"], v["pos"], v["pos"]).result
    assert seq == "G", "the forward reference base must be the VCF REF"


@pytest.mark.live
def test_live_spliceai_on_a_canonical_acceptor_site():
    out = s2f.predict("17-7674291-C-T", models=["spliceai"])
    row = out.result["models"][0]
    assert row["status"] == "ran"
    assert row["headline"]["score"] == "DS_AL" and row["headline"]["value"] > 0.9
    assert row["transcript"]["refseq"] == "NM_000546.6"


@pytest.mark.live
def test_live_grch37_is_served_too():
    out = s2f.predict("7-117199644-ATCT-A", assembly="GRCh37", models=["spliceai"])
    row = out.result["models"][0]
    assert row["status"] == "ran" and row["transcript"]["gene"] == "CFTR"


def test_a_cached_rate_limit_is_refetched_not_served(monkeypatch):
    """A rate limit is about this minute; an input error is about this variant."""
    from zebra import http as zhttp

    calls = []

    class Resp:
        def __init__(self, status, body, cached):
            self.status, self.text, self.url, self.cached = status, json.dumps(body), "https://x.invalid/", cached
            self.retrieved_at = "2026-10-05T00:00:00+00:00"

        def json(self):
            return json.loads(self.text)

    def fake_get_json(url, source=None, params=None, cache_ttl=None, **kw):
        calls.append(cache_ttl)
        if len(calls) == 1:
            return Resp(429, {"error": "Rate limit exceeded. This server only supports interactive use."}, True)
        return Resp(200, {"scores": [], "pos": 1}, False)

    monkeypatch.setattr(s2f, "get_json", fake_get_json)
    monkeypatch.setattr(s2f, "LOOKUP_INTERVAL", 0.0)
    payload, _rec = s2f.lookup("spliceai", s2f.normalise_variant("2-166042334-G-A").result)
    assert calls == [30 * 86400, 0], "the cached rate limit must be re-fetched live"
    assert "error" not in payload

    calls.clear()

    def input_error(url, source=None, params=None, cache_ttl=None, **kw):
        calls.append(cache_ttl)
        return Resp(200, {"error": "7-1-G has an unexpected reference allele.", "inputError": True}, True)

    monkeypatch.setattr(s2f, "get_json", input_error)
    payload, _rec = s2f.lookup("spliceai", s2f.normalise_variant("2-166042334-G-A").result)
    assert calls == [30 * 86400], "a deterministic input error may be served from the cache"
    assert payload["inputError"] is True
    assert zhttp  # the module the real call goes through


def test_lookup_refuses_an_assembly_it_has_no_service_for():
    with pytest.raises(UsageError):
        s2f.lookup("spliceai", s2f.normalise_variant("2-166042334-G-A").result, assembly="T2T")


def test_render_states_the_service_terms(monkeypatch):
    monkeypatch.setattr(s2f, "lookup", fake_lookup(load("spliceai/spliceai_tp53_acceptor_loss.json")))
    monkeypatch.setattr(s2f, "s2f_bin", lambda: None)
    text = s2f.render(s2f.predict("17-7674291-C-T", models=["spliceai"]).result)
    assert "CC BY-NC" in text and "research use" in text


# ===================================================================== v0.2 (W6)
# C-P2-9: the deadline-sharing code had no test (mutants M6 and M7 survived), and a run cut by the
# deadline told the caller to "raise --timeout".

import time  # noqa: E402


def _fake_s2f(tmp_path, seconds=30):
    marker = tmp_path / "ran"
    exe = tmp_path / "s2f"
    exe.write_text(f"#!/bin/sh\ntouch '{marker}'\nsleep {seconds}\n")
    exe.chmod(0o755)
    return str(exe), marker


def _deadline(monkeypatch, total):
    start = time.monotonic()
    monkeypatch.setattr(s2f, "deadline_seconds", lambda: total - (time.monotonic() - start))


def test_c_p2_9_a_model_is_cut_to_what_is_left_of_the_deadline(monkeypatch, tmp_path):
    exe, marker = _fake_s2f(tmp_path)
    monkeypatch.setenv("S2F_BIN", exe)
    monkeypatch.setenv("NVCF_RUN_KEY", "fake")
    _deadline(monkeypatch, 5.0)
    t0 = time.monotonic()
    out = s2f.predict("2-166042334-G-A", models=["evo2"], timeout=240, workspace=str(tmp_path / "ws"))
    wall = time.monotonic() - t0
    row = out.result["models"][0]
    assert marker.exists(), "the model must have started"
    assert wall < 10, f"--timeout 240 must be cut to the deadline, took {wall:.1f}s"
    assert row["status"] == "error" and "overall time budget" in row["reason"]
    assert "raise --timeout" not in row["reason"]


def test_c_p2_9_a_model_that_cannot_fit_is_not_started(monkeypatch, tmp_path):
    exe, marker = _fake_s2f(tmp_path)
    monkeypatch.setenv("S2F_BIN", exe)
    monkeypatch.setenv("NVCF_RUN_KEY", "fake")
    _deadline(monkeypatch, 3.0)
    out = s2f.predict("2-166042334-G-A", models=["evo2"], timeout=240, workspace=str(tmp_path / "ws"))
    row = out.result["models"][0]
    assert row["status"] == "not_run" and "ran out before this model started" in row["reason"]
    assert not marker.exists()


def test_c_p2_9_a_plain_timeout_still_says_raise_timeout(monkeypatch, tmp_path):
    exe, marker = _fake_s2f(tmp_path)
    monkeypatch.setattr(s2f, "deadline_seconds", lambda: None)
    call = s2f._s2f_call(exe, ["evo2"], tmp_path / "ws", 1.0)
    assert call["timed_out"] is True and call["deadline_limited"] is False
    row = s2f._from_receipt("evo2", call, tmp_path / "ws", 1.0)
    assert "raise --timeout" in row["reason"]


def test_review_b_deadline_label_follows_whether_the_deadline_cut_the_timeout(monkeypatch, tmp_path):
    exe, marker = _fake_s2f(tmp_path)
    # 2 s of deadline left after a run that --timeout (1 s) stopped: --timeout, not the deadline
    monkeypatch.setattr(s2f, "deadline_seconds", lambda: 2.0)
    call = s2f._s2f_call(exe, ["evo2"], tmp_path / "ws", 1.0)
    assert call["deadline_limited"] is False
    call = s2f._s2f_call(exe, ["evo2"], tmp_path / "ws", 1.0, deadline_cut=True)
    assert call["deadline_limited"] is True
