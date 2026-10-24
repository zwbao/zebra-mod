"""zebra.s2f: parsers on real captured responses, model gating, and live calls.

Fixtures under tests/fixtures/spliceai/ are real lookup responses, and those under
tests/fixtures/s2f/ real run summaries of the heavy models, trimmed to the fields read.
"""

from __future__ import annotations

import json
import os
import struct
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
    out = s2f.predict("17-7674291-C-T", models=["spliceai", "pangolin"])
    assert calls == ["spliceai", "pangolin"]
    rows = {r["model"]: r for r in out.result["models"]}
    assert rows["spliceai"]["headline"]["value"] != rows["pangolin"]["headline"]["value"]
    text = s2f.render(out.result)
    assert "never average" in out.result["how_to_combine"]
    assert "spliceai" in text and "pangolin" in text
    assert "consensus" not in json.dumps(out.result).lower()


# ------------------------------------------------------- summaries -> rows


def test_headline_value_top_level_and_per_head():
    assert s2f._headline_value({"score": -6.77}, "score") == -6.77
    got = s2f._headline_value({"head_stats": {"RNA_SEQ": {"log2fc": -4e-05}}}, "log2fc")
    assert got == {"head": "RNA_SEQ", "value": -4e-05}
    assert s2f._headline_value({}, None) is None
    assert s2f._headline_value({"a": 1}, "missing") is None


@pytest.mark.parametrize(
    "fixture,model,name,value",
    [("s2f/gpn_msa_variant_summary.json", "gpn_msa", "score", -6.77),
     ("s2f/evo2_score_summary.json", "evo2", "delta_loglik", -40.10391114001686)],
)
def test_summary_becomes_a_row_with_a_receipt(fixture, model, name, value, tmp_path):
    row = s2f._from_summary(model, load(fixture), tmp_path, 1.0, 240.0)
    assert row["status"] == "ran"
    assert row["headline"] == {"name": name, "value": value}
    receipt = Path(row["receipt"])
    assert receipt.name == "receipt.json" and receipt.is_file() and row["url"].startswith("file://")
    body = json.loads(receipt.read_text())
    assert body["run_id"] == row["run_id"] and body["summary"]["headline_score"] == name
    assert row["units"] and row["definitions"] and name in row["definitions"]
    assert row["declared_claim_level"] in ("molecular", "cellular")
    assert row["model_id"]


def test_alphagenome_row_reports_the_head_it_answered_for(tmp_path):
    row = s2f._from_summary("alphagenome", load("s2f/alphagenome_variant_summary.json"), tmp_path, 8.0, 240.0)
    assert row["status"] == "ran"
    assert row["headline"]["name"] == "log2fc"
    assert row["headline"]["value"]["head"] == "RNA_SEQ"
    assert row["declared_claim_level"] == "cellular"
    assert row["inputs"]["ontology"] == ["UBERON:0000955"]


def test_lookup_status_other_than_scored_is_flagged(tmp_path):
    summary = load("s2f/gpn_msa_variant_summary.json")
    summary["lookup_status"] = "absent"
    summary.pop("score", None)
    row = s2f._from_summary("gpn_msa", summary, tmp_path, 1.0, 240.0)
    assert row["lookup_status"] == "absent"
    assert "never read as 'no effect'" in row["lookup_status_note"]
    assert "absent" in row["reason"]


def test_timeout_is_reported_as_an_error_without_a_receipt(tmp_path):
    row = s2f._from_summary("alphagenome", {"timed_out": True}, tmp_path, 240.0, 240.0)
    assert row["status"] == "error" and "timed out after 240s" in row["reason"]
    assert "receipt" not in row and not any(tmp_path.iterdir())


def test_backend_failure_is_an_error_and_keeps_its_receipt(tmp_path):
    summary = load("s2f/evo2_score_summary.json")
    summary["status"] = "error"
    summary["error"] = "RuntimeError: endpoint refused the key"
    row = s2f._from_summary("evo2", summary, tmp_path, 3.0, 240.0)
    assert row["status"] == "error" and "endpoint refused" in row["reason"]
    assert Path(row["receipt"]).is_file()


# ------------------------------------------------------------ GPN-MSA table


def _bgzf_block(payload: bytes) -> bytes:
    """One BGZF block (gzip member with the BC extra field), as bgzip writes it."""
    import zlib as _z

    comp = _z.compressobj(6, _z.DEFLATED, -15)
    deflated = comp.compress(payload) + comp.flush()
    bsize = 12 + 6 + len(deflated) + 8 - 1
    head = b"\x1f\x8b\x08\x04" + b"\x00" * 4 + b"\x00\xff" + struct.pack("<H", 6) + b"BC" + struct.pack("<HH", 2, bsize)
    return head + deflated + struct.pack("<II", _z.crc32(payload) & 0xFFFFFFFF, len(payload))


BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def _table_and_index(lines):
    """A bgzipped table (chrom pos ref alt score) in one block, and its parsed tabix index."""
    payload = "".join(l + "\n" for l in lines).encode()
    table = _bgzf_block(payload) + BGZF_EOF
    names = b"12\x00"
    pos0 = int(lines[0].split("\t")[1]) - 1
    leaf = 4681 + (pos0 >> 14)
    raw = b"TBI\x01" + struct.pack("<8i", 1, 0, 1, 2, 2, ord("#"), 0, len(names)) + names
    raw += struct.pack("<i", 1) + struct.pack("<Ii", leaf, 1) + struct.pack("<QQ", 0, len(payload))
    n_intv = (pos0 >> 14) + 1
    raw += struct.pack("<i", n_intv) + struct.pack(f"<{n_intv}Q", *([0] * n_intv))
    return table, raw


def test_bgzf_and_tbi_readers_find_the_row(monkeypatch):
    from zebra import s2f_backends as B

    lines = ["12\t25245350\tC\tA\t-10.1", "12\t25245350\tC\tT\t-12.05", "12\t25245351\tC\tG\t-3.2"]
    table, raw = _table_and_index(lines)
    text, blocks = B.bgzf_decompress(table)
    assert text.decode().splitlines() == lines and blocks[0] == (0, 0)
    idx = B.parse_tbi(raw)
    assert idx["names"] == ["12"] and idx["col_beg"] == 2 and idx["meta"] == "#"
    monkeypatch.setattr(B, "_tbi", lambda url, timeout: idx)
    monkeypatch.setattr(B, "http_range", lambda url, a, b, timeout: table[a:b + 1])
    rows = B.tabix_rows("https://x.invalid/t.bgz", "12", 25245350)
    assert [r[3] for r in rows] == ["A", "T"]
    assert B.tabix_rows("https://x.invalid/t.bgz", "3", 1) == []
    got = B.gpn_msa_score("12", 25245350, "C", "T")
    assert got["status"] == "ok" and got["lookup_status"] == "scored" and got["score"] == -12.05
    assert got["self_check"]["variant"] == "KRAS G12D" and got["self_check"]["agrees"] is True
    assert B.gpn_msa_score("12", 25245350, "C", "G")["lookup_status"] == "present_but_null"
    assert B.gpn_msa_score("12", 25245352, "C", "G")["lookup_status"] == "absent"


def test_gpn_transport_failure_is_an_error_not_absent(monkeypatch):
    from zebra import s2f_backends as B

    def boom(url, timeout):
        raise B.RangeError("connection reset")

    monkeypatch.setattr(B, "_tbi", boom)
    got = B.gpn_msa_score("12", 25245350, "C", "T")
    assert got["status"] == "error" and got["lookup_status"] == "error"
    assert "not evidence that the variant is absent" in got["error"]


def test_reg2bins_covers_root_and_leaf():
    from zebra import s2f_backends as B

    bins = B.reg2bins(25245349, 25245350)
    assert bins[0] == 0 and 4681 + (25245349 >> 14) in bins and len(bins) == 6


# ------------------------------------------------------------------ Evo 2


def _npy(shape, values, descr="<f4"):
    header = repr({"descr": descr, "fortran_order": False, "shape": tuple(shape)}).encode()
    header += b" " * (-(len(header) + 11) % 16) + b"\n"
    if descr == "<f4":
        body = struct.pack(f"<{len(values)}f", *values)
    else:  # bfloat16 as 2-byte void: the upper half of each float32
        body = struct.pack(f"<{len(values)}H", *[struct.unpack("<I", struct.pack("<f", v))[0] >> 16 for v in values])
    return b"\x93NUMPY\x01\x00" + struct.pack("<H", len(header)) + header + body


def test_read_npy_float32_and_bfloat16():
    from zebra import s2f_backends as B

    assert B.read_npy(_npy((1, 2, 2), [1.0, -2.0, 0.5, 4.0])) == ((1, 2, 2), [1.0, -2.0, 0.5, 4.0])
    shape, vals = B.read_npy(_npy((2, 2), [1.0, -2.0, 0.5, 4.0], descr="|V2"))
    assert shape == (2, 2) and vals == [1.0, -2.0, 0.5, 4.0]


def test_token_logprobs_is_a_causal_log_softmax():
    import math as _m
    from zebra import s2f_backends as B

    seq = "ACG"
    flat = [0.0] * (3 * 512)
    lp = B.token_logprobs((1, 3, 512), flat, seq)
    assert lp[0] == 0.0 and lp[1] == pytest.approx(-_m.log(512)) and lp[2] == pytest.approx(-_m.log(512))


def _evo2_logits(seq):
    """Logits that favour the reference: the next base of REF scores 5, everything else 0."""
    ref = "ACGT" * 16
    flat = [0.0] * (len(seq) * 512)
    for i in range(len(seq) - 1):
        flat[i * 512 + ord(ref[i + 1])] = 5.0
    return (1, len(seq), 512), flat


def test_evo2_score_alt_is_less_likely_and_the_answering_model_is_reported(monkeypatch):
    from zebra import s2f_backends as B

    monkeypatch.setenv("NVCF_RUN_KEY", "k")
    monkeypatch.delenv("ZEBRA_EVO2_MODEL", raising=False)
    asked = []

    def forward(seq, model, key, timeout):
        asked.append(model)
        if model == "evo2-7b":
            raise RuntimeError("HTTP 503: degraded")
        return _evo2_logits(seq)

    monkeypatch.setattr(B, "evo2_forward", forward)
    ref = "ACGT" * 16
    # window 64 centred so that the variant sits at index 31: pos - 32 + 1 = start
    got = B.evo2_score("1", 1031, ref[31], "A" if ref[31] != "A" else "C", window=64,
                       reference=lambda c, a, b: ref)
    assert got["status"] == "ok" and got["delta_loglik"] < 0
    assert got["model"]["id"] == "evo2-40b" and got["fallback_from"] == "evo2-7b"
    assert asked[0] == "evo2-7b"
    step = next(s for s in got["preflight_steps"] if s["step"] == "model")
    assert step["ok"] is False and "answered by evo2-40b" in step["detail"]


def test_evo2_ref_mismatch_is_an_error(monkeypatch):
    from zebra import s2f_backends as B

    monkeypatch.setenv("NVCF_RUN_KEY", "k")
    ref = "ACGT" * 16
    wrong = "A" if ref[31] != "A" else "C"
    got = B.evo2_score("1", 1031, wrong, "G", window=64, reference=lambda c, a, b: ref)
    assert got["status"] == "error" and "REF mismatch" in got["error"]


# ------------------------------------------------------------- AlphaGenome


def _fake_python(tmp_path, body, sleep=0):
    exe = tmp_path / "py"
    exe.write_text("#!/bin/sh\n"
                   'if [ "$1" = "-c" ]; then exit 0; fi\n'
                   f"sleep {sleep}\n"
                   f"echo '{body}'\n")
    exe.chmod(0o755)
    return str(exe)


def test_alphagenome_runner_output_becomes_a_row(monkeypatch, tmp_path):
    from zebra import s2f_backends as B

    body = json.dumps({"status": "ok", "model": {"id": "alphagenome hosted API", "version": "0.9.0"},
                       "head_stats": {"RNA_SEQ": {"num_tracks": 2, "log2fc": -0.01, "local_max_abs_log2fc": -0.8,
                                                  "local_track": "brain", "local_offset_bp": 12}}})
    monkeypatch.setenv("ZEBRA_ALPHAGENOME_PYTHON", _fake_python(tmp_path, body))
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "k")
    monkeypatch.setattr(s2f, "_ref_mismatch", lambda v: None)
    assert B.alphagenome_importable()
    row = s2f.run_alphagenome(s2f.normalise_variant("2-166042334-G-A").result, "GRCh38", tmp_path / "ws", 30.0,
                              ontology="UBERON:0000955")
    assert row["status"] == "ran" and row["headline"] == {"name": "log2fc", "value": {"head": "RNA_SEQ", "value": -0.01}}
    assert row["head_stats"]["RNA_SEQ"]["local_max_abs_log2fc"] == -0.8
    assert row["inputs"]["ontology"] == ["UBERON:0000955"] and "clinical decision-making" in row["terms"]


def test_alphagenome_ref_mismatch_never_calls_the_api(monkeypatch, tmp_path):
    monkeypatch.setenv("ZEBRA_ALPHAGENOME_PYTHON", _fake_python(tmp_path, "{}"))
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "k")
    monkeypatch.setattr(s2f, "_ref_mismatch", lambda v: "C")
    row = s2f.run_alphagenome(s2f.normalise_variant("2-166042334-G-A").result, "GRCh38", tmp_path / "ws", 30.0,
                              ontology="UBERON:0000955")
    assert row["status"] == "error" and "genome has C" in row["reason"]


def test_alphagenome_without_the_client_is_not_run(monkeypatch, tmp_path):
    exe = tmp_path / "py"
    exe.write_text("#!/bin/sh\nexit 1\n")
    exe.chmod(0o755)
    monkeypatch.setenv("ZEBRA_ALPHAGENOME_PYTHON", str(exe))
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "k")
    row = s2f.run_alphagenome(s2f.normalise_variant("2-166042334-G-A").result, "GRCh38", tmp_path, 10.0,
                              ontology="UBERON:0000955")
    assert row["status"] == "not_run" and "pip install alphagenome" in row["reason"]


# -------------------------------------------------------- gating / defaults


def test_alphagenome_without_a_key_never_runs(monkeypatch, tmp_path):
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    row = s2f.run_alphagenome(s2f.normalise_variant("2-166042334-G-A").result, "GRCh38",
                              tmp_path, 10.0, ontology="UBERON:0000955")
    assert row["status"] == "not_run" and "ALPHAGENOME_API_KEY" in row["reason"]


def test_alphagenome_without_an_ontology_never_runs(monkeypatch, tmp_path):
    """Without a tissue the answer would be for an arbitrary one; silence would be wrong science."""
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "x")
    row = s2f.run_alphagenome(s2f.normalise_variant("2-166042334-G-A").result, "GRCh38",
                              tmp_path, 10.0, ontology=None)
    assert row["status"] == "not_run"
    assert "--ontology" in row["reason"] and "UBERON:0000955" in row["reason"]


def test_evo2_without_a_key_never_runs(monkeypatch, tmp_path):
    monkeypatch.delenv("NVCF_RUN_KEY", raising=False)
    monkeypatch.delenv("EVO2_API_KEY", raising=False)
    row = s2f.run_evo2(s2f.normalise_variant("2-166042334-G-A").result, "GRCh38", tmp_path, 10.0)
    assert row["status"] == "not_run" and "NVCF_RUN_KEY" in row["reason"]


def test_gpn_msa_refuses_grch37_and_indels(tmp_path):
    snv = s2f.normalise_variant("2-166042334-G-A").result
    assert "hg38" in s2f.run_gpn_msa(snv, "GRCh37", tmp_path, 10.0)["reason"]
    indel = s2f.normalise_variant("7-117559590-ATCT-A").result
    row = s2f.run_gpn_msa(indel, "GRCh38", tmp_path, 10.0)
    assert row["status"] == "not_run" and "SNVs only" in row["reason"]


def test_default_models_without_keys_are_the_splice_models_and_gpn_msa(monkeypatch):
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    monkeypatch.delenv("NVCF_RUN_KEY", raising=False)
    monkeypatch.delenv("EVO2_API_KEY", raising=False)
    v = s2f.normalise_variant("2-166042334-G-A").result
    chosen, skipped = s2f._select(None, v, "GRCh38", None)
    assert chosen == ["spliceai", "pangolin", "gpn_msa"]
    assert {s["model"] for s in skipped} == {"evo2", "alphagenome"}
    assert all(s["reason"] for s in skipped)
    indel = s2f.normalise_variant("7-117559590-ATCT-A").result
    assert s2f._select(None, indel, "GRCh38", None)[0] == ["spliceai", "pangolin"]


def test_explicit_models_are_honoured_and_validated():
    v = s2f.normalise_variant("2-166042334-G-A").result
    chosen, skipped = s2f._select(["evo2", "evo2", "spliceai"], v, "GRCh38", None)
    assert chosen == ["evo2", "spliceai"] and skipped == []
    with pytest.raises(UsageError):
        s2f._select(["alphafold"], v, "GRCh38", None)


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


def test_render_reports_not_run_rows_and_ceilings(monkeypatch):
    monkeypatch.setattr(s2f, "lookup", fake_lookup(load("spliceai/spliceai_tp53_acceptor_loss.json")))
    monkeypatch.delenv("NVCF_RUN_KEY", raising=False)
    monkeypatch.delenv("EVO2_API_KEY", raising=False)
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
    text = s2f.render(s2f.predict("17-7674291-C-T", models=["spliceai"]).result)
    assert "CC BY-NC" in text and "research use" in text


# ===================================================================== v0.2 (W6)
# C-P2-9: the deadline-sharing code had no test (mutants M6 and M7 survived), and a run cut by the
# deadline told the caller to "raise --timeout".

import time  # noqa: E402


def _slow_alphagenome(monkeypatch, tmp_path, seconds=30):
    """A child interpreter that answers the import check at once and then sleeps through the run."""
    marker = tmp_path / "ran"
    exe = tmp_path / "py"
    exe.write_text("#!/bin/sh\n"
                   'if [ "$1" = "-c" ]; then exit 0; fi\n'
                   f"touch '{marker}'\nsleep {seconds}\n")
    exe.chmod(0o755)
    monkeypatch.setenv("ZEBRA_ALPHAGENOME_PYTHON", str(exe))
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "fake")
    monkeypatch.setattr(s2f, "_ref_mismatch", lambda v: None)
    return marker


def _deadline(monkeypatch, total):
    start = time.monotonic()
    monkeypatch.setattr(s2f, "deadline_seconds", lambda: total - (time.monotonic() - start))


def test_c_p2_9_a_model_is_cut_to_what_is_left_of_the_deadline(monkeypatch, tmp_path):
    marker = _slow_alphagenome(monkeypatch, tmp_path)
    _deadline(monkeypatch, 6.0)
    t0 = time.monotonic()
    out = s2f.predict("2-166042334-G-A", models=["alphagenome"], ontology="UBERON:0000955", timeout=240,
                      workspace=str(tmp_path / "ws"))
    wall = time.monotonic() - t0
    row = out.result["models"][0]
    assert marker.exists(), "the model must have started"
    assert wall < 12, f"--timeout 240 must be cut to the deadline, took {wall:.1f}s"
    assert row["status"] == "error" and "overall time budget" in row["reason"]
    assert "raise --timeout" not in row["reason"]


def test_c_p2_9_a_model_that_cannot_fit_is_not_started(monkeypatch, tmp_path):
    marker = _slow_alphagenome(monkeypatch, tmp_path)
    _deadline(monkeypatch, 3.0)
    out = s2f.predict("2-166042334-G-A", models=["alphagenome"], ontology="UBERON:0000955", timeout=240,
                      workspace=str(tmp_path / "ws"))
    row = out.result["models"][0]
    assert row["status"] == "not_run" and "ran out before this model started" in row["reason"]
    assert not marker.exists()


def test_c_p2_9_a_plain_timeout_still_says_raise_timeout(monkeypatch, tmp_path):
    marker = _slow_alphagenome(monkeypatch, tmp_path)
    monkeypatch.setattr(s2f, "deadline_seconds", lambda: None)
    out = s2f.predict("2-166042334-G-A", models=["alphagenome"], ontology="UBERON:0000955", timeout=1,
                      workspace=str(tmp_path / "ws"))
    row = out.result["models"][0]
    assert marker.exists() and row["status"] == "error" and "raise --timeout" in row["reason"]


@pytest.mark.live
def test_live_gpn_msa_published_table_matches_the_published_value():
    from zebra import s2f_backends as B

    got = B.gpn_msa_score("12", 25245350, "C", "T")
    assert got["lookup_status"] == "scored" and got["self_check"]["agrees"] is True
