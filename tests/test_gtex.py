"""GTEx median expression per tissue (CP1-17: no expression layer, no way to pick an RNA-test tissue).

Offline tests parse captured GTEx v2 responses (tests/fixtures/gtex/, built by
tools/aso/fetch_fixtures.py). The live tests check the two calls still answer.
"""

from __future__ import annotations

import json
import os

import pytest

from zebra.core import UsageError
from zebra.http import Response, SourceError
from zebra.sources import gtex

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "gtex")


def text(name):
    path = os.path.join(FIX, name)
    if not os.path.exists(path):  # pragma: no cover - the fixtures are committed
        pytest.skip(f"missing fixture {path}; rebuild with tools/aso/fetch_fixtures.py")
    with open(path) as fh:
        return fh.read()


def fake_http(monkeypatch, mapping):
    """Route GTEx calls to captured bodies by a substring of the URL."""
    calls = []

    def fake(url, source, **kw):
        calls.append(url)
        params = kw.get("params") or {}
        if "medianGeneExpression" in url:
            # the module docstring makes this load-bearing: without datasetId the service
            # answers an empty list with HTTP 200
            assert params.get("datasetId") in gtex.DATASETS, params
            assert params.get("gencodeId"), params
        if "reference/gene" in url:
            assert params.get("gencodeVersion") in set(gtex.DATASETS.values()), params
        query = "&".join(f"{k}={v}" for k, v in params.items())
        for key, name in mapping.items():
            if key in url or key in query:
                body = name if isinstance(name, str) and name.startswith("{") else text(name)
                return Response(f"{url}?{query}", 200, body, "2026-10-06T00:00:00+00:00", True)
        raise AssertionError(f"unexpected GTEx call {url} {query}")

    monkeypatch.setattr(gtex, "get_json", fake)
    return calls


SETUP = {"geneId=CFTR": "reference_gene_cftr_v26.json",
         "medianGeneExpression": "median_expression_cftr_v8.json"}


# ---------------------------------------------------------------- gene resolution


def test_resolve_gene_reads_the_gencode_id_of_the_datasets_release(monkeypatch):
    fake_http(monkeypatch, SETUP)
    out = gtex.resolve_gene("CFTR", dataset="gtex_v8")
    assert out.result["found"] is True
    assert out.result["gencode_id"] == "ENSG00000001626.14"   # GENCODE v26, which gtex_v8 is built on
    assert out.result["gencode_version"] == "v26"
    assert out.sources and out.sources[0]["db"] == "GTEx reference/gene"


def test_resolve_gene_rejects_an_unknown_dataset():
    with pytest.raises(ValueError):
        gtex.resolve_gene("CFTR", dataset="gtex_v99")
    with pytest.raises(ValueError):
        gtex.resolve_gene("")


def test_a_gene_gtex_does_not_know_is_not_found_not_absent(monkeypatch):
    fake_http(monkeypatch, {"geneId=NOTAGENE123": "reference_gene_missing.json"})
    out = gtex.resolve_gene("NOTAGENE123")
    assert out.result["found"] is False and out.result["gencode_id"] is None
    assert any("not found in GTEx" in w and "not 'not expressed'" in w for w in out.warnings)


# ---------------------------------------------------------------- expression rows


def test_median_expression_carries_the_ontology_id_for_every_tissue(monkeypatch):
    fake_http(monkeypatch, SETUP)
    out = gtex.median_expression("ENSG00000001626.14", dataset="gtex_v8")
    rows = out.result
    assert len(rows) == 54
    by_name = {r["tissue"]: r for r in rows}
    assert by_name["Pancreas"]["uberon"] == "UBERON:0001150"
    assert by_name["Pancreas"]["ontology_is_uberon"] is True
    assert by_name["Pancreas"]["median_tpm"] == pytest.approx(65.79, abs=0.1)
    # the two cell lines carry EFO ids, and the flag says so instead of pretending
    assert by_name["Cells_Cultured_fibroblasts"]["uberon"] == "EFO:0002009"
    assert by_name["Cells_Cultured_fibroblasts"]["ontology_is_uberon"] is False
    for row in rows:
        assert row["unit"] == "TPM"
        assert row["band"] in {"expressed", "expressed_low", "trace", "not_detected"}


def test_expression_bands_are_monotonic():
    assert gtex._band(100.0)["band"] == "expressed"
    assert gtex._band(10.0)["band"] == "expressed"
    assert gtex._band(9.99)["band"] == "expressed_low"
    assert gtex._band(1.0)["band"] == "expressed_low"
    assert gtex._band(0.5)["band"] == "trace"
    assert gtex._band(0.0)["band"] == "not_detected"
    assert gtex._band(None)["band"] is None


def test_a_gencode_id_from_the_wrong_release_is_a_named_warning(monkeypatch):
    empty = json.dumps({"data": [], "paging_info": {"numberOfPages": 0, "totalNumberOfItems": 0}})
    fake_http(monkeypatch, {"medianGeneExpression": empty})
    out = gtex.median_expression("ENSG00000001626.17", dataset="gtex_v8")
    assert out.result == []
    assert any("does not belong to that dataset's release" in w for w in out.warnings)
    assert any("not an absence of expression" in w for w in out.warnings)


# ---------------------------------------------------------------- top_tissues contract


def test_top_tissues_matches_the_contract_and_always_carries_the_rna_test_tissues(monkeypatch):
    """CP1-17: top_tissues(symbol, n) -> {gene, gencode_id, tissues:[{tissue, uberon, median_tpm}]}."""
    fake_http(monkeypatch, SETUP)
    out = gtex.top_tissues("CFTR", n=5)
    r = out.result
    assert r["gene"] == "CFTR" and r["gencode_id"] == "ENSG00000001626.14"
    assert len(r["tissues"]) == 5 and r["tissues_total"] == 54
    for row in r["tissues"]:
        assert {"tissue", "uberon", "median_tpm"} <= set(row)
    # ordered by median TPM, highest first
    values = [row["median_tpm"] for row in r["tissues"]]
    assert values == sorted(values, reverse=True)
    assert r["tissues"][0]["tissue"] == "Pancreas"
    # blood, lymphoblastoid, fibroblasts, skin and muscle are there whether or not they rank
    names = [row["tissue"] for row in r["rna_test_tissues"]]
    assert names == [n for n, _ in gtex.RNA_TEST_TISSUES]
    for row in r["rna_test_tissues"]:
        assert row["why_it_matters"] or row.get("note")


def test_top_tissues_hint_names_an_accessible_tissue_or_says_there_is_none(monkeypatch):
    fake_http(monkeypatch, SETUP)
    out = gtex.top_tissues("CFTR", n=3)
    hint = out.result["rna_test_hint"]
    # CFTR is absent from blood, fibroblasts and muscle; skin is the only accessible tissue above 1 TPM
    assert "Skin" in hint
    blood = next(t for t in out.result["rna_test_tissues"] if t["tissue"] == "Whole_Blood")
    assert blood["band"] == "not_detected"


def test_top_tissues_says_nothing_was_measured_for_an_unknown_gene(monkeypatch):
    fake_http(monkeypatch, {"geneId=NOTAGENE123": "reference_gene_missing.json"})
    out = gtex.top_tissues("NOTAGENE123")
    r = out.result
    assert r["tissues"] == [] and r["rna_test_tissues"] == []
    assert "nothing was measured, nothing is implied" in r["note"]


def test_caveats_name_the_nmd_trap_and_the_bulk_tissue_limit(monkeypatch):
    fake_http(monkeypatch, SETUP)
    caveats = " ".join(gtex.top_tissues("CFTR", n=2).result["caveats"])
    assert "nonsense-mediated decay" in caveats
    assert "bulk tissue" in caveats
    assert "not a detection limit" in caveats


def test_a_failed_expression_call_leaves_the_gene_row_and_a_warning(monkeypatch):
    def fake(url, source, **kw):
        if "reference/gene" in url:
            return Response(url, 200, text("reference_gene_cftr_v26.json"), "2026-10-06T00:00:00+00:00", True)
        raise SourceError("GTEx medianGeneExpression", url, 503, "service unavailable")

    monkeypatch.setattr(gtex, "get_json", fake)
    out = gtex.top_tissues("CFTR")
    r = out.result
    # the gene row survives, and nothing negative is claimed about expression
    assert r["gencode_id"] == "ENSG00000001626.14"
    assert r["expression_retrieved"] is False
    assert r["tissues"] == [] and r["rna_test_tissues"] == []
    assert "Nothing is implied" in r["note"]
    assert "cannot be said from GTEx" in r["rna_test_hint"]
    assert any("GTEx medianGeneExpression unavailable" in w for w in out.warnings)


def test_a_lookup_mismatch_never_becomes_a_negative_finding(monkeypatch):
    """P0: an empty table means the GENCODE id did not match the release, not "not expressed"."""
    empty = json.dumps({"data": [], "paging_info": {"numberOfPages": 0, "totalNumberOfItems": 0}})
    fake_http(monkeypatch, {"geneId=CFTR": "reference_gene_cftr_v26.json", "medianGeneExpression": empty})
    out = gtex.top_tissues("CFTR", n=3)
    r = out.result
    assert r["tissues"] == [] and r["tissues_total"] == 0
    assert r["expression_retrieved"] is False
    hint = r["rna_test_hint"]
    assert "cannot be said from GTEx" in hint
    assert "reaches 1 TPM" not in hint
    assert "would need the disease tissue" not in hint
    for row in r["rna_test_tissues"]:
        assert "not in gtex_v8" not in (row.get("note") or "")
        assert "not read at all" in (row.get("note") or "")
    from zebra.commands import expression as command
    rendered = command.render(r)
    assert "reaches 1 TPM" not in rendered


def test_a_tissue_with_no_usable_median_is_not_counted_as_zero(monkeypatch):
    rows = {"data": [{"median": None, "tissueSiteDetailId": "Whole_Blood", "ontologyId": "UBERON:0013756",
                      "gencodeId": "ENSG00000001626.14", "unit": "TPM"},
                     {"median": "NaN", "tissueSiteDetailId": "Muscle_Skeletal", "ontologyId": "UBERON:0011907",
                      "gencodeId": "ENSG00000001626.14", "unit": "TPM"},
                     {"median": 0.5, "tissueSiteDetailId": "Lung", "ontologyId": "UBERON:0008952",
                      "gencodeId": "ENSG00000001626.14", "unit": "TPM"}],
            "paging_info": {"numberOfPages": 1, "totalNumberOfItems": 3}}
    fake_http(monkeypatch, {"geneId=CFTR": "reference_gene_cftr_v26.json",
                            "medianGeneExpression": json.dumps(rows)})
    out = gtex.top_tissues("CFTR", n=3)
    blood = next(t for t in out.result["rna_test_tissues"] if t["tissue"] == "Whole_Blood")
    assert blood["median_tpm"] is None and blood["band"] is None
    # a NaN must not sort to the top as if it were the highest-expressing tissue
    assert out.result["tissues"][0]["tissue"] == "Lung"
    assert any("no usable median" in w for w in out.warnings)
    assert "had no usable median" in out.result["rna_test_hint"]


def test_rows_for_another_gene_are_dropped_with_a_warning(monkeypatch):
    rows = {"data": [{"median": 9.0, "tissueSiteDetailId": "Lung", "ontologyId": "UBERON:0008952",
                      "gencodeId": "ENSG00000999999.1", "unit": "TPM"},
                     {"median": 1.0, "tissueSiteDetailId": "Pancreas", "ontologyId": "UBERON:0001150",
                      "gencodeId": "ENSG00000001626.14", "unit": "TPM"}],
            "paging_info": {"numberOfPages": 1, "totalNumberOfItems": 2}}
    fake_http(monkeypatch, {"geneId=CFTR": "reference_gene_cftr_v26.json",
                            "medianGeneExpression": json.dumps(rows)})
    out = gtex.top_tissues("CFTR", n=3)
    assert [t["tissue"] for t in out.result["tissues"]] == ["Pancreas"]
    assert any("ENSG00000999999.1" in w for w in out.warnings)


def test_a_data_field_that_is_not_a_list_is_an_error_not_a_release_mismatch(monkeypatch):
    for payload in ('{"data": "oops"}', '{"data": {"a": 1}}'):
        fake_http(monkeypatch, {"medianGeneExpression": payload})
        with pytest.raises(SourceError) as err:
            gtex.median_expression("ENSG00000001626.14")
        assert "not a list" in str(err.value)


def test_top_is_validated_rather_than_silently_clamped():
    for bad in (0, -1, 61, "x", None):
        with pytest.raises(UsageError):
            gtex.top_tissues("CFTR", n=bad)


# ---------------------------------------------------------------- the command


def test_the_expression_command_is_registered():
    """CP1-17: `zebra expression <gene>` exists and takes the flags the skill needs."""
    import argparse

    from zebra.commands import expression as command

    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command")
    command.register(sub)
    args = parser.parse_args(["expression", "CFTR", "--top", "4", "--dataset", "gtex_v10"])
    assert args.gene == "CFTR" and args.top == 4 and args.dataset == "gtex_v10"
    assert callable(args.func)
    assert parser.parse_args(["expression", "CFTR"]).dataset == gtex.DEFAULT_DATASET


def test_the_expression_command_refuses_an_empty_gene():
    import argparse

    from zebra.commands import expression as command

    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command")
    command.register(sub)
    args = parser.parse_args(["expression", "  "])
    with pytest.raises(UsageError):
        args.func(args)


def test_render_shows_the_rna_test_tissues_and_the_caveats(monkeypatch):
    fake_http(monkeypatch, SETUP)
    from zebra.commands import expression as command

    out = gtex.top_tissues("CFTR", n=3)
    rendered = command.render(out.result)
    assert "Whole_Blood" in rendered and "Cells_Cultured_fibroblasts" in rendered
    assert "Muscle_Skeletal" in rendered
    assert "which RNA test:" in rendered
    assert "nonsense-mediated decay" in rendered


# ---------------------------------------------------------------- live


@pytest.mark.live
def test_live_gtex_cftr():
    out = gtex.top_tissues("CFTR", n=5)
    r = out.result
    assert r["gencode_id"] and r["tissues_total"] >= 50
    assert r["expression_retrieved"] is True
    # CFTR is a secretory-epithelium gene: pancreas and gut, not blood or fibroblasts
    assert r["tissues"][0]["tissue"] == "Pancreas" and r["tissues"][0]["median_tpm"] > 10
    assert all(t["uberon"] for t in r["tissues"])
    blood = next(t for t in r["rna_test_tissues"] if t["tissue"] == "Whole_Blood")
    fibro = next(t for t in r["rna_test_tissues"] if t["tissue"] == "Cells_Cultured_fibroblasts")
    assert blood["band"] == "not_detected" and fibro["band"] == "not_detected"


@pytest.mark.live
def test_live_gtex_both_releases_resolve_their_own_gencode_id():
    ids = {d: gtex.resolve_gene("CFTR", dataset=d).result["gencode_id"] for d in gtex.DATASETS}
    assert len(set(ids.values())) == len(ids), ids
    for dataset, gencode_id in ids.items():
        rows = gtex.median_expression(gencode_id, dataset=dataset).result
        assert rows, (dataset, gencode_id)


@pytest.mark.live
def test_live_gtex_hint_points_at_blood_for_a_blood_expressed_gene():
    # HBB is a blood gene: the hint must name blood, not send anyone for a biopsy
    out = gtex.top_tissues("HBB", n=3)
    blood = next(t for t in out.result["rna_test_tissues"] if t["tissue"] == "Whole_Blood")
    assert blood["median_tpm"] > 100 and blood["band"] == "expressed"
    assert "Whole_Blood" in out.result["rna_test_hint"]


@pytest.mark.live
def test_live_gtex_hint_refuses_blood_for_a_brain_restricted_gene():
    # SCN1A is brain-restricted: no accessible tissue reaches 1 TPM, and the hint says so
    out = gtex.top_tissues("SCN1A", n=3)
    assert out.result["tissues"][0]["tissue"].startswith("Brain")
    blood = next(t for t in out.result["rna_test_tissues"] if t["tissue"] == "Whole_Blood")
    assert (blood["median_tpm"] or 0.0) < 1.0
    assert "none of the accessible tissues" in out.result["rna_test_hint"]
