"""zebra cnv: CNV/CMA intervals, exon-level del/dup, copy-number and repeat results.

Offline tests answer Ensembl and ClinGen from captured responses (tests/fixtures/cnv,
see make_fixtures.py and provenance.json there; the 15q fixture of the original tests
is read from tests/fixtures/vcf). Moved here from tests/test_vcf.py in v0.2; tests
named test_<id>_… are the regression tests for that review finding.
"""

from __future__ import annotations

import gzip
import json
import urllib.error
from pathlib import Path

import pytest

from zebra.core import Outcome, UsageError
from zebra.http import Response, SourceError, source_record

FIX = Path(__file__).parent / "fixtures" / "vcf"
CNV_FIX = Path(__file__).parent / "fixtures" / "cnv"
STAMP = "2026-10-06T00:00:00+00:00"


def _clingen_text(name: str) -> str:
    return gzip.decompress((CNV_FIX / (name + ".gz")).read_bytes()).decode("utf-8")


@pytest.fixture
def fake_cnv_sources(monkeypatch):
    """Ensembl overlap/lookup and the ClinGen bulk files answered from captured responses."""
    from zebra import cnv as C
    from zebra.sources import clingen
    from zebra.sources import ensembl as ens

    overlap = json.loads(gzip.decompress((CNV_FIX / "overlap.json.gz").read_bytes()))
    old_overlap = json.loads((FIX / "overlap_cnv.json").read_text())["response"]  # 15:25200000-25500000
    lookups = json.loads((CNV_FIX / "lookups.json").read_text())
    calls = {"overlap": [], "dosage": [], "lookup": [], "validate": []}

    def get_json(url, source, params=None, **kw):
        calls["overlap"].append(url)
        build = "GRCh37" if "grch37." in url else "GRCh38"
        region = url.rsplit("/", 1)[1]
        key = f"{build}:{region}"
        if key in overlap:
            body = overlap[key]
        elif region == "15:25200000-25500000":
            body = old_overlap
        else:
            body = []  # a window no capture covers: no genes (tests that use it say so)
        return Response(url, 200, json.dumps(body), STAMP, False)

    def request(url, source, validate=None, **kw):
        """The ClinGen bulk TSVs, served like zebra.http.request: a body `validate` rejects raises."""
        calls["dosage"].append(url)
        text = _clingen_text(url.rsplit("/", 1)[1])
        calls["validate"].append(validate is not None)
        if validate is not None and validate(text):
            raise SourceError(source, url, 200, validate(text))
        return Response(url, 200, text, STAMP, False)

    def lookup_symbol(symbol, assembly="GRCh38", expand=False):
        calls["lookup"].append((symbol, expand))
        data = lookups.get(f"{assembly}:{symbol.upper()}")
        if data is None:
            if symbol.upper() == "NOSUCHGENE1":  # what Ensembl answers for an unknown symbol
                raise SourceError("Ensembl lookup", f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}",
                                  400, "No valid lookup found for symbol NOSUCHGENE1")
            raise SourceError("Ensembl lookup", f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}", 503,
                              "not captured in the offline fixture")
        if not expand:
            data = {k: v for k, v in data.items() if k != "Transcript"}
        return Outcome(data, sources=[source_record("Ensembl lookup", symbol,
                                                    url="https://rest.ensembl.org/lookup/symbol")])

    monkeypatch.setattr(C, "get_json", get_json)
    monkeypatch.setattr(clingen, "request", request)
    monkeypatch.setattr(ens, "lookup_symbol", lookup_symbol)
    assert clingen.DOSAGE_TSV["GRCh38"] and clingen.REGION_TSV["GRCh38"]  # the real urls the bulk fetch uses
    return calls


# ------------------------------------------------------- moved from test_vcf.py

def test_p1d_input_forms_are_read_or_refused():
    from zebra import cnv as C

    iscn = C.parse("arr[GRCh38] 15q11.2q13.1(23123715_28193120)x1")
    assert iscn["kind"] == "cnv" and iscn["chrom"] == "15" and iscn["copy_number"] == 1
    assert (iscn["start"], iscn["end"]) == (23123715, 28193120) and iscn["cnv_type"] == "loss"
    assert iscn["assembly"] == "GRCh38"
    assert C.parse("seq[hg19] 22q11.21(18,648,855_21,800,471)x3")["assembly"] == "GRCh37"
    assert C.parse("seq[hg19] 22q11.21(18,648,855_21,800,471)x3")["cnv_type"] == "gain"

    coords = C.parse("chr15:23123715-28193120 loss")
    assert coords["kind"] == "cnv" and coords["cnv_type"] == "loss" and coords["copy_number"] is None

    exon = C.parse("DMD exon 45-50 deletion")
    assert {k: exon[k] for k in ("kind", "gene", "transcript", "first", "last", "cds_start", "cds_end", "cnv_type")} \
        == {"kind": "exon_cnv", "gene": "DMD", "transcript": None, "first": 45, "last": 50, "cds_start": None,
            "cds_end": None, "cnv_type": "loss"}
    assert C.parse("MLPA: DMD exon 8 duplication".split(": ")[1])["cnv_type"] == "gain"

    hgvs = C.parse("NM_004006.3:c.6439-?_7309+?del")
    assert hgvs["kind"] == "exon_cnv" and (hgvs["cds_start"], hgvs["cds_end"]) == (6439, 7309)
    assert hgvs["transcript"] == "NM_004006.3" and hgvs["cnv_type"] == "loss"

    cn = C.parse("SMN1 exon 7 copy number 0")
    assert {k: cn[k] for k in ("kind", "gene", "exon", "copy_number")} == \
        {"kind": "copy_number", "gene": "SMN1", "exon": "7", "copy_number": 0}
    assert C.parse("SMN2 copy number 2")["copy_number"] == 2
    assert C.parse("SMN1 0 copies")["copy_number"] == 0

    rep = C.parse("FMR1 CGG 230")
    assert {k: rep[k] for k in ("kind", "gene", "motif", "repeat_count")} == \
        {"kind": "repeat_expansion", "gene": "FMR1", "motif": "CGG", "repeat_count": 230}
    assert C.parse("HTT (CAG)n 42 repeats")["repeat_count"] == 42
    assert C.parse("FMR1 CGG 55-200")["repeat_count"] == "55-200"

    # nothing is invented: a band without coordinates is refused
    with pytest.raises(UsageError, match="names a band but no coordinates"):
        C.parse("del(15)(q11.2q13.1)")
    with pytest.raises(UsageError, match="cannot read"):
        C.parse("the array was abnormal")
    with pytest.raises(UsageError, match="Accepted forms"):
        C.parse("")


def test_p1d_cnv_reports_genes_dosage_and_acmg_inputs(fake_cnv_sources):
    from zebra import cnv as C

    got = C.card("chr15:25200000-25500000 loss", copies=1, inheritance="de_novo")
    r = got.result
    assert r["kind"] == "cnv" and r["region"] == {"chrom": "15", "start": 25200000, "end": 25500000,
                                                  "length_bp": 300001, "band_as_reported": None}
    assert r["cnv_type"] == "loss" and r["copy_number"] == 1
    assert r["genes"]["total"] == 38 and r["genes"]["protein_coding"] == 1
    assert r["genes"]["protein_coding_symbols"] == ["UBE3A"]
    # one bulk TSV per table (genes, regions), never one request per gene
    assert sorted(u.rsplit("/", 1)[1] for u in fake_cnv_sources["dosage"]) == \
        ["ClinGen_gene_curation_list_GRCh38.tsv", "ClinGen_region_curation_list_GRCh38.tsv"]
    inputs = r["acmg_cnv_inputs"]
    assert inputs["framework"].startswith("ACMG/ClinGen technical standard")
    assert inputs["section_1_variant_type"]["contains_protein_coding_genes"] is True
    hi = inputs["section_2_overlap_with_established_regions_or_genes"]["clingen_established_haploinsufficient_genes"]
    assert hi == [{"gene": "UBE3A", "score": "3", "overlap": "whole gene",
                   "description": "Sufficient evidence for dosage pathogenicity"}]
    assert inputs["section_3_gene_number"]["protein_coding_genes"] == 1
    assert inputs["section_4_case_level_evidence"] is None
    assert inputs["section_5_inheritance_and_family_history"] == {"reported_inheritance": "de_novo"}
    # the interpretation itself is never invented
    assert inputs["classification"] is None and "does not score or classify" in inputs["classification_note"]
    assert any(s["db"] == "Ensembl overlap" for s in got.sources)
    assert any(s["db"] == "ClinGen dosage sensitivity" for s in got.sources)
    assert any(s["db"] == "ClinGen dosage sensitivity (regions)" for s in got.sources)
    assert r["caveats"] and any("not a breakpoint" in c for c in r["caveats"])


def test_p1d_cnv_without_a_copy_number_says_so(fake_cnv_sources):
    from zebra import cnv as C

    got = C.card("15:25200000-25500000")
    assert got.result["copy_number"] is None and got.result["cnv_type"] == "unknown"
    assert any("no copy number in the input" in w for w in got.warnings)
    with pytest.raises(UsageError, match="copy number 1 and --copies"):
        C.card("arr[GRCh38] 15q11.2(25200000_25500000)x1", copies=3)


def test_p1d_exon_deletion_resolves_coordinates_and_frame(fake_cnv_sources):
    """DMD exon 45-50 is 871 bp: out of frame, and exon 51's removal restores it."""
    from zebra import cnv as C

    got = C.card("DMD exon 45-50 deletion")
    r = got.result
    assert r["kind"] == "exon_cnv" and r["gene"] == "DMD" and r["cnv_type"] == "loss"
    assert r["transcript"]["resolved"] == "ENST00000357033" and r["transcript"]["exon_total"] == 79
    assert r["exons"]["count"] == 6
    assert [e["length_bp"] for e in r["exons"]["per_exon"]] == [176, 148, 150, 186, 102, 109]
    assert r["coordinates"] == {"chrom": "X", "start": 31819975, "end": 31968514,
                                "note": "the exon boundaries; the real breakpoints lie in the flanking introns"}
    assert r["frame"]["bases"] == 871 and r["frame"]["modulo_3"] == 1
    assert r["frame"]["consequence"] == "out of frame"
    restoration = {c["exon"]: c["restores_frame"] for c in r["frame_restoration"]["candidates"]}
    assert restoration == {44: False, 51: True}
    assert "separate question" in r["frame_restoration"]["note"]
    assert fake_cnv_sources["lookup"] == [("DMD", True)]


def test_p1d_exon_hgvs_form_gives_the_same_frame(fake_cnv_sources):
    from zebra import cnv as C

    got = C.card("NM_004006.3:c.6439-?_7309+?del", gene="DMD")
    r = got.result
    assert r["cds_span"] == {"from": 6439, "to": 7309, "length_bp": 871}
    assert r["frame"]["modulo_3"] == 1 and r["frame"]["consequence"] == "out of frame"
    assert r["transcript"]["as_reported"] == "NM_004006.3"
    assert any("not sequenced" in w for w in got.warnings)
    # an in-frame example: a 3-base multiple
    inframe = C.card("NM_004006.3:c.6439-?_7308+?del", gene="DMD").result
    assert inframe["frame"]["bases"] == 870 and inframe["frame"]["consequence"] == "in frame"


def test_p1d_exon_number_beyond_the_transcript_is_refused(fake_cnv_sources):
    from zebra import cnv as C

    with pytest.raises(UsageError, match="has 79 exons"):
        C.card("DMD exon 80 deletion")


def test_p1d_copy_number_and_repeats_are_structure_not_prediction(fake_cnv_sources):
    """No classification, ever; since CP1-5 the cited meaning is given beside the structure."""
    from zebra import cnv as C

    got = C.card("SMN1 exon 7 copy number 0", related=["SMN2 copy number 2"], method="MLPA")
    r = got.result
    assert r["kind"] == "copy_number" and r["gene"] == "SMN1" and r["copy_number"] == 0 and r["exon"] == "7"
    assert r["related_results"] == [{"gene": "SMN2", "copy_number": 2, "exon": None}]
    assert r["method"] == "MLPA" and r["classification"] is None
    assert any("2+0" in m for m in r["mechanism"])  # the silent-carrier caveat
    assert any("modifier" in m for m in r["mechanism"])
    assert r["gene_location"]["seq_region_name"] == "5"

    rep = C.card("FMR1 CGG 230").result
    assert rep["kind"] == "repeat_expansion" and rep["motif"] == "CGG" and rep["repeat_count"] == 230
    assert rep["classification"] is None
    assert rep["thresholds"]["categories_for_this_result"][0]["categories"] == ["full mutation"]
    assert "research reference" in rep["thresholds"]["label"]
    assert any("length measurement" in m for m in rep["mechanism"])
    with pytest.raises(UsageError, match="--related takes further copy-number results"):
        C.card("SMN1 copy number 1", related=["FMR1 CGG 230"])


def test_p1d_findings_are_recorded_in_the_case(fake_cnv_sources, tmp_path, capsys):
    """Every new form reaches case.json as its own variant kind, and case.json stays readable."""
    from zebra import case as case_mod
    from zebra.cli import main

    case_dir = tmp_path / "cnvcase"
    case_mod.init(str(case_dir))
    for text, extra in (("DMD exon 45-50 deletion", []),
                        ("SMN1 exon 7 copy number 0", ["--related", "SMN2 copy number 2"]),
                        ("FMR1 CGG 230", []),
                        ("15:25200000-25500000 loss", ["--copies", "1"])):
        code = main(["--json", "cnv", text, "--case", str(case_dir), "--record", *extra])
        env = json.loads(capsys.readouterr().out)
        assert code == 0 and env["ok"], env
        assert env["result"]["recorded_in_case"]["kind"] in case_mod.VARIANT_KINDS
    data = case_mod.load(str(case_dir))
    assert [v["kind"] for v in data["variants"]] == ["exon_cnv", "copy_number", "repeat_expansion", "cnv"]
    summary = case_mod.summary(str(case_dir))
    labels = {v["kind"]: v["label"] for v in summary["variants"]}
    assert labels["exon_cnv"] == "DMD exon 45-50 loss"
    assert labels["copy_number"] == "SMN1 copy number 0"
    assert labels["repeat_expansion"] == "FMR1 CGG 230 repeats"
    assert labels["cnv"] == "loss 15:25200000-25500000 CN1"
    # a case written before these kinds existed still reads
    data["variants"].append({"id": "v9", "gene": "SCN1A", "hgvs_c": "NM_001165963.4:c.2134C>T"})
    (case_dir / "case.json").write_text(json.dumps(data), "utf-8")
    old = next(v for v in case_mod.summary(str(case_dir))["variants"] if v["id"] == "v9")
    assert old["kind"] == "small" and old["label"] == "NM_001165963.4:c.2134C>T"


def test_p1d_record_without_a_case_is_refused(fake_cnv_sources, capsys):
    from zebra.cli import main

    assert main(["--json", "cnv", "FMR1 CGG 230", "--record"]) == 2
    env = json.loads(capsys.readouterr().out)
    assert env["error"]["type"] == "UsageError" and "--record needs a case" in env["error"]["message"]


def test_p1d_a_small_indel_is_sent_to_zebra_variant(fake_cnv_sources):
    """The exon-boundary HGVS form must not swallow an ordinary 3 bp deletion."""
    from zebra import cnv as C

    with pytest.raises(UsageError, match="that is a small variant"):
        C.card("NM_000492.4:c.1521_1523del")
    assert C.parse("NM_004006.3:c.6439-?_7309+?del")["cds_start"] == 6439  # '?' marks the exon form
    assert C.parse("NM_004006.3:c.6439_7309del")["cds_start"] == 6439  # 871 bp: large enough to be an exon event


def test_p1d_impossible_coordinates_are_refused(fake_cnv_sources):
    """A CNV interval is checked against the real chromosome lengths before anything is looked up."""
    from zebra import cnv as C

    with pytest.raises(UsageError, match="runs past the end of chromosome 15"):
        C.card("chr15:23123715-999000000 loss")
    with pytest.raises(UsageError, match="not a chromosome zebra can place"):
        C.card("chr23:100-200 loss")
    with pytest.raises(UsageError, match="ends before it starts"):
        C.card("15:25500000-25200000 loss")
    assert fake_cnv_sources["overlap"] == []  # refused before any request
    # an ISCN string's own build is honoured, and checked against that build's lengths
    assert C.parse("arr[GRCh37] 15q11.2(23123715_28193120)x1")["assembly"] == "GRCh37"


def test_p1d_a_copy_number_on_x_is_not_a_gain_by_itself(fake_cnv_sources):
    """One copy of X is normal in a male and a loss in a female: without the sex zebra must not pick one."""
    from zebra import cnv as C

    for text, copies in (("arr[GRCh38] Xq28(154021812_154137257)x2", 2),
                         ("arr[GRCh38] Xp21.1(31000000_31200000)x1", 1),
                         ("chrX:31000000-31200000 x1", 1)):
        got = C.parse(text)
        assert got["chrom"] == "X" and got["copy_number"] == copies
        assert got["cnv_type"] == "unknown"
        assert "normal in a male" in got["cnv_type_note"]
    # zero copies is a loss in either sex, and the report's own word always wins
    assert C.parse("arr[GRCh38] Xq28(154021812_154137257)x0")["cnv_type"] == "loss"
    assert C.parse("chrX:31000000-31200000 deletion")["cnv_type"] == "loss"
    # autosomes are unambiguous
    assert C.parse("chr15:23123715-25193120 x1")["cnv_type"] == "loss"
    assert C.parse("chr15:23123715-25193120 x3")["cnv_type"] == "gain"
    got = C.card("chrX:31000000-31200000 x1")
    assert got.result["cnv_type"] == "unknown"
    assert any("normal in a male" in w for w in got.warnings)
    assert got.result["acmg_cnv_inputs"]["scope"].startswith("unknown")


def test_p1d_a_transcript_in_the_gene_slot_is_handled(fake_cnv_sources):
    from zebra import cnv as C

    got = C.parse("NM_004006.3 exon 45-50 deletion")
    assert got["gene"] is None and got["transcript"] == "NM_004006.3"
    with pytest.raises(UsageError, match="which gene"):
        C.card("NM_004006.3 exon 45-50 deletion")
    with_gene = C.card("NM_004006.3 exon 45-50 deletion", gene="DMD").result
    assert with_gene["frame"]["bases"] == 871 and with_gene["gene"] == "DMD"
    assert with_gene["transcript"]["as_reported"] == "NM_004006.3"
    with pytest.raises(UsageError, match="exons are numbered from 1"):
        C.parse("DMD exon 0 deletion")


def test_p1d_a_band_pair_is_never_read_as_coordinates(fake_cnv_sources):
    """"chr22:11.21-11.23" is 22q11.21-q11.23, not an 11 bp interval."""
    from zebra import cnv as C

    with pytest.raises(UsageError, match="names a band but no coordinates|does not read as one interval"):
        C.card("chr22:11.21-11.23 deletion")
    with pytest.raises(UsageError, match="names a band but no coordinates|does not read as one interval"):
        C.card("chr7:11.23-11.25 del")
    assert fake_cnv_sources["overlap"] == []  # nothing was looked up
    # a real interval with a dotted separator still reads
    assert C.parse("15:25200000..25500000 loss")["end"] == 25500000


def test_p1d_an_oversized_or_unreadable_input_is_refused(fake_cnv_sources):
    from zebra import cnv as C

    with pytest.raises(UsageError, match="runs past the end of chromosome"):
        C.card("chr1:1-999999999 loss")
    assert C.MAX_REGION_BP == 50_000_000
    with pytest.raises(UsageError, match="will not walk more than"):
        C.genes_in_region("1", 1, 60_000_000)
    with pytest.raises(UsageError, match="'~' \\(x1~2\\)|copy number followed by"):
        C.card("arr[GRCh38] 15q11.2(23123715_23200000)x1-2")
    with pytest.raises(UsageError, match="has text before the ISCN result"):
        C.parse("junk arr[GRCh38] 15q11.2(1000_2000)x1 junk")
    with pytest.raises(UsageError, match="at most 400 characters"):
        C.parse("chr15:1-2 loss " + "x" * 500)


def test_p1d_a_copy_number_that_contradicts_the_report_is_refused(fake_cnv_sources):
    from zebra import cnv as C

    with pytest.raises(UsageError, match="the report says loss and copy number 3 on 15 cannot be a loss"):
        C.card("chr15:23123715-23200000 loss", copies=3)
    with pytest.raises(UsageError, match="--copies applies to a CNV or a copy-number result"):
        C.card("FMR1 CGG 230", copies=3)
    with pytest.raises(UsageError, match="--copies applies to a CNV or a copy-number result"):
        C.card("DMD exon 45-50 deletion", copies=1)
    # agreeing values are fine
    assert C.card("chr15:25200000-25500000 loss", copies=1).result["copy_number"] == 1


def test_p1d_the_record_keeps_the_exon_and_every_related_result(fake_cnv_sources, tmp_path):
    """SMN2 copy number is the modifier the mechanism note names: it must persist."""
    from zebra import case as case_mod
    from zebra import cnv as C

    case_dir = tmp_path / "smn"
    case_mod.init(str(case_dir))
    got = C.card("SMN1 exon 7 copy number 0", related=["SMN2 copy number 2"], method="MLPA")
    fields = C.case_fields(got.result, method="MLPA")
    entry = case_mod.add_variant(str(case_dir), **fields)
    assert entry["copy_number"] == 0 and entry["exons"] == "7" and entry["method"] == "MLPA"
    assert "SMN2 copy number 2" in entry["note"]
    assert entry["genes"] == ["SMN1", "SMN2"]


# ---------------------------------------------------------- CP0-3: regions

ACCEPTANCE = [
    # (id, ISCN as a laboratory would print it, ISCA id, region name fragment, HI score, relation)
    ("22q11.2 A-D (DiGeorge)", "arr[GRCh38] 22q11.21(18648855_21800471)x1", "ISCA-37446",
     "22q11.2 recurrent (DGS) region (proximal, A-B, A-C, or A-D)", "3", "cnv_contains_it"),
    ("1p36 terminal", "arr[GRCh38] 1p36.33p36.31(849466_6823542)x1", "ISCA-37434",
     "1p36 terminal region", "3", "cnv_contains_it"),
    ("Williams-Beuren 7q11.23", "arr[GRCh38] 7q11.23(73330452_74778226)x1", "ISCA-37392",
     "7q11.23 recurrent (Williams-Beuren syndrome) region", "3", "cnv_contains_it"),
    ("PWS/AS 15q11-q13", "arr[GRCh38] 15q11.2q13.1(23123715_28193120)x1", "ISCA-37478",
     "15q11.2q13 recurrent (PWS/AS) region", "3", "cnv_contains_it"),
    ("Smith-Magenis 17p11.2", "arr[GRCh38] 17p11.2(16810000_20450000)x1", "ISCA-37418",
     "17p11.2 recurrent (SMS/PLS) region", "3", "cnv_contains_it"),
    # a typical array call stops a few kb short of the curated BP4-BP5 boundary: partial, said so
    ("16p11.2 BP4-BP5", "arr[GRCh38] 16p11.2(29580020_30180020)x1", "ISCA-37400",
     "16p11.2 recurrent region (proximal, BP4-BP5)", "3", "partial"),
]


@pytest.mark.parametrize("label,iscn,isca,name,hi,relation", ACCEPTANCE, ids=[a[0] for a in ACCEPTANCE])
def test_cp0_3_each_recurrent_deletion_reports_its_curated_region(fake_cnv_sources, label, iscn, isca, name, hi,
                                                                   relation):
    from zebra import cnv as C

    got = C.card(iscn)
    r = got.result
    regions = r["clingen_regions"]
    assert regions["status"] == "checked" and regions["regions_in_file"] == 513
    entry = next(e for e in regions["dosage_curated"] if e["isca_id"] == isca)
    assert name in entry["name"]
    assert entry["haploinsufficiency"]["score"] == hi
    assert entry["haploinsufficiency"]["description"] == "Sufficient evidence for dosage pathogenicity"
    assert entry["relation"] == relation
    assert entry["url"] == f"https://search.clinicalgenome.org/kb/gene-dosage/region/{isca}"
    assert entry["location"].startswith("chr")
    sec2 = r["acmg_cnv_inputs"]["section_2_overlap_with_established_regions_or_genes"]
    est = next(e for e in sec2["clingen_established_regions"] if e["isca_id"] == isca)
    if relation == "cnv_contains_it":
        assert entry["region_fraction_covered"] == 1.0
        assert est["acmg_section_2"].startswith("2A row")
    else:
        assert 0.9 < entry["region_fraction_covered"] < 1.0
        # partial, but the named gene (TBX6) is inside: the curator decides 2A vs 2B, which the hint says
        assert est["acmg_section_2"].startswith("2A or 2B row") and "TBX6" in est["acmg_section_2"]
        assert entry["region_bp_outside_cnv"]["after"] > 0
    assert isca in got.text and "HI 3" in got.text
    assert not any("could not be read" in w for w in got.warnings)


def test_cp0_3_the_22q11_deletion_names_tbx1_and_no_longer_reads_none(fake_cnv_sources):
    from zebra import cnv as C

    r = C.card("arr[GRCh38] 22q11.21(18648855_21800471)x1").result
    dgs = next(e for e in r["clingen_regions"]["dosage_curated"] if e["isca_id"] == "ISCA-37446")
    assert dgs["named_genes"] == [{"gene": "TBX1", "in_cnv": "whole gene"}]
    assert dgs["triplosensitivity"]["score"] == "3" and dgs["haploinsufficiency"]["pmids"]
    # the distal D-E/D-F region is only partly covered: reported with its fraction, as 2B
    distal = next(e for e in r["clingen_regions"]["dosage_curated"] if e["isca_id"] == "ISCA-37397")
    assert distal["relation"] == "partial" and 0.2 < distal["region_fraction_covered"] < 0.3
    # population regions are kept apart from the curated ones
    assert r["clingen_regions"]["benign_or_population"]
    assert all(e["haploinsufficiency_score"] == "40" for e in r["clingen_regions"]["benign_or_population"])
    assert r["clingen_regions"]["benign_or_population_total"] == len(r["clingen_regions"]["benign_or_population"])


def test_cp0_3_every_spanned_gene_is_checked_without_a_40_gene_cut(fake_cnv_sources):
    from zebra import cnv as C

    got = C.card("arr[GRCh38] 22q11.21(18648855_21800471)x1")
    dosage = got.result["clingen_dosage"]
    assert dosage["status"] == "checked"
    assert dosage["genes_checked"] == len(got.result["genes"]["protein_coding_symbols"]) == 58  # not 40
    assert not any("first 40" in w for w in got.warnings)
    assert {"TBX1", "CRKL", "PRODH"} <= {g["gene"] for g in dosage["records"]}


def test_cp0_3_a_gene_desert_reports_no_curated_region(fake_cnv_sources):
    """13q21 (chr13:55-56 Mb): no ClinGen region, no ClinGen gene, no protein-coding gene (captured live)."""
    from zebra import cnv as C

    got = C.card("chr13:55000000-56000000 loss", copies=1)
    r = got.result
    assert r["clingen_regions"]["status"] == "checked" and r["clingen_regions"]["regions_in_file"] == 513
    assert r["clingen_regions"]["dosage_curated"] == [] and r["clingen_regions"]["benign_or_population"] == []
    assert r["clingen_dosage"]["status"] == "checked" and r["clingen_dosage"]["records"] == []
    assert r["genes"]["protein_coding"] == 0 and r["genes"]["total"] == 8
    assert "no curated (dosage-scored) region overlaps this interval" in got.text
    assert "NOT CHECKED" not in got.text


def test_cp0_3_grch37_coordinates_use_the_grch37_region_file(fake_cnv_sources):
    from zebra import cnv as C

    got = C.card("chr22:18648855-21800471 hg19 deletion", copies=1)
    r = got.result
    assert r["assembly"] == "GRCh37"
    dgs = next(e for e in r["clingen_regions"]["dosage_curated"] if e["isca_id"] == "ISCA-37446")
    assert dgs["location"] == "chr22:18912231-21465672" and dgs["relation"] == "cnv_contains_it"
    assert any(u.endswith("ClinGen_region_curation_list_GRCh37.tsv") for u in fake_cnv_sources["dosage"])


def test_cp0_3_region_table_parsing_and_overlap_arithmetic():
    from zebra.sources import clingen

    regions = clingen.parse_region_table(_clingen_text("ClinGen_region_curation_list_GRCh38.tsv"))
    assert len(regions) == 513
    assert sum(1 for r in regions if r["start"] is None) == 8  # 'tbd' locations
    wbs = next(r for r in regions if r["isca_id"] == "ISCA-37392")
    assert (wbs["chrom"], wbs["start"], wbs["end"]) == ("7", 73330452, 74728173)
    assert clingen.region_named_genes("1q21.1 recurrent region (includes RBM8A and GJA5)") == ["RBM8A", "GJA5"]
    assert clingen.region_named_genes("15q13.3 recurrent region (D-CHRNA7 to BP5) (includes CHRNA7, OTUD7A)") == \
        ["CHRNA7", "OTUD7A"]
    ov = clingen.overlap("7", 100, 199, 150, 300)
    assert ov["relation"] == "partial" and ov["overlap_bp"] == 50 and ov["item_fraction_covered"] == round(50 / 151, 4)
    assert clingen.overlap("7", 100, 400, 150, 300)["relation"] == "cnv_contains_it"
    assert clingen.overlap("7", 160, 170, 150, 300)["relation"] == "cnv_within_it"
    assert clingen.overlap("7", 150, 300, 150, 300)["relation"] == "identical"
    assert clingen.file_date(_clingen_text("ClinGen_region_curation_list_GRCh38.tsv")) == "28 Apr,2026"


# --------------------------------------------------- B-P0-1 / B-P2-2: HGVS

def test_b_p0_1_intronic_offsets_keep_their_sign(fake_cnv_sources):
    """c.6438+1234_7310-567del is DMD exons 45-50 (871 coding bp, out of frame), not 873 bp in frame."""
    from zebra import cnv as C

    p = C.parse("NM_004006.3:c.6438+1234_7310-567del")
    assert (p["cds_start"], p["cds_end"]) == (6439, 7309)
    got = C.card("NM_004006.3:c.6438+1234_7310-567del")
    assert got.result["frame"]["bases"] == 871 and got.result["frame"]["consequence"] == "out of frame"
    assert not any("'?'" in w for w in got.warnings)  # no '?' in the input, no '?' warning
    one = C.card("NM_004006.3:c.6438+100_6614+50del").result  # exon 45 alone
    assert one["frame"]["bases"] == 176 and one["frame"]["consequence"] == "out of frame"
    assert C.parse("NM_004006.3:c.6439-200_7309+300del")["cds_start"] == 6439  # already the coding bases
    # the HGVS uncertain-breakpoint notation gives the same span
    rng = C.parse("NM_004006.3:c.(6438+1_6439-1)_(7309+1_7310-1)del")
    assert (rng["cds_start"], rng["cds_end"]) == (6439, 7309)
    with pytest.raises(UsageError, match="removes no coding base"):
        C.parse("NM_004006.3:c.6438+5_6439-3del")


def test_b_p0_1_with_the_gene_the_hgvs_form_gets_exons_and_restoration(fake_cnv_sources):
    from zebra import cnv as C

    r = C.card("NM_004006.3(DMD):c.6438+1234_7310-567del").result
    assert r["gene"] == "DMD" and r["exons"]["first"] == 45 and r["exons"]["last"] == 50
    assert r["exons"]["whole_exons"] is True and r["frame"]["bases"] == 871
    assert {c["exon"]: c["restores_frame"] for c in r["frame_restoration"]["candidates"]} == {44: False, 51: True}


def test_b_p2_2_the_gene_in_parentheses_is_kept(fake_cnv_sources):
    from zebra import cnv as C

    assert C.parse("NM_004006.3(DMD):c.6439-?_7309+?del")["gene"] == "DMD"
    got = C.card("NM_004006.3(DMD):c.6439-?_7309+?del")
    assert C.case_fields(got.result)["gene"] == "DMD"
    with pytest.raises(UsageError, match="names DMD and --gene says"):
        C.card("NM_004006.3(DMD):c.6439-?_7309+?del", gene="UTRN")


# --------------------------------------------------------- B-P1-1: UTR

def test_b_p1_1_untranslated_exons_are_not_counted_as_coding(fake_cnv_sources):
    from zebra import cnv as C

    utr = C.card("UBE3A exon 2 deletion").result  # ENST00000648336 exon 2 is 5'UTR only
    assert utr["frame"]["bases"] == 0 and utr["frame"]["consequence"] == "no coding sequence"
    assert utr["frame_restoration"] is None
    start = C.card("UBE3A exon 3-4 deletion").result  # exon 3 holds the start codon
    assert start["frame"]["consequence"] == "not assessed" and start["frame"]["modulo_3"] is None
    assert "start codon" in start["frame"]["reading"] and start["frame_restoration"] is None
    assert start["frame"]["bases"] == 62 and start["frame"]["genomic_bases"] == 162
    # the coding length of the DMD canonical transcript is the 3685-codon CDS plus the stop codon
    dmd = C.card("DMD exon 45-50 deletion").result
    assert dmd["transcript"]["coding_length_bp"] == 11058


# ------------------------------------------------------- B-P1-2: SMN1 exons

def test_b_p1_2_smn1_exon_7_is_legacy_numbering_and_a_dosage_result(fake_cnv_sources):
    from zebra import cnv as C

    got = C.card("SMN1 exon 7 deletion")
    r = got.result
    assert r["kind"] == "copy_number" and r["copy_number"] is None
    assert "frame" not in r and "frame_restoration" not in r  # no exon-skipping arithmetic for SMN1
    mapped = r["exon_numbering"]["ensembl"][0]
    assert mapped == {"legacy_exon": "7", "ensembl_ordinal_exon": 8, "transcript": "ENST00000380707", "chrom": "5",
                      "start": 70951941, "end": 70951994, "length_bp": 54}
    assert any("clinical exon 7 is Ensembl exon 8" in w for w in got.warnings)
    assert any("homozygous deletion is 0 copies" in w for w in got.warnings)
    assert C.card("SMN1 exon 7 homozygous deletion").result["copy_number"] == 0
    assert C.card("SMN1 exon 7 heterozygous deletion").result["copy_number"] == 1
    assert C.parse("SMN1 exon 2a deletion")["exon"] == "2a"
    with pytest.raises(UsageError, match="split into 2a and 2b"):
        C.parse("SMN2 exon 2 deletion")
    with pytest.raises(UsageError, match="legacy numbering zebra knows only for SMN1 and SMN2"):
        C.parse("DMD exon 2a deletion")


# --------------------------------------------------- B-P1-3: in-frame loss

def test_b_p1_3_an_in_frame_deletion_gets_no_skip_target(fake_cnv_sources):
    from zebra import cnv as C

    got = C.card("DMD exon 45-47 deletion")
    r = got.result
    assert r["frame"]["consequence"] == "in frame" and r["frame"]["bases"] == 474
    assert r["frame_restoration"] is None and "already in frame" in r["frame_restoration_note"]
    assert "restores the frame" not in got.text and "already in frame" in got.text


# ----------------------------------------- B-P1-4 / B-P2-1: words after it

def test_b_p1_4_the_coordinate_form_reads_or_refuses_what_follows(fake_cnv_sources):
    from zebra import cnv as C

    hg19 = C.parse("chr22:18631365-21800471 hg19 deletion")
    assert hg19["assembly"] == "GRCh37" and hg19["cnv_type"] == "loss"
    got = C.card("chr22:18631365-21800471 hg19 deletion", copies=1)
    assert got.result["assembly"] == "GRCh37"
    mos = C.parse("chr15:23123715-28193120 x1~2")
    assert mos["mosaic"] and mos["copy_number"] is None and mos["copy_number_range"] == "1~2"
    assert mos["cnv_type"] == "loss"
    pct = C.parse("chr15:23123715-28193120 x1 mosaic 30%")
    assert pct["mosaic"] and pct["mosaic_fraction"] == "30%" and pct["copy_number"] == 1
    with pytest.raises(UsageError, match="more than one finding"):
        C.parse("chr15:23123715-28193120 x1, chr16:29580020-30180020 x3")
    with pytest.raises(UsageError, match="which zebra cannot read"):
        C.parse("chr15:23123715-28193120 loss see comment")
    with pytest.raises(UsageError, match="NCBI36/hg18"):
        C.parse("chr15:23123715-28193120 hg18 loss")
    with pytest.raises(UsageError, match="more than one copy number"):
        C.parse("chr15:23123715-28193120 x1 cn=3")
    # the ISCN path reads its mosaic range the same way
    iscn = C.parse("arr[GRCh38] 15q11.2(23123715_23200000)x2~3")
    assert iscn["mosaic"] and iscn["cnv_type"] == "gain" and iscn["copy_number_range"] == "2~3"
    mos_card = C.card("chr15:25200000-25500000 x1~2").result
    assert mos_card["mosaic"] is True and any("mosaic" in c for c in mos_card["caveats"])


def test_b_p2_1_report_typography_and_tags_are_read(fake_cnv_sources):
    from zebra import cnv as C

    assert C.parse("arr[GRCh38] 22q11.21(18648855_21800471)×1")["copy_number"] == 1  # ×
    assert C.parse("chr15:23123715–28193120 loss")["end"] == 28193120  # en dash
    assert C.parse("chr15：23123715-28193120 loss")["start"] == 23123715  # full-width colon
    tagged = C.parse("arr[GRCh38] 22q11.21(18648855_21800471)x1 dn")
    assert tagged["inheritance"] == "de_novo"
    assert C.parse("arr[GRCh38] 22q11.21(18648855_21800471)x1 mat")["inheritance"] == "maternal"
    assert C.card("chr15:25200000-25500000 x1 pat").query["inheritance"] == "paternal"
    with pytest.raises(UsageError, match="tag says paternal and --inheritance says maternal"):
        C.card("chr15:25200000-25500000 x1 pat", inheritance="maternal")
    with pytest.raises(UsageError, match="NCBI36/hg18"):
        C.parse("arr[NCBI36] 22q11.21(17000000_19000000)x1")
    assert C.parse("SMN1 exon 7 homozygous deletion")["copy_number"] == 0
    two = C.parse("FMR1 CGG 30/230")
    assert two["repeat_count"] == "30/230" and [a["low"] for a in two["alleles"]] == [30, 230]
    assert C.parse("FMR1 CGG >200")["alleles"][0]["low"] == 201
    c9 = C.parse("C9orf72 G4C2 800")
    assert c9["motif"] == "GGGGCC" and c9["repeat_count"] == 800


# ----------------------------------------- B-P1-5: validated dosage tables

def test_b_p1_5_a_maintenance_page_is_one_named_failure_not_no_hi_genes(fake_cnv_sources, monkeypatch):
    from zebra import cnv as C
    from zebra.sources import clingen

    def request(url, source, validate=None, **kw):
        body = "<html><body>Scheduled maintenance</body></html>"
        problem = validate(body) if validate else None
        if problem:
            raise SourceError(source, url, 200, problem)
        return Response(url, 200, body, STAMP, False)

    monkeypatch.setattr(clingen, "request", request)
    got = C.card("arr[GRCh38] 7q11.23(73330452_74778226)x1")
    r = got.result
    assert r["clingen_dosage"]["status"] == "unavailable" and r["clingen_regions"]["status"] == "unavailable"
    failures = [w for w in got.warnings if "unavailable" in w and "ClinGen dosage map" in w]
    assert len(failures) == 2  # one per table, not one per gene
    assert "NOT CHECKED" in got.text and "established HI: none" not in got.text
    assert any("means 'not checked', not 'none'" in w for w in got.warnings)


def test_b_p1_5_a_bad_body_is_never_cached_or_served(tmp_path, monkeypatch):
    """Through the real zebra.http.request: a poisoned cache entry is not served, and a bad answer is not cached."""
    from zebra import http
    from zebra.sources import clingen

    monkeypatch.setenv("ZEBRA_CACHE_DIR", str(tmp_path / "c"))
    url = clingen.DOSAGE_TSV["GRCh38"]
    path = http._cache_path(f"GET {url} ")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"stored": __import__("time").time(), "status": 200,
                                "text": "Service temporarily unavailable", "retrieved_at": STAMP}), "utf-8")
    good = gzip.decompress((CNV_FIX / "ClinGen_gene_curation_list_GRCh38.tsv.gz").read_bytes())
    bodies = [b"<html>maintenance</html>"] * 4  # the live server answers an error page on every attempt

    class FakeResp:
        status = 200

        def __init__(self, body):
            self.body = body

        def read(self, *a):
            return self.body

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    class Opener:
        def open(self, req, timeout=None):
            return FakeResp(bodies.pop(0) if bodies else good)

    monkeypatch.setattr(http, "_opener", lambda: Opener())
    monkeypatch.setattr(http, "_wait_for_retry", lambda delay: True)
    with pytest.raises(SourceError, match="HTML page"):
        clingen.fetch_table("gene", "GRCh38")
    assert "Service temporarily unavailable" in path.read_text()  # the bad answer did not replace it either
    assert clingen.table_problem("Service temporarily unavailable", "gene")
    assert "only 3 data rows" in clingen.table_problem(
        "#Gene Symbol\tGene ID\n" + "\n".join("\t".join(["X"] * 23) for _ in range(3)), "gene")
    resp = clingen.fetch_table("gene", "GRCh38")  # the server recovers: the real table, now cached
    assert resp.text.startswith("#ClinGen Gene Curation Results")
    assert json.loads(path.read_text())["text"].startswith("#ClinGen Gene Curation Results")
    assert clingen.dosage("UBE3A", "GRCh38").result["haploinsufficiency"]["score"] == "3"


# ----------------------------------------------- B-P1-6: section 3 bands

def test_b_p1_6_gain_bands_are_not_the_loss_bands(fake_cnv_sources):
    from zebra import cnv as C

    gain = C.card("arr[GRCh38] 7q11.23(73330452_74778226)x3").result["acmg_cnv_inputs"]["section_3_gene_number"]
    assert "0-34, 35-49 and 50+" in gain["clingen_bands"] and "Table 2" in gain["clingen_bands"]
    assert gain["band_for_this_count"]["gain"].startswith("3A")
    assert gain["source"]["pmid"] == "31690835"
    loss = C.card("chr15:25200000-25500000 loss", copies=1).result["acmg_cnv_inputs"]["section_3_gene_number"]
    assert "0-24, 25-34 and 35+" in loss["clingen_bands"] and "Table 1" in loss["clingen_bands"]
    assert C._section3("gain", 40)["band_for_this_count"] == {"gain": "3B (35-49 genes; 0.45 points in the table)"}
    assert C._section3("loss", 40)["band_for_this_count"] == {"loss": "3C (35+ genes; 0.90 points in the table)"}
    assert C._section3("loss", 24)["band_for_this_count"]["loss"].startswith("3A")
    assert C._section3("gain", 50)["band_for_this_count"]["gain"].startswith("3C")
    both = C._section3("unknown", 30)
    assert set(both["band_for_this_count"]) == {"loss", "gain"} and "does not say which" in both["clingen_bands"]
    assert "RefSeq" in both["counted_as"]


# ------------------------------------------------------------ the other P2s

def test_b_p2_3_score_meanings_are_printed(fake_cnv_sources):
    from zebra import cnv as C

    text = C.card("arr[GRCh38] 22q11.21(18648855_21800471)x1").text
    assert "30 autosomal recessive gene; 40 dosage sensitivity unlikely" in text
    assert "ISCA-37446" in text and "HI 3 (Sufficient evidence for dosage pathogenicity)" in text


def test_b_p2_7_cli_takes_related_and_smn2_copies(fake_cnv_sources, capsys):
    from zebra.cli import main

    assert main(["--json", "cnv", "SMN1 copy number 0", "--smn2-copies", "3"]) == 0
    env = json.loads(capsys.readouterr().out)
    assert env["result"]["related_results"] == [{"gene": "SMN2", "copy_number": 3, "exon": None}]
    assert env["result"]["interpretation"]["smn2"]["phenotype_distribution_supportive_care_only"]["SMA II"] == "54%"


def test_b_p2_8_a_repeated_record_is_not_added_twice(fake_cnv_sources, tmp_path, capsys):
    from zebra import case as case_mod
    from zebra.cli import main

    case_dir = tmp_path / "dup"
    case_mod.init(str(case_dir))
    seen = []
    for _ in range(2):
        assert main(["--json", "cnv", "FMR1 CGG 230", "--case", str(case_dir), "--record"]) == 0
        seen.append(json.loads(capsys.readouterr().out)["result"]["recorded_in_case"])
    assert [s["already_recorded"] for s in seen] == [False, True] and seen[0]["id"] == seen[1]["id"]
    assert len(case_mod.load(str(case_dir))["variants"]) == 1
    assert main(["--json", "cnv", "arr[GRCh38] 22q11.21(18648855_21800471)x1", "--case", str(case_dir),
                 "--record"]) == 0
    capsys.readouterr()
    cnv = [v for v in case_mod.load(str(case_dir))["variants"] if v["kind"] == "cnv"][0]
    assert len(cnv["genes"]) == 50 and "62 protein-coding genes spanned" in cnv["note"]
    assert "ISCA-37446" in cnv["note"]


def test_b_p2_9_x_copy_number_uses_the_proband_sex(fake_cnv_sources, tmp_path, capsys):
    from zebra import case as case_mod
    from zebra import cnv as C
    from zebra.cli import main

    assert C.parse("chrX:31000000-31200000 x2", sex="male")["cnv_type"] == "gain"
    assert C.parse("chrX:31000000-31200000 x1", sex="female")["cnv_type"] == "loss"
    assert C.parse("chrX:31000000-31200000 x1", sex="male")["cnv_type"] == "unknown"  # normal male dosage
    assert C.parse("chrX:100000-200000 x1", sex="male")["cnv_type"] == "loss"  # PAR1: diploid in both sexes
    case_dir = tmp_path / "boy"
    case_mod.init(str(case_dir))
    case_mod.set_profile(str(case_dir), sex="male")
    assert main(["--json", "cnv", "chrX:31000000-31200000 x2", "--case", str(case_dir)]) == 0
    env = json.loads(capsys.readouterr().out)
    assert env["result"]["cnv_type"] == "gain" and env["result"]["sex_used"] == "male"
    assert any("taken from the case profile" in w for w in env["warnings"])


def test_b_p2_10_large_and_mtdna_intervals_are_read(fake_cnv_sources):
    """Large intervals are read (summary only) and stay inside the 60,000-character envelope."""
    from zebra import cnv as C

    big = C.card("chr1:1000000-70000000 loss", copies=1)
    r = big.result
    assert r["genes"]["total"] is None and r["clingen_regions"]["status"] == "checked"
    assert r["acmg_cnv_inputs"]["section_3_gene_number"]["protein_coding_genes"] is None
    assert any(e["isca_id"] == "ISCA-37434" for e in r["clingen_regions"]["dosage_curated"])
    assert fake_cnv_sources["overlap"] == []  # no 70 Mb Ensembl walk
    assert len(r["clingen_regions"]["benign_or_population"]) <= C.BENIGN_LIMIT
    assert len(json.dumps(r)) < 45_000
    mt = C.card("chrM:8470-13447 deletion")
    assert mt.result["region"]["chrom"] == "MT" and mt.result["clingen_regions"]["status"] == "not applicable"
    assert any("heteroplasmy" in c for c in mt.result["caveats"])
    with pytest.raises(UsageError, match="runs past the end of chromosome MT"):
        C.card("chrM:8470-20000 deletion")


# ------------------------------------------------- CP1-5: SMN1/SMN2 meaning

def test_cp1_5_smn1_and_smn2_combinations_have_a_cited_meaning(fake_cnv_sources):
    from zebra import cnv as C

    both = C.card("SMN1 0 copies, SMN2 3 copies")
    r = both.result
    assert r["gene"] == "SMN1" and r["copy_number"] == 0
    assert r["related_results"] == [{"gene": "SMN2", "copy_number": 3, "exon": None}]
    interp = r["interpretation"]
    assert "research reference" in interp["label"] and interp["source"]["url"].endswith("NBK1352/")
    assert interp["source"]["revision"] == "Last Revision: February 12, 2026"
    assert "5q spinal muscular atrophy" in interp["smn1"]["reading"]
    dist = interp["smn2"]["phenotype_distribution_supportive_care_only"]
    assert (dist["SMA I"], dist["SMA II"], dist["SMA III/IV"]) == ("15%", "54%", "31%")
    assert "never a prognosis for one person" in dist["not_a_prediction"]
    assert "two, three, or four copies of SMN2" in interp["smn2"]["treatment_timing_as_genereviews_summarises_it"]
    assert r["classification"] is None
    # the flag form gives the same reading; a contradiction is refused
    flag = C.card("SMN1 copy number 0", smn2_copies=3).result
    assert flag["interpretation"]["smn2"]["phenotype_distribution_supportive_care_only"] == dist
    with pytest.raises(UsageError, match="SMN2 copy number 3 and --smn2-copies says 2"):
        C.card("SMN1 0 copies, SMN2 3 copies", smn2_copies=2)
    with pytest.raises(UsageError, match="goes with an SMN1 result"):
        C.card("SMN2 copy number 3", smn2_copies=3)
    with pytest.raises(UsageError, match="goes with an SMN1 copy-number result"):
        C.card("FMR1 CGG 230", smn2_copies=3)
    carrier = C.card("SMN1 copy number 1").result["interpretation"]["smn1"]
    assert "carrier" in carrier["reading"] and "2%-5%" in carrier["reading"]
    two = C.card("SMN1 copy number 2", smn2_copies=2).result["interpretation"]
    assert "1/670" in two["smn1"]["reading"] and "[2+0]" in two["smn1"]["reading"]
    assert "no diagnostic or carrier meaning" in two["smn2"]["reading"]
    four = C.card("SMN1 copy number 0", smn2_copies=4).result["interpretation"]["smn2"]
    assert four["phenotype_distribution_supportive_care_only"]["smn2_copies"] == ">=4"
    five = C.card("SMN1 copy number 0", smn2_copies=5).result["interpretation"]["smn2"]
    assert "five copies of SMN2" in five["five_or_more"] and "deferred until symptom onset" in \
        five["treatment_timing_as_genereviews_summarises_it"]
    alone = C.card("SMN2 copy number 3").result["interpretation"]
    assert "does not diagnose SMA" in alone["reading"]


# ---------------------------------------------- CP1-5: repeat categories

REPEAT_CASES = [
    ("FMR1 CGG 44", ["normal"]), ("FMR1 CGG 45", ["intermediate (gray zone)"]),
    ("FMR1 CGG 54", ["intermediate (gray zone)"]), ("FMR1 CGG 55", ["premutation"]),
    ("FMR1 CGG 200", ["premutation"]), ("FMR1 CGG 201", ["full mutation"]), ("FMR1 CGG >200", ["full mutation"]),
    ("HTT CAG 26", ["normal"]), ("HTT CAG 27", ["intermediate"]), ("HTT CAG 35", ["intermediate"]),
    ("HTT CAG 36", ["reduced penetrance (HD-causing)"]), ("HTT CAG 39", ["reduced penetrance (HD-causing)"]),
    ("HTT CAG 40", ["full penetrance (HD-causing)"]),
    ("DMPK CTG 34", ["normal"]), ("DMPK CTG 35", ["mutable normal (premutation)"]),
    ("DMPK CTG 49", ["mutable normal (premutation)"]), ("DMPK CTG 50", ["full penetrance"]),
    ("DMPK CAG 1200", ["full penetrance"]),  # the same tract read on the other strand
    ("FXN GAA 33", ["normal"]), ("FXN GAA 34", ["intermediate (mutable normal)"]),
    ("FXN GAA 44", ["intermediate (mutable normal)", "borderline"]), ("FXN GAA 66", ["pathogenic (full penetrance)"]),
    ("C9orf72 G4C2 24", ["normal"]), ("C9orf72 G4C2 25", ["uncertain significance (intermediate)"]),
    ("C9orf72 GGGGCC 61", ["pathogenic"]),
]


@pytest.mark.parametrize("text,expected", REPEAT_CASES, ids=[c[0] for c in REPEAT_CASES])
def test_cp1_5_repeat_sizes_fall_in_the_cited_categories(fake_cnv_sources, text, expected):
    from zebra import cnv as C

    r = C.card(text).result
    ref = r["thresholds"]
    assert ref["categories_for_this_result"][0]["categories"] == expected
    assert "research reference" in ref["label"] and ref["source"]["url"].startswith("https://www.ncbi.nlm.nih.gov/books/")
    assert ref["source"]["pmid"] and ref["source"]["revision"] and ref["source"]["retrieved_at"].startswith("2026-10-06")
    assert r["classification"] is None


def test_cp1_5_repeat_edges_disagreements_and_mismatches_are_said(fake_cnv_sources):
    from zebra import cnv as C

    dm35 = C.card("DMPK CTG 35").result["thresholds"]["categories_for_this_result"][0]
    assert "EMQN 2012 lists 35 as normal" in dm35["sources_disagree"]
    assert C.card("DMPK CTG 35").result["thresholds"]["second_source"]["pmid"] == "22643181"
    assert C.card("FMR1 CGG 55").result["thresholds"]["categories_for_this_result"][0]["near_a_boundary"]
    assert "near_a_boundary" not in C.card("FMR1 CGG 30").result["thresholds"]["categories_for_this_result"][0]
    assert C.card("HTT CAG 38").result["thresholds"]["categories_for_this_result"][0]["near_a_boundary"]
    span = C.card("FMR1 CGG 150-300").result["thresholds"]["categories_for_this_result"][0]
    assert span["categories"] == ["premutation", "full mutation"]
    assert "size mosaicism" in span["meanings"][0] and "not as a premutation" in span["meanings"][0]
    pair = C.card("FMR1 CGG 30/230").result["thresholds"]["categories_for_this_result"]
    assert [a["categories"] for a in pair] == [["normal"], ["full mutation"]]
    wrong = C.card("FMR1 CAG 230")
    assert wrong.result["thresholds"]["categories_for_this_result"] is None
    assert any("FMR1's repeat is CGG" in w for w in wrong.warnings)
    other = C.card("ATXN3 CAG 70").result
    assert other["thresholds"] is None and "no retrieved threshold table for ATXN3" in other["thresholds_note"]
    assert C.card("C9orf72 G4C2 800").result["gene"] == "C9orf72"


def test_cp1_5_every_category_quote_carries_its_own_numbers():
    """The quoted source sentence must contain the numbers the category uses (no number without its sentence)."""
    from zebra import cnv as C

    for gene, locus in C.REPEAT_LOCI.items():
        assert locus["source"] in C.SOURCES
        for cat in locus["categories"]:
            nums = {str(cat["min"])} if cat["max"] is None else {str(cat["max"])}
            if gene == "HTT" and cat["name"] == "normal":
                nums = {"26"}
            if gene == "FMR1" and cat["name"] == "full mutation":
                nums = {"200"}
            if gene == "DMPK" and cat["name"] == "full penetrance":
                nums = {"50"}
            if gene == "C9ORF72" and cat["name"] == "pathogenic":
                nums = {"61"}
            if gene == "FXN" and cat["max"] is None:
                nums = {"66"}
            assert any(n in cat["quote"] for n in nums), (gene, cat["name"])
    for src in C.SOURCES.values():
        assert src["url"].startswith("https://") and src["pmid"] and src["retrieved_at"].startswith("2026-10-06")


# ------------------------------------------------------------------- live

@pytest.mark.live
@pytest.mark.parametrize("label,iscn,isca,name,hi,relation", ACCEPTANCE, ids=[a[0] for a in ACCEPTANCE])
def test_live_cp0_3_acceptance_regions(label, iscn, isca, name, hi, relation):
    from zebra import cnv as C

    r = C.card(iscn).result
    entry = next(e for e in r["clingen_regions"]["dosage_curated"] if e["isca_id"] == isca)
    assert name in entry["name"] and entry["haploinsufficiency"]["score"] == hi and entry["relation"] == relation
    assert r["clingen_dosage"]["status"] == "checked"


@pytest.mark.live
def test_live_cp0_3_gene_desert_and_grch37():
    from zebra import cnv as C

    desert = C.card("chr13:55000000-56000000 loss", copies=1).result
    assert desert["clingen_regions"]["status"] == "checked" and desert["clingen_regions"]["dosage_curated"] == []
    assert desert["genes"]["protein_coding"] == 0
    hg19 = C.card("arr[hg19] 22q11.21(18648855_21800471)x1").result
    assert any(e["isca_id"] == "ISCA-37446" and e["relation"] == "cnv_contains_it"
               for e in hg19["clingen_regions"]["dosage_curated"])


@pytest.mark.live
def test_live_clingen_tables_are_valid_for_both_builds():
    from zebra.sources import clingen

    for kind in ("gene", "region"):
        for build in ("GRCh38", "GRCh37"):
            resp = clingen.fetch_table(kind, build)
            assert clingen.table_problem(resp.text, kind) is None


@pytest.mark.live
def test_live_b_p0_1_and_b_p1_2_against_ensembl():
    from zebra import cnv as C

    dmd = C.card("NM_004006.3(DMD):c.6438+1234_7310-567del").result
    assert dmd["frame"]["bases"] == 871 and dmd["exons"]["first"] == 45 and dmd["exons"]["last"] == 50
    smn = C.card("SMN1 exon 7 homozygous deletion", smn2_copies=2).result
    assert smn["exon_numbering"]["ensembl"][0]["ensembl_ordinal_exon"] == 8
    assert smn["exon_numbering"]["ensembl"][0]["length_bp"] == 54
    assert C.card("UBE3A exon 2 deletion").result["frame"]["consequence"] == "no coding sequence"


# ------------------------------- adversarial review (round 1, engineering)

def test_adv_a_p0_1_a_thousands_comma_in_a_repeat_count_is_never_two_alleles(fake_cnv_sources):
    from zebra import cnv as C

    for text in ("C9orf72 G4C2 1,000", "DMPK CTG 1,200", "FXN GAA 1,100"):
        with pytest.raises(UsageError, match="one number with a thousands comma or two alleles"):
            C.parse(text)
    assert C.parse("DMPK CTG 1200")["repeat_count"] == 1200
    assert C.parse("FMR1 CGG 30, 230")["repeat_count"] == "30/230"
    assert C.parse("FMR1 CGG 30/230")["repeat_count"] == "30/230"


def test_adv_a_p0_2_the_x_pseudoautosomal_bounds_follow_the_build(fake_cnv_sources):
    from zebra import cnv as C

    # GRCh37 PAR1 ends at 2,699,520: this interval is X-specific there (one copy is normal in a male)
    p37 = C.parse("chrX:2700000-2780000 x1", sex="male", assembly="GRCh37")
    assert p37["in_par"] is False and p37["cnv_type"] == "unknown"
    assert C.parse("chrX:2700000-2780000 x1", sex="male", assembly="GRCh38")["cnv_type"] == "loss"  # PAR1 in GRCh38
    # GRCh37 PAR2 starts at 154,931,044: inside it two copies are normal in a male too
    assert C.parse("chrX:155000000-155100000 x1", sex="male", assembly="GRCh37")["cnv_type"] == "loss"
    got = C.card("chrX:2700000-2780000 x1", assembly="GRCh37", sex="male").result
    assert got["cnv_type"] == "unknown" and not any("pseudoautosomal" in c for c in got["caveats"])


def test_adv_a_p1_1_a_report_word_the_count_cannot_mean_is_refused(fake_cnv_sources):
    from zebra import cnv as C

    for text in ("chr15:25200000-25500000 loss x2", "chrX:31000000-31200000 x1 dup", "chrY:3000000-4000000 x2 loss"):
        with pytest.raises(UsageError, match="cannot be a"):
            C.parse(text)
    with pytest.raises(UsageError, match="cannot be a loss"):
        C.card("chr15:25200000-25500000 loss", copies=2)
    # a word the count allows in one of the two sexes stays (X x1 is a loss in a female)
    assert C.parse("chrX:31000000-31200000 x1 del")["cnv_type"] == "loss"


def test_adv_a_p1_2_an_smn_deletion_word_and_a_copy_number_must_agree(fake_cnv_sources):
    from zebra import cnv as C

    with pytest.raises(UsageError, match="a deletion leaves 0 or 1 copies"):
        C.card("SMN1 exon 7 deletion", copies=3)
    with pytest.raises(UsageError, match="a duplication 3 or more"):
        C.card("SMN1 exon 7 duplication", copies=1)
    ok = C.card("SMN1 exon 7 deletion", copies=0)
    assert ok.result["copy_number"] == 0 and not any("without saying" in w for w in ok.warnings)


def test_adv_a_p1_3_an_hgvs_bracket_spanning_exons_is_refused(fake_cnv_sources):
    from zebra import cnv as C

    with pytest.raises(UsageError, match="each bracket must name one intron"):
        C.parse("NM_004006.3(DMD):c.(6438+1_6615-1)_(6912+1_6913-1)del")


def test_adv_a_p1_4_an_enst_of_another_gene_is_named_not_swapped_silently(fake_cnv_sources):
    from zebra import cnv as C

    got = C.card("ENST00000380707(DMD):c.6439-?_7309+?del")
    assert any("not one of Ensembl's DMD transcripts" in w for w in got.warnings)
    exon = C.card("ENST00000380707 exon 45-50 deletion", gene="DMD")
    assert any("not one of Ensembl's DMD transcripts" in w for w in exon.warnings)


def test_adv_a_p1_5_a_chromosome_or_unknown_symbol_is_not_a_gene(fake_cnv_sources):
    from zebra import cnv as C

    for text in ("chrX x1", "chr21 x3", "X copy number 1", "chr7 exon 2 deletion"):
        with pytest.raises(UsageError, match="names a chromosome where a gene belongs"):
            C.parse(text)
    with pytest.raises(UsageError, match="Ensembl has no gene with the symbol 'NOSUCHGENE1'"):
        C.card("NOSUCHGENE1 copy number 2")
    with pytest.raises(UsageError, match="Ensembl has no gene with the symbol 'NOSUCHGENE1'"):
        C.card("NOSUCHGENE1 exon 2 deletion")


def test_adv_a_p1_6_an_ensembl_overlap_body_is_validated_before_caching(fake_cnv_sources, monkeypatch):
    from zebra import cnv as C

    assert C._overlap_problem('[{"id": "ENSG1"') == "a truncated or malformed JSON list"
    assert "error" in C._overlap_problem('{"error": "Too many requests"}')
    assert "error page" in C._overlap_problem("<html>busy</html>")
    assert C._overlap_problem("[]") is None
    seen = {}
    real = C.get_json

    def spy(url, source, **kw):
        seen["validate"] = kw.get("validate")
        return real(url, source, **kw)

    monkeypatch.setattr(C, "get_json", spy)
    C.genes_in_region("13", 55000000, 56000000)
    assert seen["validate"] is C._overlap_problem


def test_adv_a_p1_7_an_ensembl_failure_keeps_the_clingen_answer(fake_cnv_sources, monkeypatch):
    from zebra import cnv as C

    def down(url, source, **kw):
        raise SourceError(source, url, None, "deadline reached")

    monkeypatch.setattr(C, "get_json", down)
    got = C.card("arr[GRCh38] 7q11.23(73330452_74778226)x1")
    r = got.result
    assert r["genes"]["total"] is None and any("Ensembl overlap unavailable" in w for w in got.warnings)
    assert r["clingen_regions"]["status"] == "checked"
    assert any(e["isca_id"] == "ISCA-37392" for e in r["clingen_regions"]["dosage_curated"])


def test_adv_a_p1_7_a_partial_gene_walk_is_a_lower_bound(fake_cnv_sources, monkeypatch):
    from zebra import cnv as C

    real = C.get_json

    def second_window_fails(url, source, **kw):
        if url.endswith("1:5349466-6823542"):
            raise SourceError(source, url, None, "deadline reached")
        return real(url, source, **kw)

    monkeypatch.setattr(C, "get_json", second_window_fails)
    got = C.card("arr[GRCh38] 1p36.33p36.31(849466_6823542)x1")
    r = got.result
    assert r["genes"]["complete"] is False and any("lower bound" in w for w in got.warnings)
    sec3 = r["acmg_cnv_inputs"]["section_3_gene_number"]
    assert sec3["protein_coding_genes"] is None and sec3["protein_coding_genes_at_least"] >= 35
    assert sec3["band_for_this_count"] == {"loss": "3C (35+ genes; 0.90 points in the table)"}  # already certain
    assert r["clingen_dosage"]["genes_checked"] is None


def test_adv_a_p1_8_a_different_finding_is_not_taken_for_one_already_recorded(fake_cnv_sources, tmp_path, capsys):
    from zebra import case as case_mod
    from zebra.cli import main

    case_dir = tmp_path / "ident"
    case_mod.init(str(case_dir))
    for text in ("chr15:25200000-25500000 loss", "chr15:25200000-25500000 x1~2", "chr15:25200000-25500000 loss dn",
                 "DMD exon 45 deletion", "NM_004006.3(DMD):c.6450_6600del", "FMR1 CGG 30/230", "FMR1 CGG 230/30"):
        assert main(["--json", "cnv", text, "--case", str(case_dir), "--record"]) == 0
        capsys.readouterr()
    kinds = [v["kind"] for v in case_mod.load(str(case_dir))["variants"]]
    assert kinds.count("cnv") == 3 and kinds.count("exon_cnv") == 2 and kinds.count("repeat_expansion") == 1


def test_adv_a_p1_9_what_was_not_checked_is_null_not_empty(fake_cnv_sources, monkeypatch):
    from zebra import cnv as C
    from zebra.sources import clingen

    def request(url, source, validate=None, **kw):
        raise SourceError(source, url, 200, validate("<html>maintenance</html>"))

    monkeypatch.setattr(clingen, "request", request)
    r = C.card("arr[GRCh38] 7q11.23(73330452_74778226)x1").result
    sec2 = r["acmg_cnv_inputs"]["section_2_overlap_with_established_regions_or_genes"]
    assert sec2["clingen_established_regions"] is None and sec2["clingen_established_haploinsufficient_genes"] is None
    assert r["clingen_regions"]["dosage_curated"] is None and r["clingen_dosage"]["records"] is None
    assert r["clingen_regions"]["other_overlapping_without_dosage_evidence"] is None


def test_adv_a_p1_10_malformed_numbers_are_usage_errors(fake_cnv_sources, capsys):
    from zebra import cnv as C
    from zebra.cli import main

    for text in ("chr15:,-28193120 loss", "chr1:1,000,000-2,000,0000 loss", "chr1:1,000,000-20,00,000 loss"):
        with pytest.raises(UsageError):
            C.parse(text)
    assert main(["--json", "cnv", "arr[GRCh38] 15q11.2(,_,)x1"]) == 2
    assert json.loads(capsys.readouterr().out)["error"]["type"] == "UsageError"
    assert C.parse("chr1:1,000,000-2,000,000 loss")["end"] == 2_000_000


def test_adv_a_p2_smaller_input_defects(fake_cnv_sources, tmp_path, capsys):
    from zebra import case as case_mod
    from zebra import cnv as C
    from zebra.cli import main

    assert C.parse("arr[GRCh38] 1Q21.1(146000000_147000000)x1")["chrom"] == "1"
    with pytest.raises(UsageError):
        C.parse("DMD exons 3 8 deletion")
    with pytest.raises(UsageError):
        C.parse("DMD exon 0012 deletion")
    with pytest.raises(UsageError, match="there is no exon 9"):
        C.parse("SMN1 exon 9 copy number 0")
    with pytest.raises(UsageError, match="cannot read"):
        C.parse("1 45-50 deletion")
    with pytest.raises(UsageError, match="more than one finding"):
        C.parse("arr[GRCh38] 15q11.2(1000_2000)x1, arr[GRCh38] 16p11.2(1000_2000)x3")
    with pytest.raises(UsageError, match="is not a range"):
        C.parse("chr15:25200000-25500000 x1~1")
    with pytest.raises(UsageError, match="not a fraction of cells"):
        C.parse("chr15:25200000-25500000 x1 mosaic 300%")
    with pytest.raises(UsageError, match="contradict"):
        C.parse("chr15:25200000-25500000 x1 hom")
    with pytest.raises(UsageError, match="--related takes further copy-number results"):
        C.card("chr15:25200000-25500000 loss", related=["SMN2 copy number 2"])
    with pytest.raises(UsageError, match="--gene names the gene of an exon-level result"):
        C.card("FMR1 CGG 230", gene="DMD")
    assert C.card("chrM:100-200 x1~2").result["cnv_type"] == "unknown"
    # GABRD lies outside this part of 1p36: said, not "?"
    r = C.card("chr1:5349466-6823542 loss", copies=1).result
    gabrd = next(g for e in r["clingen_regions"]["dosage_curated"] if e["isca_id"] == "ISCA-37434"
                 for g in e["named_genes"] if g["gene"] == "GABRD")
    assert gabrd["in_cnv"] == "no"
    # a case profile written "M" is a male proband
    case_dir = tmp_path / "m"
    case_mod.init(str(case_dir))
    case_mod.set_profile(str(case_dir), sex="M")
    assert main(["--json", "cnv", "chrX:31000000-31200000 x2", "--case", str(case_dir)]) == 0
    assert json.loads(capsys.readouterr().out)["result"]["cnv_type"] == "gain"
    # a failed write keeps the analysis and says it was not recorded
    bad = tmp_path / "bad"
    case_mod.init(str(bad))
    (bad / "case.json").write_text("{not json", "utf-8")
    assert main(["--json", "cnv", "FMR1 CGG 230", "--case", str(bad), "--record"]) == 0
    env = json.loads(capsys.readouterr().out)
    assert env["result"]["recorded_in_case"]["recorded"] is False and env["result"]["thresholds"]
    assert any("NOT recorded" in w for w in env["warnings"])


def test_e_6_clingen_validity_writes_the_good_answer_over_a_bad_cached_body(tmp_path, monkeypatch):
    """A bad cached validity CSV is replaced in the cache, not just read past (W4's E-6, in zebra/sources/clingen.py)."""
    import time

    from zebra import http
    from zebra.sources import clingen

    monkeypatch.setenv("ZEBRA_CACHE_DIR", str(tmp_path / "c"))
    path = http._cache_path(f"GET {clingen.VALIDITY_CSV} ")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"stored": time.time(), "status": 200, "text": "<html>Service unavailable</html>",
                                "retrieved_at": STAMP}), "utf-8")
    good = (Path(__file__).parent / "fixtures" / "clingen" / "validity_excerpt.csv").read_bytes()
    requests = []

    class FakeResp:
        status = 200

        def read(self, *a):
            return good

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    class Opener:
        def open(self, req, timeout=None):
            requests.append(req.full_url)
            return FakeResp()

    monkeypatch.setattr(http, "_opener", lambda: Opener())
    rows = clingen.validity(symbol="SCN1A").result
    assert rows and rows[0]["gene"] == "SCN1A"
    assert "GENE SYMBOL" in json.loads(path.read_text())["text"]  # the good body replaced the bad one
    assert clingen.validity(symbol="SCN1A").result == rows and len(requests) == 1  # now served from the cache
    assert clingen.validity_problem("<html>x</html>") and clingen.validity_problem("GENE SYMBOL,x") is None


# ------------------------------- adversarial review (round 1, medical fidelity)

def test_adv_b_p0_1_an_fmr1_full_mutation_size_is_not_read_without_methylation(fake_cnv_sources):
    from zebra import cnv as C

    got = C.card("FMR1 CGG 210")
    allele = got.result["thresholds"]["categories_for_this_result"][0]
    joined = " ".join(allele["meanings"])
    assert "methylation result is part of the answer" in joined
    assert "in a male: in GeneReviews Table 3, males with a completely methylated full mutation" in joined
    assert "in a female" in joined and any("about 230" in n for n in allele["range_notes"])
    assert "methylation" in got.text and "about 230" in got.text  # printed, not only in the JSON
    male = C.card("FMR1 CGG 210", sex="male").result["thresholds"]["categories_for_this_result"][0]["meanings"]
    assert not any("in a female" in m for m in male)


def test_adv_b_p0_2_an_smn1_count_of_another_exon_gets_no_exon_7_reading(fake_cnv_sources):
    from zebra import cnv as C

    r = C.card("SMN1 exon 8 copy number 0").result["interpretation"]
    assert "no diagnostic or carrier reading" in r["smn1"]["reading"]
    assert "5q spinal muscular atrophy" not in r["smn1"]["reading"]
    assert "molecular criterion" in C.card("SMN1 exon 7 copy number 0").result["interpretation"]["smn1"]["reading"]
    assert "molecular criterion" in C.card("SMN1 exon 7-8 copy number 0").result["interpretation"]["smn1"]["reading"]


def test_adv_b_p0_3_fxn_is_read_across_both_alleles(fake_cnv_sources):
    from zebra import cnv as C

    def both(text):
        return C.card(text).result["thresholds"]["reading_of_both_alleles"]

    assert "biallelic" in both("FXN GAA 800/800") and "96%" in both("FXN GAA 800/800")
    assert "LOFA/VLOFA" in both("FXN GAA 50/800")
    assert "sequence analysis of FXN" in both("FXN GAA 8/800")  # P1-4: the next step for one expansion
    assert "whether it is one allele" in both("FXN GAA 800")
    per_allele = C.card("FXN GAA 800/800").result["thresholds"]["categories_for_this_result"]
    assert not any("carrier" in m for a in per_allele for m in a["meanings"])


def test_adv_b_p1_5_notes_reach_the_text_output(fake_cnv_sources):
    from zebra import cnv as C

    assert "underestimate HTT CAG length by two repeats" in C.card("HTT CAG 35").text
    assert "EMQN" in C.card("DMPK CTG 50").text


def test_adv_b_p1_6_fmr1_readings_follow_the_sex(fake_cnv_sources):
    from zebra import cnv as C

    male = C.card("FMR1 CGG 80", sex="male")
    text = " ".join(male.result["thresholds"]["categories_for_this_result"][0]["meanings"])
    assert "all of his daughters and to none of his sons" in text and "older than 50" in text
    assert "in a female" not in text
    female = " ".join(C.card("FMR1 CGG 30/80", sex="female").result["thresholds"]["categories_for_this_result"][1]
                      ["meanings"])
    assert "FXPOI" in female and "children with FXS" in female and "in a male" not in female
    two = C.card("FMR1 CGG 30/80", sex="male")
    assert any("a male has one allele" in w for w in two.warnings)


def test_adv_b_p1_7_an_fmr1_range_across_premutation_and_full_is_not_a_premutation(fake_cnv_sources):
    from zebra import cnv as C

    allele = C.card("FMR1 CGG 200-800").result["thresholds"]["categories_for_this_result"][0]
    assert "size mosaicism" in allele["meanings"][0] and len(allele["meanings"]) == 1
    assert "not associated with FXS" not in " ".join(allele["meanings"])


def test_adv_b_p1_8_to_p1_10_smn_readings_keep_their_scope(fake_cnv_sources):
    from zebra import cnv as C

    two = C.card("SMN1 copy number 2").result["interpretation"]["smn1"]["reading"]
    assert "NOT known to have a family history" in two and "approximately 6% of parents" in two
    assert any("two-copy result (apparently not a carrier)" in m
               for m in C.card("SMN1 copy number 2").result["mechanism"])
    timing = C.card("SMN1 copy number 0", smn2_copies=5).result["interpretation"]["smn2"]
    assert timing["treatment_timing_as_genereviews_summarises_it"].startswith(
        "for an infant detected by newborn screening, after confirmatory SMN1 testing")
    one = C.card("SMN1 copy number 1").result["interpretation"]["smn1"]["reading"]
    assert "consistent with carrier status" in one and "the result of a carrier" not in one


def test_adv_b_p1_11_to_p1_13_apparent_homozygosity_and_low_penetrance_are_said(fake_cnv_sources):
    from zebra import cnv as C

    assert "apparently homozygous" in C.card("HTT CAG 17/17").result["thresholds"]["apparent_homozygosity"]
    assert "apparent homozygosity" in C.card("C9orf72 G4C2 5").result["thresholds"]["apparent_homozygosity"]
    assert "apparent_homozygosity" not in C.card("HTT CAG 17/20").result["thresholds"]
    reduced = " ".join(C.card("HTT CAG 37").result["thresholds"]["categories_for_this_result"][0]["meanings"])
    assert "0.2%-2%" in reduced and "significantly higher in families with HD" in reduced
    c9 = " ".join(C.card("C9orf72 G4C2 47").result["thresholds"]["categories_for_this_result"][0]["meanings"])
    assert "cosegregate with the disorder in a family was 47" in c9


def test_adv_b_p1_14_a_region_is_read_by_the_score_of_this_cnv_type(fake_cnv_sources, monkeypatch):
    """15q11.2 BP1-BP2 (ISCA-37448): HI 3 but TS 40 — for a gain it is a population (benign) region."""
    from zebra import cnv as C

    gain = C.card("chr15:22782170-23040134 dup", copies=3).result
    assert not any(e["isca_id"] == "ISCA-37448" for e in gain["clingen_regions"]["dosage_curated"])
    bp = next(e for e in gain["clingen_regions"]["benign_or_population"] if e["isca_id"] == "ISCA-37448")
    assert bp["triplosensitivity_score"] == "40" and bp["acmg_section_2"].startswith("2C-2G")
    loss = C.card("chr15:22782170-23040134 del", copies=1).result
    assert any(e["isca_id"] == "ISCA-37448" for e in loss["clingen_regions"]["dosage_curated"])


def test_adv_b_p1_15_a_homozygous_loss_names_the_recessive_genes(fake_cnv_sources):
    from zebra import cnv as C

    r = C.card("chr15:25200000-25500000 x0").result  # UBE3A only; no score-30 gene here: no note
    assert "recessive_genes_lost" not in r
    pws = C.card("arr[GRCh38] 15q11.2q13.1(23123715_28193120)x0").result
    rec = pws["recessive_genes_lost"]
    assert rec["copy_number"] == 0 and "OCA2" in {g["gene"] for g in rec["genes"]}
    assert "lost on both alleles" in rec["reading"]
    het = C.card("arr[GRCh38] 15q11.2q13.1(23123715_28193120)x1").result["recessive_genes_lost"]
    assert "carrier status" in het["reading"]
    benign = [e for e in pws["clingen_regions"]["benign_or_population"]]
    assert all("do not speak to a homozygous loss" in e["acmg_section_2"] for e in benign if e["acmg_section_2"])


def test_adv_b_p1_16_a_female_x_loss_carries_the_x_linked_note(fake_cnv_sources, monkeypatch):
    from zebra import cnv as C
    from zebra.sources import clingen

    # DMD (HI 3) by position, with no Ensembl capture for the window: the note rests on ClinGen alone
    female = C.card("chrX:31800000-31990000 deletion", copies=1, sex="female")
    assert any("females may manifest symptoms" in c for c in female.result["caveats"])
    assert "females may manifest symptoms" in female.text
    male = C.card("chrX:31800000-31990000 deletion", copies=0, sex="male")
    assert not any("females may manifest" in c for c in male.result["caveats"])


def test_adv_b_p1_17_the_missing_gene_list_names_its_real_reason(fake_cnv_sources, monkeypatch):
    from zebra import cnv as C

    def down(url, source, **kw):
        raise SourceError(source, url, None, "deadline reached")

    monkeypatch.setattr(C, "get_json", down)
    r = C.card("arr[GRCh38] 7q11.23(73330452_74778226)x1").result
    assert r["genes"]["note"] == "not listed: Ensembl did not answer (see warnings)"


def test_adv_b_p1_18_the_exon_text_carries_the_arithmetic_only_caveats(fake_cnv_sources):
    from zebra import cnv as C

    text = C.card("DMD exon 45-50 deletion").text
    assert "arithmetic only" in text and "separate question" in text and "caveat: frame arithmetic only" in text
    assert C.card("DMD exon 52 deletion").text.splitlines()[0] == "DMD loss exon 52"
