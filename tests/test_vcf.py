"""zebra vcf: reading, genotypes, inheritance classes, restriction, ranking (offline with captured Ensembl
responses), plus one live end-to-end triage."""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path

import pytest

from zebra import vcf as V
from zebra.core import Outcome, UsageError
from zebra.http import Response, SourceError, source_record
from zebra.sources import ensembl

FIX = Path(__file__).parent / "fixtures" / "vcf"
GZ = str(FIX / "trio.vcf.gz")
PLAIN = str(FIX / "trio.vcf")
FAM2 = str(FIX / "family2.vcf")
FAM2_GZ = str(FIX / "family2.vcf.gz")
BAD = str(FIX / "bad.vcf")
REAL_HPO = Path(os.path.expanduser("~/.cache/zebra-mod/hpo"))


# ------------------------------------------------------------------ mocks

@pytest.fixture(autouse=True)
def _moi_rows(monkeypatch):
    """Gene inheritance modes from a trimmed copy of HPO genes_to_phenotype.txt (release 2026-09-01)."""
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(FIX / "hpo"))


@pytest.fixture
def fake_ensembl(monkeypatch):
    """VEP and lookup answered from responses captured from rest.ensembl.org (see make_fixtures.py).

    Both families' captures are served: the original trio and the second family
    added for the triage defects (vep_family2.json).
    """
    by_input = {}
    for name in ("vep_trio.json", "vep_family2.json"):
        for r in json.loads((FIX / name).read_text())["records"]:
            by_input[r["input"]] = r
    lookup = json.loads((FIX / "lookup_genes.json").read_text())["response"]
    calls = {"vep": [], "lookup": []}

    def vep_batch(variants, assembly="GRCh38"):
        calls["vep"].append(list(variants))
        lines = [ensembl.vcf_line(*v) for v in variants]
        return Outcome([by_input[l] for l in lines if l in by_input],
                       sources=[source_record("Ensembl VEP", f"fixture ({len(lines)} variants)", url="https://rest.ensembl.org/vep/human/region")])

    def post_json(url, payload, source, **kw):
        calls["lookup"].append(list(payload["symbols"]))
        upper = {k.upper(): v for k, v in lookup.items()}  # Ensembl matches case-insensitively, keys = input
        data = {s: upper[s.upper()] for s in payload["symbols"] if s.upper() in upper}
        return Response(url, 200, json.dumps(data), "2026-10-05T15:24:00+00:00", False)

    monkeypatch.setattr(ensembl, "vep_batch", vep_batch)
    monkeypatch.setattr(V, "post_json", post_json)
    return calls


def _genes():
    return V.read_gene_file(str(FIX / "genes.txt"))


def _by_gene(res):
    return {c["gene"]: c for c in res["candidates"]}


# ------------------------------------------------------------- genotypes

def test_parse_gt_forms():
    assert V.parse_gt("0/1") == ((0, 1), False)
    assert V.parse_gt("0|1") == ((0, 1), True)
    assert V.parse_gt("1") == ((1,), False)
    assert V.parse_gt("./.") == ((None, None), False)
    assert V.parse_gt(".") == (None, False)
    assert V.parse_gt("1/2") == ((1, 2), False)


@pytest.mark.parametrize("raw,k,expected", [
    ("0/1", 1, "het"), ("0|1", 1, "het"), ("1/1", 1, "hom_alt"), ("1", 1, "hemi"), ("0", 1, "hom_ref"),
    ("0/0", 1, "hom_ref"), ("./.", 1, "missing"), (".", 1, "missing"), ("./1", 1, "het"), ("./0", 1, "missing"),
    ("0/2", 1, "other"), ("0/2", 2, "het"), ("1/2", 1, "het"), ("1/2", 2, "het"), ("2/2", 2, "hom_alt"),
])
def test_zygosity(raw, k, expected):
    assert V.zygosity(V.parse_call(["GT"], raw), k) == expected


def test_call_depth_quality_and_allele_balance():
    c = V.parse_call(["GT", "AD", "DP", "GQ"], "0/1:18,17:35:99")
    assert (c.dp, c.gq, c.ad) == (35, 99.0, [18, 17])
    assert c.ab(1) == pytest.approx(17 / 35)
    c2 = V.parse_call(["GT", "AD"], "0/1:20,18,2")  # no DP: depth from AD; multi-allelic AD
    assert c2.depth() == 40 and c2.reads(2) == (20, 2)
    c3 = V.parse_call(["GT", "AD", "DP", "GQ"], "0/1:.:.:.")
    assert c3.depth() is None and c3.gq is None and c3.ab(1) is None


def test_normalisation():
    assert V.norm_chrom("chr2") == "2" and V.norm_chrom("chrM") == "MT" and V.norm_chrom("chrx") == "X"
    assert V.norm_chrom("MT") == "MT" and V.norm_chrom("17") == "17"
    assert V.normalize_allele(100, "ATT", "AT") == (100, "AT", "A")  # split multi-allelic leftovers
    assert V.normalize_allele(100, "CTG", "CAG") == (101, "T", "A")
    assert V.normalize_allele(100, "g", "a") == (100, "G", "A")
    assert V.is_symbolic("*") and V.is_symbolic("<DEL>") and not V.is_symbolic("A")


# ---------------------------------------------------------------- inspect

def test_inspect_bgzip_and_plain():
    gz = V.inspect(GZ).result
    plain = V.inspect(PLAIN).result
    assert V.compression(GZ) == "bgzip" and gz["compression"] == "bgzip" and plain["compression"] == "plain"
    for r in (gz, plain):
        assert r["samples"] == ["M", "F", "P", "S"]
        assert r["variants"] == 13 and r["scan_complete"] and r["multiallelic_records"] == 1
        assert r["build"]["guess"] == "GRCh38" and not r["build"]["conflict"]
        assert r["genotypes"] is True
        assert r["filters"] == {"PASS": 12, "LowQual": 1}
        assert r["annotations"]["gnomad_af_info_keys"] == []  # INFO/AF is the call set's own frequency
        assert r["per_sample"]["S"]["sex_hint"].startswith("male-like")
        assert r["per_sample"]["F"]["sex_hint"].startswith("male-like")  # haploid ref call on X
    assert gz["contig_naming"].startswith("plain") and plain["contig_naming"].startswith("chr")
    assert set(gz["chromosomes"]) == set(plain["chromosomes"]) == {"2", "3", "5", "7", "X"}


def test_gzip_without_bgzf_blocks(tmp_path):
    p = tmp_path / "x.vcf.gz"
    with gzip.open(p, "wb") as fh:
        fh.write(Path(PLAIN).read_bytes())
    assert V.compression(str(p)) == "gzip"
    assert V.inspect(str(p)).result["variants"] == 13


def _header_with(contigs, reference=None):
    lines = ["##fileformat=VCFv4.2"] + [f"##contig=<ID={c},length={n}>" for c, n in contigs]
    if reference:
        lines.append(f"##reference={reference}")
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP")
    import io

    h, _ = V.read_header(io.StringIO("\n".join(lines) + "\n"))
    return h


def test_build_guess():
    assert V.guess_build(_header_with([("1", 249250621), ("X", 155270560)]))["guess"] == "GRCh37"
    assert V.guess_build(_header_with([("chr1", 248956422)]))["guess"] == "GRCh38"
    assert V.guess_build(_header_with([], reference="hs37d5.fa"))["guess"] == "GRCh37"
    assert V.guess_build(_header_with([]))["guess"] is None
    conflict = V.guess_build(_header_with([("1", 248956422)], reference="human_g1k_v37.fasta"))
    assert conflict["guess"] is None and conflict["conflict"]


def test_header_parses_csq_ann_and_gnomad_keys():
    import io

    text = "\n".join([
        "##fileformat=VCFv4.2",
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|gnomADe_AF|gnomADe_ASJ_AF|MAX_AF">',
        "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID'\">",
        '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="gnomAD v4 AF">',
        '##INFO=<ID=gnomAD_AF_asj,Number=A,Type=Float,Description="gnomAD v4 AF, Ashkenazi">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="cohort AF">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]) + "\n"
    h, first = V.read_header(io.StringIO(text))
    assert first is None and h.samples == [] and not h.has_format
    assert h.csq_fields[:4] == ["Allele", "Consequence", "IMPACT", "SYMBOL"]
    assert h.ann_fields[:4] == ["Allele", "Annotation", "Annotation_Impact", "Gene_Name"]
    assert h.gnomad_info_keys() == ["gnomAD_AF"]  # bottlenecked group and cohort AF excluded
    assert V._is_af_field("gnomADe_AF") and not V._is_af_field("gnomADe_ASJ_AF") and not V._is_af_field("MAX_AF")


# ----------------------------------------------------------- inheritance

def _p(z, adequate=False, why=()):
    return {"z": z, "adequate_ref": adequate, "why": list(why)}


def test_classify_trio_classes():
    ok = _p("hom_ref", True)
    assert V._classify("auto", "female", "het", ok, ok)[:2] == ("de_novo", "de_novo")
    low = _p("hom_ref", False, ["DP 4 < 10"])
    cls, _, flags = V._classify("auto", "female", "het", low, ok)
    assert cls == "possible_de_novo" and "mother: DP 4 < 10" in flags[0]
    assert V._classify("auto", None, "het", _p("het"), ok)[:2] == ("inherited_het", "maternal")
    assert V._classify("auto", None, "het", _p("het"), _p("het"))[:2] == ("inherited_het", "both")
    assert V._classify("auto", None, "hom_alt", _p("het"), _p("missing"))[0] == "hom_recessive"
    cls, _, flags = V._classify("auto", None, "hom_alt", _p("het"), ok)
    assert cls == "hom_recessive" and any("Mendelian conflict: father" in f for f in flags)


def test_classify_x_mt_and_singletons():
    ok = _p("hom_ref", True)
    assert V._classify("x_nonpar", "male", "hemi", _p("het"), ok)[:2] == ("x_hemizygous", "maternal")
    assert V._classify("x_nonpar", "male", "hemi", ok, ok)[:2] == ("x_hemizygous", "de_novo")
    cls, _, flags = V._classify("x_nonpar", "male", "het", None, None)
    assert cls == "x_hemizygous" and "heterozygous call on male X" in flags[0]
    cls, _, flags = V._classify("x_nonpar", "female", "hemi", None, None)
    assert cls == "hom" and "female" in flags[0]
    assert V._classify("x_par", "male", "het", None, None)[0] == "het"  # PAR behaves as autosome
    assert V._classify("x_nonpar", None, "hemi", None, None)[0] == "x_hemizygous"
    assert V._classify("mt", None, "hom_alt", None, None)[0] == "mitochondrial"
    assert V._classify("auto", None, "het", None, None)[0] == "het"
    assert V._classify("auto", None, "hom_alt", None, None)[0] == "hom"
    cls, _, flags = V._classify("auto", None, "het", _p("hom_ref", True), None)  # duo
    assert cls == "het" and "other parent was not tested" in flags[0]


def test_parent_adequacy_needs_depth_quality_and_no_alt_reads():
    fmt = ["GT", "AD", "DP", "GQ"]
    assert V._parent(V.parse_call(fmt, "0/0:30,0:30:90"), 1, 10, 20)["adequate_ref"]
    assert not V._parent(V.parse_call(fmt, "0/0:27,3:30:90"), 1, 10, 20)["adequate_ref"]  # 3 ALT reads
    assert not V._parent(V.parse_call(fmt, "0/0:4,0:4:12"), 1, 10, 20)["adequate_ref"]
    assert not V._parent(V.parse_call(["GT"], "0/0"), 1, 10, 20)["adequate_ref"]  # no depth at all
    assert V._parent(V.parse_call(fmt, "0:20,0:20:60"), 1, 10, 20)["adequate_ref"]  # haploid father on X


def test_ctype_pars():
    assert V._ctype("X", 2_000_000, "GRCh38") == "x_par" and V._ctype("X", 2_600_000, "GRCh37") == "x_par" and V._ctype("X", 2_700_000, "GRCh37") == "x_nonpar"
    assert V._ctype("X", 31_178_721, "GRCh38") == "x_nonpar"
    assert V._ctype("X", 155_800_000, "GRCh38") == "x_par" and V._ctype("X", 155_800_000, "GRCh37") == "x_nonpar"


def test_moi_fit():
    assert V.moi_fit("de_novo", ["AD"], "auto") == "fits"
    assert V.moi_fit("de_novo", ["AR"], "auto") == "against"
    assert V.moi_fit("hom_recessive", ["AR"], "auto") == "fits"
    assert V.moi_fit("comphet", ["AD"], "auto") == "against"
    assert V.moi_fit("inherited_het", ["XL", "XLR"], "x_nonpar") == "against"  # female carrier of DMD
    assert V.moi_fit("x_hemizygous", ["XLR"], "x_nonpar") == "fits"
    assert V.moi_fit("het", [], "auto") == "unknown"


# ---------------------------------------------------- annotation & scores

def _vep_record(inp):
    for r in json.loads((FIX / "vep_trio.json").read_text())["records"]:
        if r["input"] == inp:
            return r
    raise KeyError(inp)


def test_annotate_vep_captured_scn1a():
    a = V.annotate_vep(_vep_record("2 166042334 . G A . . ."))
    assert a["gene"] == "SCN1A" and a["mane"] == "NM_001165963.4"
    assert a["hgvsc"].endswith(":c.2134C>T") and a["hgvsp"].endswith(":p.Arg712Ter")
    assert a["impact"] == "HIGH" and a["af"]["filter_af"] is None and "pathogenic" in a["clinvar"]


def test_af_summary_grpmax_excludes_bottlenecked():
    g551d = V.annotate_vep(_vep_record("7 117587806 . G A . . ."))["af"]
    assert g551d["exome"] == pytest.approx(0.0004041) and g551d["grpmax_group"] == "gnomadg_nfe"
    assert g551d["filter_af"] == pytest.approx(0.0005734)
    n1303k = V.annotate_vep(_vep_record("7 117652877 . C G . . ."))["af"]
    # asj 0.000865 and remaining 0.000945 are larger but excluded, as in gnomAD's grpmax
    assert n1303k["grpmax_group"] == "gnomadg_nfe" and n1303k["filter_af"] == pytest.approx(0.0002354)
    zero = V.annotate_vep(_vep_record("2 166053034 . C A . . ."))["af"]
    assert zero["filter_af"] == 0  # present in gnomAD with zero count, not absent
    assert V.af_summary({"af": 0.2, "eas": 0.3})["source"].startswith("1000 Genomes")


def test_predictors_never_outrank_lof():
    a = V.annotate_vep(_vep_record("7 117587806 . G A . . ."))
    score, support = V.predictor_support(a, "GRCh38")
    assert score == V.PREDICTOR_CAP and len(support) == 3
    assert score < V.SEVERITY["HIGH"]
    assert V.predictor_support({"spliceai_max": 0.25}, "GRCh38")[0] == 0.65


def test_region_index_overlaps():
    idx = V.RegionIndex({"A": {"symbol": "A", "chrom": "2", "start": 100, "end": 5000},
                         "B": {"symbol": "B", "chrom": "2", "start": 1000, "end": 1200},
                         "C": {"symbol": "C", "chrom": "3", "start": 10, "end": 20}}, pad=50)
    assert sorted(idx.hits("2", 1100, 1100)) == ["A", "B"]
    assert idx.hits("2", 5049, 5049) == ["A"] and idx.hits("2", 5051, 5051) == []
    assert idx.hits("3", 70, 70) == ["C"] and idx.hits("X", 10, 10) == []


def test_gene_regions_reports_unknown_symbols(fake_ensembl):
    got = V.gene_regions(["SCN1A", "NOTAGENE1", "cftr"], "GRCh38")
    assert set(got.result["regions"]) == {"SCN1A", "cftr"}
    assert got.result["regions"]["SCN1A"] == {"symbol": "SCN1A", "chrom": "2", "start": 165984641, "end": 166182806,
                                              "id": "ENSG00000144285", "biotype": "protein_coding"}
    assert got.result["missing"] == ["NOTAGENE1"] and "NOTAGENE1" in got.warnings[0]
    assert len(fake_ensembl["lookup"]) == 1  # one batched POST


# ---------------------------------------------------------------- triage

def test_triage_trio_genes_ranks_scn1a_de_novo_first(fake_ensembl, tmp_path):
    out = tmp_path / "t.tsv"
    got = V.triage(GZ, "P", mother="M", father="F", sex="female", genes=_genes(), out=str(out))
    res = got.result
    assert res["candidates"][0]["gene"] == "SCN1A" and res["candidates"][0]["class"] == "de_novo"
    assert res["candidates"][0]["hgvsp"] == "p.Arg712Ter"
    c = res["counts"]
    assert (c["records"], c["alleles"], c["filter_pass"], c["proband_carries_alt"], c["quality_pass"]) == (13, 14, 13, 12, 11)
    assert (c["in_gene_regions"], c["annotated"], c["af_pass_recessive"], c["comphet_genes"], c["af_pass_dominant"]) == (9, 9, 6, 1, 6)
    g = _by_gene(res)
    assert g["ALDH7A1"]["class"] == "hom_recessive" and g["ALDH7A1"]["moi_fit"] == "fits"
    cftr = [x for x in res["candidates"] if x["gene"] == "CFTR"]
    assert {x["class"] for x in cftr} == {"comphet"} and {x["origin"] for x in cftr} == {"maternal", "paternal"}
    assert cftr[0]["partners"] == [cftr[1]["variant"]]
    assert g["NGLY1"]["class"] == "possible_de_novo"
    assert g["DMD"]["class"] == "inherited_het" and g["DMD"]["moi_fit"] == "against"
    variants = {x["variant"] for x in res["candidates"]}
    # low-depth and LowQual SCN1A stops, common SNPs and out-of-list TTN never reach the list
    for gone in ("2-166036097-A-C", "2-166036116-C-A", "7-117559479-G-A", "2-166053034-C-T", "2-178560163-T-A"):
        assert gone not in variants
    assert res["thresholds"]["max_af_dominant"] == 0.0001 and res["phenotype"] is None
    assert set(res["weights"]) == {"variant", "inheritance"}
    comp = res["candidates"][0]["components"]
    assert comp["phenotype"] is None and comp["variant"] == 1.0 and comp["inheritance"] == 1.0
    # only candidate variants (no sample data) were sent, in one VEP batch
    assert len(fake_ensembl["vep"]) == 1 and len(fake_ensembl["vep"][0]) == 9
    assert all(len(v) == 4 for v in fake_ensembl["vep"][0])
    lines = out.read_text().splitlines()
    assert lines[0].startswith("## zebra vcf triage") and "dominant classes" in lines[1]
    assert lines[2].split("\t")[:3] == ["rank", "score", "phenotype_score"]
    assert len(lines) == 3 + 6 and lines[3].split("\t")[8] == "SCN1A"
    assert "SCN1A" in got.text.splitlines()[4]


def test_triage_plain_chr_vcf_gives_same_calls(fake_ensembl):
    a = V.triage(GZ, "P", mother="M", father="F", sex="female", genes=_genes()).result
    b = V.triage(PLAIN, "P", mother="M", father="F", sex="female", genes=_genes()).result
    assert [(x["variant"], x["class"], x["score"]) for x in a["candidates"]] == \
        [(x["variant"], x["class"], x["score"]) for x in b["candidates"]]


def test_triage_male_sibling_hemizygous(fake_ensembl):
    res = V.triage(GZ, "S", mother="M", father="F", sex="male", genes=_genes()).result
    top = res["candidates"][0]
    assert top["gene"] == "DMD" and top["class"] == "x_hemizygous" and top["origin"] == "maternal"
    assert top["moi_fit"] == "fits"
    assert "SCN1A" not in _by_gene(res)  # S does not carry the de novo


def test_triage_singleton_comphet_unphased_and_no_de_novo(fake_ensembl):
    res = V.triage(GZ, "P", genes=_genes()).result
    assert res["mode"] == "singleton" and any("de novo" in n for n in res["notes"])
    classes = {x["variant"]: x["class"] for x in res["candidates"]}
    assert classes["2-166042334-G-A"] == "het"
    assert classes["7-117587806-G-A"] == classes["7-117652877-C-G"] == "comphet_unphased"


def test_triage_duo_infers_trans(fake_ensembl):
    res = V.triage(GZ, "P", mother="M", genes=_genes()).result
    cftr = [x for x in res["candidates"] if x["gene"] == "CFTR"]
    assert len(cftr) == 2 and {x["class"] for x in cftr} == {"comphet"}
    assert any("in trans inferred from one parent" in f for f in cftr[0]["flags"])


def test_triage_budget_prioritises_de_novo_and_hom(fake_ensembl):
    got = V.triage(GZ, "P", mother="M", father="F", sex="female", max_annotate=2)
    res = got.result
    assert res["counts"]["sent_to_vep"] == 2 and res["counts"]["within_budget"] == 2
    assert res["restriction"]["mode"].startswith("prioritised")
    assert {x["class"] for x in res["candidates"]} <= {"de_novo", "hom_recessive"}
    assert any("NOT annotated" in w for w in got.warnings)


def test_triage_vep_failure_keeps_variants_unscored(fake_ensembl, monkeypatch):
    def broken(variants, assembly="GRCh38"):
        raise SourceError("Ensembl VEP", "https://rest.ensembl.org/vep/human/region", 503, "No server is available")

    monkeypatch.setattr(ensembl, "vep_batch", broken)
    got = V.triage(GZ, "P", mother="M", father="F", sex="female", genes=_genes())
    assert any(w.startswith("Ensembl VEP batch 1/1") and "503" in w for w in got.warnings)
    assert any("could not be annotated" in w for w in got.warnings)
    assert got.result["total_candidates"] == 9 and all(c["score"] is None for c in got.result["candidates"])


def test_triage_uses_vcf_info_annotations(fake_ensembl, tmp_path):
    vcf = tmp_path / "ann.vcf"
    vcf.write_text("\n".join([
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=248956422>",
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|gnomADe_AF">',
        '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="gnomAD AF">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP",
        "chr1\t1000\t.\tA\tG\t50\tPASS\tCSQ=G|missense_variant|MODERATE|GENEA|0.00001;gnomAD_AF=0.00001\tGT:DP:GQ\t1/1:30:90",
        "chr1\t2000\t.\tC\tT\t50\tPASS\tCSQ=T|intron_variant|MODIFIER|GENEA|\tGT:DP:GQ\t1/1:30:90",
        "chr1\t3000\t.\tG\tA\t50\tPASS\tCSQ=A|stop_gained|HIGH|GENEB|0.2;gnomAD_AF=0.2\tGT:DP:GQ\t0/1:30:90",
        "chr1\t4000\t.\tGT\tG\t50\tPASS\tCSQ=-|frameshift_variant|HIGH|GENEC|\tGT:DP:GQ\t1/1:30:90",
    ]) + "\n")
    got = V.triage(str(vcf), "P")
    r = got.result
    assert r["restriction"]["mode"] == "VCF INFO annotations"
    pre = r["restriction"]["info_prefilter"]
    assert pre["dropped_af_gt_max"] == 1 and pre["dropped_modifier_only"] == 1
    assert r["counts"]["info_prefilter_pass"] == 2
    # VEP (fixture) does not know these positions: the VCF's own annotation is used, and said so
    genes = {c["gene"]: c for c in r["candidates"]}
    assert set(genes) == {"GENEA", "GENEC"}
    assert genes["GENEC"]["impact"] == "HIGH" and genes["GENEC"]["annotation"].startswith("VCF INFO")


def test_triage_usage_errors(fake_ensembl, tmp_path):
    with pytest.raises(UsageError, match="not a sample"):
        V.triage(GZ, "X1")
    with pytest.raises(UsageError, match="contradicts"):
        V.triage(GZ, "P", assembly="GRCh37", genes=_genes())
    with pytest.raises(UsageError, match="must be different"):
        V.triage(GZ, "P", mother="P")
    bare = tmp_path / "bare.vcf"
    bare.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP\n"
                    "1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:90\n")
    with pytest.raises(UsageError, match="build unknown"):
        V.triage(str(bare), "P")
    with pytest.raises(UsageError, match="phenotypes"):
        V.triage(GZ, "P", hpo_genes=True)


@pytest.mark.skipif(not (REAL_HPO / "hp.json").exists(), reason="local HPO release not fetched (zebra hpo fetch)")
def test_triage_hpo_genes_from_case(fake_ensembl, monkeypatch, tmp_path):
    from zebra import case as case_mod

    monkeypatch.setenv("ZEBRA_HPO_DIR", str(REAL_HPO))
    case_dir = tmp_path / "dravet"
    case_mod.init(str(case_dir), title="Dravet-like")
    for hpo, label in (("HP:0002373", "Febrile seizure"), ("HP:0002133", "Status epilepticus"),
                       ("HP:0001263", "Global developmental delay"), ("HP:0007359", "Focal-onset seizure"),
                       ("HP:0002123", "Generalized myoclonic seizure")):
        case_mod.add_phenotype(str(case_dir), hpo, label)
    got = V.triage(GZ, "P", mother="M", father="F", sex="female", hpo_genes=True, case_dir=str(case_dir),
                   today="2026-10-05")
    res = got.result
    assert res["restriction"]["from_hpo"] >= 150 and "SCN1A" in fake_ensembl["lookup"][0]
    top = res["candidates"][0]
    assert top["gene"] == "SCN1A" and top["components"]["phenotype"] > 0.9
    assert "Dravet" in (top["phenotype_via"] or "")
    assert res["phenotype"]["origin"].startswith("case ") and len(res["phenotype"]["present"]) == 5
    assert "CFTR" not in _by_gene(res)  # not a gene for this phenotype: outside the restriction
    tsv = case_dir / "reports" / "triage-2026-10-05.tsv"
    assert tsv.exists() and res["tsv"] == str(tsv.resolve())
    assert any(s["db"] == "HPO annotations (local)" for s in got.sources)


def test_cli_inspect_json(capsys):
    from zebra.cli import main

    assert main(["--json", "vcf", "inspect", GZ]) == 0
    env = json.loads(capsys.readouterr().out)
    assert env["ok"] and env["command"] == "vcf inspect" and env["result"]["build"]["guess"] == "GRCh38"


# ------------------------------------------------------------------ live

@pytest.mark.live
def test_live_triage_end_to_end(tmp_path):
    os.environ.pop("ZEBRA_NO_CACHE", None)
    got = V.triage(GZ, "P", mother="M", father="F", sex="female", genes=_genes(), out=str(tmp_path / "live.tsv"))
    top = got.result["candidates"][0]
    assert top["gene"] == "SCN1A" and top["class"] == "de_novo" and top["mane"] == "NM_001165963.4"
    assert top["hgvsc"] == "c.2134C>T" and top["impact"] == "HIGH"
    assert any(s["db"] == "Ensembl VEP" and s["url"].startswith("https://rest.ensembl.org") for s in got.sources)
    assert any(s["db"] == "Ensembl lookup" for s in got.sources)
    assert got.result["comphet_genes"] == ["CFTR"]


def test_not_a_vcf_and_missing_file(tmp_path):
    junk = tmp_path / "notes.txt"
    junk.write_text("patient notes\nnot a vcf\n")
    with pytest.raises(UsageError, match="no #CHROM"):
        V.inspect(str(junk))
    with pytest.raises(UsageError, match="no such file"):
        V.inspect(str(tmp_path / "absent.vcf"))


def test_genes_path_and_no_tsv_outside_a_case(fake_ensembl, tmp_path):
    got = V.triage(GZ, "P", mother="M", father="F", sex="female", genes=str(FIX / "genes.txt"),
                   case_dir=str(tmp_path / "not-a-case"))
    assert got.result["tsv"] is None and any("no TSV written" in w for w in got.warnings)
    assert not (tmp_path / "not-a-case").exists()
    assert got.result["candidates"][0]["gene"] == "SCN1A"


# ===================================================================
# Regression tests for the defects found in review (docs/ROADMAP.md).
# Each is named after its backlog id.
# ===================================================================

def _info_vcf(tmp_path, info_lines, data_lines, name="info.vcf"):
    """A tiny VCF with the INFO header lines and data lines given."""
    p = tmp_path / name
    p.write_text("\n".join([
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr7,length=159345973>",
        *info_lines,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP",
        *data_lines,
    ]) + "\n")
    return str(p)


# ------------------------------------------------- S3: gnomAD AF parsing

def test_s3_frequency_parsing(tmp_path, fake_ensembl):
    """Allele counts must never be read as frequencies ("af" is a substring of "afr")."""
    import io

    header_text = "\n".join([
        "##fileformat=VCFv4.2",
        '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="x">',
        '##INFO=<ID=gnomAD_AF_afr,Number=A,Type=Float,Description="x">',
        '##INFO=<ID=gnomADe_AF,Number=A,Type=Float,Description="x">',
        '##INFO=<ID=gnomAD_faf95_grpmax,Number=A,Type=Float,Description="x">',
        '##INFO=<ID=gnomAD_AN_afr,Number=1,Type=Integer,Description="x">',
        '##INFO=<ID=gnomAD_AC_afr,Number=A,Type=Integer,Description="x">',
        '##INFO=<ID=gnomAD_nhomalt_afr,Number=A,Type=Integer,Description="x">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="this call set">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]) + "\n"
    h, _ = V.read_header(io.StringIO(header_text))
    assert h.gnomad_info_keys() == ["gnomAD_AF", "gnomAD_AF_afr", "gnomADe_AF", "gnomAD_faf95_grpmax"]
    for counted in ("gnomAD_AN_afr", "gnomAD_AC_afr", "gnomAD_nhomalt_afr", "AF"):
        assert counted not in h.gnomad_info_keys()

    # end to end: a real nonsense variant with AN_afr = 41000 in the INFO must survive
    vcf = _info_vcf(
        tmp_path,
        ['##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="x">',
         '##INFO=<ID=gnomAD_AN_afr,Number=1,Type=Integer,Description="x">',
         '##INFO=<ID=gnomAD_AC_afr,Number=A,Type=Integer,Description="x">',
         '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
         'Format: Allele|Consequence|IMPACT|SYMBOL">'],
        ["chr7\t117587806\t.\tG\tA\t500\tPASS\tCSQ=A|missense_variant|MODERATE|CFTR;gnomAD_AF=0;"
         "gnomAD_AN_afr=41000;gnomAD_AC_afr=0\tGT:AD:DP:GQ\t1/1:0,30:30:90"],
    )
    res = V.triage(vcf, "P").result
    assert res["restriction"]["info_prefilter"]["dropped_af_gt_max"] == 0
    assert res["counts"]["info_prefilter_pass"] == 1
    assert res["total_candidates"] == 1 and res["candidates"][0]["gene"] == "CFTR"


def test_s3_an_sized_value_in_an_af_field_is_refused(tmp_path, fake_ensembl):
    """A value above 1 in an AF-named field is a count: it must not drop the variant as 'too common'."""
    vcf = _info_vcf(
        tmp_path,
        ['##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="x">'],
        ["chr7\t117587806\t.\tG\tA\t500\tPASS\tgnomAD_AF=41000\tGT:AD:DP:GQ\t1/1:0,30:30:90"],
    )
    got = V.triage(vcf, "P")
    assert got.result["total_candidates"] == 1
    assert any("gnomAD_AF" in w and "above 1" in w for w in got.warnings)


# ------------------------------- E11: common pathogenic recessive alleles

def test_e11_clinvar_pathogenic_allele_survives_the_af_filter(fake_ensembl):
    """CFTR F508del (gnomAD grpmax 0.0149) must not be deleted by the 1 % default."""
    got = V.triage(FAM2, "P2", mother="M2", father="F2", genes=_genes())
    res = got.result
    f508 = next(c for c in res["candidates"] if c["variant"] == "7-117559590-ATCT-A")
    assert f508["gnomad"]["filter_af"] == pytest.approx(0.01498)
    assert f508["af_exempt"].startswith("ClinVar")
    assert any("ClinVar reports pathogenic" in f for f in f508["flags"])
    assert any("kept because ClinVar reports them pathogenic" in w for w in got.warnings)
    assert res["counts"]["af_pass_recessive"] == 3


def test_e11_af_removed_variants_are_named_in_a_warning(fake_ensembl):
    """A variant removed on frequency alone is reported, not dropped in silence."""
    got = V.triage(GZ, "P", mother="M", father="F", sex="female", genes=_genes())
    dropped = [w for w in got.warnings if "removed by gnomAD AF" in w]
    assert dropped and "7-117559479-G-A" in dropped[0]  # CFTR V470M, a common benign SNP


# ---------------------------------------- F16: ClinVar significance, indels

def test_f16_clinvar_significance_on_indels(fake_ensembl):
    """VEP writes an empty allele in clin_sig_allele for a deletion; "" and "-" are the same allele."""
    rec = next(r for r in json.loads((FIX / "vep_family2.json").read_text())["records"]
               if r["input"].startswith("7 117559590"))
    assert rec["allele_string"] == "TCT/-"  # the VEP spelling that used not to match
    a = V.annotate_vep(rec)
    assert "pathogenic" in a["clinvar"] and "likely_pathogenic" in a["clinvar"]
    assert V.clinvar_pathogenic(a) == ["pathogenic", "likely_pathogenic"]
    assert V._allele_keys("-") == {"-", ""} and V._allele_keys("") == {"-", ""}
    assert V._allele_keys("A") == {"A"}
    got = V.triage(FAM2, "P2", mother="M2", father="F2", genes=_genes())
    f508 = next(c for c in got.result["candidates"] if c["variant"] == "7-117559590-ATCT-A")
    assert "pathogenic" in (f508["clinvar"] or [])


# ------------------------------------------------- F21: truncated gzip

def test_f21_truncated_gzip_gives_a_usage_error_not_a_traceback(tmp_path, capsys):
    """A half-downloaded .vcf.gz must produce a JSON error envelope."""
    whole = Path(FAM2_GZ).read_bytes()
    trunc = tmp_path / "trunc.vcf.gz"
    trunc.write_bytes(whole[: len(whole) // 2])
    with pytest.raises(UsageError, match="truncated or corrupt"):
        V.inspect(str(trunc))

    from zebra.cli import main

    code = main(["--json", "vcf", "inspect", str(trunc)])
    env = json.loads(capsys.readouterr().out)
    assert code == 2 and env["ok"] is False and env["error"]["type"] == "UsageError"
    assert "truncated or corrupt" in env["error"]["message"]


def test_f21_garbage_gzip_body_is_also_an_envelope(tmp_path, capsys):
    p = tmp_path / "junk.vcf.gz"
    p.write_bytes(b"\x1f\x8b\x08\x00" + b"\x00" * 60)
    from zebra.cli import main

    assert main(["--json", "vcf", "inspect", str(p)]) == 2
    env = json.loads(capsys.readouterr().out)
    assert env["ok"] is False and env["error"]["type"] == "UsageError"


# ----------------------------------------------- F22: malformed data lines

def test_f22_malformed_lines_are_counted_and_reported(fake_ensembl):
    """bad.vcf holds 7 data lines, 5 of which cannot be read: none may vanish silently."""
    got = V.inspect(BAD)
    res = got.result
    assert res["variants"] == 2  # the space-separated, 3e2, 7-column, POS 0 and empty-ALT lines are not variants
    assert res["malformed_lines"] == 5
    note = next(n for n in res["notes"] if "could not be read" in n)
    assert "5 of 7 data line(s)" in note
    assert "space-separated" in note and "not an integer" in note.replace("POS '3e2' is", "POS is")
    assert "line " in note

    rec, why = V.parse_line_why("chr1\t0\t.\tA\tG\t.\t.\t.", 0)
    assert rec is None and "not a 1-based position" in why
    rec, why = V.parse_line_why("chr1\t10\t.\tA\t\t.\t.\t.", 0)
    assert rec is None and "empty REF or ALT" in why


def test_f22_triage_refuses_a_file_that_is_mostly_unreadable(tmp_path, fake_ensembl):
    lines = ["chr7 117587806 . G A 500 PASS . GT:DP:GQ 0/1:30:90"] * 30
    lines.append("chr7\t117587806\t.\tG\tA\t500\tPASS\t.\tGT:DP:GQ\t0/1:30:90")
    vcf = _info_vcf(tmp_path, [], lines, name="mostly_bad.vcf")
    with pytest.raises(UsageError, match="of the data lines; triage refuses"):
        V.triage(vcf, "P")


def test_f22_triage_warns_about_a_few_bad_lines(tmp_path, fake_ensembl):
    lines = ["chr7 117587806 . G A 500 PASS . GT:DP:GQ 0/1:30:90"]
    lines += ["chr7\t%d\t.\tG\tA\t500\tPASS\t.\tGT:DP:GQ\t0/1:30:90" % (117587806 + i) for i in range(20)]
    vcf = _info_vcf(tmp_path, [], lines, name="one_bad.vcf")
    got = V.triage(vcf, "P")
    assert got.result["counts"]["malformed_lines"] == 1
    assert any("could not be read" in w for w in got.warnings)


# ------------------------------------------------------ F31: inheritance

def test_f31_comphet_partner_is_not_held_to_the_dominant_cutoff(fake_ensembl):
    """A de novo plus an inherited allele in one recessive gene is a comp-het candidate, not two dominant ones."""
    got = V.triage(FAM2, "P2", mother="M2", father="F2", genes=_genes())
    res = got.result
    cftr = [c for c in res["candidates"] if c["gene"] == "CFTR"]
    assert len(cftr) == 2
    # the de novo keeps its class (the strongest label in triage); its partner,
    # an ordinary inherited het, becomes the comp-het
    assert {c["class"] for c in cftr} == {"de_novo", "comphet_unphased"}
    assert {c["comphet"] for c in cftr} == {"unphased"}
    assert {c["origin"] for c in cftr} == {"de_novo", "maternal"}
    for c in cftr:
        assert c["partners"] and any("the recessive AF threshold is used" in f for f in c["flags"])
        assert c["af_exempt"]  # both are held to the recessive threshold, not the dominant one
    # the de novo partner's AF (5.7e-4) is above the dominant cut-off and it still survives
    g551d = next(c for c in cftr if c["variant"] == "7-117587806-G-A")
    assert g551d["gnomad"]["filter_af"] > res["thresholds"]["max_af_dominant"]
    assert res["counts"]["comphet_genes"] == 1 and res["counts"]["af_pass_dominant"] == 3


def test_f31_missing_parental_genotype_is_unknown_not_absent(fake_ensembl, tmp_path):
    """A `./.` parental call must not be read as "not carried", and must not prove trans."""
    cls, _, flags = V._classify("auto", None, "het", _p("missing"), None)
    assert cls == "het" and any("carrier status unknown" in f for f in flags)
    cls, _, flags = V._classify("auto", None, "het", _p("hom_ref", False, ["DP 4 < 10"]), None)
    assert cls == "het" and any("not usable" in f and "carrier status unknown" in f for f in flags)
    assert V._parent_ref_ok(_p("hom_ref", True), None) is True
    assert V._parent_ref_ok(_p("missing"), None) is False
    assert V._parent_ref_ok(_p("het"), None) is True  # a carrier is known, not unknown

    # duo: the mother is ./. at one CFTR allele and 0/1 at the other
    vcf = tmp_path / "duo_missing.vcf"
    vcf.write_text("\n".join([
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr7,length=159345973>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tM\tP",
        "chr7\t117587806\t.\tG\tA\t500\tPASS\t.\tGT:AD:DP:GQ\t./.:0,0:0:0\t0/1:15,17:32:99",
        "chr7\t117652877\t.\tC\tG\t500\tPASS\t.\tGT:AD:DP:GQ\t0/1:14,14:28:99\t0/1:16,14:30:99",
    ]) + "\n")
    res = V.triage(str(vcf), "P", mother="M", genes=_genes()).result
    cftr = [c for c in res["candidates"] if c["gene"] == "CFTR"]
    assert len(cftr) == 2 and {c["class"] for c in cftr} == {"comphet_unphased"}
    assert any("carrier status unknown, so trans was NOT inferred" in f for c in cftr for f in c["flags"])
    assert not any("in trans inferred from one parent" in f for c in cftr for f in c["flags"])


def test_f31_diploid_male_x_is_inferred_from_the_genotypes(fake_ensembl):
    """A `1/1` call on male X is hemizygous, not a homozygote with a Mendelian conflict."""
    inferred = V.infer_sex(FAM2, "P2", "GRCh38")
    assert inferred["sex"] == "male" and "called Y non-PAR sites carry an ALT allele" in inferred["basis"]
    assert inferred["counts"]["y_alt"] == 2

    got = V.triage(FAM2, "P2", mother="M2", father="F2", genes=_genes())
    res = got.result
    assert res["sex"] == "male" and res["sex_inferred"]["sex"] == "male"
    dmd = next(c for c in res["candidates"] if c["gene"] == "DMD")
    assert dmd["class"] == "x_hemizygous" and dmd["origin"] == "maternal"
    assert not any("Mendelian conflict" in f for f in dmd["flags"])
    assert any("inferred as male" in w for w in got.warnings)
    # --sex still wins, and the explicit call agrees
    explicit = V.triage(FAM2, "P2", mother="M2", father="F2", sex="male", genes=_genes()).result
    assert explicit["sex_inferred"] is None
    assert next(c for c in explicit["candidates"] if c["gene"] == "DMD")["class"] == "x_hemizygous"


def test_f31_unknown_sex_never_blames_the_father_on_x(fake_ensembl, tmp_path):
    """When the sex cannot be inferred, an X hom-alt call says so instead of pointing at UPD or a sample swap."""
    cls, _, flags = V._classify("x_nonpar", None, "hom_alt", _p("het"), _p("hom_ref", True))
    assert cls == "hom_recessive"
    assert any("sex unknown" in f and "hemizygous" in f for f in flags)
    assert not any("Mendelian conflict" in f for f in flags)


def test_f31_phase_set_resolves_cis_and_trans(fake_ensembl, tmp_path):
    """FORMAT/PS and a 1/2 call settle the phase without the parents."""
    vcf = tmp_path / "phased.vcf"
    vcf.write_text("\n".join([
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr7,length=159345973>",
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP",
        "chr7\t117587806\t.\tG\tA\t500\tPASS\t.\tGT:AD:DP:GQ:PS\t0|1:15,17:32:99:117587806",
        "chr7\t117652877\t.\tC\tG\t500\tPASS\t.\tGT:AD:DP:GQ:PS\t0|1:16,14:30:99:117587806",
    ]) + "\n")
    res = V.triage(str(vcf), "P", genes=_genes(), max_af_dominant=0.01).result
    cftr = [c for c in res["candidates"] if c["gene"] == "CFTR"]
    assert len(cftr) == 2
    assert {c["class"] for c in cftr} == {"het"}  # same haplotype: not a compound heterozygote
    assert all(any("in cis" in f for f in c["flags"]) for c in cftr)
    assert res["counts"]["comphet_genes"] == 0

    trans = tmp_path / "trans.vcf"
    trans.write_text("\n".join([
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr7,length=159345973>",
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP",
        "chr7\t117587806\t.\tG\tA\t500\tPASS\t.\tGT:AD:DP:GQ:PS\t0|1:15,17:32:99:1",
        "chr7\t117652877\t.\tC\tG\t500\tPASS\t.\tGT:AD:DP:GQ:PS\t1|0:16,14:30:99:1",
    ]) + "\n")
    res = V.triage(str(trans), "P", genes=_genes()).result
    cftr = [c for c in res["candidates"] if c["gene"] == "CFTR"]
    assert {c["class"] for c in cftr} == {"comphet"}
    assert all(any("read-backed phase, in trans" in f for f in c["flags"]) for c in cftr)


def test_f31_two_alts_at_one_site_are_in_trans(tmp_path, fake_ensembl):
    """A 1/2 call puts the two ALT alleles on different chromosomes by definition."""
    a = {"site": "7:100:G", "variant": "7-100-G-A", "pos": 100, "trans_at_site": True, "phase": {"gt": "1/2"}}
    b = {"site": "7:100:G", "variant": "7-100-G-C", "pos": 100, "trans_at_site": True, "phase": {"gt": "1/2"}}
    rel, why = V._phase_relation(a, b)
    assert rel == "trans" and "necessarily in trans" in why
    cis = {"site": "7:200:G", "variant": "7-200-G-A", "pos": 200, "phase": {"ps": 7, "hap": 1}}
    cis2 = {"site": "7:300:G", "variant": "7-300-G-A", "pos": 300, "phase": {"ps": 7, "hap": 1}}
    assert V._phase_relation(cis, cis2)[0] == "cis"
    unknown = {"site": "7:400:G", "variant": "7-400-G-A", "pos": 400, "phase": {"ps": None, "hap": None}}
    assert V._phase_relation(cis, unknown)[0] is None
    # the same allele written twice is not a pair, and a phase set that spans
    # more than read-backed phasing can reach settles nothing
    same = dict(cis, site="7:200:G")
    assert V._phase_relation(cis, same)[0] is None
    far = {"site": "7:9000000:G", "variant": "7-9000000-G-A", "pos": 9_000_000, "phase": {"ps": 7, "hap": 1}}
    rel, why = V._phase_relation(cis, far)
    assert rel is None and "too far apart for read-backed phasing" in why

    vcf = _info_vcf(
        tmp_path,
        ['##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
         'Format: Allele|Consequence|IMPACT|SYMBOL">'],
        ["chr7\t900000\t.\tG\tA,C\t500\tPASS\tCSQ=A|missense_variant|MODERATE|GENEZ,C|stop_gained|HIGH|GENEZ"
         "\tGT:AD:DP:GQ\t1/2:0,15,14:29:99"],
        name="multiallelic.vcf",
    )
    res = V.triage(vcf, "P").result
    assert {c["class"] for c in res["candidates"]} == {"comphet"}
    assert all(any("necessarily in trans" in f for f in c["flags"]) for c in res["candidates"])


def test_f31_dominant_cutoff_uses_a_bound_not_a_per_group_point_estimate(fake_ensembl, tmp_path):
    """One observation in a small gnomAD group must not delete a de novo candidate."""
    # a single allele in gnomAD v4 genomes EAS (AN ~ 5,208) is a 1.9e-4 point estimate
    summary = V.af_summary({"gnomadg": 1.3e-5, "gnomadg_eas": 1 / 5208})
    assert summary["filter_af"] == pytest.approx(1 / 5208)
    assert summary["bound_af"] == pytest.approx(1.3e-5)
    assert "point estimate" in summary["bound_basis"]
    # when a filtering AF is served, that is the bound
    with_faf = V.af_summary({"gnomadg": 1.3e-5, "gnomadg_eas": 1 / 5208, "gnomadg_faf95": 3e-6})
    assert with_faf["bound_af"] == pytest.approx(3e-6) and with_faf["faf95"] == pytest.approx(3e-6)
    assert "filtering AF" in with_faf["bound_basis"]

    rec = _vep_record("2 166042334 . G A . . .")  # SCN1A de novo, absent from gnomAD
    patched = json.loads(json.dumps(rec))
    for cv in patched["colocated_variants"]:
        cv["frequencies"] = {"A": {"gnomadg": 1.3e-5, "gnomadg_eas": 1 / 5208}}
    a = V.annotate_vep(patched)
    assert a["af"]["grpmax"] == pytest.approx(1 / 5208) and a["af"]["bound_af"] == pytest.approx(1.3e-5)

    records = {r["input"]: r for r in json.loads((FIX / "vep_trio.json").read_text())["records"]}
    records[patched["input"]] = patched  # only SCN1A's frequencies are changed

    def patched_batch(variants, assembly="GRCh38"):
        lines = [ensembl.vcf_line(*v) for v in variants]
        return Outcome([records[l] for l in lines if l in records],
                       sources=[source_record("Ensembl VEP", "fixture",
                                              url="https://rest.ensembl.org/vep/human/region")])

    import zebra.sources.ensembl as ens

    original = ens.vep_batch
    ens.vep_batch = patched_batch
    try:
        got = V.triage(GZ, "P", mother="M", father="F", sex="female", genes=_genes())
    finally:
        ens.vep_batch = original
    scn1a = next(c for c in got.result["candidates"] if c["variant"] == "2-166042334-G-A")
    assert scn1a["class"] == "de_novo"  # kept: the bound, not the per-group estimate, is compared
    assert any("is not evidence against pathogenicity" in f for f in scn1a["flags"])


# --------------------------------------------- C-P2-7: in-region memory cap

def test_c_p2_7_in_region_candidates_are_capped(fake_ensembl):
    """The per-class cap applies with --genes too, and says how many were skipped."""
    got = V.triage(GZ, "P", mother="M", father="F", sex="female", genes=_genes(), max_annotate=1)
    res = got.result
    assert res["counts"]["sent_to_vep"] == 1
    assert res["restriction"]["skipped_over_budget"]
    assert any("NOT annotated" in w and "narrow the gene list" in w for w in got.warnings)


# -------------------------------------------- P1d: the other report forms

@pytest.fixture
def fake_cnv_sources(monkeypatch):
    """Ensembl overlap/lookup and ClinGen dosage answered from captured responses."""
    from zebra import cnv as C
    from zebra.sources import clingen
    from zebra.sources import ensembl as ens

    overlap = json.loads((FIX / "overlap_cnv.json").read_text())
    dosage_tsv = json.loads((FIX / "dosage_cnv.json").read_text())["tsv"]
    dmd = json.loads((FIX / "lookup_dmd.json").read_text())["response"]
    calls = {"overlap": [], "dosage": [], "lookup": []}

    def get_json(url, source, params=None, **kw):
        calls["overlap"].append(url)
        return Response(url, 200, json.dumps(overlap["response"]), "2026-10-06T00:00:00+00:00", False)

    def request(url, source, **kw):
        """The ClinGen bulk TSV, fetched once and parsed per gene (its real header and rows)."""
        calls["dosage"].append(url)
        return Response(url, 200, dosage_tsv, "2026-10-06T00:00:00+00:00", False)

    def lookup_symbol(symbol, assembly="GRCh38", expand=False):
        calls["lookup"].append((symbol, expand))
        if symbol.upper() == "DMD":
            data = dmd if expand else {k: v for k, v in dmd.items() if k != "Transcript"}
            return Outcome(data, sources=[source_record("Ensembl lookup", symbol,
                                                        url="https://rest.ensembl.org/lookup/symbol")])
        return Outcome({"seq_region_name": "5", "start": 70925030, "end": 70953942, "strand": 1,
                        "id": "ENSG00000172062", "biotype": "protein_coding"},
                       sources=[source_record("Ensembl lookup", symbol,
                                              url="https://rest.ensembl.org/lookup/symbol")])

    monkeypatch.setattr(C, "get_json", get_json)
    monkeypatch.setattr(C, "request", request)
    monkeypatch.setattr(ens, "lookup_symbol", lookup_symbol)
    assert clingen.DOSAGE_TSV["GRCh38"]  # the real url the bulk fetch uses
    return calls


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
    assert exon == {"kind": "exon_cnv", "input": "DMD exon 45-50 deletion", "gene": "DMD", "transcript": None,
                    "first": 45, "last": 50, "cds_start": None, "cds_end": None, "cnv_type": "loss"}
    assert C.parse("SMN1 exon 7 deletion")["first"] == 7
    assert C.parse("MLPA: DMD exon 8 duplication".split(": ")[1])["cnv_type"] == "gain"

    hgvs = C.parse("NM_004006.3:c.6439-?_7309+?del")
    assert hgvs["kind"] == "exon_cnv" and (hgvs["cds_start"], hgvs["cds_end"]) == (6439, 7309)
    assert hgvs["transcript"] == "NM_004006.3" and hgvs["cnv_type"] == "loss"

    cn = C.parse("SMN1 exon 7 copy number 0")
    assert cn == {"kind": "copy_number", "input": "SMN1 exon 7 copy number 0", "gene": "SMN1", "exon": "7",
                  "copy_number": 0}
    assert C.parse("SMN2 copy number 2")["copy_number"] == 2
    assert C.parse("SMN1 0 copies")["copy_number"] == 0

    rep = C.parse("FMR1 CGG 230")
    assert rep == {"kind": "repeat_expansion", "input": "FMR1 CGG 230", "gene": "FMR1", "motif": "CGG",
                   "repeat_count": 230}
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
    assert len(fake_cnv_sources["dosage"]) == 1  # one bulk TSV, not one request per gene
    assert "ClinGen_gene_curation_list" in fake_cnv_sources["dosage"][0]
    inputs = r["acmg_cnv_inputs"]
    assert inputs["framework"].startswith("ACMG/ClinGen technical standard")
    assert inputs["section_1_variant_type"]["contains_protein_coding_genes"] is True
    hi = inputs["section_2_overlap_with_established_regions_or_genes"]["clingen_established_haploinsufficient_genes"]
    assert hi == [{"gene": "UBE3A", "score": "3", "description": "Sufficient evidence for dosage pathogenicity"}]
    assert inputs["section_3_gene_number"]["protein_coding_genes"] == 1
    assert inputs["section_4_case_level_evidence"] is None
    assert inputs["section_5_inheritance_and_family_history"] == {"reported_inheritance": "de_novo"}
    # the interpretation itself is never invented
    assert inputs["classification"] is None and "does not score or classify" in inputs["classification_note"]
    assert any(s["db"] == "Ensembl overlap" for s in got.sources)
    assert any(s["db"] == "ClinGen dosage sensitivity" for s in got.sources)
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
    assert rep["classification"] is None and rep["thresholds"] is None
    assert "gene-specific" in rep["thresholds_note"]
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


def test_e11_the_clinvar_exemption_stops_at_the_ba1_threshold(fake_ensembl, tmp_path):
    """A ClinVar assertion must not smuggle a >5 % allele past the frequency filter."""
    rec = json.loads(json.dumps(_vep_record("7 117587806 . G A . . .")))
    for cv in rec["colocated_variants"]:
        cv["frequencies"] = {"A": {"gnomade": 0.21, "gnomadg": 0.19}}
    a = V.annotate_vep(rec)
    assert V.clinvar_pathogenic(a) and a["af"]["filter_af"] == pytest.approx(0.21)

    records = {r["input"]: r for r in json.loads((FIX / "vep_trio.json").read_text())["records"]}
    records[rec["input"]] = rec

    def patched_batch(variants, assembly="GRCh38"):
        lines = [ensembl.vcf_line(*v) for v in variants]
        return Outcome([records[l] for l in lines if l in records],
                       sources=[source_record("Ensembl VEP", "fixture", url="https://rest.ensembl.org/vep")])

    import zebra.sources.ensembl as ens

    original = ens.vep_batch
    ens.vep_batch = patched_batch
    try:
        got = V.triage(GZ, "P", mother="M", father="F", sex="female", genes=_genes())
    finally:
        ens.vep_batch = original
    assert "7-117587806-G-A" not in {c["variant"] for c in got.result["candidates"]}
    assert any("above the ClinGen BA1 threshold" in w for w in got.warnings)


def test_f31_sex_inference_states_exactly_what_it_measured(fake_ensembl, tmp_path):
    """The basis must never assert a count that was not taken, and Y noise must not make a female male."""
    got = V.infer_sex(FAM2, "P2", "GRCh38")
    assert got["sex"] == "male"
    assert got["counts"]["scan_complete"] is True  # the whole file was read: no order-dependent answer
    assert "2 of 2 called Y non-PAR sites" in got["basis"]
    assert "X non-PAR: 1 ALT call(s)" in got["basis"]
    # no X or Y data at all: no guess, and the counts say why
    assert V.infer_sex(PLAIN, "M", "GRCh38")["sex"] is None
    missing = V.infer_sex(FAM2, "nobody", "GRCh38")
    assert missing["sex"] is None and "not a sample" in missing["basis"]

    def write(name, rows):
        p = tmp_path / name
        head = ["##fileformat=VCFv4.2", "##contig=<ID=chrX,length=156040895>",
                "##contig=<ID=chrY,length=57227415>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP"]
        p.write_text("\n".join(head + rows) + "\n")
        return str(p)

    # a female with two stray Y ALT calls among many called Y sites, X clearly het
    rows = [f"chrX\t{31178721 + i * 97}\t.\tG\tA\t500\tPASS\t.\tGT:DP:GQ\t0/1:30:99" for i in range(40)]
    rows += [f"chrY\t{2786042 + i * 53}\t.\tT\tG\t500\tPASS\t.\tGT:DP:GQ\t"
             + ("1/1:20:60" if i < 2 else "0/0:20:60") for i in range(40)]
    female = V.infer_sex(write("female_y_noise.vcf", rows), "P", "GRCh38")
    assert female["sex"] == "female", female["basis"]
    assert "2 of 40 called Y non-PAR sites" not in female["basis"]  # the male branch was not taken
    assert "no consistent ALT calls" in female["basis"]
    # the same genotypes with Y written before X give the same answer (no early stop)
    reordered = V.infer_sex(write("female_reordered.vcf", rows[40:] + rows[:40]), "P", "GRCh38")
    assert reordered["sex"] == "female" and reordered["counts"]["y_alt"] == female["counts"]["y_alt"]
    # a real male: most called Y sites carry the ALT allele
    male_rows = [f"chrX\t{31178721 + i * 97}\t.\tG\tA\t500\tPASS\t.\tGT:DP:GQ\t1/1:30:99" for i in range(40)]
    male_rows += [f"chrY\t{2786042 + i * 53}\t.\tT\tG\t500\tPASS\t.\tGT:DP:GQ\t1/1:20:60" for i in range(20)]
    male = V.infer_sex(write("male.vcf", male_rows), "P", "GRCh38")
    assert male["sex"] == "male" and "20 of 20 called Y non-PAR sites" in male["basis"]


def test_f31_the_parents_outrank_a_phase_block_that_disagrees(fake_ensembl, tmp_path):
    """A trio-proven trans pair must survive a phase set that calls it cis, with the conflict named."""
    vcf = tmp_path / "conflict.vcf"
    vcf.write_text("\n".join([
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr7,length=159345973>",
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tM\tF\tP",
        # maternal G551D and paternal N1303K, yet both phased onto haplotype 1 of one PS block
        "chr7\t117587806\t.\tG\tA\t500\tPASS\t.\tGT:AD:DP:GQ:PS\t0/1:16,15:31:99\t0/0:30,0:30:90"
        "\t0|1:15,17:32:99:117587806",
        "chr7\t117652877\t.\tC\tG\t500\tPASS\t.\tGT:AD:DP:GQ:PS\t0/0:29,0:29:87\t0/1:14,14:28:99"
        "\t0|1:16,14:30:99:117587806",
    ]) + "\n")
    res = V.triage(str(vcf), "P", mother="M", father="F", genes=_genes()).result
    cftr = [c for c in res["candidates"] if c["gene"] == "CFTR"]
    assert len(cftr) == 2 and {c["class"] for c in cftr} == {"comphet"}
    assert {c["origin"] for c in cftr} == {"maternal", "paternal"}
    assert all(any("the parents' genotypes put them in trans" in f for f in c["flags"]) for c in cftr)
    assert all(not any("not a compound heterozygote" in f for f in c["flags"]) for c in cftr)


def test_f31_a_female_hom_alt_x_call_still_flags_the_father(fake_ensembl):
    """She has her father's X: a hom-ref father IS a conflict. Only the male/unknown case is exempt."""
    cls, _, flags = V._classify("x_nonpar", "female", "hom_alt", _p("het"), _p("hom_ref", True))
    assert cls == "hom_recessive"
    assert any("Mendelian conflict: father is hom-ref" in f for f in flags)
    # a male proband never reaches this branch; an unknown sex says so instead of blaming a parent
    assert V._classify("x_nonpar", "male", "hom_alt", _p("het"), _p("hom_ref", True))[0] == "x_hemizygous"
    cls, _, flags = V._classify("x_nonpar", None, "hom_alt", _p("het"), _p("hom_ref", True))
    assert not any("Mendelian conflict" in f for f in flags)


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
    """One copy of X is normal in a male and a loss in a female: zebra must not pick one."""
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


# -------------- second-round review findings (same ids, deeper defects)

def _one_sample_vcf(tmp_path, rows, name, contigs=("chr7\t159345973",)):
    p = tmp_path / name
    head = ["##fileformat=VCFv4.2"] + [f"##contig=<ID={c.split()[0]},length={c.split()[1]}>" for c in contigs]
    head.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP")
    p.write_text("\n".join(head + rows) + "\n")
    return str(p)


def test_f31_one_allele_on_two_lines_is_not_a_compound_heterozygote(fake_ensembl, tmp_path):
    """A merged call set can write one deletion twice; it must not pair with itself."""
    vcf = _one_sample_vcf(tmp_path, [
        "chr7\t117559589\tcallerA\tCATCT\tCA\t500\tPASS\t.\tGT:AD:DP:GQ\t0/1:18,17:35:99",
        "chr7\t117559590\tcallerB\tATCT\tA\t500\tPASS\t.\tGT:AD:DP:GQ\t0/1:18,17:35:99",
    ], "dup.vcf")
    got = V.triage(vcf, "P", genes=_genes())
    res = got.result
    assert res["comphet_genes"] == [] and res["total_candidates"] == 1
    only = res["candidates"][0]
    assert only["variant"] == "7-117559590-ATCT-A" and only["class"] == "het"
    assert only["partners"] is None and only["comphet"] is None
    assert any("repeated an allele already counted" in w for w in got.warnings)
    assert any("written on more than one line" in f for f in only["flags"])
    # and the helper refuses to pair an allele with itself
    same = {"variant": "7-1-A-G", "site": "7:1:A", "pos": 1, "trans_at_site": True, "phase": {"gt": "1/2"}}
    assert V._phase_relation(same, dict(same))[0] is None


def test_f31_the_dominant_threshold_names_what_it_deletes(fake_ensembl, tmp_path):
    """A silent deletion at the dominant cut-off is what made triage answer 'no candidate'."""
    vcf = _one_sample_vcf(tmp_path, [
        "chr7\t117587806\t.\tG\tA\t500\tPASS\t.\tGT:AD:DP:GQ\t0/1:16,15:31:99",
    ], "dom.vcf")
    got = V.triage(vcf, "P", genes=_genes())
    assert got.result["total_candidates"] == 0
    named = [w for w in got.warnings if "the dominant cut-off" in w]
    assert named and "7-117587806-G-A CFTR" in named[0]
    assert "--max-af-dominant" in named[0]
    assert V.triage(vcf, "P", genes=_genes(), max_af_dominant=0.01).result["total_candidates"] == 1


def test_f31_a_cis_call_that_suppresses_a_comphet_says_so(fake_ensembl, tmp_path):
    """Suppressing a comp-het on phase alone must leave a trace the reader can act on."""
    vcf = _one_sample_vcf(tmp_path, [
        "chr7\t117587806\t.\tG\tA\t500\tPASS\t.\tGT:AD:DP:GQ:PS\t0|1:15,17:32:99:117587806",
        "chr7\t117652877\t.\tC\tG\t500\tPASS\t.\tGT:AD:DP:GQ:PS\t0|1:16,14:30:99:117587806",
    ], "cis.vcf")
    got = V.triage(vcf, "P", genes=_genes())
    assert got.result["comphet_genes"] == []
    assert any("are on one haplotype according to the phase set" in w for w in got.warnings)
    assert any("the dominant cut-off" in w for w in got.warnings)  # and what that then cost


def test_f31_a_phase_set_wider_than_read_backed_phasing_settles_nothing(fake_ensembl, tmp_path):
    """Population phasing written back as PS must not be read as read-backed phase."""
    assert V.PHASE_TRUSTED_SPAN == 500_000
    # DMD spans 2.2 Mb, so two alleles in it can sit further apart than reads reach
    vcf = _one_sample_vcf(tmp_path, [
        "chrX\t31178721\t.\tG\tA\t500\tPASS\t.\tGT:AD:DP:GQ:PS\t0|1:13,14:27:99:1",
        "chrX\t32000000\t.\tC\tG\t500\tPASS\t.\tGT:AD:DP:GQ:PS\t0|1:16,14:30:99:1",
    ], "wide_ps.vcf", contigs=("chrX\t156040895",))
    got = V.triage(vcf, "P", sex="female", genes=_genes(), max_af_dominant=0.01)
    dmd = [c for c in got.result["candidates"] if c["gene"] == "DMD"]
    assert len(dmd) == 2, [c["variant"] for c in got.result["candidates"]]
    assert {c["class"] for c in dmd} == {"comphet_unphased"}
    assert all(any("too far apart for read-backed phasing" in f for f in c["flags"]) for c in dmd)
    # within a read-backed distance the same phase set still settles cis
    close = _one_sample_vcf(tmp_path, [
        "chrX\t31178721\t.\tG\tA\t500\tPASS\t.\tGT:AD:DP:GQ:PS\t0|1:13,14:27:99:1",
        "chrX\t31180000\t.\tC\tG\t500\tPASS\t.\tGT:AD:DP:GQ:PS\t0|1:16,14:30:99:1",
    ], "close_ps.vcf", contigs=("chrX\t156040895",))
    assert V.triage(close, "P", sex="female", genes=_genes(), max_af_dominant=0.01).result["comphet_genes"] == []


def test_s3_a_csq_frequency_above_one_is_refused_too(fake_ensembl, tmp_path):
    """The CSQ path had no plausibility guard, so a count there deleted the variant."""
    vcf = _info_vcf(
        tmp_path,
        ['##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
         'Format: Allele|Consequence|IMPACT|SYMBOL|gnomADe_AF">'],
        ["chr7\t117587806\t.\tG\tA\t500\tPASS\tCSQ=A|missense_variant|MODERATE|CFTR|4321"
         "\tGT:AD:DP:GQ\t1/1:0,30:30:90"],
        name="csq_count.vcf",
    )
    got = V.triage(vcf, "P")
    assert got.result["total_candidates"] == 1
    assert any("gnomADe_AF" in w and "above 1" in w for w in got.warnings)


def test_s3_frequency_field_names_and_the_type_guard(fake_ensembl):
    import io

    lines = ["##fileformat=VCFv4.2"]
    for key, typ in (("gnomad_exomes_AF", "Float"), ("gnomAD_genomes_AF_nfe", "Float"),
                     ("gnomAD_AFR_AF", "Float"), ("AF_grpmax", "Float"),
                     ("AF_popmax", "Integer"), ("gnomAD_AC_afr", "Integer")):
        lines.append(f'##INFO=<ID={key},Number=A,Type={typ},Description="x">')
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    h, _ = V.read_header(io.StringIO("\n".join(lines) + "\n"))
    assert h.gnomad_info_keys() == ["gnomad_exomes_AF", "gnomAD_genomes_AF_nfe", "gnomAD_AFR_AF", "AF_grpmax"]


def test_e11_a_contested_clinvar_record_does_not_exempt_a_common_allele(fake_ensembl, tmp_path):
    """"Conflicting classifications" is not an assertion a frequency filter should yield to."""
    rec = json.loads(json.dumps(_vep_record("7 117587806 . G A . . .")))
    for cv in rec["colocated_variants"]:
        cv["clin_sig"] = ["benign", "likely_benign", "pathogenic"]
        cv["clin_sig_allele"] = "A:benign;A:likely_benign;A:pathogenic"
        cv["frequencies"] = {"A": {"gnomade": 0.04, "gnomadg": 0.04}}
    a = V.annotate_vep(rec)
    assert V.clinvar_conflicting(a) and V.clinvar_pathogenic(a) == []

    records = {r["input"]: r for r in json.loads((FIX / "vep_trio.json").read_text())["records"]}
    records[rec["input"]] = rec

    def patched(variants, assembly="GRCh38"):
        lines = [ensembl.vcf_line(*v) for v in variants]
        return Outcome([records[l] for l in lines if l in records],
                       sources=[source_record("Ensembl VEP", "fixture", url="https://rest.ensembl.org/vep")])

    import zebra.sources.ensembl as ens

    original = ens.vep_batch
    ens.vep_batch = patched
    try:
        got = V.triage(GZ, "P", mother="M", father="F", sex="female", genes=_genes())
    finally:
        ens.vep_batch = original
    assert "7-117587806-G-A" not in {c["variant"] for c in got.result["candidates"]}
    assert any("conflicting ClinVar classifications" in w for w in got.warnings)


def test_s3_a_bottlenecked_group_is_still_a_frequency(fake_ensembl):
    """grpmax excludes bottlenecked groups; a frequency *filter* must not treat them as absent."""
    only_fin = V.af_summary({"gnomadg_fin": 0.3})
    assert only_fin["filter_af"] == pytest.approx(0.3)
    assert only_fin["bottlenecked_only"] == "gnomadg_fin" and "bottlenecked" in only_fin["source"]
    # when a non-bottlenecked frequency exists, grpmax semantics are unchanged
    both = V.af_summary({"gnomadg": 0.0001, "gnomadg_nfe": 0.0002, "gnomadg_fin": 0.3})
    assert both["grpmax_group"] == "gnomadg_nfe" and both.get("bottlenecked_only") is None


def test_f22_records_without_gt_and_positions_off_the_contig_are_reported(fake_ensembl, tmp_path):
    vcf = _one_sample_vcf(tmp_path, [
        "chr7\t117587806\t.\tG\tA\t500\tPASS\t.\tGT:AD:DP:GQ\t0/1:16,15:31:99",
        "chr7\t117587900\t.\tG\tA\t500\tPASS\t.",                       # no FORMAT/GT at all
        "chr7\t900000000\t.\tG\tA\t500\tPASS\t.\tGT:AD:DP:GQ\t0/1:16,15:31:99",
    ], "oddities.vcf")
    got = V.triage(vcf, "P", genes=_genes(), max_af_dominant=0.01)
    assert got.result["counts"]["records_without_gt"] == 1
    assert any("no GT in FORMAT" in n for n in got.result["notes"])
    assert any("past the end of their chromosome" in w for w in got.warnings)


def test_f31_a_de_novo_call_without_parental_gq_says_so(fake_ensembl):
    ok_no_gq = V._parent(V.parse_call(["GT", "AD", "DP"], "0/0:30,0:30"), 1, 10, 20)
    assert ok_no_gq["adequate_ref"] and ok_no_gq["no_gq"]
    cls, origin, flags = V._classify("auto", "female", "het", ok_no_gq, ok_no_gq)
    assert (cls, origin) == ("de_novo", "de_novo")
    assert any("no genotype quality (GQ)" in f for f in flags)
    with_gq = V._parent(V.parse_call(["GT", "AD", "DP", "GQ"], "0/0:30,0:30:99"), 1, 10, 20)
    assert not with_gq["no_gq"]
    assert not any("GQ" in f for f in V._classify("auto", "female", "het", with_gq, with_gq)[2])


def test_f22_a_small_file_that_is_wholly_unreadable_is_refused(fake_ensembl, tmp_path):
    """Seven garbage lines used to produce a confident 'no candidate'."""
    rows = ["chr7 117587806 . G A 500 PASS . GT:DP:GQ 0/1:30:90"] * 7
    with pytest.raises(UsageError, match="100% of the data lines"):
        V.triage(_one_sample_vcf(tmp_path, rows, "all_bad.vcf"), "P")
    # a single bad line among many is a warning, not a refusal
    rows = ["chr7 117587806 . G A 500 PASS . GT:DP:GQ 0/1:30:90"]
    rows += [f"chr7\t{117587806 + i}\t.\tG\tA\t500\tPASS\t.\tGT:DP:GQ\t0/1:30:90" for i in range(20)]
    got = V.triage(_one_sample_vcf(tmp_path, rows, "one_bad2.vcf"), "P")
    assert got.result["counts"]["malformed_lines"] == 1


def test_f21_triage_names_a_corrupt_stream_instead_of_blaming_the_samples(tmp_path):
    whole = Path(FAM2_GZ).read_bytes()
    trunc = tmp_path / "half.vcf.gz"
    trunc.write_bytes(whole[: int(len(whole) * 0.6)])
    with pytest.raises(UsageError, match="truncated or corrupt"):
        V.triage(str(trunc), "P2", mother="M2", father="F2")


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


def test_p1d_an_oversized_or_mosaic_input_is_refused(fake_cnv_sources):
    from zebra import cnv as C

    with pytest.raises(UsageError, match="runs past the end of chromosome"):
        C.card("chr1:1-999999999 loss")
    assert C.MAX_REGION_BP == 50_000_000
    with pytest.raises(UsageError, match="will not walk more than"):
        C.genes_in_region("1", 1, 60_000_000)
    with pytest.raises(UsageError, match="copy-number range"):
        C.card("arr[GRCh38] 15q11.2(23123715_23200000)x1~2")
    with pytest.raises(UsageError, match="after the copy number"):
        C.card("arr[GRCh38] 15q11.2(23123715_23200000)x3 mos")
    with pytest.raises(UsageError, match="names a band but no coordinates"):
        C.parse("junk arr[GRCh38] 15q11.2(1000_2000)x1 junk")


def test_p1d_a_copy_number_that_contradicts_the_report_is_refused(fake_cnv_sources):
    from zebra import cnv as C

    with pytest.raises(UsageError, match="the report says loss and the copy number says gain"):
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


def test_p1d_the_gene_count_bands_are_labelled_section_3(fake_cnv_sources):
    from zebra import cnv as C

    inputs = C.card("chr15:25200000-25500000 loss", copies=1).result["acmg_cnv_inputs"]
    bands = inputs["section_3_gene_number"]["clingen_bands"]
    assert "section 3 bands (3A/3B/3C)" in bands and "2A/2B" not in bands
    assert "<25, 25-34 and >=35" in bands


def test_f22_a_nonsense_genotype_or_ref_is_not_quietly_accepted():
    assert V.parse_gt("-1/1") == ((None, 1), False)  # a negative allele index is not an allele
    assert V.zygosity(V.parse_call(["GT"], "-1/1"), 1) == "het"  # one real ALT copy, one unknown
    rec, why = V.parse_line_why("chr1\t100\t.\tn\tG\t.\t.\t.", 0)
    assert rec is not None and rec.ref == "N"  # N is a reference base VEP understands
    rec, why = V.parse_line_why("chr1\t100\t.\tRY\tG\t.\t.\t.", 0)
    assert rec is None and "is not a DNA sequence" in why
