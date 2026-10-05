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
REAL_HPO = Path(os.path.expanduser("~/.cache/zebra-mod/hpo"))


# ------------------------------------------------------------------ mocks

@pytest.fixture(autouse=True)
def _moi_rows(monkeypatch):
    """Gene inheritance modes from a trimmed copy of HPO genes_to_phenotype.txt (release 2026-09-01)."""
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(FIX / "hpo"))


@pytest.fixture
def fake_ensembl(monkeypatch):
    """VEP and lookup answered from responses captured from rest.ensembl.org (see make_fixtures.py)."""
    vep = json.loads((FIX / "vep_trio.json").read_text())
    by_input = {r["input"]: r for r in vep["records"]}
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
