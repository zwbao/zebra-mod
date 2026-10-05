"""Variant card: normalisation, ClinVar matching, gnomAD/ClinVar parsing (offline, real captured responses)
and the acceptance variants live (pytest -m live)."""

import json
import os

import pytest

from zebra import acmg
from zebra.core import Outcome
from zebra.sources import clinvar, ensembl, gnomad
from zebra.sources import variant as V

FIX = os.path.join(os.path.dirname(__file__), "fixtures")


def load(*parts):
    with open(os.path.join(FIX, *parts)) as fh:
        return json.load(fh)


def window(name):
    s = load("vep", f"{name}.seq.json")
    return V.Window(s["chrom"], s["start"], s["seq"])


def vcf_of(name):
    rec = load("vep", f"{name}.json")
    w = window(name)
    _, alts = V.vep_alleles(rec)
    pos, ref, alt, mismatch = V.raw_vcf(rec, alts[0], w.base)
    assert mismatch is None
    return (rec["seq_region_name"],) + V.left_normalize(pos, ref, alt, w.base)


# ---------------------------------------------------------------- normalisation

def test_revcomp():
    assert V.revcomp("ACGTN") == "NACGT"


def test_scn1a_minus_strand_hgvs_gives_forward_vcf():
    # VEP: start 166042334, strand -1, allele_string C/T (transcript strand); reference base is G
    rec = load("vep", "scn1a_r712x_grch38.json")
    assert rec["strand"] == -1 and rec["allele_string"] == "C/T"
    assert vcf_of("scn1a_r712x_grch38") == ("2", 166042334, "G", "A")


def test_scn1a_grch37():
    assert vcf_of("scn1a_r712x_grch37") == ("2", 166898844, "G", "A")


@pytest.mark.parametrize("name", ["cftr_f508del_hgvs_grch38", "cftr_f508del_vcf_grch38", "cftr_rs113993960_grch38"])
def test_f508del_left_aligned_from_every_input(name):
    # HGVS gives the 3'-shifted CTT/- at 117559592; VCF input TCT/- at 117559591; rsID TCTT/T/TCTTCTT
    assert vcf_of(name) == ("7", 117559590, "ATCT", "A")


def test_duplication_left_aligned():
    w = window("cftr_rs113993960_grch38")
    rec = load("vep", "cftr_rs113993960_grch38.json")
    _, alts = V.vep_alleles(rec)
    assert alts == ["T", "TCTTCTT"]
    pos, ref, alt, _ = V.raw_vcf(rec, alts[1], w.base)
    assert V.left_normalize(pos, ref, alt, w.base) == (117559590, "A", "ATCT")


def test_ref_mismatch_is_reported():
    rec = dict(load("vep", "cftr_g551d_grch38.json"))
    rec["allele_string"] = "C/A"  # reference at 117587806 is G
    pos, ref, alt, mismatch = V.raw_vcf(rec, "A", window("cftr_g551d_grch38").base)
    assert mismatch and "reference has G" in mismatch


def test_out_of_window_raises():
    w = V.Window("7", 100, "ACGT")
    with pytest.raises(V.OutOfWindow):
        w.base(99)


def test_spdi_matches_vcf():
    w = window("cftr_f508del_hgvs_grch38")
    _, pos0, d, i = V.parse_spdi("NC_000007.14:117559590:TCTT:T")
    assert V.left_normalize(*V.spdi_to_vcf(pos0, d, i, w.base), w.base) == (117559590, "ATCT", "A")


# ---------------------------------------------------------------- VEP-derived fields

def test_predictors_g551d():
    rec = load("vep", "cftr_g551d_grch38.json")
    tc = V.pick_transcript(rec["transcript_consequences"])
    assert tc["mane_select"] == "NM_000492.4"
    p = V.predictors(tc)
    assert p["revel"] == pytest.approx(0.99)
    assert p["alphamissense"] == {"score": pytest.approx(0.9897), "class": "likely_pathogenic"}
    assert p["spliceai"]["max"] == 0 and p["cadd_phred"] > 20
    assert p["sift"]["prediction"].startswith("deleterious")


def test_pick_transcript_prefers_given_gene_then_mane():
    tcs = [{"gene_symbol": "CFTR-AS1", "canonical": 1, "biotype": "lncRNA"},
           {"gene_symbol": "CFTR", "mane_select": "NM_000492.4"}]
    assert V.pick_transcript(tcs)["gene_symbol"] == "CFTR"
    assert V.pick_transcript(tcs, gene="CFTR-AS1")["gene_symbol"] == "CFTR-AS1"


def test_vep_frequency_fallback():
    rec = load("vep", "cftr_f508del_hgvs_grch38.json")
    fb = V.vep_frequencies(rec, "-")
    assert fb["rsid"] == "rs113993960"
    assert fb["grpmax_group"].endswith("nfe") and 0.01 < fb["grpmax_af"] < 0.02


def test_clinvar_candidates_from_var_synonyms():
    vcvs, rsid = V._clinvar_candidates(load("vep", "scn1a_r712x_grch38.json"))
    assert "VCV000189886" in vcvs and rsid == "rs794726730"


# ---------------------------------------------------------------- gnomAD parsing

def test_gnomad_present_v4():
    body = load("gnomad", "f508del_r4.json")
    r = gnomad.parse_variant(body["data"], "gnomad_r4", "7-117559590-ATCT-A")
    assert r["found"] and not r["absent"] and r["caid"] == "CA118639"
    assert r["total"]["ac"] == 19237 and r["total"]["an"] == 1612320
    assert r["grpmax"]["group"] == "nfe" and r["grpmax"]["an"] == 1178514
    assert r["grpmax"]["af"] == pytest.approx(17610 / 1178514)
    assert r["faf95"]["value"] == pytest.approx(0.01475724) and r["faf95"]["group"] == "nfe"
    assert "discrepant_frequencies" in r["filters"]
    ids = {p["id"] for p in r["populations"]}
    assert "nfe" in ids and not any("_" in i or ":" in i or i in ("XX", "XY", "") for i in ids)
    assert r["liftover"]["variant_id"] == "7-117199644-ATCT-A"
    assert r["covered"] is True


def test_gnomad_absent_with_coverage():
    body = load("gnomad", "scn1a_r712x_absent_r4.json")
    assert body["errors"][0]["message"] == "Variant not found"
    r = gnomad.parse_variant(body["data"], "gnomad_r4", "2-166042334-G-A")
    assert r["absent"] and not r["found"]
    assert r["coverage"]["exome"]["pos"] == 166042334 and r["coverage"]["exome"]["mean"] > 20
    assert r["covered"] is True and r["total"]["ac"] == 0


def test_gnomad_v2_grpmax_sums_exome_and_genome():
    body = load("gnomad", "f508del_r2_1.json")
    r = gnomad.parse_variant(body["data"], "gnomad_r2_1", "7-117199644-ATCT-A")
    assert r["joint"] is None and r["total"]["ac"] == 1776 + 251
    assert r["grpmax"]["group"] == "nfe" and r["grpmax"]["an"] == 113626 + 15408
    assert r["grpmax"]["af"] == pytest.approx((1394 + 204) / (113626 + 15408))
    assert r["faf95"]["value"] == pytest.approx(0.01175284) and "genome" in r["faf95"]["basis"]
    assert {"asj", "fin", "oth"} <= {p["id"] for p in r["populations"]}  # listed, but never grpmax


def test_gnomad_constraint():
    c = gnomad.parse_constraint(load("gnomad", "constraint_scn1a.json")["data"]["gene"], "GRCh38")
    assert c["loeuf"] == pytest.approx(0.107, abs=1e-3) and c["pLI"] == pytest.approx(1)
    assert c["mis_z"] > 5
    s = gnomad.parse_constraint(load("gnomad", "constraint_smn1.json")["data"]["gene"], "GRCh38")
    assert s["pLI"] is None and "no_exp_lof" in s["flags"]


def test_gnomad_query_shapes():
    assert "joint" in gnomad.build_variant_query("gnomad_r4")
    assert "joint" not in gnomad.build_variant_query("gnomad_r2_1")
    assert gnomad.variant_id("chr7", 117559590, "atct", "a") == "7-117559590-ATCT-A"


# ---------------------------------------------------------------- ClinVar parsing and matching

def test_clinvar_parse_and_stars():
    res = load("clinvar", "esummary_f508del_locus.json")["result"]
    r = clinvar.parse_summary(res["7105"])
    assert r["vcv"] == "VCV000007105" and r["classification"] == "Pathogenic"
    assert r["review_status"] == "practice guideline" and r["stars"] == 4
    assert r["canonical_spdi"] == "NC_000007.14:117559590:TCTT:T" and r["caid"] == "CA118639"
    assert r["rsid"] == "rs113993960" and r["last_evaluated"] == "2004-03-03"
    assert any("MONDO:0009061" in c["ids"] for c in r["conditions"])
    assert r["submissions_scv"] == len(res["7105"]["supporting_submissions"]["scv"])
    assert clinvar.parse_summary(res["634837"])["compound"] is True
    assert clinvar.stars("criteria provided, multiple submitters, no conflicts") == 2
    assert clinvar.stars("criteria provided, conflicting classifications") == 1
    assert clinvar.stars("no assertion criteria provided") == 0


def test_clinvar_term():
    t = clinvar.term_for(vcv=["VCV000189886"], rsid="rs794726730", chrom="2", pos=166042334)
    assert t == "189886[VID] OR rs794726730[VRID] OR (2[CHR] AND 166042334[CPOS])"
    assert "[C37]" in clinvar.term_for(chrom="2", pos=1, assembly="GRCh37")
    assert clinvar.vcv_uid("VCV000007105.213") == "7105"


def test_clinvar_match_exact_allele_rejects_haplotype():
    res = load("clinvar", "esummary_f508del_locus.json")["result"]
    recs = [clinvar.parse_summary(res[u]) for u in res["uids"]]
    exact, others = V.match_clinvar(recs, ("7", 117559590, "ATCT", "A"), "GRCh38", window("cftr_f508del_hgvs_grch38"), None)
    assert [r["vcv"] for r in exact] == ["VCV000007105"]
    assert [r["vcv"] for r in others] == ["VCV004073507"]  # haplotype 634837 dropped


def test_clinvar_match_grch37_shifts_spdi():
    res = load("clinvar", "esummary_189886.json")["result"]
    recs = [clinvar.parse_summary(res[u]) for u in res["uids"]]
    exact, _ = V.match_clinvar(recs, ("2", 166898844, "G", "A"), "GRCh37", window("scn1a_r712x_grch37"), None)
    assert [r["vcv"] for r in exact] == ["VCV000189886"]
    exact, _ = V.match_clinvar(recs, ("2", 166898844, "G", "T"), "GRCh37", window("scn1a_r712x_grch37"), None)
    assert exact == []


# ---------------------------------------------------------------- card, offline end to end

@pytest.fixture
def offline(monkeypatch):
    """card() with every network call replaced by captured responses."""
    cases = {
        "NM_001165963.4:c.2134C>T": ("scn1a_r712x_grch38", "scn1a_r712x_absent_r4", "esummary_189886"),
        "NM_000492.4:c.1521_1523del": ("cftr_f508del_hgvs_grch38", "f508del_r4", "esummary_f508del_locus"),
        "NM_000492.4:c.1652G>A": ("cftr_g551d_grch38", "g551d_r4", None),
    }
    state = {}

    def vep(text, assembly="GRCh38"):
        state["case"] = cases[text]
        return Outcome(load("vep", f"{cases[text][0]}.json"), sources=[{"db": "Ensembl VEP"}])

    def sequence(chrom, start, end, assembly="GRCh38", strand=1):
        s = load("vep", f"{state['case'][0]}.seq.json")
        w = V.Window(s["chrom"], s["start"], s["seq"])
        return Outcome(w.slice(max(start, w.start), min(end, w.end)), sources=[{"db": "Ensembl sequence"}])

    def gvariant(chrom, pos, ref, alt, assembly="GRCh38"):
        body = load("gnomad", f"{state['case'][1]}.json")
        return Outcome(gnomad.parse_variant(body["data"], "gnomad_r4", gnomad.variant_id(chrom, pos, ref, alt)),
                       sources=[{"db": "gnomAD"}])

    def search(term, retmax=40):
        return Outcome({"term": term, "count": 1, "ids": []}, sources=[{"db": "ClinVar esearch"}])

    def summaries(uids):
        name = state["case"][2]
        if not name:
            return Outcome([])
        res = load("clinvar", f"{name}.json")["result"]
        return Outcome([clinvar.parse_summary(res[u]) for u in res["uids"]], sources=[{"db": "ClinVar"}])

    def no_litvar(text, gene):
        raise ImportError("litvar not here")

    monkeypatch.setattr(ensembl, "vep", vep)
    monkeypatch.setattr(ensembl, "sequence", sequence)
    monkeypatch.setattr(gnomad, "variant", gvariant)
    monkeypatch.setattr(clinvar, "search", search)
    monkeypatch.setattr(clinvar, "summaries", summaries)
    monkeypatch.setattr(V, "_litvar", no_litvar)


def test_card_offline_scn1a(offline):
    out = V.card("NM_001165963.4:c.2134C>T")
    r = out.result
    assert r["vcf"]["id"] == "2-166042334-G-A"
    assert r["hgvs"]["c"] == "NM_001165963.4:c.2134C>T" and r["hgvs"]["p"] == "p.Arg712Ter"
    assert r["consequence"]["terms"] == ["stop_gained"] and r["consequence"]["exon"] == "15/29"
    assert r["clinvar"]["vcv"] == "VCV000189886" and r["clinvar"]["stars"] >= 1
    a = r["acmg_inputs"]
    assert a["gnomad_ac"] == 0 and a["grpmax_af"] is None and a["is_missense"] is False
    assert any("LitVar unavailable" in w for w in out.warnings)
    assert any(s["code"] == "PM2_Supporting" for s in acmg.suggest_from_data(a))


def test_card_offline_f508del(offline):
    r = V.card("NM_000492.4:c.1521_1523del").result
    assert r["vcf"]["id"] == "7-117559590-ATCT-A"
    assert r["clinvar"]["vcv"] == "VCV000007105" and r["clinvar"]["stars"] == 4
    assert r["consequence"]["terms"] == ["inframe_deletion"]
    a = r["acmg_inputs"]
    assert a["gnomad_ac"] == 19237 and 0.005 < a["grpmax_af"] < 0.02 and a["grpmax_an"] > 2000
    assert a["faf95"] == pytest.approx(0.01475724)


def test_card_offline_g551d_suggest(offline):
    r = V.card("NM_000492.4:c.1652G>A", gene="SCN1A")
    a = r.result["acmg_inputs"]
    assert a["revel"] == pytest.approx(0.99) and a["is_missense"]
    codes = {s["code"] for s in acmg.suggest_from_data(a, inheritance="AR")}
    assert {"PP3_Strong", "PM2_Supporting"} <= codes
    assert any("gene mismatch" in w for w in r.warnings)
    assert r.result["clinvar"] is None and "no ClinVar record" in r.result["clinvar_note"]


def test_render_text(offline):
    from zebra.commands.variant import render

    text = render(V.card("NM_000492.4:c.1521_1523del").result)
    assert "7-117559590-ATCT-A" in text and "ClinVar: Pathogenic" in text and "acmg_inputs:" in text


# ---------------------------------------------------------------- live

@pytest.mark.live
def test_live_f508del_hgvs():
    r = V.card("NM_000492.4:c.1521_1523del").result
    assert r["vcf"]["id"] == "7-117559590-ATCT-A"
    assert "inframe_deletion" in r["consequence"]["terms"]
    assert r["clinvar"]["classification"].startswith("Pathogenic") and r["clinvar"]["stars"] >= 2
    assert 0.005 < r["population"]["total"]["af"] < 0.02


@pytest.mark.live
def test_live_f508del_vcf():
    r = V.card("7-117559590-ATCT-A").result
    assert r["hgvs"]["p"] == "p.Phe508del" and r["clinvar"]["vcv"] == "VCV000007105"


@pytest.mark.live
def test_live_scn1a_stop_absent():
    r = V.card("NM_001165963.4:c.2134C>T").result
    assert r["vcf"]["id"] == "2-166042334-G-A"
    assert r["consequence"]["terms"] == ["stop_gained"]
    assert r["population"]["absent"] and r["acmg_inputs"]["gnomad_ac"] == 0
    assert "athogenic" in (r["clinvar"] or {}).get("classification", "")


@pytest.mark.live
def test_live_g551d_predictors_and_suggest():
    r = V.card("NM_000492.4:c.1652G>A").result
    assert r["predictors"]["revel"] >= 0.9 and r["predictors"]["alphamissense"]["class"] == "likely_pathogenic"
    codes = {s["code"] for s in acmg.suggest_from_data(r["acmg_inputs"], inheritance="AR")}
    assert {"PP3_Strong", "PM2_Supporting"} <= codes


@pytest.mark.live
def test_live_grch37():
    r = V.card("2-166898844-G-A", assembly="GRCh37").result
    assert r["gene"] == "SCN1A" and r["consequence"]["terms"] == ["stop_gained"]
    assert r["population"]["dataset"] == "gnomad_r2_1"
    assert r["clinvar"]["vcv"] == "VCV000189886"
