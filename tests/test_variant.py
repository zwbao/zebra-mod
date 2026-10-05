"""Variant card: normalisation, ClinVar matching, gnomAD/ClinVar parsing (offline, real captured responses)
and the acceptance variants live (pytest -m live)."""

import json
import os

import pytest

from zebra import acmg
from zebra.core import Outcome
from zebra.http import Response, SourceError
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


# ---------------------------------------------------------------- E7: never fabricate an allele

def test_E7_mismatch_message_names_both_bases_and_the_matching_build():
    msg = V.ref_mismatch_message("7", 117199644, "ATCT", "CAGC", "GRCh38", "ATCT")
    assert "REF ATCT does not match the GRCh38 reference at 7:117199644" in msg
    assert "which has CAGC" in msg
    assert "GRCh37 has ATCT there, which matches: rerun with --assembly GRCh37" in msg


def test_E7_mismatch_message_when_the_other_build_does_not_match_either():
    msg = V.ref_mismatch_message("7", 117199644, "ATCT", "CAGC", "GRCh38", "GGGG")
    assert "GRCh37 has GGGG there, which does not match either" in msg
    assert "rerun" not in msg


def test_E7_mismatch_message_when_the_other_build_could_not_be_checked():
    msg = V.ref_mismatch_message("7", 1, "A", "C", "GRCh37", None)
    assert "GRCh38 could not be checked" in msg


def test_E7_vcf_input_with_wrong_ref_is_refused_before_vep(monkeypatch):
    """A GRCh37 F508del given as GRCh38 must raise, not reach VEP or gnomAD."""
    reference = {("GRCh38", 117199644): "CAGC", ("GRCh37", 117199644): "ATCT"}

    def sequence(chrom, start, end, assembly="GRCh38", strand=1):
        return Outcome(reference[(assembly, start)][:end - start + 1], sources=[{"db": "Ensembl sequence"}])

    def boom(*a, **kw):
        raise AssertionError("VEP must not be called once the REF check has failed")

    monkeypatch.setattr(ensembl, "sequence", sequence)
    monkeypatch.setattr(ensembl, "vep", boom)
    with pytest.raises(V.RefMismatch) as err:
        V.card("7-117199644-ATCT-A")
    text = str(err.value)
    assert "REF ATCT does not match the GRCh38 reference" in text and "CAGC" in text
    assert "--assembly GRCh37" in text
    # the refusal must not leave room for an "absent from gnomAD" reading
    assert "absence from gnomAD or ClinVar is not evidence" in text


def test_E7_correct_ref_passes_the_check(monkeypatch):
    seen = []

    def sequence(chrom, start, end, assembly="GRCh38", strand=1):
        seen.append((assembly, start, end))
        return Outcome("ATCT"[:end - start + 1], sources=[{"db": "Ensembl sequence"}])

    monkeypatch.setattr(ensembl, "sequence", sequence)
    warnings, sources = [], []
    V.verify_input_ref("7", 117199644, "ATCT", "GRCh38", sources, warnings)
    assert warnings == [] and seen == [("GRCh38", 117199644, 117199647)]


def test_E7_unreachable_reference_warns_instead_of_refusing(monkeypatch):
    def sequence(*a, **kw):
        raise SourceError("Ensembl sequence", "u", 503, "down")

    monkeypatch.setattr(ensembl, "sequence", sequence)
    warnings, sources = [], []
    V.verify_input_ref("7", 1, "A", "GRCh38", sources, warnings)
    assert any("not checked against the GRCh38 reference" in w for w in warnings)


def test_E7_post_vep_mismatch_raises(monkeypatch):
    """An HGVS whose VEP coordinates contradict the reference is refused, not warned about."""
    rec = dict(load("vep", "cftr_g551d_grch38.json"))
    rec["allele_string"] = "C/A"  # the reference at 117587806 is G

    def vep(text, assembly="GRCh38"):
        return Outcome(rec, sources=[])

    def sequence(chrom, start, end, assembly="GRCh38", strand=1):
        s = load("vep", "cftr_g551d_grch38.seq.json")
        w = V.Window(s["chrom"], s["start"], s["seq"])
        return Outcome(w.slice(max(start, w.start), min(end, w.end)), sources=[])

    monkeypatch.setattr(ensembl, "vep", vep)
    monkeypatch.setattr(ensembl, "sequence", sequence)
    with pytest.raises(V.RefMismatch) as err:
        V.card("NM_000492.4:c.1652G>A")
    assert "does not match the GRCh38 reference" in str(err.value)


# ------------------------------------------------- P1e: the transcript must overlap the variant

def test_P1e_m3243ag_annotates_MT_TL1_not_MT_ND1():
    rec = load("vep", "mt_tl1_m3243ag_grch38.json")
    assert rec["most_severe_consequence"] == "non_coding_transcript_exon_variant"
    # the whole list still offers MT-ND1 as a canonical protein-coding transcript 64 bp away
    assert any(t["gene_symbol"] == "MT-ND1" and t.get("canonical") for t in rec["transcript_consequences"])
    tc = V.pick_transcript(rec["transcript_consequences"], None, rec["most_severe_consequence"])
    assert tc["gene_symbol"] == "MT-TL1" and tc["biotype"] == "Mt_tRNA"
    assert V.overlaps(tc) and tc.get("distance") is None
    # MT-TL1 is the only transcript the variant is inside
    assert [t["gene_symbol"] for t in V.overlapping_transcripts(rec["transcript_consequences"])] == ["MT-TL1"]


def test_P1e_non_coding_rna_gene_wins_over_an_overlapping_protein_coding_MANE():
    """RNU4ATAC n.51G>A sits in a CLASP1 intron; both genes have a MANE transcript."""
    rec = load("vep", "rnu4atac_n51ga_grch38.json")
    genes = {t["gene_symbol"] for t in V.overlapping_transcripts(rec["transcript_consequences"])}
    assert {"RNU4ATAC", "CLASP1"} <= genes
    mane = [t for t in rec["transcript_consequences"] if t.get("mane_select") and V.overlaps(t)]
    assert {t["gene_symbol"] for t in mane} == {"RNU4ATAC", "CLASP1"}
    tc = V.pick_transcript(rec["transcript_consequences"], None, rec["most_severe_consequence"])
    assert tc["gene_symbol"] == "RNU4ATAC" and tc["mane_select"] == "NR_023343.3"
    assert "non_coding_transcript_exon_variant" in tc["consequence_terms"]


def test_P1e_intergenic_variant_has_no_gene_and_no_transcript():
    rec = load("vep", "intergenic_2_36000000_grch38.json")
    assert rec["most_severe_consequence"] == "intergenic_variant"
    assert not rec.get("transcript_consequences")
    assert V.pick_transcript(rec.get("transcript_consequences") or []) is None
    assert V.overlapping_transcripts(rec.get("transcript_consequences") or []) == []


def test_P1e_overlaps_rejects_upstream_and_downstream():
    assert not V.overlaps({"distance": 64, "consequence_terms": ["upstream_gene_variant"]})
    assert not V.overlaps({"consequence_terms": ["downstream_gene_variant"]})
    assert V.overlaps({"consequence_terms": ["intron_variant"]})
    assert not V.overlaps({"consequence_terms": []})


def test_P1e_nearest_transcript_is_the_closest_one():
    rec = load("vep", "mt_tl1_m3243ag_grch38.json")
    near = V.nearest_transcript(rec["transcript_consequences"])
    assert near["gene_symbol"] == "MT-RNR2" and near["distance"] == 14


# ------------------------------------------------- F37: transcript versions and bare GENE:c.

def test_F37_hgvs_regex_is_anchored_at_both_ends():
    assert ensembl.classify_input("NM_000492.4:c.1521_1523del") == "hgvs"
    # one argv string carrying a flag must not pass as HGVS and be annotated on the default build
    assert ensembl.classify_input("NM_000492.4:c.1521_1523del --assembly GRCh37") == "unknown"
    assert ensembl.classify_input("SCN1A:c.2134C>T") == "gene_hgvs"
    assert ensembl.classify_input("SCN1A:c.2134C>T extra") == "unknown"


def test_F37_strip_version():
    assert ensembl.strip_version("NM_001165963.1") == "NM_001165963"
    assert ensembl.strip_version("NM_001165963") == "NM_001165963"


def test_F37_legacy_transcript_version_is_retried_without_the_version(monkeypatch):
    calls = []

    def fake_get(url, source, **kw):
        calls.append(url)
        if "NM_001165963.1" in url:
            raise SourceError("Ensembl VEP", url, 400,
                              "Unable to parse HGVS notation 'NM_001165963.1:c.2134C>T': "
                              "Could not get a Transcript object for 'NM_001165963.1'")
        return Response(url, 200, json.dumps([load("vep", "scn1a_r712x_grch38.json")]),
                        "2026-10-06T00:00:00+00:00", False)

    monkeypatch.setattr(ensembl, "get_json", fake_get)
    out = ensembl.vep("NM_001165963.1:c.2134C>T")
    assert len(calls) == 2
    assert calls[0].endswith("NM_001165963.1:c.2134C>T") and calls[1].endswith("NM_001165963:c.2134C>T")
    assert any("not in this Ensembl release" in w and "c. numbering can differ" in w for w in out.warnings)


def test_F37_other_vep_errors_are_not_retried(monkeypatch):
    def fake_get(url, source, **kw):
        raise SourceError("Ensembl VEP", url, 400, "Unable to parse HGVS notation")

    monkeypatch.setattr(ensembl, "get_json", fake_get)
    with pytest.raises(SourceError):
        ensembl.vep("NM_001165963.1:c.2134C>T")


def test_F37_bare_gene_hgvs_warns_which_transcript_was_used(monkeypatch):
    def vep(text, assembly="GRCh38"):
        return Outcome(load("vep", "scn1a_r712x_grch38.json"), sources=[])

    def sequence(chrom, start, end, assembly="GRCh38", strand=1):
        s = load("vep", "scn1a_r712x_grch38.seq.json")
        w = V.Window(s["chrom"], s["start"], s["seq"])
        return Outcome(w.slice(max(start, w.start), min(end, w.end)), sources=[])

    monkeypatch.setattr(ensembl, "vep", vep)
    monkeypatch.setattr(ensembl, "sequence", sequence)
    monkeypatch.setattr(gnomad, "variant", lambda *a, **k: None)
    monkeypatch.setattr(clinvar, "search", lambda *a, **k: Outcome({"term": "x", "count": 0, "ids": []}))
    monkeypatch.setattr(clinvar, "summaries", lambda uids: Outcome([]))
    monkeypatch.setattr(V, "_litvar", lambda text, gene: None)
    out = V.card("SCN1A:c.2134C>T")
    assert any("you gave only a gene name" in w and "NM_001165963.4" in w and
               "c. numbering is transcript-specific" in w for w in out.warnings)


# ---------------------------------------------------------------- live

@pytest.mark.live
def test_live_E7_grch37_coordinates_sent_as_grch38_are_refused():
    with pytest.raises(V.RefMismatch) as err:
        V.card("7-117199644-ATCT-A")
    assert "GRCh37 has ATCT there, which matches" in str(err.value)
    # and the same input with the right build works
    r = V.card("7-117199644-ATCT-A", assembly="GRCh37").result
    assert r["vcf"]["id"] == "7-117199644-ATCT-A" and r["gene"] == "CFTR"


@pytest.mark.live
def test_live_P1e_m3243ag_is_MT_TL1():
    r = V.card("NC_012920.1:m.3243A>G").result
    assert r["gene"] == "MT-TL1"
    assert r["consequence"]["terms"] == ["non_coding_transcript_exon_variant"]
    assert r["transcript"]["overlaps_variant"] is True
    assert r["clinvar"]["classification"].startswith("Pathogenic")


@pytest.mark.live
def test_live_P1e_rnu4atac_and_intergenic():
    r = V.card("2-121530930-G-A").result
    assert r["gene"] == "RNU4ATAC" and r["transcript"]["refseq"] == "NR_023343.3"
    assert any(a["gene"] == "CLASP1" for a in r["transcript"]["also_overlapping"])
    i = V.card("2-36000000-A-G")
    assert i.result["gene"] is None and i.result["transcript"] is None
    assert any("lies in no transcript" in w for w in i.warnings)


@pytest.mark.live
def test_live_F37_legacy_version_resolves_with_a_warning():
    out = V.card("NM_001165963.1:c.2134C>T")
    assert out.result["vcf"]["id"] == "2-166042334-G-A"
    assert any("not in this Ensembl release" in w for w in out.warnings)


# ------------------------------------- S7: grpmax follows gnomAD's own per-dataset definition

def test_S7_grpmax_excluded_groups_per_dataset():
    assert "mid" not in gnomad.grpmax_excluded("joint")
    assert "mid" not in gnomad.grpmax_excluded("exome")
    assert "mid" in gnomad.grpmax_excluded("genome")
    for kind in ("joint", "exome", "genome"):
        assert {"ami", "asj", "fin", "rmi", "remaining"} <= set(gnomad.grpmax_excluded(kind))


def test_S7_mefv_m694v_grpmax_is_the_middle_eastern_group():
    """MEFV p.Met694Val: the founder allele's grpmax group is `mid`, ~20x the `amr` AF."""
    fx = load("gnomad", "s7_mefv_m694v_r4.json")
    r = gnomad.parse_variant(fx["data"], fx["dataset"], fx["variant_id"])
    pops = {p["id"]: p for p in r["populations"]}
    assert pops["mid"]["af"] > 15 * pops["amr"]["af"]  # measured 19.8x on 2026-10-06
    assert r["grpmax"]["group"] == "mid"
    assert r["grpmax"]["af"] == pytest.approx(pops["mid"]["af"])
    # 4.6e-03 is 6.6x above the 7e-04 recessive PM2 threshold, so PM2 must not fire
    assert r["grpmax"]["af"] > 0.0007
    assert "mid" in r["grpmax"]["groups_considered"]
    assert "mid" not in r["grpmax"]["groups_excluded"]
    assert "gnomAD help topic 'grpmax'" in r["grpmax"]["basis"]


def test_S7_grpmax_is_computed_by_exclusion_so_a_new_group_counts():
    groups = {"afr": {"ac": 1, "an": 1000}, "newgrp": {"ac": 50, "an": 1000},
              "fin": {"ac": 900, "an": 1000}}
    best = gnomad.grpmax(groups, "joint")
    assert best["group"] == "newgrp"  # not dropped for being unlisted
    assert best["groups_excluded"] == ["fin"]


def test_S7_faf95_basis_says_it_is_the_group_with_the_highest_faf():
    fx = load("gnomad", "s7_mefv_m694v_r4.json")
    r = gnomad.parse_variant(fx["data"], fx["dataset"], fx["variant_id"])
    basis = r["faf95"]["basis"]
    assert "GroupMax FAF" in basis and "highest FAF" in basis
    assert "not necessarily the group with the highest AF" in basis
    assert r["faf95"]["datasets"] == ["joint"]


def test_S7_v2_faf95_reports_exome_and_genome_separately():
    data = {"variant": {"variant_id": "1-1-A-G", "exome": {"ac": 5, "an": 1000, "faf95": {"popmax": 0.004,
            "popmax_population": "nfe"}, "populations": [{"id": "nfe", "ac": 5, "an": 1000}]},
            "genome": {"ac": 1, "an": 500, "faf95": {"popmax": 0.001, "popmax_population": "afr"},
                       "populations": [{"id": "afr", "ac": 1, "an": 500}]}}}
    r = gnomad.parse_variant(data, "gnomad_r2_1", "1-1-A-G")
    assert r["faf95"]["value"] == 0.004 and r["faf95"]["from"] == "exome"
    assert set(r["faf95"]["per_dataset"]) == {"exome", "genome"}
    assert "NOT a filtering AF computed on the pooled sample" in r["faf95"]["basis"]


# ------------------------------- F33: a covered site is judged on the fraction over 20x

def test_F33_covered_uses_the_fraction_over_20x_not_mean_depth():
    # the real DMD site: exome mean 2.2 / over_20 0.01, genome mean 24.1 / over_20 0.61
    cov = {"exome": {"mean": 2.206, "median": 0, "over_20": 0.01},
           "genome": {"mean": 24.076, "median": 23, "over_20": 0.61}}
    assert gnomad.covered(cov) is False
    assert gnomad.covered({"genome": {"mean": 30.0, "median": 30, "over_20": 0.97}}) is True
    assert gnomad.covered({}) is None


def test_F33_dmd_absent_site_is_not_called_covered():
    fx = load("gnomad", "f33_dmd_x31121491_absent_r4.json")
    r = gnomad.parse_variant(fx["data"], fx["dataset"], fx["variant_id"])
    assert r["found"] is False
    assert r["coverage"]["genome"]["mean"] >= 20          # mean alone said "covered"
    assert r["coverage"]["genome"]["over_20"] < 0.8        # the fraction says otherwise
    assert r["covered"] is False
    assert r["total"]["ac"] is None                        # no AC 0 to feed PM2
    assert "61%" in r["note"] or "0.61" in str(r["coverage_detail"]["fraction_over_20"])
    assert "mean depth is not used" in r["coverage_rule"]
    assert "allele number" in r["total"]["an_note"]


def test_F33_coverage_detail_reports_the_fractions_and_the_threshold():
    cov = {"exome": {"mean": 2.2, "median": 0, "over_20": 0.01},
           "genome": {"mean": 24.1, "median": 23, "over_20": 0.61}}
    d = gnomad.covered_detail(cov)
    assert d["threshold_over_20"] == 0.8
    assert d["fraction_over_20"] == {"exome": 0.01, "genome": 0.61}
    assert d["covered_by"] is None


def test_F33_median_is_the_fallback_when_over_20_is_absent():
    assert gnomad.covered({"genome": {"mean": 30.0, "median": 31}}) is True
    assert gnomad.covered({"genome": {"mean": 30.0, "median": 4}}) is False
    assert "median depth was used" in gnomad.covered_detail({"genome": {"median": 31}})["note"]


def test_F33_vep_fallback_frequencies_split_the_group_sets():
    rec = {"colocated_variants": [{"id": "rs1", "frequencies": {"G": {
        "gnomade_mid": 0.02, "gnomadg_mid": 0.03, "gnomade_nfe": 0.001, "gnomadg_nfe": 0.001}}}]}
    fb = V.vep_frequencies(rec, "G")
    # mid counts for the exome set and not for the genome set
    assert fb["grpmax_group"] == "e:mid" and fb["grpmax_af"] == 0.02
    assert "mid" in fb["grpmax_basis"]


@pytest.mark.live
def test_live_S7_mefv_grpmax_is_mid():
    r = V.card("NM_000243.3:c.2080A>G").result
    g = r["population"]["grpmax"]
    assert g["group"] == "mid" and g["af"] > 0.004
    assert r["acmg_inputs"]["grpmax_af"] > 0.0007


@pytest.mark.live
def test_live_F33_dmd_site_is_not_covered():
    r = V.card("X-31121491-C-T").result
    p = r["population"]
    assert p["absent"] is True and p["covered"] is False
    assert p["coverage"]["genome"]["over_20"] < 0.8
    assert r["acmg_inputs"]["gnomad_ac"] is None


def test_E7_refusal_is_a_usage_error_so_every_command_reports_it_as_bad_input():
    """`zebra acmg suggest` does not wrap card(); a plain ValueError became an InternalError."""
    from zebra.core import UsageError

    assert issubclass(V.RefMismatch, UsageError) and issubclass(V.RefMismatch, ValueError)


@pytest.mark.live
def test_live_E7_acmg_suggest_refuses_the_same_input(capsys):
    from zebra import cli

    code = cli.main(["--json", "acmg", "suggest", "7-117199644-ATCT-A"])
    env = json.loads(capsys.readouterr().out)
    assert code != 0 and env["ok"] is False
    assert env["error"]["type"] == "UsageError"
    assert "does not match the GRCh38 reference" in env["error"]["message"]
    assert "traceback_tail" not in env["error"]
