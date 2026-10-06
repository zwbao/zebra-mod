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
    _no_contract_sources(monkeypatch)


# The Chinese-frequency (W4) and MaveDB (W7) modules call the network; offline card tests replace
# them, and the contract tests below put the real dispatchers back with fake modules behind them.
_REAL_CHINA, _REAL_MAVE = getattr(V, "china_frequencies", None), getattr(V, "mavedb_scores", None)


def _no_contract_sources(monkeypatch):
    monkeypatch.setattr(V, "china_frequencies", lambda vcf, warnings, assembly="GRCh38": None, raising=False)
    monkeypatch.setattr(V, "mavedb_scores", lambda gene, p, c, warnings: None, raising=False)


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
    _no_contract_sources(monkeypatch)  # the W4/W7 modules would go to the network
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


# ===================================================================== v0.2 (W6)
# One test per finding id in the v0.2 reviews; each was checked to fail on the unfixed code.

import sys  # noqa: E402
import types  # noqa: E402

from zebra.core import UsageError  # noqa: E402


def _text_fixture(*parts):
    with open(os.path.join(FIX, *parts), encoding="utf-8") as fh:
        return fh.read()


# ---- E-5: MT variants went to the nuclear query as "MT" (HTTP 500 after ~25 s of retries)

def test_e_5_mt_variant_is_sent_to_the_mitochondrial_query_as_M(monkeypatch):
    seen = {}

    def fake_graphql(query, variables, label):
        seen["query"], seen["vars"] = query, variables
        return Response("https://gnomad.broadinstitute.org/api", 200, "", "t", False), \
            load("gnomad", "mito_m3243ag_r3.json")

    monkeypatch.setattr(gnomad, "_graphql", fake_graphql)
    out = gnomad.variant("MT", 3243, "A", "G", "GRCh38")
    assert "mitochondrial_variant" in seen["query"] and "variant(variantId" not in seen["query"]
    assert seen["vars"]["id"] == "M-3243-A-G"
    r = out.result
    assert r["mitochondrial"] is True and r["found"] is True
    assert (r["ac_het"], r["ac_hom"], r["an"]) == (6, 0, 56383)
    assert r["af_het"] == pytest.approx(6 / 56383) and r["max_heteroplasmy"] == pytest.approx(0.464)
    assert "95-100%" in r["definitions"]
    assert gnomad.variant_id("chrM", 3243, "a", "g") == "M-3243-A-G"


def test_e_5_an_mt_allele_missing_from_the_callset_is_not_called_absent():
    r = gnomad.parse_mito(None, "M-3243-A-C", "u")
    assert r["found"] is False and r["absent"] is None and "not the same as AC = 0" in r["note"]


@pytest.fixture
def offline_mt(monkeypatch):
    """card('m.3243A>G') with VEP, the reference, gnomAD mtDNA, MITOMASTER and ClinVar replayed."""
    seq = "N" * 3300
    seq = seq[:3242] + "A" + seq[3243:]

    monkeypatch.setattr(ensembl, "vep", lambda text, assembly="GRCh38": Outcome(
        load("vep", "mt_tl1_m3243ag_grch38.json"), sources=[{"db": "Ensembl VEP"}]))
    monkeypatch.setattr(ensembl, "sequence", lambda c, s, e, a="GRCh38", strand=1: Outcome(
        seq[s - 1:e], sources=[{"db": "Ensembl sequence"}]))
    monkeypatch.setattr(gnomad, "_graphql", lambda q, v, l: (
        Response("https://gnomad.broadinstitute.org/api", 200, "", "t", False), load("gnomad", "mito_m3243ag_r3.json")))
    monkeypatch.setattr(clinvar, "search", lambda term, retmax=40: Outcome({"term": term, "count": 0, "ids": []}))
    monkeypatch.setattr(clinvar, "summaries", lambda uids: Outcome([]))
    monkeypatch.setattr(V, "_litvar", lambda text, gene: (_ for _ in ()).throw(ImportError("offline")))
    monkeypatch.setattr(V, "request", lambda *a, **k: Response(V.MITOMASTER, 200,
                                                                _text_fixture("mitomap", "mitomaster_3243G.tsv"), "t", False))
    _no_contract_sources(monkeypatch)


def test_e_5_mt_card_reports_mtdna_frequencies_and_no_pm2_inputs(offline_mt):
    out = V.card("m.3243A>G")
    r = out.result
    assert r["gene"] == "MT-TL1"
    pop = r["population"]
    assert pop["mitochondrial"] is True and pop["ac_het"] == 6
    a = r["acmg_inputs"]
    assert a["mitochondrial"] is True and a["gnomad_ac"] is None and a["grpmax_af"] is None
    assert a["mt_af_het"] == pytest.approx(6 / 56383)
    assert not any("gnomAD unavailable" in w or "Invalid chromosome" in w for w in out.warnings)
    assert "PM2_Supporting" not in [s["code"] for s in acmg.suggest(a, inheritance="AD")["suggested"]]


# ---- CP1-4: heteroplasmy on input, MITOMAP

def test_cp1_4_heteroplasmy_is_read_from_the_variant_text(offline_mt):
    r = V.card("m.3243A>G 35%").result
    h = r["mitochondrial"]["heteroplasmy"]
    assert h["level"] == pytest.approx(0.35) and h["class"] == "heteroplasmic"
    assert h["gnomad_max_heteroplasmy"] == pytest.approx(0.464)
    assert r["acmg_inputs"]["heteroplasmy"] == pytest.approx(0.35)
    assert V.card("m.3243A>G", heteroplasmy=0.35).result["mitochondrial"]["heteroplasmy"]["percent"] == 35.0


def test_cp1_4_a_percentage_on_a_nuclear_variant_is_refused():
    assert V.split_heteroplasmy("m.3243A>G heteroplasmy 35%") == ("m.3243A>G", pytest.approx(0.35))
    with pytest.raises(UsageError, match="mtDNA heteroplasmy"):
        V.split_heteroplasmy("7-117559590-ATCT-A 35%")
    with pytest.raises(UsageError):
        V.split_heteroplasmy("m.3243A>G 135%")


def test_cp1_4_mitomap_row_is_parsed_with_its_status_caveat(offline_mt):
    mm = V.card("m.3243A>G").result["mitochondrial"]["mitomap"]
    assert mm["found"] is True and "MELAS" in mm["disease_reported"]
    assert mm["listed_as_disease_mutation"] is True and mm["genbank_sequences_with_variant"] == 10
    assert "Cfrm" in mm["status_note"]


def test_cp1_4_mitomap_refuses_an_indel_without_a_request():
    with pytest.raises(ValueError, match="single-base"):
        V.mitomap(3243, "AT", "A")


# ---- CP1-9: LitVar merged CFTR p.Gly542Arg (VUS) with p.Gly542Ter (842 PMIDs)

def test_cp1_9_litvar_record_for_another_allele_is_excluded_and_counted(monkeypatch):
    from zebra.sources import litvar

    monkeypatch.setattr(litvar, "get_json", lambda url, source, params=None, **kw: Response(
        url, 200, _text_fixture("litvar", "auto_cftr_g542r.json"), "t", False))
    res = litvar.lookup("p.Gly542Arg", gene="CFTR").result
    lit = V._litvar_compact(res, ["rs113993959"], "p.Gly542Arg", "NM_000492.4:c.1624G>A")
    assert lit["records"] == [] and lit["pmid_count_max"] is None
    assert lit["excluded"][0]["litvar_id"] == "litvar@rs113993959##"
    assert lit["excluded"][0]["pmid_count"] == 842 and lit["pmids_excluded_max"] == 842
    assert "another allele" in lit["excluded"][0]["reason"]
    same = V._litvar_compact(res, ["rs113993959"], "p.Gly542Ter", None)
    assert same["records"][0]["pmid_count"] == 842 and same["excluded"] == []


# ---- probe P2: ClinVar review status as "n of 4 stars"

def test_probe_clinvar_review_status_is_rendered_as_n_of_4_stars():
    assert clinvar.stars_text(3).startswith("3 of 4 stars")
    assert "4 = practice guideline" in clinvar.stars_text(3)
    res = load("clinvar", "esummary_f508del_locus.json")["result"]
    rec = next(clinvar.parse_summary(res[u]) for u in res["uids"] if res[u].get("accession") == "VCV000007105")
    assert clinvar.compact(rec)["stars_text"].startswith("4 of 4 stars")


# ---- E-13: survivors of the stream-E mutation run

def test_e_13_clinvar_expert_panel_is_three_stars():
    assert clinvar.stars("reviewed by expert panel") == 3


def test_e_13_a_non_overlapping_mane_transcript_loses_without_most_severe():
    tcs = [{"gene_symbol": "NEAR", "mane_select": "NM_1.1", "distance": 64,
            "consequence_terms": ["upstream_gene_variant"], "biotype": "protein_coding"},
           {"gene_symbol": "HERE", "consequence_terms": ["non_coding_transcript_exon_variant"], "biotype": "Mt_tRNA"}]
    assert V.pick_transcript(tcs)["gene_symbol"] == "HERE"


# ---- E-11: a wrong REF in a gene-only HGVS was reported as a source failure

def test_e_11_vep_ref_mismatch_is_bad_input_not_a_source_failure(monkeypatch):
    def boom(base, text, params):
        raise SourceError("Ensembl VEP", "https://rest.ensembl.org/vep/human/hgvs/DMD:c.6439C>T", 400,
                          '{"error":"Reference allele extracted from DMD:c.6439 (G) does not match reference allele '
                          'given by HGVS notation DMD:c.6439C>T (C)"}')

    monkeypatch.setattr(ensembl, "_vep_hgvs", boom)
    with pytest.raises(UsageError, match="reference base in this HGVS"):
        ensembl.vep("DMD:c.6439C>T")


# ---- C-P2-8 / E-6: refetch paths never repaired the cache

def test_c_p2_8_gnomad_refetch_writes_the_good_answer_back(monkeypatch):
    calls = []

    def fake_post(url, payload, source, **kw):
        calls.append(kw)
        if len(calls) == 1:
            return Response(url, 200, json.dumps({"errors": [{"message": "Internal error"}]}), "t", True)
        return Response(url, 200, json.dumps({"data": {}}), "t", False)

    monkeypatch.setattr(gnomad, "post_json", fake_post)
    gnomad._graphql("query {}", {}, "x")
    assert calls[1].get("refresh") is True and calls[1].get("cache_ttl", 0) > 0


def test_c_p2_8_clinvar_refresh_writes_the_good_answer_back(monkeypatch):
    seen = {}

    def fake_get(url, source, params=None, **kw):
        seen.update(kw)
        return Response(url, 200, "{}", "t", False)

    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    monkeypatch.setattr(clinvar, "get_json", fake_get)
    clinvar._eutils("esearch.fcgi", {"term": "x"}, refresh=True)
    assert seen["refresh"] is True and seen["cache_ttl"] == clinvar.CACHE_TTL


# ---- probe P2-3: suggest said "allele number not reported" for a site the card had called covered

def test_sp_p2_3_pm2_basis_carries_the_card_coverage(offline):
    a = V.card("NM_001165963.4:c.2134C>T").result["acmg_inputs"]
    assert a["site_covered"] is True and "of samples at >=20x" in a["coverage_text"]
    pm2 = next(s for s in acmg.suggest(a, inheritance="AD")["suggested"] if s["code"] == "PM2_Supporting")
    assert "site covered" in pm2["basis"] and "allele number not reported" not in pm2["basis"]


# ---- CP1-2: GRCh37 input annotated on the mapped GRCh38 coordinates too

@pytest.fixture
def offline_37(monkeypatch):
    def vep(text, assembly="GRCh38"):
        name = "scn1a_r712x_grch37" if assembly == "GRCh37" else "scn1a_r712x_grch38"
        return Outcome(load("vep", f"{name}.json"), sources=[{"db": f"Ensembl VEP {assembly}"}])

    def sequence(chrom, start, end, assembly="GRCh38", strand=1):
        s = load("vep", f"scn1a_r712x_{assembly.lower()}.seq.json")
        w = V.Window(s["chrom"], s["start"], s["seq"])
        return Outcome(w.slice(max(start, w.start), min(end, w.end)), sources=[{"db": "Ensembl sequence"}])

    def gvariant(chrom, pos, ref, alt, assembly="GRCh38"):
        ds = "gnomad_r4" if assembly == "GRCh38" else "gnomad_r2_1"
        body = load("gnomad", "scn1a_r712x_absent_r4.json")
        return Outcome(gnomad.parse_variant(body["data"], ds, gnomad.variant_id(chrom, pos, ref, alt)),
                       sources=[{"db": f"gnomAD {ds}"}])

    monkeypatch.setattr(ensembl, "vep", vep)
    monkeypatch.setattr(ensembl, "sequence", sequence)
    monkeypatch.setattr(ensembl, "map_assembly", lambda c, s, e, a, b: Outcome(
        {"chrom": "2", "start": s + (166042334 - 166898844), "end": e + (166042334 - 166898844), "strand": 1,
         "assembly": b}, sources=[{"db": "Ensembl assembly map"}]))
    monkeypatch.setattr(gnomad, "variant", gvariant)
    monkeypatch.setattr(clinvar, "search", lambda term, retmax=40: Outcome({"term": term, "count": 0, "ids": []}))
    monkeypatch.setattr(clinvar, "summaries", lambda uids: Outcome([]))
    monkeypatch.setattr(V, "_litvar", lambda text, gene: (_ for _ in ()).throw(ImportError("offline")))
    _no_contract_sources(monkeypatch)


def test_cp1_2_grch37_input_is_also_annotated_on_grch38(offline_37):
    out = V.card("NM_001165963.4:c.2134C>T", assembly="GRCh37")
    r = out.result
    assert r["builds"]["GRCh37"]["id"] == "2-166898844-G-A"
    assert r["builds"]["GRCh38"]["id"] == "2-166042334-G-A"
    assert r["grch38"]["population"]["dataset"] == "gnomad_r4"
    assert r["acmg_inputs"]["frequency_source"].startswith("gnomAD v4 on the GRCh38 coordinates 2-166042334-G-A")
    assert "spliceai" in (r["predictors"].get("filled_from_grch38") or [])
    assert r["predictors"]["spliceai"] is not None
    assert not any("GRCh37 VEP serves no AlphaMissense or SpliceAI" in w for w in out.warnings)


def test_cp1_2_a_reference_that_changed_between_builds_is_not_carried_over(offline_37, monkeypatch):
    real = ensembl.sequence

    def seq(chrom, start, end, assembly="GRCh38", strand=1):
        got = real(chrom, start, end, assembly, strand)
        return Outcome("T" * len(got.result), sources=got.sources) if assembly == "GRCh38" else got

    monkeypatch.setattr(ensembl, "sequence", seq)
    out = V.card("NM_001165963.4:c.2134C>T", assembly="GRCh37")
    assert "grch38" not in out.result
    assert any("reference changed between builds" in w for w in out.warnings)
    assert any("GRCh38 view could not be built" in w for w in out.warnings)


# ---- contracts with W4 (Chinese frequencies) and W7 (MaveDB), behind guarded imports

def _real_contracts(monkeypatch):
    monkeypatch.setattr(V, "china_frequencies", _REAL_CHINA)
    monkeypatch.setattr(V, "mavedb_scores", _REAL_MAVE)


def test_contracts_missing_modules_are_named_gaps(offline, monkeypatch):
    _real_contracts(monkeypatch)
    # None in sys.modules makes `import` raise ImportError, as in a build without the module
    monkeypatch.setitem(sys.modules, "zebra.sources.china_freq", None)
    monkeypatch.setitem(sys.modules, "zebra.sources.mavedb", None)
    out = V.card("NM_000492.4:c.1652G>A")
    assert any(w.startswith("Chinese population frequencies: not checked") for w in out.warnings)
    assert any(w.startswith("MaveDB functional scores: not checked") for w in out.warnings)
    assert "population_china" not in out.result and "functional_scores" not in out.result


def test_contracts_china_freq_and_mavedb_are_called_with_the_agreed_signatures(offline, monkeypatch):
    calls = {}

    china = types.ModuleType("zebra.sources.china_freq")

    def lookup(chrom, pos, ref, alt, assembly="GRCh38"):
        calls["china"] = (chrom, pos, ref, alt, assembly)
        return Outcome({"datasets": [{"name": "ChinaMAP", "population": "Chinese", "af": 0.001, "ac": 20,
                                      "an": 20000, "url": "https://example.invalid/chinamap"}],
                        "checked": ["ChinaMAP", "NyuWa"]}, sources=[{"db": "ChinaMAP"}])

    china.lookup = lookup
    mave = types.ModuleType("zebra.sources.mavedb")

    def mlookup(gene, hgvs_p=None, hgvs_c=None):
        calls["mave"] = (gene, hgvs_p, hgvs_c)
        return Outcome({"datasets": [{"urn": "urn:mavedb:00000001-a", "title": "CFTR DMS", "target": "CFTR"}],
                        "scores": [{"urn": "urn:mavedb:00000001-a#1", "hgvs": "p.Gly551Asp", "score": -1.2,
                                    "interpretation": "loss of function"}]}, sources=[{"db": "MaveDB"}])

    mave.lookup = mlookup
    _real_contracts(monkeypatch)
    monkeypatch.setitem(sys.modules, "zebra.sources.china_freq", china)
    monkeypatch.setitem(sys.modules, "zebra.sources.mavedb", mave)
    out = V.card("NM_000492.4:c.1652G>A")
    r = out.result
    assert calls["china"] == ("7", 117587806, "G", "A", "GRCh38")
    assert calls["mave"][0] == "CFTR" and calls["mave"][1] == "p.Gly551Asp"
    assert r["population_china"]["datasets"][0]["name"] == "ChinaMAP"
    assert r["functional_scores"]["mavedb"]["scores"][0]["score"] == -1.2
    assert any(s.get("db") == "ChinaMAP" for s in out.sources) and any(s.get("db") == "MaveDB" for s in out.sources)
    from zebra.commands.variant import render

    text = render(r)
    assert "Chinese cohorts: ChinaMAP" in text and "MaveDB: urn:mavedb:00000001-a#1" in text


def test_contracts_a_failing_contract_module_is_a_warning_not_a_crash(offline, monkeypatch):
    china = types.ModuleType("zebra.sources.china_freq")

    def lookup(*a, **k):
        raise SourceError("ChinaMAP", "https://example.invalid", 503, "down")

    china.lookup = lookup
    _real_contracts(monkeypatch)
    monkeypatch.setitem(sys.modules, "zebra.sources.china_freq", china)
    monkeypatch.setitem(sys.modules, "zebra.sources.mavedb", None)
    out = V.card("NM_000492.4:c.1652G>A")
    assert "population_china" not in out.result
    assert any(w.startswith("Chinese population frequencies unavailable") for w in out.warnings)


# ---- CP1-3: PVS1 / PS1 / PM5 inputs, end to end offline

def test_cp1_3_judgement_inputs_for_scn1a_r712x(offline, monkeypatch):
    from zebra.sources import clingen

    t = load("ensembl", "transcript_ENST00000674923_scn1a.json")
    monkeypatch.setattr(ensembl, "transcript", lambda tid, assembly="GRCh38": Outcome(
        {"id": t["id"], "strand": t["strand"], "exons": [{"start": e["start"], "end": e["end"]} for e in t["Exon"]],
         "translation": {"start": t["Translation"]["start"], "end": t["Translation"]["end"],
                         "length": t["Translation"]["length"]}}, sources=[{"db": "Ensembl lookup"}]))
    monkeypatch.setattr(clingen, "dosage", lambda s, assembly="GRCh38": Outcome(
        {"haploinsufficiency": {"score": "3", "description": "Sufficient Evidence for Haploinsufficiency"},
         "url": "https://example.invalid/scn1a"}))
    monkeypatch.setattr(gnomad, "gene_constraint", lambda s, a="GRCh38", gene_id=None: Outcome(
        gnomad.parse_constraint(load("gnomad", "constraint_scn1a.json")["data"]["gene"], "GRCh38")))
    r = V.card("NM_001165963.4:c.2134C>T").result
    ji = V.judgement_inputs(r).result
    pv = ji["pvs1_inputs"]
    assert pv["nmd"]["nmd_predicted"] is True and pv["nmd"]["ptc_exon"] == "15/29"
    assert pv["lof_mechanism"]["clingen_hi_score"] == "3" and pv["lof_mechanism"]["gnomad_loeuf"] < 0.2
    assert "not a PVS1 call" in pv["decision_tree"]
    assert "ps1_pm5_inputs" not in ji


def test_cp1_3_judgement_inputs_for_a_missense_list_the_codon_records(offline, monkeypatch):
    res = load("clinvar", "esummary_cftr_codon551.json")["result"]
    recs = [clinvar.parse_summary(res[u]) for u in res["uids"]]
    seen = {}

    def at_positions(chrom, start, end, assembly="GRCh38"):
        seen["range"] = (chrom, start, end, assembly)
        return Outcome(recs, sources=[{"db": "ClinVar"}])

    monkeypatch.setattr(clinvar, "at_positions", at_positions)
    r = V.card("NM_000492.4:c.1652G>A").result
    ji = V.judgement_inputs(r).result
    assert seen["range"] == ("7", 117587804, 117587808, "GRCh38")
    cr = ji["ps1_pm5_inputs"]
    # the card has no ClinVar match offline: G551D's own record is still not its own PS1 input
    assert cr["ps1_inputs"] == []
    assert "VCV000007142" in [x["vcv"] for x in cr["pm5_inputs"]]
    assert "pvs1_inputs" not in ji


# ---- live (v0.2)

@pytest.mark.live
def test_live_e_5_m3243ag_gets_mtdna_frequencies():
    out = V.card("m.3243A>G 35%")
    pop = out.result["population"]
    assert pop["mitochondrial"] is True and pop["found"] is True and pop["ac_het"] >= 1 and pop["an"] > 50000
    assert not any("Invalid chromosome" in w for w in out.warnings)
    assert out.result["mitochondrial"]["heteroplasmy"]["level"] == pytest.approx(0.35)


@pytest.mark.live
def test_live_cp1_2_grch37_g542r_is_annotated_on_grch38():
    r = V.card("7-117227832-G-A", assembly="GRCh37").result
    assert r["builds"]["GRCh38"]["id"] == "7-117587778-G-A"
    assert r["predictors"]["alphamissense"] is not None and r["grch38"]["population"]["dataset"] == "gnomad_r4"


@pytest.mark.live
def test_live_cp1_9_g542r_does_not_inherit_g542x_literature():
    lit = V.card("NM_000492.4:c.1624G>A").result["literature"]["litvar"]
    assert all(x.get("litvar_id") != "litvar@rs113993959##" for x in lit["records"])
    assert any(x.get("litvar_id") == "litvar@rs113993959##" for x in lit["excluded"])


# ---- adversarial review (round 1) findings on the code above

def test_review_b_p1_2_mtdna_flags_hold_when_gnomad_does_not_answer(offline_mt, monkeypatch):
    def down(q, v, l):
        raise SourceError("gnomAD", "https://gnomad.broadinstitute.org/api", 503, "down")

    monkeypatch.setattr(gnomad, "_graphql", down)
    out = V.card("m.3243A>G 35%")
    a = out.result["acmg_inputs"]
    assert a["mitochondrial"] is True and a["heteroplasmy"] == pytest.approx(0.35) and a["gnomad_ac"] is None
    s = acmg.suggest(a, inheritance="AD")
    assert any("mtDNA" in n for n in s["not_assessed"]) and "PM2_Supporting" not in [x["code"] for x in s["suggested"]]
    assert "not computed for mtDNA" in V.judgement_inputs(out.result).result["note"]
    assert any("not absence" in w for w in out.warnings)


def test_review_b_p1_3_an_rsid_query_is_not_gene_filtered(monkeypatch):
    from zebra.sources import litvar

    row = [{"_id": "litvar@rs199474657##", "rsid": "rs199474657", "gene": ["MT-ND1", "MT-ND2"], "hgvs": "p.A3243G",
            "name": "p.A3243G", "pmids_count": 96}]
    monkeypatch.setattr(litvar, "get_json", lambda url, source, params=None, **kw: Response(url, 200, json.dumps(row), "t", False))
    res = litvar.lookup("rs199474657", gene="MT-TL1").result
    assert [m["litvar_id"] for m in res["matches"]] == ["litvar@rs199474657##"]


@pytest.mark.parametrize("behaviour", ["raise_runtime", "raise_usage", "return_tuple", "fail_import"])
def test_review_b_p1_4_any_failure_of_a_contract_module_is_a_named_gap(offline, monkeypatch, behaviour):
    _real_contracts(monkeypatch)
    monkeypatch.setitem(sys.modules, "zebra.sources.mavedb", None)
    if behaviour == "fail_import":
        real_import_module = V.importlib.import_module

        def import_module(name):
            if name == "zebra.sources.china_freq":
                raise RuntimeError("broken at import")
            return real_import_module(name)

        monkeypatch.setattr(V.importlib, "import_module", import_module)
    else:
        china = types.ModuleType("zebra.sources.china_freq")

        def lookup(*a, **k):
            if behaviour == "raise_runtime":
                raise RuntimeError("boom")
            if behaviour == "raise_usage":
                raise UsageError("REF differs from the reference")
            return ("not", "an outcome")

        china.lookup = lookup
        monkeypatch.setitem(sys.modules, "zebra.sources.china_freq", china)
    out = V.card("NM_000492.4:c.1652G>A")
    assert "population_china" not in out.result
    assert any(w.startswith("Chinese population frequencies: not") for w in out.warnings), out.warnings


def test_review_b_p1_4_gtex_failure_is_a_named_gap_in_the_gene_card(monkeypatch):
    from zebra.sources import gene as gene_src

    gtex = types.ModuleType("zebra.sources.gtex")
    gtex.top_tissues = lambda symbol, n=10: (_ for _ in ()).throw(RuntimeError("GTEx parser broke"))
    monkeypatch.setitem(sys.modules, "zebra.sources.gtex", gtex)
    monkeypatch.setattr(gene_src, "hgnc", lambda s: Outcome({"symbol": "SCN1A", "hgnc_id": None}))
    monkeypatch.setattr(gene_src, "_run", lambda tasks: {"gtex": (tasks["gtex"][1](), [])})
    out = gene_src.card("SCN1A")
    assert "expression" not in out.result
    assert any("GTEx expression: not checked — RuntimeError" in w for w in out.warnings)


def test_review_b_p1_5_grch37_judgement_inputs_use_one_transcript_throughout(monkeypatch):
    seen = {}
    t = load("ensembl", "transcript_ENST00000674923_scn1a.json")

    def transcript(tid, assembly="GRCh38"):
        seen["tid"], seen["asm"] = tid, assembly
        return Outcome({"id": t["id"], "strand": t["strand"],
                        "exons": [{"start": e["start"], "end": e["end"]} for e in t["Exon"]],
                        "translation": {"start": t["Translation"]["start"], "end": t["Translation"]["end"],
                                        "length": t["Translation"]["length"]}})

    monkeypatch.setattr(ensembl, "transcript", transcript)
    from zebra.sources import clingen
    monkeypatch.setattr(clingen, "dosage", lambda s, assembly="GRCh38": Outcome(None))
    monkeypatch.setattr(gnomad, "gene_constraint", lambda s, a="GRCh38", gene_id=None: Outcome(None))
    r = {"assembly": "GRCh37", "gene": "SCN1A", "vcf": {"chrom": "2", "pos": 166898844, "ref": "G", "alt": "A"},
         "transcript": {"ensembl": "ENST00000303395.4"}, "hgvs": {"p": "p.Arg712Ter"},
         "consequence": {"terms": ["splice_donor_variant"], "intron": "12/25"}, "acmg_inputs": {},
         "grch38": {"vcf": {"chrom": "2", "pos": 166042334, "ref": "G", "alt": "A"},
                    "transcript": {"ensembl": "ENST00000674923"}, "consequence_terms": ["splice_donor_variant"],
                    "intron": "15/28", "hgvsp": None}}
    pv = V.judgement_inputs(r).result["pvs1_inputs"]
    assert seen == {"tid": "ENST00000674923", "asm": "GRCh38"}
    assert pv["canonical_splice_site"]["intron"] == "15/28" and pv["canonical_splice_site"]["exon_affected"] == "15/29"


def test_review_b_p1_6_mavedb_interpretation_dict_is_rendered_not_dumped():
    from zebra.commands.variant import render

    r = {"input": "x", "assembly": "GRCh38", "consequence": {}, "predictors": {}, "acmg_inputs": {},
         "functional_scores": {"mavedb": {"datasets": [{}] * 7, "datasets_searched": ["a"] * 6, "scores": [
             {"urn": "urn:mavedb:1", "hgvs": "p.Gly551Asp", "score": -1.2,
              "interpretation": {"label": "Functionally abnormal", "acmg_criterion": "PS3",
                                 "acmg_evidence_strength": "STRONG", "source": "MaveDB primary score calibration"}}]}}}
    text = render(r)
    assert "(Functionally abnormal; calibration says PS3 STRONG)" in text and "{'" not in text
    r["functional_scores"]["mavedb"]["scores"] = []
    assert "in the 6 score set(s) searched" in render(r)


def test_review_b_p2_heteroplasmy_on_an_rsid_or_mt_hgvs_is_deferred_to_the_card():
    assert V.split_heteroplasmy("rs199474657 35%") == ("rs199474657", pytest.approx(0.35))
    assert V.split_heteroplasmy("NC_012920.1:m.3243A>G 35%")[1] == pytest.approx(0.35)


def test_review_b_p2_mitomap_checks_position_allele_and_reference(monkeypatch):
    tsv = _text_fixture("mitomap", "mitomaster_3243G.tsv")
    monkeypatch.setattr(V, "request", lambda *a, **k: Response(V.MITOMASTER, 200, tsv, "t", False))
    assert V.mitomap(3243, "A", "G").result["found"] is True
    assert V.mitomap(3243, "A", "C").result["found"] is False  # another allele: the row is not used
    assert V.mitomap(3244, "A", "G").result["found"] is False  # another position
    with pytest.raises(ValueError, match="reference base"):
        V.mitomap(3243, "C", "G")
    out = V.mitomap(3243, "A", "G")
    assert out.sources[0]["url"] == V.MITOMASTER
    header, row = tsv.splitlines()[0].split("\t"), tsv.splitlines()[1].split("\t")
    drop = [i for i, h in enumerate(header) if h not in ("is_mmut", "is_rtmut")]
    trimmed = "\t".join(header[i] for i in drop) + "\n" + "\t".join(row[i] for i in drop) + "\n"
    monkeypatch.setattr(V, "request", lambda *a, **k: Response(V.MITOMASTER, 200, trimmed, "t", False))
    assert V.mitomap(3243, "A", "G").result["listed_as_disease_mutation"] is None


def test_review_b_p2_mtdna_skips_of_china_and_mavedb_are_named(offline_mt, monkeypatch):
    _real_contracts(monkeypatch)
    out = V.card("m.3243A>G")
    assert any(w.startswith("Chinese population frequencies: not checked for an mtDNA variant") for w in out.warnings)
    assert any(w.startswith("MaveDB functional scores: not checked for an mtDNA variant") for w in out.warnings)


def test_review_b_p2_hg19_chrm_on_grch37_is_warned_about(offline_mt):
    out = V.card("chrM-3243-A-G", assembly="GRCh37")
    assert any("Yoruba" in w and "NC_012920" in w for w in out.warnings)


def test_review_b_gaps_parse_mito_absent_flag_and_population_order():
    v = load("gnomad", "mito_m3243ag_r3.json")["data"]["mitochondrial_variant"]
    r = gnomad.parse_mito(v, "M-3243-A-G", "u")
    assert r["absent"] is False and [p["id"] for p in r["populations"][:2]] == ["afr", "nfe"]
    zero = dict(v, ac_het=0, ac_hom=0)
    assert gnomad.parse_mito(zero, "M-3243-A-G", "u")["absent"] is True


def test_review_b_gaps_map_assembly_refuses_split_reverse_and_resized_mappings(monkeypatch):
    def answer(maps):
        return lambda url, source, **kw: Response(url, 200, json.dumps({"mappings": maps}), "t", False)

    one = {"original": {}, "mapped": {"seq_region_name": "7", "start": 100, "end": 100, "strand": 1}}
    monkeypatch.setattr(ensembl, "get_json", answer([one]))
    assert ensembl.map_assembly("7", 50, 50).result["start"] == 100
    for maps in ([one, one], [dict(one, mapped=dict(one["mapped"], strand=-1))],
                 [dict(one, mapped=dict(one["mapped"], end=105))]):
        monkeypatch.setattr(ensembl, "get_json", answer(maps))
        with pytest.raises(ValueError):
            ensembl.map_assembly("7", 50, 50)


def test_review_b_gaps_clinvar_range_search_uses_the_build_field(monkeypatch):
    terms = []
    monkeypatch.setattr(clinvar, "search", lambda term, retmax=40: (terms.append(term),
                        Outcome({"term": term, "count": 0, "ids": []}))[1])
    monkeypatch.setattr(clinvar, "summaries", lambda uids: Outcome([]))
    clinvar.at_positions("7", 117227834, 117227830, "GRCh37")
    clinvar.at_positions("7", 117587804, 117587808, "GRCh38")
    assert terms == ["7[CHR] AND 117227830:117227834[C37]", "7[CHR] AND 117587804:117587808[CPOS]"]


def test_review_b_gaps_ensembl_transcript_parses_the_real_lookup(monkeypatch):
    text = _text_fixture("ensembl", "transcript_ENST00000674923_scn1a.json")
    seen = {}

    def get(url, source, params=None, **kw):
        seen["url"], seen["params"] = url, params
        return Response(url, 200, text, "t", False)

    monkeypatch.setattr(ensembl, "get_json", get)
    t = ensembl.transcript("ENST00000674923.1").result
    assert seen["url"].endswith("/lookup/id/ENST00000674923") and seen["params"] == {"expand": 1}
    assert t["strand"] == -1 and len(t["exons"]) == 29 and t["translation"]["length"] == 2009


def test_eval_bare_c_change_with_a_gene_is_read_as_gene_hgvs(monkeypatch):
    """A report's "SCN1A c.2134C>T" arrives as variant "c.2134C>T" plus gene "SCN1A": it is annotated, not refused."""
    from zebra.sources import ensembl
    from zebra.sources import variant as V

    seen = []

    def fake_vep(text, assembly):
        seen.append(text)
        raise RuntimeError("stop here")

    monkeypatch.setattr(ensembl, "vep", fake_vep)
    try:
        V.card("c.2134C>T", gene="SCN1A")
    except RuntimeError:
        pass
    assert seen == ["SCN1A:c.2134C>T"]
