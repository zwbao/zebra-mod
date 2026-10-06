import pytest

from zebra import acmg


@pytest.mark.parametrize(
    "codes,points,cls",
    [
        (["PVS1", "PS2", "PM2_Supporting"], 13, "Pathogenic"),
        (["PVS1", "PM2_Supporting"], 9, "Likely pathogenic"),
        (["PS3", "PM2_Supporting", "PP3"], 6, "Likely pathogenic"),
        (["PM2_Supporting", "PP3"], 2, "Uncertain significance"),
        (["BS1", "BP4"], -5, "Likely benign"),
        (["BS1", "BS2"], -8, "Benign"),
        (["PS3", "BS1"], 0, "Uncertain significance"),
    ],
)
def test_points(codes, points, cls):
    r = acmg.classify(codes)
    assert r["points"] == points
    assert r["classification"] == cls


def test_ba1_standalone():
    r = acmg.classify(["BA1"])
    assert r["classification"] == "Benign" and r["points"] is None


def test_richards_disagreement_is_reported():
    r = acmg.classify(["PVS1", "PM2_Supporting"])
    assert r["classification_richards_2015"] == "Uncertain significance"
    assert r["agree"] is False and "note" in r


def test_richards_rules():
    assert acmg.classify(["PS1", "PS3"])["classification_richards_2015"] == "Pathogenic"
    assert acmg.classify(["PM1", "PM2", "PM5"])["classification_richards_2015"] == "Likely pathogenic"
    # a lone BS1 meets no benign rule, so the 2015 rules see no contradiction; points weigh it
    r = acmg.classify(["PVS1", "PS3", "BS1"])
    assert r["classification_richards_2015"] == "Pathogenic" and r["classification"] == "Likely pathogenic"
    assert acmg.classify(["PS1", "PM1", "BS1", "BP4"])["classification_richards_2015"] == "Uncertain significance"


def test_pm1_pp3_cap():
    r = acmg.classify(["PM1", "PP3_Strong"])
    assert r["points"] == 4
    assert any("cap" in w for w in r["warnings"])


def test_warnings():
    r = acmg.classify(["PVS1", "PP3", "PP5", "PM2"])
    text = " ".join(r["warnings"])
    assert "double counting" in text and "PP5" in text and "PM2_Supporting" in text


def test_duplicate_counts_once():
    r = acmg.classify(["PM2_Supporting", "PM2"])
    assert r["points"] == 1


@pytest.mark.parametrize("bad", ["PX1", "PVS2", "PM2_Huge", "BA1_Strong", "PP3_StandAlone"])
def test_bad_codes(bad):
    with pytest.raises(ValueError):
        acmg.classify([bad])


def test_suggest_revel_and_frequency():
    s = acmg.suggest_from_data({"revel": 0.95, "consequence": ["missense_variant"], "gnomad_ac": 0, "spliceai_max": 0.0})
    codes = {x["code"] for x in s}
    assert "PP3_Strong" in codes and "PM2_Supporting" in codes and "BP4" not in codes


def test_suggest_benign_side():
    s = acmg.suggest_from_data({"revel": 0.01, "consequence": "missense_variant", "grpmax_af": 0.08, "grpmax_an": 30000})
    codes = {x["code"] for x in s}
    assert "BP4_Strong" in codes and "BA1" in codes


def test_suggest_bs1_from_max_af():
    s = acmg.suggest_from_data({"faf95": 0.001, "max_credible_af": 0.00004, "gnomad_ac": 30})
    assert [x["code"] for x in s] == ["BS1"]


def test_spliceai_not_with_truncating():
    s = acmg.suggest_from_data({"spliceai_max": 0.6, "consequence": ["stop_gained"]})
    assert all(x["code"] != "PP3" for x in s)


# --------------------------------------------------------------------- regressions
# One test per backlog id in docs/ROADMAP.md.


def _codes(res):
    return [s["code"] for s in res]


# ---- S1: BA1 fired on the gnomAD point estimate, with no filtering AF and no exception list

def test_s1_ba1_needs_filtering_af():
    """A 60x lower faf95 must win over the grpmax point estimate."""
    s = acmg.suggest({"grpmax_af": 0.06, "grpmax_an": 2500, "faf95": 0.001, "max_credible_af": 0.0001})
    assert "BA1" not in _codes(s["suggested"])
    # the calibrated criterion in the middle is what the numbers actually support
    assert "BS1" in _codes(s["suggested"])


def test_s1_ba1_from_faf95_above_five_percent():
    s = acmg.suggest({"faf95": 0.07, "grpmax_af": 0.08, "grpmax_an": 30000})
    ba1 = next(x for x in s["suggested"] if x["code"] == "BA1")
    assert "faf95" in ba1["basis"]
    # the exception list is not bundled, so the output must not claim it was applied
    joined = " ".join(ba1["caveats"])
    assert "has NOT been checked" in joined and "Ghosh" in joined
    assert "outside the BA1 exception list" not in ba1["rule"]


def test_s1_ba1_point_estimate_still_allowed_with_an_but_labelled():
    s = acmg.suggest({"grpmax_af": 0.07102, "grpmax_an": 68000})
    ba1 = next(x for x in s["suggested"] if x["code"] == "BA1")
    assert "point estimate" in ba1["basis"] and "68000" in ba1["basis"]


def test_s1_ba1_with_clinvar_plp_is_loudly_flagged():
    """HFE p.Cys282Tyr: ClinVar Pathogenic at 2 stars, grpmax AF 7.1%."""
    s = acmg.suggest({"grpmax_af": 0.07102, "grpmax_an": 68000, "revel": 0.872,
                      "consequence": ["missense_variant"],
                      "clinvar_classification": "Pathogenic/Pathogenic, low penetrance; risk factor"})
    ba1 = next(x for x in s["suggested"] if x["code"] == "BA1")
    assert ba1["requires_review"] is True
    assert "LOW-PENETRANCE" in ba1["conflict"] and "DO NOT apply BA1" in ba1["conflict"]
    assert any("DO NOT apply BA1" in c for c in s["caveats"])


def test_s1_ba1_without_clinvar_says_to_check_it():
    s = acmg.suggest({"grpmax_af": 0.08, "grpmax_an": 30000})
    ba1 = next(x for x in s["suggested"] if x["code"] == "BA1")
    assert any("no ClinVar classification was supplied" in c for c in ba1["caveats"])


# ---- S2: BP4 for a missense from SpliceAI alone

def test_s2_bp4_not_from_spliceai_alone_on_a_missense():
    s = acmg.suggest({"spliceai_max": 0.05, "consequence": ["missense_variant"]})
    assert _codes(s["suggested"]) == []
    assert any("protein functional impact has been excluded" in n for n in s["not_assessed"])


def test_s2_bp4_from_spliceai_on_a_missense_once_protein_impact_excluded():
    s = acmg.suggest({"spliceai_max": 0.05, "revel": 0.10, "consequence": ["missense_variant"]})
    # one BP4, from the protein axis, not two
    assert _codes(s["suggested"]) == ["BP4_Moderate"]


def test_s2_bp4_from_spliceai_on_an_intronic_variant_is_fine():
    s = acmg.suggest({"spliceai_max": 0.05, "consequence": ["intron_variant"]})
    assert _codes(s["suggested"]) == ["BP4"]


def test_s2_alphamissense_can_exclude_protein_impact():
    s = acmg.suggest({"spliceai_max": 0.02, "alphamissense": 0.08, "consequence": ["missense_variant"]})
    assert _codes(s["suggested"]) == ["BP4"]
    assert any("AlphaMissense" in c and "Cheng" in c for c in s["caveats"])


def test_s2_bp7_follows_walker_2023_scope_not_every_variant():
    """Walker et al. 2023 replaced the old "never BP7" stance: BP7 after BP4, synonymous outside the
    splice region; a synonymous variant in the first/last 3 exonic bases gets BP4 but no BP7."""
    s = acmg.suggest({"spliceai_max": 0.02, "consequence": ["synonymous_variant"]})
    assert _codes(s["suggested"]) == ["BP4", "BP7"]
    assert "Walker et al. 2023" in s["suggested"][1]["rule"]
    edge = acmg.suggest({"spliceai_max": 0.02, "consequence": ["synonymous_variant", "splice_region_variant"]})
    assert _codes(edge["suggested"]) == ["BP4"]
    assert any(n.startswith("BP7: not offered") and "last 3 bases" in n for n in edge["not_assessed"])


# ---- S4: the BS1 branch was unreachable

def test_s4_bs1_reachable_from_an_argument():
    """CFTR p.Arg117His inputs: BS1 once the maximum credible AF is supplied."""
    d = {"revel": 0.807, "grpmax_af": 0.00271, "grpmax_an": 1179756, "gnomad_ac": 3445,
         "faf95": 0.00263119, "consequence": ["missense_variant"], "is_missense": True}
    assert _codes(acmg.suggest(d)["suggested"]) == ["PP3_Moderate"]
    with_max = acmg.suggest(d, max_credible_af=0.002)
    assert _codes(with_max["suggested"]) == ["PP3_Moderate", "BS1"]


def test_s4_bs1_absence_is_stated_not_silent():
    s = acmg.suggest({"faf95": 0.001, "revel": 0.8, "consequence": ["missense_variant"]})
    assert any("BS1:" in n and "maximum credible AF" in n for n in s["not_assessed"])
    assert any("not evidence" in n for n in s["not_assessed"])


# ---- S5: computational evidence alone reached Benign

def test_s5_computational_only_cannot_reach_benign():
    r = acmg.classify(["BP4_VeryStrong"])
    assert r["points"] == -8
    assert r["classification"] == "Uncertain significance"
    assert r["capped_at_uncertain"] == "Benign"
    assert any("Pejaver" in w for w in r["warnings"])


def test_s5_computational_only_cannot_reach_likely_benign():
    assert acmg.classify(["BP4_Strong", "BP7"])["classification"] == "Uncertain significance"


def test_s5_one_real_benign_code_lifts_the_cap():
    r = acmg.classify(["BS1", "BP4"])
    assert r["computational_only"] is False and r["classification"] == "Likely benign"


def test_s5_note_says_the_points_rest_on_one_criterion():
    # PM3_VeryStrong is 8 points (Likely pathogenic) but satisfies no 2015 clause
    r = acmg.classify(["PM3_VeryStrong"])
    assert r["classification"] == "Likely pathogenic"
    assert r["classification_richards_2015"] == "Uncertain significance"
    assert "single criterion" in r["note"] and "PM3_VeryStrong" in r["note"]


# ---- F23: contradictory or duplicate suggestions

def test_f23_no_bp4_alongside_pp3():
    s = acmg.suggest({"revel": 0.10, "spliceai_max": 0.6, "consequence": ["missense_variant"]})
    assert _codes(s["suggested"]) == ["PP3"]
    assert any("withheld because SpliceAI" in n for n in s["not_assessed"])


def test_f23_bp4_without_a_splice_prediction_says_so():
    s = acmg.suggest({"revel": 0.10, "spliceai_max": None, "consequence": ["missense_variant"]})
    bp4 = next(x for x in s["suggested"] if x["code"].startswith("BP4"))
    assert "splice impact not assessed" in bp4["basis"]


def test_f23_only_one_pp3_is_offered():
    s = acmg.suggest({"revel": 0.95, "spliceai_max": 0.6, "consequence": ["missense_variant"]})
    assert _codes(s["suggested"]) == ["PP3_Strong"]
    assert any("PP3 counts once" in c for c in s["caveats"])


def test_f23_classify_warns_on_pp3_with_bp4():
    r = acmg.classify(["PP3_Moderate", "BP4_Moderate", "PM2_Supporting"])
    assert any("PP3 with BP4" in w for w in r["warnings"])


def test_f23_classify_warns_on_pp3_with_bp7():
    assert any("PP3 with BP7" in w for w in acmg.classify(["PP3", "BP7"])["warnings"])


def test_f23_duplicate_pp3_still_counts_once():
    r = acmg.classify(["PP3_Strong", "PP3"])
    assert r["points"] == 4 and any("given twice" in w for w in r["warnings"])


def test_f23_spliceai_uninformative_band_is_named():
    s = acmg.suggest({"spliceai_max": 0.162, "consequence": ["intron_variant"]})
    assert _codes(s["suggested"]) == []
    assert any("uninformative band" in c for c in s["caveats"])


# ---- F24: frequency codes went missing silently

def test_f24_ba1_with_unknown_an_is_reported_not_dropped():
    s = acmg.suggest({"grpmax_af": 0.30, "grpmax_an": None, "gnomad_ac": 5000,
                      "consequence": ["missense_variant"]})
    assert _codes(s["suggested"]) == []
    assert any("BA1:" in n and "not reported" in n for n in s["not_assessed"])


def test_f24_gnomad_an_key_is_accepted_as_an_alternative():
    s = acmg.suggest({"grpmax_af": 0.30, "gnomad_an": 40000})
    assert "BA1" in _codes(s["suggested"])


# ---- F25: PM2 from AC = 0 with no coverage check; none for dominant unless AC == 0

def test_f25_pm2_refused_when_the_site_was_not_callable():
    s = acmg.suggest({"gnomad_ac": 0, "grpmax_an": 0}, inheritance="AD")
    assert _codes(s["suggested"]) == []
    assert any("not callable" in n for n in s["not_assessed"])


def test_f25_pm2_refused_below_the_allele_number_floor():
    s = acmg.suggest({"gnomad_ac": 0, "grpmax_an": 300}, inheritance="AD")
    assert _codes(s["suggested"]) == []
    assert any("only 300 alleles" in n for n in s["not_assessed"])


def test_f25_pm2_from_absence_names_the_allele_number():
    s = acmg.suggest({"gnomad_ac": 0, "grpmax_an": 150000}, inheritance="AD")
    pm2 = next(x for x in s["suggested"] if x["code"] == "PM2_Supporting")
    assert "AN=150000" in pm2["basis"]


def test_f25_no_pm2_for_a_dominant_disorder_just_for_being_below_the_max_credible_af():
    """D-P0-2: this test used to lock the bug in. Below the ceiling is not rarity: a dominant
    variant present in gnomAD gets no PM2 from the maximum credible AF (ClinGen SVI 2020)."""
    s = acmg.suggest({"grpmax_af": 1.2e-6, "gnomad_ac": 1, "grpmax_an": 800000, "faf95": 1.0e-6},
                     inheritance="AD", max_credible_af=4e-5)
    assert _codes(s["suggested"]) == []
    assert any(n.startswith("PM2: not offered") and "present in gnomAD" in n for n in s["not_assessed"])
    assert any("BS1 does not apply" in c and "not evidence that the variant is rare" in c for c in s["caveats"])


# ---- F26: any code could be modified to any strength

@pytest.mark.parametrize("bad,ceiling", [("PP1_VeryStrong", "Strong"), ("PP3_VeryStrong", "Strong"),
                                         ("PM2_Strong", "Moderate"), ("PM2_VeryStrong", "Moderate")])
def test_f26_strength_ceilings_are_refused(bad, ceiling):
    with pytest.raises(ValueError) as err:
        acmg.classify([bad])
    assert f"cannot be applied above {ceiling}" in str(err.value)
    # the error names the source of the ceiling, so it can be checked
    assert any(name in str(err.value) for name in ("Biesecker", "Pejaver", "Richards"))


@pytest.mark.parametrize("ok", ["PP1_Strong", "PP3_Strong", "PM2", "PM2_Supporting", "BP4_VeryStrong",
                                "PS2_VeryStrong", "PM3_Strong"])
def test_f26_everything_at_or_below_its_ceiling_is_accepted(ok):
    assert acmg.classify([ok])["codes"][0]["code"] == ok.split("_")[0]


def test_f26_pp1_ceiling_agrees_with_the_stats_module():
    from zebra import stats

    assert stats.segregation(ad_meioses=50)["pp1_code"] == "PP1_Strong"
    acmg.classify(["PP1_Strong"])  # the strongest PP1 stats can produce is accepted


# ---- F27: the PM1 + PP3 cap was not applied to the 2015 reading

def test_f27_pm1_pp3_cap_applies_to_both_readings():
    r = acmg.classify(["PM1", "PP3_Strong"])
    assert r["points"] == 4
    assert r["pm1_pp3_cap_applied"] is True
    assert r["classification_richards_2015"] == "Uncertain significance"
    assert r["agree"] is True
    assert any("Capped in both readings" in w for w in r["warnings"])


def test_f27_cap_not_applied_when_the_sum_is_within_strong():
    r = acmg.classify(["PM1", "PP3"])
    assert r["points"] == 3 and r["pm1_pp3_cap_applied"] is False


def test_s1_ba1_withheld_when_faf95_is_low_says_why():
    s = acmg.suggest({"grpmax_af": 0.06, "grpmax_an": 2500, "faf95": 0.001})
    assert "BA1" not in _codes(s["suggested"])
    assert any("BA1 is not offered" in c and "few alleles" in c for c in s["caveats"])


# ---- S6 / F32: the s2f interpretation text claimed BP7 for every variant class and hid
# the uninformative band, the supporting cap, the run-settings/calibration mismatch, the
# max-over-tissues headline and SpliceAI/Pangolin's correlation.
#
# These assertions live here, with the ACMG threshold contract they mirror, because
# tests/test_s2f.py covers s2f's runtime and is owned elsewhere; only the interpretation
# strings in zebra/s2f.py are changed by this work.

def _s2f():
    from zebra import s2f

    return s2f


def test_s6_spliceai_reading_does_not_hand_bp7_to_every_variant():
    text = _s2f().SPLICEAI_READING
    assert "BP4/BP7" not in text, "BP7 must not be bundled with BP4 for every class"
    assert "BP7 is NOT established by this number alone" in text
    # Walker et al. 2023's scope, which replaced the conservation requirement of Richards et al. 2015
    assert "only after BP4" in text and "+7/-21" in text and "last 3 bases of the exon" in text
    assert "canonical +/-1 and +/-2" in text and "PVS1 covers that mechanism" in text


def test_f32_spliceai_reading_names_the_uninformative_band():
    text = _s2f().SPLICEAI_READING
    assert "0.1 < delta < 0.2 is UNINFORMATIVE" in text
    assert "not a lean" in text
    assert "0.162" in text, "the CFTR 3849+10kbC>T reproduction belongs in the text"


def test_f32_spliceai_reading_states_the_pp3_supporting_cap():
    text = _s2f().SPLICEAI_READING
    assert "SUPPORTING weight only" in text and "do not upgrade" in text
    assert "Walker et al. 2023" in text and "AJHG 110:1046" in text
    assert ">= 0.5 is high confidence" not in text


def test_f32_spliceai_reading_discloses_the_run_settings_mismatch():
    text = _s2f().SPLICEAI_READING
    assert "+/-500" in text and "+/-4,999" in text
    assert "mask=1" in text and "voids them" in text


def test_f32_spliceai_thresholds_match_the_acmg_module():
    s2f = _s2f()
    assert f">= {acmg.SPLICEAI_PP3}" in s2f.SPLICEAI_READING
    assert f"<= {acmg.SPLICEAI_BP4}" in s2f.SPLICEAI_READING


def test_f32_pangolin_headline_is_disclosed_as_a_max_over_tissues():
    s2f = _s2f()
    for key in ("DS_SG", "DS_SL"):
        d = s2f.PANGOLIN_DEFS[key]
        assert "heart, liver, brain, testis" in d
        assert "LARGEST" in d
    assert "maximum over tissues" in s2f.PANGOLIN_DEFS["DS_SG"]
    assert "maximum (gain) or minimum (loss) over heart, liver, brain and testis" in s2f.PANGOLIN_READING


def test_f32_pangolin_and_spliceai_are_not_called_independent():
    text = _s2f().PANGOLIN_READING
    assert "NOT independent" in text
    assert "Zeng & Li" in text and "Genome Biol 23:103" in text
    assert "largely expected" in text and "weak corroboration" in text
    assert "agreement between the two is evidence" not in text.lower()


def test_f32_pangolin_carries_no_clingen_thresholds():
    assert "no ClinGen-calibrated PP3/BP4 thresholds" in _s2f().PANGOLIN_READING


def test_f27_cap_weakens_pm1_too_when_pp3_alone_cannot_satisfy_it():
    """PM1_VeryStrong + PP3 is 9 points; the pair's ceiling is 4, so both must give."""
    r = acmg.classify(["PM1_VeryStrong", "PP3"])
    assert r["points"] == 4 and r["pm1_pp3_cap_applied"] is True
    assert r["classification"] == "Uncertain significance"
    assert r["classification_richards_2015"] == "Uncertain significance"


# ===================================================================== v0.2 (W6)
# One test per finding id in the v0.2 reviews; each was checked to fail on the unfixed code.

import argparse  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402

from zebra.core import Outcome, UsageError  # noqa: E402

FIX = os.path.join(os.path.dirname(__file__), "fixtures")


def _fx(*parts):
    with open(os.path.join(FIX, *parts)) as fh:
        return json.load(fh)


def _suggest_args(**kw):
    base = dict(prevalence=None, allelic=None, genetic=1.0, penetrance=1.0, inheritance="unknown",
                inheritance_mode=None)
    base.update(kw)
    return argparse.Namespace(**base)


# F508del, gnomAD v4 joint: AC 19,237, grpmax nfe AF 0.0149, faf95 0.01476 (tests/fixtures/gnomad/f508del_r4.json)
F508DEL = {"grpmax_af": 0.0149, "grpmax_an": 1179756, "gnomad_ac": 19237, "faf95": 0.01475724,
           "consequence": ["inframe_deletion"], "clinvar_classification": "Pathogenic"}


# ---- D-P0-1: --inheritance AR computed the maximum credible AF with the monoallelic formula

def test_d_p0_1_ar_inheritance_uses_the_biallelic_formula():
    from zebra.commands import acmg as cmd

    res = cmd._max_credible_af(_suggest_args(prevalence=0.0004, allelic=0.9, inheritance="AR"))
    assert res["max_credible_af"] == pytest.approx(0.9 * 0.0004 ** 0.5)  # 0.018, not 0.00018
    assert "biallelic" in res["model"] and "--inheritance AR" in res["formula_basis"]
    ad = cmd._max_credible_af(_suggest_args(prevalence=0.0004, allelic=0.9, inheritance="AD"))
    assert ad["max_credible_af"] == pytest.approx(0.0004 * 0.9 / 2)


def test_d_p0_1_a_formula_contradicting_the_inheritance_is_refused():
    from zebra.commands import acmg as cmd

    with pytest.raises(UsageError, match="contradicts"):
        cmd._max_credible_af(_suggest_args(prevalence=0.0004, allelic=0.9, inheritance="AR",
                                           inheritance_mode="monoallelic"))


def test_d_p0_1_x_linked_or_unknown_inheritance_must_choose_the_formula():
    from zebra.commands import acmg as cmd

    for inh in ("XLR", "XLD", "unknown"):
        with pytest.raises(UsageError, match="--inheritance-mode"):
            cmd._max_credible_af(_suggest_args(prevalence=0.0004, allelic=0.9, inheritance=inh))
    res = cmd._max_credible_af(_suggest_args(prevalence=0.0004, allelic=0.9, inheritance="XLR",
                                             inheritance_mode="biallelic"))
    assert "no X-linked formula" in res["formula_basis"]


def test_d_p0_1_bs1_against_a_clinvar_pathogenic_assertion_requires_review():
    """The monoallelic ceiling (0.00018) put F508del above it; BS1 then carried no conflict at all."""
    s = acmg.suggest(F508DEL, inheritance="AR", max_credible_af=0.00018)
    bs1 = next(x for x in s["suggested"] if x["code"] == "BS1")
    assert bs1["requires_review"] is True and "DO NOT apply BS1" in bs1["conflict"]
    assert any("DO NOT apply BS1" in c for c in s["caveats"])


# ---- D-P0-2: "faf95 <= maximum credible AF" emitted PM2_Supporting (F508del at 1.5 %)

def test_d_p0_2_f508del_below_the_disease_ceiling_gets_no_pm2():
    s = acmg.suggest(F508DEL, inheritance="AR", max_credible_af=0.018)
    assert _codes(s["suggested"]) == []
    assert any("not evidence that the variant is rare" in c for c in s["caveats"])
    assert any(n.startswith("PM2: not offered") and "0.07%" in n for n in s["not_assessed"])


def test_d_p0_2_pm2_still_given_for_absence_and_for_extremely_rare_recessive_alleles():
    assert _codes(acmg.suggest({"gnomad_ac": 0, "grpmax_an": 150000}, inheritance="AD")["suggested"]) == ["PM2_Supporting"]
    rare = acmg.suggest({"gnomad_ac": 3, "grpmax_af": 1e-5, "grpmax_an": 500000, "faf95": 4e-6},
                        inheritance="AR", max_credible_af=0.018)
    assert _codes(rare["suggested"]) == ["PM2_Supporting"]
    assert "extremely low frequency if recessive" in rare["suggested"][0]["rule"]


def test_d_p0_2_bs1_and_pm2_are_never_offered_together():
    s = acmg.suggest({"gnomad_ac": 30, "grpmax_af": 5e-4, "grpmax_an": 60000, "faf95": 3e-4},
                     inheritance="AR", max_credible_af=1e-4)
    assert _codes(s["suggested"]) == ["BS1"]
    assert any(n.startswith("PM2: not assessed — BS1 is offered") for n in s["not_assessed"])


def test_d_p0_2_a_gene_specific_pm2_ceiling_is_honoured_and_attributed():
    d = {"gnomad_ac": 2, "grpmax_af": 1.5e-5, "grpmax_an": 130000}
    assert _codes(acmg.suggest(d, inheritance="AD")["suggested"]) == []
    s = acmg.suggest(d, inheritance="AD", pm2_max_af=2e-5)
    assert _codes(s["suggested"]) == ["PM2_Supporting"]
    assert "caller" in s["suggested"][0]["basis"] and "VCEP" in s["suggested"][0]["rule"]


# ---- D-P0-3: the PM1 + PP3 cap deleted PP3 when the excess was odd

def test_d_p0_3_pp3_strong_with_pm1_supporting_counts_four_points():
    r = acmg.classify(["PP3_Strong", "PM1_Supporting"])
    assert r["points"] == 4
    assert [c["label"] for c in r["codes_after_cap"]] == ["PP3_Strong"]
    assert any("PM1 not counted" in w for w in r["warnings"])


def test_d_p0_3_the_review_example_reaches_likely_pathogenic():
    r = acmg.classify(["PP3_Strong", "PM1_Supporting", "PM2_Supporting", "PS4_Moderate"])
    assert r["points"] == 7
    assert r["classification"] == "Likely pathogenic"
    assert r["classification_richards_2015"] == "Likely pathogenic"


def test_d_p0_3_pair_never_exceeds_strong_and_keeps_the_largest_total():
    for pm1 in ("PM1_Supporting", "PM1", "PM1_Strong", "PM1_VeryStrong"):
        for pp3 in ("PP3", "PP3_Moderate", "PP3_Strong"):
            r = acmg.classify([pm1, pp3])
            raw = sum(c["points"] for c in r["codes"])
            assert r["points"] == min(raw, 4), (pm1, pp3, r["points"])


# ---- D-P2-6: the 2015 column ignored the SVI 2020 PM2 combining addition

def test_d_p2_6_svi_2020_reading_gives_pvs1_plus_pm2_supporting_likely_pathogenic():
    r = acmg.classify(["PVS1", "PM2_Supporting"])
    assert r["classification_richards_2015"] == "Uncertain significance"
    assert r["classification_richards_2015_svi_2020"] == "Likely pathogenic"
    assert "ClinGen SVI's PM2 recommendation (2020)" in r["note"]


# ---- D-P2-8: mutants that survived the offline suite

def test_d_p2_8_revel_pp3_band_edges():
    def pp3(revel):
        return [c for c in _codes(acmg.suggest({"revel": revel, "consequence": ["missense_variant"]})["suggested"])
                if c.startswith("PP3")]
    assert pp3(0.932) == ["PP3_Strong"] and pp3(0.931) == ["PP3_Moderate"]
    assert pp3(0.773) == ["PP3_Moderate"] and pp3(0.772) == ["PP3"]
    assert pp3(0.644) == ["PP3"] and pp3(0.643) == []


def test_d_p2_8_ba1_from_faf95_is_at_five_percent_exactly():
    assert "BA1" in _codes(acmg.suggest({"faf95": 0.051, "grpmax_af": 0.06, "grpmax_an": 30000})["suggested"])
    assert "BA1" not in _codes(acmg.suggest({"faf95": 0.03, "grpmax_af": 0.04, "grpmax_an": 30000})["suggested"])
    assert "BA1" not in _codes(acmg.suggest({"faf95": 0.05, "grpmax_af": 0.05, "grpmax_an": 30000})["suggested"])


def test_d_p2_8_recessive_pm2_cut_is_0_07_percent():
    def pm2(af):
        return "PM2_Supporting" in _codes(acmg.suggest({"gnomad_ac": 10, "grpmax_af": af, "grpmax_an": 100000},
                                                       inheritance="AR")["suggested"])
    assert pm2(0.00069) and not pm2(0.0008) and not pm2(0.0007)


def test_d_p2_8_likely_benign_benign_boundary_is_minus_seven():
    assert acmg.classify(["BS1", "BP1", "BP5"])["points"] == -6
    assert acmg.classify(["BS1", "BP1", "BP5"])["classification"] == "Likely benign"
    seven = acmg.classify(["BS1", "BP1", "BP5", "BP2"])
    assert seven["points"] == -7 and seven["classification"] == "Benign"


def test_d_p2_8_richards_likely_benign_needs_two_supporting():
    r = acmg.classify(["BP1", "BP5"])
    assert r["classification_richards_2015"] == "Likely benign"
    assert acmg.classify(["BP1"])["classification_richards_2015"] == "Uncertain significance"


# ---- CP1-1: Walker et al. 2023 scope for the splicing codes, and the model named

def test_cp1_1_bp7_for_deep_intronic_but_not_inside_the_splice_region():
    def run(c):
        return acmg.suggest({"spliceai_max": 0.03, "consequence": ["intron_variant"], "hgvs_c": c})
    assert _codes(run("NM_000492.4:c.3718-2477C>T")["suggested"]) == ["BP4", "BP7"]
    assert _codes(run("NM_000492.4:c.1584+7A>G")["suggested"]) == ["BP4", "BP7"]
    assert _codes(run("NM_000492.4:c.1584+6A>G")["suggested"]) == ["BP4"]
    assert _codes(run("NM_000492.4:c.1585-20A>G")["suggested"]) == ["BP4"]
    assert _codes(run("NM_000492.4:c.1585-21A>G")["suggested"]) == ["BP4", "BP7"]
    unread = run(None)
    assert _codes(unread["suggested"]) == ["BP4"]
    assert any("could not be read" in n for n in unread["not_assessed"])


def test_cp1_1_no_bp7_without_bp4_and_never_for_a_missense():
    assert "BP7" not in _codes(acmg.suggest({"spliceai_max": 0.15, "consequence": ["synonymous_variant"]})["suggested"])
    s = acmg.suggest({"spliceai_max": 0.02, "revel": 0.1, "consequence": ["missense_variant"]})
    assert "BP7" not in _codes(s["suggested"])


def test_cp1_1_the_basis_names_the_model_and_run_that_drove_the_code():
    src = "SpliceAI via the Broad SpliceAI-lookup, raw scores, +/-4999 nt"
    s = acmg.suggest({"spliceai_max": 0.31, "spliceai_source": src, "consequence": ["intron_variant"]})
    assert _codes(s["suggested"]) == ["PP3"]
    assert src in s["suggested"][0]["basis"]


def test_cp1_1_classify_accepts_bp4_with_bp7_only_in_walker_scope_and_says_so():
    r = acmg.classify(["BP4", "BP7", "BS2"])
    assert any("Walker et al. 2023" in w and "+7/-21" in w for w in r["warnings"])
    assert not any("cannot both apply" in w and w.startswith("BP4 with BP7") for w in r["warnings"])


# ---- E-5 / CP1-4: an mtDNA variant never gets nuclear frequency codes

def test_e_5_mitochondrial_inputs_get_no_nuclear_frequency_codes():
    s = acmg.suggest({"mitochondrial": True, "gnomad_ac": None, "mt_af_hom": 0.0, "mt_af_het": 1e-4},
                     inheritance="AD", max_credible_af=1e-3)
    assert not any(c in ("PM2_Supporting", "BS1", "BA1") for c in _codes(s["suggested"]))
    assert any("McCormick et al. 2020" in n for n in s["not_assessed"])


# ---- CP1-3: PVS1 / PS1 / PM5 inputs are structure, from real exon structures

def _scn1a_structure():
    d = _fx("ensembl", "transcript_ENST00000674923_scn1a.json")
    return acmg.transcript_exons(d["Exon"], d["strand"]), d["strand"], d["Translation"]


def test_cp1_3_scn1a_r712x_is_predicted_to_undergo_nmd():
    exons, strand, tl = _scn1a_structure()
    cds_start = tl["start"] if strand == 1 else tl["end"]
    n = acmg.nmd_inputs(exons, strand, cds_start, 712, tl["length"])
    assert n["ptc_exon"] == "15/29"  # VEP: exon 15/29 (tests/fixtures/vep/scn1a_r712x_grch38.json)
    assert n["nmd_predicted"] is True
    assert n["fraction_of_protein_after_ptc"] == pytest.approx((2009 - 712 + 1) / 2009, abs=1e-3)


def test_cp1_3_last_exon_and_last_50_nt_of_the_penultimate_exon_escape_nmd():
    # synthetic plus-strand transcript: exons of 100, 100, 100 nt, translation starting at base 1
    exons = acmg.transcript_exons([{"start": 1, "end": 100}, {"start": 201, "end": 300},
                                   {"start": 401, "end": 500}], 1)
    # last junction is after cDNA 200; codon 51 starts at cDNA 151 (50 nt from the end), codon 50 at 148
    assert acmg.nmd_inputs(exons, 1, 1, 51, 100)["nmd_predicted"] is False
    assert acmg.nmd_inputs(exons, 1, 1, 50, 100)["nmd_predicted"] is True
    assert acmg.nmd_inputs(exons, 1, 1, 70, 100)["position"] == "in the last exon"
    single = acmg.transcript_exons([{"start": 1, "end": 900}], 1)
    assert acmg.nmd_inputs(single, 1, 1, 10, 299)["nmd_predicted"] is False


def test_cp1_3_premature_stop_is_read_from_the_protein_hgvs():
    assert acmg.ptc_codon_from_hgvsp("ENSP1:p.Arg712Ter", ["stop_gained"])[:2] == (712, 712)
    assert acmg.ptc_codon_from_hgvsp("p.Glu123GlyfsTer45", ["frameshift_variant"])[:2] == (123, 167)
    assert acmg.ptc_codon_from_hgvsp("p.Gly551fs", ["frameshift_variant"])[:2] == (551, None)


def test_cp1_3_codon_records_split_into_ps1_and_pm5_inputs():
    """CFTR codon 551 (live ClinVar capture): G551S (Pathogenic, 3 stars) and G551V (Likely pathogenic)
    are PM5 inputs for G551D; G551D itself is excluded; the frameshifts and the synonymous change are not
    missense at the residue."""
    from zebra.sources import clinvar

    res = _fx("clinvar", "esummary_cftr_codon551.json")["result"]
    recs = [clinvar.parse_summary(res[u]) for u in res["uids"]]
    out = acmg.codon_records(recs, "p.Gly551Asp", exclude_vcv="VCV000007120", gene="CFTR")
    assert out["ps1_inputs"] == []
    pm5 = {r["vcv"]: r for r in out["pm5_inputs"]}
    assert "VCV000007142" in pm5 and pm5["VCV000007142"]["stars"] == 3
    assert "VCV001333303" in pm5
    assert "VCV000007120" not in pm5
    as_ps1 = acmg.codon_records(recs, "p.Gly551Ser", exclude_vcv="VCV000007142", gene="CFTR")
    assert [r["vcv"] for r in as_ps1["ps1_inputs"]] == []  # G551S is the variant itself here
    assert {r["vcv"] for r in as_ps1["pm5_inputs"]} >= {"VCV000007120", "VCV001333303"}
    # a different nucleotide change for the same amino-acid change is a PS1 input...
    other = acmg.codon_records(recs, "p.Gly551Asp", exclude_vcv=None, gene="CFTR", exclude_c="NM_000492.4:c.1652G>C")
    assert [r["vcv"] for r in other["ps1_inputs"]] == ["VCV000007120"]
    assert "never" not in other["rule"] and "the skill decides" in other["rule"]
    # ...but the variant itself never is, even when the card found no ClinVar match for it
    # (reviewer: G551D was listed as its own PS1 input when exclude_vcv was missing)
    for kw in ({"exclude_c": "NM_000492.4:c.1652G>A"}, {"exclude_spdi": ":117587805:G:A"}):
        own = acmg.codon_records(recs, "p.Gly551Asp", exclude_vcv=None, gene="CFTR", **kw)
        assert own["ps1_inputs"] == [], kw


def test_cp1_3_codon_records_filters_non_missense_other_genes_and_conflicting():
    from zebra.sources import clinvar

    res = _fx("clinvar", "esummary_cftr_codon551.json")["result"]
    recs = [clinvar.parse_summary(res[u]) for u in res["uids"]]
    out = acmg.codon_records(recs, "p.Gly551Glu", gene="CFTR")
    listed = {r["vcv"] for r in out["ps1_inputs"] + out["pm5_inputs"]}
    assert "VCV004818480" not in listed and "VCV000053320" not in listed  # frameshifts at the codon
    assert "VCV003491706" not in listed  # the synonymous p.Gly551=
    assert acmg.codon_records(recs, "p.Gly551Glu", gene="OTHER")["pm5_inputs"] == []
    fake = [dict(recs[0], vcv="VCV9", title="NM_000492.4(CFTR):c.1652G>C (p.Gly551Ala)",
                 classification="Conflicting classifications of pathogenicity", compound=False),
            dict(recs[0], vcv="VCV8", title="NM_000492.4(CFTR):c.1652G>T (p.Gly551Val)",
                 classification="not pathogenic", compound=False)]
    out = acmg.codon_records(fake, "p.Gly551Glu", gene="CFTR")
    assert out["pm5_inputs"] == [] and {o["vcv"] for o in out["other_missense_at_codon"]} == {"VCV8", "VCV9"}
    with pytest.raises(ValueError):
        acmg.codon_records(recs, "p.Arg117Ter")


# ---- CP1-1 / CP1-18: `acmg suggest` calls the splice models and says which one drove the codes

def _card(terms, vcf=("7", 117639961, "C", "T"), assembly="GRCh38", spliceai=None, gene="CFTR", hgvs_c=None,
          grch38=None):
    data = {"consequence": terms, "spliceai_max": spliceai, "is_missense": any("missense" in t for t in terms),
            "hgvs_c": hgvs_c, "spliceai_source": "SpliceAI precomputed scores served by Ensembl VEP" if spliceai is not None else None}
    r = {"variant": "x", "assembly": assembly, "gene": gene, "transcript": {"ensembl": "ENST00000003084.11"},
         "vcf": {"chrom": vcf[0], "pos": vcf[1], "ref": vcf[2], "alt": vcf[3], "id": "-".join(map(str, vcf))},
         "consequence": {"terms": terms}, "acmg_inputs": data}
    if grch38:
        r["grch38"] = grch38
    return r, dict(data)


def _fake_predict(spliceai_rows, pangolin_rows=None, calls=None, raises=None):
    def predict(variant, assembly="GRCh38", models=None, distance=500, mask=0, **kw):
        if calls is not None:
            calls.append({"variant": variant, "assembly": assembly, "models": models, "distance": distance, "mask": mask})
        if raises:
            raise raises
        models_out = [{"model": "spliceai", "status": "ran", "transcripts": spliceai_rows}]
        if pangolin_rows is not None:
            models_out.append({"model": "pangolin", "status": "ran", "transcripts": pangolin_rows})
        return Outcome({"models": models_out}, sources=[{"db": "spliceai via Broad SpliceAI-lookup"}])
    return predict


SPLICE_ROW = [{"gene": "CFTR", "transcript": "ENST00000003084.11", "refseq": "NM_000492.4",
               "max_delta": {"score": "DS_DG", "delta": 0.03, "position": 117639963}}]


def test_cp1_1_splice_relevant_variant_is_scored_at_walker_settings(monkeypatch):
    from zebra import s2f
    from zebra.commands import acmg as cmd

    calls = []
    monkeypatch.setattr(s2f, "predict", _fake_predict(SPLICE_ROW, [{"gene": "CFTR", "refseq": "NM_000492.4",
                                                                     "max_delta": {"score": "DS_SG", "delta": 0.05}}],
                                                      calls=calls))
    r, data = _card(["intron_variant"], hgvs_c="NM_000492.4:c.3718-2477C>T")
    out = cmd.splice_evidence(r, data)
    assert calls == [{"variant": "7-117639961-C-T", "assembly": "GRCh38", "models": ["spliceai", "pangolin"],
                      "distance": 4999, "mask": 0}]
    assert out.result["used"] == "SpliceAI-lookup" and data["spliceai_max"] == pytest.approx(0.03)
    assert "+/-4999 nt" in data["spliceai_source"] and "DS_DG" in data["spliceai_source"]
    assert "not used for a code" in out.result["pangolin"]["note"]
    s = acmg.suggest(data)
    assert _codes(s["suggested"]) == ["BP4", "BP7"]
    assert "SpliceAI-lookup" in s["suggested"][0]["basis"]


def test_cp1_1_an_indel_without_a_vep_score_is_scored_and_a_covered_missense_is_not(monkeypatch):
    from zebra import s2f
    from zebra.commands import acmg as cmd

    calls = []
    monkeypatch.setattr(s2f, "predict", _fake_predict(SPLICE_ROW, calls=calls))
    r, data = _card(["inframe_deletion"], vcf=("7", 117559590, "ATCT", "A"), spliceai=None)
    assert cmd.splice_relevance(r, data) == "an indel with no precomputed SpliceAI score from VEP"
    cmd.splice_evidence(r, data)
    assert len(calls) == 1
    r, data = _card(["missense_variant"], vcf=("7", 117587806, "G", "A"), spliceai=0.0)
    assert cmd.splice_relevance(r, data) is None
    cmd.splice_evidence(r, data)
    assert len(calls) == 1
    r, data = _card(["stop_gained"], spliceai=None)
    assert cmd.splice_relevance(r, data) is None  # PVS1 covers a null variant's mechanism


def test_cp1_1_grch37_input_is_scored_on_the_mapped_grch38_coordinates(monkeypatch):
    from zebra import s2f
    from zebra.commands import acmg as cmd

    calls = []
    monkeypatch.setattr(s2f, "predict", _fake_predict(SPLICE_ROW, calls=calls))
    r, data = _card(["synonymous_variant"], vcf=("7", 117279915, "C", "T"), assembly="GRCh37",
                    grch38={"vcf": {"id": "7-117639961-C-T"}})
    cmd.splice_evidence(r, data)
    assert (calls[0]["variant"], calls[0]["assembly"]) == ("7-117639961-C-T", "GRCh38")


def test_cp1_18_unreachable_lookup_falls_back_to_vep_and_says_so(monkeypatch):
    from zebra import s2f
    from zebra.commands import acmg as cmd
    from zebra.http import SourceError

    monkeypatch.setattr(s2f, "predict", _fake_predict(SPLICE_ROW, raises=SourceError(
        "spliceai (SpliceAI-lookup)", "https://spliceai-38-xwkwwwxdwq-uc.a.run.app", None, "network error: timed out")))
    r, data = _card(["synonymous_variant"], spliceai=0.01)
    out = cmd.splice_evidence(r, data)
    assert out.result["used"] == "VEP precomputed SpliceAI"
    assert data["spliceai_max"] == 0.01 and "Ensembl VEP" in data["spliceai_source"]
    assert any("SpliceAI-lookup unavailable" in w and "timed out" in w for w in out.warnings)
    assert "Ensembl VEP" in acmg.suggest(data)["suggested"][0]["basis"]
    # and with no VEP score either, the axis is "not assessed", never "no effect"
    r, data = _card(["intron_variant"], spliceai=None)
    out = cmd.splice_evidence(r, data)
    assert out.result["used"] is None and data["spliceai_max"] is None
    assert "not assessed" in out.result["note"]
    assert not any(c.startswith("BP") for c in _codes(acmg.suggest(data)["suggested"]))


def test_cp1_1_suggest_command_wires_the_splice_and_judgement_inputs(monkeypatch):
    from zebra import s2f
    from zebra.commands import acmg as cmd
    from zebra.sources import variant as V

    r, _ = _card(["intron_variant"], hgvs_c="NM_000492.4:c.3718-2477C>T")
    r["acmg_inputs"].update({"gnomad_ac": 0, "site_covered": True, "coverage_text": "genome 97% of samples at >=20x"})
    seen = {}

    def card(variant, assembly="GRCh38", heteroplasmy=None, contracts=True):
        seen["contracts"] = contracts
        return Outcome(r, sources=[{"db": "VEP"}])

    monkeypatch.setattr(V, "card", card)
    monkeypatch.setattr(V, "judgement_inputs", lambda res: Outcome({"note": "no PVS1 or PS1/PM5 inputs"}))
    monkeypatch.setattr(s2f, "predict", _fake_predict([dict(SPLICE_ROW[0], max_delta={"score": "DS_DG", "delta": 0.33,
                                                                                         "position": 117639963})]))
    args = argparse.Namespace(variant="7-117639961-C-T", assembly="GRCh38", inheritance="AR", prevalence=None,
                              allelic=None, genetic=1.0, penetrance=1.0, inheritance_mode=None, pm2_max_af=None,
                              heteroplasmy=None, no_splice_lookup=False)
    out = cmd._suggest(args)
    assert _codes(out.result["suggested"]) == ["PP3", "PM2_Supporting"]
    assert "SpliceAI-lookup" in out.result["suggested"][0]["basis"]
    assert "splicing score used: SpliceAI via the Broad SpliceAI-lookup" in out.text
    assert seen["contracts"] is False  # the Chinese-cohort and MaveDB lookups feed no code: not run here


@pytest.mark.live
def test_live_d_p0_1_and_2_f508del_with_ar_prevalence_gets_neither_bs1_nor_pm2(capsys):
    from zebra import cli

    code = cli.main(["--json", "acmg", "suggest", "NM_000492.4:c.1521_1523del", "--inheritance", "AR",
                     "--prevalence", "0.0004", "--allelic", "0.9", "--no-splice-lookup"])
    env = json.loads(capsys.readouterr().out)
    assert code == 0, env
    assert env["result"]["max_credible_af"]["max_credible_af"] == pytest.approx(0.018)
    assert _codes(env["result"]["suggested"]) == []


@pytest.mark.live
def test_live_cp1_1_deep_intronic_cftr_is_scored_by_the_lookup(capsys):
    from zebra import cli

    code = cli.main(["--json", "acmg", "suggest", "7-117639961-C-T", "--inheritance", "AR"])
    env = json.loads(capsys.readouterr().out)
    assert code == 0
    sp = env["result"]["splicing"]
    assert sp["splice_relevant"].startswith("splice-relevant class")
    # either the lookup answered at the calibrated settings, or the fallback is said
    assert sp["used"] == "SpliceAI-lookup" or "SpliceAI-lookup unavailable" in (sp.get("note") or "")


@pytest.mark.live
def test_live_cp1_3_scn1a_r712x_pvs1_inputs(capsys):
    from zebra import cli

    code = cli.main(["--json", "acmg", "suggest", "NM_001165963.4:c.2134C>T", "--inheritance", "AD", "--no-splice-lookup"])
    env = json.loads(capsys.readouterr().out)
    assert code == 0
    pv = env["result"]["pvs1_inputs"]
    assert pv["nmd"]["nmd_predicted"] is True and pv["nmd"]["ptc_exon"] == "15/29"
    assert pv["lof_mechanism"]["clingen_hi_score"] == "3"


# ---- adversarial review (round 1) findings on the code above

def test_review_a_p0_1_a_caller_pm2_ceiling_is_not_overridden_by_the_recessive_convention():
    d = {"gnomad_ac": 40, "grpmax_af": 5e-4, "grpmax_an": 80000}
    assert _codes(acmg.suggest(d, inheritance="AR")["suggested"]) == ["PM2_Supporting"]  # generic 0.07%
    s = acmg.suggest(d, inheritance="AR", pm2_max_af=2e-4)
    assert _codes(s["suggested"]) == []
    assert any("above the gene-specific ceiling" in n for n in s["not_assessed"])
    assert _codes(acmg.suggest(d, inheritance="AR", pm2_max_af=5e-4)["suggested"]) == ["PM2_Supporting"]  # <= is in


def test_review_a_p1_2_frequency_flags_refuse_nan_inf_and_out_of_range(capsys):
    from zebra import cli

    for flag, value in (("--pm2-max-af", "inf"), ("--pm2-max-af", "nan"), ("--pm2-max-af", "2"),
                        ("--prevalence", "0"), ("--allelic", "nan")):
        assert cli.main(["--json", "acmg", "suggest", "7-117559590-ATCT-A", flag, value]) == 2, (flag, value)
        capsys.readouterr()


def test_review_a_p1_2_allelic_without_prevalence_is_refused():
    from zebra.commands import acmg as cmd

    with pytest.raises(UsageError, match="--allelic needs --prevalence"):
        cmd._max_credible_af(_suggest_args(allelic=0.9))
    with pytest.raises(UsageError, match="contradicts"):
        cmd._max_credible_af(_suggest_args(inheritance="AD", inheritance_mode="biallelic"))


def test_review_a_p1_3_bp7_is_judged_by_the_exon_side_end_of_an_intronic_range():
    for c in ("NM_1.1:c.124-30_124-5del", "NM_1.1:c.124-25_124-3del", "NM_1.1:c.123+3_123+40del"):
        s = acmg.suggest({"spliceai_max": 0.02, "consequence": ["intron_variant"], "hgvs_c": c})
        assert _codes(s["suggested"]) == ["BP4"], c
    s = acmg.suggest({"spliceai_max": 0.02, "consequence": ["intron_variant"], "hgvs_c": "NM_1.1:c.124-60_124-22del"})
    assert _codes(s["suggested"]) == ["BP4", "BP7"]
    # an intronic term on a variant that also touches the exon is not "intronic" for BP7
    s = acmg.suggest({"spliceai_max": 0.02, "consequence": ["intron_variant", "coding_sequence_variant"],
                      "hgvs_c": "NM_1.1:c.123+30A>G"})
    assert "BP7" not in _codes(s["suggested"])


def _lookup_row(gene, delta, refseq, tid, priority="MS"):
    return {"gene": gene, "transcript": tid, "refseq": refseq, "priority": priority,
            "max_delta": {"score": "DS_AG", "delta": delta, "position": 1}}


def test_review_a_p1_4_a_score_for_another_gene_never_drives_the_codes(monkeypatch):
    from zebra import s2f
    from zebra.commands import acmg as cmd

    monkeypatch.setattr(s2f, "predict", _fake_predict([_lookup_row("CFTR-AS2", 0.0, None, "ENST00000456270.1", "C")]))
    r, data = _card(["intron_variant"], spliceai=0.5, hgvs_c="NM_000492.4:c.3718-2477C>T")
    out = cmd.splice_evidence(r, data)
    assert out.result["used"] == "VEP precomputed SpliceAI" and data["spliceai_max"] == 0.5
    assert "scored no transcript of CFTR" in out.result["note"]
    assert _codes(acmg.suggest(data)["suggested"]) == ["PP3"]


def test_review_a_p1_4_the_cards_own_transcript_is_preferred(monkeypatch):
    from zebra import s2f
    from zebra.commands import acmg as cmd

    rows = [_lookup_row("CFTR", 0.05, "NM_999.1", "ENST00000999999.1", "MS"),
            _lookup_row("CFTR", 0.35, "NM_000492.4", "ENST00000003084.11", "C")]
    monkeypatch.setattr(s2f, "predict", _fake_predict(rows))
    r, data = _card(["intron_variant"], hgvs_c="NM_000492.4:c.3718-2477C>T")
    cmd.splice_evidence(r, data)
    assert data["spliceai_max"] == pytest.approx(0.35) and "NM_000492.4" in data["spliceai_source"]


def test_review_a_p1_4_benign_splicing_is_withheld_when_another_gene_transcript_predicts_an_effect(monkeypatch):
    from zebra import s2f
    from zebra.commands import acmg as cmd

    rows = [_lookup_row("CFTR", 0.02, "NM_000492.4", "ENST00000003084.11"),
            _lookup_row("CFTR", 0.6, None, "ENST00000426809.5", "N")]
    monkeypatch.setattr(s2f, "predict", _fake_predict(rows))
    r, data = _card(["intron_variant"], hgvs_c="NM_000492.4:c.3718-2477C>T")
    cmd.splice_evidence(r, data)
    s = acmg.suggest(data)
    assert _codes(s["suggested"]) == []
    assert any("withheld" in n and "ENST00000426809.5" in n for n in s["not_assessed"])


def test_review_a_p1_5_protein_bp4_is_withheld_while_splicing_is_uninformative():
    s = acmg.suggest({"revel": 0.1, "spliceai_max": 0.15, "consequence": ["missense_variant"]})
    assert _codes(s["suggested"]) == []
    assert any("uninformative band" in n and n.startswith("BP4 (protein axis") for n in s["not_assessed"])
    assert _codes(acmg.suggest({"revel": 0.1, "spliceai_max": 0.1, "consequence": ["missense_variant"]})["suggested"]) == ["BP4_Moderate"]


def test_review_a_p1_7_x_linked_recessive_does_not_get_the_autosomal_recessive_pm2_convention():
    s = acmg.suggest({"gnomad_ac": 40, "grpmax_af": 5e-4, "grpmax_an": 80000}, inheritance="XLR")
    assert _codes(s["suggested"]) == []
    assert any("hemizygous male" in n for n in s["not_assessed"])
    assert _codes(acmg.suggest({"gnomad_ac": 0, "grpmax_an": 80000}, inheritance="XLR")["suggested"]) == ["PM2_Supporting"]


def test_review_a_p2_inputs_are_coerced_or_refused_never_compared_as_strings():
    assert "BA1" in _codes(acmg.suggest({"grpmax_af": "0.06", "grpmax_an": "5000"})["suggested"])
    with pytest.raises(ValueError):
        acmg.suggest({"grpmax_af": 0.06, "grpmax_an": "many"})


def test_review_a_p2_clinvar_not_pathogenic_and_conflicting_are_not_plp():
    assert acmg.is_plp("Pathogenic/Likely pathogenic") and acmg.is_plp("Likely pathogenic")
    assert not acmg.is_plp("not pathogenic") and not acmg.is_plp("probable-non-pathogenic")
    assert not acmg.is_plp("Conflicting classifications of pathogenicity")


def test_review_a_p2_recessive_pm2_needs_enough_alleles_and_silence_is_named():
    s = acmg.suggest({"gnomad_ac": 1, "grpmax_af": 0.00069, "grpmax_an": 1400}, inheritance="AR")
    assert _codes(s["suggested"]) == [] and any("only 1400 alleles" in n for n in s["not_assessed"])
    s = acmg.suggest({}, inheritance="AR")
    assert any(n.startswith("PM2: not assessed — no gnomAD") for n in s["not_assessed"])


def test_review_a_p2_premature_stop_guards():
    assert acmg.ptc_codon_from_hgvsp("p.Gly117AlafsTer0", ["frameshift_variant"])[1] is None
    exons = acmg.transcript_exons([{"start": 1, "end": 300}, {"start": 401, "end": 500}], 1)
    with pytest.raises(ValueError):
        acmg.nmd_inputs(exons, 1, 1, 0, 100)


def test_review_a_p2_stop_lost_is_not_a_pvs1_class_and_string_consequence_is_read():
    from zebra.commands import acmg as cmd

    assert "stop_lost" not in acmg.NULL_CONSEQUENCES
    r, data = _card(["intron_variant"])
    data["consequence"] = "intron_variant"
    assert cmd.splice_relevance(r, data) is not None


def test_review_a_p2_alphamissense_does_not_overrule_an_indeterminate_revel():
    s = acmg.suggest({"revel": 0.5, "alphamissense": 0.1, "spliceai_max": 0.02, "consequence": ["missense_variant"]})
    assert "BP4" not in _codes(s["suggested"])
    s = acmg.suggest({"alphamissense": 0.1, "spliceai_max": 0.02, "consequence": ["missense_variant"]})
    assert _codes(s["suggested"]) == ["BP4"]


def test_review_a_gaps_spliceai_exact_boundaries():
    def codes(sp):
        return _codes(acmg.suggest({"spliceai_max": sp, "consequence": ["intron_variant"],
                                    "hgvs_c": "NM_1.1:c.100+50A>G"})["suggested"])
    assert codes(0.2) == ["PP3"] and codes(0.1999) == [] and codes(0.1) == ["BP4", "BP7"] and codes(0.1001) == []


def test_review_a_gaps_nmd_boundary_with_a_translation_start_inside_exon_1():
    # exon 1 = cDNA 1-100 (start codon at cDNA 11), exon 2 = 101-200, exon 3 = 201-300; last junction after 200.
    exons = acmg.transcript_exons([{"start": 1, "end": 100}, {"start": 201, "end": 300}, {"start": 401, "end": 500}], 1)
    # codon k starts at cDNA 11 + 3(k-1): k=47 -> 149 (52 nt before the junction), k=48 -> 152 (49 nt)
    assert acmg.nmd_inputs(exons, 1, 11, 47, 90)["nmd_predicted"] is True
    assert acmg.nmd_inputs(exons, 1, 11, 48, 90)["nmd_predicted"] is False
    # cDNA 151 is exactly the 50th nt from the end of the penultimate exon: escapes
    exons2 = acmg.transcript_exons([{"start": 1, "end": 100}, {"start": 201, "end": 300}, {"start": 401, "end": 500}], 1)
    assert acmg.nmd_inputs(exons2, 1, 10, 48, 90)["ptc_cdna"] == 151
    assert acmg.nmd_inputs(exons2, 1, 10, 48, 90)["nmd_predicted"] is False
    assert acmg.nmd_inputs(exons2, 1, 9, 48, 90)["nmd_predicted"] is True  # cDNA 150: 51 nt before


def test_review_a_gaps_cap_keeps_pm1_and_weakens_pp3_on_a_tie():
    r = acmg.classify(["PP3_Strong", "PM1"])
    assert [c["label"] for c in r["codes_after_cap"]] == ["PM1", "PP3_Moderate"]


def test_review_a_gaps_svi_rule_needs_a_supporting_code_and_respects_the_computational_cap():
    assert acmg.classify(["PVS1"])["classification_richards_2015_svi_2020"] == "Uncertain significance"
    r = acmg.classify(["BP4_Strong", "BP7"])
    assert r["classification_richards_2015_svi_2020"] == "Uncertain significance"


def test_review_a_gaps_absence_at_a_poorly_covered_site_is_not_pm2():
    s = acmg.suggest({"gnomad_ac": 0, "site_covered": False, "coverage_text": "genome 41% of samples at >=20x"},
                     inheritance="AD")
    assert _codes(s["suggested"]) == [] and any("not adequately covered" in n for n in s["not_assessed"])
