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


def test_s2_bp7_is_never_claimed_in_the_rule_text():
    s = acmg.suggest({"spliceai_max": 0.02, "consequence": ["synonymous_variant"]})
    assert all("BP7" not in x["rule"] for x in s["suggested"])
    assert any("BP7 is not offered" in c for c in s["caveats"])


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


def test_f25_pm2_for_a_dominant_disorder_below_the_max_credible_af():
    s = acmg.suggest({"grpmax_af": 1.2e-6, "gnomad_ac": 1, "grpmax_an": 800000, "faf95": 1.0e-6},
                     inheritance="AD", max_credible_af=4e-5)
    assert _codes(s["suggested"]) == ["PM2_Supporting"]
    assert "maximum credible AF" in s["suggested"][0]["rule"]


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
    assert "BP7 is NOT established by this number" in text
    assert "not to be highly conserved" in text and "Richards et al. 2015" in text
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
