import math

import pytest

from zebra import stats as S


def test_fisher_matches_reference():
    # reference values from scipy.stats.fisher_exact
    assert S.fisher_exact(3, 1, 1, 3)["p_two_sided"] == pytest.approx(0.485714, abs=1e-6)
    assert S.fisher_exact(8, 2, 1, 5)["p_two_sided"] == pytest.approx(0.034965, abs=1e-6)
    assert S.fisher_exact(8, 2, 1, 5)["p_greater"] == pytest.approx(0.024476, abs=1e-6)
    assert S.fisher_exact(0, 10, 5, 20)["p_two_sided"] == pytest.approx(0.291505, abs=1e-6)


def test_distributions():
    assert S.poisson_sf(5, 1.2) == pytest.approx(0.0077457883, rel=1e-8)
    assert S.binom_sf(3, 100, 0.01) == pytest.approx(0.0793732023, rel=1e-8)
    assert S.t_sf_two_sided(2.3, 7) == pytest.approx(0.0549910952, rel=1e-8)
    assert S.chi2_sf_1df(3.84) == pytest.approx(0.0500435212, rel=1e-8)
    assert S._t_quantile(0.975, 5) == pytest.approx(2.5705818, rel=1e-6)


def test_whiffin_hcm_example():
    # Whiffin 2017: HCM prevalence 1/500, allelic 0.02, penetrance 0.5 -> 4e-5
    r = S.max_credible_af(1 / 500, 0.02, 1.0, 0.5, "monoallelic")
    assert r["max_credible_af"] == pytest.approx(4e-5)


def test_biallelic_max_af():
    r = S.max_credible_af(1 / 10000, 0.5, 1.0, 1.0, "biallelic")
    assert r["max_credible_af"] == pytest.approx(0.5 * 0.01)


def test_segregation_strengths():
    """Strengths now come from Biesecker et al. 2024 Table 3 points, not the Tavtigian odds floor.

    Reconciled under F28: the old route read the Jarvik & Browning likelihood ratio against
    the Tavtigian 2018 odds (supporting 2.08), which is 1-2 points short of ClinGen at every
    count because Biesecker 2024 shifts the supporting odds path to 2.0.
    """
    assert S.segregation(ad_meioses=1)["pp1_code"] == "PP1"
    assert S.segregation(ad_meioses=2)["pp1_code"] == "PP1_Moderate"
    assert S.segregation(ad_meioses=4)["pp1_code"] == "PP1_Strong"
    assert S.segregation(ad_meioses=10)["pp1_code"] == "PP1_Strong"
    assert S.segregation(ar_affected_sibs=1)["pp1_code"] == "PP1_Moderate"
    assert S.segregation(ar_affected_sibs=2)["pp1_code"] == "PP1_Strong"


def test_carrier():
    r = S.recessive_from_prevalence(1 / 2500)
    assert r["pathogenic_allele_freq"] == pytest.approx(0.02)
    assert r["carrier_freq"] == pytest.approx(2 * 0.98 * 0.02)


def test_xlinked_bayes():
    r = S.xlinked_carrier_posterior(0.5, 3)
    assert r["posterior_carrier"] == pytest.approx(1 / 9)


def test_km_and_logrank():
    curve = S.kaplan_meier([1, 2, 2, 3, 4, 5], [1, 1, 0, 1, 0, 1])
    assert [round(r["survival"], 4) for r in curve] == [0.8333, 0.6667, 0.4444, 0.0]
    lr = S.logrank([1, 2, 3, 4], [1, 1, 1, 1], [5, 6, 7, 8], [1, 1, 1, 1])
    assert lr["p"] < 0.05


def test_nof1():
    a = S.nof1_analyze([5.1, 6.2, 5.8, 6.9], [4.0, 5.5, 5.9, 5.2])
    assert a["p_two_sided"] == pytest.approx(0.10976473, rel=1e-6)
    d = S.nof1_pairs_needed(1.0, 1.0)
    assert d["pairs"] == math.ceil((1.959964 + 0.841621) ** 2)


def test_denovo():
    r = S.denovo_enrichment(3, 1000, 1e-5)
    assert r["expected"] == pytest.approx(0.02)
    assert r["p"] < 2e-6


# --------------------------------------------------------------------- regressions
# One test per backlog id in docs/ROADMAP.md.

import zebra.commands.stats as CS
from zebra.core import UsageError


class _Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


# ---- E10: stats km hung on a NaN time or a non-0/1 event code, and miscounted ties

def test_e10_km_rejects_non_finite_times():
    with pytest.raises(ValueError) as err:
        S.kaplan_meier([1, 2, float("nan"), 3], [1, 1, 0, 1])
    assert "finite" in str(err.value) and "row 2" in str(err.value)
    with pytest.raises(ValueError):
        S.kaplan_meier([1, float("inf")], [1, 0])


def test_e10_km_rejects_unknown_event_codes():
    with pytest.raises(ValueError) as err:
        S.kaplan_meier([1, 2, 3, 4], [1, 2, 1, 0])
    assert "event code 2" in str(err.value) and "1/2" in str(err.value)


def test_e10_km_accepts_r_style_coding_explicitly():
    """1 = censored, 2 = event (R survival). The same rows read differently under each coding."""
    rows_t, rows_e = [1, 2, 3, 4], [2, 1, 2, 1]
    curve = S.kaplan_meier(rows_t, rows_e, event_coding="1/2")
    assert [round(r["survival"], 4) for r in curve] == [0.75, 0.375]
    assert [r["time"] for r in curve] == [1.0, 3.0]
    with pytest.raises(ValueError):
        S.kaplan_meier(rows_t, rows_e, event_coding="0/1")


def test_e10_km_counts_tied_events_once():
    """Four events, two of them tied at t=2: one row per distinct time, each counted once."""
    curve = S.kaplan_meier([1, 2, 2, 3, 4], [1, 1, 1, 1, 0])
    assert [r["time"] for r in curve] == [1.0, 2.0, 3.0]
    assert [r["events"] for r in curve] == [1, 2, 1]
    assert [round(r["survival"], 4) for r in curve] == [0.8, 0.4, 0.2]
    assert sum(r["events"] for r in curve) == 4


def test_e10_km_terminates_on_every_tie_pattern():
    # the old loop advanced by the rows it recognised, so any unrecognised row looped forever
    for events in ([1, 1, 1], [0, 0, 1], [1, 0, 1], [0, 1, 0]):
        assert isinstance(S.kaplan_meier([2, 2, 2], events), list)


def test_e10_logrank_validates_too():
    with pytest.raises(ValueError):
        S.logrank([1, 2], [1, 2], [3, 4], [1, 0])
    lr = S.logrank([1, 2, 3, 4], [2, 2, 2, 2], [5, 6, 7, 8], [2, 2, 2, 2], event_coding="1/2")
    assert lr["p"] < 0.05


def test_e10_cli_km_turns_a_bad_event_code_into_a_usage_error(tmp_path):
    csv_path = tmp_path / "km.csv"
    csv_path.write_text("time,event\n1,1\n2,2\n3,1\n4,0\n", encoding="utf-8")
    args = _Args(csv=str(csv_path), time_col="time", event_col="event", group_col=None, event_coding="0/1")
    with pytest.raises(UsageError) as err:
        CS._run(CS._km)(args)
    assert "event code 2" in str(err.value)
    # the same file read under 1/2 is still rejected, because 0 is not a code there either
    args.event_coding = "1/2"
    with pytest.raises(UsageError) as err2:
        CS._run(CS._km)(args)
    assert "event code 0 is not valid under coding '1/2'" in str(err2.value)
    # a genuine R-style file goes through, and the event count follows the coding
    r_path = tmp_path / "km_r.csv"
    r_path.write_text("time,event\n1,2\n2,1\n3,2\n4,1\n", encoding="utf-8")
    args.csv = str(r_path)
    out = CS._run(CS._km)(args)
    assert out.result["groups"]["all"]["events"] == 2
    assert out.result["groups"]["all"]["censored"] == 2
    assert "R `survival`" in out.result["event_coding"]


def test_e10_cli_km_rejects_a_nan_time(tmp_path):
    csv_path = tmp_path / "km.csv"
    csv_path.write_text("time,event\n1,1\n2,1\nnan,0\n3,1\n", encoding="utf-8")
    args = _Args(csv=str(csv_path), time_col="time", event_col="event", group_col=None, event_coding="0/1")
    with pytest.raises(UsageError) as err:
        CS._run(CS._km)(args)
    assert "finite" in str(err.value)


# ---- F17 (stats half): nan/inf reached the formulas and produced a confident wrong verdict

@pytest.mark.parametrize("bad", ["nan", "inf", "-inf", "Infinity", "NaN"])
def test_f17_finite_argparse_type_refuses_non_finite(bad):
    import argparse

    with pytest.raises(argparse.ArgumentTypeError):
        CS.finite(bad)


def test_f17_finite_accepts_ordinary_numbers():
    assert CS.finite("0.0004") == pytest.approx(0.0004)
    assert CS.finite("1e-5") == pytest.approx(1e-5)


def test_f17_maxaf_cannot_be_asked_about_a_nan_faf95():
    # the wrong verdict was "filtering AF is below the maximum credible AF" for faf95 = nan
    args = _Args(prevalence=0.0002, allelic=0.1, genetic=1.0, penetrance=1.0,
                 inheritance="monoallelic", an=None, faf95=CS.finite("0.001"))
    res = CS._run(CS._maxaf)(args).result
    assert res["too_common"] is True and "supports BS1" in res["reading"]


# ---- F18: out-of-range inputs were accepted or crashed with a traceback

@pytest.mark.parametrize("call,needle", [
    (lambda: S.recessive_from_prevalence(2), "prevalence must be in (0, 1)"),
    (lambda: S.recessive_from_prevalence(1), "prevalence must be in (0, 1)"),
    (lambda: S.recessive_from_prevalence(-0.1), "prevalence must be in (0, 1)"),
    (lambda: S.recessive_from_alleles([-0.5, 0.2]), "allele_freqs[0] must be in [0, 1]"),
    (lambda: S.recessive_from_alleles([0.6, 0.6]), "summed allele frequency must be in (0, 1)"),
    (lambda: S.recurrence("AD-inherited", penetrance=1.5), "penetrance must be in [0, 1]"),
    (lambda: S.recurrence("AD-de-novo", germline_mosaic_risk=2), "germline_mosaic_risk must be in [0, 1]"),
    (lambda: S.xlinked_carrier_posterior(2, 1), "prior must be in (0, 1)"),
    (lambda: S.xlinked_carrier_posterior(1, 0), "prior must be in (0, 1)"),
    (lambda: S.xlinked_carrier_posterior(0.5, -1), "unaffected_sons must be >= 0"),
    (lambda: S.max_credible_af(2, 0.1, 1.0, 1.0), "prevalence must be in (0, 1]"),
    (lambda: S.max_credible_af(0.1, 0.1, 1.0, 0), "penetrance must be in (0, 1]"),
    (lambda: S.max_tolerated_ac(1.5, 1000), "max_af must be in [0, 1]"),
    (lambda: S.denovo_enrichment(-1, 1000, 1e-5), "observed must be >= 0"),
    (lambda: S.nof1_pairs_needed(1.0, 1.0, alpha=0), "alpha must be in (0, 1)"),
    (lambda: S.nof1_pairs_needed(1.0, 1.0, power=1), "power must be in (0, 1)"),
    (lambda: S.segregation(ad_meioses=-1), "ad_meioses must be >= 0"),
    (lambda: S.wilson_ci(5, 2), "k must be <= 2"),
])
def test_f18_ranges_are_validated_with_the_bound_named(call, needle):
    with pytest.raises(ValueError) as err:
        call()
    assert needle in str(err.value)


def test_f18_segregation_does_not_overflow():
    r = S.segregation(ad_meioses=2000)
    assert r["likelihood_ratio"] is None
    assert r["lod"] == pytest.approx(2000 * math.log10(2))
    assert r["pp1_code"] == "PP1_Strong"


def test_f18_nan_is_refused_as_a_probability():
    with pytest.raises(ValueError) as err:
        S.recessive_from_prevalence(float("nan"))
    assert "finite" in str(err.value)


@pytest.mark.parametrize("args,needle", [
    (_Args(prevalence=1.0, allele_freqs=None), "prevalence must be in (0, 1)"),
    (_Args(prevalence=None, allele_freqs=[-0.5, 0.2]), "must be in [0, 1]"),
])
def test_f18_cli_stats_raises_usage_error_not_a_traceback(args, needle):
    with pytest.raises(UsageError) as err:
        CS._run(CS._carrier)(args)
    assert needle in str(err.value)


def test_f18_cli_recurrence_zero_division_is_a_usage_error():
    args = _Args(mode="XLR-bayes", prior=2.0, unaffected_sons=1, affected_sons=0,
                 penetrance=1.0, mosaic=0.01)
    with pytest.raises(UsageError) as err:
        CS._run(CS._recurrence)(args)
    assert "prior must be in (0, 1)" in str(err.value)


def test_f18_cli_segregation_overflow_is_not_a_traceback():
    args = _Args(ad_meioses=2000, ar_affected_sibs=0, ar_unaffected_sibs=0, xlr_male_meioses=0,
                 nonsegregations=0, unaffected_carriers=0, full_penetrance=False)
    assert CS._run(CS._segregation)(args).result["lod"] > 600


# ---- F19: the Poisson tail lost all precision and underflowed to 1.0

def test_f19_poisson_tail_keeps_its_precision():
    # reference values computed in exact arithmetic (mpmath) and quoted in the review
    assert S.poisson_sf(25, 2 * 1000 * 1e-5) == pytest.approx(2.122e-68, rel=1e-3)
    assert S.poisson_sf(2000, 2 * 1000000 * 0.0005) == pytest.approx(3.058e-170, rel=1e-3)
    # and the small-mu case the old code already got right is unchanged
    assert S.poisson_sf(5, 1.2) == pytest.approx(0.0077457883, rel=1e-8)


def test_f19_poisson_does_not_underflow_to_one():
    r = S.denovo_enrichment(2000, 1000000, 0.0005)
    assert 0 < r["p"] < 1e-100, r["p"]
    assert S.denovo_enrichment(25, 1000, 1e-5)["p"] < 1e-60


def test_f19_poisson_edges():
    assert S.poisson_sf(0, 5.0) == 1.0
    assert S.poisson_sf(1, 0.0) == 0.0
    # P(X >= k) for k well below the mean is close to, but never above, 1
    assert 0.99 < S.poisson_sf(1, 1000.0) <= 1.0
    # cross-checked against an exact high-precision summation of the same distribution:
    # P(X >= 1000) for Poisson(1000) = 0.5042052441802155
    assert S.poisson_sf(1000, 1000.0) == pytest.approx(0.50420524418, rel=1e-10)


def test_f19_binomial_tail_is_computed_in_log_space():
    assert S.binom_sf(3, 100, 0.01) == pytest.approx(0.0793732023, rel=1e-8)
    # a tail that exp(lbeta) alone would have turned into a hard zero
    assert 0 < S.binom_sf(200, 100000, 0.0005) < 1e-40


# ---- F28: segregation under-counted vs ClinGen, and there was no BS4

@pytest.mark.parametrize("kw,points,strength", [
    (dict(ad_meioses=1), 1.0, "Supporting"),
    (dict(ad_meioses=2), 2.0, "Moderate"),
    (dict(ad_meioses=3), 3.0, "Moderate"),
    (dict(ad_meioses=4), 4.0, "Strong"),
    (dict(ar_affected_sibs=1), 2.0, "Moderate"),
    (dict(ar_affected_sibs=2), 4.0, "Strong"),
    (dict(xlr_male_meioses=3), 3.0, "Moderate"),
])
def test_f28_segregation_points_match_clingen_table_3(kw, points, strength):
    r = S.segregation(**kw)
    assert r["points"] == pytest.approx(points)
    assert r["pp1_strength"] == strength


def test_f28_model_string_cites_the_exact_source():
    m = S.segregation(ad_meioses=2)["model"]
    assert "Biesecker et al. 2024" in m and "AJHG 111:24-38" in m and "Table 3" in m
    assert "38103548" in m


def test_f28_locus_points_are_capped_at_five():
    r = S.segregation(ar_affected_sibs=5)
    assert r["points_uncapped"] == pytest.approx(10.0) and r["points"] == pytest.approx(5.0)
    assert any("capped" in n for n in r["notes"])
    assert r["pp1_strength"] == "Strong", "PP1 is not applied above Strong"


def test_f28_unaffected_relatives_need_full_penetrance():
    with pytest.raises(ValueError) as err:
        S.segregation(ar_unaffected_sibs=3)
    assert "fully penetrant" in str(err.value)
    r = S.segregation(ar_unaffected_sibs=3, full_penetrance=True)
    assert r["points"] == pytest.approx(1.2)


def test_f28_bs4_is_available():
    r = S.segregation(nonsegregations=2)
    assert r["bs4_code"] == "BS4_Moderate" and r["bs4_points"] == pytest.approx(2.0)
    assert r["pp1_code"] is None
    assert S.segregation(nonsegregations=4)["bs4_code"] == "BS4_Strong"
    assert "not a published table" in r["bs4_basis"]


def test_f28_unaffected_carriers_are_bs4_under_full_penetrance():
    with pytest.raises(ValueError):
        S.segregation(unaffected_carriers=1)
    assert S.segregation(unaffected_carriers=1, full_penetrance=True)["bs4_code"] == "BS4"


def test_f28_both_directions_at_once_is_flagged():
    r = S.segregation(ad_meioses=3, nonsegregations=2)
    assert r["pp1_code"] and r["bs4_code"]
    assert any("cannot both be applied" in n for n in r["notes"])


def test_f28_jarvik_browning_lr_is_still_reported():
    r = S.segregation(ad_meioses=3)
    assert r["likelihood_ratio"] == pytest.approx(8.0)
    assert "Jarvik & Browning 2016" in r["likelihood_ratio_model"]


def test_f28_cli_exposes_the_new_counts():
    args = _Args(ad_meioses=0, ar_affected_sibs=0, ar_unaffected_sibs=0, xlr_male_meioses=2,
                 nonsegregations=3, unaffected_carriers=0, full_penetrance=False)
    res = CS._run(CS._segregation)(args).result
    assert res["pp1_code"] == "PP1_Moderate" and res["bs4_code"] == "BS4_Moderate"
