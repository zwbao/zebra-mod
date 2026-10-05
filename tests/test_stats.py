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
    assert S.segregation(ad_meioses=1)["pp1_code"] is None or S.segregation(ad_meioses=1)["pp1_code"] == "PP1"
    assert S.segregation(ad_meioses=2)["pp1_code"] == "PP1"
    assert S.segregation(ad_meioses=3)["pp1_code"] == "PP1_Moderate"
    assert S.segregation(ad_meioses=5)["pp1_code"] == "PP1_Strong"
    assert S.segregation(ad_meioses=10)["pp1_code"] == "PP1_Strong"


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
