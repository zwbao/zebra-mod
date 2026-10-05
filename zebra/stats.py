"""Classical statistics for rare-disease genetics, standard library only.

Each function states its model in its docstring; the CLI (`zebra stats`)
prints the model with the number so nobody reads a figure without its
assumptions. Nothing here needs network access.
"""

from __future__ import annotations

import math
from statistics import NormalDist
from typing import Any, Dict, List, Optional, Sequence, Tuple

_N = NormalDist()

# ClinGen/Tavtigian 2018 odds of pathogenicity per evidence strength (prior 0.10, exponent 2).
ODDS_PATH = {"Supporting": 2.08, "Moderate": 4.33, "Strong": 18.7, "VeryStrong": 350.0}
# Tavtigian et al. 2020 (Hum Mutat 41:1734) point scale, used to read Bayesian points
# back as an ACMG strength.
POINTS_TO_STRENGTH = ((8, "VeryStrong"), (4, "Strong"), (2, "Moderate"), (1, "Supporting"))


# ---------------------------------------------------------------- argument checks
# Every probability and count argument is range-checked here, and the error names the
# bound it violated, so the CLI can turn it into a usage error instead of a traceback.

def _prob(name: str, value: Any, lo: float = 0.0, hi: float = 1.0,
          lo_open: bool = False, hi_open: bool = False) -> float:
    try:
        v = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a number in {'(' if lo_open else '['}{lo}, {hi}{')' if hi_open else ']'}") from None
    if not math.isfinite(v):
        raise ValueError(f"{name} must be a finite number, not {value!r}")
    if v < lo or (lo_open and v == lo) or v > hi or (hi_open and v == hi):
        raise ValueError(
            f"{name} must be in {'(' if lo_open else '['}{lo}, {hi}{')' if hi_open else ']'}, got {v!r}"
        )
    return v


def _count(name: str, value: Any, maximum: Optional[int] = None) -> int:
    try:
        v = int(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a whole number >= 0") from None
    if v < 0:
        raise ValueError(f"{name} must be >= 0, got {v}")
    if maximum is not None and v > maximum:
        raise ValueError(f"{name} must be <= {maximum}, got {v}")
    return v


# ---------------------------------------------------------------- distributions

def _betacf(a: float, b: float, x: float) -> float:
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1.0 + aa / c
        c = c if abs(c) > 1e-300 else 1e-300
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1.0 + aa / c
        c = c if abs(c) > 1e-300 else 1e-300
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 3e-14:
            break
    return h


def betainc(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta I_x(a, b)."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x)
    if x < (a + 1) / (a + b + 2):
        # exp(lbeta) alone underflows to 0.0 once lbeta < -745, which silently turns a
        # representable tail probability (say 1e-170) into a hard zero. Multiply inside
        # the exponent instead; the continued fraction is positive by construction.
        cf = _betacf(a, b, x) / a
        return math.exp(lbeta + math.log(cf)) if cf > 0 else 0.0
    # Here the result is O(1), so the subtraction loses nothing worth keeping.
    return 1.0 - math.exp(lbeta) * _betacf(b, a, 1.0 - x) / b


def t_sf_two_sided(t: float, df: float) -> float:
    return betainc(df / 2.0, 0.5, df / (df + t * t))


def chi2_sf_1df(x: float) -> float:
    return math.erfc(math.sqrt(max(x, 0.0) / 2.0))


def _poisson_tail(k: int, mu: float, upper: bool) -> float:
    """Poisson(mu) mass over i >= k (upper) or 0 <= i <= k-1 (lower).

    Summed outwards from the largest term in the range and scaled by it, so nothing
    underflows: `exp(-mu)` on its own is 0.0 in double precision once mu > 745, which
    is why the naive `1 - cdf` form returns p = 1.0 for every large expectation.
    """
    lo = k if upper else 0
    hi: Optional[int] = None if upper else k - 1
    star = max(lo, min(int(mu), hi)) if hi is not None else max(lo, int(mu))
    log_star = -mu + star * math.log(mu) - math.lgamma(star + 1.0)
    total = 1.0
    term, i = 1.0, star
    while i > lo:  # downwards: P(i-1)/P(i) = i / mu
        term *= i / mu
        i -= 1
        if term <= 0.0:
            break
        total += term
        if term < 1e-18 * total:
            break
    term, i = 1.0, star
    while hi is None or i < hi:  # upwards: P(i+1)/P(i) = mu / (i+1)
        i += 1
        term *= mu / i
        if term <= 0.0:
            break
        total += term
        if term < 1e-18 * total:
            break
    return math.exp(log_star + math.log(total))


def poisson_sf(k: int, mu: float) -> float:
    """P(X >= k) for X ~ Poisson(mu), computed in log space from the shorter tail."""
    if k <= 0:
        return 1.0
    if mu <= 0:
        return 0.0
    if k > mu:
        # The answer is small: sum the upper tail directly, so no cancellation.
        return min(1.0, _poisson_tail(k, mu, upper=True))
    # The answer is O(1) and the lower tail is the shorter sum; the subtraction is safe
    # because the result is not small.
    return max(0.0, min(1.0, 1.0 - _poisson_tail(k, mu, upper=False)))


def binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p)."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    if not 0.0 <= p <= 1.0:
        raise ValueError("p must be in [0, 1]")
    return betainc(k, n - k + 1, p)


def wilson_ci(k: int, n: int, conf: float = 0.95) -> Tuple[float, float]:
    n = _count("n", n)
    k = _count("k", k, maximum=n)
    conf = _prob("conf", conf, 0.0, 1.0, lo_open=True, hi_open=True)
    if n == 0:
        return (0.0, 1.0)
    z = _N.inv_cdf(1 - (1 - conf) / 2)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


# ---------------------------------------------------------------- Fisher / burden

def _log_hyper(a: int, b: int, c: int, d: int) -> float:
    def lc(n: int, k: int) -> float:
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    return lc(a + b, a) + lc(c + d, c) - lc(a + b + c + d, a + c)


def fisher_exact(a: int, b: int, c: int, d: int) -> Dict[str, float]:
    """Two-sided Fisher exact test on [[a, b], [c, d]] (sum of tables no more likely than observed).

    Odds ratio with Haldane-Anscombe correction (+0.5) when a cell is zero, and
    its Woolf 95% confidence interval.
    """
    if min(a, b, c, d) < 0:
        raise ValueError("counts must be non-negative")
    r1, c1, n = a + b, a + c, a + b + c + d
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    p_obs = _log_hyper(a, b, c, d)
    p = 0.0
    p_less = 0.0
    p_greater = 0.0
    for x in range(lo, hi + 1):
        lp = _log_hyper(x, r1 - x, c1 - x, n - r1 - c1 + x)
        prob = math.exp(lp)
        if lp <= p_obs + 1e-7:
            p += prob
        if x <= a:
            p_less += prob
        if x >= a:
            p_greater += prob
    aa, bb, cc, dd = (a, b, c, d) if min(a, b, c, d) > 0 else (a + 0.5, b + 0.5, c + 0.5, d + 0.5)
    odds = (aa * dd) / (bb * cc)
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return {
        "p_two_sided": min(1.0, p),
        "p_greater": min(1.0, p_greater),
        "p_less": min(1.0, p_less),
        "odds_ratio": odds,
        "or_ci95": (math.exp(math.log(odds) - 1.96 * se), math.exp(math.log(odds) + 1.96 * se)),
        "haldane_corrected": min(a, b, c, d) == 0,
    }


def burden(case_carriers: int, case_n: int, control_carriers: int, control_n: int) -> Dict[str, Any]:
    """Gene-level carrier burden (collapsing test): carriers vs non-carriers in cases and controls.

    One-sided interest is enrichment in cases (p_greater). Controls from a
    public resource (gnomAD) differ in sequencing and calling; say so.
    """
    if case_carriers > case_n or control_carriers > control_n:
        raise ValueError("carriers cannot exceed sample size")
    f = fisher_exact(case_carriers, case_n - case_carriers, control_carriers, control_n - control_carriers)
    return {
        "model": "collapsing burden: Fisher exact on carriers vs non-carriers (one-sided p_greater = enrichment in cases)",
        "case_rate": case_carriers / case_n if case_n else None,
        "control_rate": control_carriers / control_n if control_n else None,
        **f,
    }


# ---------------------------------------------------------------- de novo enrichment

def denovo_enrichment(observed: int, trios: int, mu: float) -> Dict[str, Any]:
    """Poisson test for excess de novo variants in one gene (Samocha et al. 2014 framework).

    Expected = 2 * trios * mu, where mu is the per-haploid-genome, per-generation
    mutation rate for the variant class in this gene (e.g. from the Samocha
    model or gnomAD mutation-rate tables). One-sided p = P(X >= observed).
    """
    observed = _count("observed", observed)
    if trios <= 0 or mu <= 0:
        raise ValueError("trios and mu must be positive")
    mu = _prob("mu", mu, 0.0, 1.0, lo_open=True)
    if not math.isfinite(float(trios)):
        raise ValueError("trios must be a finite whole number")
    expected = 2.0 * trios * mu
    return {
        "model": "Poisson(2 * trios * mu); one-sided P(X >= observed)",
        "observed": observed,
        "expected": expected,
        "ratio": observed / expected if expected else None,
        "p": poisson_sf(observed, expected),
        "exome_wide_alpha": "0.05 / (19000 genes * classes tested) ~ 2.6e-6 for one class",
    }


# ---------------------------------------------------------------- frequencies

def max_credible_af(prevalence: float, allelic: float, genetic: float, penetrance: float,
                    inheritance: str = "monoallelic") -> Dict[str, Any]:
    """Maximum credible population allele frequency (Whiffin et al. 2017, Genet Med 19:1151).

    monoallelic: q_max = prevalence * genetic * allelic / (2 * penetrance)
    biallelic:   q_max = allelic * sqrt(prevalence * genetic / penetrance)
    prevalence: affected fraction of the population; genetic: share of cases due
    to this gene; allelic: largest share of those cases one allele can explain;
    penetrance: chance a carrier (monoallelic) or a biallelic genotype is affected.
    A variant whose filtering AF (gnomAD faf95) exceeds q_max is too common to be
    a fully penetrant cause (supports BS1).
    """
    prevalence = _prob("prevalence", prevalence, 0.0, 1.0, lo_open=True)
    allelic = _prob("allelic", allelic, 0.0, 1.0, lo_open=True)
    genetic = _prob("genetic", genetic, 0.0, 1.0, lo_open=True)
    penetrance = _prob("penetrance", penetrance, 0.0, 1.0, lo_open=True)
    if inheritance == "monoallelic":
        q = prevalence * genetic * allelic / (2 * penetrance)
        formula = "prevalence * genetic * allelic / (2 * penetrance)"
    elif inheritance == "biallelic":
        q = allelic * math.sqrt(prevalence * genetic / penetrance)
        formula = "allelic * sqrt(prevalence * genetic / penetrance)"
    else:
        raise ValueError("inheritance must be monoallelic or biallelic")
    return {"model": f"Whiffin 2017 maximum credible AF, {inheritance}: {formula}", "max_credible_af": q,
            "inputs": {"prevalence": prevalence, "allelic": allelic, "genetic": genetic, "penetrance": penetrance}}


def max_tolerated_ac(max_af: float, an: int, conf: float = 0.95) -> int:
    """Largest allele count still consistent (one-sided, `conf`) with AF <= max_af among `an` alleles."""
    max_af = _prob("max_af", max_af)
    an = _count("an", an)
    conf = _prob("conf", conf, 0.0, 1.0, lo_open=True, hi_open=True)
    k = 0
    while binom_sf(k + 1, an, max_af) > 1 - conf:
        k += 1
        if k > an:
            break
    return k


def recessive_from_prevalence(prevalence: float) -> Dict[str, float]:
    """Hardy-Weinberg for a fully penetrant recessive disorder: q = sqrt(prevalence), carriers = 2pq."""
    # prevalence must be strictly inside (0, 1): at 1 every allele is pathogenic and 2pq is 0.
    prevalence = _prob("prevalence", prevalence, 0.0, 1.0, lo_open=True, hi_open=True)
    q = math.sqrt(prevalence)
    p = 1 - q
    return {"model": "Hardy-Weinberg, full penetrance, panmixia", "pathogenic_allele_freq": q,
            "carrier_freq": 2 * p * q, "one_in_carriers": 1 / (2 * p * q) if q > 0 else None}


def recessive_from_alleles(allele_freqs: Sequence[float]) -> Dict[str, Any]:
    """Genetic prevalence from summed pathogenic allele frequencies: (sum q)^2 (homozygotes + compound hets)."""
    if not allele_freqs:
        raise ValueError("give at least one allele frequency in [0, 1]")
    q = sum(_prob(f"allele_freqs[{i}]", v) for i, v in enumerate(allele_freqs))
    if not 0 < q < 1:
        raise ValueError(f"summed allele frequency must be in (0, 1), got {q!r}")
    prev = q * q
    return {"model": "genetic prevalence = (sum of P/LP allele frequencies)^2; carriers = 2q(1-q)",
            "summed_q": q, "genetic_prevalence": prev, "one_in": 1 / prev if prev > 0 else None,
            "carrier_freq": 2 * q * (1 - q)}


# ---------------------------------------------------------------- segregation

def segregation(ad_meioses: int = 0, ar_affected_sibs: int = 0, ar_unaffected_sibs: int = 0) -> Dict[str, Any]:
    """Cosegregation likelihood ratio by counting meioses (Jarvik & Browning 2016, AJHG 98:1077), full penetrance.

    Dominant: each informative meiosis that transmits variant with disease
    (affected relative beyond the proband who carries it, or an obligate
    carrier link) halves the chance of chance cosegregation: LR = 2^m.
    Recessive: each additional affected sibling with the same two alleles adds
    LR 4; each unaffected sibling without that genotype adds 4/3.
    Mapped to ACMG PP1 strength with Tavtigian 2018 odds (Supporting >= 2.08,
    Moderate >= 4.33, Strong >= 18.7). Reduced penetrance or phenocopies lower
    the LR: run a full likelihood model before leaning on Strong.
    """
    lr = (2.0 ** ad_meioses) * (4.0 ** ar_affected_sibs) * ((4.0 / 3.0) ** ar_unaffected_sibs)
    strength = None
    for name in ("VeryStrong", "Strong", "Moderate", "Supporting"):
        if lr >= ODDS_PATH[name]:
            strength = name
            break
    if strength == "VeryStrong":
        strength = "Strong"  # PP1 is not applied above Strong
    return {
        "model": "counting meioses, full penetrance, no phenocopies (Jarvik & Browning 2016)",
        "likelihood_ratio": lr,
        "lod": math.log10(lr) if lr > 0 else None,
        "pp1_strength": strength,
        "pp1_code": None if strength is None else ("PP1" if strength == "Supporting" else f"PP1_{strength}"),
        "thresholds": {k: v for k, v in ODDS_PATH.items() if k != "VeryStrong"},
    }


# ---------------------------------------------------------------- recurrence risk

def xlinked_carrier_posterior(prior: float, unaffected_sons: int, affected_sons: int = 0) -> Dict[str, Any]:
    """Bayesian carrier probability for a woman at risk of an X-linked recessive allele.

    Each unaffected son multiplies the carrier likelihood by 1/2 (vs 1 if not a
    carrier); an affected son makes her an obligate carrier (germline mosaicism aside).
    """
    if affected_sons > 0:
        return {"model": "Bayes, X-linked recessive", "posterior_carrier": 1.0, "note": "affected son: obligate carrier unless germline mosaic or de novo"}
    like_carrier = 0.5 ** unaffected_sons
    post = prior * like_carrier / (prior * like_carrier + (1 - prior))
    return {"model": "Bayes, X-linked recessive", "prior": prior, "unaffected_sons": unaffected_sons,
            "posterior_carrier": post, "risk_next_son_affected": post * 0.5}


def recurrence(mode: str, **kw: Any) -> Dict[str, Any]:
    """Recurrence risk for the next pregnancy of the same parents, by inheritance pattern."""
    if mode == "AR":
        return {"mode": mode, "risk": 0.25, "note": "both parents confirmed heterozygous carriers"}
    if mode == "AD-inherited":
        pen = float(kw.get("penetrance", 1.0))
        return {"mode": mode, "risk": 0.5 * pen, "note": "affected/carrier parent; risk = 0.5 x penetrance"}
    if mode == "AD-de-novo":
        mosaic = float(kw.get("germline_mosaic_risk", 0.01))
        return {"mode": mode, "risk": mosaic, "note": "variant absent in both parents' blood; residual risk is germline mosaicism (commonly quoted ~1%, gene-dependent)"}
    if mode == "XLR-carrier-mother":
        return {"mode": mode, "risk_son_affected": 0.5, "risk_daughter_carrier": 0.5, "risk_any_child_affected": 0.25}
    raise ValueError("mode must be AR, AD-inherited, AD-de-novo or XLR-carrier-mother")


# ---------------------------------------------------------------- natural history

def kaplan_meier(times: Sequence[float], events: Sequence[int]) -> List[Dict[str, float]]:
    """Kaplan-Meier survival estimate with Greenwood variance; event=1, censored=0."""
    data = sorted(zip(times, events), key=lambda x: (x[0], -x[1]))
    at_risk = len(data)
    s = 1.0
    var_sum = 0.0
    out = []
    i = 0
    while i < len(data):
        t = data[i][0]
        d = sum(1 for tt, e in data if tt == t and e == 1)
        c = sum(1 for tt, e in data if tt == t and e == 0)
        if d > 0:
            s *= 1 - d / at_risk
            if at_risk - d > 0:
                var_sum += d / (at_risk * (at_risk - d))
            se = s * math.sqrt(var_sum)
            out.append({"time": t, "at_risk": at_risk, "events": d, "survival": s,
                        "ci95_low": max(0.0, s - 1.96 * se), "ci95_high": min(1.0, s + 1.96 * se)})
        at_risk -= d + c
        i += d + c
    return out


def median_survival(curve: List[Dict[str, float]]) -> Optional[float]:
    for row in curve:
        if row["survival"] <= 0.5:
            return row["time"]
    return None


def logrank(times_a: Sequence[float], events_a: Sequence[int], times_b: Sequence[float], events_b: Sequence[int]) -> Dict[str, float]:
    """Two-group log-rank test (chi-square, 1 df)."""
    pooled = sorted(set(t for t, e in zip(times_a, events_a) if e) | set(t for t, e in zip(times_b, events_b) if e))
    o_minus_e = 0.0
    var = 0.0
    for t in pooled:
        na = sum(1 for x in times_a if x >= t)
        nb = sum(1 for x in times_b if x >= t)
        da = sum(1 for x, e in zip(times_a, events_a) if x == t and e)
        db = sum(1 for x, e in zip(times_b, events_b) if x == t and e)
        n, d = na + nb, da + db
        if n < 2:
            continue
        o_minus_e += da - d * na / n
        var += d * (na / n) * (nb / n) * (n - d) / (n - 1)
    chi2 = o_minus_e ** 2 / var if var > 0 else 0.0
    return {"chi2": chi2, "p": chi2_sf_1df(chi2), "observed_minus_expected_a": o_minus_e}


# ---------------------------------------------------------------- N-of-1 trials

def nof1_pairs_needed(effect: float, sd_diff: float, alpha: float = 0.05, power: float = 0.8) -> Dict[str, Any]:
    """Treatment/control period pairs needed in a single-patient crossover (normal approximation, paired).

    effect: smallest worthwhile mean difference; sd_diff: SD of within-pair
    differences (period-to-period noise). Add washout periods; carryover and
    a progressive disease break the exchangeability this assumes.
    """
    if effect <= 0 or sd_diff <= 0:
        raise ValueError("effect and sd_diff must be positive")
    z = _N.inv_cdf(1 - alpha / 2) + _N.inv_cdf(power)
    n = math.ceil((z * sd_diff / effect) ** 2)
    return {"model": "paired normal approximation", "pairs": max(n, 2), "standardized_effect": effect / sd_diff}


def nof1_analyze(treatment: Sequence[float], control: Sequence[float]) -> Dict[str, Any]:
    """Paired t-test on matched treatment/control periods of one patient."""
    if len(treatment) != len(control) or len(treatment) < 2:
        raise ValueError("need >= 2 matched treatment/control pairs")
    diffs = [a - b for a, b in zip(treatment, control)]
    n = len(diffs)
    mean = sum(diffs) / n
    sd = math.sqrt(sum((d - mean) ** 2 for d in diffs) / (n - 1))
    if sd == 0:
        return {"model": "paired t-test", "pairs": n, "mean_difference": mean, "sd": 0.0, "t": None, "p_two_sided": None,
                "note": "no variation between pairs; the test is undefined"}
    t = mean / (sd / math.sqrt(n))
    half = _t_quantile(0.975, n - 1) * sd / math.sqrt(n)
    return {"model": "paired t-test", "pairs": n, "mean_difference": mean, "sd": sd, "t": t, "df": n - 1,
            "p_two_sided": t_sf_two_sided(t, n - 1), "ci95": (mean - half, mean + half)}


def _t_quantile(p: float, df: float) -> float:
    lo, hi = 0.0, 1000.0
    target = 2 * (1 - p)
    for _ in range(200):
        mid = (lo + hi) / 2
        if t_sf_two_sided(mid, df) > target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2
