"""zebra stats: classical statistics for rare-disease genetics (offline)."""

from __future__ import annotations

import argparse
import csv
import math
from typing import Any, Dict, List

from zebra import stats as S
from zebra.core import Outcome, UsageError


def finite(raw: str) -> float:
    """argparse float that refuses nan/inf.

    argparse's own `type=float` accepts "nan" and "inf", which then travel through the
    formulas into the JSON envelope as the non-standard literals NaN and Infinity (and,
    for `--faf95 nan`, produce a confidently wrong "not too common" verdict, because
    every comparison with NaN is False).
    """
    try:
        v = float(raw)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError(f"{raw!r} is not a number") from None
    if not math.isfinite(v):
        raise argparse.ArgumentTypeError(f"{raw!r} is not a finite number")
    return v


def _fmt(d: Dict[str, Any]) -> str:
    lines = []
    for k, v in d.items():
        if isinstance(v, float):
            v = f"{v:.4g}"
        elif isinstance(v, tuple):
            v = "(" + ", ".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in v) + ")"
        lines.append(f"{k}: {v}")
    return "\n".join(lines)


def _run(fn):
    def handler(args: argparse.Namespace) -> Outcome:
        try:
            res = fn(args)
        # Every out-of-range argument should read as a usage error, never as a traceback:
        # ZeroDivisionError and OverflowError are what the un-validated formulas used to
        # raise, and are kept here as a backstop for any path that still slips through.
        except (ValueError, ZeroDivisionError, OverflowError, ArithmeticError) as err:
            raise UsageError(str(err) or type(err).__name__) from None
        return Outcome(res, text=_fmt(res) if isinstance(res, dict) else None, query={k: v for k, v in vars(args).items() if k not in ("func",) and not callable(v)})

    return handler


def _segregation(a):
    return S.segregation(a.ad_meioses, a.ar_affected_sibs, a.ar_unaffected_sibs,
                         xlr_male_meioses=a.xlr_male_meioses, nonsegregations=a.nonsegregations,
                         unaffected_carriers=a.unaffected_carriers, full_penetrance=a.full_penetrance)


def _maxaf(a):
    res = S.max_credible_af(a.prevalence, a.allelic, a.genetic, a.penetrance, a.inheritance)
    if a.an:
        res["max_tolerated_ac"] = S.max_tolerated_ac(res["max_credible_af"], a.an)
        res["an"] = a.an
    if a.faf95 is not None:
        res["faf95"] = a.faf95
        res["too_common"] = a.faf95 > res["max_credible_af"]
        res["reading"] = ("filtering AF exceeds the maximum credible AF: too common to be a fully penetrant cause (supports BS1)"
                          if res["too_common"] else "filtering AF is below the maximum credible AF: frequency does not argue against causality")
    return res


def _carrier(a):
    if a.allele_freqs:
        return S.recessive_from_alleles(a.allele_freqs)
    if a.prevalence:
        return S.recessive_from_prevalence(a.prevalence)
    raise ValueError("give --prevalence (affected fraction) or --allele-freqs (P/LP allele frequencies)")


def _recurrence(a):
    if a.mode == "XLR-bayes":
        return S.xlinked_carrier_posterior(a.prior, a.unaffected_sons, a.affected_sons)
    return S.recurrence(a.mode, penetrance=a.penetrance, germline_mosaic_risk=a.mosaic)


def _fisher(a):
    return S.fisher_exact(a.a, a.b, a.c, a.d)


def _burden(a):
    return S.burden(a.case_carriers, a.case_n, a.control_carriers, a.control_n)


def _denovo(a):
    return S.denovo_enrichment(a.observed, a.trios, a.mu)


def _km(a):
    groups: Dict[str, List[tuple]] = {}
    try:
        fh = open(a.csv, newline="", encoding="utf-8")
    except OSError as err:
        raise ValueError(f"{a.csv}: cannot be read ({err})") from None
    with fh:
        for lineno, row in enumerate(csv.DictReader(fh), start=2):
            for col in (a.time_col, a.event_col):
                if col not in row or row[col] in (None, ""):
                    raise ValueError(f"{a.csv} line {lineno}: column {col!r} is missing or empty "
                                     f"(columns found: {', '.join(k for k in row if k)})")
            try:
                t = float(row[a.time_col])
                e = int(float(row[a.event_col]))
            except (TypeError, ValueError) as err:
                raise ValueError(f"{a.csv} line {lineno}: need numeric {a.time_col!r} and {a.event_col!r} ({err})") from None
            g = row.get(a.group_col, "all") if a.group_col else "all"
            groups.setdefault(g, []).append((t, e))
    if not groups:
        raise ValueError(f"{a.csv}: no data rows")
    coding = a.event_coding
    out: Dict[str, Any] = {
        "model": "Kaplan-Meier with Greenwood 95% CI (linear, clipped to [0,1]); log-rank for two groups",
        "event_coding": f"{coding}: {S.EVENT_CODINGS[coding][2]}",
        "groups": {},
    }
    for g, rows in groups.items():
        # the validated rows are what gets counted, so `events` cannot disagree with the curve
        checked = S.check_survival_input([r[0] for r in rows], [r[1] for r in rows], coding)
        curve = S.kaplan_meier([r[0] for r in checked], [r[1] for r in checked])
        out["groups"][g] = {"n": len(checked), "events": sum(r[1] for r in checked),
                            "censored": sum(1 for r in checked if r[1] == 0),
                            "median": S.median_survival(curve), "curve": curve}
    if len(groups) == 2:
        (ga, ra), (gb, rb) = list(groups.items())
        out["logrank"] = {"groups": [ga, gb],
                          **S.logrank([r[0] for r in ra], [r[1] for r in ra],
                                      [r[0] for r in rb], [r[1] for r in rb], event_coding=coding)}
    return out


def _nof1(a):
    if a.treatment and a.control:
        return S.nof1_analyze(a.treatment, a.control)
    if a.effect and a.sd_diff:
        return S.nof1_pairs_needed(a.effect, a.sd_diff, a.alpha, a.power)
    raise ValueError("design: --effect and --sd-diff; analysis: --treatment and --control (matched period values)")


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("stats", help="classical statistics: segregation, maxaf, carrier, recurrence, fisher, burden, denovo, km, nof1")
    ss = p.add_subparsers(dest="action", metavar="<method>")

    q = ss.add_parser("segregation", help="co-segregation → PP1 and non-segregation → BS4, as ClinGen points (Biesecker et al. 2024)")
    q.add_argument("--ad-meioses", type=int, default=0, help="informative meioses in a dominant pedigree (1.0 point each)")
    q.add_argument("--xlr-male-meioses", type=int, default=0, help="informative male meioses, X-linked recessive (1.0 point each)")
    q.add_argument("--ar-affected-sibs", type=int, default=0, help="additional affected relatives with the same biallelic genotype (2.0 points each)")
    q.add_argument("--ar-unaffected-sibs", type=int, default=0, help="unaffected relatives without the biallelic genotype (0.4 points each; needs --full-penetrance)")
    q.add_argument("--nonsegregations", type=int, default=0, help="BS4: affected relatives who do NOT carry the variant")
    q.add_argument("--unaffected-carriers", type=int, default=0, help="BS4: unaffected carriers in a fully penetrant dominant family (needs --full-penetrance)")
    q.add_argument("--full-penetrance", action="store_true",
                   help="assert full penetrance; ClinGen counts unaffected individuals only under it")
    q.set_defaults(func=_run(_segregation))

    q = ss.add_parser("maxaf", help="maximum credible allele frequency (Whiffin 2017) → BS1 check")
    q.add_argument("--prevalence", type=finite, required=True, help="affected fraction, e.g. 0.0000625 for 1 in 16,000")
    q.add_argument("--allelic", type=finite, required=True, help="max share of cases from one allele")
    q.add_argument("--genetic", type=finite, default=1.0, help="share of cases due to this gene")
    q.add_argument("--penetrance", type=finite, default=1.0)
    q.add_argument("--inheritance", choices=("monoallelic", "biallelic"), default="monoallelic")
    q.add_argument("--an", type=int, help="allele number of the reference sample (for max tolerated AC)")
    q.add_argument("--faf95", type=finite, help="the variant's filtering AF (gnomAD faf95) to compare")
    q.set_defaults(func=_run(_maxaf))

    q = ss.add_parser("carrier", help="carrier frequency / genetic prevalence (Hardy-Weinberg)")
    q.add_argument("--prevalence", type=finite)
    q.add_argument("--allele-freqs", type=finite, nargs="*")
    q.set_defaults(func=_run(_carrier))

    q = ss.add_parser("recurrence", help="recurrence risk by inheritance; XLR-bayes for carrier posterior")
    q.add_argument("--mode", required=True, choices=("AR", "AD-inherited", "AD-de-novo", "XLR-carrier-mother", "XLR-bayes"))
    q.add_argument("--penetrance", type=finite, default=1.0)
    q.add_argument("--mosaic", type=finite, default=0.01, help="germline mosaicism risk for apparently de novo")
    q.add_argument("--prior", type=finite, default=0.5, help="XLR-bayes: prior carrier probability")
    q.add_argument("--unaffected-sons", type=int, default=0)
    q.add_argument("--affected-sons", type=int, default=0)
    q.set_defaults(func=_run(_recurrence))

    q = ss.add_parser("fisher", help="2x2 Fisher exact test [[a, b], [c, d]]")
    for k in ("a", "b", "c", "d"):
        q.add_argument(f"--{k}", type=int, required=True)
    q.set_defaults(func=_run(_fisher))

    q = ss.add_parser("burden", help="gene burden: carriers in cases vs controls")
    q.add_argument("--case-carriers", type=int, required=True)
    q.add_argument("--case-n", type=int, required=True)
    q.add_argument("--control-carriers", type=int, required=True)
    q.add_argument("--control-n", type=int, required=True)
    q.set_defaults(func=_run(_burden))

    q = ss.add_parser("denovo", help="de novo enrichment in one gene (Poisson)")
    q.add_argument("--observed", type=int, required=True)
    q.add_argument("--trios", type=int, required=True)
    q.add_argument("--mu", type=finite, required=True, help="per-haploid mutation rate for the class in this gene")
    q.set_defaults(func=_run(_denovo))

    q = ss.add_parser("km", help="Kaplan-Meier (natural history) from a CSV; log-rank for two groups")
    q.add_argument("--csv", required=True)
    q.add_argument("--time-col", default="time")
    q.add_argument("--event-col", default="event")
    q.add_argument("--group-col")
    q.add_argument("--event-coding", choices=tuple(sorted(S.EVENT_CODINGS)), default="0/1",
                   help="status column convention: 0/1 (0 censored, 1 event; default) or 1/2 "
                        "(the R survival package's coding, 1 censored, 2 event). Any other code is refused.")
    q.set_defaults(func=_run(_km))

    q = ss.add_parser("nof1", help="N-of-1 trial: pairs needed (design) or paired analysis")
    q.add_argument("--effect", type=finite)
    q.add_argument("--sd-diff", type=finite)
    q.add_argument("--alpha", type=finite, default=0.05)
    q.add_argument("--power", type=finite, default=0.8)
    q.add_argument("--treatment", type=finite, nargs="*")
    q.add_argument("--control", type=finite, nargs="*")
    q.set_defaults(func=_run(_nof1))
