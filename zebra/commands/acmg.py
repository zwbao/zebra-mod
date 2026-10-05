"""zebra acmg: ACMG/AMP points and combining rules; code suggestions from data."""

from __future__ import annotations

import argparse
from typing import Any, Dict, Optional

from zebra import acmg, stats
from zebra.core import Outcome, UsageError


def _classify(args: argparse.Namespace) -> Outcome:
    if not args.codes:
        raise UsageError("give evidence codes, e.g. PVS1 PM2_Supporting PP3")
    try:
        res = acmg.classify(args.codes)
    except ValueError as err:
        raise UsageError(str(err)) from None
    lines = [" + ".join(c["label"] for c in res["codes"]),
             f"points: {res['points'] if res['points'] is not None else 'BA1 (stand-alone)'} → {res['classification']}",
             f"2015 combining rules → {res['classification_richards_2015']}"]
    if res.get("capped_at_uncertain"):
        lines.append(f"capped: the points total alone would read {res['capped_at_uncertain']}; computational evidence "
                     f"cannot classify on its own")
    if res.get("note"):
        lines.append(f"note: {res['note']}")
    for w in res["warnings"]:
        lines.append(f"warning: {w}")
    lines.append("research-grade: confirm with an accredited laboratory before clinical use")
    return Outcome(res, text="\n".join(lines), query={"codes": args.codes})


def _max_credible_af(args: argparse.Namespace) -> Dict[str, Any]:
    """The Whiffin 2017 maximum credible AF from the flags, or why it is absent."""
    if args.prevalence is None:
        return {
            "max_credible_af": None,
            "why_absent": ("not given: pass --prevalence (and --allelic, --genetic, --penetrance, "
                           "--inheritance-mode) so the Whiffin et al. 2017 maximum credible AF can be computed. "
                           "Without it BS1 cannot be assessed."),
        }
    if args.allelic is None:
        raise UsageError("--prevalence needs --allelic (the largest share of cases one allele can explain)")
    try:
        res = stats.max_credible_af(args.prevalence, args.allelic, args.genetic, args.penetrance,
                                    args.inheritance_mode)
    except ValueError as err:
        raise UsageError(str(err)) from None
    return {"max_credible_af": res["max_credible_af"], "model": res["model"], "inputs": res["inputs"]}


def _suggest(args: argparse.Namespace) -> Outcome:
    from zebra.sources import variant as variant_src

    card = variant_src.card(args.variant, assembly=args.assembly)
    data: Dict[str, Any] = dict(card.result.get("acmg_inputs", {}))
    # The caller may know ClinVar's aggregate classification already; the variant card
    # carries it. BA1 must never be offered against a P/LP assertion without saying so.
    clinvar = card.result.get("clinvar") or {}
    if isinstance(clinvar, dict) and not data.get("clinvar_classification"):
        assertion = clinvar.get("germline_classification") or clinvar.get("classification") or clinvar.get("clinical_significance")
        if assertion:
            data["clinvar_classification"] = assertion

    maxaf = _max_credible_af(args)
    computed: Optional[float] = maxaf["max_credible_af"]
    res_s = acmg.suggest(data, inheritance=None if args.inheritance == "unknown" else args.inheritance,
                         max_credible_af=computed)
    suggestions = res_s["suggested"]
    res = {
        "variant": card.result.get("variant"),
        "inputs": data,
        "max_credible_af": maxaf,
        "suggested": suggestions,
        "not_assessed": res_s["not_assessed"],
        "caveats": res_s["caveats"],
        "not_from_data": [
            "PVS1 (null variant in a gene where loss of function causes disease: ClinGen PVS1 decision tree, Abou Tayoun 2018)",
            "PS1/PM5 (same amino-acid change / another missense at the residue known pathogenic)",
            "PS2/PM6 (de novo, confirmed or assumed)", "PS3/BS3 (functional studies, calibrated per Brnich 2019)",
            "PS4 (case-control enrichment)", "PM1 (hotspot / critical domain)", "PM3/BP2 (in trans / in cis)",
            "PM4/BP3 (in-frame length change)", "PP1/BS4 (segregation: zebra stats segregation)", "PP4 (phenotype specific to the gene)",
        ],
        "note": "Suggestions cover only codes that follow from numbers with ClinGen-calibrated thresholds; each still needs a look at coverage, transcript and gene context.",
    }
    lines = [f"variant: {res['variant']}"]
    for s in suggestions:
        lines.append(f"  {s['code']:<16} {s['basis']}\n                   rule: {s['rule']}")
        if s.get("conflict"):
            lines.append(f"                   CONFLICT: {s['conflict']}")
        for c in s.get("caveats") or []:
            lines.append(f"                   caveat: {c}")
    if not suggestions:
        lines.append("  no code follows from the numbers alone")
    if computed is not None:
        lines.append(f"maximum credible AF: {computed:.4g} ({maxaf['model']})")
    else:
        lines.append(f"maximum credible AF: {maxaf['why_absent']}")
    for n in res["not_assessed"]:
        lines.append(f"not assessed: {n}")
    for c in res["caveats"]:
        lines.append(f"caveat: {c}")
    lines.append("judgement codes to consider: " + "; ".join(res["not_from_data"]))
    out = Outcome(res, sources=card.sources, warnings=card.warnings, text="\n".join(lines),
                  query={"variant": args.variant, "inheritance": args.inheritance,
                         "prevalence": args.prevalence, "allelic": args.allelic, "genetic": args.genetic,
                         "penetrance": args.penetrance, "inheritance_mode": args.inheritance_mode})
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("acmg", help="ACMG/AMP: points + combining rules (classify); data-driven code suggestions (suggest)")
    a = p.add_subparsers(dest="action", metavar="<action>")
    q = a.add_parser("classify", help="classify from evidence codes")
    q.add_argument("codes", nargs="*")
    q.set_defaults(func=_classify, no_ledger=True)
    q = a.add_parser("suggest", help="codes that follow from frequency and calibrated predictors")
    q.add_argument("variant")
    q.add_argument("--assembly", choices=("GRCh38", "GRCh37"), default="GRCh38")
    q.add_argument("--inheritance", choices=("AD", "AR", "XLD", "XLR", "unknown"), default="unknown")
    # BS1 needs the Whiffin et al. 2017 maximum credible AF; these are its inputs.
    q.add_argument("--prevalence", type=float,
                   help="disease prevalence (affected fraction), e.g. 0.0004 for 1 in 2,500; needed for BS1")
    q.add_argument("--allelic", type=float,
                   help="largest share of cases one allele can explain (maximum allelic contribution)")
    q.add_argument("--genetic", type=float, default=1.0, help="share of cases due to this gene (default 1.0)")
    q.add_argument("--penetrance", type=float, default=1.0, help="penetrance of the genotype (default 1.0)")
    q.add_argument("--inheritance-mode", choices=("monoallelic", "biallelic"), default="monoallelic",
                   help="which Whiffin 2017 formula to use for the maximum credible AF (default monoallelic)")
    q.set_defaults(func=_suggest)
