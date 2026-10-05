"""zebra acmg: ACMG/AMP points and combining rules; code suggestions from data."""

from __future__ import annotations

import argparse
from typing import Any, Dict

from zebra import acmg
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
    if res.get("note"):
        lines.append(f"note: {res['note']}")
    for w in res["warnings"]:
        lines.append(f"warning: {w}")
    lines.append("research-grade: confirm with an accredited laboratory before clinical use")
    return Outcome(res, text="\n".join(lines), query={"codes": args.codes})


def _suggest(args: argparse.Namespace) -> Outcome:
    from zebra.sources import variant as variant_src

    card = variant_src.card(args.variant, assembly=args.assembly)
    data: Dict[str, Any] = card.result.get("acmg_inputs", {})
    suggestions = acmg.suggest_from_data(data, inheritance=None if args.inheritance == "unknown" else args.inheritance)
    res = {
        "variant": card.result.get("variant"),
        "inputs": data,
        "suggested": suggestions,
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
    if not suggestions:
        lines.append("  no code follows from the numbers alone")
    lines.append("judgement codes to consider: " + "; ".join(res["not_from_data"]))
    out = Outcome(res, sources=card.sources, warnings=card.warnings, text="\n".join(lines),
                  query={"variant": args.variant, "inheritance": args.inheritance})
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
    q.set_defaults(func=_suggest)
