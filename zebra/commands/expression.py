"""zebra expression: where a gene is expressed (GTEx), and which RNA test that makes possible."""

from __future__ import annotations

import argparse
from typing import Any, Dict, List

from zebra.core import Outcome, UsageError
from zebra.sources import gtex


def render(r: Dict[str, Any]) -> str:
    lines: List[str] = [f"{r['gene']}  {r.get('gencode_id') or '-'}  ({r['dataset']}, GENCODE "
                        f"{r.get('gencode_version') or '-'}; {r['unit']})"]
    if r.get("note") and not r.get("tissues"):
        lines.append(r["note"])
        return "\n".join(lines)
    lines.append(f"highest-expressing tissues ({len(r['tissues'])} of {r['tissues_total']}):")
    for t in r["tissues"]:
        tpm = t.get("median_tpm")
        lines.append(f"  {t['tissue']:<38} {('-' if tpm is None else format(tpm, '>10.3f'))}  "
                     f"{str(t.get('uberon') or '-'):<16} {t.get('band') or '-'}")
    lines.append("")
    lines.append("tissues an RNA test could use:")
    for t in r["rna_test_tissues"]:
        tpm = t.get("median_tpm")
        lines.append(f"  {t['tissue']:<38} {('-' if tpm is None else format(tpm, '>10.3f'))}  "
                     f"{str(t.get('uberon') or '-'):<16} {t.get('band') or 'not measured here'}")
        lines.append(f"  {'':<38} {t.get('why_it_matters') or t.get('note') or ''}")
    lines.append("")
    lines.append(f"which RNA test: {r.get('rna_test_hint')}")
    if r.get("ontology_hint"):
        lines.append(f"ontology term: {r['ontology_hint']}")
    lines.append("")
    for c in r.get("caveats") or []:
        lines.append(f"! {c}")
    return "\n".join(lines)


def _expression(args: argparse.Namespace) -> Outcome:
    gene = (args.gene or "").strip()
    if not gene:
        raise UsageError("give a gene symbol, e.g. zebra expression CFTR")
    out = gtex.top_tissues(gene, n=args.top, dataset=args.dataset)
    out.text = render(out.result)
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "expression",
        help="median expression per tissue for a gene (GTEx), with the blood / fibroblast / muscle values that "
             "decide which RNA test is feasible and the tissue's ontology id",
    )
    p.add_argument("gene", help="gene symbol, e.g. CFTR")
    p.add_argument("--top", type=int, default=10, help="how many highest-expressing tissues to list (default 10)")
    p.add_argument("--dataset", choices=tuple(sorted(gtex.DATASETS)), default=gtex.DEFAULT_DATASET,
                   help=f"GTEx release (default {gtex.DEFAULT_DATASET}); the GENCODE id is resolved per release, so "
                        f"the two give slightly different medians")
    p.set_defaults(func=_expression)
