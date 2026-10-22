"""zebra repeat: repeat sizes read by rule (offline) — D4Z4 units from a 4q35 fragment, DMPK CTG classes."""

from __future__ import annotations

import argparse

from zebra import repeats
from zebra.commands.stats import finite


def _d4z4(args: argparse.Namespace):
    return repeats.d4z4(args.fragment_kb, digest=args.digest, haplotype=args.haplotype,
                        methylation_pct=args.methylation_pct, locus=args.locus)


def _dmpk(args: argparse.Namespace):
    return repeats.dmpk(args.ctg)


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("repeat", help="repeat sizes read by rule: D4Z4 units from a 4q35 fragment (FSHD1), "
                                      "DMPK CTG classes (DM1)")
    rs = p.add_subparsers(dest="repeat_cmd", required=True)
    d = rs.add_parser("d4z4", help="D4Z4 units from the fragment size a laboratory reports")
    d.add_argument("--fragment-kb", type=finite, required=True, help="the contracted fragment, kb")
    d.add_argument("--digest", default="ecori", help="EcoRI (default), EcoRI/BlnI or EcoRI/HindIII "
                                                     "(a HindIII-only or BlnI-only digest is refused)")
    d.add_argument("--haplotype", help="4qA or 4qB, as reported")
    d.add_argument("--locus", help="4q or 10q, only when the laboratory states which chromosome the allele is on")
    d.add_argument("--methylation-pct", type=finite, help="D4Z4 methylation %%, as reported (recorded, not interpreted)")
    d.set_defaults(func=_d4z4)
    m = rs.add_parser("dmpk", help="DM1 class from a DMPK CTG repeat count")
    m.add_argument("--ctg", type=int, required=True, help="CTG repeats, as reported (the larger allele)")
    m.set_defaults(func=_dmpk)
