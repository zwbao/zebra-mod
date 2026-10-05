"""zebra edit: base-editing feasibility screen for reverting a variant to the reference."""

from __future__ import annotations

import argparse
import re
from typing import Optional, Tuple

from zebra import editing
from zebra.core import Outcome, UsageError

WINDOW_RE = re.compile(r"^(\d{1,2})\s*[-:]\s*(\d{1,2})$")


def _window(raw: Optional[str]) -> Optional[Tuple[int, int]]:
    """None means "use the window of the chemistry this correction needs"."""
    if raw is None or not str(raw).strip():
        return None
    m = WINDOW_RE.match(str(raw).strip())
    if not m:
        raise UsageError("--window takes LO-HI protospacer positions, e.g. 4-8")
    return int(m.group(1)), int(m.group(2))


def _edit(args: argparse.Namespace) -> Outcome:
    return editing.screen(
        args.variant,
        assembly=args.assembly,
        pam=args.pam,
        window=_window(args.window),
        annotate_bystanders=args.annotate_bystanders,
    )


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "edit",
        help="base-editing feasibility screen: could an ABE/CBE revert this SNV, and with which protospacers?",
    )
    p.add_argument("variant", help="chrom-pos-ref-alt (SNV), transcript HGVS or an rsID")
    p.add_argument("--assembly", choices=("GRCh38", "GRCh37"), default="GRCh38")
    p.add_argument("--pam", choices=tuple(editing.PAMS), default="NGG",
                   help="NGG = SpCas9 (default); NG = engineered relaxed-PAM variants, lower efficiency")
    p.add_argument("--window", default=None,
                   help="editing window as protospacer positions LO-HI, counted from the PAM-distal end. "
                        "Default: the window of the chemistry the correction needs — CBE 4-8 (BE3/BE4; Komor "
                        "et al. 2016) or ABE 4-7 (ABE7.10; Gaudelli et al. 2017). Use 3-9 for ABE8e "
                        "(Richter et al. 2020).")
    p.add_argument("--annotate-bystanders", action="store_true",
                   help="put each bystander edit through VEP to see whether it would change coding")
    p.set_defaults(func=_edit)
