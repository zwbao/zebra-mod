"""zebra aso: splice-switching antisense feasibility screen for one variant."""

from __future__ import annotations

import argparse

from zebra import aso
from zebra.core import Outcome


def _aso(args: argparse.Namespace) -> Outcome:
    return aso.screen(
        args.variant,
        assembly=args.assembly,
        lengths=aso.parse_lengths(args.lengths),
        event=args.event,
        min_delta=args.min_delta,
        top=args.top,
        uniqueness=args.uniqueness,
        distance=args.distance,
        stride=args.stride,
        body_stride=args.body_stride,
        timeout=args.timeout,
        with_precedents=not args.no_precedents,
    )


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "aso",
        help="splice-switching antisense feasibility screen: which windows on the pre-mRNA would cover the "
             "aberrant splice site this variant is predicted to create?",
    )
    p.add_argument("variant", help="chrom-pos-ref-alt, transcript HGVS (NM_...:c....) or an rsID")
    p.add_argument("--assembly", choices=("GRCh38", "GRCh37"), default="GRCh38")
    p.add_argument("--lengths", default=f"{aso.DEFAULT_LENGTHS[0]}-{aso.DEFAULT_LENGTHS[1]}",
                   help=f"target lengths as LO-HI or a comma-separated list, between {aso.LENGTH_LIMITS[0]} and "
                        f"{aso.LENGTH_LIMITS[1]} (default "
                        f"{aso.DEFAULT_LENGTHS[0]}-{aso.DEFAULT_LENGTHS[1]})")
    p.add_argument("--event", choices=aso.EVENTS, default="auto",
                   help="auto (default) reads the event off the SpliceAI gains: a gained acceptor AND donor a "
                        "plausible exon apart is a cryptic exon, one alone is a cryptic acceptor or donor. "
                        "pseudoexon / cryptic_acceptor / cryptic_donor force that reading and stop when the "
                        "prediction does not support it. exon_skip ignores the gains and targets the authentic "
                        "splice sites of the exon the variant lies in (or nearest)")
    p.add_argument("--min-delta", type=float, default=aso.MIN_DELTA,
                   help=f"SpliceAI delta that counts as supporting a splice effect (default {aso.MIN_DELTA}, the "
                        f"ClinGen SVI threshold; Walker et al. 2023). A gain above {aso.DELTA_FLOOR} but below "
                        f"this is still screened, flagged as being in the uninformative band")
    p.add_argument("--top", type=int, default=aso.DEFAULT_TOP,
                   help=f"how many candidates to report per screen (default {aso.DEFAULT_TOP}); the full count is "
                        f"always reported")
    p.add_argument("--uniqueness", action="store_true",
                   help="check each shortlisted target against the human genome with one NCBI BLAST search "
                        "(asynchronous; adds roughly 30-120 s). Without it every candidate says uniqueness was "
                        "not checked and the ranking leaves that component out")
    p.add_argument("--distance", type=int, default=500,
                   help="SpliceAI/Pangolin window around the variant, 1-10000 (default 500)")
    p.add_argument("--stride", type=int, default=aso.SITE_STRIDE,
                   help=f"step in bases between successive windows over a splice site (default {aso.SITE_STRIDE}; "
                        f"1 enumerates every offset)")
    p.add_argument("--body-stride", type=int, default=aso.BODY_STRIDE,
                   help=f"step in bases between windows tiling a cryptic exon's body (default {aso.BODY_STRIDE})")
    p.add_argument("--timeout", type=float, default=240.0,
                   help="per-model wall-clock limit for the splicing prediction, in seconds (default 240)")
    p.add_argument("--no-precedents", action="store_true",
                   help="skip the Europe PMC lookup of the individualized-antisense precedent and eligibility "
                        "records (they are retrieved live; nothing is cited from memory either way)")
    p.set_defaults(func=_aso)
