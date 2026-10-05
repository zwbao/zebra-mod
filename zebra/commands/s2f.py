"""zebra s2f: sequence-to-function predictions for one variant, model by model."""

from __future__ import annotations

import argparse
from typing import List, Optional

from zebra import s2f
from zebra.core import Outcome, UsageError


def _models(raw: Optional[str]) -> Optional[List[str]]:
    if not raw:
        return None
    names = [m.strip().lower().replace("-", "_") for m in raw.split(",") if m.strip()]
    if not names:
        raise UsageError(f"--models needs at least one of {', '.join(s2f.MODELS)}")
    return names


def _predict(args: argparse.Namespace) -> Outcome:
    out = s2f.predict(
        args.variant,
        assembly=args.assembly,
        models=_models(args.models),
        distance=args.distance,
        workspace=getattr(args, "workspace", None),
        ontology=args.ontology,
        timeout=args.timeout,
        case=getattr(args, "case", None),
        mask=args.mask,
        evo2_window=args.evo2_window,
    )
    out.text = s2f.render(out.result)
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("s2f", help="sequence-to-function models for a variant (splicing, regulation, likelihood, constraint)")
    a = p.add_subparsers(dest="action", metavar="<action>")
    q = a.add_parser(
        "predict",
        help="run SpliceAI/Pangolin (Broad lookup) and, when available, AlphaGenome/Evo 2/GPN-MSA via the s2f CLI",
    )
    q.add_argument("variant", help="chrom-pos-ref-alt, transcript HGVS (NM_...:c....) or an rsID")
    q.add_argument("--assembly", choices=("GRCh38", "GRCh37"), default="GRCh38")
    q.add_argument("--models", help=f"comma-separated subset of {', '.join(s2f.MODELS)} (default: the splice models plus every s2f model ready here)")
    q.add_argument("--distance", type=int, default=500, help="SpliceAI/Pangolin window around the variant, 1-10000 (default 500)")
    q.add_argument("--mask", type=int, choices=(0, 1), default=0,
                   help="0 (default) = raw delta scores, the scale the ClinGen thresholds (Walker et al. 2023) were "
                        "fitted on, so keep it for variant interpretation. 1 = masked: gains at an existing site and "
                        "losses at a non-site are zeroed, which reads more cleanly but voids those thresholds. "
                        "Pangolin is sent the same setting, where 1 means its own --mask.")
    q.add_argument("--ontology", help="UBERON/CL CURIE for AlphaGenome's tissue or cell type, e.g. UBERON:0000955 (brain)")
    q.add_argument("--timeout", type=float, default=s2f.DEFAULT_TIMEOUT, help=f"per-model wall-clock limit in seconds (default {s2f.DEFAULT_TIMEOUT:g})")
    q.add_argument("--evo2-window", type=int, default=2048, help="Evo 2 variant window in bases (default 2048)")
    q.add_argument("--workspace", help="where s2f run directories go (default: <case>/evidence/s2f, else ~/.cache/zebra-mod/s2f-runs)")
    q.set_defaults(func=_predict)
