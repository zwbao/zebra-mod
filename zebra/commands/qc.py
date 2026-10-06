"""zebra qc: genomic QC of a local (family) VCF — sex, relatedness, ROH, Mendelian errors / UPD, mosaic de novo.

Computed on this machine; nothing is sent anywhere. The analysis is in `zebra.qc`.
"""

from __future__ import annotations

import argparse

from zebra import qc as qc_mod
from zebra.core import Outcome, UsageError


def _qc(args: argparse.Namespace) -> Outcome:
    if args.min_dp < 0 or args.min_gq < 0:
        raise UsageError("--min-dp and --min-gq must be >= 0")
    if args.max_samples < 1:
        raise UsageError("--max-samples must be at least 1")
    if not args.max_seconds > 0:
        raise UsageError("--max-seconds must be positive")
    return qc_mod.run(args.vcf, ped=args.ped, proband=args.proband, mother=args.mother, father=args.father,
                      sex=args.sex, assembly=args.assembly, min_dp=args.min_dp, min_gq=args.min_gq,
                      max_samples=args.max_samples, max_seconds=args.max_seconds,
                      progress=qc_mod.stderr_progress("zebra qc"))


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("qc", help="family VCF QC on this machine: sex check, KING-robust relatedness (sample swap / "
                                  "non-paternity), runs of homozygosity, Mendelian errors and UPD, mosaic de novo")
    p.add_argument("vcf", help="VCF file (plain, gzip or bgzip), joint-called for a family")
    p.add_argument("--ped", metavar="FILE", help="PED file (family, id, father, mother, sex 1/2/0, phenotype 1/2/0): "
                                                 "names the trios and the expected relationships to check")
    p.add_argument("--proband", help="proband sample (without --ped)")
    p.add_argument("--mother", help="mother sample (without --ped)")
    p.add_argument("--father", help="father sample (without --ped)")
    p.add_argument("--sex", choices=("male", "female"), help="the proband's stated sex, checked against genotypes")
    p.add_argument("--assembly", choices=("GRCh38", "GRCh37"), help="genome build (default: from the header)")
    p.add_argument("--min-dp", type=int, default=10, help="min depth per genotype (default 10)")
    p.add_argument("--min-gq", type=int, default=20, help="min genotype quality (default 20)")
    p.add_argument("--max-samples", type=int, default=qc_mod.MAX_SAMPLES,
                   help=f"samples assessed at most (default {qc_mod.MAX_SAMPLES}; kinship is pairwise)")
    p.add_argument("--max-seconds", type=float, default=qc_mod.QC_MAX_SECONDS,
                   help=f"wall-time cap on the scan (default {qc_mod.QC_MAX_SECONDS:g}); a capped scan says so")
    p.set_defaults(func=_qc)
