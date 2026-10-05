"""zebra vcf: inspect and triage a local VCF (singleton, duo or trio). The VCF is never uploaded."""

from __future__ import annotations

import argparse

from zebra import vcf as vcf_mod
from zebra.core import Outcome, UsageError


def _inspect(args: argparse.Namespace) -> Outcome:
    return vcf_mod.inspect(args.vcf)


def _triage(args: argparse.Namespace) -> Outcome:
    if args.max_annotate < 1:
        raise UsageError("--max-annotate must be at least 1")
    genes = vcf_mod.read_gene_file(args.genes) if args.genes else None
    return vcf_mod.triage(
        args.vcf, proband=args.proband, mother=args.mother, father=args.father, sex=args.sex,
        max_af=args.max_af, min_dp=args.min_dp, min_gq=args.min_gq, assembly=args.assembly, genes=genes,
        hpo_genes=args.hpo_genes, case_dir=getattr(args, "case", None), max_annotate=args.max_annotate,
        out=args.out, hpo_terms=args.hpo, max_af_dominant=args.max_af_dominant,
    )


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("vcf", help="local VCF reanalysis: inspect, triage (only candidate variants go to VEP)")
    vs = p.add_subparsers(dest="action", metavar="<action>")

    q = vs.add_parser("inspect", help="samples, variant count, build guess, genotypes and annotations present")
    q.add_argument("vcf", help="VCF file (plain, gzip or bgzip)")
    q.set_defaults(func=_inspect)

    q = vs.add_parser("triage", help="filter by quality and inheritance, annotate survivors (VEP), rank by phenotype fit")
    q.add_argument("vcf", help="VCF file (plain, gzip or bgzip)")
    q.add_argument("--proband", required=True, help="proband sample name")
    q.add_argument("--mother", help="mother sample name")
    q.add_argument("--father", help="father sample name")
    q.add_argument("--sex", choices=("male", "female"), help="proband sex (X hemizygosity)")
    q.add_argument("--max-af", type=float, default=0.01,
                   help="max gnomAD AF for recessive/hemizygous classes (default 0.01)")
    q.add_argument("--max-af-dominant", type=float, default=None,
                   help="max gnomAD AF for de novo / heterozygous classes (default min(--max-af, 0.0001))")
    q.add_argument("--min-dp", type=int, default=10, help="min proband depth (default 10); parents need it for de novo")
    q.add_argument("--min-gq", type=int, default=20, help="min genotype quality (default 20)")
    q.add_argument("--assembly", choices=("GRCh38", "GRCh37"), help="genome build (default: from the VCF header)")
    q.add_argument("--genes", metavar="FILE", help="gene symbols to restrict to (one per line or comma-separated)")
    q.add_argument("--hpo-genes", action="store_true",
                   help="restrict to the top 200 genes for the case's HPO profile (offline ranking)")
    q.add_argument("--hpo", nargs="*", metavar="HP:0000000", help="HPO terms to use instead of the case's")
    q.add_argument("--max-annotate", type=int, default=1500,
                   help="max variants sent to VEP when there is no gene set (default 1500)")
    q.add_argument("--out", metavar="TSV", help="ranked TSV (default <case>/reports/triage-<date>.tsv)")
    q.set_defaults(func=_triage)
