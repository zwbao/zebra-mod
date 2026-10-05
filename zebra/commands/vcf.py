"""zebra vcf: inspect and triage a local VCF (singleton, duo or trio). The VCF is never uploaded.

This module also registers `zebra cnv`, the sibling entry for the result forms
that are not an SNV/indel list: a CMA/CNV-seq interval, an exon-level deletion
or duplication, a copy-number count (SMN1/SMN2) and a repeat expansion. They
live next to `vcf` because they answer the same question from a different
report: what did the laboratory actually find. The analysis is in `zebra.cnv`.
"""

from __future__ import annotations

import argparse
from typing import Any, Dict

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


def _cnv(args: argparse.Namespace) -> Outcome:
    from zebra import cnv as cnv_mod

    out = cnv_mod.card(args.result, assembly=args.assembly, gene=args.gene, copies=args.copies,
                       inheritance=args.inheritance, method=args.method, related=args.related)
    if args.record:
        from zebra import case as case_mod

        target = getattr(args, "case", None)
        if not target:
            raise UsageError("--record needs a case: pass --case <dir> (or set ZEBRA_CASE)")
        try:
            fields: Dict[str, Any] = cnv_mod.case_fields(out.result, method=args.method)
            entry = case_mod.add_variant(target, **fields)
        except case_mod.CaseError as err:
            raise UsageError(str(err)) from None
        out.result = dict(out.result)
        out.result["recorded_in_case"] = {"id": entry["id"], "kind": entry["kind"], "case": target}
        if out.text:
            out.text += f"\n[recorded in the case as {entry['id']} ({entry['kind']})]"
    return out


def register(sub: argparse._SubParsersAction) -> None:
    c = sub.add_parser("cnv", help="a CNV/CMA interval, an exon-level del/dup, a copy-number or repeat-expansion "
                                   "result: genes spanned, dosage sensitivity, frame, ACMG CNV inputs")
    c.add_argument("result", help='e.g. "chr15:23123715-28193120 loss" | "arr[GRCh38] 22q11.21(18648855_21800471)x1" '
                                 '| "DMD exon 45-50 deletion" | "NM_004006.3:c.6439-?_7309+?del" '
                                 '| "SMN1 exon 7 copy number 0" | "FMR1 CGG 230"')
    c.add_argument("--assembly", choices=("GRCh38", "GRCh37"), default="GRCh38",
                   help="build of the coordinates (default GRCh38; an ISCN string's own build wins)")
    c.add_argument("--gene", help="gene symbol when the input does not name one (exon forms)")
    c.add_argument("--copies", type=int, help="copy number the report gives, when it is not in the string")
    c.add_argument("--inheritance", choices=("de_novo", "maternal", "paternal", "biparental", "unknown"),
                   help="ACMG CNV section 5 input, as the family study found it")
    c.add_argument("--method", help="how it was measured (CMA, CNV-seq, MLPA, ddPCR, repeat-primed PCR, …)")
    c.add_argument("--related", nargs="*", metavar="RESULT",
                   help='further copy-number results from the same report, e.g. "SMN2 copy number 2"')
    c.add_argument("--record", action="store_true", help="also record the finding in the case (--case)")
    c.set_defaults(func=_cnv)

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
