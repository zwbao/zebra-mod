"""zebra cnv: a CNV/CMA interval, an exon-level deletion or duplication, a copy-number count
(SMN1/SMN2) or a repeat expansion — the result forms that are not an SNV/indel list.

They answer the same question as `zebra vcf` from a different report: what did
the laboratory actually find. The analysis is in `zebra.cnv`.
"""

from __future__ import annotations

import argparse
from typing import Any, Dict, Optional

from zebra.core import Outcome, UsageError


def _case_sex(target: Optional[str]) -> Optional[str]:
    """The proband's sex recorded in the case profile, when it is male or female (B-P2-9)."""
    if not target:
        return None
    from zebra import case as case_mod

    try:
        data = case_mod.load(target)
    except (case_mod.CaseError, OSError, ValueError):
        return None
    from zebra import cnv as cnv_mod

    return cnv_mod._sex_norm(str((data.get("proband") or {}).get("sex") or ""))


def _cnv(args: argparse.Namespace) -> Outcome:
    from zebra import cnv as cnv_mod

    target = getattr(args, "case", None)
    sex = args.sex
    case_sex = _case_sex(target)
    sex_note = None
    if sex and case_sex and sex != case_sex:
        sex_note = f"--sex says {sex} and the case profile says {case_sex}: --sex was used"
    if not sex and case_sex:
        sex = case_sex
        sex_note = f"the proband's sex ({case_sex}) was taken from the case profile"
    out = cnv_mod.card(args.result, assembly=args.assembly, gene=args.gene, copies=args.copies,
                       inheritance=args.inheritance, method=args.method, related=args.related,
                       smn2_copies=args.smn2_copies, sex=sex)
    if sex_note and out.result.get("kind") == "cnv" and out.result.get("region", {}).get("chrom") in ("X", "Y"):
        out.warnings.append(sex_note)
    if args.record:
        from zebra import case as case_mod

        if not target:
            raise UsageError("--record needs a case: pass --case <dir> (or set ZEBRA_CASE)")
        fields: Dict[str, Any] = cnv_mod.case_fields(out.result, method=args.method,
                                                     inheritance=args.inheritance or out.result.get("inheritance")
                                                     or (out.query or {}).get("inheritance"))
        try:
            existing = next((v for v in case_mod.load(target).get("variants") or []
                             if isinstance(v, dict) and v.get("kind") == fields.get("kind")
                             and cnv_mod.record_identity(v) == cnv_mod.record_identity(fields)), None)
            entry = existing or case_mod.add_variant(target, **fields)
        except case_mod.CaseError as err:
            # the analysis stands; only the write failed, and the result says so
            out.result = dict(out.result)
            out.result["recorded_in_case"] = {"recorded": False, "case": target, "error": str(err)}
            out.warnings.append(f"NOT recorded in the case: {err}")
            if out.text:
                out.text += f"\n[NOT recorded in the case: {err}]"
            return out
        out.result = dict(out.result)
        out.result["recorded_in_case"] = {"recorded": True, "id": entry["id"], "kind": entry["kind"],
                                          "case": target, "already_recorded": existing is not None}
        if out.text:
            out.text += (f"\n[already in the case as {entry['id']} ({entry['kind']}): not added again]" if existing
                         else f"\n[recorded in the case as {entry['id']} ({entry['kind']})]")
    return out


def register(sub: argparse._SubParsersAction) -> None:
    c = sub.add_parser("cnv", help="a CNV/CMA interval, an exon-level del/dup, a copy-number or repeat-expansion "
                                   "result: genes spanned, ClinGen region and gene dosage, frame, SMN1/SMN2 and "
                                   "repeat-size meaning (cited), ACMG CNV inputs")
    c.add_argument("result", help='e.g. "chr15:23123715-28193120 loss" | "arr[GRCh38] 22q11.21(18648855_21800471)x1" '
                                 '| "DMD exon 45-50 deletion" | "NM_004006.3:c.6439-?_7309+?del" '
                                 '| "SMN1 0 copies, SMN2 3 copies" | "FMR1 CGG 230"')
    c.add_argument("--assembly", choices=("GRCh38", "GRCh37"), default="GRCh38",
                   help="build of the coordinates (default GRCh38; a build written in the result wins)")
    c.add_argument("--gene", help="gene symbol when the input does not name one (exon forms)")
    c.add_argument("--copies", type=int, help="copy number the report gives, when it is not in the string")
    c.add_argument("--inheritance", choices=("de_novo", "maternal", "paternal", "biparental", "unknown"),
                   help="ACMG CNV section 5 input, as the family study found it")
    c.add_argument("--method", help="how it was measured (CMA, CNV-seq, MLPA, ddPCR, repeat-primed PCR, …)")
    c.add_argument("--related", nargs="*", metavar="RESULT",
                   help='further copy-number results from the same report, e.g. "SMN2 copy number 2"')
    c.add_argument("--smn2-copies", dest="smn2_copies", type=int,
                   help="SMN2 copy number reported with an SMN1 result (the main modifier)")
    c.add_argument("--sex", choices=("male", "female"),
                   help="proband's sex, which decides whether an X/Y copy number is a loss or a gain "
                        "(default: the case profile's)")
    c.add_argument("--record", action="store_true", help="also record the finding in the case (--case)")
    c.set_defaults(func=_cnv)
