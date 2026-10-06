"""zebra report export: a Markdown report as HTML, Word and PDF that a family or a clinic can open.

For the person who runs zebra-mod on a family's behalf (a clinician, a genetic
counsellor, a volunteer): the family letter and the visit-preparation sheet go
to the family as files, with no terminal on their side.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import List, Optional

from zebra import case as case_mod
from zebra import report_export
from zebra.core import Outcome, UsageError

FORMATS = ("html", "docx", "pdf")


def _formats(text: str) -> List[str]:
    out = [f.strip().lower() for f in text.split(",") if f.strip()]
    out = ["docx" if f in ("word", "doc") else f for f in out]
    bad = [f for f in out if f not in FORMATS]
    if bad or not out:
        raise argparse.ArgumentTypeError(f"formats are {', '.join(FORMATS)} (comma-separated), not {text!r}")
    return out


def _cases_for(md: Path, explicit: Optional[str]) -> List[str]:
    """Every case whose identifiers apply: the case folder the report sits in, and --case."""
    found: List[str] = []
    for parent in [md.resolve().parent, *md.resolve().parents][:4]:
        if (parent / "case.json").is_file():
            found.append(str(parent))
            break
    if explicit and os.path.realpath(explicit) not in {os.path.realpath(f) for f in found}:
        found.append(explicit)
    return found


def _export(args: argparse.Namespace) -> Outcome:
    md = Path(args.report).expanduser()
    if md.suffix.lower() not in (".md", ".markdown", ".txt"):
        raise UsageError(f"{md.name}: export takes the Markdown report (.md) the report skill wrote")
    cases = _cases_for(md, args.case)
    identifiers: List[str] = []
    warnings: List[str] = []
    for case_dir in cases:
        try:
            identifiers += list(case_mod.load(case_dir)["privacy"].get("identifiers") or [])
        except case_mod.CaseError as err:
            raise UsageError(f"cannot read the case's protected identifiers ({err}); nothing was exported") from None
    if not cases:
        warnings.append("no case found for this report: checked for ID numbers only, not for the patient's name")
    try:
        got = report_export.export(str(md), args.to, out_dir=args.out, title=args.title,
                                   identifiers=identifiers, allow_identifiers=args.allow_identifiers)
    except FileNotFoundError as err:
        raise UsageError(str(err)) from None
    except PermissionError as err:
        raise UsageError(str(err)) from None
    except ValueError as err:
        raise UsageError(str(err)) from None
    warnings.extend(got.pop("warnings"))
    lines = [f"{f['format']}: {f['path']} ({f['bytes']:,} bytes)" for f in got["files"]]
    lines += [f"{s['format']}: not made — {s['reason']}" for s in got["skipped"]]
    return Outcome(got, warnings=warnings, text="\n".join(lines),
                   query={"report": os.path.basename(str(md)), "to": args.to})


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("report", help="export a Markdown report as HTML, Word and PDF (Chinese typography included)")
    rs = p.add_subparsers(dest="action", metavar="<action>")
    q = rs.add_parser("export", help="write <report>.html/.docx/.pdf next to the report (or in --out)")
    q.add_argument("report", help="the Markdown report, e.g. <case>/reports/family-letter-2026-10-06.md")
    q.add_argument("--to", type=_formats, default=["docx", "pdf"],
                   help="html, docx (Word), pdf — comma-separated (default docx,pdf)")
    q.add_argument("--out", help="directory for the files (default: next to the report)")
    q.add_argument("--title", help="document title (default: the report's first heading)")
    q.add_argument("--allow-identifiers", action="store_true",
                   help="export even if the report contains the case's protected identifiers (an internal copy only)")
    q.set_defaults(func=_export, no_ledger=True)
    p.set_defaults(func=lambda a: (_ for _ in ()).throw(UsageError("zebra report export <report.md> [--to docx,pdf]")))
