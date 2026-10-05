"""`zebra` command line: one subcommand per module in zebra/commands.

Each module in zebra/commands defines `register(subparsers)`, adding its
parsers with `set_defaults(func=handler)`; a handler takes the parsed args and
returns an Outcome. `--json` prints the envelope the mod's tools read;
otherwise the Outcome's text (or the JSON, pretty) is printed.
"""

from __future__ import annotations

import argparse
import importlib
import json
import os
import pkgutil
import sys
from typing import List, Optional

from zebra import __version__
from zebra.core import Outcome, UsageError
from zebra.http import SourceError


class _SubParser(argparse.ArgumentParser):
    """Every subcommand also takes --json and --case, wherever they are typed."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.add_argument("--json", action="store_true", default=argparse.SUPPRESS, help=argparse.SUPPRESS)
        self.add_argument("--case", default=argparse.SUPPRESS, help=argparse.SUPPRESS)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="zebra",
        description="Rare-disease research toolkit (zebra-mod). Every answer carries its sources.",
    )
    parser.add_argument("--version", action="version", version=f"zebra {__version__}")
    parser.add_argument("--json", action="store_true", help="print the JSON envelope (what the mod's tools read)")
    parser.add_argument(
        "--case",
        default=os.environ.get("ZEBRA_CASE"),
        help="case directory: sources are appended to its evidence ledger (default $ZEBRA_CASE)",
    )
    sub = parser.add_subparsers(dest="command", metavar="<command>", parser_class=_SubParser)
    import zebra.commands as commands_pkg

    for info in sorted(pkgutil.iter_modules(commands_pkg.__path__), key=lambda m: m.name):
        if info.name.startswith("_"):
            continue
        module = importlib.import_module(f"zebra.commands.{info.name}")
        register = getattr(module, "register", None)
        if register:
            register(sub)
    return parser


def _envelope(command: str, outcome: Outcome, ledger_ids: List[str]) -> dict:
    return {
        "ok": True,
        "zebra": __version__,
        "command": command,
        "query": outcome.query,
        "result": outcome.result,
        "sources": outcome.sources,
        "warnings": outcome.warnings,
        "ledger": ledger_ids,
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 1
    command = " ".join(p for p in [args.command, getattr(args, "action", None)] if p)
    try:
        outcome: Outcome = args.func(args)
    except UsageError as err:
        return _fail(args, command, "UsageError", str(err), 2)
    except SourceError as err:
        return _fail(args, command, "SourceError", f"{err.source}: {err.message}", 3, source=err.source, status=err.status)
    except KeyboardInterrupt:
        return 130

    ledger_ids: List[str] = []
    if args.case and outcome.sources and not getattr(args, "no_ledger", False):
        from zebra import case as case_mod

        try:
            ledger_ids = case_mod.append_ledger(args.case, command, outcome.query, outcome.sources)
        except (OSError, case_mod.CaseError) as err:
            outcome.warnings.append(f"evidence ledger not written: {err}")

    if args.json:
        print(json.dumps(_envelope(command, outcome, ledger_ids), ensure_ascii=False, indent=1))
    else:
        if outcome.text is not None:
            print(outcome.text)
        else:
            print(json.dumps(outcome.result, ensure_ascii=False, indent=2))
        if outcome.warnings:
            print("\nwarnings:", file=sys.stderr)
            for w in outcome.warnings:
                print(f"  - {w}", file=sys.stderr)
        if ledger_ids:
            print(f"\n[ledger {ledger_ids[0]}..{ledger_ids[-1]} → {args.case}/evidence/ledger.jsonl]", file=sys.stderr)
    return 0


def _fail(args: argparse.Namespace, command: str, kind: str, message: str, code: int, **extra) -> int:
    if getattr(args, "json", False):
        print(json.dumps({"ok": False, "zebra": __version__, "command": command, "error": {"type": kind, "message": message, **extra}}, ensure_ascii=False))
    else:
        print(f"zebra {command}: {message}", file=sys.stderr)
    return code


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
