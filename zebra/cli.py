"""`zebra` command line: one subcommand per module in zebra/commands.

Each module in zebra/commands defines `register(subparsers)`, adding its
parsers with `set_defaults(func=handler)`; a handler takes the parsed args and
returns an Outcome. `--json` prints the envelope the mod's tools read;
otherwise the Outcome's text (or the JSON, pretty) is printed.

Envelope contract (what the mod's tools rely on): with `--json`, stdout is
always one valid JSON object — for results, for usage and argparse errors, for
`--help`/`--version`, and for unexpected exceptions. Non-finite numbers become
null plus a warning (JSON has no NaN), and an oversized result is trimmed with
a warning that names how many items were dropped, so the envelope stays valid.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib
import io
import json
import math
import os
import pkgutil
import sys
import traceback
from typing import Any, Dict, List, Optional, Tuple

from zebra import __version__
from zebra.core import Outcome, UsageError
from zebra.http import SourceError

# The mod cuts a tool result at 60,000 characters; trimming here instead keeps
# the JSON valid and the warnings, sources and ledger ids intact.
DEFAULT_MAX_BYTES = 60_000
TRACEBACK_TAIL = 2_000


class _ArgvError(Exception):
    """argparse wanted to exit with a usage error; raised so `main` can shape the envelope."""

    def __init__(self, message: str, usage: str):
        super().__init__(message)
        self.message = message
        self.usage = usage


class _Parser(argparse.ArgumentParser):
    """argparse that raises instead of writing usage to stderr and calling sys.exit."""

    def error(self, message: str):  # noqa: D102 - argparse hook
        raise _ArgvError(message, self.format_usage().strip())


class _SubParser(_Parser):
    """Every subcommand also takes --json and --case, wherever they are typed."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.add_argument("--json", action="store_true", default=argparse.SUPPRESS, help=argparse.SUPPRESS)
        self.add_argument("--case", default=argparse.SUPPRESS, help=argparse.SUPPRESS)


def _no_flag_value(inner):
    """Reject a value that starts with '-' where a value is expected (flag smuggling).

    `zebra variant -- --help` then returns a UsageError envelope instead of
    printing help, and `--` still works for values that legitimately begin with
    a dash-less token.
    """

    def convert(text: str):
        if isinstance(text, str) and len(text) > 1 and text[0] == "-" and not (text[1].isdigit() or text[1] == "."):
            raise argparse.ArgumentTypeError(
                f"{text!r} looks like a flag, not a value; if it really is the value, zebra still refuses it here "
                "(pass the value without a leading '-')"
            )
        return inner(text) if inner is not None else text

    convert.__name__ = getattr(inner, "__name__", "str")
    return convert


def _guard_positionals(parser: argparse.ArgumentParser) -> None:
    for act in parser._actions:  # noqa: SLF001 - argparse has no public walk
        if isinstance(act, argparse._SubParsersAction):  # noqa: SLF001
            for sub in act.choices.values():
                _guard_positionals(sub)
            continue
        if act.option_strings or act.type in (int, float) or act.choices:
            continue
        act.type = _no_flag_value(act.type)


def build_parser() -> argparse.ArgumentParser:
    parser = _Parser(
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
    _guard_positionals(parser)
    return parser


# ------------------------------------------------------------- JSON hygiene

def _finite(value: Any, path: str, bad: List[str]) -> Any:
    """Replace NaN/Infinity (which JSON cannot hold) with null, naming where they were."""
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            bad.append(f"{path} was {value!r}")
            return None
        return value
    if isinstance(value, dict):
        return {k: _finite(v, f"{path}.{k}" if path else str(k), bad) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite(v, f"{path}[{i}]", bad) for i, v in enumerate(value)]
    return value


def _longest_list(node: Any, path: str = "result") -> Tuple[Optional[list], str, int, int]:
    """The longest list anywhere under `node` (by serialised length), with its dotted path."""
    best: Tuple[Optional[list], str, int, int] = (None, "", 0, 0)
    if isinstance(node, list):
        if node:
            best = (node, path, len(node), len(_dump(node, None)))
        for i, item in enumerate(node):
            got = _longest_list(item, f"{path}[{i}]")
            if got[3] > best[3]:
                best = got
    elif isinstance(node, dict):
        for key, item in node.items():
            got = _longest_list(item, f"{path}.{key}")
            if got[3] > best[3]:
                best = got
    return best


def _dump(envelope: Dict[str, Any], indent: Optional[int] = 1) -> str:
    return json.dumps(envelope, ensure_ascii=False, indent=indent, allow_nan=False)


def _trim(envelope: Dict[str, Any], budget: int) -> str:
    """Serialise, trimming the longest list in `result` until it fits. The envelope stays valid JSON."""
    text = _dump(envelope)
    if budget <= 0 or len(text) <= budget:
        return text
    for _ in range(40):
        target, path, size, _ = _longest_list(envelope.get("result"))
        if target is None or size == 0:
            break
        # trim proportionally; the loop corrects itself if the first cut is not enough
        keep = min(size - 1, max(1 if size > 1 else 0, int(size * budget / len(text)) - 1))
        drop = size - keep
        del target[keep:]
        envelope.setdefault("warnings", []).append(
            f"result trimmed to stay under {budget} characters: {drop} of {size} items dropped from {path}"
            + ("" if keep else " (the list is now empty)")
        )
        text = _dump(envelope)
        if len(text) <= budget:
            return text
    if len(text) > budget:
        kept = {k: v for k, v in envelope.items() if k != "result"}
        kept["warnings"] = list(kept.get("warnings") or []) + [
            f"result dropped entirely: it does not fit in {budget} characters even after trimming its lists"
        ]
        kept["result"] = None
        text = _dump(kept)
    return text


def _envelope(command: str, outcome: Outcome, ledger_ids: List[str]) -> Dict[str, Any]:
    """Warnings, ledger ids and sources come before `result`: they are what the doctrine cites."""
    bad: List[str] = []
    result = _finite(outcome.result, "result", bad)
    query = _finite(outcome.query, "query", bad)
    sources = _finite(outcome.sources, "sources", bad)
    warnings = list(outcome.warnings)
    if bad:
        warnings.append("not a finite number, reported as null: " + "; ".join(bad[:5])
                        + (f" (+{len(bad) - 5} more)" if len(bad) > 5 else ""))
    return {
        "ok": True,
        "zebra": __version__,
        "command": command,
        "query": query,
        "warnings": warnings,
        "ledger": ledger_ids,
        "sources": sources,
        "result": result,
    }


def _max_bytes() -> int:
    raw = os.environ.get("ZEBRA_MAX_BYTES")
    if raw is None:
        return DEFAULT_MAX_BYTES
    try:
        return max(0, int(raw))
    except ValueError:
        return DEFAULT_MAX_BYTES


# ------------------------------------------------------------------- main

def main(argv: Optional[List[str]] = None) -> int:
    tokens = list(sys.argv[1:] if argv is None else argv)
    want_json = "--json" in tokens  # known before parsing, so parse errors can be JSON too
    parser = build_parser()
    buffer = io.StringIO()
    try:
        with contextlib.redirect_stdout(buffer) if want_json else contextlib.nullcontext():
            args = parser.parse_args(argv)
    except _ArgvError as err:
        return _fail(want_json, "zebra", "UsageError", f"{err.message}\n{err.usage}", 2)
    except SystemExit as err:  # --help / --version / an argparse exit we did not intercept
        code = 0 if err.code is None else int(err.code)
        text = buffer.getvalue()
        if want_json:
            print(_dump({"ok": code == 0, "zebra": __version__, "command": "zebra", "warnings": [],
                         "result": {"stdout": text}} if code == 0 else
                        {"ok": False, "zebra": __version__, "command": "zebra",
                         "error": {"type": "UsageError", "message": text.strip() or f"exit {code}"}}))
        else:
            sys.stdout.write(text)
        return code
    except Exception as err:  # noqa: BLE001 - the envelope must survive anything
        return _internal(want_json, "zebra", err)

    if not getattr(args, "func", None):
        if want_json:
            return _fail(True, "zebra", "UsageError", parser.format_help().strip(), 1)
        parser.print_help()
        return 1
    command = " ".join(p for p in [args.command, getattr(args, "action", None)] if p)
    try:
        outcome: Outcome = args.func(args)
    except UsageError as err:
        return _fail(getattr(args, "json", False), command, "UsageError", str(err), 2)
    except SourceError as err:
        return _fail(getattr(args, "json", False), command, "SourceError", f"{err.source}: {err.message}", 3,
                     source=err.source, status=err.status)
    except KeyboardInterrupt:
        return 130
    except Exception as err:  # noqa: BLE001
        return _internal(getattr(args, "json", False), command, err)

    ledger_ids: List[str] = []
    if args.case and outcome.sources and not getattr(args, "no_ledger", False):
        from zebra import case as case_mod

        try:
            ledger_ids = case_mod.append_ledger(args.case, command, outcome.query, outcome.sources)
        except (OSError, case_mod.CaseError) as err:
            outcome.warnings.append(f"evidence ledger not written: {err}")

    try:
        if args.json:
            print(_trim(_envelope(command, outcome, ledger_ids), _max_bytes()))
        else:
            if outcome.text is not None:
                print(outcome.text)
            else:
                print(json.dumps(_finite(outcome.result, "result", []), ensure_ascii=False, indent=2, allow_nan=False))
            if outcome.warnings:
                print("\nwarnings:", file=sys.stderr)
                for w in outcome.warnings:
                    print(f"  - {w}", file=sys.stderr)
            if ledger_ids:
                print(f"\n[ledger {ledger_ids[0]}..{ledger_ids[-1]} → {args.case}/evidence/ledger.jsonl]", file=sys.stderr)
    except Exception as err:  # noqa: BLE001 - rendering must not produce a traceback either
        return _internal(getattr(args, "json", False), command, err)
    return 0


def _fail(want_json: bool, command: str, kind: str, message: str, code: int, **extra) -> int:
    if want_json:
        print(_dump({"ok": False, "zebra": __version__, "command": command,
                     "error": {"type": kind, "message": message, **extra}}, indent=None))
    else:
        print(f"zebra {command}: {message}", file=sys.stderr)
    return code


def _internal(want_json: bool, command: str, err: BaseException) -> int:
    tail = "".join(traceback.format_exception(type(err), err, err.__traceback__))[-TRACEBACK_TAIL:]
    if want_json:
        print(_dump({"ok": False, "zebra": __version__, "command": command,
                     "error": {"type": "InternalError", "exception": type(err).__name__, "message": str(err),
                               "traceback_tail": tail}}, indent=None))
    else:
        print(f"zebra {command}: internal error ({type(err).__name__}: {err})", file=sys.stderr)
        print(tail, file=sys.stderr)
    return 70


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
