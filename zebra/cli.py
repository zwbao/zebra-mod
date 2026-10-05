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
MESSAGE_CAP = 4_000  # an error message never carries a whole oversized input back


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
    """The longest list anywhere under `node`, by item count, with its dotted path.

    Item count, not serialised size: a list's serialisation always contains its
    children's, so ranking by size could only ever pick the outermost list, and
    trimming that drops whole records (a variant, a disease) when trimming an
    inner list would have freed the same bytes.
    """
    best: Tuple[Optional[list], str, int, int] = (None, "", 0, 0)
    if isinstance(node, list):
        if node:
            best = (node, path, len(node), len(_dump(node, None)))
        for i, item in enumerate(node):
            got = _longest_list(item, f"{path}[{i}]")
            if (got[2], got[3]) > (best[2], best[3]):
                best = got
    elif isinstance(node, dict):
        for key, item in node.items():
            got = _longest_list(item, f"{path}.{key}")
            if (got[2], got[3]) > (best[2], best[3]):
                best = got
    return best


def _subparser(parser: argparse.ArgumentParser, name: str) -> Optional[argparse.ArgumentParser]:
    for act in parser._actions:  # noqa: SLF001
        if isinstance(act, argparse._SubParsersAction):  # noqa: SLF001
            return act.choices.get(name)
    return None


def _longest_string(node: Any) -> Optional[Tuple[Any, Any, int]]:
    """The longest string value under `node`, as (holder, key, length), so it can be cut in place."""
    best: Optional[Tuple[Any, Any, int]] = None
    stack = [node]
    while stack:
        current = stack.pop()
        items = current.items() if isinstance(current, dict) else (
            enumerate(current) if isinstance(current, list) else ())
        for key, value in items:
            if isinstance(value, str):
                if best is None or len(value) > best[2]:
                    best = (current, key, len(value))
            elif isinstance(value, (dict, list)):
                stack.append(value)
    return best


def _dump(envelope: Any, indent: Optional[int] = 1, ascii_only: bool = False) -> str:
    return json.dumps(envelope, ensure_ascii=ascii_only, indent=indent, allow_nan=False)


def _emit(text: str, stream=None) -> None:
    """Write one line to stdout, surviving a stdout that cannot encode it.

    A non-UTF-8 stdout (PYTHONIOENCODING, a cp1252 console, C-locale coercion
    disabled) would otherwise raise UnicodeEncodeError *inside* the error
    handler and leave the caller with empty stdout and a traceback — which is
    exactly what the envelope contract exists to prevent. This toolkit answers
    queries in Chinese, so it is first-class input.
    """
    out = stream or sys.stdout
    try:
        out.write(text + "\n")
        out.flush()
        return
    except UnicodeEncodeError:
        pass
    encoding = getattr(out, "encoding", None) or "ascii"
    out.write(text.encode(encoding, "backslashreplace").decode(encoding, "replace") + "\n")
    out.flush()


def _emit_json(envelope: Dict[str, Any], indent: Optional[int] = 1) -> None:
    """Write a JSON envelope, re-encoding non-ASCII as \\uXXXX if stdout cannot take it."""
    try:
        _emit(_dump(envelope, indent))
    except UnicodeEncodeError:  # pragma: no cover - _emit handles it, this is the belt
        _emit(_dump(envelope, indent, ascii_only=True))


def _trim(envelope: Dict[str, Any], budget: int) -> str:
    """Serialise, trimming the longest list in `result` until it fits. The envelope stays valid JSON.

    The number of items kept is found by bisection, so as much of the result
    survives as the budget allows, and the warning names exactly how many
    items were dropped.
    """
    text = _dump(envelope)
    if budget <= 0 or len(text) <= budget:
        return text
    warnings = envelope.setdefault("warnings", [])
    # the bulk is normally `result`, but an oversized query string or a long
    # list of source records must be cut too: a consumer that hard-cuts the
    # body at its own limit would otherwise be handed invalid JSON
    for key in ("query", "result"):
        node = envelope.get(key)
        big = _longest_string(node)
        if big is not None and big[2] > max(200, budget // 4):
            holder, field, length = big
            keep = max(100, budget // 4)
            holder[field] = str(holder[field])[:keep] + f"… [{length - keep} more characters, cut to stay inside " \
                                                        f"the {budget}-character limit]"
            warnings.append(f"{key}.{field} was {length} characters and was cut to {keep} to stay under "
                            f"{budget} characters")
            text = _dump(envelope)
            if len(text) <= budget:
                return text
    for _ in range(40):
        target, path, size, bytes_ = _longest_list(envelope.get("result"))
        sources = envelope.get("sources")
        if isinstance(sources, list) and sources:
            source_bytes = len(_dump(sources, None))
            if target is None or source_bytes > bytes_:
                target, path, size = sources, "sources", len(sources)
        if target is None or size == 0:
            break
        whole = list(target)
        slot = len(warnings)
        what = "provenance record(s)" if path == "sources" else "items"
        warnings.append(f"result trimmed to stay under {budget} characters: {size} of {size} {what} dropped "
                        f"from {path} (the list is now empty)")  # the longest wording, so the budget holds
        low, high = 0, size - 1  # at least one item must go, or there is no progress
        while low < high:
            mid = (low + high + 1) // 2
            target[:] = whole[:mid]
            if len(_dump(envelope)) <= budget:
                low = mid
            else:
                high = mid - 1
        target[:] = whole[:low]
        dropped = size - low
        warnings[slot] = (f"{'sources' if path == 'sources' else 'result'} trimmed to stay under {budget} "
                          f"characters: {dropped} of {size} {what} dropped from {path}"
                          + ("" if low else " (the list is now empty)")
                          + (" — the ledger ids still name every source" if path == "sources" else ""))
        text = _dump(envelope)
        if len(text) <= budget:
            return text
    if len(text) > budget:
        kept = {k: v for k, v in envelope.items() if k != "result"}
        had_lists = _longest_list(envelope.get("result"))[0] is not None
        kept["warnings"] = list(kept.get("warnings") or []) + [
            f"result dropped entirely: it does not fit in {budget} characters"
            + (" even after trimming its lists" if had_lists else " and holds no list that could be trimmed")
        ]
        kept["result"] = None
        text = _dump(kept)
    return text


def _nonfinite_warning(bad: List[str]) -> str:
    return ("not a finite number, reported as null: " + "; ".join(bad[:5])
            + (f" (+{len(bad) - 5} more)" if len(bad) > 5 else ""))


def _envelope(command: str, outcome: Outcome, ledger_ids: List[str]) -> Dict[str, Any]:
    """Warnings, ledger ids and sources come before `result`: they are what the doctrine cites."""
    bad: List[str] = []
    result = _finite(outcome.result, "result", bad)
    query = _finite(outcome.query, "query", bad)
    sources = _finite(outcome.sources, "sources", bad)
    warnings = list(outcome.warnings)
    if bad:
        warnings.append(_nonfinite_warning(bad))
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


def _cap(message: str, limit: Optional[int] = None) -> str:
    """Keep an error message inside the size budget: a 200 kB query must not be echoed whole."""
    if limit is None:
        budget = _max_bytes()
        limit = MESSAGE_CAP if budget <= 0 else max(500, min(MESSAGE_CAP, budget // 2))
    text = str(message)
    if len(text) <= limit:
        return text
    return text[:limit] + f"… [{len(text) - limit} more characters, cut to stay inside the size budget]"


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
    # known before parsing, so parse errors can be JSON too; argparse accepts any
    # unambiguous prefix of --json, and the two must never disagree
    want_json = any(t.startswith("--j") and "--json".startswith(t) for t in tokens)
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
            _emit_json({"ok": code == 0, "zebra": __version__, "command": "zebra", "query": {}, "warnings": [],
                        "ledger": [], "sources": [], "result": {"stdout": text}} if code == 0 else
                       {"ok": False, "zebra": __version__, "command": "zebra",
                        "error": {"type": "UsageError", "message": text.strip() or f"exit {code}"}})
        else:
            _emit(text.rstrip("\n"))
        return code
    except KeyboardInterrupt:
        return _fail(want_json, "zebra", "Interrupted", "interrupted while reading the arguments", 130)
    except Exception as err:  # noqa: BLE001 - the envelope must survive anything
        return _internal(want_json, "zebra", err)

    if not getattr(args, "func", None):
        named = getattr(args, "command", None)
        sub = _subparser(parser, named) if named else None
        help_text = (sub or parser).format_help().strip()
        if want_json:
            return _fail(True, named or "zebra", "UsageError", help_text, 2 if sub else 1)
        _emit(help_text)
        return 2 if sub else 1
    command = " ".join(p for p in [args.command, getattr(args, "action", None)] if p)
    try:
        outcome: Outcome = args.func(args)
    except UsageError as err:
        return _fail(getattr(args, "json", False), command, "UsageError", str(err), 2)
    except SourceError as err:
        return _fail(getattr(args, "json", False), command, "SourceError", f"{err.source}: {err.message}", 3,
                     source=err.source, status=err.status)
    except KeyboardInterrupt:
        return _fail(getattr(args, "json", False), command, "Interrupted",
                     "interrupted before the answer was complete", 130)
    except SystemExit as err:  # a handler (or a library) called sys.exit
        return _internal(getattr(args, "json", False), command,
                         RuntimeError(f"the command exited with status {err.code!r} instead of returning a result"))
    except Exception as err:  # noqa: BLE001
        return _internal(getattr(args, "json", False), command, err)

    ledger_ids: List[str] = []
    if args.case and outcome.sources and not getattr(args, "no_ledger", False):
        from zebra import case as case_mod

        try:
            # the ledger is a durable artifact: it gets the sanitised copies, so
            # a non-finite number cannot make a row unparseable by strict JSON
            ledger_ids = case_mod.append_ledger(args.case, command, _finite(outcome.query, "query", []),
                                                _finite(outcome.sources, "sources", []))
        except (OSError, case_mod.CaseError) as err:
            outcome.warnings.append(f"evidence ledger not written: {err}")

    try:
        if args.json:
            _emit(_trim(_envelope(command, outcome, ledger_ids), _max_bytes()))
        else:
            if outcome.text is not None:
                _emit(outcome.text)
            else:
                nonfinite: List[str] = []
                _emit(_dump(_finite(outcome.result, "result", nonfinite), indent=2))
                if nonfinite:
                    outcome.warnings.append(_nonfinite_warning(nonfinite))
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
        _emit_json({"ok": False, "zebra": __version__, "command": command,
                    "error": {"type": kind, "message": _cap(message), **extra}}, indent=None)
    else:
        _emit(f"zebra {command}: {message}", sys.stderr)
    return code


def _internal(want_json: bool, command: str, err: BaseException) -> int:
    tail = "".join(traceback.format_exception(type(err), err, err.__traceback__))[-TRACEBACK_TAIL:]
    if want_json:
        _emit_json({"ok": False, "zebra": __version__, "command": command,
                    "error": {"type": "InternalError", "exception": type(err).__name__, "message": _cap(str(err)),
                              "traceback_tail": tail}}, indent=None)
    else:
        _emit(f"zebra {command}: internal error ({type(err).__name__}: {err})", sys.stderr)
        _emit(tail, sys.stderr)
    return 70


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
