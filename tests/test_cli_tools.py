"""The mod's tool schemas (hooks/tools.ts) must say what the CLI accepts (A-P1-4, A-P1-5).

tools.ts is TypeScript and the hook test kit cannot run the CLI, so the comparison is made here:
the declarations are read from the source text and checked against the real argparse parser.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from zebra import case as case_mod
from zebra.cli import build_parser

TOOLS_TS = (Path(__file__).resolve().parent.parent / "hooks" / "tools.ts").read_text("utf-8")


def _subparsers(parser: argparse.ArgumentParser) -> dict:
    return next(a for a in parser._actions if isinstance(a, argparse._SubParsersAction)).choices


def _flags(parser: argparse.ArgumentParser) -> set:
    return {o[2:] for a in parser._actions for o in a.option_strings
            if o.startswith("--") and o not in ("--help", "--json", "--case")}


def test_a_p1_4_stats_flags_match_the_cli_parser_exactly():
    block = re.search(r"const STATS_FLAGS[^{]*\{(.*?)\n\}", TOOLS_TS, re.S).group(1)
    declared = {m.group(1): set(re.findall(r"'([\w-]+)'", m.group(2)))
                for m in re.finditer(r"(\w+):\s*\[(.*?)\]", block, re.S)}
    methods = _subparsers(_subparsers(build_parser())["stats"])
    assert set(declared) == set(methods), "rare_stats methods differ from `zebra stats`"
    for name, sub in methods.items():
        assert declared[name] == _flags(sub), f"stats {name}: tool {sorted(declared[name])} vs CLI {sorted(_flags(sub))}"


def test_a_p1_5_every_variant_kind_the_tool_offers_is_accepted_by_the_case(tmp_path):
    m = re.search(r"kind: \{\s*type: 'string',\s*enum: \[([^\]]*)\]", TOOLS_TS)
    kinds = re.findall(r"'([\w-]+)'", m.group(1))
    assert set(kinds) == set(case_mod.VARIANT_KINDS)
    for alias, kind in case_mod.KIND_ALIASES.items():
        assert case_mod.check_variant(kind=alias, description="x")["kind"] == kind
