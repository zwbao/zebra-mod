#!/usr/bin/env python3
"""Where the agent's answers come from, on a batch of cases run with and without the plugin.

    python3 -I tools/bench/agent/decomposition.py --cases cases.jsonl --batch tools/bench/results/agent/batch1_ids.txt \
        --ranker-run runs/main --arm "Opus 5.5" runs/main runs/noplugin \
        --arm "DeepSeek-V4.1-Flash" runs/deepseek runs/deepseek_noplugin --out decomposition.json

Each batch case is classified by whether the local phenotype ranker alone (local list of the first phenotype_rank
result in --ranker-run) and the model alone (the arm's run without the plugin; "no_answer" when that run never
answered) ranked the correct diagnosis first, and we count how often the agent (the arm's run with the plugin) did,
all with the MONDO rule (score.py's lenient column). Also: how representative the batch is of each full dataset,
whether the identifiers in the final answers appear verbatim in a tool result of the same run, and process measures
(turns, zebra-mod tool calls by tool, skills loaded, time). Local only, no model calls.
"""

from __future__ import annotations

import argparse
import json
import re
import statistics
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import tools.bench.agent.tool_vs_final as tvf  # noqa: E402
from tools.bench.agent.score import ancestors_of, mondo_of, norm_id  # noqa: E402

ZEBRA = "mcp__zebra-mod__"


def per_case(run: Path) -> dict:
    return {r["id"]: r for r in (json.loads(l) for l in open(run / "per_case.jsonl", encoding="utf-8"))}


def transcript(run: Path, cid: str) -> list:
    return json.loads((run / "cases" / cid / "transcript.json").read_text("utf-8"))


def tool_text(t: list) -> str:
    out = []
    for m in t:
        if m.get("type") == "user" and isinstance(m["message"].get("content"), list):
            for b in m["message"]["content"]:
                if b.get("type") == "tool_result":
                    c = b.get("content")
                    out.append(c if isinstance(c, str) else " ".join(x.get("text", "") for x in c if isinstance(x, dict)))
    return "\n".join(out)


def ranker_first(cases: dict, run: Path, ids: list, m) -> dict:
    out = {}
    for cid in ids:
        truth = {norm_id(t) for t in cases[cid]["truth"]} - {None}
        tm = set().union(*(m(t) for t in truth)) if truth else set()
        tl = tvf.tool_lists(transcript(run, cid))
        out[cid] = bool(tl) and tvf.first_rank(tl["local"], truth, tm, m) == 1
    return out


def provenance(run: Path, ids: list) -> dict:
    total = verbatim = 0
    for cid in ids:
        res = json.loads((run / "cases" / cid / "result.json").read_text("utf-8"))
        txt = tool_text(transcript(run, cid))
        for e in res.get("differential") or []:
            x = (e.get("id") or "").strip()
            mm = re.match(r"(OMIM|ORPHA|MONDO)[:_](\d+)$", x, re.I)
            if not x:
                continue
            total += 1
            verbatim += bool(mm and re.search(rf"\b{mm.group(1)}[:_]{mm.group(2)}\b", txt, re.I))
    return {"identifiers": total, "verbatim_in_tool_output": verbatim, "pct": round(100 * verbatim / max(1, total), 2)}


def process(run: Path, ids: list, per: dict) -> dict:
    calls, cases_with, skills = Counter(), Counter(), Counter()
    for cid in ids:
        seen = set()
        for m in transcript(run, cid):
            if m.get("type") != "assistant":
                continue
            for b in m["message"].get("content", []):
                if b.get("type") != "tool_use":
                    continue
                if b["name"].startswith(ZEBRA):
                    n = b["name"][len(ZEBRA):]
                    calls[n] += 1
                    seen.add(n)
                elif b["name"] == "Skill":
                    skills[b["input"].get("skill", "?")] += 1
        cases_with.update(seen)
    answered = [cid for cid in ids if cid in per]
    return {"answered": len(answered), "mean_turns": round(statistics.mean(per[c]["turns"] or 0 for c in answered), 1),
            "zebra_calls_per_case": round(sum(calls.values()) / len(ids), 2), "calls_by_tool": dict(calls.most_common()),
            "cases_with_tool": dict(cases_with.most_common()), "skills_loaded": dict(skills.most_common()),
            "median_seconds": round(statistics.median(per[c]["wall"] for c in answered))}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--batch", type=Path, required=True, help="whitespace-separated case ids")
    ap.add_argument("--ranker-run", type=Path, required=True)
    ap.add_argument("--arm", nargs=3, action="append", metavar=("LABEL", "WITH_PLUGIN", "WITHOUT_PLUGIN"), required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    cases = {c["id"]: c for c in (json.loads(l) for l in open(a.cases, encoding="utf-8"))}
    ids = sorted(a.batch.read_text().split())
    m = mondo_of()
    tvf.ANC = ancestors_of()
    rank = ranker_first(cases, a.ranker_run, ids, m)
    out = {"n": len(ids), "ranker_alone_R1": sum(rank.values()), "arms": {}}
    for label, with_run, without_run in a.arm:
        agent, bare = per_case(Path(with_run)), per_case(Path(without_run))
        state = {c: ("right" if bare[c]["lenient"] == 1 else "wrong") if c in bare else "no_answer" for c in ids}
        hit = {c: agent[c]["lenient"] == 1 for c in ids}
        cells = {}
        for r in (True, False):
            for s in ("right", "wrong", "no_answer"):
                cell = [c for c in ids if rank[c] == r and state[c] == s]
                if cell:
                    cells[f"ranker_{'right' if r else 'wrong'}|model_{s}"] = {"n": len(cell), "agent_right": sum(hit[c] for c in cell)}
        misses = [c for c in ids if not hit[c]]
        out["arms"][label] = {
            "model_alone_R1": sum(s == "right" for s in state.values()), "agent_R1": sum(hit.values()),
            "either_source_R1": sum(rank[c] or state[c] == "right" for c in ids), "cells": cells,
            "agent_misses": len(misses), "misses_where_neither_source_right": sum(not rank[c] and state[c] != "right" for c in misses),
            "identifiers": provenance(Path(with_run), ids), "process_with_plugin": process(Path(with_run), ids, agent)}
    # how representative the batch is: the ranker run's agent (DeepRare judge) and ranker (MONDO), batch vs the rest
    main_per = per_case(a.ranker_run)
    rest = sorted(set(main_per) - set(ids))
    rank_rest = ranker_first(cases, a.ranker_run, rest, m)
    rank_all = {**rank, **rank_rest}
    rep = {}
    for ds in sorted({cases[c]["dataset"] for c in rest}):
        for g, sel in (("batch", ids), ("rest", rest)):
            s = [c for c in sel if cases[c]["dataset"] == ds]
            rep[f"{ds}|{g}"] = {"n": len(s), "agent_R1_deeprare": round(100 * sum(main_per[c]["deeprare"] == 1 for c in s) / len(s), 1),
                                "ranker_R1_mondo": round(100 * sum(rank_all[c] for c in s) / len(s), 1)}
    out["batch_vs_rest"] = rep
    out["ranker_run_identifiers_all_cases"] = provenance(a.ranker_run, sorted(main_per))
    a.out.write_text(json.dumps(out, indent=1) + "\n", "utf-8")
    for label, r in out["arms"].items():
        print(f"{label}: ranker {out['ranker_alone_R1']}, model alone {r['model_alone_R1']}, agent {r['agent_R1']} of {len(ids)}; "
              f"misses {r['agent_misses']} ({r['misses_where_neither_source_right']} where neither source was right)")
        for k, v in r["cells"].items():
            print(f"  {k:32s} n={v['n']:3d}  agent right {v['agent_right']:3d} ({100 * v['agent_right'] / v['n']:.0f}%)")


if __name__ == "__main__":
    main()
