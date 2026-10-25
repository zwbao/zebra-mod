#!/usr/bin/env python3
"""Score agent-benchmark runs: id match (strict) and an LLM judge on names (as DeepRare did).

    python3 -I tools/bench/agent/score.py --cases cases.jsonl --run RUN_DIR [--judge] [--judge-model M]

strict  the first rank whose id is a truth id, or is joined to one by a MONDO skos:exactMatch
        (tools/bench/mapping, the pinned SSSOM); a broader or narrower class is a miss.
judge   an LLM reads the truth (names + ids) and the ranked list and returns the first rank that
        names the same disease ("same disease or a subtype/synonym of it"), as DeepRare scored
        with GPT-4o. Batches of 20 cases, blind to which system produced the list; verdicts are
        cached in RUN_DIR/judge.jsonl.
Writes RUN_DIR/scores.json: Recall@1/3/5/10 per dataset and method, Wilson 95% intervals,
cost and wall time.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

KS = (1, 3, 5, 10)


def norm_id(x):
    if not x or not isinstance(x, str):
        return None
    x = x.strip().upper().replace("ORPHANET:", "ORPHA:").replace("ORPHA_", "ORPHA:")
    m = re.match(r"(OMIM|ORPHA|MONDO):\s*0*(\d+)$", x)
    if not m:
        return None
    ns, n = m.groups()
    return f"MONDO:{int(n):07d}" if ns == "MONDO" else f"{ns}:{int(n)}"


def mondo_of():
    from tools.bench.mapping import rows_by_object
    rows = rows_by_object()
    return lambda i: {r["subject"] for r in rows.get(i, [])} | ({i} if i and i.startswith("MONDO:") else set())


def wilson(k, n, z=1.96):
    if n == 0:
        return (None, None)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (round(100 * (c - h), 1), round(100 * (c + h), 1))


JUDGE = """You grade rare-disease differential diagnoses. For each case you get the true diagnosis (one or more equivalent names/ids) and a ranked list of predictions. Return, for each case, the rank (1-based) of the FIRST prediction that refers to the same disease as the truth — the same entity, a synonym, or a more specific subtype of it (e.g. "Bardet-Biedl syndrome 1" for truth "Bardet-Biedl syndrome"); a broader group, a different subtype, or a merely related disease does not count. If none matches, return 0. Judge by meaning, not by ids alone; ids may be wrong.

Answer only with a JSON object {{"<case id>": rank, ...}}.

{cases}"""


def run_judge(items, model, cache: Path):
    done = {}
    if cache.exists():
        for l in open(cache, encoding="utf-8"):
            r = json.loads(l)
            done[r["id"]] = r["rank"]
    todo = [it for it in items if it["id"] not in done]
    for b in range(0, len(todo), 20):
        batch = todo[b:b + 20]
        text = "\n\n".join(
            f"Case {it['id']}\nTruth: " + "; ".join(it["truth_names"]) + f" ({', '.join(it['truth'])})\nPredictions:\n"
            + "\n".join(f"{i+1}. {p.get('name')} ({p.get('id')})" for i, p in enumerate(it["pred"][:10]))
            for it in batch)
        out = subprocess.run([str(Path.home() / ".local/bin/claude"), "-p", JUDGE.format(cases=text), "--model", model,
                              "--tools", "", "--no-session-persistence", "--output-format", "json"],
                             capture_output=True, text=True, cwd=str(cache.parent), timeout=600)
        res = ""
        for m in json.loads(out.stdout):
            if m.get("type") == "result":
                res = m.get("result") or ""
        got = json.loads(re.search(r"\{.*\}", res, re.S).group(0))
        with open(cache, "a", encoding="utf-8") as fh:
            for it in batch:
                r = int(got.get(it["id"], 0) or 0)
                done[it["id"]] = r
                fh.write(json.dumps({"id": it["id"], "rank": r}) + "\n")
        print(f"judged {min(b+20, len(todo))}/{len(todo)}", flush=True)
    return done


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--judge", action="store_true")
    ap.add_argument("--judge-model", default="claude-sonnet-5-5")
    a = ap.parse_args()
    cases = {c["id"]: c for c in (json.loads(l) for l in open(a.cases, encoding="utf-8"))}
    m = mondo_of()
    results, failed = [], []
    for d in sorted((a.run / "cases").iterdir()):
        if (d / "result.json").exists():
            results.append(json.loads((d / "result.json").read_text("utf-8")))
        elif (d / "failed.json").exists():
            failed.append(json.loads((d / "failed.json").read_text("utf-8")))
    per = []
    for r in results:
        c = cases[r["id"]]
        truth = {norm_id(t) for t in c["truth"]} - {None}
        tm = set().union(*(m(t) for t in truth)) if truth else set()
        strict, bad_ids = 0, 0
        for i, p in enumerate(r["differential"][:10]):
            pid = norm_id(p.get("id"))
            if p.get("id") and not pid:
                bad_ids += 1
            if pid and (pid in truth or (m(pid) & tm)) and not strict:
                strict = i + 1
        per.append({"id": r["id"], "dataset": r["dataset"], "strict": strict, "pred": r["differential"],
                    "truth": c["truth"], "truth_names": c["truth_names"], "cost": r.get("cost_usd") or 0,
                    "wall": r.get("wall_s") or 0, "turns": r.get("turns"), "malformed_ids": bad_ids,
                    "n_tools": len(r.get("tools") or [])})
    if a.judge:
        j = run_judge(per, a.judge_model, a.run / "judge.jsonl")
        for p in per:
            p["judge"] = j.get(p["id"], 0)
    out = {"n_results": len(results), "n_failed": len(failed), "failed_ids": [f["id"] for f in failed], "datasets": {}}
    for ds in sorted({p["dataset"] for p in per}) + ["ALL"]:
        rows = [p for p in per if ds == "ALL" or p["dataset"] == ds]
        n = len(rows)
        e = {"n": n, "cost_usd": round(sum(p["cost"] for p in rows), 2),
             "mean_cost_usd": round(sum(p["cost"] for p in rows) / n, 3),
             "mean_wall_s": round(sum(p["wall"] for p in rows) / n, 1),
             "mean_tool_calls": round(sum(p["n_tools"] for p in rows) / n, 1)}
        for meth in ("strict", "judge") if a.judge else ("strict",):
            for k in KS:
                hit = sum(1 for p in rows if 0 < p[meth] <= k)
                e[f"{meth}_R@{k}"] = round(100 * hit / n, 1)
                e[f"{meth}_R@{k}_ci"] = wilson(hit, n)
        out["datasets"][ds] = e
    (a.run / "scores.json").write_text(json.dumps(out, indent=1), "utf-8")
    with open(a.run / "per_case.jsonl", "w", encoding="utf-8") as fh:
        for p in per:
            fh.write(json.dumps({k: v for k, v in p.items() if k != "pred"}, ensure_ascii=False) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
