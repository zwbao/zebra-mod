#!/usr/bin/env python3
"""Where the truth was in the tools' lists vs in the agent's final answer, per case.

    python3 -I tools/bench/agent/tool_vs_final.py --cases cases.jsonl --run RUN_DIR

Reads the first phenotype_rank result in each transcript (local, Monarch, PubCaseFinder lists and the
consensus) and scores each list with the same strict id join as score.py. Prints a table of
"tool had it in top-k, final answer did not" (the agent dropped it) and the reverse (the agent added it).
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from tools.bench.agent.score import ancestors_of, mondo_of, norm_id  # noqa: E402

ANC = None


def first_rank(ids, truth, tm, m):
    """Lenient match (truth or a MONDO descendant of it), as score.py's lenient column."""
    for i, x in enumerate(ids):
        pid = norm_id(x)
        if pid and (pid in truth or any(ANC(y) & tm for y in m(pid))):
            return i + 1
    return 0


def tool_lists(transcript):
    for x in transcript:
        if x.get("type") != "user":
            continue
        for c in x["message"]["content"]:
            if not (isinstance(c, dict) and c.get("type") == "tool_result"):
                continue
            s = c["content"] if isinstance(c["content"], str) else c["content"][0].get("text", "")
            try:
                r = json.loads(s)["result"]
            except Exception:
                continue
            if not isinstance(r, dict) or "per_source" not in r:
                continue
            ps = r["per_source"]
            out = {"local": [d["id"] for d in ps.get("local", {}).get("diseases", [])],
                   "monarch": [d["id"] for d in ps.get("monarch", {}).get("diseases", [])],
                   "pcf": [d["id"] for d in ps.get("pubcasefinder", {}).get("diseases_omim", [])],
                   "consensus": [h["key"] for h in r.get("consensus", [])]}
            # local list may carry OMIM and ORPHA lists separately
            for k in ("diseases_orpha", "orphanet"):
                if ps.get("local", {}).get(k):
                    out["local"] += [d["id"] for d in ps["local"][k]]
            return out
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--k", type=int, default=5)
    a = ap.parse_args()
    cases = {c["id"]: c for c in (json.loads(l) for l in open(a.cases, encoding="utf-8"))}
    per = {json.loads(l)["id"]: json.loads(l) for l in open(a.run / "per_case.jsonl", encoding="utf-8")}
    global ANC
    m, ANC = mondo_of(), ancestors_of()
    stats = defaultdict(Counter)
    rows = []
    for cid, p in per.items():
        c = cases[cid]
        truth = {norm_id(t) for t in c["truth"]} - {None}
        tm = set().union(*(m(t) for t in truth)) if truth else set()
        res = json.loads((a.run / "cases" / cid / "result.json").read_text("utf-8"))
        final = first_rank([d.get("id") for d in res["differential"]], truth, tm, m)
        tl = tool_lists(json.loads((a.run / "cases" / cid / "transcript.json").read_text("utf-8")))
        ds = c["dataset"]
        st = stats[ds]
        st["n"] += 1
        if tl is None:
            st["no_rank_call"] += 1
            continue
        loc = first_rank(tl["local"], truth, tm, m)
        st["local_top1"] += loc == 1
        st["final_top1"] += final == 1
        st[f"local_top{a.k}"] += 0 < loc <= a.k
        st[f"final_top{a.k}"] += 0 < final <= a.k
        st["dropped_from_local_top1"] += loc == 1 and final != 1
        st["promoted_to_top1"] += final == 1 and loc != 1
        st[f"lost_from_local_top{a.k}"] += 0 < loc <= a.k and not (0 < final <= a.k)
        rows.append({"id": cid, "dataset": ds, "local": loc, "final": final,
                     "monarch": first_rank(tl["monarch"], truth, tm, m), "pcf": first_rank(tl["pcf"], truth, tm, m),
                     "consensus": first_rank(tl["consensus"], truth, tm, m), "judge": p.get("judge"), "top1": (res["differential"][0] or {}).get("name"),
                     "truth_name": c["truth_names"][0], "n_present": len(c["present"])})
    (a.run / "tool_vs_final.jsonl").write_text("\n".join(json.dumps(r) for r in rows) + "\n", "utf-8")
    keys = ["n", "no_rank_call", "local_top1", "final_top1", "dropped_from_local_top1", "promoted_to_top1",
            f"local_top{a.k}", f"final_top{a.k}", f"lost_from_local_top{a.k}"]
    print("| dataset | " + " | ".join(keys) + " |")
    tot = Counter()
    for ds, st in sorted(stats.items()):
        tot.update(st)
        print(f"| {ds} | " + " | ".join(str(st[k]) for k in keys) + " |")
    print("| ALL | " + " | ".join(str(tot[k]) for k in keys) + " |")


if __name__ == "__main__":
    main()
