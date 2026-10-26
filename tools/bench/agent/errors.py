#!/usr/bin/env python3
"""Error analysis of an agent-benchmark run: which cases were missed, and the measurable reasons.

    python3 -I tools/bench/agent/errors.py --cases cases.jsonl --run RUN_DIR [--packets N]

Needs score.py (per_case.jsonl) and tool_vs_final.py (tool_vs_final.jsonl) run first. A miss is a case whose
truth is not first in the final answer (lenient: truth or a MONDO descendant). For each miss: where the truth was
in the final top 10 and in the tools' lists, whether the truth disease has HPO annotations at all, how many
present terms the case had, and what was ranked first instead. Writes RUN_DIR/errors.json (the tables) and
RUN_DIR/miss_packets.jsonl (one compact record per miss, for a reader to classify by failure mode).
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from tools.bench.agent.score import mondo_of, norm_id  # noqa: E402
from tools.bench.common import hpo_files  # noqa: E402


def annotated_ids():
    ids = set()
    for line in open(Path(hpo_files()["dir"]) / "phenotype.hpoa", encoding="utf-8"):
        if line[:1] != "#" and not line.startswith("database_id"):
            ids.add(norm_id(line.split("\t", 1)[0].replace("ORPHA:", "ORPHA:")))
    return ids - {None}


def nbucket(n):
    return "1-3" if n <= 3 else "4-6" if n <= 6 else "7-10" if n <= 10 else ">10"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--run", type=Path, required=True)
    a = ap.parse_args()
    cases = {c["id"]: c for c in (json.loads(l) for l in open(a.cases, encoding="utf-8"))}
    per = {json.loads(l)["id"]: json.loads(l) for l in open(a.run / "per_case.jsonl", encoding="utf-8")}
    tvf = {json.loads(l)["id"]: json.loads(l) for l in open(a.run / "tool_vs_final.jsonl", encoding="utf-8")}
    m, ann = mondo_of(), annotated_ids()
    mondo_ann = set().union(*(m(i) for i in ann))

    def has_annotation(truth):
        for t in truth:
            t = norm_id(t)
            if t and (t in ann or (m(t) & mondo_ann)):
                return True
        return False

    out = {"by_dataset": {}, "by_terms": {}, "confusions": {}, "diseases": {}}
    packets = []
    strata = defaultdict(Counter)
    for ds in sorted({p["dataset"] for p in per.values()}):
        rows = [p for p in per.values() if p["dataset"] == ds]
        st = Counter(n=len(rows))
        dis = defaultdict(lambda: [0, 0])
        conf = Counter()
        for p in rows:
            c, t = cases[p["id"]], tvf.get(p["id"], {})
            final = t.get("final", p["lenient"])
            hit = final == 1
            name = c["truth_names"][0].split("/")[-1].strip()
            dis[name][0] += 1
            dis[name][1] += hit
            nb = nbucket(len(c["present"]))
            strata[(ds, nb)]["n"] += 1
            strata[(ds, nb)]["hit"] += hit
            if hit:
                continue
            st["miss"] += 1
            best = min([x for x in (t.get("local", 0), t.get("monarch", 0), t.get("pcf", 0)) if x > 0] or [0])
            st["final_2-10" if 0 < final <= 10 else "final_absent"] += 1
            st["tool_had_it_first" if best == 1 else "tool_2-10" if 0 < best <= 10 else
               "tool_11-50" if best else "no_tool_had_it"] += 1
            annotated = has_annotation(c["truth"])
            st["truth_unannotated"] += not annotated
            st["terms_1-3"] += len(c["present"]) <= 3
            conf[(name, (t.get("top1") or "")[:80])] += 1
            res = json.loads((a.run / "cases" / p["id"] / "result.json").read_text("utf-8"))
            tr = json.loads((a.run / "cases" / p["id"] / "transcript.json").read_text("utf-8"))
            final_text = next((x.get("result") or "" for x in tr if x.get("type") == "result"), "")
            packets.append({"id": p["id"], "dataset": ds, "truth": c["truth"], "truth_names": c["truth_names"],
                            "present": [h["label"] for h in c["present"]],
                            "excluded": [h["label"] for h in c["excluded"]],
                            "final_rank_of_truth": final, "tool_ranks": {k: t.get(k, 0) for k in ("local", "monarch", "pcf")},
                            "truth_has_hpo_annotation": annotated,
                            "differential": [f"{d.get('name')} ({d.get('id')})" for d in res["differential"][:10]],
                            "reasoning_excerpt": final_text.split("```json")[0][-2500:]})
        out["by_dataset"][ds] = dict(st)
        out["diseases"][ds] = sorted(([k, v[0], v[1]] for k, v in dis.items() if v[0] - v[1] > 0),
                                     key=lambda r: -(r[1] - r[2]))[:12]
        out["confusions"][ds] = [[k[0], k[1], n] for k, n in conf.most_common(10)]
    for (ds, nb), s in sorted(strata.items()):
        out["by_terms"].setdefault(ds, {})[nb] = [s["n"], round(100 * s["hit"] / s["n"], 1)]
    (a.run / "errors.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), "utf-8")
    with open(a.run / "miss_packets.jsonl", "w", encoding="utf-8") as fh:
        for pk in packets:
            fh.write(json.dumps(pk, ensure_ascii=False) + "\n")
    print(json.dumps(out["by_dataset"], indent=1))
    print(json.dumps(out["by_terms"]))
    print(len(packets), "packets")


if __name__ == "__main__":
    main()
