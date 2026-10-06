#!/usr/bin/env python3
"""Choose ranking parameters on the development half only, then score the choice once on held-out.

    python3 -I tools/bench/tune.py dev        # every candidate on the dev half (writes results/tuning.json)
    python3 -I tools/bench/tune.py heldout NAME [NAME ...]   # the named candidates on the held-out half

Candidates are in CANDIDATES below (the last ones were added after the first development pass,
still on the development half only). Selection rule: the candidate with the best development
top-10 (OMIM-only) overall; it is then run on the held-out half, once, and kept as zebra's
default only if its held-out top-10 beats the shipped code's (docs/BENCHMARK.md reports both).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.bench.common import RESULTS, read_jsonl, work_dir  # noqa: E402
from tools.bench import run_local  # noqa: E402

CANDIDATES = {
    "resnik": {"method": "resnik"},
    "resnik_excl0": {"method": "resnik", "params": {"excluded_weight": 0.0}},
    "resnik_excl25": {"method": "resnik", "params": {"excluded_weight": 0.25}},
    "resnik_excl50": {"method": "resnik", "params": {"excluded_weight": 0.5}},
    "lr": {"method": "lr"},
    "lr_fu75": {"method": "lr", "params": {"f_unknown": 0.75}},
    "lr_fu50": {"method": "lr", "params": {"f_unknown": 0.5}},
    "lr_cap50": {"method": "lr", "params": {"excluded_cap": 0.5}},
    "lr_cap99": {"method": "lr", "params": {"excluded_cap": 0.99}},
    "lr_unexp-1": {"method": "lr", "params": {"unexplained": -1.0}},
    "lr_unexp-2": {"method": "lr", "params": {"unexplained": -2.0}},
    "lr_floor10": {"method": "lr", "params": {"f_floor": 0.1}},
    "lr_noexcl": {"method": "lr", "params": {"excluded_cap": 0.0}},
    # added after the first dev pass showed every excluded-term penalty costing accuracy: only an
    # absent feature the disease has (very) frequently, at least 80% of patients, is evidence against it
    "lr_minf80": {"method": "lr", "params": {"excluded_min_f": 0.8}},
    "lr_minf80_cap50": {"method": "lr", "params": {"excluded_min_f": 0.8, "excluded_cap": 0.5}},
    "lr_minf95": {"method": "lr", "params": {"excluded_min_f": 0.95}},
    "lr_noexcl_fu50": {"method": "lr", "params": {"excluded_cap": 0.0, "f_unknown": 0.5}},
    "resnik_excl10": {"method": "resnik", "params": {"excluded_weight": 0.1}},
    "resnik_minf80": {"method": "resnik", "params": {"resnik_excluded_min_f": 0.8}},
    "resnik_minf80_w25": {"method": "resnik", "params": {"resnik_excluded_min_f": 0.8, "excluded_weight": 0.25}},
    # CP1-16: ties broken by how frequent the matched annotations are, before specificity
    "resnik_excl0_tbf": {"method": "resnik", "params": {"excluded_weight": 0.0, "tiebreak_frequency": 1.0}},
}


def top10(label: str) -> dict:
    rows = [r for r in read_jsonl(work_dir() / f"local_{label}.jsonl.gz") if not r.get("error")]
    n = len(rows)
    out = {"n": n}
    for key in ("rank_omim", "rank_all", "gene_rank"):
        sub = [r for r in rows if key != "gene_rank" or r.get("has_gene")]
        for k in (1, 10):
            out[f"{key}_top{k}"] = round(sum(1 for r in sub if r.get(key) and r[key] <= k) / len(sub), 4) if sub else None
    return out


def run(name: str, split: str) -> dict:
    cand = CANDIDATES[name]
    opts = {"method": cand["method"]}
    if cand.get("params"):
        opts["params"] = cand["params"]
    label = f"{split}_{name}"
    run_local.main(["--label", label, "--split", split, "--options", json.dumps(opts), "--workers", "12"])
    return {"candidate": name, "split": split, "options": opts, **top10(label)}


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    split = argv[0]
    names = argv[1:] or list(CANDIDATES)
    path = RESULTS / "tuning.json"
    log = json.loads(path.read_text("utf-8")) if path.exists() else {"rule": __doc__.strip().splitlines()[-3:], "runs": []}
    for name in names:
        row = run(name, split)
        log["runs"] = [r for r in log["runs"] if not (r["candidate"] == name and r["split"] == split)] + [row]
        path.write_text(json.dumps(log, indent=1) + "\n", "utf-8")
        print(json.dumps(row), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
