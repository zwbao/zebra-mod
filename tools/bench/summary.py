#!/usr/bin/env python3
"""Collect every benchmark number into tools/bench/results/benchmark.json (what docs/BENCHMARK.md cites).

    python3 -I tools/bench/summary.py --before base010 --after LABEL [--remote-before FILE --remote-after FILE]

--before / --after are run_local labels over ALL cases (the shipped 0.1.0 ranker and the new
default); the remote files are remote_eval.py --json outputs. Also stores the tuning log, the
provenance of every input (phenopacket-store tag and sha256, MONDO SSSOM tag and sha256, HPO
release) and a compact per-case table of the --after run (results/local_cases.jsonl.gz).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.bench import mapping, metrics  # noqa: E402
from tools.bench.common import RESULTS, git_sha, read_jsonl, work_dir, write_jsonl  # noqa: E402
from zebra import __version__  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", required=True)
    ap.add_argument("--after", required=True)
    ap.add_argument("--remote-before")
    ap.add_argument("--remote-after")
    ap.add_argument("--reps", type=int, default=1000)
    args = ap.parse_args(argv)
    store = json.loads((work_dir() / "store.json").read_text("utf-8"))
    out = {
        "zebra_version": __version__, "git": git_sha(),
        "dataset": {k: store["store"].get(k) for k in ("repo", "tag", "published_at", "asset_url", "bytes", "sha256",
                                                       "digest_matches", "retrieved_at", "members")},
        "cases": {k: store[k] for k in ("cases", "skipped", "diseases", "split_cases", "split_diseases", "leaked_cases",
                                        "truth_not_in_hpoa", "cases_with_unresolved_terms", "seed")},
        "hpo_version": store["hpo_version"],
        "mondo_sssom": {k: v for k, v in mapping.manifest().items() if k != "path"},
        "local": {"before": metrics.summarise(args.before, args.reps), "after": metrics.summarise(args.after, args.reps)},
    }
    tuning = RESULTS / "tuning.json"
    if tuning.exists():
        out["tuning"] = json.loads(tuning.read_text("utf-8"))
    for key, path in (("remote_before", args.remote_before), ("remote_after", args.remote_after)):
        if path:
            data = json.loads(Path(path).read_text("utf-8"))
            out[key] = {k: v for k, v in data.items() if k != "rows"}
    sample = RESULTS / "sample.json"
    if sample.exists():
        out["remote_sample"] = json.loads(sample.read_text("utf-8"))
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "benchmark.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n", "utf-8")
    keep = ("case", "disease", "split", "bucket", "leaked", "in_hpoa", "n_present", "n_excluded", "rank_omim",
            "rank_all", "gene_rank", "truth_tied")
    before = {r["case"]: r for r in read_jsonl(work_dir() / f"local_{args.before}.jsonl.gz")}
    rows = []
    for r in read_jsonl(work_dir() / f"local_{args.after}.jsonl.gz"):
        row = {k: r.get(k) for k in keep}
        b = before.get(r["case"]) or {}
        row.update({"before_rank_omim": b.get("rank_omim"), "before_rank_all": b.get("rank_all"),
                    "before_gene_rank": b.get("gene_rank")})
        rows.append(row)
    write_jsonl(RESULTS / "local_cases.jsonl.gz", rows)
    print(json.dumps({"written": str(RESULTS / "benchmark.json"), "cases": len(rows)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
