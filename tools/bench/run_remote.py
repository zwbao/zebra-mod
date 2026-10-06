#!/usr/bin/env python3
"""Run the web rankers (Monarch semsim, PubCaseFinder) on a stratified sample of benchmark cases.

    python3 -I tools/bench/run_remote.py sample [--per-bucket 60]
    python3 -I tools/bench/run_remote.py monarch [--per-bucket 60]
    python3 -I tools/bench/run_remote.py pubcasefinder [--per-bucket 25]

`sample` draws the sample once and writes tools/bench/results/sample.json: held-out
cases only (the local ranker is tuned on the development half, so every consensus
number stays held-out), one case per disease, the same number of cases in each
present-term stratum (1-3, 4-6, 7-10, >10). The rankers then take the first N cases
of each stratum, so the PubCaseFinder cases are a subset of the Monarch cases and the
two are compared on the same patients.

Rate limits. PubCaseFinder publishes "no more than 10 requests per minute, 100 per
hour, 1,000 per day" (swagger-pubcasefinder.json, checked 2026-10-06). zebra.http
already spaces requests to the host 6 s apart; this runner adds rolling windows
(default at most 90 requests an hour and 900 a day, leaving room for anything else on
this machine) and a minimum gap of 40 s, in a state file shared under a lock, and it
counts each request before sending it, with retries off (Window).
Monarch publishes no limit; requests are sequential with a 1 s gap. Two requests per
case: diseases + genes (PubCaseFinder: the OMIM list and the gene list; the truth in
phenopacket-store is always an OMIM id, so the Orphanet list is not requested).
Every answer goes through zebra's cache ($ZEBRA_CACHE_DIR, 7 days), so a rerun within
the week costs no requests.

Output: work_dir()/remote_<source>.jsonl, appended case by case (a rerun resumes).
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.bench.common import BUCKETS, RESULTS, SEED, read_jsonl, work_dir  # noqa: E402
from zebra import hpo_local  # noqa: E402
from zebra.http import SourceError  # noqa: E402

LIMIT = 50


def draw_sample(per_bucket: int) -> dict:
    cases = [c for c in read_jsonl(work_dir() / "cases.jsonl") if c["split"] == "heldout"]
    by_disease = defaultdict(list)
    for c in cases:
        by_disease[c["disease"]].append(c)
    rng = random.Random(SEED + 1)
    diseases = sorted(by_disease)
    rng.shuffle(diseases)
    chosen = {name: [] for name, _lo, _hi in BUCKETS}
    for d in diseases:
        pool = sorted(by_disease[d], key=lambda c: c["case"])
        rng.shuffle(pool)
        for c in pool:
            if len(chosen[c["bucket"]]) < per_bucket:
                chosen[c["bucket"]].append(c["case"])
                break
        if all(len(v) >= per_bucket for v in chosen.values()):
            break
    out = {"seed": SEED + 1, "per_bucket": per_bucket, "split": "heldout", "rule": "one case per disease",
           "cases": chosen, "counts": {k: len(v) for k, v in chosen.items()}}
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "sample.json").write_text(json.dumps(out, indent=1) + "\n", "utf-8")
    return out


def sample_cases(per_bucket: int):
    """The first `per_bucket` cases of each stratum, interleaved (one per stratum in turn), so a run
    that is stopped early still covers every stratum evenly."""
    sample = json.loads((RESULTS / "sample.json").read_text("utf-8"))
    wanted = []
    for i in range(per_bucket):
        for name, _lo, _hi in BUCKETS:
            if i < len(sample["cases"][name]):
                wanted.append(sample["cases"][name][i])
    rows = {c["case"]: c for c in read_jsonl(work_dir() / "cases.jsonl")}
    return [rows[c] for c in wanted]


class Window:
    """Rolling limits on network requests to one host, shared through a locked state file.

    `reserve()` blocks until one more request fits -- at least `gap` seconds after the last,
    fewer than `per_hour` in the last hour and `per_day` in the last 24 h -- and records it
    BEFORE the request is sent, so a request that fails, times out or is retried inside a
    call still counts, and two runners on one machine share one budget (fcntl lock).
    """

    def __init__(self, path: Path, per_hour: int, gap: float, per_day: int = 900):
        self.path, self.per_hour, self.gap, self.per_day = path, per_hour, gap, per_day
        self.lock_path = path.with_suffix(path.suffix + ".lock")

    def _load(self, now: float) -> list:
        try:
            stamps = json.loads(self.path.read_text("utf-8"))
        except (OSError, ValueError):
            stamps = []
        return sorted(t for t in stamps if isinstance(t, (int, float)) and now - t < 86400)

    @property
    def stamps(self) -> list:
        return [t for t in self._load(time.time()) if time.time() - t < 3600]

    def reserve(self) -> None:
        import fcntl

        while True:
            with open(self.lock_path, "a+") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                now = time.time()
                stamps = self._load(now)
                hour = [t for t in stamps if now - t < 3600]
                waits = []
                if stamps:
                    waits.append(stamps[-1] + self.gap - now)
                if len(hour) >= self.per_hour:
                    waits.append(hour[0] + 3600 - now + 1)
                if len(stamps) >= self.per_day:
                    waits.append(stamps[0] + 86400 - now + 1)
                w = max(waits) if waits else 0
                if w <= 0:
                    stamps.append(now)
                    self.path.write_text(json.dumps(stamps), "utf-8")
                    return
            time.sleep(min(w, 60))


def remote_terms(idx, present):
    out = []
    for t in present:
        pid, _note = idx.primary(t)
        if pid and pid not in out:
            out.append(pid)
    return sorted(out)


def run_monarch(case, terms):
    from zebra.sources import monarch

    rec = {"case": case["case"], "source": "monarch", "terms": terms}
    t0 = time.time()
    for group, key in (("Human Diseases", "diseases"), ("Human Genes", "genes")):
        try:
            got = monarch.semsim_rank(terms, limit=LIMIT, group=group)
            rec[key] = [{"id": h["id"], "name": h["name"], "symbol": h.get("symbol"), "xref": h["xref"],
                         "score": h.get("score")} for h in got.result["hits"]]
            rec[f"{key}_cached"] = bool(got.sources and got.sources[0].get("cached"))
        except (SourceError, ValueError, KeyError, TypeError, AttributeError) as err:
            rec[f"{key}_error"] = f"{type(err).__name__}: {err}"[:300]
        time.sleep(1.0)
    rec["seconds"] = round(time.time() - t0, 2)
    return rec


def run_pcf(case, terms, window: Window):
    from zebra.sources import pubcasefinder as pcf

    rec = {"case": case["case"], "source": "pubcasefinder", "terms": terms}
    t0 = time.time()
    for target, key in (("omim", "diseases"), ("gene", "genes")):
        window.reserve()  # counted before it is sent: a failure or a cache re-request still counts
        try:
            # retries=0: one call is one request (two only when a bad body came from the cache,
            # which the reserve before it already allows for by the 40 s gap and the margin to 100/h)
            got = pcf.ranked(terms, target, LIMIT, timeout=90, retries=0)
            cached = bool(got.sources and got.sources[0].get("cached"))
            rec[key] = [({"id": h["id"], "name": h["name"], "served_rank": h["rank"]} if target != "gene"
                         else {"symbol": h["symbol"], "gene_id": h["gene_id"], "served_rank": h["rank"]})
                        for h in got.result["hits"]]
            rec[f"{key}_cached"] = cached
        except (SourceError, ValueError, KeyError, TypeError, AttributeError) as err:
            rec[f"{key}_error"] = f"{type(err).__name__}: {err}"[:300]
    rec["seconds"] = round(time.time() - t0, 2)
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=("sample", "monarch", "pubcasefinder"))
    ap.add_argument("--per-bucket", type=int)
    ap.add_argument("--per-hour", type=int, default=90, help="PubCaseFinder network requests per rolling hour (max 100)")
    ap.add_argument("--gap", type=float, default=40.0, help="seconds between PubCaseFinder network requests (min 6)")
    ap.add_argument("--per-day", type=int, default=900, help="PubCaseFinder requests per rolling 24 h (max 1000)")
    args = ap.parse_args(argv)
    if args.what == "sample":
        print(json.dumps(draw_sample(args.per_bucket or 60)["counts"]))
        return 0
    if args.per_hour > 100 or args.gap < 6 or args.per_day > 1000:
        raise SystemExit("PubCaseFinder's published limits are 10/minute, 100/hour and 1,000/day: refusing "
                         "--per-hour > 100, --gap < 6 or --per-day > 1000")
    per_bucket = args.per_bucket or (60 if args.what == "monarch" else 25)
    idx = hpo_local.load()
    out_path = work_dir() / f"remote_{args.what}.jsonl"
    # a case whose lists came back with an error is run again (its old row is superseded: remote_eval
    # keeps the last row per case)
    done = {r["case"] for r in read_jsonl(out_path) if not any(k.endswith("_error") for k in r)} \
        if out_path.exists() else set()
    window = Window(work_dir() / "pcf_requests.json", args.per_hour, args.gap, args.per_day)
    todo = [c for c in sample_cases(per_bucket) if c["case"] not in done]
    print(f"{args.what}: {len(todo)} cases to run, {len(done)} already done", flush=True)
    for i, case in enumerate(todo, 1):
        terms = remote_terms(idx, case["present"])
        rec = run_monarch(case, terms) if args.what == "monarch" else run_pcf(case, terms, window)
        with open(out_path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
        errs = [k for k in rec if k.endswith("_error")]
        print(f"[{i}/{len(todo)}] {case['case']} {rec['seconds']}s {'ERR ' + ','.join(errs) if errs else 'ok'}",
              flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
