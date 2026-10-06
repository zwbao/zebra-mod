"""Regression gate for the phenotype ranker against the recorded benchmark (docs/BENCHMARK.md).

Marked live because it needs the installed HPO release and the phenopacket-store case list
(tools/bench/fetch_store.py + cases.py); it skips when either is missing or when the installed
HPO release differs from the one the benchmark recorded (the positions would legitimately move).
It re-scores a fixed slice of held-out cases with the configuration `zebra phenotype rank` ships
and requires the very positions recorded in tools/bench/results/local_cases.jsonl.gz: any change
to the ranking is caught here, and must come with a re-run of the benchmark.
"""

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "tools" / "bench" / "results"


@pytest.fixture
def bench(monkeypatch):
    monkeypatch.delenv("ZEBRA_CACHE_DIR", raising=False)  # the real cache holds the case list
    monkeypatch.delenv("ZEBRA_HPO_DIR", raising=False)
    from tools.bench.common import read_jsonl, work_dir
    from zebra import hpo_local

    summary = RESULTS / "benchmark.json"
    recorded = RESULTS / "local_cases.jsonl.gz"
    cases_path = work_dir() / "cases.jsonl"
    if not (summary.exists() and recorded.exists() and cases_path.exists()):
        pytest.skip("benchmark results or the phenopacket-store case list are not present")
    try:
        idx = hpo_local.load()
    except hpo_local.HpoDataMissing as err:
        pytest.skip(str(err))
    meta = json.loads(summary.read_text("utf-8"))
    if idx.version != meta["hpo_version"]:
        pytest.skip(f"installed HPO {idx.version} differs from the benchmark's {meta['hpo_version']}")
    cases = {c["case"]: c for c in read_jsonl(cases_path)}
    rows = [r for r in read_jsonl(recorded) if r["split"] == "heldout" and r["case"] in cases]
    rows.sort(key=lambda r: r["case"])
    return idx, cases, rows[::25][:200], meta


@pytest.mark.live
def test_shipped_configuration_reproduces_the_recorded_positions(bench):
    from zebra import hpo_local
    from zebra.commands.phenotype import LOCAL_METHOD, LOCAL_PARAMS

    idx, cases, rows, _meta = bench
    assert rows, "no held-out rows to check"
    moved = []
    for r in rows:
        c = cases[r["case"]]
        res = hpo_local.rank(idx, c["present"], c["excluded"], top=50, db=("OMIM",), method=LOCAL_METHOD,
                             params=dict(LOCAL_PARAMS))
        ids = [d["disease"] for d in res["diseases"]]
        pos = next((i for i, d in enumerate(ids, 1) if d in c["truth"]), None)
        if pos != r["rank_omim"]:
            moved.append((r["case"], r["rank_omim"], pos))
    assert not moved, (f"{len(moved)} of {len(rows)} positions moved (case, recorded, now): {moved[:5]}. The shipped "
                       f"configuration is {LOCAL_METHOD} {LOCAL_PARAMS}; the benchmark recorded "
                       f"{_meta['local']['after'].get('run', {}).get('options')}. Re-run tools/bench (docs/BENCHMARK.md).")


@pytest.mark.live
def test_the_recorded_benchmark_is_of_the_shipped_configuration(bench):
    from zebra.commands.phenotype import LOCAL_METHOD, LOCAL_PARAMS

    _idx, _cases, _rows, meta = bench
    recorded = meta["local"]["after"]["run"]["options"]
    assert recorded.get("method", "resnik") == LOCAL_METHOD
    assert recorded.get("params", {}) == LOCAL_PARAMS, (
        f"shipped {LOCAL_PARAMS} differs from the benchmarked {recorded.get('params')}: re-run tools/bench")
    assert meta["local"]["after"]["subsets"]["heldout"]["disease_omim"]["top10"]["rate"] >= 0.65
