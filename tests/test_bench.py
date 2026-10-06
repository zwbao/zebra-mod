"""tools/bench: the phenopacket-store benchmark harness (offline parts)."""

import json
import os
import shutil
import zipfile

import pytest

from tools.bench import common, mapping
from tools.bench import cases as bench_cases
from tools.bench import run_remote
from zebra import hpo_local

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "hpo")


def test_wilson_interval_matches_the_textbook_value():
    lo, hi = common.wilson(5, 10)
    assert lo == pytest.approx(0.2366, abs=1e-4) and hi == pytest.approx(0.7634, abs=1e-4)
    assert common.wilson(0, 0) == (None, None)


def test_topk_counts_positions_and_resamples_whole_diseases():
    rows = ([{"disease": "D1", "pos": 1}] * 4 + [{"disease": "D2", "pos": 12}] * 4 + [{"disease": "D3", "pos": None}] * 2)
    out = common.topk_with_ci(rows, "pos", ks=(1, 10, 25), reps=300)
    assert out["n"] == 10 and out["diseases"] == 3
    assert out["top1"]["rate"] == 0.4 and out["top10"]["rate"] == 0.4 and out["top25"]["rate"] == 0.8
    lo, hi = out["top1"]["ci95_cluster"]
    assert lo <= 0.4 <= hi
    # three clusters: the bootstrap interval is far wider than the case-level Wilson interval would suggest
    assert hi - lo > out["top1"]["ci95_wilson"][1] - out["top1"]["ci95_wilson"][0] - 0.2
    assert common.topk_with_ci([], "pos") == {"n": 0, "diseases": 0}


def test_buckets():
    assert [common.bucket(n) for n in (1, 3, 4, 6, 7, 10, 11, 40)] == ["1-3", "1-3", "4-6", "4-6", "7-10", "7-10",
                                                                     ">10", ">10"]


@pytest.fixture
def fake_sssom(monkeypatch):
    rows = {
        "OMIM:607208": [{"subject": "MONDO:0011803", "subject_label": "DEE6B", "predicate": "skos:exactMatch",
                         "object": "OMIM:607208"}],
        "ORPHA:999": [{"subject": "MONDO:0011803", "subject_label": "DEE6B", "predicate": "skos:exactMatch",
                       "object": "ORPHA:999"}],
        "ORPHA:33069": [{"subject": "MONDO:0100135", "subject_label": "Dravet syndrome", "predicate": "skos:exactMatch",
                         "object": "ORPHA:33069"}],
    }
    mapping.rows_by_object.cache_clear()
    monkeypatch.setattr(mapping, "rows_by_object", lambda: rows)
    return rows


def test_a_hit_counts_only_through_an_exact_join(fake_sssom):
    truth = ["OMIM:607208"]
    # an Orphanet id exactly mapped to the same MONDO class is the same disease
    assert mapping.first_match(["OMIM:1", "ORPHA:999", "OMIM:607208"], truth) == 2
    # a different (broader) class is a miss, however related
    assert mapping.first_match(["ORPHA:33069", "OMIM:2"], truth) is None
    # a Monarch MONDO hit counts when it is the truth's exactMatch or carries the truth as an xref
    items = [{"source": "monarch", "id": "MONDO:0100135", "xref": ["ORPHA:33069"]},
             {"source": "monarch", "id": "MONDO:0011803", "xref": []}]
    assert mapping.first_match(None, truth, items) == 2
    items = [{"source": "monarch", "id": "MONDO:0000001", "xref": ["OMIM:607208"]}]
    assert mapping.first_match(None, truth, items) == 1


def _store(tmp_path, packets):
    store = tmp_path / "store" / "0.0.1"
    store.mkdir(parents=True)
    with zipfile.ZipFile(store / "all_phenopackets.zip", "w") as zf:
        for name, pp in packets.items():
            zf.writestr(f"0.0.1/{name}", json.dumps(pp))
    (store / "manifest.json").write_text(json.dumps({"tag": "0.0.1", "sha256": "x"}), "utf-8")
    return store


def _pp(pmid, disease, present, excluded=(), gene=None):
    gi = [{"interpretationStatus": "CAUSATIVE",
           "variantInterpretation": {"variationDescriptor": {"id": "v", "geneContext": {"symbol": gene, "valueId": "x"}}}}] \
        if gene else []
    return {"id": pmid, "phenotypicFeatures": [{"type": {"id": t}} for t in present]
            + [{"type": {"id": t}, "excluded": True} for t in excluded],
            "interpretations": [{"progressStatus": "SOLVED", "diagnosis": {"disease": {"id": disease},
                                                                           "genomicInterpretations": gi}}],
            "metaData": {"externalReferences": [{"id": pmid}]}}


def test_cases_split_by_disease_and_flag_leakage(tmp_path, monkeypatch):
    hpo = tmp_path / "hpo"
    hpo.mkdir()
    for f in hpo_local.FILES:
        shutil.copy(os.path.join(FIX, f), hpo / f)
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(hpo))
    monkeypatch.setenv("ZEBRA_BENCH_DIR", str(tmp_path / "work"))
    packets = {}
    for i in range(6):
        packets[f"G/PMID_1_{i}.json"] = _pp("PMID:1", "OMIM:1", ["HP:0002373", "HP:0001263"], ["HP:0001252"], "SCN1A")
        packets[f"H/PMID_2_{i}.json"] = _pp("PMID:2", "OMIM:2", ["HP:0001252"])
        packets[f"I/PMID_3_{i}.json"] = _pp("PMID:3", f"OMIM:{10 + i}", ["HP:0000518"])
    packets["J/none.json"] = {"id": "x", "phenotypicFeatures": [], "interpretations": []}
    info = bench_cases.build(_store(tmp_path, packets))
    rows = common.read_jsonl(tmp_path / "work" / "cases.jsonl")
    assert info["cases"] == 18 and info["skipped"] == {"no diagnosis": 1}
    halves = {}
    for r in rows:
        halves.setdefault(r["disease"], set()).add(r["split"])
    assert all(len(h) == 1 for h in halves.values()), "a disease must never be in both halves"
    assert {r["split"] for r in rows} == {"dev", "heldout"}
    by = {r["case"]: r for r in rows}
    # the fixture HPOA cites PMID:1 for OMIM:1, so those cases were part of their disease's own annotations
    assert by["G/PMID_1_0.json"]["leaked"] is True and by["H/PMID_2_0.json"]["leaked"] is False
    assert by["G/PMID_1_0.json"]["genes"] == ["SCN1A"] and by["G/PMID_1_0.json"]["excluded"] == ["HP:0001252"]
    assert by["I/PMID_3_0.json"]["in_hpoa"] is False  # OMIM:10 has no annotations: local can never find it
    again = bench_cases.build(_store(tmp_path / "again", packets))
    assert again["split_cases"] == info["split_cases"], "the split is deterministic"


def test_remote_sample_is_heldout_one_case_per_disease(tmp_path, monkeypatch):
    work = tmp_path / "work"
    monkeypatch.setenv("ZEBRA_BENCH_DIR", str(work))
    monkeypatch.setattr(run_remote, "RESULTS", tmp_path / "results")
    rows = []
    for d in range(30):
        for j, n in enumerate((2, 5, 8, 12)):
            rows.append({"case": f"c{d}_{j}", "disease": f"OMIM:{d}", "bucket": common.bucket(n),
                         "split": "heldout" if d % 3 else "dev", "present": ["HP:0001250"] * n})
    common.write_jsonl(work / "cases.jsonl", rows)
    out = run_remote.draw_sample(4)
    chosen = [c for v in out["cases"].values() for c in v]
    assert out["counts"] == {"1-3": 4, "4-6": 4, "7-10": 4, ">10": 4}
    diseases = [c.split("_")[0] for c in chosen]
    assert len(set(diseases)) == len(diseases), "one case per disease"
    assert all(int(d[1:]) % 3 for d in diseases), "held-out diseases only"
    picked = run_remote.sample_cases(2)
    assert len(picked) == 8
    # interleaved by stratum, so a run stopped early still covers every stratum
    assert [c["case"] for c in picked][:4] == [out["cases"][b][0] for b in ("1-3", "4-6", "7-10", ">10")]


def test_pubcasefinder_window_keeps_the_published_limits(tmp_path, monkeypatch):
    now = [100_000.0]
    slept = []
    monkeypatch.setattr(run_remote.time, "time", lambda: now[0])

    def fake_sleep(s):
        slept.append(s)
        now[0] += s

    monkeypatch.setattr(run_remote.time, "sleep", fake_sleep)
    w = run_remote.Window(tmp_path / "w.json", per_hour=3, gap=40, per_day=5)
    for _ in range(3):
        w.reserve()
    assert sum(slept) == pytest.approx(80, abs=1e-6)  # two 40 s gaps
    w.reserve()  # the fourth request in the hour waits for the first to age out
    assert now[0] >= 100_000.0 + 3600
    # a second runner on the same state file shares the budget (and the file survives restarts)
    w2 = run_remote.Window(tmp_path / "w.json", per_hour=3, gap=40, per_day=5)
    w2.reserve()
    assert len(json.loads((tmp_path / "w.json").read_text())) == 5
    start = now[0]
    w2.reserve()  # the sixth request in 24 h waits for the first to be a day old
    assert now[0] >= 100_000.0 + 86400 and now[0] > start


def test_pubcasefinder_requests_are_counted_before_they_are_sent(tmp_path, monkeypatch):
    """A failing call still uses budget: the reservation is made before the request."""
    from zebra.http import SourceError
    from zebra.sources import pubcasefinder

    calls = []

    def failing(terms, target, limit, timeout=60, retries=1):
        calls.append(retries)
        raise SourceError("PubCaseFinder", "u", 503, "busy")

    monkeypatch.setattr(pubcasefinder, "ranked", failing)
    w = run_remote.Window(tmp_path / "w.json", per_hour=90, gap=0.0)
    rec = run_remote.run_pcf({"case": "c"}, ["HP:0001250"], w)
    assert calls == [0, 0], "retries are off: one call is one request"
    assert len(json.loads((tmp_path / "w.json").read_text())) == 2
    assert "diseases_error" in rec and "genes_error" in rec


def test_remote_runner_refuses_rates_above_the_published_limit():
    with pytest.raises(SystemExit, match="published limits"):
        run_remote.main(["pubcasefinder", "--per-hour", "150"])
    with pytest.raises(SystemExit, match="published limits"):
        run_remote.main(["pubcasefinder", "--gap", "1"])
    with pytest.raises(SystemExit, match="published limits"):
        run_remote.main(["pubcasefinder", "--per-day", "2000"])
