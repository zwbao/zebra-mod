#!/usr/bin/env python3
"""Build the benchmark case list from a phenopacket-store release.

    python3 -I tools/bench/cases.py [--store DIR/<tag>]

Reads every phenopacket in all_phenopackets.zip (without extracting it) through
zebra.phenopacket.read, and writes work_dir()/cases.jsonl, one row per case:

  case          zip member name (gene folder / PMID_..._individual.json)
  pmid          the case report's PMID (metaData.externalReferences)
  truth         the diagnosed disease ids (interpretations + non-excluded diseases)
  genes         causal gene symbols (interpretationStatus CAUSATIVE or CONTRIBUTORY)
  present / excluded   HPO ids as recorded
  n_present, bucket    number of present terms and its stratum (1-3, 4-6, 7-10, >10)
  split         "dev" or "heldout": diseases are shuffled with a fixed seed and each
                disease goes, with all its cases, to the half that has fewer cases so far
  in_hpoa       the truth disease has phenotype annotations in the installed HPO release
                (without them the local ranker cannot find it at all)
  leaked        the truth disease's HPOA rows cite this case's PMID: the case's own paper
                was one of the sources of the annotations it is scored against
  unresolved    present ids the installed HPO release does not resolve (newer/obsolete)

Also writes work_dir()/store.json (the release manifest plus counts).
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.bench.common import SEED, bucket, hpo_files, version_key, work_dir, write_jsonl  # noqa: E402
from tools.bench.fetch_store import default_dir  # noqa: E402
from zebra import hpo_local, phenopacket  # noqa: E402


def hpoa_references() -> dict:
    refs = defaultdict(set)
    path = hpo_local.data_dir() / "phenotype.hpoa"
    with open(path, encoding="utf-8") as fh:
        rows = (line for line in fh if not line.startswith("#"))
        for row in csv.DictReader(rows, delimiter="\t"):
            for ref in (row.get("reference") or "").split(";"):
                ref = ref.strip()
                if ref.startswith("PMID:"):
                    refs[row["database_id"]].add(ref)
    return refs


def build(store: Path) -> dict:
    manifest = json.loads((store / "manifest.json").read_text("utf-8"))
    idx = hpo_local.load()
    refs = hpoa_references()
    rows = []
    skipped = defaultdict(int)
    with zipfile.ZipFile(store / "all_phenopackets.zip") as zf:
        for name in sorted(n for n in zf.namelist() if n.endswith(".json")):
            try:
                pp = phenopacket.read(json.loads(zf.read(name)))
            except (ValueError, phenopacket.PhenopacketError) as err:
                skipped[f"unreadable: {type(err).__name__}"] += 1
                continue
            truth = [d["id"] for d in pp["diseases"] if not d["excluded"]]
            if not truth:
                skipped["no diagnosis"] += 1
                continue
            present = [t["id"] for t in pp["present"]]
            if not present:
                skipped["no present HPO term"] += 1
                continue
            genes = sorted({g["symbol"] for g in pp["genes"] if g["symbol"]
                            and (g["status"] in (None, "CAUSATIVE", "CONTRIBUTORY"))})
            pmid = next((r for r in pp["references"] if r.startswith("PMID:")), None)
            primary = truth[0]
            rows.append({
                "case": name.split("/", 1)[-1], "pmid": pmid, "truth": truth, "disease": primary, "genes": genes,
                "present": present, "excluded": [t["id"] for t in pp["excluded"]],
                "n_present": len(present), "bucket": bucket(len(present)),
                "in_hpoa": primary in idx.disease_terms,
                "leaked": bool(pmid and pmid in refs.get(primary, ())),
                "unresolved": [t for t in present if idx.primary(t)[0] is None],
            })
    # split by disease: shuffle diseases, then fill the smaller half
    per = defaultdict(int)
    for r in rows:
        per[r["disease"]] += 1
    diseases = sorted(per)
    random.Random(SEED).shuffle(diseases)
    halves = {"dev": 0, "heldout": 0}
    assign = {}
    for d in diseases:
        half = "dev" if halves["dev"] <= halves["heldout"] else "heldout"
        assign[d] = half
        halves[half] += per[d]
    for r in rows:
        r["split"] = assign[r["disease"]]
    out = work_dir() / "cases.jsonl"
    write_jsonl(out, rows)
    info = {"store": manifest, "hpo_version": idx.version, "hpo_files": hpo_files(), "cases": len(rows),
            "skipped": dict(skipped),
            "diseases": len(per), "split_cases": halves,
            "split_diseases": {h: sum(1 for d in assign.values() if d == h) for h in halves},
            "leaked_cases": sum(r["leaked"] for r in rows), "truth_not_in_hpoa": sum(not r["in_hpoa"] for r in rows),
            "cases_with_unresolved_terms": sum(bool(r["unresolved"]) for r in rows), "seed": SEED,
            "path": str(out)}
    (work_dir() / "store.json").write_text(json.dumps(info, indent=2) + "\n", "utf-8")
    return info


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", type=Path, help="directory holding all_phenopackets.zip and manifest.json")
    args = ap.parse_args(argv)
    store = args.store
    if store is None:
        tags = sorted((p for p in default_dir().iterdir() if (p / "manifest.json").exists()),
                      key=lambda p: version_key(p.name))
        if not tags:
            raise SystemExit("no phenopacket-store release found: run tools/bench/fetch_store.py first")
        store = tags[-1]
    print(json.dumps(build(store), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
