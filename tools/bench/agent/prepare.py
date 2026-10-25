#!/usr/bin/env python3
"""Build the case list for the agent benchmark (zebra-mod + Claude, end to end).

    python3 -I tools/bench/agent/prepare.py --rarebench DIR [--mygene2 FILE] --out cases.jsonl

One row per case: id, dataset, present / excluded HPO ids with labels (from the pinned HPO
release), truth ids and truth names. Nothing else about the case (no PMID, no gene, no text)
goes to the model.

  RareBench (Chen et al., KDD 2024; HuggingFace chenxz/RareBench): MME, HMS, LIRICAL,
    RAMEDIS jsonl, truth = RareDisease ids (OMIM / ORPHA / CCRD), names from
    mapping/disease_mapping.json.
  phenopacket-store: the 240-case stratified held-out sample of tools/bench/results/sample.json
    (one case per disease, release 0.1.27), present and excluded terms as curated.
  RDguru: the 238-case test set released with RDguru (present terms only, as the paper ran it).
  MyGene2: optional, a jsonl with the same fields as RareBench.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from tools.bench.common import hpo_files, work_dir  # noqa: E402


def hpo_labels() -> dict:
    hp = json.loads((Path(hpo_files()["dir"]) / "hp.json").read_text("utf-8"))
    out = {}
    for n in hp["graphs"][0]["nodes"]:
        cid = n.get("id", "").rsplit("/", 1)[-1].replace("_", ":")
        if cid.startswith("HP:") and n.get("lbl"):
            out[cid] = n["lbl"]
    return out


def omim_names() -> dict:
    names = {}
    for line in open(Path(hpo_files()["dir"]) / "phenotype.hpoa", encoding="utf-8"):
        if line.startswith("#") or line.startswith("database_id"):
            continue
        f = line.split("\t")
        names.setdefault(f[0], f[1])
    return names


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rarebench", type=Path, required=True)
    ap.add_argument("--mygene2", type=Path)
    ap.add_argument("--rdguru", type=Path, help="checkout of github.com/1678492347/rdguru-test")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    labels, dnames = hpo_labels(), omim_names()
    dmap = json.loads((a.rarebench / "disease_mapping.json").read_text("utf-8"))
    rows = []

    def term(h):
        return {"id": h, "label": labels.get(h)}

    def add(ds, i, present, excluded, truth, tnames):
        rows.append({"id": f"{ds}-{i:04d}", "dataset": ds, "present": [term(h) for h in present],
                     "excluded": [term(h) for h in excluded], "truth": truth, "truth_names": tnames})

    sets = [(ds, a.rarebench / "x" / "data" / f"{ds}.jsonl") for ds in ("MME", "HMS", "LIRICAL", "RAMEDIS")]
    if a.mygene2:
        sets.append(("MyGene2", a.mygene2))
    for ds, path in sets:
        for i, line in enumerate(open(path, encoding="utf-8")):
            r = json.loads(line)
            truth = r["RareDisease"]
            add(ds, i, r["Phenotype"], [], truth, sorted({dmap.get(t) or dnames.get(t) or t for t in truth}))

    if a.rdguru:  # RDguru (Yang et al., JBHI 2024): 238 case reports, truth ORPHA, present terms only
        import csv
        for i, r in enumerate(csv.DictReader(open(a.rdguru / "summary.csv", encoding="utf-8", errors="replace"))):
            pp = json.loads((a.rdguru / "case_info" / f"{r['File Name']}.json").read_text("utf-8"))
            present = [f["type"]["id"] for f in pp.get("phenotypicFeatures", []) if not f.get("excluded")]
            add("RDguru", i, present, [], [f"ORPHA:{r['ORPHA']}"], [r["Disease Name"]])
            rows[-1]["source_case"] = r["File Name"]

    sample = json.loads((ROOT / "tools/bench/results/sample.json").read_text("utf-8"))
    keep = set()
    for v in sample["cases"].values():
        keep.update(v if isinstance(v, list) else [v])
    by_case = {}
    for line in open(work_dir() / "cases.jsonl", encoding="utf-8"):
        r = json.loads(line)
        if r["case"] in keep:
            by_case[r["case"]] = r
    for i, c in enumerate(sorted(by_case)):
        r = by_case[c]
        truth = r["truth"] if isinstance(r["truth"], list) else [r["truth"]]
        add("phenopacket", i, r["present"], r["excluded"], truth, sorted({dnames.get(t, t) for t in truth}))
        rows[-1]["source_case"] = c
        rows[-1]["bucket"] = r.get("bucket")
        rows[-1]["leaked"] = r.get("leaked")

    missing = sum(1 for r in rows for t in r["present"] + r["excluded"] if not t["label"])
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(json.dumps(r, ensure_ascii=False) + "\n")
    from collections import Counter
    print(json.dumps({"cases": len(rows), "by_dataset": Counter(r["dataset"] for r in rows),
                      "terms_without_label": missing}, default=dict))


if __name__ == "__main__":
    main()
