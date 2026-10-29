#!/usr/bin/env python3
"""Do RareBench literature cases reappear in phenopacket-store (which records each case's PMID and whether that
paper is an HPO annotation source)? Match on HPO present-term sets: exact, or Jaccard >= 0.8 with the same truth."""
import json, sys
from collections import defaultdict
from pathlib import Path
import argparse
ap = argparse.ArgumentParser()
ap.add_argument("--cases", type=Path, required=True)
ap.add_argument("--out", type=Path, required=True)
a = ap.parse_args()
cases = [json.loads(l) for l in open(a.cases)]
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from tools.bench.common import work_dir  # noqa
store = [json.loads(l) for l in open(work_dir() / "cases.jsonl")]
idx = defaultdict(list)
for s in store:
    S = frozenset(s["present"]); idx[len(S)].append((S, s))
out = {}
for ds in ["MME", "HMS", "LIRICAL", "RAMEDIS", "MyGene2", "RDguru"]:
    cs = [c for c in cases if c["dataset"] == ds]
    exact = near = leaked = 0
    for c in cs:
        A = frozenset(t["id"] for t in c["present"]); truth = set(c["truth"])
        best = None
        for L in range(max(1, int(len(A) * 0.8)), int(len(A) / 0.8) + 2):
            for S, s in idx.get(L, []):
                j = len(A & S) / len(A | S)
                if j >= 0.8 and (best is None or j > best[0]):
                    st = set(s["truth"] if isinstance(s["truth"], list) else [s["truth"]])
                    if j == 1.0 or st & truth:
                        best = (j, s)
        if best:
            near += 1; exact += best[0] == 1.0; leaked += bool(best[1].get("leaked"))
    out[ds] = {"n": len(cs), "matched": near, "exact": exact, "matched_paper_is_hpoa_source": leaked}
    print(ds, out[ds], flush=True)
json.dump(out, open(a.out, "w"), indent=1)
