#!/usr/bin/env python3
"""Why zebra-mod gains or loses against DeepRare per dataset: how retrievable each set is from itself.

For every case, the most similar OTHER case in the same dataset (Jaccard over HPO sets); a hit when that neighbour's
diagnosis matches (MONDO exact or descendant either way). Also disease concentration. Local only, no model calls.
"""
import json, sys, sys
from collections import Counter, defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from tools.bench.agent.score import mondo_of, ancestors_of, norm_id  # noqa

import argparse
ap = argparse.ArgumentParser()
ap.add_argument("--cases", type=Path, required=True)
ap.add_argument("--out", type=Path, required=True)
a = ap.parse_args()
cases = [json.loads(l) for l in open(a.cases)]
m, anc = mondo_of(), ancestors_of()

def mset(c):
    s = set()
    for t in c["truth"]:
        t = norm_id(t)
        if t:
            s |= m(t) | ({t} if t.startswith(("OMIM", "ORPHA")) else set())
    return s

def same(a, b):
    """The neighbour b's truth ids, read as predictions for case a, scored exactly as score.py scores the agent
    (lenient: a prediction is the truth or a MONDO descendant of it; one-way)."""
    truth = {norm_id(t) for t in a["truth"]} - {None}
    tm = set().union(*(m(t) for t in truth)) if truth else set()
    for pid in ({norm_id(t) for t in b["truth"]} - {None}):
        if pid in truth or any(anc(x) & tm for x in m(pid)):
            return True
    return False

def same_exact(a, b):
    return bool(({norm_id(t) for t in a["truth"]} - {None}) & ({norm_id(t) for t in b["truth"]} - {None}))

out = {}
for ds in ["MME", "HMS", "LIRICAL", "RAMEDIS", "MyGene2", "RDguru", "phenopacket"]:
    cs = [c for c in cases if c["dataset"] == ds]
    sets = [frozenset(t["id"] for t in c["present"]) for c in cs]
    dis = Counter(frozenset(c["truth"]) for c in cs)
    hit1 = hit5 = dup = dupsame = ex1 = 0
    for i, c in enumerate(cs):
        sims = sorted(((len(sets[i] & sets[j]) / max(1, len(sets[i] | sets[j])), j) for j in range(len(cs)) if j != i), reverse=True)
        if not sims:
            continue
        hit1 += same(c, cs[sims[0][1]])
        ex1 += same_exact(c, cs[sims[0][1]])
        hit5 += any(same(c, cs[j]) for _, j in sims[:5])
        if sims[0][0] >= 0.8:
            dup += 1
            dupsame += same(c, cs[sims[0][1]])
    n = len(cs)
    out[ds] = {"n": n, "diseases": len(dis), "cases_per_disease": round(n / len(dis), 2),
               "share_in_diseases_with_5plus": round(100 * sum(v for v in dis.values() if v >= 5) / n, 1),
               "top_disease": [sorted(dis.most_common(1)[0][0]), dis.most_common(1)[0][1]],
               "loo_1nn_R1": round(100 * hit1 / n, 1), "loo_1nn_R1_exact_id": round(100 * ex1 / n, 1), "loo_5nn_R5": round(100 * hit5 / n, 1),
               "near_dup_jacc80": round(100 * dup / n, 1), "near_dup_same_dx": round(100 * dupsame / n, 1),
               "median_terms": sorted(len(s) for s in sets)[n // 2]}
    print(ds, out[ds], flush=True)
json.dump(out, open(a.out, "w"), indent=1)
