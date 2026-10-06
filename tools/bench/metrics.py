#!/usr/bin/env python3
"""Top-k metrics for the local ranker runs, by split, stratum and leakage, with 95% intervals.

    python3 -I tools/bench/metrics.py LABEL [LABEL ...] [--json OUT]

Reads work_dir()/local_<label>.jsonl.gz (tools/bench/run_local.py) and prints, for each
label: top-1/3/10/25 of the OMIM-only disease rank (primary), of the OMIM+Orphanet rank
(exact MONDO joins), and of the gene rank, for all cases and for each subset below.
Subsets: dev / heldout; leaked / not leaked (the case's own paper is a source of its
disease's HPOA annotations); present-term strata 1-3, 4-6, 7-10, >10.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.bench.common import BUCKETS, read_jsonl, topk_with_ci, work_dir  # noqa: E402

INHERITANCE = {"HP:0000006": "autosomal dominant", "HP:0000007": "autosomal recessive",
               "HP:0001417": "X-linked", "HP:0001419": "X-linked", "HP:0001423": "X-linked",
               "HP:0001427": "mitochondrial"}


def disease_categories() -> dict:
    """disease -> (inheritance, dominant organ system). Cheap proxies for "disease category":
    the inheritance modes HPOA annotates (aspect I; several -> "mixed"), and the top-level HPO
    branch (a child of Phenotypic abnormality) holding most of the disease's annotations."""
    import csv
    from collections import Counter, defaultdict

    from zebra import hpo_local

    idx = hpo_local.load()
    modes = defaultdict(set)
    with open(hpo_local.data_dir() / "phenotype.hpoa", encoding="utf-8") as fh:
        for row in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            if row.get("aspect") == "I" and row["hpo_id"] in INHERITANCE:
                modes[row["database_id"]].add(INHERITANCE[row["hpo_id"]])
    top = {t for t, ps in idx.parents.items() if hpo_local.ROOT_PHENO in ps}
    out = {}
    for d, terms in idx.disease_terms.items():
        c = Counter()
        for t in terms:
            for a in idx.ancestors(idx.alt.get(t, t)) & top:
                c[a] += 1
        # ties broken by term id, so the category does not depend on set iteration order
        organ = idx.names.get(min(c.items(), key=lambda kv: (-kv[1], kv[0]))[0], "?") if c else "?"
        m = modes.get(d) or set()
        out[d] = ("mixed" if len(m) > 1 else next(iter(m)) if m else "not annotated", organ)
    return out

SUBSETS = [("all", lambda r: True), ("dev", lambda r: r["split"] == "dev"),
           ("heldout", lambda r: r["split"] == "heldout"),
           ("heldout, not leaked", lambda r: r["split"] == "heldout" and not r["leaked"]),
           ("leaked", lambda r: r["leaked"]), ("not leaked", lambda r: not r["leaked"])] + \
          [(f"heldout, {b} present terms", (lambda b: lambda r: r["split"] == "heldout" and r["bucket"] == b)(b))
           for b, _lo, _hi in BUCKETS] + \
          [(f"all, {b} present terms", (lambda b: lambda r: r["bucket"] == b)(b)) for b, _lo, _hi in BUCKETS]


_CATS = None


def summarise(label: str, reps: int = 1000) -> dict:
    global _CATS
    if _CATS is None:
        _CATS = disease_categories()
    rows = [r for r in read_jsonl(work_dir() / f"local_{label}.jsonl.gz")]
    for r in rows:
        r["inheritance"], r["organ"] = _CATS.get(r["disease"], ("not annotated", "?"))
    errors = [r for r in rows if r.get("error")]
    ok = [r for r in rows if not r.get("error")]
    out = {"label": label, "cases": len(rows), "errors": len(errors), "subsets": {}}
    try:
        out["run"] = json.loads((work_dir() / f"local_{label}.meta.json").read_text("utf-8"))
    except OSError:
        pass
    extra = []
    for inh in sorted({r["inheritance"] for r in ok}):
        extra.append((f"heldout, inheritance: {inh}", (lambda v: lambda r: r["split"] == "heldout" and r["inheritance"] == v)(inh)))
    organs = {}
    for r in ok:
        organs[r["organ"]] = organs.get(r["organ"], 0) + 1
    for org, n in sorted(organs.items(), key=lambda kv: -kv[1])[:8]:
        extra.append((f"heldout, organ system: {org}", (lambda v: lambda r: r["split"] == "heldout" and r["organ"] == v)(org)))
    for name, pred in SUBSETS + extra:
        sub = [r for r in ok if pred(r)]
        genes = [r for r in sub if r.get("has_gene")]
        out["subsets"][name] = {
            "disease_omim": topk_with_ci(sub, "rank_omim", reps=reps),
            "disease_omim_orpha": topk_with_ci(sub, "rank_all", reps=reps),
            "gene": topk_with_ci(genes, "gene_rank", reps=reps),
        }
    return out


def fmt(cell: dict, k: int, ci: bool = True) -> str:
    t = cell.get(f"top{k}")
    if not t:
        return "-"
    lo, hi = t["ci95_cluster"]
    return f"{t['rate'] * 100:.1f}" + (f" ({lo * 100:.1f}-{hi * 100:.1f})" if ci else "")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("labels", nargs="+")
    ap.add_argument("--json")
    ap.add_argument("--reps", type=int, default=1000)
    args = ap.parse_args(argv)
    res = {lab: summarise(lab, args.reps) for lab in args.labels}
    for lab, r in res.items():
        print(f"## {lab}: {r['cases']} cases, {r['errors']} errors")
        print("| subset | n | OMIM top-1 | top-3 | top-10 | top-25 | +ORPHA top-10 | gene n | gene top-1 | gene top-10 |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for name, s in r["subsets"].items():
            d, a, g = s["disease_omim"], s["disease_omim_orpha"], s["gene"]
            print(f"| {name} | {d['n']} | {fmt(d, 1)} | {fmt(d, 3)} | {fmt(d, 10)} | {fmt(d, 25)} | {fmt(a, 10)} | "
                  f"{g['n']} | {fmt(g, 1)} | {fmt(g, 10)} |")
    if args.json:
        Path(args.json).write_text(json.dumps(res, indent=1) + "\n", "utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
