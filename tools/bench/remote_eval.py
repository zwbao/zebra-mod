#!/usr/bin/env python3
"""Score Monarch, PubCaseFinder, the local ranker and their consensus on the web sample.

    python3 -I tools/bench/remote_eval.py [--method resnik|lr] [--module FILE] [--hpo-dir DIR] [--top 15]
                                          [--label NAME] [--json OUT]

Reads work_dir()/remote_monarch.jsonl and remote_pubcasefinder.jsonl (run_remote.py),
re-ranks the same cases with the local ranker (OMIM + Orphanet, as `zebra phenotype
rank` does), and rebuilds the consensus exactly as the CLI does: the top `--top`
(default 15, the CLI default) of each source, joined with
zebra.commands.phenotype.resolve_keys / build_consensus / gene_consensus. The joins
use MONDO's SSSOM exactMatch rows (tools/bench/mapping.py) plus the xrefs Monarch
returns with each hit -- the same kind of exact links the CLI fetches live.

Differences from the CLI, kept on purpose and stated in docs/BENCHMARK.md: PubCaseFinder's
Orphanet list is not requested (budget: two requests per case, not three), so the consensus
here has one PubCaseFinder list where the CLI has two; exact links come from MONDO's SSSOM file
(obsolete MONDO subjects are not followed to their replacement, as the CLI does through OLS).

Scoring a source's list (position as listed, and `ties_fair`: the truth equally likely at any
position of the block of hits sharing its score -- PubCaseFinder serves tied ranks, 1, 1, 3):
  local      truth OMIM id, or an Orphanet/MONDO id joined to it by exactMatch
  monarch    a MONDO hit joined to the truth by exactMatch or carrying it as an xref;
             MONDO hits with no OMIM link at all are counted and reported as `unlinkable`
  pubcasefinder  the OMIM list: plain id equality
  consensus  the first consensus row whose join key is the truth's key
A list that came back with an error is a miss for that case, and the error is counted.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.bench import mapping  # noqa: E402
from tools.bench.common import read_jsonl, topk_with_ci, work_dir  # noqa: E402
from tools.bench.run_local import _module  # noqa: E402
from zebra.commands.phenotype import build_consensus, gene_consensus, _ranked_items  # noqa: E402

KS = (1, 3, 10, 25)


def pos_in(ids, truth):
    return next((i for i, x in enumerate(ids, 1) if x in truth), None)


def tie_block(scores, pos):
    """(first position, size) of the block of list entries sharing the score at `pos` (1-based)."""
    if pos is None or not scores or scores[pos - 1] is None:
        return None, None
    s0 = scores[pos - 1]
    same = [i for i, v in enumerate(scores, 1) if v is not None and round(v, 6) == round(s0, 6)]
    return min(same), len(same)


def evaluate(method, module, hpo_dir, top, params=None):
    import os

    if hpo_dir:
        os.environ["ZEBRA_HPO_DIR"] = hpo_dir
    mod = _module(module)
    idx = mod.load()
    if not isinstance(idx, mod.Index):  # a pickle written by another copy of the module
        idx = mod.load(rebuild=True)
    kw = {}
    if not module and not method and not params:
        # by default, exactly what `zebra phenotype rank` runs
        from zebra.commands.phenotype import LOCAL_METHOD, LOCAL_PARAMS

        method, params = LOCAL_METHOD, dict(LOCAL_PARAMS)
    if method:
        kw["method"] = method
    if params:
        kw["params"] = params
    cases = {c["case"]: c for c in read_jsonl(work_dir() / "cases.jsonl")}
    remote = {}
    for src in ("monarch", "pubcasefinder"):
        path = work_dir() / f"remote_{src}.jsonl"
        remote[src] = {r["case"]: r for r in read_jsonl(path)} if path.exists() else {}
    rows = []
    evaluate.config = {"method": method, "params": params, "module": module}
    sssom = mapping.rows_by_object()
    mondo_with_omim = {r["subject"] for obj, rs in sssom.items() if obj.startswith("OMIM:") for r in rs}
    for name in sorted(set(remote["monarch"]) | set(remote["pubcasefinder"])):
        case = cases[name]
        truth = case["truth"]
        tgenes = {g.upper() for g in case["genes"]}
        row = {"case": name, "disease": case["disease"], "bucket": case["bucket"], "leaked": case["leaked"],
               "truth_has_mondo": any(t in sssom for t in truth),
               "has_gene": bool(tgenes), "in_pcf_sample": name in remote["pubcasefinder"],
               "in_monarch_sample": name in remote["monarch"]}
        try:
            res = mod.rank(idx, case["present"], case["excluded"], top=50, db=("OMIM", "ORPHA"), **kw)
        except ValueError as err:  # e.g. only onset/inheritance terms: the CLI skips local with a warning
            row["local_error"] = str(err)
            res = {"diseases": [], "genes": []}
        local_d = [{"id": d["disease"], "name": d["name"], "rank": i, "score": d["score"]}
                   for i, d in enumerate(res["diseases"], 1)]
        local_g = [{"symbol": g["gene"], "rank": i} for i, g in enumerate(
            [g for g in res["genes"] if g["gene"] not in ("-", "")], 1)]
        row["local"] = mapping.first_match([d["id"] for d in local_d], truth)
        row["local_tie_first"], row["local_tied"] = tie_block([d["score"] for d in local_d], row["local"])
        row["local_gene"] = pos_in([g["symbol"].upper() for g in local_g], tgenes) if tgenes else None
        per = {"local": {"diseases": local_d[:top], "genes": local_g[:top]}}
        m = remote["monarch"].get(name)
        if m is not None:
            if "diseases" in m:
                items = [{"source": "monarch", "id": h["id"], "xref": [x for x in h["xref"] if x.split(":")[0] in ("OMIM", "ORPHA")]}
                         for h in m["diseases"]]
                row["monarch"] = mapping.first_match(None, truth, items)
                row["monarch_tie_first"], row["monarch_tied"] = tie_block([h.get("score") for h in m["diseases"]],
                                                                          row["monarch"])
                # top-10 MONDO hits that no exact link ties to any OMIM id: they can never count as the truth
                row["monarch_unlinkable_top10"] = sum(1 for h in m["diseases"][:10]
                                                      if not any(x.startswith("OMIM:") for x in h["xref"])
                                                      and h["id"] not in mondo_with_omim)
                per["monarch"] = {"diseases": [{"id": h["id"], "name": h["name"], "rank": i, "xref": it["xref"]}
                                               for i, (h, it) in enumerate(zip(m["diseases"], items), 1)][:top]}
            else:
                row["monarch_error"] = m.get("diseases_error")
            if "genes" in m:
                syms = [(h.get("symbol") or h.get("name") or "").upper() for h in m["genes"]]
                row["monarch_gene"] = pos_in(syms, tgenes) if tgenes else None
                per.setdefault("monarch", {"diseases": []})["genes"] = [{"symbol": s, "rank": i} for i, s in enumerate(syms, 1)][:top]
        p = remote["pubcasefinder"].get(name)
        if p is not None:
            if "diseases" in p:
                ids = [h["id"] for h in p["diseases"]]
                row["pubcasefinder"] = pos_in(ids, truth)
                if row["pubcasefinder"]:
                    served = p["diseases"][row["pubcasefinder"] - 1]["served_rank"]
                    row["pubcasefinder_tie_first"] = served
                    row["pubcasefinder_tied"] = sum(1 for h in p["diseases"] if h["served_rank"] == served)
                # the CLI's consensus orders by the rank PubCaseFinder serves (ties included)
                per["pubcasefinder"] = {"diseases_omim": [{"id": h["id"], "name": h["name"], "rank": h["served_rank"]}
                                                          for h in p["diseases"]][:top],
                                        "diseases_orphanet": []}
            else:
                row["pubcasefinder_error"] = p.get("diseases_error")
            if "genes" in p:
                syms = [(h.get("symbol") or "").upper() for h in p["genes"]]
                row["pubcasefinder_gene"] = pos_in(syms, tgenes) if tgenes else None
                per.setdefault("pubcasefinder", {"diseases_omim": [], "diseases_orphanet": []})["genes"] = \
                    [{"symbol": s, "rank": h["served_rank"]} for s, h in zip(syms, p["genes"])][:top]
        # consensus over the sources this case has (CLI rule: in the top N of >= 2 sources)
        for label, sources in (("consensus_local_monarch", ("local", "monarch")),
                               ("consensus_all", ("local", "monarch", "pubcasefinder"))):
            if not all(s in per and per[s] for s in sources):
                continue
            sub = {s: per[s] for s in sources}
            for s in sub.values():
                s.setdefault("diseases", [])
                s.setdefault("diseases_omim", [])
                s.setdefault("diseases_orphanet", [])
            items = _ranked_items(sub)
            keys = mapping.keys_for(items, truth)
            cons = build_consensus(items, keys)
            tkeys = {keys[t][0] for t in truth}
            row[label] = next((r["order"] for r in cons if r["key"] in tkeys), None)
            row[label + "_len"] = len(cons)
            gcons = gene_consensus(sub)
            row[label + "_gene"] = next((i for i, g in enumerate(gcons, 1) if g["symbol"] in tgenes), None) if tgenes else None
        rows.append(row)
    return rows


def table(rows, top):
    out = {}
    msample = [r for r in rows if r["in_monarch_sample"]]
    psample = [r for r in rows if r["in_pcf_sample"]]
    m_cols = ("local", "monarch", "consensus_local_monarch")
    p_cols = ("local", "monarch", "pubcasefinder", "consensus_local_monarch", "consensus_all")
    for sample_name, sample, cols in (
            ("monarch sample", msample, m_cols),
            ("monarch sample, leaked", [r for r in msample if r["leaked"]], m_cols),
            ("monarch sample, not leaked", [r for r in msample if not r["leaked"]], m_cols),
            ("pubcasefinder sample", psample, p_cols),
            ("pubcasefinder sample, not leaked", [r for r in psample if not r["leaked"]], p_cols)):
        block = {}
        for c in cols:
            block[c] = {"disease": topk_with_ci(sample, c), "gene": topk_with_ci([r for r in sample if r["has_gene"]],
                                                                                  c + "_gene")}
            errs = sum(1 for r in sample if r.get(c + "_error"))
            if errs:
                block[c]["errors"] = errs
        out[sample_name] = block
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--method")
    ap.add_argument("--module")
    ap.add_argument("--hpo-dir")
    ap.add_argument("--params", default="{}")
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--label", default="remote")
    ap.add_argument("--json")
    args = ap.parse_args(argv)
    rows = evaluate(args.method, args.module, args.hpo_dir, args.top, json.loads(args.params))
    res = {"label": args.label, "local_config": getattr(evaluate, "config", None), "module": args.module,
           "top_per_source_for_consensus": args.top,
           "cases": len(rows), "tables": table(rows, args.top),
           "truth_without_mondo_exactmatch": sum(1 for r in rows if not r["truth_has_mondo"]),
           "monarch_top10_unlinkable_hits": sum(r.get("monarch_unlinkable_top10", 0) for r in rows),
           "monarch_top10_hits": 10 * sum(1 for r in rows if "monarch" in r),
           "rows": rows}
    print(f"truth OMIM ids with no MONDO exactMatch: {res['truth_without_mondo_exactmatch']} of {len(rows)} cases; "
          f"Monarch top-10 hits with no OMIM link: {res['monarch_top10_unlinkable_hits']} of {res['monarch_top10_hits']}")
    for sample, block in res["tables"].items():
        print(f"## {sample}")
        print("| list | n | top-1 | top-3 | top-10 | top-25 | gene n | gene top-1 | gene top-10 |")
        print("|---|---|---|---|---|---|---|---|---|")
        for c, b in block.items():
            d, g = b["disease"], b["gene"]

            def f(x, k):
                t = x.get(f"top{k}")
                return "-" if not t else (f"{t['rate'] * 100:.1f} ({t['ci95_wilson'][0] * 100:.0f}-"
                                          f"{t['ci95_wilson'][1] * 100:.0f}; ties {t['ties_fair'] * 100:.1f})")
            print(f"| {c} | {d['n']} | {f(d, 1)} | {f(d, 3)} | {f(d, 10)} | {f(d, 25)} | {g['n']} | {f(g, 1)} | {f(g, 10)} |"
                  + (f" errors {b['errors']}" if b.get("errors") else ""))
    if args.json:
        Path(args.json).write_text(json.dumps(res, indent=1) + "\n", "utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
