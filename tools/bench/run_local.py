#!/usr/bin/env python3
"""Score zebra's offline ranker (zebra.hpo_local.rank) on the benchmark cases.

    python3 -I tools/bench/run_local.py --label NAME [--split all|dev|heldout] [--options JSON]
                                        [--module FILE] [--hpo-dir DIR] [--no-excluded] [--workers N]

Per case, two rankings are made, exactly as zebra makes them:
  rank_omim  the truth's position when only OMIM diseases are ranked (the truth is always an
             OMIM id, so this is plain id equality: the primary number)
  rank_all   the position of the first hit that is the truth or joined to it by MONDO exactMatch,
             in the OMIM + Orphanet list `zebra phenotype rank` shows (secondary number)
  gene_rank  the position of the first causal gene in the gene list of that same ranking
Positions beyond 50 are recorded as null (not in the top 50).
`--options` passes keyword arguments to rank() (e.g. '{"method": "lr"}'); `--module` scores
another copy of hpo_local (how the shipped 0.1.0 baseline is reproduced: `git show
54c0b4e:zebra/hpo_local.py > /tmp/hpo_local_010.py`).
Output: work_dir()/local_<label>.jsonl.gz
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
from multiprocessing import get_context
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.bench.common import read_jsonl, work_dir, write_jsonl  # noqa: E402

TOP = 50
_W = {}


def _module(path):
    if not path:
        from zebra import hpo_local

        return hpo_local
    spec = importlib.util.spec_from_file_location("hpo_local_alt", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["hpo_local_alt"] = mod  # dataclasses and pickle look the module up by name
    spec.loader.exec_module(mod)
    return mod


def _init(module_path, hpo_dir, options, no_excluded):
    if hpo_dir:
        os.environ["ZEBRA_HPO_DIR"] = hpo_dir
    mod = _module(module_path)
    idx = mod.load()
    if not isinstance(idx, mod.Index):  # a pickle written by another copy of the module: rebuild with this one
        idx = mod.load(rebuild=True)
    _W.update(mod=mod, idx=idx, options=options, no_excluded=no_excluded)
    from tools.bench import mapping

    mapping.rows_by_object()
    _W["mapping"] = mapping


def score_case(case):
    mod, idx, kw = _W["mod"], _W["idx"], _W["options"]
    excluded = [] if _W["no_excluded"] else case["excluded"]
    out = {k: case[k] for k in ("case", "disease", "split", "bucket", "leaked", "in_hpoa", "n_present")}
    out["n_excluded"] = len(case["excluded"])
    t0 = time.time()
    try:
        r_omim = mod.rank(idx, case["present"], excluded, top=TOP, db=("OMIM",), **kw)
        r_all = mod.rank(idx, case["present"], excluded, top=TOP, db=("OMIM", "ORPHA"), **kw)
    except ValueError as err:
        out["error"] = str(err)
        return out
    ids = [d["disease"] for d in r_omim["diseases"]]
    truth = case["truth"]
    out["rank_omim"] = next((i for i, d in enumerate(ids, 1) if d in truth), None)
    hit = next((d for d in r_omim["diseases"] if d["disease"] in truth), None)
    out["truth_tied"] = hit["tied_at_this_score"] if hit else None
    if hit:  # the block of positions sharing the truth's score: [first, first + tied - 1]
        out["rank_omim_tied"] = hit["tied_at_this_score"]
        out["rank_omim_tie_first"] = next(i for i, d in enumerate(r_omim["diseases"], 1)
                                          if round(d["score"], 4) == round(hit["score"], 4))
    out["top_ties"] = r_omim["ties"]
    out["rank_all"] = _W["mapping"].first_match([d["disease"] for d in r_all["diseases"]], truth)
    genes = [g["gene"] for g in r_all["genes"] if g["gene"] not in ("-", "")]
    tg = set(case["genes"])
    out["gene_rank"] = next((i for i, g in enumerate(genes, 1) if g in tg), None) if tg else None
    out["has_gene"] = bool(tg)
    out["seconds"] = round(time.time() - t0, 3)
    return out


def _sha(path):
    import hashlib

    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--split", default="all", choices=("all", "dev", "heldout"))
    ap.add_argument("--options", default="{}")
    ap.add_argument("--module")
    ap.add_argument("--hpo-dir")
    ap.add_argument("--no-excluded", action="store_true")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--cases", help="only these case names (comma list)")
    args = ap.parse_args(argv)
    cases = read_jsonl(work_dir() / "cases.jsonl")
    if args.split != "all":
        cases = [c for c in cases if c["split"] == args.split]
    if args.cases:
        want = set(args.cases.split(","))
        cases = [c for c in cases if c["case"] in want]
    if args.limit:
        cases = cases[: args.limit]
    options = json.loads(args.options)
    t0 = time.time()
    # fail here, not in the pool: a pool whose initializer raises respawns workers forever
    _init(args.module, args.hpo_dir, options, args.no_excluded)
    if cases:
        score_case(cases[0])
    ctx = get_context("spawn")
    with ctx.Pool(args.workers, initializer=_init,
                  initargs=(args.module, args.hpo_dir, options, args.no_excluded)) as pool:
        rows = list(pool.imap(score_case, cases, chunksize=8))
    path = work_dir() / f"local_{args.label}.jsonl.gz"
    write_jsonl(path, rows)
    from tools.bench.common import git_sha, hpo_files

    meta = {"label": args.label, "split": args.split, "options": options, "module": args.module,
            "module_sha256": _sha(args.module) if args.module else None, "git": git_sha(),
            "hpo": hpo_files(args.hpo_dir), "no_excluded": args.no_excluded, "cases": len(rows),
            "seconds": round(time.time() - t0, 1), "path": str(path)}
    (work_dir() / f"local_{args.label}.meta.json").write_text(json.dumps(meta, indent=2), "utf-8")
    n = len(rows)
    for key in ("rank_omim", "rank_all", "gene_rank"):
        vals = [r.get(key) for r in rows if not r.get("error") and (key != "gene_rank" or r.get("has_gene"))]
        m = len(vals)
        print(key, {k: round(sum(1 for v in vals if v and v <= k) / m, 4) for k in (1, 3, 10, 25)} if m else None, m)
    print(json.dumps(meta))
    return 0


if __name__ == "__main__":
    sys.exit(main())
