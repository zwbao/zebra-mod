"""Exact disease-id joins for scoring: MONDO's SSSOM exactMatch rows, fed to zebra's own join.

zebra's consensus joins ids from different namespaces only through exact
mappings (`zebra.commands.phenotype.resolve_keys`: Monarch exactMatch rows
MONDO -> OMIM/ORPHA, plus the MONDO xrefs Monarch returns with each hit). The CLI
fetches those rows from Monarch's /mappings endpoint per query; the benchmark
scores 10,000 cases, so it reads the same kind of rows from MONDO's own SSSOM file
(fetched and pinned by tools/bench/fetch_store.py) and passes them to the very same
`resolve_keys`. A hit counts as the right disease only when it is the truth id itself
or is joined to it by these exact links; a broader or narrower class is a miss.
"""

from __future__ import annotations

import csv
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from tools.bench.fetch_store import default_dir
from zebra.commands.phenotype import resolve_keys


def sssom_path() -> Path:
    base = default_dir().parent / "mondo"
    from tools.bench.common import version_key

    tags = sorted((p for p in base.iterdir() if (p / "manifest.json").exists()), key=lambda p: version_key(p.name)) \
        if base.exists() else []
    if not tags:
        raise SystemExit("MONDO SSSOM not fetched: run tools/bench/fetch_store.py")
    return tags[-1] / "mondo.sssom.tsv"


def manifest() -> Dict[str, Any]:
    return json.loads((sssom_path().parent / "manifest.json").read_text("utf-8"))


@lru_cache(maxsize=1)
def rows_by_object() -> Dict[str, List[Dict[str, Any]]]:
    out: Dict[str, List[Dict[str, Any]]] = {}
    with open(sssom_path(), encoding="utf-8") as fh:
        reader = csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t")
        for r in reader:
            if r.get("predicate_id") != "skos:exactMatch":
                continue
            subj, obj = r.get("subject_id") or "", r.get("object_id") or ""
            if not subj.startswith("MONDO:"):
                continue
            if obj.startswith("Orphanet:"):
                obj = "ORPHA:" + obj.split(":", 1)[1]
            elif not obj.startswith("OMIM:"):
                continue
            out.setdefault(obj, []).append({"subject": subj, "subject_label": r.get("subject_label") or "",
                                            "predicate": "skos:exactMatch", "object": obj})
    return out


def keys_for(items: Sequence[Dict[str, Any]], extra_ids: Sequence[str] = ()) -> Dict[str, Any]:
    """Group key (MONDO where an exact link exists) for every item id and every extra id."""
    allitems = list(items) + [{"source": "truth", "id": i} for i in extra_ids]
    index = rows_by_object()
    rows = []
    for it in allitems:
        rows.extend(index.get(it["id"], ()))
    return resolve_keys(allitems, rows, {})


def first_match(hit_ids: Sequence[str], truth: Sequence[str], items: Optional[Sequence[Dict[str, Any]]] = None) -> Optional[int]:
    """1-based position of the first hit joined to any truth id, or None."""
    its = list(items) if items is not None else [{"source": "local", "id": h} for h in hit_ids]
    keys = keys_for(its, truth)
    tkeys = {keys[t][0] for t in truth}
    for pos, it in enumerate(its, 1):
        if keys[it["id"]][0] in tkeys or it["id"] in truth:
            return pos
    return None
