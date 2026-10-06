"""Shared helpers for tools/bench: paths, the case list, buckets, confidence intervals."""

from __future__ import annotations

import gzip
import json
import math
import os
import random
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from zebra.http import cache_dir  # noqa: E402

SEED = 20261006
KS = (1, 3, 10, 25)
BUCKETS = (("1-3", 1, 3), ("4-6", 4, 6), ("7-10", 7, 10), (">10", 11, 10 ** 9))
RESULTS = ROOT / "tools" / "bench" / "results"


def work_dir() -> Path:
    """Large per-run files (case list, per-case ranks): $ZEBRA_BENCH_DIR or $ZEBRA_CACHE_DIR/bench/work."""
    d = Path(os.environ.get("ZEBRA_BENCH_DIR") or (cache_dir() / "bench" / "work"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def bucket(n: int) -> str:
    for name, lo, hi in BUCKETS:
        if lo <= n <= hi:
            return name
    return "0"


def read_jsonl(path: Path) -> List[Dict[str, Any]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    out = []
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def write_jsonl(path: Path, rows: Iterable[Dict[str, Any]]) -> int:
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(path) + ".part")
    with opener(tmp, "wt", encoding="utf-8") as fh:
        for r in rows:
            fh.write(json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n")
            n += 1
    os.replace(tmp, path)
    return n


def wilson(k: int, n: int, z: float = 1.959964) -> Tuple[Optional[float], Optional[float]]:
    if n == 0:
        return None, None
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def topk_with_ci(rows: Sequence[Dict[str, Any]], key: str, ks: Sequence[int] = KS, cluster_key: str = "disease",
                 reps: int = 1000, seed: int = SEED) -> Dict[str, Any]:
    """Top-k rates for the 1-based positions in rows[key] (None = not found), with 95% intervals.

    Two intervals: Wilson (cases treated as independent) and a cluster bootstrap that
    resamples whole diseases with replacement. Cases of one disease -- often one family
    in one paper -- are not independent, so the bootstrap interval is the one to quote.
    """
    clusters: Dict[str, List[float]] = {}
    fair: Dict[str, List[float]] = {}
    for r in rows:
        pos = r.get(key)
        c = clusters.setdefault(str(r[cluster_key]), [0] + [0] * len(ks))
        c[0] += 1
        f = fair.setdefault(str(r[cluster_key]), [0.0] * len(ks))
        first, tied = r.get(key + "_tie_first"), r.get(key + "_tied")
        for i, k in enumerate(ks):
            if pos is not None and pos <= k:
                c[i + 1] += 1
            # ties counted fairly: the truth equally likely at any position of its tie block
            if pos is not None and first and tied:
                f[i] += min(1.0, max(0.0, (k - first + 1) / tied))
            elif pos is not None and pos <= k:
                f[i] += 1.0
    n = sum(c[0] for c in clusters.values())
    out: Dict[str, Any] = {"n": n, "diseases": len(clusters)}
    if n == 0:
        return out
    totals = [sum(c[i + 1] for c in clusters.values()) for i in range(len(ks))]
    keys = sorted(clusters)
    rng = random.Random(seed)
    boot: List[List[float]] = [[] for _ in ks]
    for _ in range(reps):
        sn = 0
        sh = [0] * len(ks)
        for _j in range(len(keys)):
            c = clusters[keys[rng.randrange(len(keys))]]
            sn += c[0]
            for i in range(len(ks)):
                sh[i] += c[i + 1]
        for i in range(len(ks)):
            boot[i].append(sh[i] / sn if sn else 0.0)
    for i, k in enumerate(ks):
        vals = sorted(boot[i])
        lo, hi = wilson(totals[i], n)
        # macro: the mean over diseases of each disease's own rate (a 462-case family weighs as one)
        macro = sum(clusters[c][i + 1] / clusters[c][0] for c in keys) / len(keys)
        out[f"top{k}"] = {"rate": round(totals[i] / n, 4), "hits": totals[i],
                          "ci95_cluster": [round(vals[int(0.025 * (reps - 1))], 4), round(vals[int(0.975 * (reps - 1))], 4)],
                          "ci95_wilson": [round(lo, 4), round(hi, 4)],
                          "macro_by_disease": round(macro, 4),
                          "ties_fair": round(sum(f[i] for f in fair.values()) / n, 4)}
    return out


def version_key(name: str):
    """Sort release tags numerically: 0.1.3 < 0.1.27, v2026-09-01 by its numbers."""
    import re as _re

    return [int(x) for x in _re.findall(r"\d+", name)] or [0]


def hpo_files(hpo_dir: Optional[str] = None) -> Dict[str, Any]:
    """sha256, size and stated version of the HPO release files a run used."""
    import hashlib

    from zebra import hpo_local

    d = Path(hpo_dir) if hpo_dir else hpo_local.data_dir()
    out: Dict[str, Any] = {"dir": str(d)}
    for name in hpo_local.FILES:
        p = d / name
        if not p.exists():
            continue
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        out[name] = {"bytes": p.stat().st_size, "sha256": h.hexdigest()}
        if name == "phenotype.hpoa":
            with open(p, encoding="utf-8") as fh:
                for line in fh:
                    if not line.startswith("#"):
                        break
                    if line.lower().startswith(("#version", "#date")):
                        out[name]["header"] = out[name].get("header", []) + [line.strip()]
    return out


def git_sha() -> Optional[str]:
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=10)
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "zebra"], capture_output=True,
                               text=True, timeout=10)
        return out.stdout.strip() + ("+uncommitted-zebra-changes" if dirty.stdout.strip() else "")
    except (OSError, subprocess.SubprocessError):
        return None
