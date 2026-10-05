"""PubCaseFinder (DBCLS): phenotype-driven ranking of diseases (OMIM, Orphanet) and genes.

API: https://pubcasefinder.dbcls.jp/api/pcf_get_ranked_list?target=omim|orphanet|gene&format=tsv&hpo_id=HP:..,HP:..
(spec: https://pubcasefinder.dbcls.jp/static/data/api/swagger-pubcasefinder.json). No key.

Usage limits published by DBCLS: at most 10 requests/minute, 100/hour, 1,000/day;
heavy jobs off-peak. One `rank` costs three requests (omim, orphanet, gene) and is
cached for a week. The service returns the whole ranked list (~8,000 rows), so
the TSV format is used (JSON is ~10 MB and four times slower). It takes HPO ids
only: excluded phenotypes are not supported.
"""

from __future__ import annotations

from typing import Any, Dict, List, Sequence

from zebra.core import Outcome
from zebra.http import request, source_record

BASE = "https://pubcasefinder.dbcls.jp/api"
TARGETS = ("omim", "orphanet", "gene")
WEB = "https://pubcasefinder.dbcls.jp"


def parse_tsv(text: str, target: str, limit: int) -> List[Dict[str, Any]]:
    """Rows as served (ranks may tie: 1, 1, 3, ...); disease ids OMIM:/ORPHA:, genes GENEID: + symbol."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    if "Rank" not in header[0]:
        raise ValueError(f"PubCaseFinder TSV without a header row: {lines[0][:80]!r}")
    out: List[Dict[str, Any]] = []
    for ln in lines[1:]:
        cols = ln.split("\t")
        row = dict(zip(header, cols))
        try:
            rank = int(row.get("Rank", "0"))
            score = round(float(row.get("Score", "nan")), 4)
        except ValueError:
            continue
        matched = [t for t in (row.get("Matched_Phenotype") or "").split(",") if t]
        if target == "gene":
            item = {"rank": rank, "score": score, "gene_id": row.get("NCBI_Gene_ID"),
                    "symbol": row.get("HGNC_Gene_Symbol"), "matched": matched}
        else:
            did = row.get("OMIM_ID") or row.get("ORPHA_ID")
            item = {"rank": rank, "score": score, "id": did, "name": row.get("Disease_Name"), "matched": matched,
                    "genes": [g for g in (row.get("Causative_Gene") or "").split(",") if g]}
        out.append(item)
        if len(out) >= limit:
            break
    return out


def ranked(hpo_ids: Sequence[str], target: str = "omim", limit: int = 15, timeout: float = 60.0) -> Outcome:
    """PubCaseFinder's ranked list for one target (omim, orphanet or gene), top `limit` rows.

    One retry at most (usage limits); typical time 10-20 s per list.
    """
    if target not in TARGETS:
        raise ValueError(f"target must be one of {', '.join(TARGETS)}")
    hpo = ",".join(hpo_ids)
    resp = request(f"{BASE}/pcf_get_ranked_list", source="PubCaseFinder",
                   params={"target": target, "format": "tsv", "hpo_id": hpo},
                   accept="text/tab-separated-values, text/plain, */*", timeout=timeout, retries=1,
                   cache_ttl=7 * 86400)
    rows = parse_tsv(resp.text, target, limit)
    return Outcome({"target": target, "hits": rows},
                   sources=[source_record("PubCaseFinder", f"{target}: {hpo}", resp,
                                          note="ranked list (pcf_get_ranked_list, TSV)")])
