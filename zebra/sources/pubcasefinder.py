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
from zebra.http import request
from zebra.sources import record as source_record
from zebra.sources import validated_text

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
            raw_score = row.get("Score")
            if raw_score is None or not str(raw_score).strip():
                continue  # a row with no score would serialise as NaN, which is not valid JSON
            score = round(float(raw_score), 4)
            if score != score or score in (float("inf"), float("-inf")):
                continue
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


def ranked(hpo_ids: Sequence[str], target: str = "omim", limit: int = 15, timeout: float = 60.0,
           retries: int = 1) -> Outcome:
    """PubCaseFinder's ranked list for one target (omim, orphanet or gene), top `limit` rows.

    One retry at most by default (usage limits); typical time 10-20 s per list. A caller that
    meters its own requests against the published limits passes retries=0, so that one call is
    at most one request (plus one re-request when a bad body came from the cache).
    """
    if target not in TARGETS:
        raise ValueError(f"target must be one of {', '.join(TARGETS)}")
    hpo = ",".join(hpo_ids)
    if not hpo:
        raise ValueError("PubCaseFinder needs at least one HPO id")

    def fetch(refresh=False):
        return request(f"{BASE}/pcf_get_ranked_list", source="PubCaseFinder",
                       params={"target": target, "format": "tsv", "hpo_id": hpo},
                       accept="text/tab-separated-values, text/plain, */*", timeout=timeout,
                       retries=max(0, min(int(retries), 1)),
                       # E-6: a re-request after a bad cached body skips the cache READ but still
                       # WRITES the good answer back; cache_ttl=0 used to skip both, so every later
                       # call replayed the bad entry and paid another request for the whole TTL
                       cache_ttl=7 * 86400, refresh=refresh)

    used = []

    def refetch():
        used.append(fetch(refresh=True))
        return used[-1]

    resp = fetch()
    # F10: an empty or HTML 200 was cached for 7 days and parsed to `hits: []`
    # with no warning, which reads as "PubCaseFinder matched no disease for this
    # patient". `Rank` is the first column of the real TSV header
    # (`Rank\tScore\tOMIM_ID\tDisease_Name\tMatched_Phenotype\tCausative_Gene`,
    # checked live 2026-10-06).
    text = validated_text(resp, "PubCaseFinder", must_contain="Rank", refetch=refetch)
    resp = used[-1] if used else resp  # provenance names the answer actually used
    rows = parse_tsv(text, target, limit)
    return Outcome({"target": target, "hits": rows},
                   sources=[source_record("PubCaseFinder", f"{target}: {hpo}", resp,
                                          note="ranked list (pcf_get_ranked_list, TSV)")])
