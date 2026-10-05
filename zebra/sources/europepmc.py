"""Europe PMC REST: literature search, PMID metadata, abstracts.

API: https://www.ebi.ac.uk/europepmc/webservices/rest (no key). Queries use
Europe PMC syntax and are passed through unchanged, e.g.
`"Dravet syndrome" AND (case report)`, `SCN1A AND PUB_YEAR:[2020 TO 2026]`.
Results keep Europe PMC's order (relevance unless a sort is asked for).
"""

from __future__ import annotations

import html
import re
from typing import Any, Dict, List, Optional, Sequence

from zebra.core import Outcome
from zebra.http import get_json, source_record

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
SORTS = {"relevance": None, "date": "P_PDATE_D desc", "cited": "CITED desc"}
_FREE = ("OA", "F", "R")  # open access, free, free to read (registration may apply)
_TAG_RE = re.compile(r"<[^>]+>")
PMID_RE = re.compile(r"^\d{1,9}$")
PMCID_RE = re.compile(r"^PMC\d+$", re.I)


def _clean(text: Optional[str]) -> Optional[str]:
    if not text:
        return text
    return " ".join(html.unescape(_TAG_RE.sub("", text)).split())


def _first_author(rec: Dict[str, Any]) -> Optional[str]:
    authors = ((rec.get("authorList") or {}).get("author")) or []
    names = [a.get("fullName") or a.get("collectiveName") for a in authors if isinstance(a, dict)]
    names = [n for n in names if n]
    if not names:
        s = (rec.get("authorString") or "").strip().rstrip(".")
        names = [n.strip() for n in s.split(",") if n.strip()] if s else []
    if not names:
        return None
    return names[0] + (" et al." if len(names) > 1 else "")


def _pub_types(rec: Dict[str, Any]) -> List[str]:
    lst = (rec.get("pubTypeList") or {}).get("pubType")
    if isinstance(lst, list):
        return [str(x) for x in lst]
    if isinstance(rec.get("pubType"), str):
        return [p.strip() for p in rec["pubType"].split(";") if p.strip()]
    return []


def _journal(rec: Dict[str, Any]) -> Optional[str]:
    j = ((rec.get("journalInfo") or {}).get("journal")) or {}
    return j.get("medlineAbbreviation") or j.get("isoabbreviation") or j.get("title") or rec.get("journalTitle")


def _full_text_urls(rec: Dict[str, Any]) -> List[Dict[str, str]]:
    out = []
    for u in (rec.get("fullTextUrlList") or {}).get("fullTextUrl") or []:
        if u.get("availabilityCode") in _FREE and u.get("url"):
            out.append({"url": u["url"], "site": u.get("site"), "style": u.get("documentStyle"),
                        "availability": u.get("availability") or u.get("availabilityCode")})
    return out


def _best_full_text(urls: List[Dict[str, str]]) -> Optional[str]:
    rank = {("Europe_PMC", "html"): 0, ("PubMedCentral", "html"): 1}
    best = sorted(urls, key=lambda u: (rank.get((u.get("site"), u.get("style")), 2), u.get("style") != "html"))
    return best[0]["url"] if best else None


def article_url(rec: Dict[str, Any]) -> str:
    src = rec.get("source") or "MED"
    rid = rec.get("id") or rec.get("pmid")
    return f"https://europepmc.org/article/{src}/{rid}"


def parse_hit(rec: Dict[str, Any]) -> Dict[str, Any]:
    """One Europe PMC result (lite or core) → the fields that matter for citing it."""
    urls = _full_text_urls(rec)
    hit: Dict[str, Any] = {
        "id": rec.get("id"),
        "source": rec.get("source"),
        "pmid": rec.get("pmid"),
        "pmcid": rec.get("pmcid"),
        "doi": rec.get("doi"),
        "title": _clean(rec.get("title")),
        "authors": _first_author(rec),
        "journal": _journal(rec),
        "year": rec.get("pubYear"),
        "pubType": _pub_types(rec),
        "isOpenAccess": rec.get("isOpenAccess") == "Y",
        "citedByCount": rec.get("citedByCount"),
        "fullTextUrl": _best_full_text(urls),
        "url": article_url(rec),
    }
    if not hit["fullTextUrl"] and hit["pmcid"] and rec.get("inEPMC") == "Y":
        # lite results carry no URL list; this is the page core results point at
        hit["fullTextUrl"] = f"https://europepmc.org/articles/{hit['pmcid']}"
    return hit


def _search_raw(query: str, page_size: int, sort: Optional[str], result_type: str):
    params: Dict[str, Any] = {"query": query, "format": "json", "resultType": result_type,
                              "pageSize": max(1, min(page_size, 100))}
    if sort:
        params["sort"] = SORTS.get(sort, sort)
    return get_json(f"{BASE}/search", source="Europe PMC", params=params, cache_ttl=86400, timeout=60)


def search(query: str, limit: int = 15, sort: Optional[str] = None, result_type: str = "core") -> Outcome:
    """Search Europe PMC; `sort` is relevance (default), date, cited, or a raw Europe PMC sort string."""
    query = (query or "").strip()
    if not query:
        raise ValueError("empty Europe PMC query")
    resp = _search_raw(query, limit, sort, result_type)
    data = resp.json()
    hits = [parse_hit(r) for r in (data.get("resultList") or {}).get("result") or []][:limit]
    result = {"query": query, "hitCount": data.get("hitCount"), "sort": sort or "relevance", "hits": hits}
    return Outcome(result, sources=[source_record("Europe PMC", query, resp, note=f"{len(hits)} of {data.get('hitCount')} hits")])


def by_pmids(pmids: Sequence[Any], result_type: str = "core") -> Outcome:
    """Metadata for PMIDs, returned in the order given; PMIDs Europe PMC does not know are reported."""
    ids = []
    for p in pmids:
        s = str(p).strip()
        if PMID_RE.match(s) and s not in ids:
            ids.append(s)
    found: Dict[str, Dict[str, Any]] = {}
    sources = []
    for i in range(0, len(ids), 25):
        chunk = ids[i:i + 25]
        query = "EXT_ID:(" + " OR ".join(chunk) + ") AND SRC:MED"
        resp = _search_raw(query, len(chunk), None, result_type)
        for rec in (resp.json().get("resultList") or {}).get("result") or []:
            if rec.get("pmid"):
                found[str(rec["pmid"])] = parse_hit(rec)
        sources.append(source_record("Europe PMC", f"PMIDs {chunk[0]}..{chunk[-1]} ({len(chunk)})", resp))
    missing = [p for p in ids if p not in found]
    warnings = [f"Europe PMC has no record for PMID {', '.join(missing)}"] if missing else []
    return Outcome([found[p] for p in ids if p in found], sources=sources, warnings=warnings)


def abstract(pmid: str) -> Outcome:
    """Title, abstract, keywords and full-text links for one PMID (or PMCID)."""
    pid = str(pmid).strip()
    if PMID_RE.match(pid):
        query = f"EXT_ID:{pid} AND SRC:MED"
    elif PMCID_RE.match(pid):
        query = f"PMCID:{pid.upper()}"
    else:
        raise ValueError(f"{pmid!r} is not a PMID or PMCID")
    resp = _search_raw(query, 1, None, "core")
    results = (resp.json().get("resultList") or {}).get("result") or []
    if not results:
        return Outcome({"id": pid, "found": False}, sources=[source_record("Europe PMC", pid, resp)],
                       warnings=[f"Europe PMC has no record for {pid}"])
    rec = results[0]
    out = parse_hit(rec)
    out["found"] = True
    out["abstract"] = _clean(rec.get("abstractText"))
    out["keywords"] = ((rec.get("keywordList") or {}).get("keyword") or [])[:15]
    out["fullTextUrls"] = _full_text_urls(rec)[:4]
    warnings = [] if out["abstract"] else [f"{pid}: Europe PMC holds no abstract"]
    return Outcome(out, sources=[source_record("Europe PMC", pid, resp)], warnings=warnings)
