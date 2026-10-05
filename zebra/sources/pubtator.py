"""PubTator3 (NCBI): entity autocomplete and entity-based article search.

API: https://www.ncbi.nlm.nih.gov/research/pubtator3-api (no key; NCBI pacing
applies). Entities are PubTator ids such as @GENE_SCN1A,
@VARIANT_p.R712X_SCN1A_human, @DISEASE_Epilepsies_Myoclonic. Disease names map
to MeSH concepts, which can be broader than the name typed (Dravet syndrome →
"Epilepsies, Myoclonic"): the `match` field says what was matched.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

from zebra.core import Outcome
from zebra.http import get_json, source_record

BASE = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
CONCEPTS = ("GENE", "VARIANT", "DISEASE", "CHEMICAL", "SPECIES", "CELLLINE")
PAGE_SIZE = 10  # fixed by the service
MAX_PAGES = 3


def autocomplete(text: str, concept: Optional[str] = None, limit: int = 5) -> Outcome:
    """Entity ids for a name: gene symbol, variant (e.g. 'SCN1A p.R712*', rsID) or disease."""
    text = (text or "").strip()
    if not text:
        raise ValueError("empty PubTator query")
    params: Dict[str, Any] = {"query": text, "limit": max(1, min(limit, 20))}
    if concept:
        concept = concept.upper()
        if concept not in CONCEPTS:
            raise ValueError(f"concept must be one of {', '.join(CONCEPTS)}")
        params["concept"] = concept
    resp = get_json(f"{BASE}/entity/autocomplete/", source="PubTator3", params=params, cache_ttl=14 * 86400)
    data = resp.json()
    rows = data if isinstance(data, list) else []
    entities = [{
        "id": e.get("_id"),
        "name": e.get("name"),
        "type": e.get("biotype"),
        "db": e.get("db"),
        "db_id": e.get("db_id"),
        "description": e.get("description"),
        "match": _strip_marks(e.get("match")),
    } for e in rows if isinstance(e, dict) and e.get("_id")][:limit]
    return Outcome({"query": text, "concept": concept, "entities": entities},
                   sources=[source_record("PubTator3 autocomplete", text, resp)])


def _strip_marks(s: Optional[str]) -> Optional[str]:
    return s.replace("<m>", "").replace("</m>", "") if s else s


def _hit(r: Dict[str, Any]) -> Dict[str, Any]:
    authors = r.get("authors") or []
    date = r.get("date") or ""
    return {
        "pmid": str(r["pmid"]) if r.get("pmid") is not None else r.get("_id"),
        "pmcid": r.get("pmcid"),
        "doi": r.get("doi"),
        "title": r.get("title"),
        "authors": (authors[0] + (" et al." if len(authors) > 1 else "")) if authors else None,
        "journal": r.get("journal"),
        "year": date[:4] or None,
        "url": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('pmid')}/" if r.get("pmid") else None,
    }


def search(text: str, limit: int = 10) -> Outcome:
    """Articles for an entity id (or free text / boolean of entity ids), in PubTator's order."""
    text = (text or "").strip()
    if not text:
        raise ValueError("empty PubTator search")
    limit = max(1, limit)
    pages = min(MAX_PAGES, math.ceil(limit / PAGE_SIZE))
    hits: List[Dict[str, Any]] = []
    sources = []
    count = None
    for page in range(1, pages + 1):
        resp = get_json(f"{BASE}/search/", source="PubTator3", params={"text": text, "page": page}, cache_ttl=86400)
        data = resp.json()
        count = data.get("count", count)
        sources.append(source_record("PubTator3 search", f"{text} p{page}", resp))
        rows = data.get("results") or []
        hits.extend(_hit(r) for r in rows if isinstance(r, dict))
        if len(hits) >= limit or page >= (data.get("total_pages") or 0):
            break
    warnings = []
    if limit > MAX_PAGES * PAGE_SIZE:
        warnings.append(f"PubTator3 returns {PAGE_SIZE} per page; capped at {MAX_PAGES * PAGE_SIZE}")
    return Outcome({"query": text, "count": count, "hits": hits[:limit]}, sources=sources, warnings=warnings)


def resolve_entity(text: str, concept: str) -> Outcome:
    """Best entity for a gene symbol / variant / disease: an exact name match first, else the service's first."""
    got = autocomplete(text, concept=concept, limit=8)
    ents = got.result["entities"]
    pick = None
    if concept.upper() == "GENE":
        pick = next((e for e in ents if (e.get("name") or "").upper() == text.strip().upper()
                     and e.get("db") == "ncbi_gene"), None)
    if pick is None and ents:
        pick = ents[0]
    got.result = {"query": text, "concept": concept.upper(), "entity": pick, "alternatives": [e["id"] for e in ents if e is not pick][:4]}
    return got
