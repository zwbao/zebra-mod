"""Monarch Initiative API v3: diseases, genes, mappings and phenotype similarity search.

API: https://api-v3.monarchinitiative.org/v3/api (no key). OpenAPI at
https://api-v3.monarchinitiative.org/openapi.json.

Each public function fetches and returns an Outcome; the `parse_*` functions are
pure (payload -> compact dict) so they can be tested on captured responses.
Monarch spells Orphanet ids `Orphanet:N`; zebra uses `ORPHA:N` everywhere else,
so ids are normalised on the way out (`to_monarch_curie` goes the other way).
"""

from __future__ import annotations

import re
import urllib.parse
from typing import Any, Dict, Iterable, List, Optional, Sequence

from zebra.core import Outcome
from zebra.http import get_json, post_json
from zebra.sources import record as source_record

BASE = "https://api-v3.monarchinitiative.org/v3/api"
WEB = "https://monarchinitiative.org"
SEMSIM_MAX = 50  # the API's own cap on `limit`
GENE_DISEASE_CATEGORIES = ("biolink:CausalGeneToDiseaseAssociation", "biolink:CorrelatedGeneToDiseaseAssociation")
_HGNC_RE = re.compile(r"^HGNC:\d+$", re.I)


def normalize_curie(curie: Optional[str]) -> Optional[str]:
    """Monarch/MONDO spellings -> zebra's: Orphanet:N -> ORPHA:N; OMIM stays."""
    if not curie:
        return curie
    if curie.startswith("Orphanet:"):
        return "ORPHA:" + curie.split(":", 1)[1]
    return curie


def to_monarch_curie(curie: str) -> str:
    if curie.upper().startswith("ORPHA:"):
        return "Orphanet:" + curie.split(":", 1)[1]
    return curie


def web_url(curie: str) -> str:
    return f"{WEB}/{urllib.parse.quote(curie)}"


def _xref_groups(xrefs: Iterable[str]) -> Dict[str, List[str]]:
    groups: Dict[str, List[str]] = {}
    for x in xrefs or []:
        x = normalize_curie(x) or ""
        prefix = x.split(":", 1)[0]
        groups.setdefault(prefix, []).append(x)
    return groups


# ---------------------------------------------------------------- search / entity


def parse_search(data: Dict[str, Any], limit: int) -> List[Dict[str, Any]]:
    out = []
    for item in (data.get("items") or [])[:limit]:
        out.append({
            "id": item.get("id"),
            "name": item.get("name"),
            "category": item.get("category"),
            "symbol": item.get("symbol"),
            "taxon": item.get("in_taxon_label"),
            "description": _clip(item.get("description"), 300),
            "synonyms": (item.get("exact_synonym") or item.get("synonym") or [])[:8],
            "xref": [normalize_curie(x) for x in (item.get("xref") or [])][:30],
            "deprecated": bool(item.get("deprecated")),
            "url": web_url(item["id"]) if item.get("id") else None,
        })
    return out


def search(text: str, category: Optional[str] = None, limit: int = 10, match_type: Optional[str] = None,
           scope: Optional[str] = None) -> Outcome:
    """Free-text (or id) search. category e.g. biolink:Disease, biolink:Gene; scope e.g. rare_disease."""
    params: Dict[str, Any] = {"q": text, "limit": max(1, min(limit, 50))}
    if category:
        params["category"] = category
    if match_type:
        params["match_type"] = match_type
    if scope:
        params["scope"] = scope
    resp = get_json(f"{BASE}/search", source="Monarch", params=params, cache_ttl=14 * 86400)
    data = resp.json()
    return Outcome({"query": text, "total": data.get("total"), "hits": parse_search(data, limit)},
                   sources=[source_record("Monarch search", text, resp)])


def parse_entity(data: Dict[str, Any]) -> Dict[str, Any]:
    xrefs = [normalize_curie(x) for x in (data.get("xref") or [])]
    mappings = [normalize_curie(m.get("id")) for m in (data.get("mappings") or []) if m.get("id")]
    groups = _xref_groups(xrefs + [m for m in mappings if m not in xrefs])
    return {
        "id": data.get("id"),
        "name": data.get("name"),
        "category": data.get("category"),
        "symbol": data.get("symbol"),
        "description": data.get("description"),
        "synonyms": (data.get("exact_synonym") or data.get("synonym") or [])[:15],
        "deprecated": bool(data.get("deprecated")),
        "xref": xrefs,
        "exact_mappings": mappings,
        "omim": groups.get("OMIM", []) + groups.get("OMIMPS", []),
        "orpha": groups.get("ORPHA", []),
        "icd": groups.get("ICD10CM", []) + groups.get("ICD10WHO", []) + groups.get("icd11.foundation", []),
        "other_ids": {k: v for k, v in groups.items() if k not in ("OMIM", "OMIMPS", "ORPHA")},
        "subsets": data.get("subsets") or [],
        "phenotype_count": data.get("has_phenotype_count"),
        "url": web_url(data["id"]) if data.get("id") else None,
    }


def entity(entity_id: str) -> Outcome:
    """One entity (disease, gene, phenotype): name, description, synonyms, xrefs (OMIM/ORPHA/ICD...), mappings."""
    curie = to_monarch_curie(entity_id.strip())
    resp = get_json(f"{BASE}/entity/{urllib.parse.quote(curie, safe=':')}", source="Monarch", cache_ttl=14 * 86400,
                    ok_statuses=(200, 404))
    if resp.status == 404:
        return Outcome(None, sources=[source_record("Monarch entity", entity_id, resp)],
                       warnings=[f"Monarch has no entity {entity_id}"])
    return Outcome(parse_entity(resp.json()), sources=[source_record("Monarch entity", entity_id, resp)])


# ---------------------------------------------------------------- associations


def parse_associations(data: Dict[str, Any], side: str) -> List[Dict[str, Any]]:
    """side='subject' lists the genes of a disease query; side='object' the diseases of a gene query."""
    out = []
    for a in data.get("items") or []:
        other = a.get(side)
        out.append({
            "id": normalize_curie(other),
            "label": a.get(f"{side}_label"),
            "predicate": (a.get("predicate") or "").replace("biolink:", ""),
            "category": (a.get("category") or "").replace("biolink:", ""),
            "source": (a.get("primary_knowledge_source") or "").replace("infores:", ""),
            "original_predicate": a.get("original_predicate"),
            "original_id": normalize_curie(a.get(f"original_{side}")),
            "publications": (a.get("publications") or [])[:5],
            "for_id": normalize_curie(a.get("object" if side == "subject" else "subject")),
            "for_label": a.get("object_label" if side == "subject" else "subject_label"),
        })
    return out


def disease_genes(disease_id: str, limit: int = 100) -> Outcome:
    """Causal and correlated gene associations of a disease (MONDO id preferred), with predicate and source."""
    curie = to_monarch_curie(disease_id.strip())
    resp = get_json(f"{BASE}/association", source="Monarch",
                    params={"object": curie, "category": list(GENE_DISEASE_CATEGORIES), "limit": max(1, min(limit, 500))},
                    cache_ttl=14 * 86400)
    data = resp.json()
    rows = parse_associations(data, "subject")
    return Outcome({"disease": disease_id, "total": data.get("total"), "genes": rows},
                   sources=[source_record("Monarch gene-disease associations", disease_id, resp)])


def resolve_gene(symbol_or_hgnc: str) -> Outcome:
    """HGNC id for a human gene symbol (exact symbol match among human genes)."""
    q = symbol_or_hgnc.strip()
    if _HGNC_RE.match(q):
        return Outcome({"id": q.upper(), "symbol": None})
    resp = get_json(f"{BASE}/search", source="Monarch",
                    params={"q": q, "category": "biolink:Gene", "in_taxon_label": "Homo sapiens", "limit": 10},
                    cache_ttl=14 * 86400)
    hits = parse_search(resp.json(), 10)
    for h in hits:
        if (h.get("symbol") or h.get("name") or "").upper() == q.upper() and (h.get("id") or "").startswith("HGNC:"):
            return Outcome({"id": h["id"], "symbol": h.get("symbol") or h.get("name")},
                           sources=[source_record("Monarch search", q, resp)])
    raise ValueError(f"no human gene with symbol {q!r} in Monarch")


def gene_diseases(symbol_or_hgnc: str, limit: int = 100) -> Outcome:
    """Diseases a gene causes or is correlated with (OMIM, Orphanet, ClinGen ... via Monarch)."""
    gene = resolve_gene(symbol_or_hgnc)
    hgnc = gene.result["id"]
    resp = get_json(f"{BASE}/association", source="Monarch",
                    params={"subject": hgnc, "category": list(GENE_DISEASE_CATEGORIES), "limit": max(1, min(limit, 500))},
                    cache_ttl=14 * 86400)
    data = resp.json()
    out = Outcome({"gene": hgnc, "symbol": gene.result.get("symbol") or symbol_or_hgnc, "total": data.get("total"),
                   "diseases": parse_associations(data, "object")},
                  sources=gene.sources + [source_record("Monarch gene-disease associations", hgnc, resp)])
    return out


# ---------------------------------------------------------------- mappings


def parse_mappings(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    return [{
        "subject": normalize_curie(m.get("subject_id")),
        "subject_label": m.get("subject_label"),
        "predicate": m.get("predicate_id"),
        "object": normalize_curie(m.get("object_id")),
        "object_label": m.get("object_label"),
        "source": m.get("mapping_source"),
    } for m in data.get("items") or []]


def mappings(object_ids: Sequence[str] = (), subject_ids: Sequence[str] = (),
             predicate: Optional[str] = "skos:exactMatch") -> Outcome:
    """SSSOM mappings (MONDO <-> OMIM/Orphanet/...), batched: many ids in one request."""
    params: Dict[str, Any] = {"limit": 500}
    if object_ids:
        params["object_id"] = [to_monarch_curie(i) for i in object_ids]
    if subject_ids:
        params["subject_id"] = [to_monarch_curie(i) for i in subject_ids]
    if predicate:
        params["predicate_id"] = predicate
    if not object_ids and not subject_ids:
        return Outcome([])
    resp = get_json(f"{BASE}/mappings", source="Monarch", params=params, cache_ttl=30 * 86400)
    label = ",".join(list(object_ids)[:3] + list(subject_ids)[:3]) + ("…" if len(object_ids) + len(subject_ids) > 6 else "")
    return Outcome(parse_mappings(resp.json()), sources=[source_record("Monarch mappings", label, resp)])


# ---------------------------------------------------------------- semantic similarity


def parse_semsim(data: List[Dict[str, Any]], limit: int) -> List[Dict[str, Any]]:
    """Semsim search rows: subject = the disease/gene, object termset = the query phenotypes."""
    out = []
    for i, row in enumerate(data[:limit], 1):
        subj = row.get("subject") or {}
        sim = row.get("similarity") or {}
        matched = []
        for q, m in (sim.get("object_best_matches") or {}).items():
            matched.append({
                "query": m.get("match_source") or q,
                "query_label": m.get("match_source_label"),
                "match": m.get("match_target"),
                "match_label": m.get("match_target_label"),
                "via": m.get("match_subsumer"),
                "via_label": m.get("match_subsumer_label"),
                "exact": m.get("match_target") == (m.get("match_source") or q),
                "score": round(float(m.get("score") or 0.0), 3),
            })
        out.append({
            "rank": i,
            "id": normalize_curie(subj.get("id")),
            "name": subj.get("name"),
            "symbol": subj.get("symbol"),
            "score": round(float(row.get("score") or 0.0), 3),
            "xref": [normalize_curie(x) for x in (subj.get("xref") or [])],
            "matched": matched,
            "url": web_url(subj["id"]) if subj.get("id") else None,
        })
    return out


def semsim_rank(present: Sequence[str], limit: int = 15, group: str = "Human Diseases",
                metric: str = "ancestor_information_content") -> Outcome:
    """Rank diseases (or genes, group='Human Genes') by phenotype similarity to the HPO term set.

    The service has no notion of excluded phenotypes and caps `limit` at 50.
    """
    warnings: List[str] = []
    if limit > SEMSIM_MAX:
        warnings.append(f"Monarch semsim returns at most {SEMSIM_MAX} results; asked for {limit}")
        limit = SEMSIM_MAX
    body = {"termset": list(present), "group": group, "metric": metric, "limit": max(1, limit)}
    resp = post_json(f"{BASE}/semsim/search", body, source="Monarch semsim", cache_ttl=7 * 86400, timeout=60, retries=1)
    data = resp.json()
    if not isinstance(data, list):
        raise ValueError(f"semsim search returned {type(data).__name__}, expected a list")
    return Outcome({"group": group, "metric": metric, "hits": parse_semsim(data, limit)},
                   sources=[source_record("Monarch semsim search", f"{group}: {','.join(present)}", resp,
                                          note=f"metric {metric}, bidirectional")],
                   warnings=warnings)


def _clip(text: Optional[str], n: int) -> Optional[str]:
    if not text:
        return text
    return text if len(text) <= n else text[: n - 1].rstrip() + "…"
