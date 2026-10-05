"""HPO terms: the JAX ontology API, with the local release (if fetched) answered first.

APIs: https://ontology.jax.org/api/hp/ (terms, disease annotations) and, for
search, EBI OLS4 https://www.ebi.ac.uk/ols4/api/search (both keyless).

Search uses OLS4, not the JAX search endpoint (E9). On a fresh install, before
`zebra hpo fetch` has run, the JAX endpoint's own order never surfaces the term
asked for: measured on 2026-10-06, `?q=seizure&page=0..1&limit=50` returns
`totalCount` 100 and HP:0001250 *Seizure* appears on no page, so the first hits
are `HP:0100622 Maternal seizure` and `HP:0032665 Repeated focal motor
seizures`. The tool description tells the model to search HPO before recording
any phenotype, so that answer put the wrong term in a patient's record. OLS4
returns HP:0001250 first for the same query. Its hits are re-ranked here the
way `zebra.hpo_local.search` ranks the local release — exact label, then exact
synonym, then prefix, then word, then substring, shorter wording first and a
label before a synonym — so the online and offline backends agree on the order,
and the JAX endpoint is kept as a fallback through the same re-ranking.
"""

from __future__ import annotations

import re
import urllib.parse
from typing import Any, Dict, List

from zebra import hpo_local
from zebra.core import Outcome
from zebra.http import SourceError, get_json
from zebra.sources import record as source_record

BASE = "https://ontology.jax.org/api"
HPO_RE = re.compile(r"^HP:\d{7}$")


def _local():
    try:
        return hpo_local.load()
    except hpo_local.HpoDataMissing:
        return None


def term(hpo_id: str, prefer_local: bool = True) -> Outcome:
    """Verify one HPO id: name, definition, synonyms; obsolete/replaced when the local release knows."""
    hpo_id = hpo_id.strip()
    if not HPO_RE.match(hpo_id):
        raise ValueError(f"{hpo_id!r} is not an HPO id (HP:0000000)")
    idx = _local() if prefer_local else None
    if idx is not None:
        pid, note = idx.primary(hpo_id)
        result: Dict[str, Any] = {
            "id": hpo_id,
            "name": idx.names.get(pid) if pid else None,
            "primary": pid,
            "obsolete": hpo_id in idx.obsolete,
            "replaced_by": idx.obsolete.get(hpo_id),
            "note": note,
            "synonyms": idx.synonyms.get(pid or "", [])[:10],
            "backend": f"local HPO {idx.version}",
        }
        return Outcome(result, sources=[source_record("HPO", hpo_id, url=f"https://hpo.jax.org/browse/term/{hpo_id}",
                                                      note=f"local release {idx.version}")])
    resp = get_json(f"{BASE}/hp/terms/{urllib.parse.quote(hpo_id)}", source="HPO", cache_ttl=30 * 86400,
                    ok_statuses=(200, 404))
    if resp.status == 404:
        return Outcome({"id": hpo_id, "name": None, "obsolete": None, "note": "not found in HPO"},
                       sources=[source_record("HPO", hpo_id, resp)])
    data = resp.json()
    result = {
        "id": data.get("id", hpo_id),
        "name": data.get("name"),
        "definition": data.get("definition"),
        "synonyms": (data.get("synonyms") or [])[:10],
        "xrefs": (data.get("xrefs") or [])[:10],
        "obsolete": bool(data.get("obsolete")) if "obsolete" in data else None,
        "replaced_by": data.get("replacedBy"),
        "backend": "HPO API (ontology.jax.org)",
    }
    return Outcome(result, sources=[source_record("HPO", hpo_id, resp)])


OLS_SEARCH = "https://www.ebi.ac.uk/ols4/api/search"
OLS_FIELDS = "obo_id,label,synonym,ontology_name,is_defining_ontology,is_obsolete"
SEARCH_ROWS = 100  # fetch wide, then re-rank; the service's own order is not trusted


def _norm(text: str) -> str:
    """Fold case, punctuation and simple plurals, as `zebra.hpo_local` does for the local index."""
    t = re.sub(r"[^a-z0-9]+", " ", (text or "").lower()).strip()
    return re.sub(r"\b(\w{4,})s\b", r"\1", t)


def rank_hits(query: str, hits: List[Dict[str, Any]], limit: int = 10) -> List[Dict[str, Any]]:
    """Order candidate terms the way the local index does, and say what matched.

    Rank: 0 exact label, 1 exact synonym, 2 label/synonym starting with the
    query, 3 the query as a whole word, 4 the query as a substring. Within a
    rank the wording closest in length to the query wins, and the term's own
    label beats the same match on a synonym; a term matching on nothing is
    dropped only when something else matched, so a service that ranks loosely
    cannot push the exact term off the end of `limit`.
    """
    q = _norm(query)
    scored = []
    for i, h in enumerate(hits):
        label = h.get("label") or ""
        syns = [s for s in (h.get("synonyms") or []) if s]
        best = None
        for j, cand in enumerate([label] + syns):
            c = _norm(cand)
            if not c or not q:
                continue
            if c == q:
                rank = 0 if j == 0 else 1
            elif c.startswith(q):
                rank = 2
            elif re.search(rf"\b{re.escape(q)}\b", c):
                rank = 3
            elif q in c:
                rank = 4
            else:
                continue
            key = (rank, max(0, len(c) - len(q)), 0 if j == 0 else 1, len(c))
            if best is None or key < best[0]:
                best = (key, cand, j)
        if best is None:
            scored.append(((9, 0, 1, 0, i), h, None))
        else:
            scored.append((best[0] + (i,), h, best[1]))
    scored.sort(key=lambda s: s[0])
    out = []
    for key, h, matched in scored[:limit]:
        row = dict(h)
        row["matched"] = matched
        row["matched_on"] = (None if matched is None else
                             ("label" if key[2] == 0 else "synonym"))
        out.append(row)
    return out


def _ols_search(text: str, rows: int = SEARCH_ROWS):
    resp = get_json(OLS_SEARCH, source="OLS (HPO search)",
                    params={"q": text, "ontology": "hp", "rows": rows, "queryFields": "label,synonym",
                            "fieldList": OLS_FIELDS, "exact": "false"},
                    cache_ttl=14 * 86400)
    data = resp.json()
    docs = ((data.get("response") or {}).get("docs") or [])
    hits = [{"id": d.get("obo_id"), "label": d.get("label"), "synonyms": (d.get("synonym") or [])[:5]}
            for d in docs
            if str(d.get("obo_id") or "").startswith("HP:") and not d.get("is_obsolete")]
    return resp, hits


def _jax_search(text: str, rows: int = 50):
    resp = get_json(f"{BASE}/hp/search", source="HPO", params={"q": text, "page": 0, "limit": max(1, min(rows, 100))},
                    cache_ttl=14 * 86400)
    data = resp.json()
    hits = [{"id": t.get("id"), "label": t.get("name"), "synonyms": (t.get("synonyms") or [])[:5],
             "descendants": t.get("descendantCount")} for t in data.get("terms", [])]
    return resp, hits


def search(text: str, limit: int = 10) -> Outcome:
    """Online HPO term search: OLS4, re-ranked; JAX as a fallback, re-ranked the same way."""
    limit = max(1, int(limit))
    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []
    backend = None
    hits: List[Dict[str, Any]] = []
    try:
        resp, raw = _ols_search(text)
        sources.append(source_record("HPO search (EBI OLS4)", text, resp, note="ontology=hp, label+synonym fields"))
        hits, backend = raw, "EBI OLS4 (ontology=hp), re-ranked: exact label, exact synonym, prefix, word, substring"
    except (SourceError, ValueError, KeyError, TypeError, AttributeError) as err:
        warnings.append(f"OLS4 HPO search unavailable ({type(err).__name__}: {err}); fell back to the JAX HPO API, "
                        "whose own order does not put the exact term first")
    if not hits and text.isascii():
        # The JAX search endpoint rejects a non-ASCII query with a raw regex
        # error (`query: must match "^[a-zA-Z0-9 …`), so it is only tried for
        # ASCII input; a Chinese query needs the local release.
        try:
            resp, raw = _jax_search(text)
            sources.append(source_record("HPO search", text, resp))
            hits = raw
            backend = backend or "HPO API (ontology.jax.org), re-ranked"
        except (SourceError, ValueError, KeyError, TypeError, AttributeError) as err:
            warnings.append(f"HPO API search unavailable ({type(err).__name__}: {err})")
    elif not hits:
        warnings.append(f"'{text}' is not ASCII: the online HPO search covers English labels and synonyms only "
                        "(the JAX search endpoint rejects non-ASCII queries). Run `zebra hpo fetch` once to search "
                        "the Chinese HPO labels locally, or search the English term")
    ranked = rank_hits(text, hits, limit)
    if not ranked:
        warnings.append(f"no HPO term matched '{text}' online; `zebra hpo fetch` adds the local release, which also "
                        "searches Chinese labels and curated lay phrases")
    return Outcome({"query": text, "hits": ranked, "backend": backend, "candidates_considered": len(hits)},
                   sources=sources, warnings=warnings)


def disease_annotations(disease_id: str) -> Outcome:
    """Phenotypes annotated to a disease (OMIM:/ORPHA:), grouped by category."""
    resp = get_json(f"{BASE}/network/annotation/{urllib.parse.quote(disease_id)}", source="HPO", cache_ttl=14 * 86400)
    data = resp.json()
    disease = data.get("disease") or {}
    terms = []
    for cat, rows in (data.get("categories") or {}).items():
        for r in rows:
            meta = r.get("metadata") or {}
            terms.append({"id": r.get("id"), "label": r.get("name"), "category": cat, "frequency": meta.get("frequency") or None,
                          "onset": meta.get("onset") or None})
    return Outcome({"disease": {"id": disease.get("id"), "name": disease.get("name"), "mondo": disease.get("mondoId")},
                    "phenotypes": terms}, sources=[source_record("HPO annotations", disease_id, resp)])
