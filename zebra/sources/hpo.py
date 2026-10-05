"""HPO terms: the JAX ontology API, with the local release (if fetched) answered first.

API: https://ontology.jax.org/api/hp/ (no key).
"""

from __future__ import annotations

import re
import urllib.parse
from typing import Any, Dict

from zebra import hpo_local
from zebra.core import Outcome
from zebra.http import get_json, source_record

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


def search(text: str, limit: int = 10) -> Outcome:
    resp = get_json(f"{BASE}/hp/search", source="HPO", params={"q": text, "page": 0, "limit": max(1, min(limit, 50))},
                    cache_ttl=14 * 86400)
    data = resp.json()
    hits = [{"id": t.get("id"), "label": t.get("name"), "synonyms": (t.get("synonyms") or [])[:5],
             "descendants": t.get("descendantCount")} for t in data.get("terms", [])][:limit]
    return Outcome({"query": text, "hits": hits, "backend": "HPO API (ontology.jax.org)"},
                   sources=[source_record("HPO search", text, resp)])


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
