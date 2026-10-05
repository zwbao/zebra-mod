"""EBI OLS4 (ontology terms, MONDO cross-references, obsolete -> replacement) and
EBI OxO (cross-ontology mappings), used to reconcile disease ids across sources.

OLS4: https://www.ebi.ac.uk/ols4/api/ontologies/{ontology}/terms?obo_id=MONDO:0100135
OxO:  https://www.ebi.ac.uk/spot/oxo/api/mappings?fromId=OMIM:607208
No key. OLS4 cannot look a MONDO class up by one of its xrefs (searching
database_cross_reference returns nothing); for OMIM/Orphanet -> MONDO use the
Monarch /mappings endpoint or OxO.
"""

from __future__ import annotations

import urllib.parse
from typing import Any, Dict, List, Optional, Sequence

from zebra.core import Outcome
from zebra.http import get_json, source_record

OLS = "https://www.ebi.ac.uk/ols4/api"
OXO = "https://www.ebi.ac.uk/spot/oxo/api"


def _short_to_curie(short: Optional[str]) -> Optional[str]:
    if not short:
        return None
    if short.startswith("http"):
        short = short.rstrip("/").rsplit("/", 1)[-1]
    return short.replace("_", ":", 1) if ":" not in short else short


def parse_term(data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    terms = (data.get("_embedded") or {}).get("terms") or []
    defining = [t for t in terms if t.get("is_defining_ontology")] or terms
    if not defining:
        return None
    t = defining[0]
    ann = t.get("annotation") or {}
    xrefs = []
    for x in t.get("obo_xref") or []:
        if x.get("database") and x.get("id"):
            xrefs.append({"id": f"{x['database']}:{x['id']}", "relation": (x.get("description") or "").replace("MONDO:", "") or None})
    if not xrefs:
        xrefs = [{"id": x, "relation": None} for x in ann.get("database_cross_reference") or []]
    for x in xrefs:  # zebra spells Orphanet ids ORPHA:
        if x["id"].startswith("Orphanet:"):
            x["id"] = "ORPHA:" + x["id"].split(":", 1)[1]
    desc = t.get("description") or []
    return {
        "id": t.get("obo_id"),
        "label": t.get("label"),
        "definition": next((d for d in desc if not d.startswith(("This is a distinct", "OBSOLETE."))), desc[0] if desc else None),
        "notes": [d for d in desc if d.startswith("This is a distinct")],
        "synonyms": (t.get("synonyms") or [])[:15],
        "obsolete": bool(t.get("is_obsolete")),
        "replaced_by": _short_to_curie(t.get("term_replaced_by")) or _short_to_curie((ann.get("term replaced by") or [None])[0]),
        "xrefs": xrefs,
        "exact_xrefs": [x["id"] for x in xrefs if x["relation"] in ("equivalentTo",)],
        "subsets": t.get("in_subset") or [],
        "url": f"https://www.ebi.ac.uk/ols4/ontologies/{t.get('ontology_name')}/classes?obo_id={urllib.parse.quote(t.get('obo_id') or '')}",
    }


def term(curie: str, ontology: Optional[str] = None) -> Outcome:
    """One ontology class (default ontology = the CURIE prefix, e.g. mondo, hp, ordo)."""
    curie = curie.strip()
    onto = (ontology or curie.split(":", 1)[0]).lower()
    resp = get_json(f"{OLS}/ontologies/{onto}/terms", source="OLS", params={"obo_id": curie}, cache_ttl=30 * 86400,
                    ok_statuses=(200, 404))
    if resp.status == 404:
        return Outcome(None, sources=[source_record("OLS", curie, resp)], warnings=[f"OLS has no {curie} in {onto}"])
    return Outcome(parse_term(resp.json()), sources=[source_record("OLS", curie, resp)])


def resolve_obsolete(curies: Sequence[str]) -> Outcome:
    """For each MONDO id: its current replacement if obsolete, else itself (None when OLS does not know it)."""
    out: Dict[str, Optional[str]] = {}
    srcs: List[Dict[str, Any]] = []
    warnings: List[str] = []
    for c in dict.fromkeys(curies):
        got = term(c)
        srcs.extend(got.sources)
        warnings.extend(got.warnings)
        t = got.result
        if not t:
            out[c] = None
        elif t["obsolete"]:
            out[c] = t["replaced_by"]
        else:
            out[c] = c
    return Outcome(out, sources=srcs, warnings=warnings)


def parse_oxo(data: Dict[str, Any], from_id: str, exact_only: bool = True) -> List[Dict[str, Any]]:
    rows = []
    seen = set()
    for m in (data.get("_embedded") or {}).get("mappings") or []:
        scope = m.get("scope")
        if exact_only and scope != "EXACT":
            continue
        a = (m.get("fromTerm") or {}).get("curie")
        b = (m.get("toTerm") or {}).get("curie")
        other, label = (b, (m.get("toTerm") or {}).get("label")) if a and a.upper() == from_id.upper() else (a, (m.get("fromTerm") or {}).get("label"))
        if not other:
            continue
        other = other.replace("ORPHANET:", "ORPHA:").replace("Orphanet:", "ORPHA:")
        key = (other, scope)
        if key in seen:
            continue
        seen.add(key)
        rows.append({"id": other, "label": label, "scope": scope,
                     "mapping_set": ((m.get("datasource") or {}).get("prefix") or "").rsplit("/", 1)[-1]})
    return rows


def oxo_mappings(curie: str, exact_only: bool = True) -> Outcome:
    """Cross-ontology mappings of one id from EBI OxO (EXACT scope only by default)."""
    query = curie.replace("ORPHA:", "Orphanet:") if curie.upper().startswith("ORPHA:") else curie
    resp = get_json(f"{OXO}/mappings", source="OxO", params={"fromId": query, "size": 200}, cache_ttl=30 * 86400)
    return Outcome(parse_oxo(resp.json(), query, exact_only), sources=[source_record("EBI OxO", curie, resp)])
