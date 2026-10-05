"""Orphanet via the Orphadata API: names, definitions, cross-references, epidemiology,
natural history and genes of rare diseases, by ORPHAcode.

API: https://api.orphadata.com (no key; CC BY 4.0). Datasets used:
  rd-cross-referencing   name, synonyms, definition, OMIM/ICD-10/ICD-11/MONDO/UMLS/MeSH/GARD refs
                         (lang=zh serves Chinese names; that dataset is older, dated 2020-06-01)
  rd-epidemiology        prevalence (type, class, mean value, region, validation, source)
  rd-natural_history     age of onset, type of inheritance
  rd-associated-genes    genes with association type and status
A miss is HTTP 404 with a JSON error body. The name endpoint returns ONE closest
match, not necessarily the name asked for ("Gaucher's Disease" -> Alexander
disease), so callers must check the name before trusting it (`name_matches`).
"""

from __future__ import annotations

import re
import unicodedata
import urllib.parse
from typing import Any, Dict, List, Optional

from zebra.core import Outcome
from zebra.http import get_json, source_record

BASE = "https://api.orphadata.com"
CODE_RE = re.compile(r"^(?:ORPHA|ORPHANET)?[:_ ]?(\d+)$", re.I)


def orpha_code(value: Any) -> str:
    m = CODE_RE.match(str(value).strip())
    if not m:
        raise ValueError(f"{value!r} is not an ORPHAcode (ORPHA:33069)")
    return m.group(1)


def page_url(code: str, lang: str = "en") -> str:
    return f"https://www.orpha.net/{lang}/disease/detail/{code}"


def _results(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    res = (data.get("data") or {}).get("results")
    if res is None:
        return []
    return res if isinstance(res, list) else [res]


def _get(path: str, record: str, label: str, lang: Optional[str] = None):
    params = {"lang": lang} if lang else None
    resp = get_json(f"{BASE}/{path}", source="Orphadata", params=params, cache_ttl=30 * 86400, ok_statuses=(200, 404))
    return resp, (resp.json() if resp.status == 200 else None), source_record(label, record, resp)


# ---------------------------------------------------------------- cross-referencing

_RELATION = re.compile(r"^\s*([A-Z]+)\s*\(")


def parse_disorder(r: Dict[str, Any]) -> Dict[str, Any]:
    refs = []
    for e in r.get("ExternalReference") or []:
        rel_text = e.get("DisorderMappingRelation") or ""
        m = _RELATION.match(rel_text)
        refs.append({
            "source": e.get("Source"),
            "id": str(e.get("Reference")),
            "relation": m.group(1) if m else rel_text or None,  # E exact, NTBT narrower, BTNT broader, ND ...
            "relation_text": rel_text or None,
            "validated": e.get("DisorderMappingValidationStatus"),
            "icd_relation": e.get("DisorderMappingICDRelation"),
        })
    summary = r.get("SummaryInformation") or []
    definition = None
    for s in summary if isinstance(summary, list) else [summary]:
        if isinstance(s, dict) and s.get("Definition"):
            definition = s["Definition"]
            break
    code = str(r.get("ORPHAcode"))
    return {
        "id": f"ORPHA:{code}",
        "orphacode": code,
        "name": r.get("Preferred term"),
        "synonyms": r.get("Synonym") or [],
        "definition": definition,
        "typology": r.get("Typology"),
        "group": r.get("DisorderGroup"),
        "flags": [f.get("Label") or f.get("Value") for f in (r.get("DisorderFlag") or []) if isinstance(f, dict)],
        "refs": refs,
        "dataset_date": r.get("Date"),
        "url": page_url(code),
    }


def exact_refs(disorder: Dict[str, Any], source: str) -> List[str]:
    """Ids in `source` (OMIM, ICD-10, MONDO, ...) that Orphanet maps as exact (E)."""
    return [r["id"] for r in disorder.get("refs", []) if r["source"] == source and r["relation"] == "E"]


def curie_refs(disorder: Dict[str, Any]) -> Dict[str, List[Dict[str, str]]]:
    """refs grouped by source with CURIE spelling (OMIM:607208, MONDO:0011794, ICD-10:G40.4 ...)."""
    prefix = {"OMIM": "OMIM", "MONDO": "MONDO", "UMLS": "UMLS", "MeSH": "MESH", "GARD": "GARD", "MedDRA": "MEDDRA",
              "ICD-10": "ICD10", "ICD-11": "ICD11"}
    out: Dict[str, List[Dict[str, str]]] = {}
    for r in disorder.get("refs", []):
        p = prefix.get(r["source"], r["source"])
        out.setdefault(r["source"], []).append({"id": f"{p}:{r['id']}", "relation": r["relation"]})
    return out


def disorder(code: Any, lang: str = "en") -> Outcome:
    """Name, synonyms, definition and cross-references of one ORPHAcode (lang en, zh, fr, ...)."""
    c = orpha_code(code)
    resp, data, src = _get(f"rd-cross-referencing/orphacodes/{c}", f"ORPHA:{c}", "Orphanet cross-referencing", lang)
    if data is None:
        return Outcome(None, sources=[src], warnings=[f"Orphanet has no ORPHA:{c} ({lang})"])
    rows = _results(data)
    return Outcome(parse_disorder(rows[0]) if rows else None, sources=[src])


def by_omim(omim: Any, lang: str = "en") -> Outcome:
    """Orphanet disorders cross-referenced to an OMIM number (with the mapping relation)."""
    num = str(omim).upper().replace("OMIM:", "").strip()
    if not num.isdigit():
        raise ValueError(f"{omim!r} is not an OMIM number")
    resp, data, src = _get(f"rd-cross-referencing/omims/{num}", f"OMIM:{num}", "Orphanet cross-referencing", lang)
    if data is None:
        return Outcome([], sources=[src])
    rows = []
    for r in _results(data):
        d = parse_disorder(r)
        d["omim_relation"] = next((x["relation"] for x in d["refs"] if x["source"] == "OMIM" and x["id"] == num), None)
        rows.append(d)
    return Outcome(rows, sources=[src])


def by_name(name: str, lang: str = "en") -> Outcome:
    """Orphadata's single closest-name match. Check `name_matches(result, name)` before trusting it."""
    resp, data, src = _get(f"rd-cross-referencing/orphacodes/names/{urllib.parse.quote(name.strip(), safe='')}",
                           name, "Orphanet name search", lang)
    if data is None:
        return Outcome(None, sources=[src])
    rows = _results(data)
    return Outcome(parse_disorder(rows[0]) if rows else None, sources=[src])


def norm_name(text: str) -> str:
    t = unicodedata.normalize("NFKC", text or "").lower()
    t = t.replace("’", "'").replace("‘", "'")
    t = re.sub(r"'s\b", "", t)
    t = re.sub(r"[\s\-_'\",.;:()\[\]{}/]+", " ", t)
    return t.strip()


def name_matches(disorder_row: Optional[Dict[str, Any]], name: str) -> bool:
    if not disorder_row:
        return False
    q = norm_name(name)
    names = [disorder_row.get("name") or ""] + list(disorder_row.get("synonyms") or [])
    return any(norm_name(n) == q for n in names)


# ---------------------------------------------------------------- epidemiology / natural history / genes


def parse_epidemiology(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    out = []
    for r in _results(data):
        for p in r.get("Prevalence") or []:
            out.append({
                "type": p.get("PrevalenceType"),
                "class": p.get("PrevalenceClass"),
                "value": p.get("ValMoy") if p.get("ValMoy") not in (None, "0.0", "0") else None,
                "qualification": p.get("PrevalenceQualification"),
                "region": p.get("PrevalenceGeographic"),
                "validation": p.get("PrevalenceValidationStatus"),
                "source": p.get("Source"),
            })
    return out


def epidemiology(code: Any) -> Outcome:
    c = orpha_code(code)
    resp, data, src = _get(f"rd-epidemiology/orphacodes/{c}", f"ORPHA:{c}", "Orphanet epidemiology")
    return Outcome(parse_epidemiology(data) if data else [], sources=[src])


def parse_natural_history(data: Dict[str, Any]) -> Dict[str, Any]:
    rows = _results(data)
    if not rows:
        return {"onset": [], "inheritance": []}
    r = rows[0]
    return {"onset": r.get("AverageAgeOfOnset") or [], "inheritance": r.get("TypeOfInheritance") or []}


def natural_history(code: Any) -> Outcome:
    c = orpha_code(code)
    resp, data, src = _get(f"rd-natural_history/orphacodes/{c}", f"ORPHA:{c}", "Orphanet natural history")
    return Outcome(parse_natural_history(data) if data else {"onset": [], "inheritance": []}, sources=[src])


def parse_genes(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    out = []
    for r in _results(data):
        for a in r.get("DisorderGeneAssociation") or []:
            g = a.get("Gene") or {}
            refs = {e.get("Source"): e.get("Reference") for e in (g.get("ExternalReference") or [])}
            pmids = re.findall(r"(\d+)\[PMID\]", a.get("SourceOfValidation") or "")
            out.append({
                "symbol": g.get("Symbol") or refs.get("Genatlas"),
                "name": g.get("Name") or g.get("name"),
                "hgnc": f"HGNC:{refs['HGNC']}" if refs.get("HGNC") else None,
                "locus": ",".join(x.get("GeneLocus", "") for x in (g.get("Locus") or []) if x.get("GeneLocus")) or None,
                "association": a.get("DisorderGeneAssociationType"),
                "status": a.get("DisorderGeneAssociationStatus"),
                "pmids": pmids[:8],
            })
    return out


def genes(code: Any) -> Outcome:
    c = orpha_code(code)
    resp, data, src = _get(f"rd-associated-genes/orphacodes/{c}", f"ORPHA:{c}", "Orphanet genes")
    return Outcome(parse_genes(data) if data else [], sources=[src])
