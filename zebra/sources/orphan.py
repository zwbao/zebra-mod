"""Orphan designations: the EMA (EU) orphan designation register.

Data: EMA's published JSON export of all orphan designations
(https://www.ema.europa.eu/en/documents/report/medicines-output-orphan_designations-json-report_en.json,
~2 MB, refreshed by EMA; cached here for 7 days). A designation means the
medicine was accepted as an orphan for that condition in the EU, not that it is
authorised. FDA's orphan designation database (accessdata.fda.gov OOPD) answers
this client with its abuse-detection page, so US designations are not covered.
"""

from __future__ import annotations

import re
from typing import Any, Dict, Iterable, List

from zebra.core import Outcome
from zebra.http import get_json, source_record

EMA_URL = "https://www.ema.europa.eu/en/documents/report/medicines-output-orphan_designations-json-report_en.json"
_STATUS_ORDER = {"positive": 0}


def match_terms(name: str, synonyms: Iterable[str] = ()) -> List[str]:
    """The disease name plus exact synonyms long enough to match without collisions (no bare abbreviations)."""
    terms: List[str] = []
    for t in [name, *synonyms]:
        t = (t or "").strip()
        if len(t) < 6 or t.isupper() or t.lower() in (x.lower() for x in terms):
            continue
        terms.append(t)
    return terms


def ema_designations(terms: Iterable[str], limit: int = 25) -> Outcome:
    """EU orphan designations whose intended use names one of `terms` (word-bounded, case-insensitive)."""
    terms = [t for t in terms if t]
    if not terms:
        raise ValueError("no disease terms to match")
    resp = get_json(EMA_URL, source="EMA orphan designations", cache_ttl=7 * 86400, timeout=150)
    data = resp.json()
    rows = data.get("data") or []
    pats = [(t, re.compile(r"(?<![\w-])" + re.escape(t) + r"(?![\w-])", re.I)) for t in terms]
    hits: List[Dict[str, Any]] = []
    for r in rows:
        use = r.get("intended_use") or ""
        matched = [t for t, p in pats if p.search(use)]
        if not matched:
            continue
        hits.append({
            "substance": r.get("active_substance"),
            "medicine": r.get("medicine_name") or None,
            "intended_use": use,
            "status": r.get("status"),
            "eu_number": r.get("eu_designation_number"),
            "date": r.get("date_of_designation_or_refusal"),
            "url": r.get("orphan_designation_url"),
            "matched": matched[0],
        })
    hits.sort(key=lambda h: (_STATUS_ORDER.get((h["status"] or "").lower(), 1), _date_key(h["date"])))
    meta = data.get("meta") or {}
    result = {"register": "EMA orphan designations (EU)", "terms": terms, "count": len(hits),
              "dataset_timestamp": meta.get("timestamp"), "records_scanned": len(rows),
              "order": "positive first, then newest designation",
              "note": "matched on the wording of the intended use; a designation worded differently "
                      "(e.g. 'severe myoclonic epilepsy in infancy') is missed",
              "designations": hits[:limit]}
    return Outcome(result, sources=[source_record("EMA orphan designations", ", ".join(terms), resp,
                                                  note=f"dataset {meta.get('timestamp')}; {len(hits)} matched")])


def _date_key(d: Any) -> str:
    # dd/mm/yyyy → reversed for newest-first sorting
    m = re.match(r"^(\d{2})/(\d{2})/(\d{4})$", str(d or ""))
    if not m:
        return "~"
    inv = "".join(chr(ord("9") - ord(c) + ord("0")) for c in m.group(3) + m.group(2) + m.group(1))
    return inv
