"""MedlinePlus Genetics: a plain-language description of a genetic condition (F39).

API: the JSON data files MedlinePlus publishes for its Genetics summaries,
`https://medlineplus.gov/download/genetics/condition/<slug>.json` (no key;
documented at https://medlineplus.gov/about/data-files-api). The slug is the
condition's own page name — `spinal-muscular-atrophy`,
`duchenne-and-becker-muscular-dystrophy`. MedlinePlus publishes no index of
slugs that this client could retrieve (`/download/genetics/condition/` and
`/download/genetics/index.json` both answer the site's HTML 404 page, checked
2026-10-06), and the site has no search API, so a page is found by trying slugs
built from the names a disease card already holds — its Orphanet and Mondo
preferred terms and synonyms.

A miss is a clean HTTP 404 with an HTML body (checked for `dravet-syndrome` and
`no-such-condition-xyz`), so it is returned as "not found", never retried.

A name-built slug can in principle land on a different condition, so a hit is
cross-checked against the OMIM numbers the file lists in its `db-key-list`
(`spinal-muscular-atrophy.json` lists OMIM 253300, 253400, 253550, 271150).
`id_confirmed` says whether such an overlap was found; when it is False the
match rests on the name alone and the caller is expected to say so.

Written for lay readers: the text is MedlinePlus's own wording, kept as plain
text with its paragraph breaks, and is a description of the condition, never
advice about a particular patient.
"""

from __future__ import annotations

import html
import re
import urllib.parse
from typing import Any, Dict, Iterable, List, Sequence

from zebra.core import Outcome
from zebra.http import get_json
from zebra.sources import record as source_record
from zebra.sources import validated_json

BASE = "https://medlineplus.gov/download/genetics/condition"
PAGE = "https://medlineplus.gov/genetics/condition"
MAX_SLUGS = 6  # one request per candidate name; the card holds a dozen names


def slug(name: str) -> str:
    """A MedlinePlus Genetics page name from a disease name: lower case, words joined by `-`."""
    t = html.unescape(name or "").lower().replace("'", "").replace("’", "")
    t = re.sub(r"[^a-z0-9]+", "-", t).strip("-")
    return t


def slug_candidates(names: Iterable[str]) -> List[str]:
    """Distinct slugs worth trying, in the order the names were given."""
    out: List[str] = []
    for n in names:
        s = slug(n)
        if s and 3 < len(s) <= 100 and s not in out:
            out.append(s)
    return out[:MAX_SLUGS]


def _text_blocks(data: Dict[str, Any]) -> List[Dict[str, str]]:
    blocks = []
    for item in data.get("text-list") or []:
        t = (item or {}).get("text") or {}
        raw = t.get("html") or ""
        paras = [" ".join(html.unescape(re.sub(r"<[^>]+>", " ", p)).split())
                 for p in re.split(r"</p\s*>", raw)]
        body = "\n\n".join(p for p in paras if p)
        if body:
            blocks.append({"role": t.get("text-role") or "description", "text": body})
    return blocks


def parse_condition(data: Dict[str, Any], asked_slug: str, omim_ids: Sequence[str] = ()) -> Dict[str, Any]:
    keys: Dict[str, List[str]] = {}
    for item in data.get("db-key-list") or []:
        k = (item or {}).get("db-key") or {}
        if k.get("db") and k.get("key"):
            keys.setdefault(k["db"], []).append(str(k["key"]))
    want = {str(o).upper().replace("OMIM:", "").strip() for o in omim_ids}
    overlap = sorted(want & set(keys.get("OMIM") or []))
    return {
        "slug": asked_slug,
        "name": data.get("name"),
        "url": data.get("ghr_page") or f"{PAGE}/{asked_slug}",
        "text": _text_blocks(data),
        "synonyms": [(s or {}).get("synonym") for s in data.get("synonym-list") or []][:10],
        "genes": [((g or {}).get("related-gene") or {}).get("gene-symbol")
                  for g in data.get("related-gene-list") or []],
        "inheritance": [((i or {}).get("inheritance-pattern") or {}).get("memo")
                        for i in data.get("inheritance-pattern-list") or []],
        "xrefs": keys,
        "omim_overlap": overlap,
        "id_confirmed": bool(overlap),
        "reviewed": data.get("reviewed"),
        "published": data.get("published"),
        "register": "MedlinePlus Genetics (US National Library of Medicine), written for patients and families",
    }


def condition(names: Iterable[str], omim_ids: Sequence[str] = ()) -> Outcome:
    """The MedlinePlus Genetics summary for the first of `names` that has a page."""
    cands = slug_candidates(names)
    if not cands:
        return Outcome(None, warnings=["MedlinePlus Genetics not searched: no disease name to build a page name from"])
    sources: List[Dict[str, Any]] = []
    tried: List[str] = []
    for s in cands:
        resp = get_json(f"{BASE}/{urllib.parse.quote(s, safe='')}.json", source="MedlinePlus Genetics",
                        accept="application/json, */*", cache_ttl=30 * 86400, ok_statuses=(200, 404))
        sources.append(source_record("MedlinePlus Genetics", s, resp, url=f"{PAGE}/{s}"))
        tried.append(s)
        if resp.status != 200:
            continue
        data = validated_json(resp, "MedlinePlus Genetics")
        if not isinstance(data, dict) or not data.get("name"):
            continue
        res = parse_condition(data, s, omim_ids)
        res["tried"] = tried
        warnings = []
        if not res["id_confirmed"]:
            warnings.append(f"MedlinePlus Genetics page '{res['name']}' was found by name ('{s}'), and its OMIM "
                            f"numbers {', '.join(res['xrefs'].get('OMIM') or []) or '(none listed)'} do not overlap "
                            "the ids on this card: confirm it is the same condition before quoting it")
        return Outcome(res, sources=sources, warnings=warnings)
    return Outcome(None, sources=sources,
                   warnings=[f"no MedlinePlus Genetics page found for this disease: tried the page names "
                             f"{', '.join(tried)} (MedlinePlus publishes no slug index to search)"])
