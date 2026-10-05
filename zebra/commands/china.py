"""zebra china: is a disease on China's national rare disease lists?

第一批罕见病目录 (2018, 121 diseases, 国卫医发〔2018〕10号) and 第二批罕见病目录
(2023, 86 diseases, 国卫医政发〔2023〕26号), bundled in zebra/data/china_rare_diseases.json
from the official gov.cn publications (provenance inside the file). Matching is on
the names as published, Chinese or English: exact (after folding case, width,
punctuation and possessive 's), on a part of a name (the bracketed or slash-separated
alternative, e.g. Dravet综合征 / Wilson Disease), containment, similar spelling, then a
shared stretch of name (Duchenne muscular dystrophy ~ Progressive Muscular Dystrophy).
Only exact and part matches count as "on the list"; the rest are offered to verify.
"""

from __future__ import annotations

import argparse
import difflib
import json
import os
import re
import unicodedata
from typing import Any, Dict, List, Optional, Sequence

from zebra.core import Outcome, UsageError
from zebra.http import source_record

DATA_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "china_rare_diseases.json")
LIST_NAMES = {1: "第一批罕见病目录 (2018)", 2: "第二批罕见病目录 (2023)"}
_CJK = re.compile(r"[㐀-鿿]")
_cache: Dict[str, Any] = {}


class ListUnavailable(Exception):
    pass


def load() -> Dict[str, Any]:
    if "data" not in _cache:
        if not os.path.exists(DATA_FILE):
            raise ListUnavailable(f"list not bundled: {DATA_FILE} is missing")
        with open(DATA_FILE, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        if not data.get("diseases"):
            raise ListUnavailable("list file holds no diseases")
        _cache["data"] = data
    return _cache["data"]


def norm(text: str) -> str:
    """Fold width/case/quotes, drop possessive 's and every separator: 'Gaucher’s Disease' -> 'gaucherdisease'."""
    t = unicodedata.normalize("NFKC", text or "").lower()
    t = t.replace("’", "'").replace("‘", "'").replace("＇", "'")
    t = re.sub(r"'s\b", "", t)
    t = t.replace("ⅱ", "ii").replace("ⅲ", "iii").replace("ⅰ", "i")
    return re.sub(r"[\s\-_'\",.;:()\[\]{}/、，。（）·]+", "", t)


def alternates(name: str) -> List[str]:
    """The whole name plus its parts split at '/', '、' and brackets (normalized)."""
    t = unicodedata.normalize("NFKC", name or "")
    parts = [t] + re.split(r"[/()\[\]（）]", t)
    out: List[str] = []
    for p in parts:
        n = norm(p)
        if n and n not in out:
            out.append(n)
    return out


def _min_len(s: str) -> int:
    return 2 if _CJK.search(s) else 4


def match(query: str, entries: Optional[Sequence[Dict[str, Any]]] = None, limit: int = 5) -> List[Dict[str, Any]]:
    """Best list entries for one name; each with match kind exact | part | contains | similar | overlap."""
    q = norm(query)
    if not q:
        return []
    entries = entries if entries is not None else load()["diseases"]
    scored = []
    for e in entries:
        best = None  # (score, kind, on)
        for field in ("name_zh", "name_en"):
            alts = alternates(e.get(field) or "")
            for i, a in enumerate(alts):
                if q == a:
                    cand = (1.0 if i == 0 else 0.95, "exact" if i == 0 else "part", e[field])
                elif len(q) >= _min_len(q) and len(a) >= _min_len(a) and (q in a or a in q):
                    cand = (0.8 * min(len(q), len(a)) / max(len(q), len(a)) + 0.1, "contains", e[field])
                else:
                    sm = difflib.SequenceMatcher(None, q, a)
                    r = sm.ratio()
                    if r >= 0.8:
                        cand = (r * 0.75, "similar", e[field])
                    else:
                        # a shared stretch such as 'musculardystrophy' (Duchenne vs Progressive Muscular Dystrophy)
                        lcs = sm.find_longest_match(0, len(q), 0, len(a)).size
                        need = 4 if _CJK.search(q) else 12
                        cand = (0.3 + 0.2 * lcs / max(len(q), len(a)), "overlap", e[field]) if lcs >= need else None
                if cand and (best is None or cand[0] > best[0]):
                    best = cand
        if best:
            scored.append((best[0], e, best[1], best[2]))
    scored.sort(key=lambda s: (-s[0], s[1]["list"], s[1]["no"]))
    out = []
    for score, e, kind, on in scored[:limit]:
        out.append({"list": e["list"], "list_name": LIST_NAMES.get(e["list"]), "no": e["no"], "name_zh": e["name_zh"],
                    "name_en": e["name_en"], "match": kind, "matched_on": on, "score": round(score, 3)})
    return out


def lookup(names: Sequence[str], limit: int = 3) -> Dict[str, Any]:
    """Match several names of one disease (English, Chinese, synonyms); strongest matches first."""
    seen: Dict[tuple, Dict[str, Any]] = {}
    for n in names:
        if not n:
            continue
        for m in match(n, limit=limit):
            key = (m["list"], m["no"])
            m = dict(m, query=n)
            if key not in seen or m["score"] > seen[key]["score"]:
                seen[key] = m
    rows = sorted(seen.values(), key=lambda m: -m["score"])
    on = [m for m in rows if m["match"] in ("exact", "part")]
    return {"on_list": bool(on), "matches": on or rows[:limit]}


def provenance_brief(data: Dict[str, Any]) -> Dict[str, Any]:
    p = data.get("provenance") or {}
    return {"retrieved_at": p.get("retrieved_at"), "counts": p.get("counts"),
            "source_urls": {k: v.get("source_url") for k, v in (p.get("lists") or {}).items()}}


def source_for(data: Dict[str, Any], hits: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    p = data.get("provenance") or {}
    lists = p.get("lists") or {}
    recs = []
    for li in sorted({h["list"] for h in hits} or {1, 2}):
        meta = lists.get(str(li)) or {}
        nos = ",".join(str(h["no"]) for h in hits if h["list"] == li)
        recs.append(source_record("China national rare disease list", f"{meta.get('name_zh', li)} {('#' + nos) if nos else ''}".strip(),
                                  url=meta.get("source_url"),
                                  note=f"{meta.get('document', '')}; bundled copy retrieved {p.get('retrieved_at')}"))
    return recs


def _run(args: argparse.Namespace) -> Outcome:
    query = " ".join(args.query).strip()
    if not query:
        raise UsageError("give a disease name (Chinese or English)")
    try:
        data = load()
    except (ListUnavailable, OSError, ValueError) as err:
        return Outcome({"query": query, "status": "unavailable", "matches": []}, query={"query": query},
                       warnings=[f"China rare disease list unavailable ({err}); source not consulted"],
                       text=f"{query}: China rare disease list not bundled; source unavailable")
    hits = match(query, data["diseases"], limit=args.limit)
    strong = [h for h in hits if h["match"] in ("exact", "part", "contains") and h["score"] >= 0.5]
    status = "on_list" if any(h["match"] in ("exact", "part") for h in hits) else ("possible" if hits else "not_found")
    shown = strong if status == "on_list" else hits
    counts = (data.get("provenance") or {}).get("counts") or {}
    lines = []
    if status == "not_found":
        lines.append(f"{query}: not on 第一批 (2018, {counts.get('1')}) or 第二批 (2023, {counts.get('2')}) 罕见病目录 "
                     "(name match; the lists name some diseases as groups, e.g. 进行性肌营养不良)")
    else:
        head = "on the national list" if status == "on_list" else "no exact entry; closest names (verify)"
        lines.append(f"{query}: {head}")
        for h in shown:
            lines.append(f"  {h['list_name']} #{h['no']}  {h['name_zh']} / {h['name_en']}  [{h['match']}]")
    result = {"query": query, "status": status, "matches": shown, "provenance": provenance_brief(data)}
    return Outcome(result, sources=source_for(data, shown), text="\n".join(lines), query={"query": query})


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("china", help="China's national rare disease lists (第一批 2018 / 第二批 2023): is a disease on them")
    p.add_argument("query", nargs="+", help="disease name, Chinese or English")
    p.add_argument("--limit", type=int, default=5)
    p.set_defaults(func=_run)
