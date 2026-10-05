"""zebra china: is a disease on China's national rare disease lists?

第一批罕见病目录 (2018, 121 diseases, 国卫医发〔2018〕10号) and 第二批罕见病目录
(2023, 86 diseases, 国卫医政发〔2023〕26号), bundled in zebra/data/china_rare_diseases.json
from the official gov.cn publications (provenance inside the file).

Names families actually use are matched through an alias layer (P1f), because the
published names are not the ones people say: `SMA`, `DMD`, `小胖威利`, `瓷娃娃`,
`渐冻症`, `快乐木偶综合征`, `德拉韦综合征` and `普拉德-威利综合征` all answered
"not on the list" for diseases that are on the 2018 list, and `瑞特综合征`'s closest
match was 白塞病. The alias layer is zebra/data/china_disease_aliases.json: per list
entry, the official Chinese and English names and their bracketed parts, Orphanet's
English and Chinese preferred terms and synonyms for the ORPHAcode the entry links
to, and a small hand-entered table of folk names. Every alias carries a `kind` and a
`source`; a folk name is marked `folk_name` and names the official entry it maps to,
because it comes from no dataset. The file's `provenance` records the Orphanet
endpoints, their retrieval date and the dataset date, the linking rules in the order
they were tried, and the 58 list entries that could not be linked to one ORPHAcode.

Matching, in order: exact alias (after folding case, width, punctuation and
possessive 's), an alias part, containment, Chinese character-bigram overlap
(Dice, threshold BIGRAM_MIN), similar spelling, then a shared stretch of name.
Only exact, part and alias matches count as "on the list"; everything else is
offered to verify, and every answer says which alias matched, of what kind and from
what source. "not on the list" is only ever said when nothing matched at all.
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
from zebra.sources import record as source_record

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
DATA_FILE = os.path.join(DATA_DIR, "china_rare_diseases.json")
ALIAS_FILE = os.path.join(DATA_DIR, "china_disease_aliases.json")
ZH_NAMES_FILE = os.path.join(DATA_DIR, "orphanet_zh_names.json")
LIST_NAMES = {1: "第一批罕见病目录 (2018)", 2: "第二批罕见病目录 (2023)"}
_CJK = re.compile(r"[㐀-鿿]")
# Dice coefficient on character bigrams. 0.6 keeps 杜氏肌营养不良 ~ 杜氏肌营养不良症
# (0.923), 脊髓性肌萎缩症 ~ 近端脊髓性肌萎缩 (0.769) and 成骨不全症 ~ 成骨不全
# (0.857), and rejects the family of one-letter syndrome names that otherwise come
# back for any `X综合征` query (德拉韦综合征 ~ W综合征 scores 0.5).
BIGRAM_MIN = 0.6
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


def load_aliases() -> Dict[str, Any]:
    """The alias layer. Missing or unreadable is not fatal: matching falls back to the published names."""
    if "aliases" not in _cache:
        data: Dict[str, Any] = {"entries": [], "provenance": {"note": f"{ALIAS_FILE} not bundled"}}
        try:
            with open(ALIAS_FILE, "r", encoding="utf-8") as fh:
                loaded = json.load(fh)
            if loaded.get("entries"):
                data = loaded
        except (OSError, ValueError):
            pass
        index: List[Dict[str, Any]] = []
        for e in data["entries"]:
            for a in e.get("aliases") or []:
                text = a.get("text") or ""
                if not text:
                    continue
                index.append({"norm": norm(text), "bigrams": bigrams(text), "alias": a, "entry": e})
        data["_index"] = index
        _cache["aliases"] = data
    return _cache["aliases"]


def load_zh_names() -> Dict[str, Any]:
    """Orphanet's Chinese preferred term per ORPHAcode, for resolving a Chinese name to one disease."""
    if "zh_names" not in _cache:
        data: Dict[str, Any] = {"names": {}, "provenance": {"note": f"{ZH_NAMES_FILE} not bundled"}}
        try:
            with open(ZH_NAMES_FILE, "r", encoding="utf-8") as fh:
                loaded = json.load(fh)
            if loaded.get("names"):
                data = loaded
        except (OSError, ValueError):
            pass
        data["_index"] = [(code, name, bigrams(name)) for code, name in data["names"].items()]
        _cache["zh_names"] = data
    return _cache["zh_names"]


def bigrams(text: str) -> frozenset:
    """Character bigrams of a normalised name; a one-character name is its own single gram."""
    t = norm(text)
    if len(t) < 2:
        return frozenset([t]) if t else frozenset()
    return frozenset(t[i:i + 2] for i in range(len(t) - 1))


def dice(a: frozenset, b: frozenset) -> float:
    """Sørensen-Dice overlap of two bigram sets, 0.0 when either is empty."""
    if not a or not b:
        return 0.0
    return 2.0 * len(a & b) / (len(a) + len(b))


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


def alias_hits(query: str, limit: int = 5) -> List[Dict[str, Any]]:
    """List entries whose alias layer matches `query`; strongest first, each saying what matched.

    Kinds, strongest first: `alias_exact` (a normalised alias equals the query),
    `alias_part` (an alias part equals it), `alias_contains`, `alias_bigram`
    (Chinese character-bigram Dice >= BIGRAM_MIN).
    """
    q = norm(query)
    if not q:
        return []
    # Character-bigram overlap is for Chinese input only. On Latin names it
    # matches on a shared word: "NGLY1 deficiency" scores 0.72 against
    # "T2 deficiency", which is not a near-miss of anything. Latin spelling
    # variants are handled by `match()` (containment, SequenceMatcher).
    cjk_query = bool(_CJK.search(query))
    qb = bigrams(query)
    best: Dict[tuple, Dict[str, Any]] = {}
    for row in load_aliases()["_index"]:
        a, e = row["alias"], row["entry"]
        kind, score = None, 0.0
        if row["norm"] == q:
            kind, score = "alias_exact", 1.0
        elif len(q) >= _min_len(q) and len(row["norm"]) >= _min_len(row["norm"]) and (
                q in row["norm"] or row["norm"] in q):
            kind = "alias_contains"
            score = 0.75 * min(len(q), len(row["norm"])) / max(len(q), len(row["norm"]))
        elif cjk_query and _CJK.search(a.get("text") or ""):
            d = dice(qb, row["bigrams"])
            if d >= BIGRAM_MIN:
                kind, score = "alias_bigram", 0.7 * d
        if kind is None:
            continue
        key = (e["list"], e["no"])
        cand = {
            "list": e["list"], "list_name": LIST_NAMES.get(e["list"]), "no": e["no"],
            "name_zh": e["name_zh"], "name_en": e["name_en"], "orpha": e.get("orpha"),
            "orpha_name_en": e.get("orpha_name_en"), "orpha_name_zh": e.get("orpha_name_zh"),
            "orpha_disorder_group": e.get("orpha_disorder_group"),
            "match": kind, "matched_on": a["text"], "alias_kind": a["kind"],
            "alias_source": a.get("source"), "maps_to": a.get("maps_to"),
            "score": round(score, 3),
        }
        if key not in best or cand["score"] > best[key]["score"]:
            best[key] = cand
    rows = sorted(best.values(), key=lambda m: (-m["score"], m["list"], m["no"]))
    return rows[:limit]


def zh_name_hits(query: str, limit: int = 5) -> List[Dict[str, Any]]:
    """ORPHAcodes whose Orphanet Chinese preferred term matches `query` (exact, then bigram Dice)."""
    q = norm(query)
    if not q or not _CJK.search(query):
        return []
    qb = bigrams(query)
    scored = []
    for code, name, nb in load_zh_names()["_index"]:
        if norm(name) == q:
            scored.append((1.0, "exact", code, name))
            continue
        d = dice(qb, nb)
        if d >= BIGRAM_MIN:
            scored.append((d, "bigram", code, name))
    scored.sort(key=lambda s: (-s[0], int(s[2])))
    return [{"orpha": f"ORPHA:{c}", "name_zh": n, "similarity": round(s, 3), "match": k,
             "source": "Orphanet Chinese preferred term (dataset "
                       f"{(load_zh_names().get('provenance') or {}).get('dataset_date', '?')})"}
            for s, k, c, n in scored[:limit]]


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


ON_LIST_KINDS = ("exact", "part", "alias_exact", "alias_part")


def lookup(names: Sequence[str], limit: int = 3) -> Dict[str, Any]:
    """Match several names of one disease (English, Chinese, synonyms, folk names); strongest first."""
    seen: Dict[tuple, Dict[str, Any]] = {}
    for n in names:
        if not n:
            continue
        for m in list(alias_hits(n, limit=limit)) + list(match(n, limit=limit)):
            key = (m["list"], m["no"])
            m = dict(m, query=n)
            if key not in seen or m["score"] > seen[key]["score"]:
                seen[key] = m
    rows = sorted(seen.values(), key=lambda m: -m["score"])
    on = [m for m in rows if m["match"] in ON_LIST_KINDS]
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
    aliases = alias_hits(query, limit=args.limit)
    by_key = {(h["list"], h["no"]): h for h in aliases}
    for h in match(query, data["diseases"], limit=args.limit):
        key = (h["list"], h["no"])
        if key not in by_key or h["score"] > by_key[key]["score"]:
            by_key[key] = h
    hits = sorted(by_key.values(), key=lambda m: (-m["score"], m["list"], m["no"]))[:args.limit]
    strong = [h for h in hits if h["match"] in ON_LIST_KINDS + ("contains", "alias_contains", "alias_bigram")
              and h["score"] >= 0.4]
    status = "on_list" if any(h["match"] in ON_LIST_KINDS for h in hits) else ("possible" if hits else "not_found")
    shown = strong if status == "on_list" else hits
    counts = (data.get("provenance") or {}).get("counts") or {}
    alias_prov = (load_aliases().get("provenance") or {})
    lines = []
    if status == "not_found":
        # P1f: never say "not on the list" for a name that simply was not matched.
        lines.append(f"{query}: no name on 第一批 (2018, {counts.get('1')}) or 第二批 (2023, {counts.get('2')}) "
                     "罕见病目录, and no alias of one, matched this text — which is not the same as the disease "
                     "being absent from the lists. Resolve the name first (`zebra disease <name>`) and query the "
                     "official Chinese or English name, or an ORPHA/OMIM id's name.")
    else:
        head = "on the national list" if status == "on_list" else "no entry matched by name; closest (verify)"
        lines.append(f"{query}: {head}")
        for h in shown:
            lines.append(f"  {h['list_name']} #{h['no']}  {h['name_zh']} / {h['name_en']}  [{h['match']}"
                         + (f", score {h['score']}" if h["match"] != "exact" else "") + "]")
            if h.get("alias_kind"):
                lines.append(f"      matched '{h['matched_on']}' — {h['alias_kind']}"
                             + (f", maps to the official name {h['maps_to']}" if h.get("maps_to") else "")
                             + f" ({h.get('alias_source')})")
            elif h.get("matched_on"):
                lines.append(f"      matched the published name '{h['matched_on']}' "
                             f"({LIST_NAMES.get(h['list'])} #{h['no']})")
    zh = zh_name_hits(query, limit=3)
    if zh:
        lines.append("Orphanet Chinese names close to this text (use with `zebra disease <ORPHA id>`):")
        for z in zh:
            lines.append(f"  {z['orpha']}  {z['name_zh']}  [{z['match']}, similarity {z['similarity']}]")
    result = {"query": query, "status": status, "matches": shown, "orphanet_zh_candidates": zh,
              "provenance": provenance_brief(data),
              "alias_provenance": {k: alias_prov.get(k) for k in ("retrieved_at", "sources", "counts", "rules")
                                   if k in alias_prov}}
    return Outcome(result, sources=source_for(data, shown) + alias_sources(shown, zh),
                   text="\n".join(lines), query={"query": query})


def alias_sources(hits: Sequence[Dict[str, Any]], zh: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Provenance rows for the alias layer and the Chinese name index, when either was used."""
    recs: List[Dict[str, Any]] = []
    prov = load_aliases().get("provenance") or {}
    if any(h.get("alias_kind") for h in hits):
        recs.append(source_record("China rare disease list alias layer",
                                  f"china_disease_aliases.json ({(prov.get('counts') or {}).get('aliases')} aliases)",
                                  url=((prov.get("sources") or {}).get("orphanet_zh") or {}).get("url"),
                                  note="bundled; every alias carries its kind and source, folk names are marked "
                                       f"folk_name (built {prov.get('retrieved_at')})"))
    if zh:
        zprov = load_zh_names().get("provenance") or {}
        recs.append(source_record("Orphanet Chinese preferred terms",
                                  f"{zprov.get('names_kept')} ORPHAcodes", url=zprov.get("source_url"),
                                  note=f"bundled copy retrieved {zprov.get('retrieved_at')}, "
                                       f"Orphanet zh dataset {zprov.get('dataset_date')}, "
                                       f"{zprov.get('licence')}"))
    return recs


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("china", help="China's national rare disease lists (第一批 2018 / 第二批 2023): is a disease on them")
    p.add_argument("query", nargs="+", help="disease name, Chinese or English")
    p.add_argument("--limit", type=int, default=5)
    p.set_defaults(func=_run)
