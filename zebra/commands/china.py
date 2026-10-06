"""zebra china: is a disease on China's national rare disease lists? And `zebra china hospitals`.

第一批罕见病目录 (2018, 121 diseases, 国卫医发〔2018〕10号) and 第二批罕见病目录
(2023, 86 diseases, 国卫医政发〔2023〕26号), bundled in zebra/data/china_rare_diseases.json
from the official gov.cn publications (provenance inside the file).

Names families actually use are matched through an alias layer
(zebra/data/china_disease_aliases.json, built by tools/china_list/build_china_aliases.py):
per list entry, the official Chinese and English names and their bracketed parts,
Orphanet's English and Chinese preferred terms and synonyms for the ORPHAcode the
entry links to, Orphanet subtypes and group members, and a small hand-entered table
of folk names. Every alias carries a `kind` and a `source`.

What a match means (statuses, strongest first):
  on_list    the text is the entry's published name, one of its synonyms, an Orphanet
             term of the linked disease, a subtype inside the entry's qualifier
             (庞贝病 / Pompe disease inside 糖原累积病（I型、Ⅱ型）), or a member of a
             group entry by Orphanet's own classification (Duchenne muscular dystrophy
             inside 进行性肌营养不良).
  qualified  the entry names only a subtype of what was asked: 帕金森病 against
             帕金森病（青年型、早发型）, 地中海贫血 against 地中海贫血（重型）. The answer
             says which subtype the list covers; it is never "on the list" (review E-1).
  possible   an acronym (CAD, HSP, PV, GSD, SMA: ambiguous across diseases), a disease
             Orphanet only files under a group entry (group_filed: Gastroschisis under
             Short bowel syndrome), a
             containment, a Chinese character-bigram overlap (Dice >= BIGRAM_MIN), a
             similar spelling or a shared stretch of name. Offered to verify.
  not_found  nothing matched; never worded as "not on the list".
The published bracket qualifiers themselves (重型, I型、Ⅱ型, 青年型、早发型) and
word-internal glosses (冷吡啉) are not names and match nothing on their own.

`zebra china hospitals [--province 浙江]` lists the 全国罕见病诊疗协作网 hospitals
bundled in zebra/data/china_rare_network_hospitals.json (the 2019 NHC notice and any
later official update, provenance inside the file).
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
HOSPITALS_FILE = os.path.join(DATA_DIR, "china_rare_network_hospitals.json")
LIST_NAMES = {1: "第一批罕见病目录 (2018)", 2: "第二批罕见病目录 (2023)"}
_CJK = re.compile(r"[㐀-鿿]")
# Dice coefficient on character bigrams. 0.6 keeps 杜氏肌营养不良 ~ 杜氏肌营养不良症
# (0.923), 脊髓性肌萎缩症 ~ 近端脊髓性肌萎缩 (0.769) and 成骨不全症 ~ 成骨不全
# (0.857), and rejects the family of one-letter syndrome names that otherwise come
# back for any `X综合征` query (德拉韦综合征 ~ W综合征 scores 0.5).
BIGRAM_MIN = 0.6
_cache: Dict[str, Any] = {}

# match kinds, by what they let the answer say
ON_LIST_KINDS = ("exact", "part", "alias_exact", "alias_part", "subtype", "group_member", "orpha_exact")
QUALIFIED_KINDS = ("wider",)
STATUS_TEXT = {
    "on_list": "on the national list",
    "qualified": "only a subtype is on the national list (see the entry's qualifier)",
    "possible": "no entry matched this name exactly; closest entries to verify",
}
HOSPITAL_ROLES = {"national_lead": "国家级牵头医院", "provincial_lead": "省级牵头医院", "member": "成员医院"}


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
        parts: Dict[tuple, Dict[str, Any]] = {}
        by_key: Dict[tuple, Dict[str, Any]] = {}
        for e in data["entries"]:
            by_key[(e["list"], e["no"])] = e
            for a in e.get("aliases") or []:
                text = a.get("text") or ""
                if not text:
                    continue
                index.append({"norm": norm(text), "bigrams": bigrams(text), "alias": a, "entry": e})
                if str(a.get("kind", "")).startswith("official_part"):
                    parts[(e["list"], e["no"], norm(text))] = a
        data["_index"] = index
        data["_parts"] = parts
        data["_by_key"] = by_key
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
    t = re.sub(r"\b([a-z]{4,})s(?=\s+(?:disease|syndrome)\b)", r"\1", t)  # parkinsons disease ~ parkinson's disease
    t = re.sub(r"([a-z])ae", r"\1e", t)  # thalassaemia / haemophilia / anaemia ~ American spelling
    t = re.sub(r"(?<=[a-z0-9])\s*(?:综合征|综合症)$", "syndrome", t)  # Pompe病 ~ Pompe disease
    t = re.sub(r"(?<=[a-z0-9])\s*(?:病|症)$", "disease", t)
    t = t.replace("ⅱ", "ii").replace("ⅲ", "iii").replace("ⅰ", "i")
    # "Type I、II" must not fold into "typeiii" (GSD type III is not on the list): a list separator
    # between two numerals is kept as "|" (same rule as tools/china_list/build_china_aliases.py)
    t = re.sub(r"(?<=\b[ivx0-9])\s*[、,，/;]\s*(?=[ivx0-9]+\b)|(?<=[^a-z][ivx])\s*[、,，/;]\s*(?=[ivx0-9])", "|", t)
    return re.sub(r"[\s\-_'\",.;:()\[\]{}/、，。（）·]+", "", t)


def alternates(name: str) -> List[str]:
    """The whole name plus its parts split at '/' and brackets (normalized)."""
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


_GENERIC_WORDS = frozenset(("deficiency", "disease", "syndrome", "disorder", "type", "related", "associated",
                            "spectrum", "complex"))


def _acronym_token(t: str) -> bool:
    if not t or len(t) > 8 or re.search(r"[^\x00-\x7f]", t):
        return False
    caps = sum(1 for c in t if c.isupper())
    return caps >= 2 and caps >= len(re.sub(r"[^A-Za-z]", "", t)) - 1


def is_acronym(text: str) -> bool:
    """An acronym, alone or with only generic words (CAD, HSP, aHUS, X-ALD, MAD deficiency, G6P deficiency, GSD type 2)."""
    t = (text or "").strip()
    if not t:
        return False
    words = t.split()
    if len(words) == 1:
        return _acronym_token(t)
    return _acronym_token(words[0]) and all(w.lower() in _GENERIC_WORDS or re.fullmatch(r"[0-9ivx]+[a-z]?", w.lower())
                                             for w in words[1:])


def _entry_fields(e: Dict[str, Any]) -> Dict[str, Any]:
    out = {"list": e["list"], "list_name": LIST_NAMES.get(e["list"]), "no": e["no"],
           "name_zh": e["name_zh"], "name_en": e["name_en"]}
    al = load_aliases()["_by_key"].get((e["list"], e["no"])) or e
    for k in ("orpha", "orpha_name_en", "orpha_name_zh", "orpha_disorder_group", "orpha_match", "qualifier",
              "orpha_scope"):
        if al.get(k) is not None:
            out[k] = al[k]
    return out


def _alias_kind(a: Dict[str, Any]) -> str:
    """What an exact match on this alias lets the answer say."""
    if a.get("acronym") or is_acronym(a.get("text") or ""):
        return "acronym"
    if a.get("scope") == "wider":
        return "wider"
    kind = str(a.get("kind") or "")
    if kind.startswith("orphanet_subtype") or (kind == "folk_name" and a.get("qualifier_part")):
        return "subtype"
    if kind.startswith("orphanet_member"):
        # P0-2: only a member whose own name carries the group's name is a form of the listed
        # disease; one Orphanet merely files under the group (Gastroschisis under Short bowel
        # syndrome, Whipple disease under Primary immunodeficiency) is offered, never claimed
        return "group_member" if a.get("member_named") else "group_filed"
    if kind.startswith("official_subtype"):
        return "subtype"
    return "alias_exact"


_EXACT_SCORE = {"alias_exact": 1.0, "subtype": 0.97, "group_member": 0.93, "wider": 0.9, "acronym": 0.6,
                "group_filed": 0.55}


def alias_hits(query: str, limit: int = 5) -> List[Dict[str, Any]]:
    """List entries whose alias layer matches `query`; strongest first, each saying what matched.

    Kinds: `alias_exact`, `subtype`, `group_member` (on the list), `wider` (the
    entry covers only a subtype: qualified), `acronym`, `alias_contains`,
    `alias_bigram` (Chinese character-bigram Dice >= BIGRAM_MIN) — offered only.
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
            kind = _alias_kind(a)
            score = _EXACT_SCORE[kind]
        elif a.get("acronym") or is_acronym(a.get("text") or ""):
            continue  # an acronym is matched whole or not at all
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
        cand = _entry_fields(e)
        cand.update({"match": kind, "matched_on": a["text"], "alias_kind": a["kind"],
                     "alias_source": a.get("source"), "maps_to": a.get("maps_to"), "score": round(score, 3)})
        if a.get("orpha"):
            cand["alias_orpha"] = a["orpha"]
        if a.get("member_of"):
            cand["member_of"] = a["member_of"]
        if a.get("qualifier_part"):
            cand["qualifier_part"] = a["qualifier_part"]
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


def _part_role(e: Dict[str, Any], part_norm: str, whole_norm: str) -> Optional[str]:
    """How a bracket/slash part of a published name may match: 'part', 'wider' or None (not a name).

    The alias layer records which parts are names: a part it kept is a synonym
    (scope 'wider' for the unqualified head of a qualified entry); a part it
    dropped is a qualifier or a gloss (重型, I型、Ⅱ型, 冷吡啉) and matches nothing.
    Without the alias layer no part is trusted to mean the disease.
    """
    if part_norm == whole_norm:
        return "exact"
    aliases = load_aliases()
    if not aliases["entries"]:
        return None
    a = aliases["_parts"].get((e["list"], e["no"], part_norm))
    if a is None:
        return None
    if a.get("acronym") or is_acronym(a.get("text") or ""):
        return "acronym"
    return "wider" if a.get("scope") == "wider" else "part"


def match(query: str, entries: Optional[Sequence[Dict[str, Any]]] = None, limit: int = 5) -> List[Dict[str, Any]]:
    """Best list entries for one name, matched on the published names.

    Kinds: exact | part (a synonym part) | wider (the unqualified head of a
    qualified entry) | acronym | contains | similar | overlap.
    """
    q = norm(query)
    if not q:
        return []
    entries = entries if entries is not None else load()["diseases"]
    scored = []
    for e in entries:
        best = None  # (score, kind, on)
        for field in ("name_zh", "name_en"):
            alts = alternates(e.get(field) or "")
            whole = alts[0] if alts else ""
            for a in alts:
                role = _part_role(e, a, whole) if q == a else None
                if q == a and role is not None:
                    cand = ({"exact": 1.0, "part": 0.95, "wider": 0.9, "acronym": 0.6}[role], role, e[field])
                elif q == a:
                    continue  # a qualifier or gloss on its own (重型, 冷吡啉): not a disease name
                elif a != whole:
                    continue  # fuzzy matching runs on whole published names only
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
        row = _entry_fields(e)
        row.update({"match": kind, "matched_on": on, "score": round(score, 3)})
        out.append(row)
    return out


def status_of(hits: Sequence[Dict[str, Any]]) -> str:
    if any(h["match"] in ON_LIST_KINDS for h in hits):
        return "on_list"
    if any(h["match"] in QUALIFIED_KINDS for h in hits):
        return "qualified"
    return "possible" if hits else "not_found"


def by_orpha(orpha_ids: Sequence[str]) -> List[Dict[str, Any]]:
    """List entries reached by an ORPHAcode: the entry's own link, a subtype or a group member (by Orphanet)."""
    want = {str(o).upper().replace("ORPHANET:", "ORPHA:") for o in orpha_ids if o}
    if not want:
        return []
    out: Dict[tuple, Dict[str, Any]] = {}
    for e in load_aliases()["entries"]:
        key = (e["list"], e["no"])
        if e.get("orpha") in want:
            row = _entry_fields(e)
            kind = "wider" if e.get("qualifier") else "orpha_exact"
            row.update({"match": kind, "matched_on": e["orpha"], "score": 1.0 if kind == "orpha_exact" else 0.9})
            out[key] = row
            continue
        for a in e.get("aliases") or []:
            if a.get("orpha") in want and not a.get("acronym"):
                kind = _alias_kind(a)
                if kind not in ("subtype", "group_member", "group_filed"):
                    continue
                row = _entry_fields(e)
                row.update({"match": kind, "matched_on": f"{a['orpha']} ({a['text']})", "alias_kind": a["kind"],
                            "alias_source": a.get("source"), "alias_orpha": a["orpha"],
                            "score": _EXACT_SCORE[kind]})
                if a.get("member_of"):
                    row["member_of"] = a["member_of"]
                if key not in out or row["score"] > out[key]["score"]:
                    out[key] = row
                break
    return sorted(out.values(), key=lambda m: (-m["score"], m["list"], m["no"]))


def lookup(names: Sequence[str], limit: int = 3, orpha_ids: Sequence[str] = ()) -> Dict[str, Any]:
    """Match one disease (ORPHAcodes first, then names); strongest first.

    `on_list` is True only for an on_list match (never for an acronym, a
    containment or a wider name); `status` is on_list | qualified | possible |
    not_found, and `qualifier` says what a qualified entry covers.
    """
    seen: Dict[tuple, Dict[str, Any]] = {}
    for m in by_orpha(orpha_ids):
        seen[(m["list"], m["no"])] = dict(m, query=m["matched_on"])
    for n in names:
        if not n:
            continue
        for m in list(alias_hits(n, limit=limit)) + list(match(n, limit=limit)):
            key = (m["list"], m["no"])
            m = dict(m, query=n)
            if key not in seen or m["score"] > seen[key]["score"]:
                seen[key] = m
    rows = sorted(seen.values(), key=lambda m: -m["score"])
    status = status_of(rows)
    on = [m for m in rows if m["match"] in ON_LIST_KINDS]
    if status == "on_list":
        shown = on
    elif status == "qualified":
        shown = [m for m in rows if m["match"] in QUALIFIED_KINDS]
    else:
        shown = rows[:limit]
    out = {"on_list": bool(on), "status": status, "matches": shown}
    if status == "qualified":
        out["qualifier"] = [{"list": m["list"], "no": m["no"], **(m.get("qualifier") or {})} for m in shown]
    return out


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


def _describe(h: Dict[str, Any]) -> List[str]:
    """Two lines per match: the entry, then exactly what matched and what that means."""
    lines = [f"  {h['list_name']} #{h['no']}  {h['name_zh']} / {h['name_en']}  [{h['match']}"
             + (f", score {h['score']}" if h["match"] not in ("exact", "orpha_exact") else "") + "]"]
    kind = h["match"]
    q = h.get("qualifier") or {}
    if kind == "wider":
        lines.append(f"      matched '{h['matched_on']}', which is wider than this entry: the list names only "
                     f"{q.get('zh') or '?'} ({q.get('covers') or q.get('en') or '?'}). It is on the list only if the "
                     "patient's diagnosis is that subtype")
    elif kind == "subtype":
        lines.append(f"      matched '{h['matched_on']}'"
                     + (f" ({h['alias_orpha']})" if h.get("alias_orpha") else "")
                     + f", a subtype inside the entry's qualifier {h.get('qualifier_part') or q.get('zh') or ''}"
                     + f" — {h.get('alias_source')}")
    elif kind == "group_member":
        lines.append(f"      matched '{h['matched_on']}'"
                     + (f" ({h['alias_orpha']})" if h.get("alias_orpha") else "")
                     + f", a member of the group this entry names — {h.get('alias_source')}")
    elif kind == "group_filed":
        lines.append(f"      '{h['matched_on']}'" + (f" ({h['alias_orpha']})" if h.get("alias_orpha") else "")
                     + " is filed under the group this entry names in Orphanet's classification, but its name does "
                       "not carry the group's name: whether the list entry covers it is not stated — verify with the "
                       "entry's 诊疗指南 or the care team")
    elif kind == "acronym":
        lines.append(f"      '{h['matched_on']}' is an acronym ({h.get('alias_kind') or 'part of the published name'}"
                     + (f", {h.get('alias_source')}" if h.get("alias_source") else "")
                     + "); acronyms are shared by unrelated diseases, so this is not a match on the disease — "
                       "give the full name")
    elif h.get("alias_kind"):
        lines.append(f"      matched '{h['matched_on']}' — {h['alias_kind']}"
                     + (f", maps to the official name {h['maps_to']}" if h.get("maps_to") else "")
                     + f" ({h.get('alias_source')})")
    elif kind in ("exact", "part"):
        lines.append(f"      matched the published name '{h['matched_on']}' "
                     f"({LIST_NAMES.get(h['list'])} #{h['no']})")
    elif h.get("matched_on"):
        how = {"contains": "shares a stretch of text with", "similar": "is spelled like",
               "overlap": "shares a stretch of text with"}.get(kind, "resembles")
        lines.append(f"      your text {how} the published name '{h['matched_on']}' — not a match on the disease")
    if q and kind not in ("wider",):
        lines.append(f"      note: this entry covers {q.get('covers') or q.get('zh')}")
    return lines


def _run(args: argparse.Namespace) -> Outcome:
    words = [w for w in args.query]
    if words and words[0].lower() in ("hospitals", "hospital", "医院", "协作网"):
        return _hospitals(args, " ".join(words[1:]).strip())
    query = " ".join(words).strip()
    if not query:
        raise UsageError("give a disease name (Chinese or English), or `hospitals [--province 浙江]`")
    if getattr(args, "province", None):
        raise UsageError("--province goes with `zebra china hospitals`")
    if args.limit < 1:
        raise UsageError("--limit must be at least 1")
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
    status = status_of(hits)
    if status == "on_list":
        shown = [h for h in hits if h["match"] in ON_LIST_KINDS]
    elif status == "qualified":
        shown = [h for h in hits if h["match"] in QUALIFIED_KINDS]
    else:
        shown = hits
    counts = (data.get("provenance") or {}).get("counts") or {}
    alias_prov = (load_aliases().get("provenance") or {})
    warnings: List[str] = []
    lines = []
    if status == "not_found":
        # P1f: never say "not on the list" for a name that simply was not matched.
        lines.append(f"{query}: no name on 第一批 (2018, {counts.get('1')}) or 第二批 (2023, {counts.get('2')}) "
                     "罕见病目录, and no alias of one, matched this text — which is not the same as the disease "
                     "being absent from the lists. Resolve the name first (`zebra disease <name>`) and query the "
                     "official Chinese or English name, or an ORPHA/OMIM id's name.")
    else:
        lines.append(f"{query}: {STATUS_TEXT[status]}")
        for h in shown:
            lines += _describe(h)
    if status == "qualified":
        for h in shown:
            q = h.get("qualifier") or {}
            warnings.append(f"'{query}' is wider than {h['list_name']} #{h['no']} {h['name_zh']}: the list covers "
                            f"{q.get('covers') or q.get('zh')}; do not say the disease is on the list unless the "
                            "diagnosis is that subtype")
    if any(h["match"] == "acronym" for h in shown):
        warnings.append(f"'{query}' matched only as an acronym; an acronym is not a disease name (it is shared by "
                        "unrelated diseases): give the full name before saying anything about the lists")
    zh = zh_name_hits(query, limit=3)
    if zh:
        lines.append("Orphanet Chinese names close to this text (use with `zebra disease <ORPHA id>`):")
        for z in zh:
            lines.append(f"  {z['orpha']}  {z['name_zh']}  [{z['match']}, similarity {z['similarity']}]")
    result = {"query": query, "status": status, "matches": shown, "orphanet_zh_candidates": zh,
              "provenance": provenance_brief(data),
              "alias_provenance": {k: alias_prov.get(k) for k in ("retrieved_at", "builder", "sources", "counts", "rules")
                                   if k in alias_prov}}
    if status == "qualified":
        result["qualifier"] = [{"list": h["list"], "no": h["no"], **(h.get("qualifier") or {})} for h in shown]
    return Outcome(result, sources=source_for(data, shown) + alias_sources(shown, zh),
                   warnings=warnings, text="\n".join(lines), query={"query": query})


def alias_sources(hits: Sequence[Dict[str, Any]], zh: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Provenance rows for the alias layer and the Chinese name index, when either was used."""
    recs: List[Dict[str, Any]] = []
    prov = load_aliases().get("provenance") or {}
    if any(h.get("alias_kind") for h in hits):
        recs.append(source_record("China rare disease list alias layer",
                                  f"china_disease_aliases.json ({(prov.get('counts') or {}).get('aliases')} aliases)",
                                  url=((prov.get("sources") or {}).get("orphanet_zh_bulk") or {}).get("url"),
                                  note="bundled; every alias carries its kind and source, folk names are marked "
                                       f"folk_name (built {prov.get('retrieved_at')} by {prov.get('builder', '?').split(' ')[0]})"))
    if zh:
        zprov = load_zh_names().get("provenance") or {}
        recs.append(source_record("Orphanet Chinese preferred terms",
                                  f"{zprov.get('names_kept')} ORPHAcodes", url=zprov.get("source_url"),
                                  note=f"bundled copy retrieved {zprov.get('retrieved_at')}, "
                                       f"Orphanet zh dataset {zprov.get('dataset_date')}, "
                                       f"{zprov.get('licence')}"))
    return recs


# ------------------------------------------------------------------ hospitals

def load_hospitals() -> Dict[str, Any]:
    if "hospitals" not in _cache:
        if not os.path.exists(HOSPITALS_FILE):
            raise ListUnavailable(f"not bundled: {os.path.basename(HOSPITALS_FILE)} is missing")
        with open(HOSPITALS_FILE, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        if not data.get("hospitals"):
            raise ListUnavailable(f"{os.path.basename(HOSPITALS_FILE)} holds no hospitals")
        _cache["hospitals"] = data
    return _cache["hospitals"]


_PROVINCE_SUFFIX = re.compile(r"(省|市|自治区|壮族自治区|回族自治区|维吾尔自治区|特别行政区)$")


def province_key(text: str) -> str:
    """浙江省 / 浙江 / 广西壮族自治区 / 广西 -> one key."""
    t = unicodedata.normalize("NFKC", text or "").strip()
    if t in ("兵团", "新疆兵团", "生产建设兵团", "新疆建设兵团"):
        return "新疆生产建设兵团"
    t = _PROVINCE_SUFFIX.sub("", t)
    for long, short in (("广西壮族", "广西"), ("宁夏回族", "宁夏"), ("新疆维吾尔", "新疆"), ("内蒙古", "内蒙古")):
        if t.startswith(long):
            t = short
    return t


def hospital_sources(data: Dict[str, Any], used: Optional[Sequence[str]] = None) -> List[Dict[str, Any]]:
    p = data.get("provenance") or {}
    recs = []
    for sid, s in (p.get("sources") or {}).items():
        if used is not None and sid not in used:
            continue
        recs.append(source_record("全国罕见病诊疗协作网 (NHC)", f"{s.get('document_no') or s.get('title') or sid}",
                                  url=s.get("url"),
                                  note=f"{s.get('title', '')}; {s.get('published', '')}; retrieved "
                                       f"{s.get('retrieved_at')}; sha256 {str(s.get('sha256'))[:12]}…; "
                                       f"{s.get('rows')} rows" + (f"; origin {s['origin']}" if s.get("origin") else "")))
    return recs


def hospitals(province: Optional[str] = None, roles: Sequence[str] = ()) -> Outcome:
    """协作网 hospitals: all of one province, or the national and provincial lead hospitals."""
    try:
        data = load_hospitals()
    except (ListUnavailable, OSError, ValueError) as err:
        return Outcome({"status": "unavailable", "province": province, "hospitals": []},
                       warnings=[f"全国罕见病诊疗协作网 hospital list {err}; the network's membership was not "
                                 "consulted — it is not evidence that a hospital is or is not a member"])
    # current membership only: a 2019 name absent from the later list is kept in the file
    # for the record (status removed_2024) but is not offered as a network hospital
    current = [h for h in data["hospitals"] if not str(h.get("status") or "").startswith("removed")]
    gone = len(data["hospitals"]) - len(current)
    rows = current
    note = (data.get("provenance") or {}).get("coverage")
    if province:
        key = province_key(province)
        sel = [h for h in rows if province_key(h.get("province") or "") == key]
        if not sel:
            known = sorted({province_key(h.get("province") or "") for h in rows})
            return Outcome({"status": "province_not_found", "province": province, "hospitals": [],
                            "provinces": known, "coverage": note},
                           sources=hospital_sources(data),
                           warnings=[f"no 协作网 hospital is listed under '{province}' in the bundled list; provinces "
                                     f"as published: {', '.join(known)}"])
    else:
        sel = [h for h in rows if h.get("role") in (roles or ("national_lead", "provincial_lead"))]
    used = sorted({h.get("source_id") for h in sel if h.get("source_id")})
    out = [{k: h.get(k) for k in ("name", "province", "role", "role_zh", "status", "name_2019", "source_id")
            if h.get(k) is not None} for h in sel]
    counts = (data.get("provenance") or {}).get("counts") or {}
    return Outcome({"status": "ok", "province": province, "count": len(out), "hospitals": out,
                    "counts": {k: counts.get(k) for k in ("total", "national_lead", "provincial_lead", "member")},
                    "not_in_current_list": gone, "coverage": note},
                   sources=hospital_sources(data, used or None))


def _hospitals(args: argparse.Namespace, rest: str) -> Outcome:
    province = (getattr(args, "province", None) or rest or "").strip() or None
    out = hospitals(province)
    r = out.result
    lines: List[str] = []
    if r["status"] == "unavailable":
        lines.append("全国罕见病诊疗协作网 hospital list: not bundled here (see warnings); not consulted")
    elif r["status"] == "province_not_found":
        lines.append(f"{province}: no 协作网 hospital listed under this name; provinces: {', '.join(r['provinces'])}")
    else:
        head = (f"全国罕见病诊疗协作网 hospitals in {province}: {r['count']}" if province else
                f"全国罕见病诊疗协作网 lead hospitals (national + provincial): {r['count']} — "
                "`zebra china hospitals --province <省>` lists every member in a province")
        lines.append(head)
        for h in r["hospitals"]:
            lines.append(f"  {h.get('province', '?')}  {h['name']}  [{h.get('role_zh') or HOSPITAL_ROLES.get(h.get('role'), h.get('role'))}]"
                         + (f" (2019 name: {h['name_2019']})" if h.get("name_2019") else "")
                         + (" (added 2024)" if h.get("status") == "added_2024" else ""))
        if r.get("coverage"):
            cov = r["coverage"]
            lines.append("coverage: " + (cov if len(cov) <= 400 else cov[:400].rsplit(" ", 1)[0]
                                         + " … (full text in result.coverage)"))
    out.text = "\n".join(lines)
    out.query = {"hospitals": True, "province": province}
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("china", help="China's national rare disease lists (第一批 2018 / 第二批 2023): is a disease on "
                                     "them; `china hospitals [--province 浙江]`: 全国罕见病诊疗协作网 hospitals")
    p.add_argument("query", nargs="+", help="disease name, Chinese or English; or `hospitals`")
    p.add_argument("--limit", type=int, default=5)
    p.add_argument("--province", help="with `hospitals`: one province (浙江 / 浙江省)")
    p.set_defaults(func=_run)
