#!/usr/bin/env python3
"""Build zebra/data/china_disease_aliases.json: the alias layer over China's national rare disease lists.

This is the builder `provenance.builder` names (review E-10: the first version
lived in a session scratchpad and was lost). It rebuilds the whole layer from
the bundled lists and Orphadata (api.orphadata.com, CC-BY-4.0):

    python3 -I tools/china_list/build_china_aliases.py --cache <dir> \
        [--out zebra/data/china_disease_aliases.json] [--compare <old json>]

`--cache` holds every Orphadata response as fetched (one JSON file per URL), so
a rerun with a warm cache makes no request. `--compare` prints how the rebuilt
links and aliases differ from an earlier file (that is how this script was
checked against the lost one: see provenance.rebuild_check).

Pass 1 — links and retrieved aliases (the original rules, in order):
  a. normalised name_en == an English Orphanet preferred term
  b. a '/'- or bracket-separated part of name_en == an English preferred term
  c. normalised name_zh (or a part) == a Chinese preferred term (zh dataset
     2020-06-01); when several codes share it, the tie is broken only by an
     exact English preferred term or synonym among them, otherwise no link
  d. name_en (or a part) == an English Orphanet SYNONYM, checked for the 4 best
     difflib candidates of each entry no earlier rule linked
  No link is ever made on similarity. Aliases: the published names, their
  bracket/slash parts (minus the reviewed fragments in DROPPED_PARTS), and the
  linked code's English and Chinese preferred terms and synonyms.

Pass 2 — what a match on each alias means (review E-1):
  * QUALIFIED: three entries name only part of a disease — 糖原累积病（I型、Ⅱ型）,
    帕金森病（青年型、早发型）, 地中海贫血（重型）. Their unqualified head
    (帕金森病, 地中海贫血, Glycogen Storage Disease) and every Orphanet alias of
    a link that is wider than the entry get `scope: "wider"`: a match on them
    means "only a subtype of this is on the list", never "on the list".
  * Subtypes inside a qualifier are added from Orphanet with their own
    ORPHAcode (`kind: orphanet_subtype_*`, `orpha`, `qualifier_part`): the
    children of ORPHA:79201 whose names say GSD type I or II (and their
    descendants: Pompe disease is GSD II), young-onset Parkinson disease,
    and the beta-thalassemia *major* records.
  * Group entries (Orphanet "Group of disorders") get their Orphanet
    descendants (`kind: orphanet_member_*`, `orpha`, `member_of`), so
    "Duchenne muscular dystrophy" is found as a member of 进行性肌营养不良
    through Orphanet's own classification rather than through a hand-entered
    folk name (review E-9). Depth and size are capped (MEMBER_DEPTH,
    MEMBER_CAP); a capped group says so in provenance.
  * `acronym: true` marks short all-capital aliases (CAD, HSP, PV, GSD, DMD,
    SMA, ...): an acronym is ambiguous across diseases and never counts as
    "on the list" by itself.

Pass 3 — folk names (FOLK): hand-entered, `kind: folk_name`, `maps_to` the
official name, and a source that says they come from no dataset.
"""

from __future__ import annotations

import argparse
import difflib
import threading
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
import re
import sys
import time
import unicodedata
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
LISTS = os.path.join(ROOT, "zebra", "data", "china_rare_diseases.json")
OUT = os.path.join(ROOT, "zebra", "data", "china_disease_aliases.json")
API = "https://api.orphadata.com"
ZH_DATASET = "2020-06-01"
PACE = 0.25  # seconds between two Orphadata requests
LICENCE = "CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0)"
MEMBER_DEPTH = 3
MEMBER_CAP = 150
WORKERS = 6
RULE_D_CANDIDATES = 4  # difflib candidates per unlinked entry whose synonyms are checked (rule d)

# Bracket/slash parts that are not alternative disease names. Reviewed by hand
# against the published lists; a length or substring heuristic is not used
# because it would also discard genuine short names (脆骨病, 白塞病, HAE, LAM).
DROPPED_PARTS = {
    (1, 35, "I型、II型"): "severity/type qualifier of 糖原累积病, not a disease name",
    (1, 35, "Type I、II"): "severity/type qualifier of Glycogen Storage Disease, not a disease name",
    (1, 87, "青年型、早发型"): "onset qualifier of 帕金森病, not a disease name",
    (1, 87, "Young-onset , Early-onset"): "onset qualifier of Parkinson Disease, not a disease name",
    (2, 18, "冷吡啉"): "protein name (cryopyrin) glossed inside the disease name, not a disease name",
    (2, 18, "冷炎素"): "the inline gloss of 冷吡啉 (cryopyrin), a protein name, not a disease name",
    (2, 18, "相关周期性综合征"): "dangling remainder left after the inline gloss is split out, not a disease name",
    (2, 78, "重型"): "severity qualifier of 地中海贫血, not a disease name",
}

# Entries whose published name restricts the disease to a subtype (E-1).
# head: the unqualified name parts; subtype_*: how an Orphanet record is
# recognised as inside the qualifier (both regexes must match one of its names).
QUALIFIED: Dict[Tuple[int, int], Dict[str, Any]] = {
    (1, 35): {"qualifier_zh": "I型、Ⅱ型", "qualifier_en": "Type I、II",
              "covers": "glycogen storage disease types I and II only (Ⅱ型 = Pompe disease)",
              "head": ["糖原累积病", "Glycogen Storage Disease"],
              "subtype_head": r"glycogen storage disease|glycogenosis|\bgsd\b",
              # type II is Pompe disease (acid maltase deficiency); "GSD IIb" (Danon disease,
              # LAMP2) is a historical label for a different disease and is not taken
              "subtype_qualifier": r"\btype\s*(?:1|i|ia|ib|ic|id|2|ii)\b|\bgsd\s*(?:1|i|ia|ib|2|ii)\b",
              "from_classification_of": 79201},
    (1, 87): {"qualifier_zh": "青年型、早发型", "qualifier_en": "Young-onset , Early-onset",
              "covers": "young-onset / early-onset Parkinson disease only",
              "head": ["帕金森病", "Parkinson Disease"],
              "subtype_head": r"parkinson(?:'s)? disease",
              "subtype_qualifier": r"young[- ]onset|early[- ]onset",
              "candidates_from_terms": r"parkinson"},
    (2, 78): {"qualifier_zh": "重型", "qualifier_en": "major",
              "covers": "thalassaemia major only (not trait / minor or intermedia)",
              "head": ["地中海贫血"],
              "subtype_head": r"thalass(?:a)?emia",
              # "thalassemia major" itself, not "... with alpha-thalassemia as a major feature"
              "subtype_qualifier": r"thalass(?:a)?emia\W{0,2}major\b",
              "candidates_from_terms": r"thalass"},
}

# Subtype names formed from the published name itself by applying its qualifier
# (糖原累积病（I型、Ⅱ型） -> 糖原累积病I型 / 糖原累积病Ⅱ型): no outside knowledge.
_DERIVED = {
    (1, 35): [("糖原累积病I型", "zh", "I型"), ("糖原累积病Ⅱ型", "zh", "Ⅱ型"), ("糖原累积病1型", "zh", "I型"),
              ("糖原累积病2型", "zh", "Ⅱ型"), ("Glycogen Storage Disease Type I", "en", "Type I"),
              ("Glycogen Storage Disease Type II", "en", "Type II")],
    (1, 87): [("青年型帕金森病", "zh", "青年型"), ("早发型帕金森病", "zh", "早发型"),
              ("Young-onset Parkinson Disease", "en", "Young-onset"), ("Early-onset Parkinson Disease", "en", "Early-onset")],
    (2, 78): [("重型地中海贫血", "zh", "重型")],
}
for _k, _v in _DERIVED.items():
    QUALIFIED[_k]["derived"] = _v

# Published English names with an evident misspelling: the corrected spelling is an extra alias;
# the published name is kept as it is.
CORRECTED = [
    ((1, 14), "Cardiac Ion Channelopathies", "misspells 'Cardiac' as 'Cardic'"),
    ((1, 71), "Methylmalonic Acidemia", "reads 'Academia' for 'Acidemia'"),
    ((1, 72), "Mitochondrial Encephalomyopathy", "misspells 'Mitochondrial' as 'Mitochodrial'"),
]

# leading words that do not name the disease family of a group (Progressive muscular dystrophy ->
# "muscular dystrophy"; Primary immunodeficiency -> "immunodeficiency")
_GROUP_LEADS = ("progressive", "hereditary", "inherited", "genetic", "primary")


def group_core(name: str) -> str:
    words = (name or "").lower().split()
    while len(words) > 1 and words[0] in _GROUP_LEADS:
        words = words[1:]
    return " ".join(words)


FOLK_SOURCE_1 = "folk name, hand-entered; not from a retrieved dataset"
FOLK_SOURCE_2 = ("folk name or Chinese transliteration, hand-entered; not from a retrieved dataset "
                 "(second builder pass, 2026-10-06)")
FOLK_SOURCE_3 = ("folk name, hand-entered; not from a retrieved dataset (third builder pass, 2026-10-06, "
                 "review CP1-7)")
FOLK: List[Tuple[Tuple[int, int], str, str, Optional[str]]] = [
    ((1, 4), "渐冻症", FOLK_SOURCE_1, None),
    ((1, 4), "渐冻人症", FOLK_SOURCE_1, None),
    ((1, 5), "快乐木偶", FOLK_SOURCE_1, None),
    ((1, 5), "快乐木偶综合征", FOLK_SOURCE_1, None),
    ((1, 86), "瓷娃娃", FOLK_SOURCE_1, None),
    ((1, 93), "小胖威利", FOLK_SOURCE_1, None),
    ((1, 93), "普拉德-威利综合征", FOLK_SOURCE_1, None),
    ((1, 98), "杜氏肌营养不良", FOLK_SOURCE_1, None),
    ((1, 98), "杜兴氏肌营养不良", FOLK_SOURCE_1, None),
    ((1, 98), "DMD", FOLK_SOURCE_1, None),
    ((1, 105), "德拉韦综合征", FOLK_SOURCE_1, None),
    ((1, 110), "SMA", FOLK_SOURCE_1, None),
    ((2, 72), "瑞特综合征", FOLK_SOURCE_2, None),
    ((2, 72), "雷特综合征", FOLK_SOURCE_2, None),
    ((1, 39), "蝴蝶宝贝", FOLK_SOURCE_2, None),
    ((1, 39), "蝴蝶宝宝", FOLK_SOURCE_2, None),
    ((1, 73), "黏宝宝", FOLK_SOURCE_2, None),
    ((1, 73), "粘宝宝", FOLK_SOURCE_2, None),
    ((1, 2), "月亮孩子", FOLK_SOURCE_2, None),
    ((1, 90), "不食人间烟火的孩子", FOLK_SOURCE_2, None),
    # Pompe disease is GSD type II: Orphanet ORPHA:365 lists both "Pompe disease" and
    # "GSD type II" as synonyms; 庞贝病 is the usual Chinese name, which Orphanet's zh
    # dataset does not carry (its zh term is 酸性麦芽糖酶缺乏所致糖原贮积病).
    ((1, 35), "庞贝病", FOLK_SOURCE_3, "Ⅱ型"),
    ((1, 35), "庞贝氏症", FOLK_SOURCE_3, "Ⅱ型"),
]
FOLK_ORPHA = {"庞贝病": "ORPHA:365", "庞贝氏症": "ORPHA:365"}

ACRONYM_RE = re.compile(r"^[A-Z][A-Za-z0-9]*(?:-[A-Z0-9][A-Za-z0-9]*)?$")


def is_acronym(text: str) -> bool:
    """Short and mostly capitals: CAD, HSP, aHUS, X-ALD, SMAX1, GEP-NEN, CADASIL."""
    t = (text or "").strip()
    if not t or len(t) > 8 or " " in t or re.search(r"[^\x00-\x7f]", t):
        return False
    caps = sum(1 for c in t if c.isupper())
    return caps >= 2 and caps >= len(re.sub(r"[^A-Za-z]", "", t)) - 1


def norm(text: str) -> str:
    """Same folding as zebra.commands.china.norm (kept in step; the test suite checks it)."""
    t = unicodedata.normalize("NFKC", text or "").lower()
    t = t.replace("’", "'").replace("‘", "'").replace("＇", "'")
    t = re.sub(r"'s\b", "", t)
    t = re.sub(r"\b([a-z]{4,})s(?=\s+(?:disease|syndrome)\b)", r"\1", t)  # parkinsons disease ~ parkinson's disease
    t = re.sub(r"([a-z])ae", r"\1e", t)  # thalassaemia / haemophilia / anaemia ~ American spelling
    t = re.sub(r"(?<=[a-z0-9])\s*(?:综合征|综合症)$", "syndrome", t)  # Pompe病 ~ Pompe disease
    t = re.sub(r"(?<=[a-z0-9])\s*(?:病|症)$", "disease", t)
    t = t.replace("ⅱ", "ii").replace("ⅲ", "iii").replace("ⅰ", "i")
    # "Type I、II" must not fold into "typeiii": a list separator between two numerals is kept as "|"
    t = re.sub(r"(?<=\b[ivx0-9])\s*[、,，/;]\s*(?=[ivx0-9]+\b)|(?<=[^a-z][ivx])\s*[、,，/;]\s*(?=[ivx0-9])", "|", t)
    return re.sub(r"[\s\-_'\",.;:()\[\]{}/、，。（）·]+", "", t)


def parts(name: str) -> List[str]:
    """Bracket/slash-separated parts of a published name, stripped, in order, whole name excluded."""
    t = unicodedata.normalize("NFKC", name or "")
    out: List[str] = []
    for p in re.split(r"[/()\[\]（）]", t):
        p = p.strip()
        if p and p != t.strip() and p not in out:
            out.append(p)
    return out


class Fetcher:
    def __init__(self, cache: str) -> None:
        self.cache = cache
        os.makedirs(cache, exist_ok=True)
        self.last = 0.0
        self.requests = 0
        self.log: Dict[str, Dict[str, Any]] = {}
        self.lock = threading.Lock()

    def prefetch(self, calls: Sequence[Tuple[str, Optional[Dict[str, str]]]]) -> None:
        """Warm the cache for many calls at once (Orphadata answers in ~3 s; requests start PACE apart)."""
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            list(pool.map(lambda c: self.get(c[0], c[1], allow_404=True, quiet=True), calls))

    def get(self, path: str, params: Optional[Dict[str, str]] = None, allow_404: bool = False,
            quiet: bool = False) -> Optional[Dict[str, Any]]:
        url = API + path + (("?" + urllib.parse.urlencode(params)) if params else "")
        key = hashlib.sha256(url.encode()).hexdigest()
        fn = os.path.join(self.cache, key + ".json")
        if os.path.exists(fn):
            with open(fn, encoding="utf-8") as fh:
                rec = json.load(fh)
        else:
            with self.lock:
                wait = PACE - (time.time() - self.last)
                self.last = max(time.time(), self.last + PACE)
                self.requests += 1
            if wait > 0:
                time.sleep(wait)
            req = urllib.request.Request(url, headers={"Accept": "application/json",
                                                       "User-Agent": "zebra-mod china alias builder"})
            status, body = 0, b""
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(req, timeout=90) as resp:
                        status, body = resp.status, resp.read()
                    break
                except urllib.error.HTTPError as err:
                    status, body = err.code, err.read()
                    if status not in (429, 500, 502, 503, 504):
                        break
                except OSError:
                    status = 0
                time.sleep(2 * (attempt + 1))
            rec = {"url": url, "status": status, "retrieved_at": now(),
                   "sha256": hashlib.sha256(body).hexdigest(), "body": body.decode("utf-8", "replace")}
            tmp = f"{fn}.{threading.get_ident()}.tmp"
            with open(tmp, "w", encoding="utf-8") as fh:
                json.dump(rec, fh, ensure_ascii=False)
            os.replace(tmp, fn)
        with self.lock:
            self.log[url] = {"status": rec["status"], "retrieved_at": rec["retrieved_at"], "sha256": rec["sha256"]}
        if quiet:
            return None
        if rec["status"] == 404 and allow_404:
            return None
        if rec["status"] != 200:
            raise SystemExit(f"{url}: HTTP {rec['status']}")
        return json.loads(rec["body"])


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).strftime("%Y-%m-%dT%H:%M:%SZ")


def detail(f: Fetcher, code: int, lang: str) -> Optional[Dict[str, Any]]:
    got = f.get(f"/rd-cross-referencing/orphacodes/{code}", {"lang": lang}, allow_404=True)
    if not got:
        return None
    r = (got.get("data") or {}).get("results") or {}
    return {"code": code, "pref": (r.get("Preferred term") or "").strip(),
            "syn": [s.strip() for s in (r.get("Synonym") or []) if s and s.strip()],
            "group": r.get("DisorderGroup"), "date": r.get("Date")}


def children(f: Fetcher, code: int) -> List[int]:
    got = f.get(f"/rd-classification/orphacodes/{code}/hchids", allow_404=True)
    kids: List[int] = []
    for r in ((got or {}).get("data") or {}).get("results") or []:
        for c in r.get("childs") or []:
            if c not in kids:
                kids.append(int(c))
    return kids


def descendants(f: Fetcher, root: int, depth: int, cap: int) -> Tuple[List[Tuple[int, int, int]], bool]:
    """[(code, parent, level)] breadth-first under `root`, at most `cap` codes and `depth` levels."""
    seen = {root}
    out: List[Tuple[int, int, int]] = []
    frontier = [root]
    capped = False
    for level in range(1, depth + 1):
        nxt: List[int] = []
        f.prefetch([(f"/rd-classification/orphacodes/{p}/hchids", None) for p in frontier])
        for parent in frontier:
            for c in children(f, parent):
                if c in seen:
                    continue
                if len(out) >= cap:
                    capped = True
                    break
                seen.add(c)
                out.append((c, parent, level))
                nxt.append(c)
        frontier = nxt
        if not frontier:
            break
    return out, capped


def link_entries(f: Fetcher, entries: List[Dict[str, Any]], en_bulk: Dict[int, str], zh_bulk: Dict[int, str],
                 notes: List[str]) -> Dict[Tuple[int, int], Tuple[int, str]]:
    en_index: Dict[str, List[int]] = {}
    for code, term in en_bulk.items():
        en_index.setdefault(norm(term), []).append(code)
    zh_index: Dict[str, List[int]] = {}
    for code, term in zh_bulk.items():
        zh_index.setdefault(norm(term), []).append(code)
    links: Dict[Tuple[int, int], Tuple[int, str]] = {}
    for e in entries:
        key = (e["list"], e["no"])
        hits = en_index.get(norm(e["name_en"])) or []
        if len(hits) == 1:
            links[key] = (hits[0], "exact English preferred term")
            continue
        found = None
        for p in parts(e["name_en"]):
            h = en_index.get(norm(p)) or []
            if len(h) == 1:
                found = (h[0], "exact English preferred term of a bracketed/slash-separated part of the list's English name")
                break
        if found:
            links[key] = found
            continue
        zh_hits: List[int] = []
        for n in [e["name_zh"]] + parts(e["name_zh"]):
            if (e["list"], e["no"], n) in DROPPED_PARTS:
                continue
            for c in zh_index.get(norm(n)) or []:
                if c not in zh_hits:
                    zh_hits.append(c)
        if len(zh_hits) == 1:
            links[key] = (zh_hits[0], f"exact Chinese preferred term (Orphanet zh dataset {ZH_DATASET})")
            continue
        if len(zh_hits) > 1:
            wanted = {norm(e["name_en"])} | {norm(p) for p in parts(e["name_en"])}
            winners = []
            for c in zh_hits:
                d = detail(f, c, "en")
                if d and ({norm(d["pref"])} | {norm(s) for s in d["syn"]}) & wanted:
                    winners.append(c)
            if len(winners) == 1:
                links[key] = (winners[0], f"exact Chinese preferred term (Orphanet zh dataset {ZH_DATASET})")
            else:
                notes.append(f"{e['list']}/{e['no']} {e['name_zh']} -> left null: "
                             + ", ".join(f"ORPHA:{c}" for c in zh_hits)
                             + f" share the Chinese preferred term '{zh_bulk[zh_hits[0]]}' and none of them carries "
                               f"an English preferred term or synonym equal to '{e['name_en']}', so no exact rule "
                               "can choose between them")
    # rule d: synonyms, for the RULE_D_CANDIDATES best difflib candidates of each entry still unlinked
    terms = list(en_bulk.items())
    lowered = [t.lower() for _, t in terms]
    todo: List[Tuple[Dict[str, Any], List[int]]] = []
    for e in entries:
        key = (e["list"], e["no"])
        if key in links:
            continue
        cands = difflib.get_close_matches(e["name_en"].lower(), lowered, n=RULE_D_CANDIDATES, cutoff=0.0)
        codes: List[int] = []
        for c in cands:
            for code, t in terms:
                if t.lower() == c and code not in codes:
                    codes.append(code)
                    break
        todo.append((e, codes))
    f.prefetch([(f"/rd-cross-referencing/orphacodes/{c}", {"lang": "en"}) for _e, cs in todo for c in cs])
    for e, codes in todo:
        key = (e["list"], e["no"])
        names = [e["name_en"]] + parts(e["name_en"])
        wanted = {norm(n): n for n in names}
        for code in codes:
            d = detail(f, code, "en")
            if not d:
                continue
            hit = next((s for s in d["syn"] if norm(s) in wanted), None)
            if hit:
                whole = norm(hit) == norm(e["name_en"])
                links[key] = (code, "exact English synonym" if whole else
                              "exact English synonym of a bracketed/slash-separated part of the list's English name")
                break
    return links


def build(f: Fetcher, compare: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    with open(LISTS, encoding="utf-8") as fh:
        lists = json.load(fh)
    entries = lists["diseases"]
    list_name = {1: "第一批罕见病目录 (2018)", 2: "第二批罕见病目录 (2023)"}
    en_raw = f.get("/rd-cross-referencing/orphacodes", {"lang": "en"})
    zh_raw = f.get("/rd-cross-referencing/orphacodes", {"lang": "zh"})
    en_bulk = {int(r["ORPHAcode"]): r["Preferred term"].strip() for r in en_raw["data"]["results"] if r.get("Preferred term")}
    zh_bulk = {int(r["ORPHAcode"]): r["Preferred term"].strip() for r in zh_raw["data"]["results"] if r.get("Preferred term")}
    # an obsolete code is never linked (2/2 获得性血友病 matched "OBSOLETE: Acquired hemophilia", ORPHA:73274,
    # by its Chinese term; Orphanet "refers" it to a broader group, which is not the same disease)
    obsolete = {c for c, t in en_bulk.items() if t.upper().startswith("OBSOLETE")}
    en_bulk = {c: t for c, t in en_bulk.items() if c not in obsolete}
    zh_bulk = {c: t for c, t in zh_bulk.items() if c not in obsolete}
    notes: List[str] = []
    links = link_entries(f, entries, en_bulk, zh_bulk, notes)
    f.prefetch([(f"/rd-cross-referencing/orphacodes/{c}", {"lang": lg}) for c, _h in links.values() for lg in ("en", "zh")])
    out_entries: List[Dict[str, Any]] = []
    missing_zh: List[str] = []
    capped_groups: List[str] = []
    counts = {"list_entries": len(entries), "linked_to_orphanet": len(links), "aliases": 0, "folk_names": 0,
              "subtype_aliases": 0, "member_aliases": 0, "member_codes": 0, "member_named": 0, "member_filed": 0,
              "acronyms": 0, "wider_aliases": 0, "derived_subtype_aliases": 0, "corrected_spellings": 0}
    for e in entries:
        key = (e["list"], e["no"])
        src = f"China list {list_name[e['list']]} #{e['no']}"
        aliases: List[Dict[str, Any]] = []
        seen = set()

        def add(text: str, kind: str, source: str, **extra: Any) -> None:
            text = (text or "").strip()
            label = re.match(r"^(NON RARE IN EUROPE)\s*:\s*", text, re.I)
            if label:  # Orphanet's label, not part of the name (1/76 Multiple sclerosis)
                text = text[label.end():].strip()
                extra = dict(extra, orphanet_label=label.group(1).upper())
            if not text or (text, kind) in seen:
                return
            seen.add((text, kind))
            row = {"text": text, "kind": kind, "source": source}
            row.update({k: v for k, v in extra.items() if v is not None})
            aliases.append(row)

        add(e["name_zh"], "official_zh", src)
        add(e["name_en"], "official_en", src)
        for lang, field in (("zh", "name_zh"), ("en", "name_en")):
            for p in parts(e[field]):
                if (e["list"], e["no"], p) in DROPPED_PARTS:
                    continue
                add(p, f"official_part_{lang}", f"{src} (part of the published {'Chinese' if lang == 'zh' else 'English'} name)")
        row = {"list": e["list"], "no": e["no"], "name_zh": e["name_zh"], "name_en": e["name_en"],
               "orpha": None, "orpha_match": None, "orpha_name_en": None, "orpha_name_zh": None,
               "orpha_disorder_group": None}
        q = QUALIFIED.get(key)
        if key in links:
            code, how = links[key]
            den = detail(f, code, "en")
            dzh = detail(f, code, "zh")
            row.update({"orpha": f"ORPHA:{code}", "orpha_match": how, "orpha_name_en": den["pref"] if den else None,
                        "orpha_name_zh": dzh["pref"] if dzh else None,
                        "orpha_disorder_group": den["group"] if den else None})
            if den:
                add(den["pref"], "orphanet_en", f"Orphanet en preferred term, ORPHA:{code}")
            if dzh:
                add(dzh["pref"], "orphanet_zh", f"Orphanet zh preferred term, ORPHA:{code} (dataset {ZH_DATASET})")
            else:
                missing_zh.append(f"{e['list']}/{e['no']} {e['name_zh']} -> ORPHA:{code}")
            for s in (den or {}).get("syn") or []:
                add(s, "orphanet_synonym_en", f"Orphanet en synonym, ORPHA:{code}")
            for s in (dzh or {}).get("syn") or []:
                add(s, "orphanet_synonym_zh", f"Orphanet zh synonym, ORPHA:{code} (dataset {ZH_DATASET})")
        if q:
            row["qualifier"] = {"zh": q["qualifier_zh"], "en": q["qualifier_en"], "covers": q["covers"]}
            if row["orpha"]:
                row["orpha_scope"] = (f"wider than the entry: {row['orpha']} is the whole disease, the entry covers "
                                      f"{q['covers']}")
            heads = {norm(h) for h in q["head"]}
            for a in aliases:
                if a["kind"] in ("official_zh", "official_en"):
                    continue
                if a["kind"].startswith("official_part") and norm(a["text"]) not in heads:
                    continue
                a["scope"] = "wider"
                a["restricted_to"] = q["qualifier_zh"]
            _add_subtypes(f, row, q, en_bulk, add)
            # a subtype record can carry the unqualified name as a synonym (Orphanet's zh
            # synonym of Beta-thalassemia major is plain 地中海贫血): that text is still wider
            wider = {norm(a["text"]) for a in aliases if a.get("scope") == "wider"}
            for a in aliases:
                if a["kind"].startswith("orphanet_subtype") and norm(a["text"]) in wider:
                    a["scope"] = "wider"
                    a["restricted_to"] = q["qualifier_zh"]
        elif row["orpha_disorder_group"] == "Group of disorders":
            members, capped = descendants(f, int(row["orpha"].split(":")[1]), MEMBER_DEPTH, MEMBER_CAP)
            if capped:
                capped_groups.append(f"{e['list']}/{e['no']} {row['orpha']} (first {MEMBER_CAP} descendants kept, "
                                     f"depth <= {MEMBER_DEPTH})")
            core_en = group_core(row["orpha_name_en"] or "")
            core_zh = [x for x in {re.split(r"[（(/]", e["name_zh"])[0].strip(), (row.get("orpha_name_zh") or "").strip()}
                       if len(x) >= 2 and re.search(r"[㐀-鿿]", x)]
            for code, parent, level in members:
                # preferred terms from the bulk files (no per-code synonym calls: ~3 s each)
                en_name = en_bulk.get(code) or ""
                if not en_name or en_name.upper().startswith("OBSOLETE"):
                    continue
                zh_name = zh_bulk.get(code) or ""
                if not re.search(r"[㐀-鿿]", zh_name):
                    zh_name = ""  # the zh dataset holds some English text; it is not a Chinese name
                # Orphanet files diseases under a group for many reasons (cause, differential, a shared
                # feature): Gastroschisis sits under Short bowel syndrome and Whipple disease under Primary
                # immunodeficiency. Only a member whose own name carries the group's name is taken as a
                # form of the listed disease ("named"); the rest are "filed" and only offered (review P0-2).
                named = bool(core_en and re.search(r"(?<![a-z0-9])" + re.escape(core_en) + r"(?![a-z0-9])",
                                                   en_name.lower())) or \
                    bool(zh_name and any(c in zh_name for c in core_zh))
                counts["member_named" if named else "member_filed"] += 1
                via = (f"Orphanet classification: ORPHA:{code} is under {row['orpha']} "
                       f"{row['orpha_name_en']} (level {level}, parent ORPHA:{parent})"
                       + (f"; its name carries the group's name '{core_en}'" if named else
                          "; its name does not carry the group's name, so it is only filed there"))
                extra = {"orpha": f"ORPHA:{code}", "member_of": row["orpha"], "member_named": named}
                add(en_name, "orphanet_member_en", via + "; Orphanet en preferred term", **extra)
                if zh_name:
                    add(zh_name, "orphanet_member_zh", via + f"; Orphanet zh preferred term ({ZH_DATASET})", **extra)
                counts["member_codes"] += 1
        if q:
            for text, lang, qpart in q.get("derived") or []:
                add(text, f"official_subtype_{lang}", f"{src} (the published name with its qualifier applied: "
                                                      f"{qpart})", qualifier_part=qpart)
                counts["derived_subtype_aliases"] += 1
        for (ck, text, why) in CORRECTED:
            if ck == key:
                add(text, "official_en_corrected", f"{src} (published English name '{e['name_en']}' {why}; "
                                                   "spelling corrected by hand)")
                counts["corrected_spellings"] += 1
        for (fk, text, source, qpart) in FOLK:
            if fk == key:
                add(text, "folk_name", source, maps_to=e["name_zh"], qualifier_part=qpart, orpha=FOLK_ORPHA.get(text))
        for a in aliases:
            if is_acronym(a["text"]):
                a["acronym"] = True
        row["aliases"] = aliases
        out_entries.append(row)
    for row in out_entries:
        for a in row["aliases"]:
            counts["aliases"] += 1
            counts["folk_names"] += a["kind"] == "folk_name"
            counts["subtype_aliases"] += a["kind"].startswith("orphanet_subtype")
            counts["member_aliases"] += a["kind"].startswith("orphanet_member")
            counts["acronyms"] += bool(a.get("acronym"))
            counts["wider_aliases"] += a.get("scope") == "wider"
    unmatched = [{"list": r["list"], "no": r["no"], "name_zh": r["name_zh"], "name_en": r["name_en"]}
                 for r in out_entries if not r["orpha"]]
    groups: Dict[str, int] = {}
    for r in out_entries:
        if r["orpha"]:
            groups[r["orpha_disorder_group"] or "?"] = groups.get(r["orpha_disorder_group"] or "?", 0) + 1
    built = now()
    prov = {
        "retrieved_at": built,
        "builder": "tools/china_list/build_china_aliases.py (run with python3 -I; --cache keeps every Orphadata "
                   "response, so a rebuild is reproducible)",
        "sources": {
            "china_lists": {"url": "bundled: zebra/data/china_rare_diseases.json",
                            "retrieved_at": lists["provenance"].get("retrieved_at"), "rows": len(entries),
                            "licence": "official PRC government publications (国卫医发〔2018〕10号, 国卫医政发〔2023〕26号); "
                                       "see that file's provenance",
                            "note": "names used exactly as published; this alias layer never alters them"},
            "orphanet_en_bulk": _src(f, "/rd-cross-referencing/orphacodes?lang=en", len(en_bulk),
                                     "ORPHAcode + English 'Preferred term' only; no synonyms, hence the per-code detail calls."),
            "orphanet_zh_bulk": _src(f, "/rd-cross-referencing/orphacodes?lang=zh", len(zh_bulk),
                                     f"Orphanet's Chinese dataset is dated {ZH_DATASET}, i.e. older than the 第二批 list "
                                     "(2023); a Chinese preferred term can be shared by several codes."),
            "orphanet_detail": {"url": API + "/rd-cross-referencing/orphacodes/{code}?lang={en|zh}",
                                "rows": sum(1 for u in f.log if "/rd-cross-referencing/orphacodes/" in u and "lang=" in u
                                            and f.log[u]["status"] == 200),
                                "licence": LICENCE,
                                "note": "per-ORPHAcode 'Preferred term', 'Synonym' and 'DisorderGroup' for linked codes, "
                                        "rule c/d candidates, subtypes and group members"},
            "orphanet_classification": {"url": API + "/rd-classification/orphacodes/{code}/hchids",
                                        "rows": sum(1 for u in f.log if "/rd-classification/" in u),
                                        "licence": LICENCE,
                                        "note": "children of each group entry's code (all Orphanet hierarchies merged), "
                                                f"walked to depth {MEMBER_DEPTH}, at most {MEMBER_CAP} codes per group; "
                                                "members carry their bulk preferred terms (en, zh), not their synonyms"},
        },
        "counts": counts,
        "linked_by_disorder_group": groups,
        "rules": [
            "a. exact: normalised name_en == an English Orphanet preferred term",
            "b. exact: a '/'- or bracket-separated part of name_en == an English Orphanet preferred term",
            f"c. exact: normalised name_zh (or one of its '/'/bracket parts) == a Chinese Orphanet preferred term (zh dataset {ZH_DATASET})",
            "c-tiebreak. when rule c hits several ORPHAcodes sharing one Chinese preferred term, the tie is broken only by an exact English preferred term or synonym among them; otherwise the entry is left null",
            "d. exact: name_en (or one of its '/'/bracket parts) == an English Orphanet SYNONYM, checked only for the RULE_D_CANDIDATES (4) best difflib candidates of entries no earlier rule linked (candidates give a reason to fetch; the match itself is exact)",
            "no link is ever made on fuzzy/approximate similarity: an entry no exact rule links keeps orpha=null and is listed in provenance.unmatched_list_entries",
            "qualified entries (1/35, 1/87, 2/78): the unqualified head and the aliases of a wider Orphanet link carry scope='wider'; subtypes inside the qualifier come from Orphanet records whose name or synonym matches both the disease and the qualifier (kind orphanet_subtype_*, with their own orpha)",
            "group entries (Orphanet 'Group of disorders'): Orphanet descendants are aliases of kind orphanet_member_* with their own orpha and member_of; member_named=true only when the member's own name carries the group's name (group_core), otherwise it is only filed under the group and is offered, never claimed",
            "qualified entries also get the subtype names formed from the published name by applying its qualifier (official_subtype_*), and three evident misspellings of published English names get a corrected alias (official_en_corrected)",
            "obsolete Orphanet codes are never linked; Orphanet's 'NON RARE IN EUROPE:' label is removed from alias text (kept as orphanet_label)",
            "acronym=true marks short all-capital aliases; an acronym alone is never 'on the list'",
            "folk names are hand-entered, carry kind=folk_name + maps_to + a source saying so, and map to a published name_zh",
        ],
        "notes": [
            "counts.aliases counts every alias object in the file, folk names, subtypes and group members included.",
            "orpha_disorder_group is Orphanet's own classification: 'Disorder' / 'Group of disorders' / 'Subtype of disorder'. A 'Group of disorders' code carries no genes in Orphanet, which matters downstream.",
            "dropped bracket/slash fragments (qualifiers and word-internal glosses, reviewed by hand): "
            + "; ".join(f"{k[0]}/{k[1]} '{k[2]}': {v}" for k, v in DROPPED_PARTS.items()),
            "folk name 玻璃人 dropped (third pass): it is used for both haemophilia and osteogenesis imperfecta",
        ] + notes + ([f"{len(missing_zh)} linked codes have no Chinese record (Orphanet zh endpoint 404; codes newer "
                      f"than the {ZH_DATASET} Chinese dataset): " + "; ".join(missing_zh)] if missing_zh else [])
          + ([f"group member lists capped: {'; '.join(capped_groups)}"] if capped_groups else []),
        "requests_made_this_run": f.requests,
        "unmatched_list_entries": unmatched,
    }
    data = {"schema": "zebra.china_aliases/2",
            "title": "中国国家罕见病目录 别名层 (alias layer for China's national rare disease lists: official names, Orphanet terms, subtypes, group members, folk names)",
            "provenance": prov, "entries": out_entries}
    if compare is not None:
        prov["rebuild_check"] = compare_with(compare, data)
    return data


def _add_subtypes(f: Fetcher, row: Dict[str, Any], q: Dict[str, Any], en_bulk: Dict[int, str], add) -> None:
    head = re.compile(q["subtype_head"], re.I)
    qual = re.compile(q["subtype_qualifier"], re.I)
    cands: List[Tuple[int, Optional[int]]] = []
    if q.get("from_classification_of"):
        root = q["from_classification_of"]
        for c in children(f, root):
            cands.append((c, root))
    if q.get("candidates_from_terms"):
        pat = re.compile(q["candidates_from_terms"], re.I)
        for code, term in sorted(en_bulk.items()):
            if pat.search(term):
                cands.append((code, None))
    f.prefetch([(f"/rd-cross-referencing/orphacodes/{c}", {"lang": "en"}) for c, _p in cands])
    chosen: List[Tuple[int, str]] = []
    for code, _parent in cands:
        den = detail(f, code, "en")
        if not den or den["pref"].upper().startswith("OBSOLETE"):
            continue
        names = [den["pref"]] + den["syn"]
        why = next((n for n in names if head.search(n) and qual.search(n)), None)
        if why:
            chosen.append((code, why))
    # a child of a chosen subtype is inside the qualifier too (Pompe's infantile / late-onset forms)
    final: List[Tuple[int, str, Optional[int]]] = [(c, w, None) for c, w in chosen]
    for code, _w in chosen:
        for kid in children(f, code):
            if kid not in {c for c, _x, _p in final}:
                final.append((kid, f"descendant of ORPHA:{code}", code))
    f.prefetch([(f"/rd-cross-referencing/orphacodes/{c}", {"lang": lg}) for c, _w, _p in final for lg in ("en", "zh")])
    for code, why, parent in final:
        den = detail(f, code, "en")
        if not den or den["pref"].upper().startswith("OBSOLETE"):
            continue
        dzh = detail(f, code, "zh")
        via = (f"Orphanet ORPHA:{code} '{den['pref']}', inside the qualifier {q['qualifier_zh']} because "
               + (f"its name/synonym '{why}' names it" if parent is None else f"it is a {why}"))
        extra = {"orpha": f"ORPHA:{code}", "qualifier_part": q["qualifier_zh"]}
        add(den["pref"], "orphanet_subtype_en", via, **extra)
        for s in den["syn"]:
            add(s, "orphanet_subtype_en", via + "; Orphanet en synonym", **extra)
        if dzh:
            add(dzh["pref"], "orphanet_subtype_zh", via + f"; Orphanet zh preferred term ({ZH_DATASET})", **extra)
            for s in dzh["syn"]:
                add(s, "orphanet_subtype_zh", via + f"; Orphanet zh synonym ({ZH_DATASET})", **extra)


def _src(f: Fetcher, path: str, rows: int, note: str) -> Dict[str, Any]:
    url = API + path
    log = f.log.get(url) or {}
    return {"url": url, "retrieved_at": log.get("retrieved_at"), "sha256": log.get("sha256"), "rows": rows,
            "licence": LICENCE, "note": note}


def compare_with(old: Dict[str, Any], new: Dict[str, Any]) -> Dict[str, Any]:
    """How the rebuilt pass-1 layer differs from an earlier file (links, and aliases of the original kinds)."""
    original_kinds = {"official_zh", "official_en", "official_part_zh", "official_part_en", "orphanet_en",
                      "orphanet_zh", "orphanet_synonym_en", "orphanet_synonym_zh", "folk_name"}
    o = {(e["list"], e["no"]): e for e in old["entries"]}
    link_diff, alias_added, alias_removed = [], [], []
    for e in new["entries"]:
        k = (e["list"], e["no"])
        oe = o.get(k) or {}
        if oe.get("orpha") != e["orpha"]:
            link_diff.append(f"{k[0]}/{k[1]}: {oe.get('orpha')} -> {e['orpha']}")
        na = {(a["text"], a["kind"]) for a in e["aliases"] if a["kind"] in original_kinds}
        oa = {(a["text"], a["kind"]) for a in oe.get("aliases") or []}
        alias_added += [f"{k[0]}/{k[1]} {t} [{kd}]" for t, kd in sorted(na - oa)]
        alias_removed += [f"{k[0]}/{k[1]} {t} [{kd}]" for t, kd in sorted(oa - na)]
    return {"compared_with": (old.get("provenance") or {}).get("retrieved_at"),
            "links_differing": link_diff, "aliases_added": alias_added, "aliases_removed": alias_removed,
            "summary": f"{len(link_diff)} link(s) differ, {len(alias_added)} alias(es) of the original kinds added, "
                       f"{len(alias_removed)} removed"}


def main(argv: Sequence[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cache", required=True, help="directory for the Orphadata responses")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--compare", help="an earlier china_disease_aliases.json to diff the rebuild against")
    a = ap.parse_args(argv)
    old = None
    if a.compare:
        with open(a.compare, encoding="utf-8") as fh:
            old = json.load(fh)
    f = Fetcher(a.cache)
    t0 = time.time()
    data = build(f, old)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=1)
        fh.write("\n")
    print(f"{a.out}: {json.dumps(data['provenance']['counts'])}; {f.requests} request(s); {time.time() - t0:.0f} s")
    if old is not None:
        print(json.dumps(data["provenance"]["rebuild_check"], ensure_ascii=False, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
