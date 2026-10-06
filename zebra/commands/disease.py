"""zebra disease: one rare disease, as a card with every claim tied to its source.

Input: a name (English or Chinese) or an id (ORPHA:33069, OMIM:607208, MONDO:0100135).
Resolution: ids are followed through exact mappings only (Orphanet cross-references
marked E, Monarch SSSOM exactMatch, obsolete MONDO classes to their OLS replacement).
A name resolves only when a source names the disease exactly (Orphanet preferred
term or synonym, Monarch exact match); otherwise the candidates are listed with ids.

Card: ids across ORPHA/OMIM/MONDO/ICD-10/ICD-11/UMLS/MeSH/GARD, English and Chinese
names (Orphanet), definition, prevalence (class, value, region), inheritance and age
of onset (Orphanet natural history; HPO annotations), genes (Orphanet association
type; Monarch predicate and source), GeneReviews chapter, top HPO annotations by
frequency, and whether it is on China's national rare disease lists.
"""

from __future__ import annotations

import argparse
import re
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Dict, List, Optional, Tuple

from zebra.core import Outcome, UsageError, attempt

ID_RE = re.compile(r"^\s*(ORPHA|ORPHANET|OMIM|MIM|MONDO)\s*[:_ ]\s*(\d+)\s*$", re.I)
_CJK = re.compile(r"[㐀-鿿]")
# A Chinese name resolves to an Orphanet code only at or above this character-bigram
# similarity (or exactly): review E-2 (糖尿病 ~ 枫糖尿病 0.8 is maple syrup urine disease).
ZH_RESOLVE_MIN = 0.85
# ... and, above it, only when the two names differ by a closing 症/病/综合征 or punctuation:
# 0.85+ still joined 甲基丙二酸血症 to CMAMMA, 酪氨酸血症Ⅲ型 to type II (equal bigram sets),
# 地中海贫血 to alpha-thalassaemia and 神经纤维瘤病 to Neurofibroma (adversarial review P0-3)
_ZH_TAIL = re.compile(r"(综合征|综合症|症|(?<!瘤)病)$")  # 神经纤维瘤病 (-tosis) is not 神经纤维瘤 (a tumour)


def zh_near_exact(a: str, b: str) -> bool:
    """Equal once a closing 症/病/综合征 and all punctuation and spaces are dropped from both."""
    from zebra.commands import china

    def core(x: str) -> str:
        t = china.norm(x)
        return _ZH_TAIL.sub("", t)

    return bool(core(a)) and core(a) == core(b)
FREQ_ORDER = {"obligate": 1.0, "very frequent": 0.9, "frequent": 0.55, "occasional": 0.17, "very rare": 0.03,
              "excluded": 0.0}


def parse_query(text: str) -> Tuple[str, str]:
    m = ID_RE.match(text)
    if not m:
        return "name", text.strip()
    prefix = m.group(1).upper()
    num = m.group(2)
    if prefix in ("ORPHA", "ORPHANET"):
        return "ORPHA", f"ORPHA:{num}"
    if prefix in ("OMIM", "MIM"):
        return "OMIM", f"OMIM:{num}"
    return "MONDO", f"MONDO:{num.zfill(7)}"


class Resolver:
    """Collects ids for one disease; every link it makes is recorded with how it was made."""

    def __init__(self) -> None:
        self.out = Outcome(None)
        self.warnings: List[str] = []
        self.orpha: Optional[str] = None
        self.orpha_row: Optional[Dict[str, Any]] = None
        self.related_orpha: List[Dict[str, Any]] = []
        self.omim: Dict[str, str] = {}  # OMIM id -> relation/how
        self.mondo: Dict[str, str] = {}  # MONDO id -> how
        self.mondo_entity: Dict[str, Dict[str, Any]] = {}
        self.obsolete_mondo: Dict[str, Optional[str]] = {}  # given obsolete MONDO -> replacement
        self.name: Optional[str] = None
        self.notes: List[str] = []
        self.query: Optional[str] = None  # what the user typed, for the China list lookup
        self.zh_candidates: List[Dict[str, Any]] = []  # Chinese near-misses, listed when nothing resolves

    # -- helpers
    def _try(self, label: str, fn: Callable[[], Outcome]) -> Optional[Outcome]:
        got = attempt(label, fn, self.warnings)
        if got is not None:
            self.out.add(got)
        return got

    def set_orpha(self, row: Dict[str, Any], how: str) -> None:
        self.orpha, self.orpha_row = row["id"], row
        self.name = self.name or row.get("name")
        for r in row.get("refs", []):
            if r["source"] == "OMIM":
                oid = f"OMIM:{r['id']}"
                rel = {"E": "exact (Orphanet E)", "NTBT": "Orphanet: ORPHA narrower than OMIM (NTBT)",
                       "BTNT": "Orphanet: ORPHA broader than OMIM (BTNT)"}.get(r["relation"], f"Orphanet relation {r['relation']}")
                self.omim.setdefault(oid, rel)
        self.notes.append(f"ORPHA: {how}")

    def map_to_mondo(self, ids: List[str]) -> None:
        """Monarch exactMatch for OMIM/ORPHA ids; obsolete MONDO classes followed to their replacement (OLS)."""
        from zebra.sources import monarch, ols

        ids = [i for i in ids if i]
        if not ids:
            return
        got = self._try("Monarch mappings", lambda: monarch.mappings(object_ids=ids))
        if not got:
            return
        rows = [m for m in got.result if (m["subject"] or "").startswith("MONDO:")]
        rows.sort(key=lambda m: ids.index(m["object"]) if m["object"] in ids else len(ids))  # input order first
        obsolete = sorted({m["subject"] for m in rows if (m.get("subject_label") or "").lower().startswith("obsolete")})
        replaced: Dict[str, Optional[str]] = {}
        if obsolete:
            rep = self._try("OLS (obsolete MONDO)", lambda: ols.resolve_obsolete(obsolete))
            replaced = rep.result if rep else {}
        for m in rows:
            mid, obj = m["subject"], m["object"]
            if mid in obsolete:
                new = replaced.get(mid)
                if new:
                    self.mondo.setdefault(new, f"{obj} exactMatch {mid} (obsolete), replaced by {new} (Monarch + OLS)")
                else:
                    self.notes.append(f"{obj} maps to obsolete {mid} with no replacement in OLS")
            else:
                self.mondo.setdefault(mid, f"exactMatch of {obj} (Monarch)")

    # -- entry points
    def from_orpha(self, orpha_id: str, how: str = "given") -> bool:
        from zebra.sources import orphanet

        got = self._try("Orphanet", lambda: orphanet.disorder(orpha_id))
        if not got or not got.result:
            return False
        self.set_orpha(got.result, how)
        exact_omim = [k for k, v in self.omim.items() if v.startswith("exact")]
        self.map_to_mondo([orpha_id] + exact_omim)
        for ref in orphanet.exact_refs(got.result, "MONDO"):
            mid = f"MONDO:{ref}"
            if mid not in self.mondo and not any(mid in h for h in self.mondo.values()):
                self.mondo.setdefault(mid, "exact per Orphanet (not checked for obsolescence)")
        return True

    def from_omim(self, omim_id: str) -> bool:
        from zebra.sources import orphanet

        self.omim[omim_id] = "given"
        got = self._try("Orphanet (by OMIM)", lambda: orphanet.by_omim(omim_id))
        rows = got.result if got else []
        exact = [r for r in rows if r.get("omim_relation") == "E"]
        self.related_orpha = [{"id": r["id"], "name": r["name"], "relation": r.get("omim_relation")} for r in rows
                              if r.get("omim_relation") != "E"]
        if exact:
            if len(exact) > 1:
                self.notes.append("Orphanet maps this OMIM entry exactly to several ORPHAcodes: "
                                  + ", ".join(r["id"] for r in exact))
            self.set_orpha(exact[0], f"exact cross-reference of {omim_id} (Orphanet)")
        self.omim[omim_id] = "given"
        self.map_to_mondo([omim_id] + ([self.orpha] if self.orpha else []))
        if not self.name:
            from zebra.sources import hpo

            ann = self._try("HPO annotations (name)", lambda: hpo.disease_annotations(omim_id))
            if ann and ann.result.get("disease", {}).get("name"):
                self.name = ann.result["disease"]["name"]
        return bool(self.orpha or self.mondo or self.name)

    def from_mondo(self, mondo_id: str) -> bool:
        from zebra.sources import monarch, ols, orphanet

        ent = self._try("Monarch entity", lambda: monarch.entity(mondo_id))
        if ent is None or ent.result is None:
            return False
        e = ent.result
        if e.get("deprecated"):
            rep = self._try("OLS (obsolete MONDO)", lambda: ols.resolve_obsolete([mondo_id]))
            new = (rep.result if rep else {}).get(mondo_id)
            self.notes.append(f"{mondo_id} is obsolete" + (f"; replaced by {new}" if new else ""))
            if new and new != mondo_id:
                ent2 = self._try("Monarch entity", lambda: monarch.entity(new))
                if ent2 and ent2.result:
                    self.obsolete_mondo[mondo_id] = new
                    self.mondo[new] = f"replacement of obsolete {mondo_id} (OLS)"
                    mondo_id, e = new, ent2.result
        self.mondo.setdefault(mondo_id, "given")
        self.mondo_entity[mondo_id] = e
        self.name = e.get("name")
        maps = self._try("Monarch mappings", lambda: monarch.mappings(subject_ids=[mondo_id]))
        linked = [m["object"] for m in (maps.result if maps else [])] + e.get("omim", []) + e.get("orpha", [])
        for x in dict.fromkeys(linked):
            if x.startswith("OMIM:"):
                self.omim.setdefault(x, f"exactMatch/xref of {mondo_id} (Monarch)")
        orpha = next((x for x in linked if x.startswith("ORPHA:")), None)
        if orpha:
            self.from_orpha(orpha, f"exactMatch/xref of {mondo_id} (Monarch)")
        elif e.get("name"):
            got = self._try("Orphanet name search", lambda: orphanet.by_name(e["name"]))
            if got and orphanet.name_matches(got.result, e["name"]):
                self.set_orpha(got.result, f"same name as {mondo_id} in Orphanet ('{got.result['name']}'); "
                                           "no exact id mapping")
                self.map_to_mondo([k for k, v in self.omim.items() if v.startswith("exact")])
        return True


def resolve_chinese(r: Resolver, name: str) -> Tuple[Resolver, str, List[str]]:
    """Resolve a Chinese disease name (P1f, E-2). Returns (resolver, status, English names to try next).

    A Chinese name resolves only on an exact match (review E-2): Orphanet's own
    Chinese preferred term (bundled as `zebra/data/orphanet_zh_names.json`), a
    near-exact one (character-bigram Dice >= ZH_RESOLVE_MIN, e.g. 杜氏肌营养不良 ~
    杜氏肌营养不良症 0.923), or an exact alias of a national-list entry (official
    name, Orphanet term, Orphanet subtype or group member, folk name). Anything
    weaker is a candidate the caller lists with its id: 糖尿病 against 枫糖尿病
    (0.8) is maple syrup urine disease, not diabetes, and 白内障 against 蔚蓝白内障
    (0.667) is one rare cataract, not cataract. Orphadata's live closest-name
    endpoint is never trusted for Chinese (asked for 渐冻症 it answers ORPHA:90280
    冻疮样狼疮, checked 2026-10-06).

    A list alias that carries its own ORPHAcode (庞贝病 -> ORPHA:365 Pompe disease
    inside 糖原累积病（I型、Ⅱ型）; Duchenne muscular dystrophy -> ORPHA:98896 inside
    进行性肌营养不良) resolves to that disease, not to the list entry's group code.
    A list entry with no ORPHAcode hands its official English name on.
    `r.zh_candidates` keeps what was close but not accepted.
    """
    from zebra.commands import china

    zh_hits: List[Dict[str, Any]] = []
    try:
        zh_hits = china.zh_name_hits(name, limit=3)
    except (OSError, ValueError):
        pass
    english: List[str] = []
    try:
        hit = china.lookup([name])
    except china.ListUnavailable as err:
        r.warnings.append(f"China rare disease list unavailable: {err}")
        hit = {"on_list": False, "status": "unavailable", "matches": []}
    alias_hit = hit["matches"][0] if hit["on_list"] else None
    near = [z for z in zh_hits if z["match"] != "exact" and z["similarity"] >= ZH_RESOLVE_MIN
            and zh_near_exact(name, z["name_zh"])]
    if len(near) > 1:
        near = []  # two equally close Chinese names: a tie is not a resolution
    r.zh_candidates = [{"id": z["orpha"], "name": z["name_zh"], "source": f"{z['source']}, similarity {z['similarity']}"}
                       for z in zh_hits if z["match"] != "exact" and z not in near]
    if hit.get("status") in ("qualified", "possible"):
        for m in hit["matches"][:3]:
            r.zh_candidates.append({"id": m.get("alias_orpha") or m.get("orpha") or f"China list {m['list']}#{m['no']}",
                                    "name": f"{m['name_zh']} / {m['name_en']}",
                                    "source": f"national rare disease list ({m['match']} match on "
                                              f"'{m.get('matched_on')}')"})

    # Candidates in order of how much the match is worth, not in source order.
    cands: List[Tuple[int, Dict[str, Any]]] = []
    exact_zh = [z for z in zh_hits if z["match"] == "exact"]
    for z in exact_zh[:1] if len(exact_zh) == 1 else []:
        cands.append((0, {"kind": "zh", "z": z}))
    if len(exact_zh) > 1:
        r.notes.append("several Orphanet codes share this exact Chinese name: "
                       + ", ".join(f"{z['orpha']} {z['name_zh']}" for z in exact_zh))
        r.zh_candidates = [{"id": z["orpha"], "name": z["name_zh"], "source": z["source"]} for z in exact_zh] \
            + r.zh_candidates
    for z in near:
        cands.append((2, {"kind": "zh", "z": z}))
    alias_orpha = (alias_hit or {}).get("alias_orpha") or (alias_hit or {}).get("orpha")
    if alias_hit and alias_orpha:
        specific = bool(alias_hit.get("alias_orpha")) or alias_hit.get("orpha_disorder_group") == "Disorder"
        cands.append((1 if specific else 3, {"kind": "alias", "m": alias_hit, "orpha": alias_orpha}))
    cands.sort(key=lambda c: c[0])

    for _tier, c in cands:
        if c["kind"] == "zh":
            z = c["z"]
            how = (f"Chinese name matched Orphanet's Chinese preferred term '{z['name_zh']}' "
                   f"({z['orpha']}, {z['match']} match, character-bigram similarity {z['similarity']}; "
                   f"{z['source']})")
            if not r.from_orpha(z["orpha"], how):
                continue
            if z["match"] != "exact":
                r.warnings.append(f"'{name}' is not an exact Orphanet Chinese name: resolved to {z['orpha']} "
                                  f"'{z['name_zh']}' by character-bigram similarity {z['similarity']} — confirm "
                                  "it is the disease meant"
                                  + ("; next closest: " + ", ".join(
                                      f"{o['orpha']} {o['name_zh']} ({o['similarity']})"
                                      for o in zh_hits if o is not z)
                                     if len(zh_hits) > 1 else ""))
            return r, "resolved", []
        m, orpha = c["m"], c["orpha"]
        r.notes.append(f"Chinese name matched the national list entry '{m['name_zh']}' / '{m['name_en']}' "
                       f"({m['list_name']} #{m['no']}, matched '{m.get('matched_on')}' "
                       f"[{m.get('alias_kind') or m['match']}])")
        if m.get("alias_orpha"):
            how = (f"your text matched the alias '{m.get('matched_on')}' [{m.get('alias_kind')}] of "
                   f"{m['list_name']} #{m['no']} '{m['name_zh']}', which names {orpha} itself "
                   f"({m.get('alias_source')})")
        else:
            how = (f"the national list entry {m['list_name']} #{m['no']} '{m['name_zh']}' links to {orpha} "
                   f"({m.get('orpha_match')}); your text matched its alias '{m.get('matched_on')}' "
                   f"[{m.get('alias_kind')}]")
        if not r.from_orpha(orpha, how):
            continue
        if not m.get("alias_orpha") and m.get("orpha_disorder_group") and m["orpha_disorder_group"] != "Disorder":
            r.warnings.append(f"{orpha} is an Orphanet '{m['orpha_disorder_group']}', not a single "
                              "clinical entity: Orphanet assigns genes to entities, so the gene list below "
                              "may come only from Monarch. Name the specific subtype if you know it")
        return r, "resolved", []

    if alias_hit:
        m = alias_hit
        r.notes.append(f"Chinese name matched the national list entry '{m['name_zh']}' / '{m['name_en']}' "
                       f"({m['list_name']} #{m['no']}, matched '{m.get('matched_on')}' "
                       f"[{m.get('alias_kind') or m['match']}]); that entry has no single ORPHAcode, so the "
                       "English name is resolved instead")
        en = m["name_en"]
        try:  # a corrected spelling of a misspelt published name (Methylmalonic "Academia") goes first
            al = china.load_aliases()["_by_key"].get((m["list"], m["no"])) or {}
            fixed = [a["text"] for a in al.get("aliases") or [] if a.get("kind") == "official_en_corrected"]
        except (OSError, ValueError, KeyError):
            fixed = []
        english = fixed + [p.strip() for p in re.split(r"[/()（）]", en) if len(p.strip()) > 3] + [en]
    if zh_hits:
        r.notes.append("Orphanet Chinese names close to this text, none accepted (a Chinese name resolves only on an "
                       f"exact or near-exact match, similarity >= {ZH_RESOLVE_MIN}): "
                       + ", ".join(f"{z['orpha']} {z['name_zh']} ({z['similarity']})" for z in zh_hits))
    elif not english:
        r.warnings.append(f"'{name}' matched no Orphanet Chinese preferred term (bundled index, "
                          "Orphanet zh dataset 2020-06-01) and no national-list name or alias. Chinese coverage is "
                          "the Orphanet zh dataset plus the list alias layer; try the English name, or "
                          "`zebra china <name>` to see the closest list entries")
    return r, "not_found", english


def resolve(query: str) -> Tuple[Resolver, str, List[Dict[str, Any]]]:
    """Returns (resolver, status, candidates). status: resolved | ambiguous | not_found."""
    from zebra.sources import monarch, orphanet

    kind, value = parse_query(query)
    r = Resolver()
    # The card's China-list line must agree with `zebra china <same text>`: the
    # list is matched on names, and the name the family used is one of them.
    r.query = query.strip() if kind == "name" else None
    if kind == "ORPHA":
        ok = r.from_orpha(value)
        if not ok:  # Orphanet down or silent: Monarch may still know the code through Mondo
            r.map_to_mondo([value])
            current = [m for m in r.mondo]
            if current:
                r.notes.append(f"Orphanet did not return {value}; card built from Monarch ({current[0]})")
                ok = r.from_mondo(current[0])
        return r, ("resolved" if ok else "not_found"), []
    if kind == "OMIM":
        ok = r.from_omim(value)
        return r, ("resolved" if ok else "not_found"), []
    if kind == "MONDO":
        ok = r.from_mondo(value)
        return r, ("resolved" if ok else "not_found"), []

    name = value
    names_to_try = [name]
    if _CJK.search(name):
        r2, status2, cands2 = resolve_chinese(r, name)
        if status2 == "resolved":
            return r2, status2, cands2
        if not cands2:
            # E-2: no English name to follow; the Chinese near-misses are candidates, never a card
            return r, ("ambiguous" if r.zh_candidates else "not_found"), list(r.zh_candidates)
        names_to_try = cands2  # the list entry's official English names
    for n in names_to_try:
        got = r._try("Orphanet name search", lambda n=n: orphanet.by_name(n))
        if got and orphanet.name_matches(got.result, n):
            ok = r.from_orpha(got.result["id"], f"Orphanet preferred term or synonym equals '{n}'")
            return r, ("resolved" if ok else "not_found"), []
    exact_hits: List[Dict[str, Any]] = []
    for n in names_to_try:
        got = r._try("Monarch search", lambda n=n: monarch.search(n, category="biolink:Disease", limit=10, match_type="exact"))
        exact_hits = [h for h in (got.result["hits"] if got else []) if not h.get("deprecated")]
        if exact_hits:
            break
    if len(exact_hits) == 1:
        ok = r.from_mondo(exact_hits[0]["id"])
        return r, ("resolved" if ok else "not_found"), []
    cands: List[Dict[str, Any]] = []
    if exact_hits:
        cands = [{"id": h["id"], "name": h["name"], "source": "Monarch (exact name)"} for h in exact_hits]
    else:
        # a Chinese name is searched by the list entry's English name, never as Chinese text
        asked = names_to_try[-1] if _CJK.search(name) else name
        got = r._try("Monarch search", lambda: monarch.search(asked, category="biolink:Disease", limit=8))
        cands = [{"id": h["id"], "name": h["name"], "source": "Monarch search"}
                 for h in (got.result["hits"] if got else []) if not h.get("deprecated")]
        close = r._try("Orphanet name search", lambda: orphanet.by_name(asked))
        if close and close.result and close.result["id"] not in {c["id"] for c in cands}:
            cands.append({"id": close.result["id"], "name": close.result["name"], "source": "Orphanet closest name"})
    cands += [c for c in r.zh_candidates if c["id"] not in {x["id"] for x in cands}]
    return r, ("ambiguous" if cands else "not_found"), cands


# ---------------------------------------------------------------- card


def _freq_weight(f: Optional[str]) -> float:
    if not f:
        return 0.5
    s = str(f).strip().lower()
    m = re.match(r"^(\d+)\s*/\s*(\d+)$", s)
    if m and int(m.group(2)):
        return int(m.group(1)) / int(m.group(2))
    m = re.match(r"^([\d.]+)\s*%$", s)
    if m:
        return float(m.group(1)) / 100
    for k, v in FREQ_ORDER.items():
        if s.startswith(k):
            return v
    return 0.5


ORPHADATA_DATASETS = ("rd-cross-referencing", "rd-epidemiology", "rd-natural_history", "rd-associated-genes",
                      "rd-phenotypes", "rd-classification", "rd-medical-specialties")


def support_pointers(orpha: Optional[str], row: Dict[str, Any], ids: Dict[str, Any],
                     mp: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Where to find patient organisations and expert centres (F39).

    These are pointers to pages, not retrieved records, and the note says so.
    Orphanet's public API serves no patient-organisation or expert-centre
    dataset: its OpenAPI document (https://api.orphadata.com/openapi.json,
    read 2026-10-06) lists only the datasets in `ORPHADATA_DATASETS`. The
    orpha.net disease page does list organisations and expert centres, and
    GARD (NIH) lists patient organisations on the page for the GARD id Orphanet
    cross-references, so both are given as links to open.
    """
    pointers: List[Dict[str, Any]] = []
    if row.get("url"):
        pointers.append({"what": "patient organisations and expert centres",
                         "source": "Orphanet disease page (orpha.net)", "url": row["url"],
                         "retrieved": False})
    for g in (ids.get("GARD") or [])[:2]:
        gid = (g.get("id") if isinstance(g, dict) else str(g)).split(":")[-1]
        if gid.isdigit():
            pointers.append({"what": "patient organisations and advocacy groups",
                             "source": f"GARD (NIH) GARD:{gid}, cross-referenced by Orphanet",
                             "url": f"https://rarediseases.info.nih.gov/diseases/{gid}/index",
                             "retrieved": False})
    if mp and mp.get("url"):
        pointers.append({"what": "plain-language description and further reading",
                         "source": "MedlinePlus Genetics", "url": mp["url"], "retrieved": True})
    return {
        "pointers": pointers,
        "note": "Links to open, not data this command retrieved. Orphanet's public API "
                "(api.orphadata.com) serves no patient-organisation or expert-centre dataset — its OpenAPI "
                "document lists only " + ", ".join(ORPHADATA_DATASETS) + " — so the organisations themselves "
                "cannot be named here. China's 全国罕见病诊疗协作网 hospitals: `zebra china hospitals "
                "[--province <省>]`; approvals, 医保 and trials in China: `zebra access <disease>`.",
    }


def build_card(r: Resolver) -> Dict[str, Any]:
    from zebra.commands import china
    from zebra.sources import genereviews, hpo, medlineplus, monarch, orphanet

    orpha = r.orpha
    exact_omim = [k for k, v in r.omim.items() if v in ("given",) or v.startswith("exact")]
    mondo_ids = [m for m in r.mondo if m not in r.obsolete_mondo and "obsolete" not in m][:3]
    jobs: Dict[str, Tuple[str, Callable[[], Outcome]]] = {}
    if orpha:
        jobs["zh"] = ("Orphanet (Chinese)", lambda: orphanet.disorder(orpha, lang="zh"))
        jobs["epi"] = ("Orphanet epidemiology", lambda: orphanet.epidemiology(orpha))
        jobs["nh"] = ("Orphanet natural history", lambda: orphanet.natural_history(orpha))
        jobs["ogenes"] = ("Orphanet genes", lambda: orphanet.genes(orpha))
        jobs["hpo_orpha"] = ("HPO annotations", lambda: hpo.disease_annotations(orpha))
    if exact_omim:
        jobs["hpo_omim"] = ("HPO annotations", lambda: hpo.disease_annotations(exact_omim[0]))
    for i, mid in enumerate(mondo_ids[:2]):
        if mid not in r.mondo_entity:
            jobs[f"ment{i}"] = ("Monarch entity", lambda mid=mid: monarch.entity(mid))
        jobs[f"mgenes{i}"] = ("Monarch gene associations", lambda mid=mid: monarch.disease_genes(mid))
    name = r.name
    jobs["gr"] = ("GeneReviews", lambda: genereviews.chapters(omim_ids=exact_omim, name=name))
    # F39: a plain-language source for the family, and the chapter text above.
    orow = r.orpha_row or {}
    mondo_names = [(r.mondo_entity.get(m) or {}).get("name") for m in mondo_ids]
    mp_names = [n for n in dict.fromkeys(
        [orow.get("name"), r.name] + mondo_names + list(orow.get("synonyms") or [])[:4]) if n]
    jobs["mp"] = ("MedlinePlus Genetics", lambda: medlineplus.condition(mp_names, exact_omim))

    results: Dict[str, Optional[Outcome]] = {}
    warn: Dict[str, List[str]] = {k: [] for k in jobs}
    with ThreadPoolExecutor(max_workers=6) as pool:
        futs = {k: pool.submit(attempt, label, fn, warn[k]) for k, (label, fn) in jobs.items()}
        for k, f in futs.items():
            results[k] = f.result()
    for k in jobs:
        r.warnings.extend(warn[k])
        if results.get(k) is not None:
            r.out.add(results[k])

    def res(k: str) -> Any:
        o = results.get(k)
        return o.result if o is not None else None

    for i, mid in enumerate(mondo_ids[:2]):
        if res(f"ment{i}"):
            r.mondo_entity[mid] = res(f"ment{i}")
    row = r.orpha_row or {}
    zh = res("zh")
    me = next((r.mondo_entity[m] for m in mondo_ids if m in r.mondo_entity), None)

    # ids
    ids: Dict[str, Any] = {"ORPHA": orpha, "OMIM": [{"id": k, "relation": v} for k, v in r.omim.items()],
                           "MONDO": [{"id": k, "how": v} for k, v in r.mondo.items()]}
    for src, refs in orphanet.curie_refs(row).items() if row else []:
        if src in ("OMIM", "MONDO"):
            continue
        ids[src] = [{"id": x["id"], "relation": x["relation"]} for x in refs]
    if me and not row:
        for k, v in (me.get("other_ids") or {}).items():
            if k in ("ICD10CM", "icd11.foundation", "UMLS", "GARD", "MESH", "MEDGEN"):
                ids[k] = v[:5]
    if r.related_orpha:
        ids["ORPHA_related"] = r.related_orpha[:8]

    # genes
    genes: List[Dict[str, Any]] = []
    for g in res("ogenes") or []:
        genes.append({"symbol": g["symbol"], "hgnc": g["hgnc"], "association": g["association"], "status": g["status"],
                      "source": "Orphanet", "pmids": g["pmids"][:5]})
    seen_m: Dict[Tuple[str, str, str], Dict[str, Any]] = {}
    for i, mid in enumerate(mondo_ids[:2]):
        for g in (res(f"mgenes{i}") or {}).get("genes", []):
            key = (g["label"], g["predicate"], g["source"])
            if key in seen_m:
                if g.get("for_label") and g["for_label"] not in seen_m[key]["diseases"]:
                    seen_m[key]["diseases"].append(g["for_label"])
                continue
            grow = {"symbol": g["label"], "hgnc": g["id"], "association": g["predicate"],
                    "status": g["original_predicate"], "source": f"Monarch ({g['source']})",
                    "diseases": [g["for_label"]] if g.get("for_label") else []}
            seen_m[key] = grow
            genes.append(grow)

    # CP1-7: a card reached through an Orphanet group (脊髓性肌萎缩症 -> ORPHA:70) has no
    # OMIM id and a name GeneReviews does not use, so no chapter matched while the
    # English entry found NBK1352. The disease-causing genes are a second, exact key:
    # NCBI's own gene table maps SMN1 to NBK1352.
    gr_genes = [g["symbol"] for g in genes if g.get("symbol") and (
        g["source"] == "Orphanet" and str(g.get("association") or "").lower().startswith("disease-causing")
        or g["source"] != "Orphanet" and str(g.get("association") or "").lower() == "causes")]
    gr_genes = list(dict.fromkeys(gr_genes))[:4]
    if not (res("gr") or {}).get("chapters") and gr_genes:
        more = attempt("GeneReviews (by gene)", lambda: genereviews.chapters(genes=gr_genes, limit=3), r.warnings)
        if more is not None:
            r.out.add(more)
            results["gr"] = Outcome(dict(more.result, by_gene=gr_genes))

    # HPO annotations (phenotypes by frequency; inheritance / onset terms kept apart)
    ann_src = "hpo_orpha" if res("hpo_orpha") and res("hpo_orpha").get("phenotypes") else "hpo_omim"
    ann = res(ann_src) or {}
    phen = [p for p in ann.get("phenotypes", []) if p["category"] not in ("Inheritance", "Clinical course")]
    phen.sort(key=lambda p: -_freq_weight(p.get("frequency")))
    hpo_inh, hpo_onset = [], []
    for k in ("hpo_omim", "hpo_orpha"):
        for p in (res(k) or {}).get("phenotypes", []):
            if p["category"] == "Inheritance" and p["label"] not in hpo_inh:
                hpo_inh.append(p["label"])
            if p["category"] == "Clinical course" and "onset" in (p["label"] or "").lower() and p["label"] not in hpo_onset:
                hpo_onset.append(p["label"])

    # China lists (E-1): the ORPHAcode first (the entry's own link, an Orphanet subtype
    # inside a qualifier, or a member of a group entry), then the names. Acronym
    # synonyms are left out: CAPS is a synonym of both catastrophic antiphospholipid
    # syndrome and cryopyrin-associated periodic syndrome.
    names = [r.query, row.get("name"), (zh or {}).get("name"), r.name] + list(row.get("synonyms") or [])[:5] \
        + list((zh or {}).get("synonyms") or [])[:3] + ([me.get("name")] + list(me.get("synonyms") or [])[:5] if me else [])
    names = [n for n in dict.fromkeys(names) if n and (_CJK.search(n) or len(china.norm(n)) > 5)
             and not china.is_acronym(n)]  # 'caps' (Monarch writes some acronyms in lower case) is still CAPS
    try:
        cn = china.lookup(names, orpha_ids=[orpha] if orpha else [])
        cn_data = china.load()
        r.out.sources.extend(china.source_for(cn_data, cn["matches"] if cn["status"] in ("on_list", "qualified") else []))
        if cn["status"] == "qualified":
            for q in cn.get("qualifier") or []:
                r.warnings.append(f"China rare disease list: entry {q['list']}#{q['no']} covers {q.get('covers') or q.get('zh')}"
                                  " — this disease is on the list only if the diagnosis is that subtype")
    except china.ListUnavailable as err:
        cn = {"on_list": None, "status": "unavailable", "matches": [], "note": str(err)}
        r.warnings.append(f"China rare disease list unavailable: {err}")

    nh = res("nh") or {}
    gr = res("gr") or {}
    mp = res("mp")
    return {
        "plain_language": mp,
        "support": support_pointers(orpha, row, ids, mp),
        "name": row.get("name") or r.name or (me or {}).get("name"),
        "name_mondo": (me or {}).get("name") if me and (me.get("name") or "").lower() != (row.get("name") or r.name or "").lower() else None,
        "name_zh": (zh or {}).get("name"),
        "name_zh_source": "Orphanet zh dataset" + (f" ({(zh or {}).get('dataset_date', '')[:10]})" if zh else "") if zh else None,
        "synonyms": (list(row.get("synonyms") or []) + [s for s in ((me or {}).get("synonyms") or []) if s not in (row.get("synonyms") or [])])[:10],
        "synonyms_zh": (zh or {}).get("synonyms") or [],
        "ids": ids,
        "definition": row.get("definition") or (me or {}).get("description"),
        "definition_source": "Orphanet" if row.get("definition") else ("Monarch/MONDO" if me and me.get("description") else None),
        "typology": row.get("typology"),
        "prevalence": (res("epi") or [])[:8],
        "inheritance": nh.get("inheritance") or [],
        "inheritance_hpo": hpo_inh,
        "onset": nh.get("onset") or [],
        "onset_hpo": hpo_onset,
        "genes": genes[:25],
        "genereviews": gr.get("chapters", [])[:3],
        "genereviews_unavailable": not gr and any(w.startswith("GeneReviews") for w in r.warnings) or None,
        "hpo_annotations": {"disease": ann.get("disease", {}).get("id"), "count": len(phen),
                            "top": [{"id": p["id"], "label": p["label"], "frequency": p.get("frequency")} for p in phen[:15]]},
        "china_rare_list": cn,
        "links": {k: v for k, v in {
            "orphanet": row.get("url"),
            "monarch": monarch.web_url(mondo_ids[0]) if mondo_ids else None,
            "omim": [f"https://omim.org/entry/{k.split(':')[1]}" for k in exact_omim][:3] or None,
        }.items() if v},
        "notes": r.notes + ([f"{len(mondo_ids)} MONDO classes reached by exact mappings ({', '.join(mondo_ids)}); Mondo treats "
                             "them as distinct concepts, the sources above treat them as this disease"] if len(mondo_ids) > 1 else []),
    }


def _prev(p: Dict[str, Any]) -> str:
    if (p.get("type") or "").lower().startswith("cases"):
        core = f"{p.get('value') or '?'} cases/families reported"
    else:
        core = f"{p.get('type')} {p.get('class')}" + (f" (mean {p['value']})" if p.get("value") else "")
    return core + f", {p.get('region')}" + (", not yet validated" if p.get("validation") and p["validation"] != "Validated" else "")


def render(card: Dict[str, Any]) -> str:
    ids = card["ids"]
    head = card["name"] or "?"
    if card.get("name_zh"):
        head += f" ({card['name_zh']})"
    if card.get("name_mondo"):
        head += f"; MONDO name: {card['name_mondo']}"
    lines = [head]
    idparts = []
    if ids.get("ORPHA"):
        idparts.append(ids["ORPHA"])
    if ids.get("OMIM"):
        idparts.append("OMIM " + "; ".join(f"{o['id']} [{o['relation']}]" for o in ids["OMIM"]))
    if ids.get("MONDO"):
        idparts.append("MONDO " + "; ".join(f"{m['id']} [{m['how']}]" for m in ids["MONDO"]))
    for k in ("ICD-10", "ICD-11", "UMLS", "MeSH", "GARD", "ICD10CM", "icd11.foundation"):
        if ids.get(k):
            vals = ids[k]
            idparts.append(f"{k} " + ", ".join(v["id"].split(":", 1)[-1] + (f" ({v['relation']})" if v.get("relation") and v["relation"] != "E" else "")
                                              if isinstance(v, dict) else str(v) for v in vals))
    lines.append("ids: " + " | ".join(idparts))
    if ids.get("ORPHA_related"):
        lines.append("related ORPHA (not exact): " + ", ".join(f"{x['id']} {x['name']} [{x['relation']}]" for x in ids["ORPHA_related"]))
    if card.get("definition"):
        lines.append(f"definition ({card['definition_source']}): {card['definition']}")
    if card["prevalence"]:
        lines.append("prevalence (Orphanet): " + "; ".join(_prev(p) for p in card["prevalence"]))
    inh = card["inheritance"] or card["inheritance_hpo"]
    if inh:
        lines.append(f"inheritance: {', '.join(inh)}" + (" (Orphanet)" if card["inheritance"] else " (HPO annotation)"))
    on = card["onset"] or card["onset_hpo"]
    if on:
        lines.append(f"onset: {', '.join(on)}" + (" (Orphanet)" if card["onset"] else " (HPO annotation)"))
    og = [g for g in card["genes"] if g["source"] == "Orphanet"]
    mg = [g for g in card["genes"] if g["source"] != "Orphanet"]
    if og:
        lines.append("genes (Orphanet): " + "; ".join(f"{g['symbol']} ({g['association']}, {g['status']})" for g in og[:12]))
    if mg:
        lines.append("genes (Monarch): " + "; ".join(f"{g['symbol']} {g['association']} [{g['source']}]" for g in mg[:12]))
    for ch in card["genereviews"]:
        lines.append(f"GeneReviews: {ch['nbk']} {ch['title']}" + (f" ({ch['dates']})" if ch.get("dates") else "")
                     + f" {ch['url']}  [matched by {', '.join(ch['matched_by'])}]")
        for head, body in (ch.get("sections") or {}).items():
            lines.append(f"  {head}: {body if len(body) <= 600 else body[:600].rsplit(' ', 1)[0] + ' …'}")
        if ch.get("sections"):
            lines.append(f"  (chapter summary as published, PMID {ch.get('pmid')}; full chapter at {ch['url']})")
    if not card["genereviews"]:
        if card.get("genereviews_unavailable"):
            lines.append("GeneReviews: not checked — the source was unavailable (see warnings); this is not "
                         "evidence that no chapter exists")
        else:
            lines.append("GeneReviews: no chapter found for these OMIM ids, this name or the disease genes")
    mp = card.get("plain_language")
    if mp:
        desc = next((b["text"] for b in mp.get("text") or [] if b["role"] == "description"),
                    (mp.get("text") or [{}])[0].get("text", ""))
        lines.append(f"plain language ({mp['register']}): {mp['name']} — {mp['url']}"
                     + (f" [OMIM match {', '.join(mp['omim_overlap'])}]" if mp.get("omim_overlap")
                        else " [matched by name only]"))
        if desc:
            lines.append("  " + (desc if len(desc) <= 700 else desc[:700].rsplit(" ", 1)[0] + " …"))
    else:
        lines.append("plain language: no MedlinePlus Genetics page found by name")
    sup = card.get("support") or {}
    if sup.get("pointers"):
        lines.append("support (links to open, not retrieved records):")
        for p in sup["pointers"]:
            lines.append(f"  {p['what']} — {p['source']}: {p['url']}")
    if sup.get("note"):
        lines.append(f"  note: {sup['note']}")
    ha = card["hpo_annotations"]
    if ha["top"]:
        lines.append(f"HPO ({ha['disease']}, {ha['count']} terms; most frequent): " + "; ".join(
            f"{p['label']} {p['id']}" + (f" [{p['frequency']}]" if p.get("frequency") else "") for p in ha["top"][:12]))
    cn = card["china_rare_list"]
    if cn.get("on_list"):
        lines.append("China rare disease list: " + "; ".join(
            f"{m['list_name']} #{m['no']} {m['name_zh']} / {m['name_en']} [{m['match']}"
            + (f" via {m['matched_on']}" if m["match"] in ("subtype", "group_member") else "") + "]"
            for m in cn["matches"]))
    elif cn.get("status") == "qualified":
        lines.append("China rare disease list: only a subtype is listed — " + "; ".join(
            f"{m['list_name']} #{m['no']} {m['name_zh']} covers {(m.get('qualifier') or {}).get('covers', '?')}"
            for m in cn["matches"]) + " (on the list only if the diagnosis is that subtype)")
    elif cn.get("on_list") is False:
        near = "; ".join(f"{m['list_name']} #{m['no']} {m['name_zh']} [{m['match']}]" for m in cn["matches"][:2])
        lines.append("China rare disease list: not matched by id or name (not the same as absent)"
                     + (f" (closest, verify: {near})" if near else ""))
    else:
        lines.append("China rare disease list: not checked — the bundled list is unavailable (see warnings)")
    for n in card["notes"]:
        lines.append(f"note: {n}")
    return "\n".join(lines)


def out_warnings(r: Resolver) -> List[str]:
    return list(r.out.warnings)


def _run(args: argparse.Namespace) -> Outcome:
    query = " ".join(args.query).strip()
    if not query:
        raise UsageError("give a disease name or an id (ORPHA:33069, OMIM:607208, MONDO:0100135)")
    r, status, cands = resolve(query)
    if status != "resolved":
        if status == "ambiguous":
            text = f"'{query}' is not an exact disease name in Orphanet or Monarch; candidates (pick an id):\n" + \
                   "\n".join(f"  {c['id']}  {c['name']}  [{c['source']}]" for c in cands)
        elif any("unavailable" in w for w in r.warnings + out_warnings(r)):
            status = "unavailable"
            text = f"'{query}': could not be resolved; a source was unavailable (see warnings)"
        else:
            text = f"'{query}': no disease found in Orphanet or Monarch"
        if status == "ambiguous" and _CJK.search(query):
            r.warnings.append(f"'{query}' is not an exact Chinese disease name in Orphanet or the national lists, so "
                              "no card is given: the candidates below only resemble it and may be different, "
                              "specific rare diseases. Pick an id, or give the English name")
        out = r.out
        out.result = {"status": status, "query": query, "candidates": cands, "notes": r.notes}
        out.warnings = r.warnings + out.warnings
        out.text = text
        out.query = {"query": query}
        return out
    card = build_card(r)
    card["status"] = "resolved"
    card["query"] = query
    out = r.out
    out.result = card
    out.warnings = r.warnings + out.warnings
    out.text = render(card)
    out.query = {"query": query}
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("disease", help="disease card: ids, names (en/zh), definition, prevalence, inheritance, genes, "
                                       "GeneReviews, HPO, China rare disease lists")
    p.add_argument("query", nargs="+", help="name (English or Chinese) or ORPHA:/OMIM:/MONDO: id")
    p.set_defaults(func=_run)
