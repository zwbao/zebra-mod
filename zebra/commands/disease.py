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


def resolve(query: str) -> Tuple[Resolver, str, List[Dict[str, Any]]]:
    """Returns (resolver, status, candidates). status: resolved | ambiguous | not_found."""
    from zebra.sources import monarch, orphanet

    kind, value = parse_query(query)
    r = Resolver()
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
        got = r._try("Orphanet name search (zh)", lambda: orphanet.by_name(name, lang="zh"))
        if got and orphanet.name_matches(got.result, name):
            ok = r.from_orpha(got.result["id"], f"Chinese name '{got.result['name']}' in Orphanet (zh, 2020 dataset)")
            return r, ("resolved" if ok else "not_found"), []
        from zebra.commands import china

        try:
            hit = china.lookup([name])
            if hit["on_list"]:
                en = hit["matches"][0]["name_en"]
                r.notes.append(f"Chinese name matched the national list entry '{hit['matches'][0]['name_zh']}' / '{en}'")
                names_to_try = [p.strip() for p in re.split(r"[/()（）]", en) if len(p.strip()) > 3] + [en]
        except china.ListUnavailable:
            pass
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
        got = r._try("Monarch search", lambda: monarch.search(name, category="biolink:Disease", limit=8))
        cands = [{"id": h["id"], "name": h["name"], "source": "Monarch search"}
                 for h in (got.result["hits"] if got else []) if not h.get("deprecated")]
        close = r._try("Orphanet name search", lambda: orphanet.by_name(name))
        if close and close.result and close.result["id"] not in {c["id"] for c in cands}:
            cands.append({"id": close.result["id"], "name": close.result["name"], "source": "Orphanet closest name"})
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


def build_card(r: Resolver) -> Dict[str, Any]:
    from zebra.commands import china
    from zebra.sources import genereviews, hpo, monarch, orphanet

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

    # China lists
    names = [row.get("name"), (zh or {}).get("name"), r.name] + list(row.get("synonyms") or [])[:5] \
        + list((zh or {}).get("synonyms") or [])[:3] + ([me.get("name")] + list(me.get("synonyms") or [])[:5] if me else [])
    try:
        cn = china.lookup([n for n in dict.fromkeys(names) if n and (len(n) >= 4 or _CJK.search(n))])
        cn_data = china.load()
        r.out.sources.extend(china.source_for(cn_data, cn["matches"] if cn["on_list"] else []))
    except china.ListUnavailable as err:
        cn = {"on_list": None, "matches": [], "note": str(err)}
        r.warnings.append(f"China rare disease list unavailable: {err}")

    nh = res("nh") or {}
    gr = res("gr") or {}
    return {
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
    if not card["genereviews"]:
        lines.append("GeneReviews: no chapter found for these OMIM ids / this name")
    ha = card["hpo_annotations"]
    if ha["top"]:
        lines.append(f"HPO ({ha['disease']}, {ha['count']} terms; most frequent): " + "; ".join(
            f"{p['label']} {p['id']}" + (f" [{p['frequency']}]" if p.get("frequency") else "") for p in ha["top"][:12]))
    cn = card["china_rare_list"]
    if cn.get("on_list"):
        lines.append("China rare disease list: " + "; ".join(f"{m['list_name']} #{m['no']} {m['name_zh']} / {m['name_en']} [{m['match']}]"
                                                           for m in cn["matches"]))
    elif cn.get("on_list") is False:
        near = "; ".join(f"{m['list_name']} #{m['no']} {m['name_zh']} [{m['match']}]" for m in cn["matches"][:2])
        lines.append("China rare disease list: not listed by name" + (f" (closest, verify: {near})" if near else ""))
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
