"""zebra access: approved where, covered by 医保, trialled in China, and which 协作网 hospital — in one answer.

A Chinese family asks four things after a diagnosis: is the drug approved in
China (and elsewhere), does 医保 pay for it and on what condition, is there a
trial in China, and which hospital of the national rare-disease network to go
to. This command answers them for a drug, a disease or a gene, every field with
its source:

  approvals   FDA (openFDA labels and Drugs@FDA) and EMA (the agency's medicines
              export) via zebra.sources.regulators; China from the bundled
              NMPA/CDE/MOF documents (zebra/data/china_rare_drug_approvals.json).
  医保        the bundled 国家医保药品目录 (zebra/data/china_nrdl.json): rows for
              the drug's Chinese name and rows whose payment restriction names
              the disease; the 商保创新药目录 is listed apart, because 医保 does not
              pay for it.
  trials      ClinicalTrials.gov sites in China (relevance-checked); ChiCTR cannot
              be searched by a script (see zebra.commands.trials).
  hospitals   the 全国罕见病诊疗协作网 lead hospitals, or every member in one
              province with --province.

Rules the output keeps: a drug missing from a bundled list is "not in the
bundled list" (each list states its coverage), never "not approved" or "not
covered"; a source that failed is a named warning; the national-list status
uses `zebra china`'s rules (qualified entries and acronyms are never "on the list").
"""

from __future__ import annotations

import argparse
import re
from typing import Any, Dict, List, Optional, Tuple

from zebra.core import Outcome, UsageError, attempt
from zebra.sources import china_access as ca

_CJK = re.compile(r"[㐀-鿿]")
GENE_RE = re.compile(r"^[A-Z0-9][A-Z0-9\-]{1,11}$")
NOT_APPROVED_WORDING = ("not in the bundled list — the list is not a complete register (see coverage), so this is "
                        "not evidence of anything about approval")
MAX_DRUGS = 10


def _short(text: Optional[str], n: int = 300) -> Optional[str]:
    if not text:
        return text
    t = " ".join(str(text).split())
    return t if len(t) <= n else t[:n].rsplit(" ", 1)[0] + " …"


# ------------------------------------------------------------------ classify

def _classify(query: str, kind: str, out: Outcome) -> Dict[str, Any]:
    """{"kind": drug|disease|gene|unresolved, ...what was matched and how}."""
    from zebra.commands import china

    cjk = bool(_CJK.search(query))
    info: Dict[str, Any] = {"kind": kind, "query": query, "how": "given with --as" if kind != "auto" else None}
    cn = None
    try:
        cn = china.lookup([query])
    except china.ListUnavailable as err:
        out.warnings.append(f"China rare disease list unavailable: {err}")
    drug_rows: List[Dict[str, Any]] = []
    nrdl_rows: List[Dict[str, Any]] = []
    # a Chinese name that ends like a disease (胰岛素瘤, 维生素B12缺乏) is never read as a drug
    drug_like = not (cjk and ca.DISEASE_TAIL.search(query.strip()))
    if drug_like:
        try:
            drug_rows = ca.approvals_for_drug([query])
        except ca.NotBundled as err:
            out.warnings.append(f"China drug approvals {err}: a drug name could not be recognised from it")
        if cjk:
            try:
                nrdl_rows = ca.nrdl_for_drug([query])
            except ca.NotBundled as err:
                out.warnings.append(f"国家医保药品目录 {err}: a drug name could not be recognised from it")
    if kind != "auto":
        info.update({"china": cn, "drug_rows": drug_rows, "nrdl_rows": nrdl_rows})
        if kind in ("gene", "disease", "drug") and not cjk:
            from zebra.sources import opentargets as ot

            got = attempt("Open Targets search", lambda: ot.search(query, (
                {"gene": "target"}.get(kind, kind),), 10), out.warnings)
            hits = out.add(got)["hits"] if got is not None else []
            exact = next((h for h in hits if (h["name"] or "").lower() == query.lower()), None)
            if exact is not None:
                info[{"gene": "target", "disease": "ot_disease", "drug": "ot_drug"}[kind]] = exact
            elif kind in ("gene", "disease") and not (kind == "disease" and cn and cn["status"] in ("on_list", "qualified")):
                # a near miss is never answered for: the family would read another disease's access as theirs
                return dict(info, kind="unresolved", how=f"--as {kind}: no exact Open Targets {kind} name",
                            search_hits=[{k: h[k] for k in ("id", "name", "entity")} for h in hits[:6]])
        return info
    if cn and cn["status"] in ("on_list", "qualified") and cjk:
        return dict(info, kind="disease", how=f"national rare-disease list ({cn['status']})", china=cn)
    if drug_rows or nrdl_rows:
        return dict(info, kind="drug", how="bundled China drug tables (Chinese or English name)",
                    drug_rows=drug_rows, nrdl_rows=nrdl_rows)
    if cjk:
        if cn and cn["status"] in ("possible",):
            info["candidates"] = [{"list": m["list"], "no": m["no"], "name_zh": m["name_zh"], "match": m["match"],
                                   "matched_on": m.get("matched_on")} for m in cn["matches"][:5]]
        return dict(info, kind="unresolved", how="no exact national-list name and no bundled drug name")
    from zebra.sources import opentargets as ot

    did = ot.normalise_disease_id(query)
    if did:
        return dict(info, kind="disease", how="disease id", ot_disease={"id": did, "name": None}, china=cn)
    got = attempt("Open Targets search", lambda: ot.search(query, ("disease", "target", "drug"), 10), out.warnings)
    hits = out.add(got)["hits"] if got is not None else []
    if GENE_RE.match(query):
        t = next((h for h in hits if h["entity"] == "target" and (h["name"] or "").upper() == query.upper()), None)
        if t:
            return dict(info, kind="gene", how="exact gene symbol (Open Targets)", target=t)
    dis = next((h for h in hits if h["entity"] == "disease" and (h["name"] or "").lower() == query.lower()), None)
    if dis:
        return dict(info, kind="disease", how="exact disease name (Open Targets)", ot_disease=dis, china=cn)
    drg = next((h for h in hits if h["entity"] == "drug" and (h["name"] or "").lower() == query.lower()), None)
    if drg:
        return dict(info, kind="drug", how="exact drug name (Open Targets)", ot_drug=drg)
    if cn and cn["status"] in ("on_list", "qualified"):
        return dict(info, kind="disease", how=f"national rare-disease list ({cn['status']})", china=cn)
    # no exact name anywhere: list the near misses, answer for none of them (a family would read another
    # disease's approvals and hospitals as its own)
    return dict(info, kind="unresolved", how="no exact name in Open Targets, the national lists or the bundled drug "
                                            "tables", search_hits=[{k: h[k] for k in ("id", "name", "entity")}
                                                                   for h in hits[:6]])


# ------------------------------------------------------------------ blocks

def _china_list_block(cn: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if cn is None:
        return {"status": "unavailable"}
    return {"status": cn["status"], "source": "China national rare disease lists (2018, 2023), bundled",
            "entries": [{k: m.get(k) for k in ("list", "no", "name_zh", "name_en", "match", "matched_on", "qualifier")
                         if m.get(k) is not None} for m in cn["matches"][:3]]}


def _zh_names_for(cn: Optional[Dict[str, Any]], query: str) -> List[str]:
    """Chinese names to search 医保 restrictions and NMPA indications with: the text typed and the entry's names."""
    from zebra.commands import china

    names = [query] if _CJK.search(query) else []
    if cn and cn["status"] in ("on_list", "qualified"):
        by_key = china.load_aliases()["_by_key"]
        for m in cn["matches"][:3]:
            names.append(re.split(r"[（(/]", m["name_zh"])[0])
            for a in (by_key.get((m["list"], m["no"])) or {}).get("aliases") or []:
                if a.get("kind") in ("official_zh", "official_part_zh", "orphanet_zh", "orphanet_synonym_zh",
                                     "folk_name") and _CJK.search(a.get("text") or "") and not a.get("acronym"):
                    names.append(a["text"])
    return [n for n in dict.fromkeys(n.strip() for n in names) if n and len(ca.zh_norm(n)) >= 3]


def _china_approvals(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    out = []
    for r in rows[:MAX_DRUGS]:
        ap = r.get("approval")  # present only when an official document says approved / marketed in China
        out.append({"drug_zh": r.get("drug_zh"), "inn": r.get("inn"), "brand": r.get("brand"),
                    "indication_zh": _short(r.get("indication_zh"), 200),
                    "approved_in_china": ap is not None,
                    "approval": ({k: ap.get(k) for k in ("date", "year", "type", "announced", "source_id")
                                  if ap.get(k)} if ap is not None else None),
                    "evidence": [_short(e.get("what"), 160) for e in (r.get("evidence") or [])][:3],
                    "note": r.get("note"), "matched_by": r.get("matched_by")})
    return out


def _china_reading(rows: List[Dict[str, Any]]) -> str:
    """One line on what the bundled China documents say about a drug: approved, named only, or absent."""
    if not rows:
        return "not in the bundled list (not evidence of non-approval)"
    ok = [r for r in rows if r.get("approval") is not None]
    if ok:
        ap = ok[0]["approval"]
        when = ap.get("date") or ap.get("year") or ap.get("announced")
        return (f"approved in China per the bundled documents ({ok[0].get('drug_zh') or ok[0].get('inn')}"
                + (f", {when}" if when else "") + f"; {ap.get('source_id')})")
    return (f"named in the bundled documents ({rows[0].get('drug_zh') or rows[0].get('inn')}) but none of them says it "
            "is approved in China — status unknown here")


def _nrdl_block(names_zh_drugs: List[str], names_zh_disease: List[str], out: Outcome) -> Dict[str, Any]:
    try:
        nrdl = ca.load_nrdl()
    except ca.NotBundled as err:
        out.warnings.append(f"国家医保药品目录 {err}; 医保 status not checked")
        return {"status": "not_bundled", "note": str(err)}
    rows: List[Dict[str, Any]] = []
    if names_zh_drugs:
        rows += [dict(r, matched_by=r.get("matched_by") or "Chinese drug name") for r in ca.nrdl_for_drug(names_zh_drugs)]
    if names_zh_disease:
        rows += ca.nrdl_for_disease(names_zh_disease)
    seen, uniq = set(), []
    for r in rows:
        key = (r.get("list"), r.get("no"), r.get("name_zh"))
        if key not in seen:
            seen.add(key)
            uniq.append(r)
    p = nrdl.get("provenance") or {}
    out.sources.extend(ca.sources_for(nrdl, [r.get("source_id") for r in uniq], "国家医保药品目录 (NHSA)"))
    paid = [ca.compact_nrdl(r, nrdl) for r in uniq if r.get("list") == "NRDL"]
    cidl = [ca.compact_nrdl(r, nrdl) for r in uniq if r.get("list") == "CIDL"]
    if paid:
        note = "listed in the national reimbursement list; 医保 pays within the restriction shown"
    elif names_zh_drugs or names_zh_disease:
        note = ("no NRDL row found for "
                + (f"the Chinese name(s) {', '.join(names_zh_drugs[:4])}" if names_zh_drugs else "")
                + (" or " if names_zh_drugs and names_zh_disease else "")
                + (f"a payment restriction naming {', '.join(names_zh_disease[:3])}" if names_zh_disease else "")
                + " — a drug listed under another Chinese name, or listed without a disease restriction, is "
                  "not found this way; this is not proof that 医保 does not cover it")
    else:
        note = ("not checked: no Chinese drug or disease name is known for this query (the 医保 list is in Chinese); "
                "give the Chinese generic name, e.g. 诺西那生钠")
    return {"status": "listed" if paid else ("not_found" if (names_zh_drugs or names_zh_disease) else "not_checked"),
            "edition": f"{p.get('edition')} ({p.get('document_no')}, in force {p.get('effective')})",
            "searched": {"drug_names_zh": names_zh_drugs, "disease_names_zh": names_zh_disease},
            "rows": paid[:MAX_DRUGS], "commercial_insurance_list": cidl[:5], "note": note,
            "coverage": _short(ca.coverage(nrdl), 500)}


def _nmpa_block(rows: List[Dict[str, Any]], what: str, out: Outcome) -> Dict[str, Any]:
    try:
        data = ca.load_approvals()
    except ca.NotBundled as err:
        out.warnings.append(f"China (NMPA) approvals: {err}; not checked — not evidence of anything about approval "
                            "in China")
        return {"status": "not_bundled", "note": str(err)}
    out.sources.extend(ca.sources_for(data, [e.get("source_id") for r in rows for e in (r.get("evidence") or [])]
                                      + [(r.get("approval") or {}).get("source_id") for r in rows],
                                      "China drug approvals (NMPA/CDE/MOF documents)"))
    approved = [r for r in rows if r.get("approval") is not None]
    status = ("approved_in_china" if approved else
              "named_not_approved" if rows else "not_in_bundled_list")
    return {"status": status,
            "status_meaning": {"approved_in_china": "at least one official document in the bundle says it is approved "
                                                    "or marketed in China",
                               "named_not_approved": "named in the bundled documents (e.g. 临床急需境外新药名单), but "
                                                     "none says it is approved in China — status unknown here",
                               "not_in_bundled_list": "not in the bundled documents — not evidence of anything about "
                                                      "approval"}[status],
            "rows": _china_approvals(rows),
            "note": (f"{len(rows)} row(s) for {what} in the bundled official documents" if rows else
                     f"{what}: {NOT_APPROVED_WORDING}"),
            "coverage": _short(ca.coverage(data), 600)}


def _trials_block(condition_en: Optional[str], args: argparse.Namespace, out: Outcome) -> Dict[str, Any]:
    from zebra.commands.trials import CHICTR_NOTE
    from zebra.sources import ctgov

    if args.trials <= 0:
        return {"status": "skipped", "note": "--trials 0"}
    if not condition_en:
        return {"status": "not_checked", "note": "no English disease name to search ClinicalTrials.gov with",
                "chictr": CHICTR_NOTE}
    got = attempt("ClinicalTrials.gov", lambda: ctgov.search(condition_en, country="China", status=args.status,
                                                             limit=args.trials), out.warnings)
    if got is None:
        return {"status": "unavailable", "chictr": CHICTR_NOTE}
    r = out.add(got)
    studies = [{"nct_id": s["nct_id"], "title": s["title"], "phase": "/".join(s.get("phases") or []) or s.get("study_type"),
                "status": s.get("status"), "sites_in_china": s.get("sites_in_country"),
                "sites": [f"{x.get('facility')} ({x.get('city')})" for x in (s.get("sites") or [])[:3]],
                "flags": s.get("status_flags") or None, "url": s["url"]} for s in r["studies"]]
    return {"status": "ok", "source": "ClinicalTrials.gov API v2 (sites in China)", "condition": condition_en,
            "trial_status": r["status"], "total": r.get("total"), "trials": studies,
            "set_aside_as_unrelated": [f["nct_id"] for f in r.get("filtered") or []],
            "total_any_status": r.get("total_any_status"), "chictr": CHICTR_NOTE}


def _hospital_block(province: Optional[str], out: Outcome) -> Dict[str, Any]:
    from zebra.commands import china

    got = china.hospitals(province)
    r = out.add(got)
    r = dict(r)
    if r.get("coverage"):
        r["coverage"] = _short(r["coverage"], 400)
    r["note"] = ("the network's member hospitals as NHC published them; which hospital treats which disease is not "
                 "published here — ask the hospital's rare-disease clinic (罕见病门诊)")
    return r


def _agency_rows(drugs: List[str], terms: List[str], out: Outcome) -> List[Dict[str, Any]]:
    """FDA/EMA reading per drug (zebra.sources.regulators), compact."""
    from concurrent.futures import ThreadPoolExecutor

    from zebra.sources import regulators

    drugs = list(dict.fromkeys(d for d in drugs if d))[:MAX_DRUGS]
    if not drugs:
        return []
    warn: Dict[int, List[str]] = {i: [] for i in range(len(drugs))}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futs = [pool.submit(attempt, f"agency check {d}", lambda d=d: regulators.check(d, terms), warn[i])
                for i, d in enumerate(drugs)]
        got = [f.result() for f in futs]
    rows = []
    for i, (d, g) in enumerate(zip(drugs, got)):
        out.warnings.extend(warn[i])
        if g is None:
            rows.append({"drug": d, "status": "unavailable"})
            continue
        c = out.add(g)
        # the label's own indication words are kept: they say which patients a drug is for
        # ("amenable to exon 51 skipping"), which a one-line reading cannot
        rows.append({"drug": d, **{k: {f: (c.get(k) or {}).get(f) for f in ("status", "reading", "indication", "approved_on", "date", "url")
                                       if (c.get(k) or {}).get(f) is not None} for k in ("FDA", "EMA")}})
    seen, kept = set(), []
    for src in out.sources:
        key = (src.get("db"), src.get("record"), src.get("url"))
        if key not in seen:
            seen.add(key)
            kept.append(src)
    out.sources[:] = kept
    return rows


# ------------------------------------------------------------------ flows

def _disease(info: Dict[str, Any], args: argparse.Namespace, out: Outcome) -> Dict[str, Any]:
    from zebra.commands import china
    from zebra.sources import opentargets as ot

    query = info["query"]
    cn = info.get("china")
    entries = [m for m in (cn or {}).get("matches") or [] if cn and cn["status"] in ("on_list", "qualified")]
    ot_dis = info.get("ot_disease")
    if ot_dis is not None:
        name_en = ot_dis.get("name")  # None for an id: filled from Open Targets below
    else:
        name_en = entries[0]["name_en"] if entries else (None if _CJK.search(query) else query)
    if entries and not ot_dis and name_en:
        got = attempt("Open Targets search", lambda: ot.search(re.split(r"[(/（]", name_en)[0].strip(), ("disease",), 5),
                      out.warnings)
        hits = out.add(got)["hits"] if got is not None else []
        want = re.split(r"[(/（]", name_en)[0].strip().lower()
        ot_dis = next((h for h in hits if (h["name"] or "").lower() == want), None)
        al = china.load_aliases()["_by_key"].get((entries[0]["list"], entries[0]["no"])) or {}
        if ot_dis is None:
            # the list's English name may be British or misspelt: the entry's Orphanet and corrected names
            alts = [a["text"] for a in al.get("aliases") or [] if a.get("kind") in (
                "orphanet_en", "official_en_corrected", "official_part_en") and not a.get("acronym")
                and a.get("scope") != "wider"]
            for alt in dict.fromkeys(alts):
                got = attempt("Open Targets search", lambda alt=alt: ot.search(alt, ("disease",), 5), out.warnings)
                more = out.add(got)["hits"] if got is not None else []
                ot_dis = next((h for h in more if (h["name"] or "").lower() == alt.lower()), None)
                if ot_dis is not None:
                    break
        if ot_dis is None and al.get("orpha"):
            # the entry's ORPHAcode, through Monarch's exact mapping to MONDO (Open Targets indexes MONDO):
            # 真性红细胞增多症 -> ORPHA:729 -> MONDO:0009891 'acquired polycythemia vera'
            from zebra.sources import monarch

            got = attempt("Monarch mappings", lambda: monarch.mappings(object_ids=[al["orpha"]]), out.warnings)
            mondo = [m["subject"] for m in (out.add(got) if got is not None else [])
                     if str(m.get("subject") or "").startswith("MONDO:")
                     and not str(m.get("subject_label") or "").lower().startswith("obsolete")]
            if len(set(mondo)) == 1:
                ot_dis = {"id": mondo[0].replace(":", "_"), "name": None,
                          "how": f"Monarch exactMatch of {al['orpha']}"}
        if ot_dis is None and hits:
            out.warnings.append(f"Open Targets has no disease named exactly '{want}': FDA/EMA approvals were not "
                                f"listed from it (closest: {hits[0]['name']} {hits[0]['id']}); run "
                                f"`zebra therapy <MONDO id>` to choose")
    if cn is None and name_en:
        try:
            cn = china.lookup([name_en])
        except china.ListUnavailable:
            cn = None
    res: Dict[str, Any] = {"kind": "disease", "name": name_en, "matched": info["how"],
                           "national_rare_disease_list": _china_list_block(cn)}
    # approvals: FDA/EMA for the drugs Open Targets lists at APPROVAL for this disease
    approved_drugs: List[str] = []
    terms = [t for t in [name_en] if t]
    if ot_dis:
        res["open_targets"] = {"id": ot_dis["id"], "name": ot_dis.get("name"), "how": ot_dis.get("how")}
        got = attempt(f"Open Targets disease {ot_dis['id']}", lambda: ot.disease(ot_dis["id"], drug_limit=60), out.warnings)
        d = out.add(got) if got is not None else None
        if d:
            res["open_targets"]["name"] = d["name"]
            if not name_en:
                name_en = d["name"]
                res["name"] = name_en
                if cn is None or cn.get("status") in ("not_found", "possible"):
                    try:
                        cn = china.lookup([name_en])
                        res["national_rare_disease_list"] = _china_list_block(cn)
                        entries = [m for m in cn["matches"] if cn["status"] in ("on_list", "qualified")]
                    except china.ListUnavailable:
                        pass
            terms = [d["name"]] + [s for s in d.get("synonyms") or [] if len(s) >= 5 and not s.isupper()]
            approved_drugs = [r["drug"] for r in d["drugs"]["rows"] if r.get("stage") == "APPROVAL"]
    names_zh = _zh_names_for(cn, query)
    agency = _agency_rows(approved_drugs, terms, out)
    try:
        list_keys = [(m["list"], m["no"]) for m in entries]
        china_rows = ca.approvals_for_disease(list_keys, names_zh, terms)
        unnamed = [r.get("drug_zh") for r in china_rows if not r.get("inn") and r.get("drug_zh")]
        for a in agency:  # does the bundled China table hold the drugs FDA/EMA approved for this disease?
            a["china_bundled"] = _china_reading(ca.approvals_for_drug([a["drug"]]))
            if a["china_bundled"].startswith("not in the bundled list") and unnamed:
                a["china_bundled"] = ("no bundled row carries this English name; "
                                      f"{len(unnamed)} bundled row(s) for this disease have only a Chinese name "
                                      f"({', '.join(unnamed[:3])}) and may be the same drug — see china_nmpa")
    except ca.NotBundled:
        china_rows = []
    res["approvals"] = {
        "fda_ema": {"source": "FDA labels / Drugs@FDA (openFDA) and the EMA medicines register, for the drugs Open "
                              "Targets lists at APPROVAL for this disease",
                    "drugs": agency,
                    "note": None if agency else "no drug at APPROVAL stage in Open Targets for this disease, or Open "
                                                "Targets did not answer (see warnings)"},
        "china_nmpa": _nmpa_block(china_rows, f"'{name_en or query}'", out),
    }
    zh_drugs = [n for r in china_rows for n in ([r.get("drug_zh")] + list(r.get("names_zh") or [])) if n]
    res["nrdl_医保"] = _nrdl_block(zh_drugs, names_zh, out)
    res["trials_china"] = _trials_block(name_en, args, out)
    res["hospitals_协作网"] = _hospital_block(args.province, out)
    return res


def _drug(info: Dict[str, Any], args: argparse.Namespace, out: Outcome) -> Dict[str, Any]:
    from zebra.sources import regulators

    query = info["query"]
    rows = list(info.get("drug_rows") or [])
    try:
        if not rows:
            rows = ca.approvals_for_drug([query])
    except ca.NotBundled:
        rows = []
    if not _CJK.search(query):
        inn = query  # the name the user gave; a bundled English field can carry a brand (Galafold（Migalastat…）)
    else:
        inn = next((ca.clean_inn(r.get("inn")) for r in rows if ca.clean_inn(r.get("inn"))), None)
    zh = [query] if _CJK.search(query) else []
    zh += [n for r in rows for n in ([r.get("drug_zh")] + list(r.get("names_zh") or [])) if n]
    res: Dict[str, Any] = {"kind": "drug", "name": query, "inn": inn, "matched": info["how"]}
    fda: Dict[str, Any] = {"status": "not_checked", "note": "no English (INN) name known for this drug"}
    ema: Dict[str, Any] = {"status": "not_checked", "note": "no English (INN) name known for this drug"}
    if inn:
        got = attempt("openFDA Drugs@FDA", lambda: regulators.fda_by_ingredient(inn), out.warnings)
        if got is not None:
            apps = out.add(got)["applications"]
            fda = {"source": "Drugs@FDA (openFDA), applications whose products contain this ingredient",
                   "status": "approved" if any(a["approved_on"] for a in apps) else "no_record",
                   "applications": apps[:5],
                   "note": None if apps else f"no FDA application found for ingredient '{regulators.moiety(inn)}' — "
                                             "not proof it is unapproved (biologics licensed by CBER may be missing)"}
        got = attempt("EMA medicines", lambda: regulators.ema(inn), out.warnings)
        if got is not None:
            recs = out.add(got)["records"]
            ema = {"source": "EMA medicines register (centrally authorised, refused, withdrawn)",
                   "status": (recs[0]["status"] if recs else "no_record"),
                   "records": [{k: r.get(k) for k in ("medicine", "status", "authorised_on", "refused_on",
                                                       "withdrawn_on", "indication", "url")} for r in recs[:4]],
                   "note": None if recs else "no centrally authorised, refused or withdrawn EU medicine with this "
                                             "INN (a national EU authorisation is not in this register)"}
    res["approvals"] = {"FDA": fda, "EMA": ema, "china_nmpa": _nmpa_block(rows, f"'{query}'", out)}
    res["nrdl_医保"] = _nrdl_block([n for n in dict.fromkeys(zh) if n], [], out)
    res["trials_china"] = {"status": "not_checked",
                           "note": "trials are searched by disease: run `zebra access <disease>` or "
                                   "`zebra trials <disease> --term <drug> --country China`"}
    return res


def _gene(info: Dict[str, Any], args: argparse.Namespace, out: Outcome) -> Dict[str, Any]:
    from zebra.sources import opentargets as ot

    t = info["target"]
    res: Dict[str, Any] = {"kind": "gene", "name": t["name"], "ensembl": t["id"], "matched": info["how"]}
    got = attempt(f"Open Targets target {t['id']}", lambda: ot.target(t["id"]), out.warnings)
    tgt = out.add(got) if got is not None else None
    if not tgt:
        return res
    approved = [r for r in tgt["drugs"]["rows"] if r.get("stage") == "APPROVAL"][:MAX_DRUGS]
    res["approved_drugs_on_target"] = []
    for r in approved:
        try:
            cnr = ca.approvals_for_drug([r["drug"]])
        except ca.NotBundled:
            cnr = []
        res["approved_drugs_on_target"].append({
            "drug": r["drug"], "indications": r.get("indications"), "open_targets_regulatory": r["regulatory"]["headline"],
            "china_bundled": _china_reading(cnr)})
    res["top_diseases"] = tgt.get("top_diseases")
    res["next"] = ("access is answered per disease or per drug: run `zebra access <disease>` for 医保, China trials "
                   "and hospitals, or `zebra access <drug>`")
    return res


# ------------------------------------------------------------------ render

def _render(res: Dict[str, Any]) -> str:
    L: List[str] = []
    k = res.get("kind")
    L.append(f"{res.get('query')}: {k}" + (f" — {res.get('name')}" if res.get("name") and res.get("name") != res.get("query") else "")
             + f" (matched: {res.get('matched')})")
    if k == "unresolved":
        for c in res.get("candidates") or []:
            L.append(f"  candidate: {c['list']}#{c['no']} {c['name_zh']} [{c['match']}]")
        for c in res.get("search_hits") or []:
            L.append(f"  Open Targets: {c['name']} ({c['entity']} {c['id']})")
        return "\n".join(L)
    nl = res.get("national_rare_disease_list")
    if nl:
        ents = "; ".join(f"{e['list']}#{e['no']} {e['name_zh']} [{e['match']}]" for e in nl.get("entries") or [])
        L.append(f"national rare-disease list: {nl['status']}" + (f" — {ents}" if ents else ""))
    ap = res.get("approvals") or {}
    if "fda_ema" in ap:
        L.append("approvals (FDA / EMA, for drugs Open Targets lists at APPROVAL):")
        for d in ap["fda_ema"]["drugs"] or []:
            parts = []
            for a in ("FDA", "EMA"):
                b = d.get(a) or {}
                if b:
                    parts.append(f"{a} {b.get('status')}" + (f" {b.get('approved_on') or b.get('date')}"
                                                             if (b.get('approved_on') or b.get('date')) else ""))
            L.append(f"  {d['drug']}: " + "; ".join(parts) + (f" · China: {d['china_bundled']}" if d.get("china_bundled") else ""))
        if ap["fda_ema"].get("note"):
            L.append(f"  {ap['fda_ema']['note']}")
    for a in ("FDA", "EMA"):
        b = ap.get(a)
        if b:
            if b.get("applications"):
                L.append(f"{a}: " + "; ".join(f"{x['application']} {', '.join(x['brands'][:2])} approved {x['approved_on']}"
                                              for x in b["applications"][:3]))
            elif b.get("records"):
                L.append(f"{a}: " + "; ".join(f"{x['medicine']} {x['status']} {x.get('authorised_on') or x.get('refused_on') or ''}"
                                              f" {x['url']}" for x in b["records"][:3]))
            else:
                L.append(f"{a}: {b.get('status')} — {b.get('note')}")
    cn = ap.get("china_nmpa")
    if cn:
        L.append(f"China (NMPA, bundled official documents): {cn['status']}"
                 + (f" — {cn['status_meaning']}" if cn.get("status_meaning") else ""))
        for r in cn.get("rows") or []:
            apv = r.get("approval") or {}
            when = apv.get("date") or apv.get("year") or (f"announced {apv['announced']}" if apv.get("announced") else None)
            L.append(f"  {r.get('drug_zh') or '(no Chinese name)'}" + (f" ({r['inn']})" if r.get("inn") else "")
                     + ((" — approved in China" + (f" {when}" if when else " (date not stated)")
                         + (f", {apv['type']}" if apv.get("type") else "") + f" [{apv.get('source_id')}]")
                        if r.get("approved_in_china") else
                        " — named in these documents, but none says it is approved in China")
                     + (f" · {r['indication_zh']}" if r.get("indication_zh") else "")
                     + (f" · {'; '.join(r['evidence'])}" if r.get("evidence") else ""))
        if cn.get("note"):
            L.append(f"  {cn['note']}")
    nr = res.get("nrdl_医保")
    if nr:
        L.append(f"医保 ({nr.get('edition', 'NRDL')}): {nr['status']}")
        for r in nr.get("rows") or []:
            L.append(f"  {r['name_zh']} [{r.get('section')} #{r.get('no')}, {r.get('class') or '-'}类] {r['restriction']}"
                     + (f" (协议期 {r['agreement_period']})" if r.get("agreement_period") else "")
                     + (f" — found because {r['matched_by']}" if r.get("matched_by") else ""))
        for r in nr.get("commercial_insurance_list") or []:
            L.append(f"  {r['name_zh']} — {r['list']}: {r.get('indication')}")
        if nr.get("note") and nr["status"] != "listed":
            L.append(f"  {nr['note']}")
    tr = res.get("trials_china")
    if tr:
        if tr.get("status") == "ok":
            L.append(f"trials with a site in China (ClinicalTrials.gov, {tr['trial_status']}): {tr['total']} — showing "
                     f"{len(tr['trials'])}" + (f"; {len(tr['set_aside_as_unrelated'])} set aside as unrelated"
                                               if tr.get("set_aside_as_unrelated") else ""))
            for s in tr["trials"]:
                L.append(f"  {s['nct_id']} · {s.get('phase')} · {s.get('status')} · {s['title']} · "
                         f"{s.get('sites_in_china')} site(s) in China: {'; '.join(s['sites'])} {s['url']}")
            if not tr["trials"] and tr.get("total_any_status"):
                L.append(f"  none {tr['trial_status']}; {tr['total_any_status']} in any status")
        else:
            L.append(f"trials in China: {tr.get('status')} — {tr.get('note', '')}")
        if tr.get("chictr"):
            L.append(f"  note: {tr['chictr']}")
    ho = res.get("hospitals_协作网")
    if ho:
        if ho.get("status") == "ok":
            head = (f"全国罕见病诊疗协作网 in {ho['province']}: {ho['count']}" if ho.get("province") else
                    f"全国罕见病诊疗协作网 lead hospitals: {ho['count']} (--province <省> for every member)")
            L.append(head)
            for h in ho["hospitals"][:40]:
                L.append(f"  {h.get('province')} {h['name']} [{h.get('role_zh')}]")
        else:
            L.append(f"协作网 hospitals: {ho.get('status')}")
    for d in res.get("approved_drugs_on_target") or []:
        L.append(f"  {d['drug']}: {d['open_targets_regulatory']} · China: {d['china_bundled']}")
    if res.get("next"):
        L.append(f"note: {res['next']}")
    return "\n".join(L)


def _run(args: argparse.Namespace) -> Outcome:
    query = " ".join(args.query).strip()
    if not query:
        raise UsageError("give a drug (诺西那生钠 / nusinersen), a disease (脊髓性肌萎缩症 / Dravet syndrome) or a gene symbol")
    status = (args.status or "RECRUITING").upper()
    from zebra.sources import ctgov

    if status not in ctgov.STATUSES:
        raise UsageError(f"--status must be one of {', '.join(ctgov.STATUSES)}")
    args.status = status
    if args.trials < 0:
        raise UsageError("--trials must be 0 or more")
    out = Outcome(None, query={"query": query, "as": args.kind, "province": args.province, "trials": args.trials,
                               "status": status})
    info = _classify(query, args.kind, out)
    kind = info["kind"]
    if kind == "disease":
        res = _disease(info, args, out)
    elif kind == "drug":
        res = _drug(info, args, out)
    elif kind == "gene":
        if not info.get("target"):
            raise UsageError("--as gene needs a gene symbol Open Targets knows")
        res = _gene(info, args, out)
    else:
        res = {"kind": "unresolved", "matched": info.get("how"), "candidates": info.get("candidates"),
               "search_hits": info.get("search_hits"),
               "note": "not answered: no exact disease, drug or gene name matched — pick one of the candidates and ask "
                       "again with its exact name or id (MONDO_…), or give the official Chinese or English name"}
        if args.province:
            res["hospitals_协作网"] = _hospital_block(args.province, out)
    res["query"] = query
    out.result = res
    out.text = _render(res)
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("access", help="China and international access in one answer: approvals (FDA/EMA/NMPA), "
                                      "医保 (NRDL), China trials, 协作网 hospitals — for a drug, disease or gene")
    p.add_argument("query", nargs="+", help="drug (诺西那生钠 / nusinersen), disease (脊髓性肌萎缩症 / Dravet syndrome) "
                                           "or gene symbol")
    p.add_argument("--as", dest="kind", choices=("auto", "drug", "disease", "gene"), default="auto",
                   help="what the query is (default: worked out from the name)")
    p.add_argument("--province", help="list every 全国罕见病诊疗协作网 hospital in this province (浙江 / 浙江省)")
    p.add_argument("--trials", type=int, default=8, help="trials with a site in China to list (0 skips the search)")
    p.add_argument("--status", default="RECRUITING", type=str.upper,
                   help="trial status for the China search (RECRUITING default; ANY for all)")
    p.set_defaults(func=_run)
