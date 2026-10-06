"""zebra therapy: drugs and clinical candidates for a disease or a gene (Open Targets), tractability, EU orphan designations.

The query is tried as both a gene symbol and a disease against Open Targets
search; the result says which matched and how. This is the data with its
sources — mechanism fit and tiering are the reader's job.
"""

from __future__ import annotations

import argparse
import re
from typing import Any, Dict, List

from zebra.core import Outcome, UsageError, attempt
from zebra.sources import opentargets as ot
from zebra.sources import orphan

XREF_ONLY_RE = re.compile(r"^(OMIM|ORPHA|ORPHANET|MIM|GARD|ICD10|ICD-10|UMLS|MESH):", re.I)
FDA_NOTE = ("FDA orphan designations: not checked — the FDA OOPD database (accessdata.fda.gov) refuses automated "
            "access from this client; search it by hand at https://www.accessdata.fda.gov/scripts/opdlisting/oopd/")


def _resolve(query: str, out: Outcome) -> List[Dict[str, Any]]:
    """[{"as": "disease"|"target", "id", "name", "how"}] for the query."""
    q = query.strip()
    if ot.ENSG_RE.match(q):
        return [{"as": "target", "id": q.upper(), "name": None, "how": "Ensembl gene id"}]
    did = ot.normalise_disease_id(q)
    if did:
        return [{"as": "disease", "id": did, "name": None, "how": "disease id"}]
    if XREF_ONLY_RE.match(q):
        raise UsageError(f"Open Targets indexes diseases by name and MONDO/EFO id, not {q.split(':')[0]} numbers: give the "
                         "disease name or its MONDO id (disease_card returns both)")
    got = attempt("Open Targets search", lambda: ot.search(q, ("disease", "target"), 10), out.warnings)
    if got is None:
        return []
    res = out.add(got)
    hits = res["hits"]
    out.result["search_hits"] = [{"id": h["id"], "name": h["name"], "entity": h["entity"]} for h in hits[:6]]
    matched: List[Dict[str, Any]] = []
    tgt = next((h for h in hits if h["entity"] == "target" and (h["name"] or "").upper() == q.upper()), None)
    if tgt:
        matched.append({"as": "target", "id": tgt["id"], "name": tgt["name"], "how": "exact gene symbol"})
    dis = next((h for h in hits if h["entity"] == "disease" and (h["name"] or "").lower() == q.lower()), None)
    if dis:
        matched.append({"as": "disease", "id": dis["id"], "name": dis["name"], "how": "exact disease name"})
    elif not tgt:
        first = next((h for h in hits if h["entity"] == "disease"), None)
        if first:
            matched.append({"as": "disease", "id": first["id"], "name": first["name"], "exact": False,
                            "how": "first disease hit of Open Targets search (not an exact name match: check it is the disease meant)"})
            # F13: the caveat used to live only in `matched[0].how`, where the
            # reader of a built landscape never saw it. `zebra therapy SMA`
            # silently became "proximal spinal muscular atrophy" (MONDO_0019079),
            # while `zebra disease SMA` refuses a non-exact name and lists
            # candidates — so the two tools disagreed without saying so.
            others = [h for h in hits if h["entity"] == "disease" and h["id"] != first["id"]][:4]
            out.warnings.append(
                f"'{q}' is not an exact Open Targets disease name: the landscape below is for "
                f"{first['id']} '{first['name']}', the first search hit. Confirm it is the disease meant"
                + ("; other hits were " + ", ".join(f"{h['name']} ({h['id']})" for h in others) if others else "")
                + ". `zebra disease` resolves a name to an id without guessing")
    for m in matched:
        m.setdefault("exact", True)
    return matched


def _drug_line(d: Dict[str, Any], with_indications: bool) -> str:
    moa = "; ".join(f"{m['mechanism']} [{', '.join(m['targets'])}]" for m in d["mechanisms"]) or "mechanism not given"
    stage = d.get("stage") or "?"
    if d.get("drug_max_stage") and d["drug_max_stage"] != stage:
        stage += f" (any indication: {d['drug_max_stage']})"
    ev = "; ".join(f"{e['source']} {e['stage']}" + (f" {e['status']}" if e.get("status") else "") + (f" {e['url']}" if e.get("url") else "")
                   for e in d["evidence"][:3])
    line = f"  {stage:<12} {d['drug']} ({d['chembl_id']}, {d.get('type')}) — {moa}"
    reg = d.get("regulatory") or {}
    if reg.get("by_jurisdiction") or reg.get("headline"):
        line += f"\n      status by jurisdiction: {reg.get('headline')}"
        for j in reg.get("by_jurisdiction") or []:
            line += f"\n        {j['jurisdiction']}: {j['stage']} — {j['reading']}" + (f" {j['url']}" if j.get("url") else "")
    chk = d.get("agency_check")
    if chk:
        for k in ("FDA", "EMA"):
            b = chk.get(k) or {}
            if b.get("status") and b["status"] != "unavailable":
                line += f"\n        {k} check: {b.get('reading') or b['status'].replace('_', ' ')}" \
                        + (f" {b['url']}" if b.get("url") else "") \
                        + (f" (approved {b['approved_on']})" if b.get("approved_on") else "")
    if d.get("stage_warning"):
        line += f"\n      ! {d['stage_warning']}"
    if with_indications and d.get("indications"):
        line += f"\n      indications: {', '.join(d['indications'])}"
    line += f"\n      {d['reports']} reports ({d['trials']} ClinicalTrials.gov) · {ev}"
    return line


# E-8: words that say what a substance is made of or delivered by, not which substance
# it is. A shared "sodium", "hydrochloride" or "vector" hid real gaps (a designated
# "Newdrugamab sodium" matched VALPROATE SODIUM).
GAP_STOPWORDS = frozenset((
    "sodium", "disodium", "potassium", "calcium", "magnesium", "hydrochloride", "dihydrochloride", "hydrobromide",
    "bromide", "chloride", "acetate", "citrate", "phosphate", "sulfate", "sulphate", "mesylate", "maleate",
    "tartrate", "succinate", "fumarate", "besylate", "tosylate", "lactate", "gluconate", "monohydrate", "hydrate",
    "acid", "salt", "base", "ester", "virus", "viral", "vector", "vectors", "adeno", "associated", "serotype",
    "recombinant", "human", "humanised", "humanized", "monoclonal", "antibody", "antibodies", "fragment",
    "oligonucleotide", "antisense", "synthetic", "single", "stranded", "phosphorothioate", "morpholino",
    "cells", "cell", "autologous", "allogeneic", "derived", "expressing", "encoding", "containing", "gene",
    "genes", "protein", "fusion", "conjugated", "conjugate", "modified", "transduced", "complementary",
    "against", "targeting", "with", "from", "into", "that", "this", "type", "base", "length", "full",
    # suffixes shared by unrelated biologics: a shared "alfa" hid cipaglucosidase alfa and olipudase alfa
    "alfa", "alpha", "beta", "gamma", "delta", "pegol", "mab", "ase",
))


def _tokens(text: str) -> set:
    return {t for t in re.split(r"[^a-z0-9]+", (text or "").lower()) if len(t) > 3 and t not in GAP_STOPWORDS}


def _gap_label(m: Dict[str, Any]) -> str:
    """A short name for a designated substance: its trade name when the EMA register gives one."""
    if m.get("medicine"):
        return str(m["medicine"])
    s = str(m.get("substance") or "?")
    return s if len(s) <= 60 else s[:57].rsplit(" ", 1)[0] + "…"


def _designations_not_in_drug_list(designations: List[Dict[str, Any]], rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """EU orphan-designated substances with no row in the Open Targets drug list (P1g).

    A gap found this way is evidence, not a guess: both lists come from this
    run. Delandistrogene moxeparvovec (Elevidys) for Duchenne muscular
    dystrophy is the case the review found — it has an EU orphan designation
    and an FDA approval, and no row in Open Targets' 62 drugs for the disease.
    Matching is on a shared word of 4+ characters between the substance name
    (or the medicine's trade name) and a drug name, so a spelling difference
    does not invent a gap.
    """
    known = set()
    for r in rows or []:
        known |= _tokens(r.get("drug") or "")
    out = []
    for d in designations:
        subj = _tokens(d.get("substance") or "") | _tokens(d.get("medicine") or "")
        if not (subj & known):
            sub = str(d.get("substance") or "")
            out.append({"substance": sub if len(sub) <= 90 else sub[:87].rsplit(" ", 1)[0] + "…",
                        "medicine": d.get("medicine"), "status": d.get("status"), "eu_number": d.get("eu_number"),
                        "url": d.get("url"),
                        **({} if subj else {"why": "described only generically (no distinctive word to compare)"})})
    return out


AGENCY_CHECK_MAX = 12  # distinct APPROVAL moieties checked against FDA labels and the EMA register per disease
DESIGNATIONS_SHOWN = 15
GAP_SHOWN = 12


def _slim_rows(rows: List[Dict[str, Any]]) -> None:
    """Drop what repeats on every row (empty lists, the stock reading of a stage): the 60,000-character
    budget then holds more drug rows (Duchenne: 40 rows, of which 23 used to be cut)."""
    for r in rows:
        reg = r.get("regulatory") or {}
        for k in [k for k, v in reg.items() if v == []]:
            reg.pop(k)
        for j in reg.get("by_jurisdiction") or []:
            j.pop("reading", None)
            j.pop("source", None)
        r.pop("report_sources", None)


def _compact_agency(b: Dict[str, Any]) -> Dict[str, Any]:
    """Status, date and link per agency; the reading only when it says something the status does not."""
    keep = ("status", "approved_on", "date", "url") if b.get("status") == "approved_for_disease" else \
        ("status", "reading", "approved_on", "date", "url")
    return {f: b.get(f) for f in keep if b.get(f) is not None}


def _agency_checks(rows: List[Dict[str, Any]], disease_names: List[str], out: Outcome) -> None:
    """E-4: read the FDA label and EMA records for every APPROVAL row before anything is called (un)approved.

    Open Targets' maxClinicalStage APPROVAL with no agency report used to be
    labelled "do not read this as an approved therapy" — for cannabidiol and
    fenfluramine in Dravet syndrome, both FDA- and EMA-approved. The agency's
    own record now decides, and finding nothing is "unknown here".
    """
    from concurrent.futures import ThreadPoolExecutor

    from zebra.sources import regulators

    terms = [n for n in dict.fromkeys(disease_names) if n and len(n) >= 5 and not n.isupper()]
    # one check per active moiety: GIVINOSTAT and GIVINOSTAT HYDROCHLORIDE are the same drug
    by_moiety: Dict[str, List[Dict[str, Any]]] = {}
    for r in rows:
        if r.get("stage") == "APPROVAL" and r.get("drug"):
            by_moiety.setdefault(regulators.moiety(r["drug"]), []).append(r)
    keys = list(by_moiety)[:AGENCY_CHECK_MAX]
    if not keys or not terms:
        return
    warn: Dict[int, List[str]] = {i: [] for i in range(len(keys))}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futs = [pool.submit(attempt, f"agency check {by_moiety[k][0]['drug']}",
                            lambda k=k: regulators.check(by_moiety[k][0]["drug"], terms), warn[i])
                for i, k in enumerate(keys)]
        results = [f.result() for f in futs]
    pairs = []
    for i, (k, got) in enumerate(zip(keys, results)):
        out.warnings.extend(warn[i])
        if got is not None:
            chk = out.add(got)
            pairs += [(row, chk) for row in by_moiety[k]]
    for row, chk in pairs:
        # compact: the reading, the date and the link per agency (a full block per row overflowed the
        # 60,000-character tool budget for Duchenne's 40 rows)
        row["agency_check"] = {k: _compact_agency(chk.get(k) or {}) for k in ("FDA", "EMA")}
        if chk["approved_in"]:
            row.pop("stage_warning", None)
        elif chk["refused_or_withdrawn_in"] and not row["regulatory"]["approved_in"]:
            row["stage_warning"] = (f"agency record: {chk['headline']} — do not read this as an approved therapy "
                                    "where the record says refused or withdrawn")
    # one provenance row per distinct record: every check reads the same EMA export
    seen, kept = set(), []
    for src in out.sources:
        key = (src.get("db"), src.get("record"), src.get("url"))
        if key not in seen:
            seen.add(key)
            kept.append(src)
    out.sources[:] = kept
    skipped = [by_moiety[k][0]["drug"] for k in list(by_moiety)[AGENCY_CHECK_MAX:]]
    if skipped:
        out.warnings.append(f"{len(skipped)} APPROVAL row(s) were not checked against FDA/EMA records (cap "
                            f"{AGENCY_CHECK_MAX}): {', '.join(skipped[:6])}")


def _disease_block(m: Dict[str, Any], out: Outcome, lines: List[str]) -> None:
    got = attempt(f"Open Targets disease {m['id']}", lambda: ot.disease(m["id"]), out.warnings)
    if got is None:
        return
    d = out.add(got)
    if d is None:
        return
    out.result["disease"] = d
    m["name"] = m.get("name") or d["name"]
    lines.append(f"== Disease {d['id']} {d['name']} (Open Targets data {d['data_version']}; matched: {m['how']}) {d['url']}")
    if d["xrefs"]:
        lines.append("xrefs: " + ", ".join(d["xrefs"]))
    dr = d["drugs"]
    lines.append(f"what Open Targets holds for this disease: {dr['count']} drugs and clinical candidates "
                 f"(showing {dr['shown']}, {dr['order']}). Not a list of approved therapies — see each row's "
                 "status by jurisdiction")
    _agency_checks(dr["rows"], [d["name"]] + list(d.get("synonyms") or []), out)
    lines += [_drug_line(x, False) for x in dr["rows"]] or ["  none in Open Targets"]
    _slim_rows(dr["rows"])
    if dr.get("withdrawn_or_suspended"):
        out.warnings.append("withdrawn or suspended in at least one jurisdiction, despite showing APPROVAL as a "
                            "highest-ever stage: " + ", ".join(dr["withdrawn_or_suspended"])
                            + " — read each row's status by jurisdiction before calling anything approved")
    if d["top_targets"]:
        lines.append("top associated targets (Open Targets association score): " +
                     ", ".join(f"{t['symbol']} {t['score']}" for t in d["top_targets"]))
    terms = orphan.match_terms(d["name"], d.get("synonyms") or [])
    # every designation is compared with the drug list; only the first DESIGNATIONS_SHOWN are listed
    od = attempt("EMA orphan designations", lambda: orphan.ema_designations(terms, limit=1000), out.warnings)
    if od is not None:
        o = out.add(od)
        all_designations = list(o.get("designations") or [])
        o["designations"] = all_designations[:DESIGNATIONS_SHOWN]
        out.result["orphan_designations"] = o
        lines.append(f"EU orphan designations (EMA register, dataset {o['dataset_timestamp']}) naming "
                     f"{' / '.join(o['terms'])}: {o['count']}" + (f" (showing {len(o['designations'])})" if o["count"] > len(o["designations"]) else ""))
        for x in o["designations"]:
            lines.append(f"  {x['status']:<9} {x['substance']}" + (f" ({x['medicine']})" if x.get("medicine") else "")
                         + f" — {x['intended_use']} — {x['eu_number']} {x['date']} {x['url']}")
        missing = _designations_not_in_drug_list(all_designations, dr["rows"])
        if missing:
            out.result["coverage_gap"] = [{k: m[k] for k in ("substance", "medicine", "eu_number", "why") if m.get(k)}
                                          for m in missing[:GAP_SHOWN]]
            out.result["coverage_gap_total"] = len(missing)
            out.warnings.append(
                f"{len(missing)} EU orphan-designated substance(s) for this disease are possibly missing from the "
                "Open Targets drug list above (no shared distinctive word with any row; a substance described "
                "chemically can still be a listed INN), so an approved or late-stage therapy can be missing from it: "
                + "; ".join(_gap_label(m) for m in missing[:6])
                + (f"; +{len(missing) - 6} more (result.coverage_gap lists {min(len(missing), GAP_SHOWN)})"
                   if len(missing) > 6 else "")
                + ". Check each against the FDA (Drugs@FDA / Purple Book) and the EMA register; for China "
                  "(NMPA approvals, 医保) run `zebra access <disease>`")


def _target_block(m: Dict[str, Any], out: Outcome, lines: List[str]) -> None:
    got = attempt(f"Open Targets target {m['id']}", lambda: ot.target(m["id"]), out.warnings)
    if got is None:
        return
    t = out.add(got)
    if t is None:
        return
    out.result["target"] = t
    m["name"] = m.get("name") or t["symbol"]
    lines.append(f"== Target {t['symbol']} {t['id']} — {t['name']} (Open Targets data {t['data_version']}; matched: {m['how']}) {t['url']}")
    names = {"small_molecule": "small molecule", "antibody": "antibody", "protac": "PROTAC", "other_modalities": "other modalities"}
    for key, label in names.items():
        lines.append(f"tractability, {label}: " + (", ".join(t["tractability"].get(key) or []) or "no positive bucket"))
    dr = t["drugs"]
    lines.append(f"drugs and clinical candidates acting on {t['symbol']}: {dr['count']} (showing {dr['shown']}, {dr['order']})")
    lines += [_drug_line(x, True) for x in dr["rows"]] or ["  none in Open Targets"]
    if t["top_diseases"]:
        lines.append("top associated diseases: " + ", ".join(f"{x['name']} ({x['id']}) {x['score']}" for x in t["top_diseases"]))


def _therapy(args: argparse.Namespace) -> Outcome:
    query = " ".join(args.query).strip()
    if not query:
        raise UsageError("give a disease name, MONDO id or gene symbol")
    out = Outcome({"query": query, "matched": [], "disease": None, "target": None, "orphan_designations": None,
                   "notes": [FDA_NOTE]}, query={"query": query})
    matched = _resolve(query, out)
    out.result["matched"] = matched
    lines: List[str] = []
    if not matched:
        out.warnings.append(f"Open Targets matched no disease or target for '{query}'")
        hits = out.result.get("search_hits") or []
        if hits:
            lines.append("Open Targets search hits: " + ", ".join(f"{h['name']} ({h['entity']} {h['id']})" for h in hits))
    for m in matched:
        if m["as"] == "target":
            _target_block(m, out, lines)
        else:
            _disease_block(m, out, lines)
    if out.result["target"] and not out.result["disease"]:
        out.result["notes"].append("orphan designations are searched by disease name: run therapy with the disease name")
    lines += [f"note: {n}" for n in out.result["notes"]]
    out.text = "\n".join(lines) if lines else f"no Open Targets match for {query}"
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("therapy", help="therapy landscape: drugs and clinical candidates for a disease or gene (Open Targets), "
                                       "tractability, EU orphan designations")
    p.add_argument("query", nargs="+", help="disease name or MONDO/EFO id, or gene symbol / Ensembl id")
    p.set_defaults(func=_therapy)
