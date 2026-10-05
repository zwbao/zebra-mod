"""zebra therapy: drugs and clinical candidates for a disease or a gene (Open Targets), tractability, EU orphan designations.

The query is tried as both a gene symbol and a disease against Open Targets
search; the result says which matched and how. This is the data with its
sources — mechanism fit and tiering are the reader's job.
"""

from __future__ import annotations

import argparse
import re
from typing import Any, Dict, List, Optional

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
            matched.append({"as": "disease", "id": first["id"], "name": first["name"],
                            "how": "first disease hit of Open Targets search (not an exact name match: check it is the disease meant)"})
    return matched


def _drug_line(d: Dict[str, Any], with_indications: bool) -> str:
    moa = "; ".join(f"{m['mechanism']} [{', '.join(m['targets'])}]" for m in d["mechanisms"]) or "mechanism not given"
    stage = d.get("stage") or "?"
    if d.get("drug_max_stage") and d["drug_max_stage"] != stage:
        stage += f" (drug max {d['drug_max_stage']})"
    ev = "; ".join(f"{e['source']} {e['stage']}" + (f" {e['status']}" if e.get("status") else "") + (f" {e['url']}" if e.get("url") else "")
                   for e in d["evidence"][:3])
    line = f"  {stage:<12} {d['drug']} ({d['chembl_id']}, {d.get('type')}) — {moa}"
    if with_indications and d.get("indications"):
        line += f"\n      indications: {', '.join(d['indications'])}"
    line += f"\n      {d['reports']} reports ({d['trials']} ClinicalTrials.gov) · {ev}"
    return line


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
    lines.append(f"drugs and clinical candidates for this disease: {dr['count']} (showing {dr['shown']}, {dr['order']})")
    lines += [_drug_line(x, False) for x in dr["rows"]] or ["  none in Open Targets"]
    if d["top_targets"]:
        lines.append("top associated targets (Open Targets association score): " +
                     ", ".join(f"{t['symbol']} {t['score']}" for t in d["top_targets"]))
    terms = orphan.match_terms(d["name"], d.get("synonyms") or [])
    od = attempt("EMA orphan designations", lambda: orphan.ema_designations(terms), out.warnings)
    if od is not None:
        o = out.add(od)
        out.result["orphan_designations"] = o
        lines.append(f"EU orphan designations (EMA register, dataset {o['dataset_timestamp']}) naming "
                     f"{' / '.join(o['terms'])}: {o['count']}" + (f" (showing {len(o['designations'])})" if o["count"] > len(o["designations"]) else ""))
        for x in o["designations"]:
            lines.append(f"  {x['status']:<9} {x['substance']}" + (f" ({x['medicine']})" if x.get("medicine") else "")
                         + f" — {x['intended_use']} — {x['eu_number']} {x['date']} {x['url']}")


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
