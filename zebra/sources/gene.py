"""Gene card: identity, diseases, validity, dosage, constraint, panels and protein for one gene.

HGNC REST (approved symbol; previous symbols and aliases resolved, never guessed)
-> in parallel: Ensembl lookup, Monarch causal gene->disease associations (with
the disease's inheritance from Monarch), ClinGen validity + dosage, gnomAD
constraint, PanelApp England + Australia, UniProt + AlphaFold DB.
Each part runs through zebra.core.attempt: a failing source becomes a warning.
"""

from __future__ import annotations

import urllib.parse
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Dict, List, Optional, Tuple

from zebra.core import Outcome, UsageError, attempt
from zebra.http import get_json, source_record
from zebra.sources import clingen, ensembl, gnomad, panelapp, uniprot

HGNC = "https://rest.genenames.org"
MONARCH = "https://api-v3.monarchinitiative.org/v3/api"
MONARCH_UI = "https://monarchinitiative.org"
MAX_INHERITANCE_LOOKUPS = 10


# ---------------------------------------------------------------- HGNC

def parse_hgnc(doc: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "symbol": doc.get("symbol"), "hgnc_id": doc.get("hgnc_id"), "name": doc.get("name"),
        "status": doc.get("status"), "locus_type": doc.get("locus_type"), "location": doc.get("location"),
        "aliases": (doc.get("alias_symbol") or [])[:10], "previous_symbols": (doc.get("prev_symbol") or [])[:10],
        "entrez_id": doc.get("entrez_id"), "ensembl_gene_id": doc.get("ensembl_gene_id"),
        "omim": [f"OMIM:{m}" for m in doc.get("omim_id") or []], "mane_select": doc.get("mane_select") or [],
        "uniprot_ids": doc.get("uniprot_ids") or [], "orphanet": doc.get("orphanet"),
        "url": f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{doc.get('hgnc_id')}",
    }


def _hgnc_get(path: str, label: str):
    resp = get_json(f"{HGNC}/{path}", source="HGNC", cache_ttl=30 * 86400)
    return resp, (resp.json().get("response") or {})


def hgnc(symbol: str) -> Outcome:
    """Approved symbol record; a previous symbol or a unique alias is followed with a warning."""
    sym = symbol.strip()
    q = urllib.parse.quote(sym)
    resp, body = _hgnc_get(f"fetch/symbol/{q}", sym)
    sources = [source_record("HGNC", sym, resp)]
    if body.get("numFound"):
        rec = parse_hgnc(body["docs"][0])
        return Outcome(rec, sources=sources)
    for field, what in (("prev_symbol", "a previous symbol"), ("alias_symbol", "an alias")):
        resp2, body2 = _hgnc_get(f"search/{field}/{q}", sym)
        sources.append(source_record("HGNC search", f"{field}:{sym}", resp2))
        docs = body2.get("docs") or []
        if len(docs) == 1:
            new = docs[0]["symbol"]
            resp3, body3 = _hgnc_get(f"fetch/symbol/{urllib.parse.quote(new)}", new)
            sources.append(source_record("HGNC", new, resp3))
            if body3.get("numFound"):
                return Outcome(parse_hgnc(body3["docs"][0]), sources=sources,
                               warnings=[f"{sym} is {what} of {new} (HGNC); showing {new}"])
        elif len(docs) > 1:
            cands = ", ".join(f"{d['symbol']} ({d['hgnc_id']})" for d in docs[:10])
            raise UsageError(f"{sym} is {what} of several genes: {cands} — give the approved symbol")
    raise UsageError(f"{sym} is not an HGNC-approved, previous or alias symbol (rest.genenames.org)")


# ---------------------------------------------------------------- Monarch

def parse_monarch_associations(items: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    by_id: Dict[str, Dict[str, Any]] = {}
    for it in items:
        oid = it.get("object")
        if not oid:
            continue
        row = by_id.setdefault(oid, {"id": oid, "name": it.get("object_label"), "sources": [], "xrefs": [],
                                     "inheritance": None, "url": f"{MONARCH_UI}/{oid}"})
        src = str(it.get("primary_knowledge_source") or "").replace("infores:", "")
        if src and src not in row["sources"]:
            row["sources"].append(src)
        orig = it.get("original_object")
        if orig and orig != oid and orig not in row["xrefs"]:
            row["xrefs"].append(orig)
    return list(by_id.values())


def monarch_diseases(hgnc_id: str) -> Outcome:
    """Causal gene->disease associations (OMIM, ClinGen, Orphanet via Monarch), deduplicated on disease."""
    resp = get_json(f"{MONARCH}/association", source="Monarch",
                    params={"subject": hgnc_id, "category": "biolink:CausalGeneToDiseaseAssociation", "limit": 100},
                    cache_ttl=14 * 86400, timeout=60)
    data = resp.json()
    rows = parse_monarch_associations(data.get("items") or [])
    warnings: List[str] = []
    sources = [source_record("Monarch associations", hgnc_id, resp, url=f"{MONARCH_UI}/{hgnc_id}")]
    todo = rows[:MAX_INHERITANCE_LOOKUPS]
    if len(rows) > MAX_INHERITANCE_LOOKUPS:
        warnings.append(f"Monarch: inheritance looked up for the first {MAX_INHERITANCE_LOOKUPS} of {len(rows)} diseases")

    def one(row: Dict[str, Any]) -> Tuple[Dict[str, Any], Optional[Any], Optional[str]]:
        try:
            r = get_json(f"{MONARCH}/entity/{urllib.parse.quote(row['id'])}", source="Monarch", cache_ttl=30 * 86400, timeout=60)
            return row, r, None
        except Exception as err:  # noqa: BLE001 - reported as a warning below
            return row, None, str(err)

    with ThreadPoolExecutor(max_workers=4) as pool:
        for row, r, err in pool.map(one, todo):
            if r is None:
                warnings.append(f"Monarch entity {row['id']} unavailable: {err}")
                continue
            ent = r.json()
            inh = ent.get("inheritance")
            if isinstance(inh, dict) and inh.get("name"):
                row["inheritance"] = inh.get("name")
                row["inheritance_id"] = inh.get("id")
            sources.append(source_record("Monarch entity", row["id"], r, url=row["url"]))
    rows.sort(key=lambda r: (0 if "clingen" in r["sources"] else 1, r["name"] or ""))
    return Outcome({"total": data.get("total"), "diseases": rows}, sources=sources, warnings=warnings)


# ---------------------------------------------------------------- card

def _run(tasks: Dict[str, Tuple[str, Callable[[], Outcome]]]) -> Dict[str, Tuple[Optional[Outcome], List[str]]]:
    """Run labelled source calls in parallel, one thread per host; keep each one's warnings apart."""
    out: Dict[str, Tuple[Optional[Outcome], List[str]]] = {}

    def go(key: str):
        label, fn = tasks[key]
        w: List[str] = []
        return key, attempt(label, fn, w), w

    with ThreadPoolExecutor(max_workers=len(tasks) or 1) as pool:
        for key, res, w in pool.map(go, list(tasks)):
            out[key] = (res, w)
    return out


def _protein(acc: Optional[str], symbol: str) -> Outcome:
    prot = uniprot.entry(acc) if acc else uniprot.by_gene(symbol)
    w: List[str] = []
    af = attempt("AlphaFold DB", lambda: uniprot.alphafold(prot.result["accession"]), w)
    if af is not None:
        prot.result["alphafold"] = af.result
        prot.sources += af.sources
        prot.warnings += af.warnings
    prot.warnings += w
    return prot


def card(symbol: str) -> Outcome:
    sym = symbol.strip()
    if not sym:
        raise UsageError("give a gene symbol, e.g. SCN1A")
    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []
    ident = attempt("HGNC", lambda: hgnc(sym), warnings)
    if ident is not None:
        sources += ident.sources
        warnings += ident.warnings
        info = ident.result
    else:
        info = {"symbol": sym.upper(), "hgnc_id": None, "note": "HGNC unavailable: symbol not verified"}
    s = info["symbol"]
    hid = info.get("hgnc_id")
    acc = (info.get("uniprot_ids") or [None])[0]

    tasks: Dict[str, Tuple[str, Callable[[], Outcome]]] = {
        "ensembl": ("Ensembl lookup", lambda: ensembl.lookup_symbol(s)),
        "validity": ("ClinGen gene-disease validity", lambda: clingen.validity(symbol=s, hgnc_id=hid)),
        "dosage": ("ClinGen dosage", lambda: clingen.dosage(s)),
        "constraint": ("gnomAD constraint", lambda: gnomad.gene_constraint(s, "GRCh38", gene_id=info.get("ensembl_gene_id"))),
        "panelapp_ge": ("PanelApp (Genomics England)", lambda: panelapp.gene_panels(s, "GE", limit=10)),
        "panelapp_au": ("PanelApp Australia", lambda: panelapp.gene_panels(s, "AU", limit=10)),
        "protein": ("UniProt", lambda: _protein(acc, s)),
    }
    if hid:
        tasks["monarch"] = ("Monarch", lambda: monarch_diseases(hid))
    else:
        warnings.append("Monarch diseases skipped: no HGNC id")
    got = _run(tasks)

    def take(key: str) -> Any:
        res, w = got.get(key, (None, []))
        warnings.extend(w)
        if res is None:
            return None
        sources.extend(res.sources)
        warnings.extend(res.warnings)
        return res.result

    ens = take("ensembl")
    ensembl_part = None
    if isinstance(ens, dict):
        ensembl_part = {"id": ens.get("id"), "chrom": ens.get("seq_region_name"), "start": ens.get("start"),
                        "end": ens.get("end"), "strand": ens.get("strand"), "biotype": ens.get("biotype"),
                        "assembly": ens.get("assembly_name"), "url": f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={ens.get('id')}"}
        if info.get("ensembl_gene_id") and ens.get("id") and info["ensembl_gene_id"] != ens.get("id"):
            warnings.append(f"HGNC gives {info['ensembl_gene_id']} but Ensembl lookup of {s} gives {ens.get('id')}")
    mon = take("monarch")
    validity = take("validity")
    dosage = take("dosage")
    constraint = take("constraint")
    pa_ge = take("panelapp_ge")
    pa_au = take("panelapp_au")
    protein = take("protein")
    diseases = (mon or {}).get("diseases") if mon else None
    for d in diseases or []:
        if d.get("inheritance"):
            d["inheritance_source"] = "Monarch (HPO annotation)"
            continue
        moi = sorted({v["moi"] for v in validity or [] if v.get("mondo") == d["id"] and v.get("moi")})
        if moi:
            d["inheritance"] = "/".join(moi)
            d["inheritance_source"] = "ClinGen MOI"
    dosage_ok = got.get("dosage", (None, []))[0] is not None
    result = {
        "symbol": s, "hgnc": info, "ensembl": ensembl_part,
        "diseases": diseases,
        "clingen": {
            "validity": validity,
            "validity_note": "no ClinGen gene-disease validity curation for this gene" if validity == [] else None,
            "dosage": dosage,
            "dosage_note": "no ClinGen dosage curation for this gene" if dosage_ok and dosage is None else None,
        },
        "constraint": constraint,
        "panelapp": {"genomics_england": pa_ge, "australia": pa_au},
        "protein": protein,
    }
    return Outcome(result, sources=sources, warnings=warnings)
