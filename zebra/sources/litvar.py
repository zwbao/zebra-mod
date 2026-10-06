"""LitVar2 (NCBI): which papers mention a variant.

API: https://www.ncbi.nlm.nih.gov/research/litvar2-api (no key; NCBI pacing).
LitVar keeps separate records for spellings it has not linked: SCN1A
p.Arg712* is both `litvar@#6323#p.R712*` (1 paper) and `litvar@rs794726730##`
(found by "SCN1A p.Arg712Ter" / "SCN1A R712X", 21 papers). So a stop-gain is
looked up under its HGVS-equivalent spellings (*, Ter, X) and the records are
merged; each match says which spellings found it and whether it was LitVar's
first suggestion for one of them (`top`). Autocomplete also suggests other
variants of the same gene after the best hit: only `top` matches are the
variant asked about, the rest are LitVar's suggestions.
"""

from __future__ import annotations

import re
import urllib.parse
from typing import Any, Dict, List, Optional

from zebra.core import Outcome
from zebra.http import SourceError, get_json, source_record

BASE = "https://www.ncbi.nlm.nih.gov/research/litvar2-api"
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H",
       "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V", "Sec": "U", "Pyl": "O"}
AA1 = {v: k for k, v in AA3.items()}
STOP_RE = re.compile(r"^(?P<prefix>p\.)?\(?(?P<ref>[A-Z][a-z]{2}|[A-Z])(?P<pos>\d+)(?P<alt>Ter|\*|X)\)?$")
RSID_RE = re.compile(r"^rs\d+$", re.I)


def spellings(variant_text: str) -> List[str]:
    """The text as given, plus HGVS-equivalent stop-codon spellings LitVar indexes separately."""
    text = variant_text.strip()
    out = [text]
    m = STOP_RE.match(text)
    if m:
        ref, pos, given = m.group("ref"), m.group("pos"), m.group("alt")
        ref3 = ref if len(ref) == 3 else AA1.get(ref)
        ref1 = AA3.get(ref) if len(ref) == 3 else ref
        candidates = {"Ter": f"p.{ref3}{pos}Ter" if ref3 else None,
                      "X": f"p.{ref1}{pos}X" if ref1 else None,
                      "*": f"p.{ref1}{pos}*" if ref1 else None}
        for stop, alt in candidates.items():
            if stop != given and alt and alt not in out:
                out.append(alt)
    return out


def _with_gene(text: str, gene: Optional[str]) -> str:
    if not gene or RSID_RE.match(text) or ":" in text or gene.upper() in text.upper().split():
        return text
    return f"{gene} {text}"


def lookup(variant_text: str, gene: Optional[str] = None) -> Outcome:
    """LitVar records for a variant (rsID, HGVS, protein change; `gene` narrows a bare change).

    Returns {"query", "gene", "queries", "matches": [{"litvar_id", "rsid", "gene", "genes", "hgvs", "name",
    "pmid_count", "clinical_significance", "spellings", "top"}]} in LitVar's order.
    """
    text = (variant_text or "").strip()
    if not text:
        raise ValueError("empty variant")
    queries = [_with_gene(s, gene) for s in spellings(text)]
    merged: Dict[str, Dict[str, Any]] = {}
    sources = []
    warnings: List[str] = []
    failures: List[SourceError] = []
    for q in queries:
        try:
            resp = get_json(f"{BASE}/variant/autocomplete/", source="LitVar2", params={"query": q}, cache_ttl=14 * 86400)
        except SourceError as err:
            failures.append(err)
            warnings.append(f"LitVar2 '{q}' unavailable: {err.message} (HTTP {err.status or '-'})")
            continue
        sources.append(source_record("LitVar2 autocomplete", q, resp))
        rows = resp.json()
        for rank, row in enumerate(rows if isinstance(rows, list) else []):
            lid = row.get("_id")
            if not lid:
                continue
            m = merged.get(lid)
            if m is None:
                genes = row.get("gene") or []
                m = merged[lid] = {
                    "litvar_id": lid,
                    "rsid": row.get("rsid"),
                    "gene": genes[0] if genes else None,
                    "genes": genes,
                    "hgvs": row.get("hgvs"),
                    "name": row.get("name"),
                    "pmid_count": row.get("pmids_count"),
                    "clinical_significance": row.get("data_clinical_significance") or None,
                    "spellings": [],
                    "top": False,
                }
            m["spellings"].append(q)
            if rank == 0:
                m["top"] = True
    if failures and len(failures) == len(queries):
        raise failures[-1]
    matches = list(merged.values())
    # An rsID names one locus already, and LitVar's gene label for it can differ from the gene the
    # card annotated (m.3243A>G, rs199474657, is labelled MT-ND1/MT-ND2): never gene-filter it.
    if gene and not RSID_RE.match(text):
        keep = [m for m in matches if gene.upper() in [g.upper() for g in m["genes"]]]
        dropped = len(matches) - len(keep)
        if dropped:
            warnings.append(f"LitVar2: {dropped} suggestion(s) in other genes than {gene} left out")
        matches = keep
    if not matches:
        warnings.append(f"LitVar2 has no record for {text}" + (f" in {gene}" if gene else ""))
    return Outcome({"query": text, "gene": gene, "queries": queries, "matches": matches},
                   sources=sources, warnings=warnings)


def _quoted(litvar_id: str) -> str:
    return urllib.parse.quote(litvar_id.strip(), safe="")


def pmids(litvar_id: str, limit: Optional[int] = None) -> Outcome:
    """PMIDs LitVar links to one record (order as LitVar returns them)."""
    resp = get_json(f"{BASE}/variant/get/{_quoted(litvar_id)}/publications", source="LitVar2", cache_ttl=7 * 86400)
    data = resp.json()
    ids = [str(p) for p in data.get("pmids") or []]
    result = {"litvar_id": litvar_id, "pmid_count": data.get("pmids_count", len(ids)),
              "pmids": ids[:limit] if limit else ids}
    return Outcome(result, sources=[source_record("LitVar2 publications", litvar_id, resp)])


def get(litvar_id: str) -> Outcome:
    """One LitVar record: rsID, ClinGen allele id, genomic position, dbSNP class."""
    resp = get_json(f"{BASE}/variant/get/{_quoted(litvar_id)}", source="LitVar2", cache_ttl=14 * 86400)
    d = resp.json()
    result = {"litvar_id": d.get("_id", litvar_id), "rsid": d.get("rsid"), "gene": d.get("gene"), "name": d.get("name"),
              "hgvs": d.get("hgvs"), "clingen_ids": d.get("clingen_ids"), "position": d.get("data_chromosome_base_position"),
              "snp_class": d.get("data_snp_class"), "clinical_significance": d.get("data_clinical_significance")}
    return Outcome(result, sources=[source_record("LitVar2", litvar_id, resp)])


# ------------------------------------------------------------ allele matching (CP1-9)
# LitVar's rsID-level record merges every allele dbSNP files under the rsID: CFTR p.Gly542Arg
# (a VUS) and p.Gly542Ter (a common CF allele) are both rs113993959, and the record, named
# "p.G542X", carries 842 PMIDs. Each record's own spelling (`hgvs`/`name`) says which allele it
# is about; a record for another allele at the same residue is excluded and counted as excluded.

_P_RE = re.compile(r"^p\.\(?(?P<ref>[A-Z][a-z]{2}|[A-Z])(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2}|[A-Z*]|Ter|del|dup|fs.*|=)\)?$")
_C_RE = re.compile(r"^(?:[A-Za-z0-9_.]+:)?c\.(?P<pos>[-*]?\d+(?:[+-]\d+)?(?:_[-*]?\d+(?:[+-]\d+)?)?)(?P<change>.+)$")
_C_SUB_RE = re.compile(r"^([ACGT])>([ACGT])$", re.I)
_SUBSTITUTION_ALTS = set(AA3.values()) | {"*"}


def _p_key(text: Optional[str]) -> Optional[tuple]:
    """(ref, pos, alt) in one-letter form, stops as '*', '=' as the reference residue, frameshifts as 'fs'."""
    if not text:
        return None
    t = text.split(":", 1)[-1].strip()
    m = _P_RE.match(t)
    if not m:
        return None
    ref, alt = m.group("ref"), m.group("alt")
    ref1 = AA3.get(ref, ref) if len(ref) == 3 else ref
    if alt in ("Ter", "*", "X"):
        alt1 = "*"
    elif alt == "=":
        alt1 = ref1  # p.Thr854= is p.T854T
    elif alt.startswith("fs"):
        alt1 = "fs"
    elif len(alt) == 3 and alt in AA3:
        alt1 = AA3[alt]
    else:
        alt1 = alt
    return ref1, int(m.group("pos")), alt1


def _c_key(text: Optional[str]) -> Optional[tuple]:
    """(position, change) with the bases after del/dup dropped: c.5266dupC == c.5266dup."""
    if not text or "c." not in text:
        return None
    m = _C_RE.match(text.split(":", 1)[-1].strip() if ":" in text and "c." in text.split(":", 1)[-1] else text.strip())
    if not m:
        return None
    change = re.sub(r"^(del|dup)[ACGTN]+$", r"\1", m.group("change").strip(), flags=re.I)
    return m.group("pos"), change.upper()


def allele_match(record: Dict[str, Any], hgvs_p: Optional[str], hgvs_c: Optional[str]) -> str:
    """'same', 'different' or 'unknown' for one LitVar record against this variant's p. and c. HGVS.

    'different' only when both sides are SUBSTITUTIONS at the same position with another alternate
    (c.1624G>A vs c.1624G>T; p.Gly542Arg vs p.G542X): the case of LitVar's rsID record merging
    another allele. Spelling variants of the same change (c.5266dupC / c.5266dup, p.Thr854= / p.T854T)
    are 'same'; anything else -- another position (possibly another transcript's numbering), an
    indel against a substitution -- is 'unknown' and never excluded.
    """
    spelled = [s for s in (record.get("hgvs"), record.get("name")) if s]
    ours_p = _p_key(hgvs_p)
    ours_c = _c_key(hgvs_c)
    verdict = "unknown"
    for s in spelled:
        k = _p_key(s)
        if k and ours_p:
            if k == ours_p:
                return "same"
            if k[:2] == ours_p[:2] and k[2] in _SUBSTITUTION_ALTS and ours_p[2] in _SUBSTITUTION_ALTS:
                verdict = "different"
            continue
        ck = _c_key(s)
        if ck and ours_c:
            if ck == ours_c:
                return "same"
            if ck[0] == ours_c[0] and _C_SUB_RE.match(ck[1]) and _C_SUB_RE.match(ours_c[1]):
                verdict = "different"
    return verdict
