"""ClinGen: gene–disease validity and dosage sensitivity (haploinsufficiency / triplosensitivity).

ClinGen's per-gene search API ignores its query string and returns every
curation (1.9 MB / 1.5 MB), so the public bulk downloads are used instead and
filtered here; both are cached for 7 days:
  validity: https://search.clinicalgenome.org/kb/gene-validity/download  (CSV)
  dosage:   https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
Terms of use: https://www.clinicalgenome.org/docs/terms-of-use/
"""

from __future__ import annotations

import csv
import io
from typing import Any, Dict, List, Optional

from zebra.core import Outcome
from zebra.http import request
from zebra.sources import record as source_record

VALIDITY_CSV = "https://search.clinicalgenome.org/kb/gene-validity/download"
DOSAGE_TSV = {
    "GRCh38": "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv",
    "GRCh37": "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv",
}
GENE_PAGE = "https://search.clinicalgenome.org/kb/genes/{}"
CACHE_TTL = 7 * 86400
# Order used to sort curations, strongest first.
RANK = {"definitive": 0, "strong": 1, "moderate": 2, "limited": 3, "supportive": 4, "animal model only": 5,
        "no known disease relationship": 6, "disputed": 7, "refuted": 8}
DOSAGE_SCORES = {
    "3": "Sufficient evidence for dosage pathogenicity", "2": "Some evidence for dosage pathogenicity",
    "1": "Little evidence for dosage pathogenicity", "0": "No evidence available",
    "30": "Gene associated with autosomal recessive phenotype", "40": "Dosage sensitivity unlikely",
}


def _match(symbol: Optional[str], hgnc_id: Optional[str], row_symbol: str, row_hgnc: str = "") -> bool:
    if hgnc_id and row_hgnc:
        return row_hgnc.strip().upper() == hgnc_id.strip().upper()
    return bool(symbol) and row_symbol.strip().upper() == symbol.strip().upper()


def parse_validity_csv(text: str, symbol: Optional[str] = None, hgnc_id: Optional[str] = None) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    reader = csv.reader(io.StringIO(text))
    header: Optional[List[str]] = None
    for rec in reader:
        if not rec:
            continue
        if header is None:
            if rec[0].strip().upper() == "GENE SYMBOL":
                header = [h.strip().upper() for h in rec]
            continue
        if rec[0].startswith("+++"):
            continue
        row = dict(zip(header, rec))
        if not _match(symbol, hgnc_id, row.get("GENE SYMBOL", ""), row.get("GENE ID (HGNC)", "")):
            continue
        rows.append({
            "gene": row.get("GENE SYMBOL", "").strip(),
            "hgnc_id": row.get("GENE ID (HGNC)", "").strip() or None,
            "disease": " ".join(row.get("DISEASE LABEL", "").split()),
            "mondo": row.get("DISEASE ID (MONDO)", "").strip() or None,
            "moi": row.get("MOI", "").strip() or None,
            "classification": row.get("CLASSIFICATION", "").strip() or None,
            "date": (row.get("CLASSIFICATION DATE", "").strip() or "")[:10] or None,
            "gcep": row.get("GCEP", "").strip() or None,
            "sop": row.get("SOP", "").strip() or None,
            "url": row.get("ONLINE REPORT", "").strip() or None,
        })
    if header is None:
        raise ValueError("ClinGen validity CSV: header row 'GENE SYMBOL' not found")
    rows.sort(key=lambda r: (RANK.get((r["classification"] or "").lower(), 9), r["disease"]))
    return rows


def validity(symbol: Optional[str] = None, hgnc_id: Optional[str] = None) -> Outcome:
    if not symbol and not hgnc_id:
        raise ValueError("give a gene symbol or HGNC id")
    resp = request(VALIDITY_CSV, source="ClinGen validity", accept="text/csv,*/*", cache_ttl=CACHE_TTL, timeout=120)
    rows = parse_validity_csv(resp.text, symbol=symbol, hgnc_id=hgnc_id)
    label = hgnc_id or symbol
    return Outcome(rows, sources=[source_record("ClinGen gene-disease validity", label, resp,
                                                url=GENE_PAGE.format(hgnc_id) if hgnc_id else VALIDITY_CSV,
                                                note="bulk CSV filtered by gene")])


def parse_dosage_tsv(text: str, symbol: str) -> Optional[Dict[str, Any]]:
    header: Optional[List[str]] = None
    for line in text.splitlines():
        if line.startswith("#"):
            if line.startswith("#Gene Symbol"):
                header = [h.strip() for h in line.lstrip("#").split("\t")]
            continue
        if header is None or not line.strip():
            continue
        cols = line.split("\t")
        if cols[0].strip().upper() != symbol.strip().upper():
            continue
        row = dict(zip(header, cols))

        def pmids(prefix: str) -> List[str]:
            return [row[k].strip() for k in header if k.startswith(prefix) and row.get(k, "").strip()]

        hi, ts = row.get("Haploinsufficiency Score", "").strip(), row.get("Triplosensitivity Score", "").strip()
        return {
            "gene": row.get("Gene Symbol"), "ncbi_gene_id": row.get("Gene ID"), "cytoband": row.get("cytoBand"),
            "location": row.get("Genomic Location"),
            "haploinsufficiency": {"score": hi or None, "description": row.get("Haploinsufficiency Description") or DOSAGE_SCORES.get(hi),
                                   "disease": row.get("Haploinsufficiency Disease ID") or None, "pmids": pmids("Haploinsufficiency PMID")},
            "triplosensitivity": {"score": ts or None, "description": row.get("Triplosensitivity Description") or DOSAGE_SCORES.get(ts),
                                  "disease": row.get("Triplosensitivity Disease ID") or None, "pmids": pmids("Triplosensitivity PMID")},
            "last_evaluated": row.get("Date Last Evaluated") or None,
        }
    if header is None:
        raise ValueError("ClinGen dosage TSV: header '#Gene Symbol' not found")
    return None


def dosage(symbol: str, assembly: str = "GRCh38") -> Outcome:
    url = DOSAGE_TSV.get(assembly)
    if not url:
        raise ValueError("assembly must be GRCh38 or GRCh37")
    resp = request(url, source="ClinGen dosage", accept="text/tab-separated-values,text/plain,*/*",
                   cache_ttl=CACHE_TTL, timeout=120)
    row = parse_dosage_tsv(resp.text, symbol)
    if row is not None:
        row["url"] = f"https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen/clingen_gene.cgi?sym={symbol}"
    return Outcome(row, sources=[source_record("ClinGen dosage sensitivity", symbol, resp, note="bulk TSV filtered by gene")])
