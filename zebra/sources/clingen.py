"""ClinGen: gene–disease validity and dosage sensitivity (haploinsufficiency / triplosensitivity).

ClinGen's per-gene search API ignores its query string and returns every
curation (1.9 MB / 1.5 MB), so the public bulk downloads are used instead and
filtered here; all are cached for 7 days:
  validity: https://search.clinicalgenome.org/kb/gene-validity/download  (CSV)
  dosage, genes:   https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv   (and _GRCh37)
  dosage, regions: https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv (and _GRCh37)
The region list is ClinGen's curation of genomic *regions* — the recurrent
microdeletion/duplication regions (22q11.2, 7q11.23, 15q11-q13 …) and the
population (benign) CNV regions — which a gene-by-gene look-up never sees.
Terms of use: https://www.clinicalgenome.org/docs/terms-of-use/

Both TSVs are validated before they are trusted or cached (B-P1-5): the header
row must be present and the body must carry a plausible number of data rows, so
a 200 maintenance page or a truncated download becomes one named failure instead
of "no dosage-sensitive genes" for a week.

Coordinates in the TSVs ("chr22:18924718-21111383") are read as 1-based and
inclusive, the UCSC position convention they are written in; ClinGen's own BED
release of the recurrent regions carries the same numbers.
"""

from __future__ import annotations

import csv
import io
import re
from typing import Any, Dict, List, Optional, Tuple

from zebra.core import Outcome
from zebra.http import Response, request
from zebra.sources import record as source_record
from zebra.sources import validated_text

VALIDITY_CSV = "https://search.clinicalgenome.org/kb/gene-validity/download"
DOSAGE_TSV = {
    "GRCh38": "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv",
    "GRCh37": "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv",
}
REGION_TSV = {
    "GRCh38": "https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv",
    "GRCh37": "https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh37.tsv",
}
GENE_PAGE = "https://search.clinicalgenome.org/kb/genes/{}"
DOSAGE_GENE_PAGE = "https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen/clingen_gene.cgi?sym={}"
REGION_PAGE = "https://search.clinicalgenome.org/kb/gene-dosage/region/{}"  # checked 2026-10-06: 200 for ISCA-37446
CACHE_TTL = 7 * 86400
# Order used to sort curations, strongest first.
RANK = {"definitive": 0, "strong": 1, "moderate": 2, "limited": 3, "supportive": 4, "animal model only": 5,
        "no known disease relationship": 6, "disputed": 7, "refuted": 8}
DOSAGE_SCORES = {
    "3": "Sufficient evidence for dosage pathogenicity", "2": "Some evidence for dosage pathogenicity",
    "1": "Little evidence for dosage pathogenicity", "0": "No evidence available",
    "30": "Gene associated with autosomal recessive phenotype", "40": "Dosage sensitivity unlikely",
}
# The two bulk dosage tables: the header row each must carry, and the fewest data
# rows a genuine file has. On 2026-10-06 the files held 1532/1534 gene rows and
# 513/513 region rows; a body far below that is an error page or a cut download,
# never "fewer curations".
TABLES = {
    "gene": {"urls": DOSAGE_TSV, "header": "#Gene Symbol", "min_rows": 1400, "label": "ClinGen dosage (genes)"},
    "region": {"urls": REGION_TSV, "header": "#ISCA ID", "min_rows": 450, "label": "ClinGen dosage (regions)"},
}
_MIN_COLUMNS = 13  # both tables carry 23 columns; the scores are in columns 5 and 13
_LOC_RE = re.compile(r"^chr([0-9]{1,2}|X|Y|M|MT):(\d+)-(\d+)$", re.I)


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


def validity_problem(text: str) -> Optional[str]:
    """Why a body is not ClinGen's gene-validity CSV, or None when it is."""
    if not text or not text.strip():
        return "empty body"
    if text.lstrip()[:1] == "<":
        return "HTML page where the ClinGen validity CSV was expected"
    if "GENE SYMBOL" not in text:
        return "no 'GENE SYMBOL' header row: not the ClinGen validity CSV"
    return None


def validity(symbol: Optional[str] = None, hgnc_id: Optional[str] = None) -> Outcome:
    if not symbol and not hgnc_id:
        raise ValueError("give a gene symbol or HGNC id")

    def fetch(refresh: bool = False) -> Response:
        # refresh=True skips the cache read and still writes the good answer back (E-6):
        # cache_ttl=0 would read past a bad entry and leave it in place for the whole TTL
        return request(VALIDITY_CSV, source="ClinGen validity", accept="text/csv,*/*", cache_ttl=CACHE_TTL,
                       refresh=refresh, timeout=120, validate=validity_problem)

    used = [fetch()]

    def refetch() -> Response:
        used.append(fetch(refresh=True))
        return used[-1]

    text = validated_text(used[0], "ClinGen validity", must_contain="GENE SYMBOL", refetch=refetch)
    resp = used[-1]
    rows = parse_validity_csv(text, symbol=symbol, hgnc_id=hgnc_id)
    label = hgnc_id or symbol
    return Outcome(rows, sources=[source_record("ClinGen gene-disease validity", label, resp,
                                                url=GENE_PAGE.format(hgnc_id) if hgnc_id else VALIDITY_CSV,
                                                note="bulk CSV filtered by gene")])


# ------------------------------------------------------- dosage bulk tables

def table_problem(text: str, kind: str) -> Optional[str]:
    """Why `text` is not ClinGen's `kind` ("gene"/"region") dosage table, or None when it is."""
    spec = TABLES[kind]
    if not text or not text.strip():
        return "empty body"
    if text.lstrip()[:1] == "<":
        return "HTML page where the ClinGen TSV was expected"
    lines = text.splitlines()
    if not any(line.startswith(spec["header"] + "\t") for line in lines[:20]):
        return f"no {spec['header']!r} header row: not the ClinGen {kind} curation list"
    rows = sum(1 for line in lines if line and not line.startswith("#") and line.count("\t") >= _MIN_COLUMNS - 1)
    last = next((line for line in reversed(lines) if line.strip()), "")
    header = next(line for line in lines[:20] if line.startswith(spec["header"] + "\t"))
    if not last.startswith("#") and last.count("\t") < header.count("\t"):
        return "the last row is cut short: a truncated download"
    if rows < spec["min_rows"]:
        return (f"only {rows} data rows (a complete file has well over {spec['min_rows']}): "
                "a truncated download or an error page")
    return None


def fetch_table(kind: str, assembly: str = "GRCh38") -> Response:
    """The bulk TSV for `kind`, fetched once, validated before it is cached or used."""
    spec = TABLES[kind]
    url = spec["urls"].get(assembly)
    if not url:
        raise ValueError("assembly must be GRCh38 or GRCh37")

    def get(refresh: bool = False) -> Response:
        return request(url, source=spec["label"], accept="text/tab-separated-values,text/plain,*/*",
                       cache_ttl=CACHE_TTL, timeout=120, refresh=refresh,
                       validate=lambda t: table_problem(t, kind))

    resp = get()
    # `request` already refuses (and never caches) a body `table_problem` rejects;
    # this second gate covers an entry written by an older version without the check
    validated_text(resp, spec["label"], must_contain=spec["header"], refetch=lambda: get(refresh=True))
    return resp


def file_date(text: str) -> Optional[str]:
    """The date line ClinGen writes at the top of each table ('#06 Oct,2026')."""
    for line in text.splitlines()[:4]:
        got = re.match(r"^#\s*(\d{1,2} \w{3},\s*\d{4})\s*$", line)
        if got:
            return got.group(1)
    return None


def _location(text: str) -> Optional[Tuple[str, int, int]]:
    got = _LOC_RE.match((text or "").strip())
    if not got:
        return None
    chrom = got.group(1).upper()
    return ("MT" if chrom == "M" else chrom), int(got.group(2)), int(got.group(3))


def _score(raw: Optional[str]) -> Optional[str]:
    """A score as ClinGen writes it ('3', '40', 'Not yet evaluated'); empty means not given."""
    raw = (raw or "").strip()
    return raw or None


def _rows(text: str, header_name: str) -> List[Dict[str, str]]:
    header: Optional[List[str]] = None
    out: List[Dict[str, str]] = []
    for line in text.splitlines():
        if line.startswith("#"):
            if line.startswith(header_name):
                header = [h.strip() for h in line.lstrip("#").split("\t")]
            continue
        if header is None or not line.strip():
            continue
        out.append(dict(zip(header, line.split("\t"))))
    if header is None:
        raise ValueError(f"ClinGen dosage TSV: header {header_name!r} not found")
    return out


def _dosage_block(row: Dict[str, str], prefix: str) -> Dict[str, Any]:
    score = _score(row.get(f"{prefix} Score"))
    return {"score": score,
            "description": (row.get(f"{prefix} Description") or "").strip() or DOSAGE_SCORES.get(score or ""),
            "disease": (row.get(f"{prefix} Disease ID") or "").strip() or None,
            "pmids": [row[k].strip() for k in row if k.startswith(f"{prefix} PMID") and (row.get(k) or "").strip()]}


def parse_gene_table(text: str) -> List[Dict[str, Any]]:
    """Every gene curation in ClinGen_gene_curation_list_<build>.tsv, with its location parsed."""
    out = []
    for row in _rows(text, "#Gene Symbol"):
        loc = _location(row.get("Genomic Location", ""))
        out.append({
            "gene": (row.get("Gene Symbol") or "").strip(), "ncbi_gene_id": (row.get("Gene ID") or "").strip() or None,
            "cytoband": (row.get("cytoBand") or "").strip() or None,
            "location": (row.get("Genomic Location") or "").strip() or None,
            "chrom": loc[0] if loc else None, "start": loc[1] if loc else None, "end": loc[2] if loc else None,
            "haploinsufficiency": _dosage_block(row, "Haploinsufficiency"),
            "triplosensitivity": _dosage_block(row, "Triplosensitivity"),
            "last_evaluated": (row.get("Date Last Evaluated") or "").strip() or None,
        })
    return out


def parse_region_table(text: str) -> List[Dict[str, Any]]:
    """Every region curation in ClinGen_region_curation_list_<build>.tsv, with its location parsed."""
    out = []
    for row in _rows(text, "#ISCA ID"):
        isca = (row.get("ISCA ID") or "").strip()
        loc = _location(row.get("Genomic Location", ""))
        name = " ".join((row.get("ISCA Region Name") or "").split())
        out.append({
            "isca_id": isca, "name": name, "cytoband": (row.get("cytoBand") or "").strip() or None,
            "location": (row.get("Genomic Location") or "").strip() or None,
            "chrom": loc[0] if loc else None, "start": loc[1] if loc else None, "end": loc[2] if loc else None,
            "named_genes": region_named_genes(name),
            "haploinsufficiency": _dosage_block(row, "Haploinsufficiency"),
            "triplosensitivity": _dosage_block(row, "Triplosensitivity"),
            "last_evaluated": (row.get("Date Last Evaluated") or "").strip() or None,
            "url": REGION_PAGE.format(isca) if isca else None,
        })
    return out


def region_named_genes(name: str) -> List[str]:
    """The genes a region's name says it includes: '(includes RBM8A and GJA5)' -> ['RBM8A', 'GJA5']."""
    got = re.search(r"\(includes ([^)]*)\)", name or "")
    if not got:
        return []
    return [g for g in re.split(r",\s*|\s+and\s+|\s+", got.group(1).strip()) if re.match(r"^[A-Za-z0-9-]+$", g)
            and g.lower() != "and"]


def parse_dosage_tsv(text: str, symbol: str) -> Optional[Dict[str, Any]]:
    """One gene's record from the gene curation list (the shape `dosage()` returns)."""
    for row in parse_gene_table(text):
        if row["gene"].upper() != symbol.strip().upper():
            continue
        return {"gene": row["gene"], "ncbi_gene_id": row["ncbi_gene_id"], "cytoband": row["cytoband"],
                "location": row["location"], "haploinsufficiency": row["haploinsufficiency"],
                "triplosensitivity": row["triplosensitivity"], "last_evaluated": row["last_evaluated"]}
    return None


def dosage(symbol: str, assembly: str = "GRCh38") -> Outcome:
    resp = fetch_table("gene", assembly)
    row = parse_dosage_tsv(resp.text, symbol)
    if row is not None:
        row["url"] = DOSAGE_GENE_PAGE.format(symbol)
    return Outcome(row, sources=[source_record("ClinGen dosage sensitivity", symbol, resp, note="bulk TSV filtered by gene")])


# ------------------------------------------------------------- overlaps

def overlap(chrom: str, start: int, end: int, item_start: int, item_end: int) -> Dict[str, Any]:
    """How the interval [start, end] and a curated item [item_start, item_end] (same chromosome) overlap."""
    lo, hi = max(start, item_start), min(end, item_end)
    bp = max(0, hi - lo + 1)
    item_len = item_end - item_start + 1
    cnv_len = end - start + 1
    contains = start <= item_start and end >= item_end
    inside = start >= item_start and end <= item_end
    relation = ("identical" if contains and inside else "cnv_contains_it" if contains
                else "cnv_within_it" if inside else "partial")
    out = {"overlap_bp": bp, "item_fraction_covered": round(bp / item_len, 4) if item_len > 0 else None,
           "cnv_fraction_in_item": round(bp / cnv_len, 4) if cnv_len > 0 else None, "relation": relation}
    if relation == "partial" or relation == "cnv_within_it":
        out["item_bp_outside_cnv"] = {"before": max(0, start - item_start), "after": max(0, item_end - end)}
    return out


def regions_overlapping(regions: List[Dict[str, Any]], chrom: str, start: int, end: int) -> List[Dict[str, Any]]:
    """The curated regions that share at least one base with chrom:start-end, each with its overlap."""
    hits = []
    for r in regions:
        if r["chrom"] != chrom or r["start"] is None:
            continue
        if r["start"] <= end and r["end"] >= start:
            hits.append(dict(r, overlap=overlap(chrom, start, end, r["start"], r["end"])))
    return hits


def genes_overlapping(genes: List[Dict[str, Any]], chrom: str, start: int, end: int) -> List[Dict[str, Any]]:
    """The curated genes whose ClinGen location shares at least one base with chrom:start-end."""
    hits = []
    for g in genes:
        if g["chrom"] != chrom or g["start"] is None:
            continue
        if g["start"] <= end and g["end"] >= start:
            hits.append(dict(g, overlap=overlap(chrom, start, end, g["start"], g["end"])))
    return hits
