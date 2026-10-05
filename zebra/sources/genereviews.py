"""GeneReviews chapters (NCBI Bookshelf) for a disease (OMIM number, name) or a gene.

Lookup tables are NCBI's own GeneReviews mapping files
(https://ftp.ncbi.nlm.nih.gov/pub/GeneReviews/):
  NBKid_shortname_OMIM.txt                    NBK id -> OMIM numbers the chapter covers
  GRshortname_NBKid_genesymbol_dzname.txt     NBK id -> gene symbols and the chapter's disease name
Chapter metadata (title, authors, last update) comes from E-utilities
(esearch/esummary, db=books); a name with no table hit is searched there too.
The bookshelf HTML is not fetched (it answers scripts with a captcha).
"""

from __future__ import annotations

import html
import re
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from zebra.core import Outcome, attempt
from zebra.http import get_json, request
from zebra.sources import record as source_record

FTP = "https://ftp.ncbi.nlm.nih.gov/pub/GeneReviews"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OMIM_FILE = "NBKid_shortname_OMIM.txt"
GENE_FILE = "GRshortname_NBKid_genesymbol_dzname.txt"


def book_url(nbk: str) -> str:
    return f"https://www.ncbi.nlm.nih.gov/books/{nbk}/"


def parse_omim_map(text: str) -> Dict[str, List[Tuple[str, str]]]:
    """OMIM number -> [(NBK id, shortname)]."""
    out: Dict[str, List[Tuple[str, str]]] = {}
    for ln in text.splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        cols = ln.split("\t")
        if len(cols) >= 3 and cols[2].strip().isdigit():
            out.setdefault(cols[2].strip(), []).append((cols[0].strip(), cols[1].strip()))
    return out


def parse_gene_map(text: str) -> List[Dict[str, str]]:
    """Rows {shortname, nbk, gene, disease} (gene 'Not applicable' for non-genic chapters)."""
    rows = []
    for ln in text.splitlines():
        cols = ln.split("|")
        if len(cols) >= 4 and cols[1].startswith("NBK"):
            rows.append({"shortname": cols[0].strip(), "nbk": cols[1].strip(), "gene": cols[2].strip(),
                         "disease": cols[3].strip()})
    return rows


def _fetch_maps():
    r1 = request(f"{FTP}/{OMIM_FILE}", source="GeneReviews (NCBI FTP)", accept="text/plain", cache_ttl=30 * 86400, timeout=60)
    r2 = request(f"{FTP}/{GENE_FILE}", source="GeneReviews (NCBI FTP)", accept="text/plain", cache_ttl=30 * 86400, timeout=60)
    return parse_omim_map(r1.text), parse_gene_map(r2.text), [
        source_record("GeneReviews OMIM map", OMIM_FILE, r1), source_record("GeneReviews gene map", GENE_FILE, r2)]


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (s or "").lower()).strip()


def match_chapters(omim_map: Dict[str, List[Tuple[str, str]]], gene_rows: List[Dict[str, str]],
                   omim_ids: Iterable[str] = (), genes: Iterable[str] = (), name: Optional[str] = None) -> Dict[str, Dict[str, Any]]:
    """NBK id -> {nbk, title, genes, matched_by[]} from the two NCBI tables."""
    titles: Dict[str, str] = {}
    gene_of: Dict[str, List[str]] = {}
    for r in gene_rows:
        titles.setdefault(r["nbk"], r["disease"])
        if r["gene"] and r["gene"] != "Not applicable":
            gene_of.setdefault(r["nbk"], [])
            if r["gene"] not in gene_of[r["nbk"]]:
                gene_of[r["nbk"]].append(r["gene"])
    hits: Dict[str, Dict[str, Any]] = {}

    def add(nbk: str, why: str, short: Optional[str] = None) -> None:
        h = hits.setdefault(nbk, {"nbk": nbk, "title": titles.get(nbk), "shortname": short, "genes": gene_of.get(nbk, [])[:12],
                                  "url": book_url(nbk), "matched_by": []})
        if why not in h["matched_by"]:
            h["matched_by"].append(why)

    for o in omim_ids:
        num = str(o).upper().replace("OMIM:", "").strip()
        for nbk, short in omim_map.get(num, []):
            add(nbk, f"OMIM:{num}", short)
    for g in genes:
        for r in gene_rows:
            if r["gene"].upper() == str(g).upper():
                add(r["nbk"], f"gene {r['gene']}", r["shortname"])
    if name:  # the chapter's disease name equals the name, or contains it as a phrase of 2+ words
        q = _norm(name)
        multiword = len(q.split()) >= 2 and len(q) >= 10
        if len(q) >= 4:
            for nbk, title in titles.items():
                t = _norm(title)
                if q == t or (multiword and re.search(rf"\b{re.escape(q)}\b", t)):
                    add(nbk, f"name '{name}'")
    return hits


def _esummary_chapters(nbks: Sequence[str]) -> Tuple[Dict[str, Dict[str, Any]], List[Dict[str, Any]]]:
    """Chapter metadata from E-utilities for the given NBK ids (one esearch + one esummary)."""
    if not nbks:
        return {}, []
    term = "gene[book] AND (" + " OR ".join(f"{n}[All Fields]" for n in nbks) + ")"
    es = get_json(f"{EUTILS}/esearch.fcgi", source="NCBI E-utilities",
                  params={"db": "books", "term": term, "retmode": "json", "retmax": 200}, cache_ttl=14 * 86400)
    ids = (es.json().get("esearchresult") or {}).get("idlist") or []
    srcs = [source_record("NCBI Bookshelf esearch", term, es)]
    if not ids:
        return {}, srcs
    sm = get_json(f"{EUTILS}/esummary.fcgi", source="NCBI E-utilities",
                  params={"db": "books", "id": ",".join(ids), "retmode": "json"}, cache_ttl=14 * 86400)
    srcs.append(source_record("NCBI Bookshelf esummary", ",".join(nbks), sm))
    return parse_esummary(sm.json(), nbks), srcs


def parse_esummary(data: Dict[str, Any], nbks: Optional[Sequence[str]] = None) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    res = data.get("result") or {}
    for uid in res.get("uids") or []:
        r = res.get(uid) or {}
        if r.get("rtype") != "chapter" or r.get("book") != "gene":
            continue
        nbk = r.get("chapteraccessionid") or r.get("accessionid")
        if nbks and nbk not in nbks:
            continue
        info = r.get("bookinfo") or ""
        date = re.search(r"<Date>(.*?)</Date>", info)
        authors = re.search(r"<MainRecord>.*?<Contributors>(.*?)</Contributors>", info, re.S)
        out[nbk] = {"nbk": nbk, "title": html.unescape(re.sub(r"<[^>]+>", "", r.get("title") or "")).strip(),
                    "dates": date.group(1) if date else r.get("pubdate"), "authors": authors.group(1) if authors else None,
                    "uid": uid, "url": book_url(nbk)}
    return out


def chapters(omim_ids: Iterable[str] = (), genes: Iterable[str] = (), name: Optional[str] = None,
             limit: int = 5) -> Outcome:
    """GeneReviews chapters covering these OMIM numbers / genes / disease name, most specific first."""
    omim_ids, genes = list(omim_ids), list(genes)
    omim_map, gene_rows, srcs = _fetch_maps()
    hits = match_chapters(omim_map, gene_rows, omim_ids, genes, name)
    warnings: List[str] = []
    if not hits and name:
        term = f'gene[book] AND "{name}"[title]'
        es = get_json(f"{EUTILS}/esearch.fcgi", source="NCBI E-utilities",
                      params={"db": "books", "term": term, "retmode": "json", "retmax": 50}, cache_ttl=14 * 86400)
        srcs.append(source_record("NCBI Bookshelf esearch", term, es))
        ids = (es.json().get("esearchresult") or {}).get("idlist") or []
        if ids:
            sm = get_json(f"{EUTILS}/esummary.fcgi", source="NCBI E-utilities",
                          params={"db": "books", "id": ",".join(ids), "retmode": "json"}, cache_ttl=14 * 86400)
            srcs.append(source_record("NCBI Bookshelf esummary", term, sm))
            q = _norm(name)
            for nbk, meta in parse_esummary(sm.json()).items():
                t = _norm(meta["title"])
                if not (q == t or (len(q.split()) >= 2 and re.search(rf"\b{re.escape(q)}\b", t))):
                    continue  # [title] also searches section/table titles; keep chapters named for it
                hits[nbk] = {"nbk": nbk, "title": meta["title"], "shortname": None, "genes": [], "url": meta["url"],
                             "matched_by": [f"Bookshelf title search '{name}'"]}

    def specificity(h: Dict[str, Any]) -> Tuple[int, int]:
        kinds = h["matched_by"]
        return (-sum(1 for k in kinds if k.startswith("OMIM:")) * 2 - sum(1 for k in kinds if k.startswith("name")),
                len(h.get("genes") or []))

    ordered = sorted(hits.values(), key=specificity)[:limit]
    if ordered:
        got = attempt("GeneReviews chapter metadata (E-utilities)",
                      lambda: _esummary_chapters([h["nbk"] for h in ordered]), warnings)
        if got:
            meta, more = got
            srcs.extend(more)
            for h in ordered:
                m = meta.get(h["nbk"])
                if m:
                    h["title"] = m["title"] or h["title"]
                    h["dates"] = m["dates"]
                    h["authors"] = m["authors"]
    return Outcome({"query": {"omim": omim_ids, "genes": genes, "name": name}, "chapters": ordered},
                   sources=srcs, warnings=warnings)
