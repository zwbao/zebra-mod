"""zebra lit: literature — Europe PMC search, and papers that mention a gene or variant (LitVar2, PubTator3).

Each source is its own section in its own order; nothing is re-ranked or scored here.
"""

from __future__ import annotations

import argparse
from typing import Any, Dict, List, Optional

from zebra.core import Outcome, UsageError, attempt
from zebra.sources import europepmc, litvar, pubtator


def _label(p: Dict[str, Any]) -> str:
    if p.get("pmid"):
        return f"PMID {p['pmid']}"
    if p.get("pmcid"):
        return p["pmcid"]
    return f"{p.get('source') or ''}:{p.get('id')}"


def _paper_line(i: int, p: Dict[str, Any]) -> str:
    bits = [_label(p), str(p.get("year") or "-"), p.get("journal") or "-", p.get("authors") or "-", p.get("title") or "(no title)"]
    line = f"{i:>2}. " + " | ".join(bits)
    types = [t for t in p.get("pubType") or []
             if t.lower() not in ("journal article", "research-article") and not t.lower().startswith("research support")]
    if types:
        line += f" [{', '.join(types[:3])}]"
    if p.get("fullTextUrl"):
        line += f"\n    full text: {p['fullTextUrl']}"
    return line


def _litvar_section(variant: str, gene: Optional[str], limit: int, out: Outcome) -> Optional[Dict[str, Any]]:
    got = attempt("LitVar2", lambda: litvar.lookup(variant, gene), out.warnings)
    if got is None:
        return None
    res = out.add(got)
    top = [m for m in res["matches"] if m.get("top")][:3]
    order: List[str] = []
    origin: Dict[str, List[str]] = {}
    for m in top:
        pg = attempt(f"LitVar2 publications {m['litvar_id']}", lambda m=m: litvar.pmids(m["litvar_id"]), out.warnings)
        if pg is None:
            continue
        ids = out.add(pg)["pmids"]
        m["pmids"] = ids[:limit]
        for p in ids:
            if p not in origin:
                order.append(p)
                origin[p] = []
            origin[p].append(m["litvar_id"])
    papers: List[Dict[str, Any]] = []
    if order:
        meta = attempt("Europe PMC (PMID metadata)", lambda: europepmc.by_pmids(order[:limit]), out.warnings)
        if meta is not None:
            for p in out.add(meta):
                p["litvar_ids"] = origin.get(str(p.get("pmid")), [])
                papers.append(p)
    return {"query": res["query"], "gene": res["gene"], "queries": res["queries"], "matches": res["matches"],
            "papers_total": len(order), "papers": papers,
            "note": "papers are those LitVar links to the top match(es); order as LitVar lists them"}


def _pubtator_section(variant: Optional[str], gene: Optional[str], lit: Optional[Dict[str, Any]], limit: int,
                      out: Outcome) -> Optional[Dict[str, Any]]:
    entity = None
    tried: List[str] = []
    if variant:
        texts: List[str] = []
        for m in (lit or {}).get("matches") or []:
            if m.get("top") and m.get("rsid") and m["rsid"] not in texts:
                texts.append(m["rsid"])
        texts.append(f"{gene} {variant}" if gene and gene.upper() not in variant.upper() and ":" not in variant else variant)
        for t in texts[:2]:
            tried.append(t)
            got = attempt("PubTator3 autocomplete", lambda t=t: pubtator.resolve_entity(t, "VARIANT"), out.warnings)
            if got is not None:
                ent = out.add(got)["entity"]
                if ent:
                    entity = ent
                    break
        if entity is None and gene:
            out.warnings.append(f"PubTator3 has no variant entity for {variant} (tried {', '.join(tried)}): gene-level papers shown")
    if entity is None and gene:
        tried.append(gene)
        got = attempt("PubTator3 autocomplete", lambda: pubtator.resolve_entity(gene, "GENE"), out.warnings)
        if got is not None:
            entity = out.add(got)["entity"]
    if entity is None:
        if tried:
            out.warnings.append(f"PubTator3: no entity found for {', '.join(tried)}")
        return None
    got = attempt("PubTator3 search", lambda: pubtator.search(entity["id"], limit), out.warnings)
    if got is None:
        return None
    res = out.add(got)
    return {"entity": entity, "count": res["count"], "hits": res["hits"], "note": "order as PubTator3 ranks them"}


def _lit(args: argparse.Namespace) -> Outcome:
    query = (args.query or "").strip()
    gene = (args.gene or "").strip() or None
    variant = (args.variant or "").strip() or None
    abstracts = [a.strip() for a in (args.abstract or []) if a.strip()]
    if not (query or gene or variant or abstracts):
        raise UsageError("give a query, --gene, --variant, or --abstract PMID")
    if args.limit < 1:
        raise UsageError("--limit must be at least 1")
    limit = min(args.limit, 100)
    out = Outcome({}, query={"query": query or None, "gene": gene, "variant": variant, "limit": limit,
                             "sort": args.sort, "abstract": abstracts or None})
    res: Dict[str, Any] = out.result
    lines: List[str] = []

    if query:
        got = attempt("Europe PMC", lambda: europepmc.search(query, limit, args.sort), out.warnings)
        if got is not None:
            sec = res["europepmc"] = out.add(got)
            lines.append(f"== Europe PMC: {query} — {sec['hitCount']} hits, showing {len(sec['hits'])} ({sec['sort']})")
            lines += [_paper_line(i, p) for i, p in enumerate(sec["hits"], 1)] or ["  no hits"]

    lit = None
    if variant:
        lit = _litvar_section(variant, gene, limit, out)
        if lit is not None:
            res["litvar"] = lit
            lines.append(f"== LitVar2: {variant}" + (f" in {gene}" if gene else "") + f" (queried: {'; '.join(lit['queries'])})")
            for m in lit["matches"]:
                tag = "match" if m.get("top") else "LitVar suggestion"
                sig = f" · dbSNP: {', '.join(m['clinical_significance'])}" if m.get("clinical_significance") else ""
                lines.append(f"  {m['litvar_id']} · {m.get('rsid') or '-'} · {m.get('gene')} {m.get('name')} · "
                             f"{m.get('pmid_count')} papers · {tag}{sig}")
            if not lit["matches"]:
                lines.append("  no LitVar record")
            if lit["papers_total"]:
                lines.append(f"  papers linked to the match(es): {lit['papers_total']} unique PMIDs, showing {len(lit['papers'])}")
                lines += ["  " + _paper_line(i, p) for i, p in enumerate(lit["papers"], 1)]

    if variant or gene:
        pt = _pubtator_section(variant, gene, lit, limit, out)
        if pt is not None:
            res["pubtator"] = pt
            e = pt["entity"]
            lines.append(f"== PubTator3: {e['id']} ({e.get('name')}, {e.get('description') or e.get('type')}) — "
                         f"{pt['count']} articles, showing {len(pt['hits'])}")
            lines += [_paper_line(i, p) for i, p in enumerate(pt["hits"], 1)] or ["  no articles"]

    if abstracts:
        res["abstracts"] = []
        for pid in abstracts:
            got = attempt(f"Europe PMC abstract {pid}", lambda pid=pid: europepmc.abstract(pid), out.warnings)
            if got is None:
                continue
            a = out.add(got)
            res["abstracts"].append(a)
            if a.get("found"):
                lines.append(f"== {_label(a)} | {a.get('year')} | {a.get('journal')} | {a.get('authors')}\n{a.get('title')}\n"
                             f"{a.get('abstract') or '(no abstract in Europe PMC)'}"
                             + (f"\nfull text: {a['fullTextUrl']}" if a.get("fullTextUrl") else ""))
            else:
                lines.append(f"== {pid}: not found in Europe PMC")

    out.text = "\n".join(lines) if lines else "no results"
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("lit", help="literature: Europe PMC search; papers that mention a gene or variant (LitVar2, PubTator3)")
    p.add_argument("query", nargs="?", help='Europe PMC query, its syntax passed through: "Dravet syndrome" AND (case report)')
    p.add_argument("--gene", help="gene symbol (HGNC)")
    p.add_argument("--variant", help="rsID, HGVS with transcript (NM_...:c.) or protein change (p.Arg712*)")
    p.add_argument("--limit", type=int, default=15, help="papers per section (default 15)")
    p.add_argument("--sort", choices=sorted(europepmc.SORTS), default="relevance", help="Europe PMC order (default relevance)")
    p.add_argument("--abstract", nargs="+", metavar="PMID", help="also fetch title + abstract for these PMIDs/PMCIDs")
    p.set_defaults(func=_lit)
