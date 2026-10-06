"""zebra lit: literature — Europe PMC search, and papers that mention a gene or variant (LitVar2, PubTator3).

Each source is its own section in its own order; nothing is re-ranked or scored here.

A PMID coming from a tool is not evidence that the paper is about what was
asked for: LitVar2 links papers through normalised variant records, which can
merge spellings numbered on different transcripts, and both services index a
mention anywhere in a paper including supplementary tables. So every paper
from LitVar2 and PubTator3 carries the handle that put it in the list
(`matched_entity`) and what the paper's own title and abstract actually say
(`mentions`). Nothing here classifies a paper's topic; the reader judges.
"""

from __future__ import annotations

import argparse
import math
import re
from typing import Any, Dict, List, Optional, Sequence

from zebra.core import Outcome, UsageError, attempt
from zebra.sources import europepmc, litvar, pubtator

ABSTRACT_EXCERPT = 240  # characters of abstract kept in the output, after the mention check has read all of it
# A LitVar paper list is called out when at least this fraction of the papers whose text could be
# checked never name the gene; the floor keeps 1-of-2 from shouting. Per-paper tags are always shown,
# so a single stray PMID is visible without the warning.
UNMATCHED_WARN_FRACTION = 1.0 / 3
UNMATCHED_WARN_FLOOR = 2

_AA3_RE = re.compile("|".join(sorted(litvar.AA3, key=len, reverse=True)))
_AA1_CLASS = "".join(sorted(litvar.AA1))
_AA1_RE = re.compile(r"(?<![A-Za-z])([" + _AA1_CLASS + r"])(?=\d)|(?<=\d)([" + _AA1_CLASS + r"])(?![A-Za-z])")
_STOPS = ("*", "Ter", "X")


def _label(p: Dict[str, Any]) -> str:
    if p.get("pmid"):
        return f"PMID {p['pmid']}"
    if p.get("pmcid"):
        return p["pmcid"]
    return f"{p.get('source') or ''}:{p.get('id')}"


def _to_one_letter(text: str) -> str:
    """Arg712Ter -> R712Ter (three-letter amino-acid codes only; Ter/X/* are left alone)."""
    return _AA3_RE.sub(lambda m: litvar.AA3[m.group(0)], text)


def _to_three_letter(text: str) -> str:
    """R712* -> Arg712* (a one-letter code touching the position number; X is not one, so it stays)."""
    return _AA1_RE.sub(lambda m: litvar.AA1.get(m.group(1) or m.group(2), m.group(1) or m.group(2)), text)


def _word_re(text: str) -> "re.Pattern":
    """`text` as a whole word, case-insensitively: SCN1A matches `SCN1A-related`, not `SCN1AB`."""
    return re.compile(r"(?<![A-Za-z0-9])" + re.escape(text) + r"(?![A-Za-z0-9])", re.I)


def _variant_forms(variant: Optional[str], rsids: Sequence[Optional[str]] = ()) -> List[str]:
    """Every spelling of the variant to look for in a paper's own words.

    Starts from `litvar.spellings` (what the LitVar query used), adds the
    one-letter/three-letter amino-acid twin of each protein change, completes
    the stop-codon triplet (`*`, `Ter`, `X`), and drops the `p.` prefix — a
    prefixed form contains the bare one, so matching bare covers both.
    p.Arg712* gives Arg712*/Arg712Ter/Arg712X/R712*/R712Ter/R712X. The rsIDs
    of the records that produced the papers are added as further spellings.
    """
    seeds: List[str] = []
    for spelling in litvar.spellings(variant) if variant else []:
        twins = [spelling]
        if spelling.startswith("p."):  # only protein changes have amino-acid codes to convert
            twins += [_to_one_letter(spelling), _to_three_letter(spelling)]
        for form in twins:
            if form not in seeds:
                seeds.append(form)
    forms: List[str] = []
    for seed in seeds:
        stop = litvar.STOP_RE.match(seed)
        spread = [seed[:stop.start("alt")] + s for s in _STOPS] if stop else [seed]
        for form in spread:
            bare = form[2:] if form.startswith("p.") else form
            if bare and bare not in forms:
                forms.append(bare)
    for rsid in rsids:
        if rsid and rsid not in forms:
            forms.append(rsid)
    return forms


def _mentions(gene: Optional[str], forms: Sequence[str], title: Optional[str],
              abstract: Optional[str]) -> Dict[str, Any]:
    """What the paper's own title and abstract say — not what the tool that returned it claims.

    `true` the name is in the title or abstract; `false` it is not and there
    was an abstract to look in; `null` there was nothing to look in (Europe
    PMC holds no abstract for this PMID) or nothing was asked about. A paper
    can still discuss the gene in its full text or supplementary tables, which
    is where both services often found it — `false` is not a verdict on topic.
    """
    text = " ".join(t for t in (title, abstract) if t)
    found = [f for f in forms if _word_re(f).search(text)] if text else []

    def verdict(hit: bool, asked: bool) -> Optional[bool]:
        if not asked:
            return None
        if hit:
            return True
        return False if abstract else None

    return {"gene": verdict(bool(gene and text and _word_re(gene).search(text)), bool(gene)),
            "variant": verdict(bool(found), bool(forms)),
            "variant_forms_found": found}


def _attach_mentions(papers: List[Dict[str, Any]], gene: Optional[str], forms: Sequence[str]) -> None:
    """Read each paper's abstract for mentions, then keep only an excerpt of it."""
    for p in papers:
        abstract = p.pop("abstract", None)
        p["mentions"] = _mentions(gene, forms, p.get("title"), abstract)
        p["abstract_excerpt"] = None if not abstract else (
            abstract[:ABSTRACT_EXCERPT] + ("..." if len(abstract) > ABSTRACT_EXCERPT else ""))


def _mention_counts(papers: Sequence[Dict[str, Any]]) -> Dict[str, int]:
    """How many of the papers name the gene / the variant, and for how many there was no abstract."""
    ms = [p.get("mentions") or {} for p in papers]
    count = lambda key, val: sum(1 for m in ms if m.get(key) is val)  # noqa: E731
    return {"papers_checked": len(ms),
            # not the same as gene_unknown: a paper with no abstract whose title names the gene is a true
            "no_abstract": sum(1 for p in papers if p.get("abstract_excerpt") is None),
            "gene_mentioned": count("gene", True), "gene_not_mentioned": count("gene", False),
            "gene_unknown": count("gene", None),
            "variant_mentioned": count("variant", True), "variant_not_mentioned": count("variant", False),
            "variant_unknown": count("variant", None)}


def _mention_tag(p: Dict[str, Any]) -> str:
    """`[gene+ variant-]` style tag for one paper, or '' for papers nothing was checked against."""
    m = p.get("mentions")
    if not isinstance(m, dict):
        return ""
    bits = [name + ("✓" if val else "✗")
            for name, val in (("gene", m.get("gene")), ("variant", m.get("variant"))) if val is not None]
    if "abstract_excerpt" in p and p.get("abstract_excerpt") is None:
        bits.append("no abstract")
    return "[" + ", ".join(bits) + "]" if bits else ""


def _paper_line(i: int, p: Dict[str, Any]) -> str:
    bits = [_label(p), str(p.get("year") or "-"), p.get("journal") or "-", p.get("authors") or "-", p.get("title") or "(no title)"]
    line = f"{i:>2}. " + " | ".join(bits)
    tag = _mention_tag(p)
    if tag:
        line += " " + tag
    entity = p.get("matched_entity") or {}
    if entity.get("id"):
        line += f" via {entity['id']}"
    types = [t for t in p.get("pubType") or []
             if t.lower() not in ("journal article", "research-article") and not t.lower().startswith("research support")]
    if types:
        line += f" [{', '.join(types[:3])}]"
    if p.get("fullTextUrl"):
        line += f"\n    full text: {p['fullTextUrl']}"
    return line


def _counts_line(counts: Dict[str, int], gene: Optional[str], forms: Sequence[str]) -> str:
    """The section's mention tally, as one line."""
    parts: List[str] = []
    if gene:
        parts.append(f"{counts['gene_mentioned']} mention {gene}, {counts['gene_not_mentioned']} do not")
    if forms:
        parts.append(f"{counts['variant_mentioned']} mention the variant ({'/'.join(forms[:6])})")
    parts.append(f"{counts['no_abstract']} with no abstract to check")
    return f"  of the {counts['papers_checked']} shown: " + "; ".join(parts) + " — title and abstract only"


def _mention_texts(pmids: Sequence[Any], out: Outcome) -> Dict[str, Optional[str]]:
    """Abstracts for PMIDs the mention check needs but has not already fetched; {} if Europe PMC is down."""
    ids = [str(p) for p in pmids if p]
    if not ids:
        return {}
    got = attempt("Europe PMC (abstracts for the mention check)", lambda: europepmc.by_pmids(ids), out.warnings)
    if got is None:
        return {}
    return {str(p["pmid"]): p.get("abstract") for p in out.add(got) if p.get("pmid")}


def _litvar_section(variant: str, gene: Optional[str], limit: int, out: Outcome) -> Optional[Dict[str, Any]]:
    got = attempt("LitVar2", lambda: litvar.lookup(variant, gene), out.warnings)
    if got is None:
        return None
    res = out.add(got)
    # CP1-9: LitVar's rsID-level record merges every allele at the rsID (CFTR p.Gly542Arg, a VUS, came
    # back as p.G542X with 842 papers). A record spelled as ANOTHER allele at the same position is
    # excluded from the papers, and said so; one whose allele cannot be told is kept and flagged.
    text = (variant or "").strip()
    v_p = text if (text.startswith("p.") or ":p." in text) else None
    v_c = text if (text.startswith("c.") or ":c." in text) else None
    excluded: List[Dict[str, Any]] = []
    top = []
    for m in [m for m in res["matches"] if m.get("top")]:
        verdict = litvar.allele_match(m, v_p, v_c)
        m["allele"] = {"same": "this allele", "different": "another allele at the same position",
                       "unknown": "not verified from LitVar's spelling"}[verdict]
        if verdict == "different":
            excluded.append({"litvar_id": m["litvar_id"], "spelled": m.get("hgvs") or m.get("name"),
                             "pmid_count": m.get("pmid_count"), "rsid": m.get("rsid")})
        else:
            top.append(m)
    top = top[:3]
    if excluded:
        out.warnings.append(
            "LitVar2: " + "; ".join(f"{e['litvar_id']} ({e['spelled']}, {e['pmid_count']} papers)" for e in excluded)
            + f" excluded — another allele at the same position as {text}"
            + ("; LitVar merges every allele filed under an rsID" if any(e.get("rsid") for e in excluded) else ""))
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
    by_id = {m["litvar_id"]: m for m in top}
    forms = _variant_forms(variant, [m.get("rsid") for m in top])
    if order:
        meta = attempt("Europe PMC (PMID metadata)", lambda: europepmc.by_pmids(order[:limit]), out.warnings)
        if meta is not None:
            for p in out.add(meta):
                lids = origin.get(str(p.get("pmid")), [])
                p["litvar_ids"] = lids
                rec = by_id.get(lids[0]) if lids else None
                p["matched_entity"] = None if rec is None else {
                    "source": "LitVar2", "id": rec["litvar_id"], "name": rec.get("name"), "rsid": rec.get("rsid"),
                    "via": (rec.get("spellings") or [None])[0]}
                papers.append(p)
    _attach_mentions(papers, gene, forms)
    counts = _mention_counts(papers)
    if gene and counts["gene_not_mentioned"] >= max(
            UNMATCHED_WARN_FLOOR, math.ceil(counts["papers_checked"] * UNMATCHED_WARN_FRACTION)):
        out.warnings.append(
            f"LitVar2: {counts['gene_not_mentioned']} of {counts['papers_checked']} papers do not mention {gene} "
            "in their title or abstract; LitVar links papers through normalised variant records that can merge "
            "different transcript numbering — check each PMID before citing it")
    return {"query": res["query"], "gene": res["gene"], "queries": res["queries"], "matches": res["matches"],
            "excluded_other_allele": excluded,
            "papers_total": len(order), "papers": papers,
            "variant_forms_checked": forms, "mention_counts": counts,
            "note": "papers are those LitVar links to the top match(es); order as LitVar lists them. "
                    "matched_entity is the LitVar record that linked each paper; mentions were checked in the "
                    "Europe PMC title and abstract only, and LitVar indexes a mention anywhere in a paper "
                    "including supplementary tables, so gene false means the title and abstract do not name it, "
                    "not that the paper is off topic"}


def _pubtator_section(variant: Optional[str], gene: Optional[str], lit: Optional[Dict[str, Any]], limit: int,
                      out: Outcome) -> Optional[Dict[str, Any]]:
    entity = None
    via: Optional[str] = None
    tried: List[str] = []
    rsids = [m.get("rsid") for m in (lit or {}).get("matches") or [] if m.get("top")]
    if variant:
        texts: List[str] = []
        for rsid in rsids:
            if rsid and rsid not in texts:
                texts.append(rsid)
        texts.append(f"{gene} {variant}" if gene and gene.upper() not in variant.upper() and ":" not in variant else variant)
        for t in texts[:2]:
            tried.append(t)
            got = attempt("PubTator3 autocomplete", lambda t=t: pubtator.resolve_entity(t, "VARIANT"), out.warnings)
            if got is not None:
                ent = out.add(got)["entity"]
                if ent:
                    entity, via = ent, t
                    break
        if entity is None and gene:
            out.warnings.append(f"PubTator3 has no variant entity for {variant} (tried {', '.join(tried)}): gene-level papers shown")
    if entity is None and gene:
        tried.append(gene)
        got = attempt("PubTator3 autocomplete", lambda: pubtator.resolve_entity(gene, "GENE"), out.warnings)
        if got is not None:
            entity = out.add(got)["entity"]
            via = gene if entity else None
    if entity is None:
        if tried:
            out.warnings.append(f"PubTator3: no entity found for {', '.join(tried)}")
        return None
    got = attempt("PubTator3 search", lambda: pubtator.search(entity["id"], limit), out.warnings)
    if got is None:
        return None
    res = out.add(got)
    hits = res["hits"]
    matched = {"source": "PubTator3", "id": entity["id"], "name": entity.get("name"),
               "rsid": via if via and litvar.RSID_RE.match(via) else None, "via": via}
    forms = _variant_forms(variant, rsids)
    # PubTator carries titles but no abstracts; the mention check needs Europe PMC's.
    texts_by_pmid = _mention_texts([h.get("pmid") for h in hits], out)
    for h in hits:
        h["matched_entity"] = dict(matched)
        h["abstract"] = texts_by_pmid.get(str(h.get("pmid")))
    _attach_mentions(hits, gene, forms)
    return {"entity": entity, "count": res["count"], "hits": hits,
            "variant_forms_checked": forms, "mention_counts": _mention_counts(hits),
            "note": "order as PubTator3 ranks them. matched_entity is the PubTator entity searched (its name can "
                    "be the same change numbered on another transcript); mentions were checked in the Europe PMC "
                    "title and abstract only, and PubTator indexes a mention anywhere in a paper"}


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
                if m.get("top") and m.get("allele") == "another allele at the same position":
                    tag = "EXCLUDED: another allele at the same position"
                elif m.get("top") and m.get("allele"):
                    tag += f" ({m['allele']})"
                sig = f" · dbSNP: {', '.join(m['clinical_significance'])}" if m.get("clinical_significance") else ""
                lines.append(f"  {m['litvar_id']} · {m.get('rsid') or '-'} · {m.get('gene')} {m.get('name')} · "
                             f"{m.get('pmid_count')} papers · {tag}{sig}")
            if not lit["matches"]:
                lines.append("  no LitVar record")
            if lit["papers_total"]:
                lines.append(f"  papers linked to the match(es): {lit['papers_total']} unique PMIDs, showing {len(lit['papers'])}")
                lines.append("  " + _counts_line(lit["mention_counts"], gene, lit["variant_forms_checked"]))
                lines += ["  " + _paper_line(i, p) for i, p in enumerate(lit["papers"], 1)]

    if variant or gene:
        pt = _pubtator_section(variant, gene, lit, limit, out)
        if pt is not None:
            res["pubtator"] = pt
            e = pt["entity"]
            lines.append(f"== PubTator3: {e['id']} ({e.get('name')}, {e.get('description') or e.get('type')}) — "
                         f"{pt['count']} articles, showing {len(pt['hits'])}")
            if pt["hits"]:
                lines.append(_counts_line(pt["mention_counts"], gene, pt["variant_forms_checked"]))
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
