"""zebra phenotype rank: phenotype-driven differential from several independent rankers.

Sources (each ranks on its own; scores are never added up across sources):
  local          offline Resnik best-match average over the HPO release (zebra.hpo_local);
                 honours excluded terms (penalty)
  monarch        Monarch semantic similarity search (diseases and genes); no excluded terms
  pubcasefinder  PubCaseFinder ranked lists (OMIM, Orphanet, genes); no excluded terms

Agreement is reported as `consensus`: diseases in the top N of two or more sources.
Ids from different namespaces are joined only through exact mappings: Monarch's
SSSOM exactMatch (OMIM/Orphanet -> MONDO), with obsolete MONDO classes followed to
their replacement in OLS, and the MONDO xrefs Monarch returns with each hit. Every
consensus row says which ids each source used and how they were joined.
"""

from __future__ import annotations

import argparse
import re
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from zebra import hpo_local
from zebra.core import Outcome, UsageError, attempt

SOURCES = ("local", "monarch", "pubcasefinder")
HPO_RE = re.compile(r"^HP:\d{7}$")
MAX_TOP = 50
PCF_BUDGET_S = 75.0  # stop asking PubCaseFinder for more lists after this (the tool times out at 180 s)


# ---------------------------------------------------------------- inputs


def _split_ids(values: Optional[Sequence[str]]) -> List[str]:
    out: List[str] = []
    for v in values or []:
        for part in re.split(r"[,\s]+", v.strip()):
            if part:
                out.append(part.upper().replace("HP_", "HP:"))
    return out


def gather_terms(args: argparse.Namespace) -> Tuple[List[str], List[str], Optional[str]]:
    present = _split_ids(getattr(args, "present", None)) + _split_ids(getattr(args, "hpo", None))
    excluded = _split_ids(getattr(args, "exclude", None))
    case_path = getattr(args, "case", None)
    used_case = None
    if args.from_case or (not present and case_path):
        if not case_path:
            raise UsageError("--from-case needs an active case (--case <dir> or $ZEBRA_CASE)")
        from zebra import case as case_mod

        try:
            data = case_mod.load(case_path)
        except case_mod.CaseError as err:
            raise UsageError(str(err)) from None
        present += [p["id"] for p in data["phenotypes"] if p.get("status", "present") == "present"]
        excluded += [p["id"] for p in data["phenotypes"] if p.get("status") == "excluded"]
        used_case = case_path
    bad = [t for t in present + excluded if not HPO_RE.match(t)]
    if bad:
        raise UsageError(f"not HPO ids (HP:0000000): {', '.join(bad)}")
    present = list(dict.fromkeys(present))
    excluded = [t for t in dict.fromkeys(excluded) if t not in present]
    if not present:
        raise UsageError("give present phenotypes (--present HP:...), or --from-case with an active case that has some")
    return present, excluded, used_case


def parse_sources(text: Optional[str]) -> List[str]:
    if not text:
        return list(SOURCES)
    chosen = [s.strip().lower() for s in text.split(",") if s.strip()]
    bad = [s for s in chosen if s not in SOURCES]
    if bad:
        raise UsageError(f"unknown source(s) {', '.join(bad)}; choose from {', '.join(SOURCES)}")
    return list(dict.fromkeys(chosen))


# ---------------------------------------------------------------- per-source runs


def run_local(idx: Any, present: Sequence[str], excluded: Sequence[str], top: int) -> Outcome:
    from zebra.sources import record as source_record

    res = hpo_local.rank(idx, present, excluded, top=top, db=("OMIM", "ORPHA"))
    diseases = []
    for i, d in enumerate(res["diseases"], 1):
        diseases.append({
            "rank": i, "id": d["disease"], "name": d["name"], "score": d["score"], "relative": d["relative"],
            "exact_matches": sum(1 for m in d["matches"] if m["exact"]), "of": len(d["matches"]),
            "matched": _compact_matches((m["query"], m["exact"], m["via"], m["via_label"], None) for m in d["matches"]),
            "excluded_hits": d["excluded_hits"], "genes": [g for g in d["genes"] if g not in ("-", "")][:8],
        })
    real = [g for g in res["genes"] if g["gene"] not in ("-", "")]  # '-' = annotation without a gene
    genes = [{"rank": i, "symbol": g["gene"], "score": g["score"], "via": g["via"]} for i, g in enumerate(real[:top], 1)]
    block = {"method": res["method"], "hpo_version": res["hpo_version"], "excluded_supported": True,
             "diseases": diseases, "genes": genes, "notes": res["notes"],
             "query_labels": {q["id"]: q["label"] for q in res["query"] + res["excluded"]}}
    return Outcome(block, sources=[source_record("HPO annotations (phenotype.hpoa)", res["hpo_version"],
                                                 url="https://hpo.jax.org/data/annotations",
                                                 note="local Resnik ranking, zebra.hpo_local")])


def _compact_matches(rows) -> Dict[str, Any]:
    """{'exact': [query ids matched as such], 'partial': {query: 'disease term via common ancestor (label)'}}."""
    exact: List[str] = []
    partial: Dict[str, str] = {}
    for query, is_exact, via, via_label, target in rows:
        if is_exact:
            exact.append(query)
        elif via:
            partial[query] = (f"{target} " if target and target != via else "") + f"via {via} {via_label or ''}".rstrip()
        else:
            partial[query] = "no match"
    return {"exact": exact, "partial": partial}


def run_monarch(present: Sequence[str], top: int, warnings: List[str]) -> Outcome:
    from zebra.sources import monarch

    out = Outcome({"method": "Monarch semsim search (ancestor information content, bidirectional); "
                             "pages: https://monarchinitiative.org/<id>",
                   "excluded_supported": False, "diseases": [], "genes": []})
    dis = attempt("Monarch semsim (diseases)", lambda: monarch.semsim_rank(present, limit=top), warnings)
    if dis is not None:
        out.add(dis)
        out.result["diseases"] = [_compact_monarch(h) for h in dis.result["hits"]]
        out.result["query_labels"] = {m["query"]: m.get("query_label") for h in dis.result["hits"] for m in h["matched"]}
    gen = attempt("Monarch semsim (genes)", lambda: monarch.semsim_rank(present, limit=top, group="Human Genes"), warnings)
    if gen is not None:
        out.add(gen)
        out.result["genes"] = [{"rank": h["rank"], "symbol": h.get("symbol") or h.get("name"), "id": h["id"],
                                "score": h["score"]} for h in gen.result["hits"]]
    if dis is None and gen is None:
        out.result = None
    return out


def _compact_monarch(h: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "rank": h["rank"], "id": h["id"], "name": h["name"], "score": h["score"],
        "exact_matches": sum(1 for m in h["matched"] if m["exact"]), "of": len(h["matched"]),
        "matched": _compact_matches((m["query"], m["exact"], m["via"], m["via_label"], m["match"]) for m in h["matched"]),
        "xref": [x for x in h["xref"] if x.split(":")[0] in ("OMIM", "ORPHA")],
    }


def run_pubcasefinder(present: Sequence[str], top: int, warnings: List[str]) -> Outcome:
    """Three requests, one after another (DBCLS asks for <= 10 requests a minute)."""
    from zebra.sources import pubcasefinder as pcf

    out = Outcome({"method": "PubCaseFinder ranked lists (pcf_get_ranked_list)", "excluded_supported": False,
                   "diseases_omim": [], "diseases_orphanet": [], "genes": []})
    got_any = False
    start = time.monotonic()
    for target, key in (("omim", "diseases_omim"), ("orphanet", "diseases_orphanet"), ("gene", "genes")):
        if time.monotonic() - start > PCF_BUDGET_S:
            warnings.append(f"PubCaseFinder ({target}) skipped: the first lists took over {PCF_BUDGET_S:.0f} s")
            continue
        got = attempt(f"PubCaseFinder ({target})", lambda t=target: pcf.ranked(present, t, top), warnings)
        if got is not None:
            got_any = True
            out.add(got)
            out.result[key] = got.result["hits"]
    if not got_any:
        out.result = None
    return out


# ---------------------------------------------------------------- consensus


def _ranked_items(per: Dict[str, Any]) -> List[Dict[str, Any]]:
    """(source, list, id, name, rank) for every disease hit in every source."""
    items = []
    if per.get("local"):
        for d in per["local"]["diseases"]:
            items.append({"source": "local", "list": "local", "id": d["id"], "name": d["name"], "rank": d["rank"]})
    if per.get("monarch"):
        for d in per["monarch"]["diseases"]:
            items.append({"source": "monarch", "list": "monarch", "id": d["id"], "name": d["name"], "rank": d["rank"],
                          "xref": d.get("xref") or []})
    if per.get("pubcasefinder"):
        for lst in ("diseases_omim", "diseases_orphanet"):
            for d in per["pubcasefinder"][lst]:
                items.append({"source": "pubcasefinder", "list": lst.replace("diseases_", "pubcasefinder "),
                              "id": d["id"], "name": d["name"], "rank": d["rank"]})
    return items


def resolve_keys(items: List[Dict[str, Any]], mapping_rows: Optional[List[Dict[str, Any]]],
                 replaced: Optional[Dict[str, Optional[str]]]) -> Dict[str, Tuple[str, str]]:
    """raw id -> (group key, how it was linked). Exact links only; MONDO preferred as the key.

    mapping_rows: Monarch exactMatch rows (subject MONDO, object OMIM/ORPHA); replaced: obsolete MONDO -> current.
    """
    replaced = replaced or {}
    to_mondo: Dict[str, List[Tuple[str, str]]] = {}
    for m in mapping_rows or []:
        if m.get("predicate") not in (None, "skos:exactMatch"):
            continue
        mondo, obj = m.get("subject"), m.get("object")
        if not (mondo and obj and mondo.startswith("MONDO:")):
            continue
        obsolete = (m.get("subject_label") or "").lower().startswith("obsolete")
        if obsolete:
            cur = replaced.get(mondo)
            if not cur:
                continue
            how = f"{obj} =exactMatch (Monarch)=> {mondo} (obsolete) =replaced by (OLS)=> {cur}"
            to_mondo.setdefault(obj, []).append((cur, how))
        else:
            to_mondo.setdefault(obj, []).append((mondo, f"{obj} =exactMatch (Monarch)=> {mondo}"))
    # the MONDO xrefs that came back with each Monarch hit are the source's own exact links
    monarch_ids = {it["id"] for it in items if it["source"] == "monarch"}
    for it in items:
        if it["source"] == "monarch":
            for x in it.get("xref") or []:
                to_mondo.setdefault(x, []).append((it["id"], f"{x} =xref of {it['id']} (Monarch)"))
    keys: Dict[str, Tuple[str, str]] = {}
    for it in items:
        rid = it["id"]
        if rid in keys:
            continue
        if rid.startswith("MONDO:"):
            keys[rid] = (rid, "")
            continue
        cands = to_mondo.get(rid) or []
        if not cands:
            keys[rid] = (rid, "")
            continue
        uniq: Dict[str, str] = {}
        for mondo, how in cands:
            uniq.setdefault(mondo, how)
        if len(uniq) == 1:
            mondo, how = next(iter(uniq.items()))
        else:  # several MONDO classes: prefer one Monarch itself ranked
            pick = sorted(uniq, key=lambda m: (m not in monarch_ids, m))[0]
            mondo, how = pick, uniq[pick] + f" (also maps to {', '.join(sorted(set(uniq) - {pick}))})"
        keys[rid] = (mondo, how)
    return keys


def build_consensus(items: List[Dict[str, Any]], keys: Dict[str, Tuple[str, str]]) -> List[Dict[str, Any]]:
    groups: Dict[str, Dict[str, Any]] = {}
    for it in items:
        key, how = keys.get(it["id"], (it["id"], ""))
        g = groups.setdefault(key, {"key": key, "names": {}, "hits": [], "links": []})
        g["names"].setdefault(it["source"], it["name"])
        g["hits"].append({"source": it["source"], "list": it["list"], "id": it["id"], "rank": it["rank"]})
        if how and how not in g["links"]:
            g["links"].append(how)
    rows = []
    for g in groups.values():
        best: Dict[str, int] = {}
        for h in g["hits"]:
            best[h["source"]] = min(best.get(h["source"], 10 ** 6), h["rank"])
        if len(best) < 2:
            continue
        name = g["names"].get("monarch") or g["names"].get("local") or g["names"].get("pubcasefinder")
        rows.append({
            "name": name, "key": g["key"], "n_sources": len(best),
            "best_rank": {s: best[s] for s in SOURCES if s in best},
            "hits": sorted(g["hits"], key=lambda h: (SOURCES.index(h["source"]), h["rank"])),
            "linked_by": g["links"] or ["same id"],
        })
    rows.sort(key=lambda r: (-r["n_sources"], max(r["best_rank"].values()), min(r["best_rank"].values()), r["key"]))
    for i, r in enumerate(rows, 1):
        r["order"] = i
    return rows


def gene_consensus(per: Dict[str, Any]) -> List[Dict[str, Any]]:
    ranks: Dict[str, Dict[str, int]] = {}
    for s in SOURCES:
        block = per.get(s)
        if not block:
            continue
        for g in block.get("genes") or []:
            sym = (g.get("symbol") or "").upper()
            if sym:
                ranks.setdefault(sym, {})
                ranks[sym][s] = min(ranks[sym].get(s, 10 ** 6), g["rank"])
    rows = [{"symbol": sym, "n_sources": len(r), "rank": {s: r[s] for s in SOURCES if s in r}}
            for sym, r in ranks.items() if len(r) >= 2]
    rows.sort(key=lambda r: (-r["n_sources"], max(r["rank"].values()), r["symbol"]))
    return rows


def link_ids(items: List[Dict[str, Any]], warnings: List[str]) -> Tuple[Dict[str, Tuple[str, str]], List[Outcome]]:
    """Fetch the exact mappings needed to join OMIM/ORPHA hits with MONDO hits."""
    from zebra.sources import monarch, ols

    used: List[Outcome] = []
    raw = sorted({it["id"] for it in items if it["id"].split(":")[0] in ("OMIM", "ORPHA")})
    has_mondo = any(it["id"].startswith("MONDO:") for it in items)
    rows: Optional[List[Dict[str, Any]]] = None
    replaced: Dict[str, Optional[str]] = {}
    if raw and has_mondo:
        rows = []
        for i in range(0, len(raw), 60):
            got = attempt("Monarch mappings (id reconciliation)", lambda c=raw[i:i + 60]: monarch.mappings(object_ids=c), warnings)
            if got is not None:
                used.append(got)
                rows.extend(got.result)
        obsolete = sorted({m["subject"] for m in rows if (m.get("subject_label") or "").lower().startswith("obsolete")})
        if obsolete:
            got = attempt("OLS (obsolete MONDO -> replacement)", lambda: ols.resolve_obsolete(obsolete), warnings)
            if got is not None:
                used.append(got)
                replaced = got.result
    return resolve_keys(items, rows, replaced), used


# ---------------------------------------------------------------- command


def _rank(args: argparse.Namespace) -> Outcome:
    present, excluded, used_case = gather_terms(args)
    sources = parse_sources(args.sources)
    top = max(1, min(int(args.top), MAX_TOP))
    warnings: List[str] = []
    notes: List[str] = []
    if args.top > MAX_TOP:
        warnings.append(f"--top capped at {MAX_TOP}")

    idx = None
    try:
        idx = hpo_local.load()
    except hpo_local.HpoDataMissing as err:
        if "local" in sources:
            warnings.append(f"local ranking skipped: {err}. Run `zebra hpo fetch` once (~80 MB) to enable it.")
    labels: Dict[str, Optional[str]] = {}
    remote_terms = list(present)
    if idx is not None:  # send current ids to the web services (alt/obsolete ids would match nothing)
        remote_terms = []
        for t in present:
            pid, note = idx.primary(t)
            if note:
                notes.append(note)
            if pid and pid not in remote_terms:
                remote_terms.append(pid)
        labels = {t: idx.names.get(idx.primary(t)[0] or t) for t in present + excluded}

    # F14: a web source cannot be queried with an empty term set. Terms the
    # local release does not know, and obsolete terms with no replacement, are
    # dropped from `remote_terms`; when every term is dropped, PubCaseFinder
    # and Monarch were still called with an empty termset and the answer came
    # back `ok: true` with `diseases: 0`, which reads as "nothing matches this
    # patient" rather than "none of these ids is usable".
    if idx is None and [s for s in sources if s in ("monarch", "pubcasefinder")]:
        # Without the local release there is nothing to resolve ids against, so
        # each one is checked once against the HPO API (cached 30 days) rather
        # than sent to a ranker that answers an unknown id with "no match".
        from zebra.sources import hpo as hpo_src

        checked = []
        for t in remote_terms:
            got = attempt(f"HPO {t}", lambda t=t: hpo_src.term(t, prefer_local=False), warnings)
            if got is None:  # the API could not be reached: keep the term, say so
                checked.append(t)
                continue
            outcome_note = (got.result or {}).get("note")
            if (got.result or {}).get("name") is None and outcome_note == "not found in HPO":
                continue
            replacement = (got.result or {}).get("replaced_by")
            checked.append(replacement or t)
            if replacement and replacement != t:
                notes.append(f"{t} is obsolete in HPO; used {replacement} (HPO API)")
        remote_terms = list(dict.fromkeys(checked))

    dropped = [t for t in present if t not in remote_terms]
    if dropped:
        warnings.append("not sent to the web sources, because this HPO release does not resolve them to a current "
                        f"term: {', '.join(dropped)}")
    remote_sources = [s for s in sources if s in ("monarch", "pubcasefinder")]
    if remote_sources and not remote_terms:
        local_possible = "local" in sources and idx is not None
        raise UsageError(
            f"none of the {len(present)} phenotype(s) given is a term this HPO release can resolve "
            f"({', '.join(present)}): {', '.join(remote_sources)} cannot be queried with an empty set, and an "
            "empty query would come back as 'no disease matches'. Check the ids with `zebra hpo term <id>`, or "
            "search for the phenotype with `zebra hpo search <phrase>`."
            + ("" if local_possible else " Run `zebra hpo fetch` once to get the HPO release."))

    per: Dict[str, Any] = {}
    outcome = Outcome(None)
    jobs: Dict[str, Callable[[], Outcome]] = {}
    job_warn: Dict[str, List[str]] = {s: [] for s in SOURCES}
    if "monarch" in sources:
        jobs["monarch"] = lambda: run_monarch(remote_terms, top, job_warn["monarch"])
    if "pubcasefinder" in sources:
        jobs["pubcasefinder"] = lambda: run_pubcasefinder(remote_terms, top, job_warn["pubcasefinder"])
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = {name: pool.submit(fn) for name, fn in jobs.items()}
        if "local" in sources and idx is not None:
            try:
                got = run_local(idx, present, excluded, top)
                per["local"] = got.result
                outcome.add(got)
            except ValueError as err:  # e.g. only onset/inheritance terms given
                job_warn["local"].append(f"local ranking skipped: {err}")
        for name in ("monarch", "pubcasefinder"):
            if name in futures:
                got = futures[name].result()
                if got.result is not None:
                    per[name] = got.result
                outcome.add(got)
    for s in SOURCES:
        warnings.extend(job_warn[s])
    if excluded:
        for s in ("monarch", "pubcasefinder"):
            if s in per:
                notes.append(f"{s} does not take excluded phenotypes; its list ignores {', '.join(excluded)}")
    if not per:
        warnings.append("no source returned a ranking")

    items = _ranked_items(per)
    keys, link_outcomes = link_ids(items, warnings) if len(per) >= 2 else (resolve_keys(items, None, None), [])
    for o in link_outcomes:
        outcome.add(o)
    consensus = build_consensus(items, keys) if len(per) >= 2 else []
    genes_agree = gene_consensus(per) if len(per) >= 2 else []

    for s in ("local", "monarch"):  # labels: local HPO release first, else as Monarch echoed them
        for k, v in ((per.get(s) or {}).pop("query_labels", None) or {}).items():
            if v and not labels.get(k):
                labels[k] = v
    result = {
        "query": {"present": [{"id": t, "label": labels.get(t)} for t in present],
                  "excluded": [{"id": t, "label": labels.get(t)} for t in excluded],
                  "from_case": used_case},
        "sources": [s for s in sources if s in per],
        "top": top,
        "per_source": per,
        "consensus": consensus,
        "consensus_rule": f"diseases in the top {top} of >= 2 sources; ids joined only by exact mappings "
                          "(Monarch exactMatch / MONDO xref, obsolete MONDO followed via OLS); ordered by number of "
                          "sources, then by the worst of their best ranks. No combined score.",
        "gene_consensus": genes_agree,
        "notes": notes,
    }
    outcome.result = result
    outcome.warnings = warnings + outcome.warnings
    outcome.query = {"present": present, "excluded": excluded, "sources": sources, "top": top, "from_case": used_case}
    outcome.text = render(result)
    return outcome


def render(r: Dict[str, Any]) -> str:
    per = r["per_source"]
    top = r["top"]
    lines = [f"phenotype rank: {len(r['query']['present'])} present, {len(r['query']['excluded'])} excluded; "
             f"sources {', '.join(r['sources']) or 'none'} (top {top} each)"]
    lines.append("present: " + ", ".join(f"{q['id']} {q['label'] or ''}".strip() for q in r["query"]["present"]))
    if r["query"]["excluded"]:
        lines.append("excluded: " + ", ".join(f"{q['id']} {q['label'] or ''}".strip() for q in r["query"]["excluded"]))
    if len(r["sources"]) >= 2:
        lines.append("")
        lines.append(f"consensus (in the top {top} of >= 2 sources; ordered by #sources, then worst rank; no combined score):")
        if not r["consensus"]:
            lines.append("  none: the sources do not agree on any disease in their top lists")
        for c in r["consensus"][:20]:
            where = "  ".join(f"{h['list']} #{h['rank']} ({h['id']})" for h in c["hits"])
            lines.append(f"{c['order']:>2}. {c['name']} [{c['key']}]  {where}")
            for how in c["linked_by"]:
                if how != "same id":
                    lines.append(f"      link: {how}")
    if per.get("local"):
        b = per["local"]
        lines.append("")
        lines.append(f"local (Resnik BMA over HPO {b['hpo_version']}; excluded terms penalised):")
        for d in b["diseases"]:
            pen = f" excl:{','.join(d['excluded_hits'])}" if d["excluded_hits"] else ""
            lines.append(f"{d['rank']:>2}. {d['id']:<13} {d['score']:.3f} ({(d['relative'] or 0):.0%}) {d['name']}  "
                         f"[{d['exact_matches']}/{d['of']} exact] {','.join(d['genes'][:5])}{pen}")
    if per.get("monarch"):
        b = per["monarch"]
        lines.append("")
        lines.append("monarch (semsim, ancestor IC; excluded terms not supported):")
        for d in b["diseases"]:
            lines.append(f"{d['rank']:>2}. {d['id']:<14} {d['score']:.2f} {d['name']}  [{d['exact_matches']}/{d['of']} exact]")
    if per.get("pubcasefinder"):
        b = per["pubcasefinder"]
        lines.append("")
        lines.append("pubcasefinder (excluded terms not supported):")
        for lst in ("diseases_omim", "diseases_orphanet"):
            if b[lst]:
                lines.append(f"  {lst.replace('diseases_', '')}:")
            for d in b[lst]:
                lines.append(f"  {d['rank']:>2}. {d['id']:<13} {d['score']:.3f} {d['name']}  "
                             f"[{len(d['matched'])} matched] {','.join(d['genes'][:5])}")
    gl = []
    for s in r["sources"]:
        g = per[s].get("genes") or []
        if g:
            gl.append(f"  {s}: " + ", ".join(x["symbol"] or "?" for x in g[:top]))
    if gl:
        lines.append("")
        lines.append("genes:")
        lines.extend(gl)
    if r["gene_consensus"]:
        lines.append("  in >= 2 sources: " + ", ".join(
            f"{g['symbol']} ({' '.join(f'{s}#{k}' for s, k in g['rank'].items())})" for g in r["gene_consensus"][:15]))
    for n in r["notes"]:
        lines.append(f"note: {n}")
    return "\n".join(lines)


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("phenotype", help="phenotype-driven differential: rank diseases and genes across sources")
    ps = p.add_subparsers(dest="action", metavar="<action>")
    q = ps.add_parser("rank", help="rank diseases/genes for HPO terms with local, Monarch and PubCaseFinder, kept separate")
    q.add_argument("hpo", nargs="*", help="present phenotypes HP:... (same as --present)")
    q.add_argument("--present", nargs="+", action="extend", default=[], help="phenotypes present HP:...")
    q.add_argument("--exclude", nargs="+", action="extend", default=[], help="phenotypes known to be absent")
    q.add_argument("--from-case", action="store_true", help="add the present/excluded phenotypes of the active case")
    q.add_argument("--sources", default="local,monarch,pubcasefinder", help="comma list: local,monarch,pubcasefinder")
    q.add_argument("--top", type=int, default=15, help=f"results per source (max {MAX_TOP})")
    q.set_defaults(func=_rank)
