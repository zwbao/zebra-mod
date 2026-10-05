"""zebra hpo: phenotype terms (search, verify) and offline phenotype-driven ranking."""

from __future__ import annotations

import argparse
import os
import shutil
import urllib.request
from typing import Any, Dict, List

from zebra import hpo_local
from zebra.core import Outcome, UsageError, attempt
from zebra.http import USER_AGENT, _opener
from zebra.sources import record as source_record


def _fetch(args: argparse.Namespace) -> Outcome:
    target = hpo_local.data_dir()
    target.mkdir(parents=True, exist_ok=True)
    done: List[Dict[str, Any]] = []
    for name in hpo_local.FILES:
        path = target / name
        if path.exists() and path.stat().st_size > 0 and not args.force:
            done.append({"file": name, "bytes": path.stat().st_size, "status": "present"})
            continue
        url = hpo_local.RELEASE_URL.format(name=name)
        tmp = path.with_suffix(path.suffix + ".part")
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with _opener().open(req, timeout=900) as resp, open(tmp, "wb") as fh:
            shutil.copyfileobj(resp, fh, length=1 << 20)
        os.replace(tmp, path)
        done.append({"file": name, "bytes": path.stat().st_size, "status": "downloaded", "url": url})
    warnings: List[str] = []
    zh_path = target / hpo_local.ZH_FILE
    if not zh_path.exists() or args.force:
        try:
            tmp = zh_path.with_suffix(".part")
            req = urllib.request.Request(hpo_local.ZH_URL, headers={"User-Agent": USER_AGENT})
            with _opener().open(req, timeout=600) as resp, open(tmp, "wb") as fh:
                shutil.copyfileobj(resp, fh, length=1 << 20)
            os.replace(tmp, zh_path)
            done.append({"file": hpo_local.ZH_FILE, "bytes": zh_path.stat().st_size, "status": "downloaded", "url": hpo_local.ZH_URL})
        except Exception as err:  # noqa: BLE001 - optional file
            warnings.append(f"Chinese HPO labels not downloaded ({err}); Chinese phenotype search stays off")
    idx = hpo_local.load(rebuild=True)
    result = {"dir": str(target), "files": done, "chinese_labels": len(idx.zh), "hpo_version": idx.version, "terms": len(idx.names),
              "diseases": len(idx.disease_list)}
    text = f"HPO {idx.version}: {len(idx.names)} terms, {len(idx.disease_list)} annotated diseases in {target}"
    return Outcome(result, text=text, warnings=warnings,
                   sources=[source_record("HPO release", idx.version, url=hpo_local.RELEASE_URL.format(name="hp.json"))])


def _local_or_none():
    try:
        return hpo_local.load()
    except hpo_local.HpoDataMissing:
        return None


def _search(args: argparse.Namespace) -> Outcome:
    text = " ".join(args.text)
    warnings: List[str] = []
    idx = _local_or_none()
    if not args.online and idx is not None:
        hits = hpo_local.search(idx, text, limit=args.limit)
        return Outcome({"query": text, "hits": hits, "backend": f"local HPO {idx.version}"},
                       sources=[source_record("HPO", idx.version, url="https://hpo.jax.org", note="local release files")],
                       text="\n".join(f"{h['id']}\t{h['label']}" + (f"\t{h['label_zh']}" if h.get('label_zh') else "")
                                       + (f"\t(matched: {h['matched']})" if h['matched'] not in (h['label'], h.get('label_zh')) else "") for h in hits)
                       or "no match",
                       query={"text": text})
    from zebra.sources import hpo as hpo_src

    out = hpo_src.search(text, limit=args.limit)
    if idx is None:
        warnings.append("local HPO files not present: searched the web (EBI OLS4 over HPO) and re-ranked the hits "
                        "exact-label-first. Run `zebra hpo fetch` once (~80 MB) for the full release, which also "
                        "searches Chinese labels and curated lay phrases")
    out.warnings.extend(warnings)
    out.text = "\n".join(
        f"{h['id']}\t{h['label']}"
        + (f"\t(matched: {h['matched']} [{h['matched_on']}])" if h.get("matched") and h["matched"] != h["label"] else "")
        for h in out.result["hits"]) or "no match"
    out.query = {"text": text}
    return out


def _term(args: argparse.Namespace) -> Outcome:
    idx = _local_or_none()
    rows = []
    sources = []
    warnings: List[str] = []
    for tid in args.ids:
        if idx is not None and not args.online:
            pid, note = idx.primary(tid)
            row = {"id": tid, "name": idx.names.get(pid) if pid else None, "name_zh": idx.zh.get(pid or ""), "primary": pid, "note": note,
                   "obsolete": tid in idx.obsolete, "replaced_by": idx.obsolete.get(tid),
                   "parents": [{"id": p, "name": idx.names.get(p)} for p in idx.parents.get(pid or "", [])],
                   "synonyms": idx.synonyms.get(pid or "", [])[:8]}
            rows.append(row)
            sources.append(source_record("HPO", tid, url=f"https://hpo.jax.org/browse/term/{tid}", note=f"local release {idx.version}"))
        else:
            from zebra.sources import hpo as hpo_src

            got = attempt(f"HPO {tid}", lambda: hpo_src.term(tid), warnings)
            if got:
                rows.append(got.result)
                sources.extend(got.sources)
    text = "\n".join(f"{r['id']}\t{r.get('name')}" + (f"\t[{r['note']}]" if r.get("note") else "") for r in rows)
    return Outcome(rows, sources=sources, warnings=warnings, text=text, query={"ids": args.ids})


def _rank(args: argparse.Namespace) -> Outcome:
    present = list(args.hpo)
    excluded = list(args.exclude or [])
    if args.from_case or (not present and args.case):
        from zebra import case as case_mod

        data = case_mod.load(args.case)
        present += [p["id"] for p in data["phenotypes"] if p.get("status", "present") == "present"]
        excluded += [p["id"] for p in data["phenotypes"] if p.get("status") == "excluded"]
    if not present:
        raise UsageError("give HPO ids, or --from-case with --case")
    try:
        idx = hpo_local.load()
    except hpo_local.HpoDataMissing as err:
        raise UsageError(str(err)) from None
    try:
        res = hpo_local.rank(idx, present, excluded, top=args.top, db=tuple(args.db.split(",")))
    except ValueError as err:
        raise UsageError(str(err)) from None
    lines = [f"{res['method']}; HPO {res['hpo_version']}; {res['diseases_scored']} diseases scored",
             "query: " + ", ".join(f"{q['id']} {q['label']}" for q in res["query"])]
    if res["excluded"]:
        lines.append("excluded: " + ", ".join(f"{q['id']} {q['label']}" for q in res["excluded"]))
    for i, d in enumerate(res["diseases"], 1):
        exact = sum(1 for m in d["matches"] if m["exact"])
        genes = ",".join(d["genes"][:6]) + ("…" if len(d["genes"]) > 6 else "")
        pen = f" penalty {d['penalty']} for {','.join(d['excluded_hits'])}" if d["excluded_hits"] else ""
        lines.append(f"{i:>2}. {d['disease']:<14} {d['score']:.3f} ({d['relative']:.0%})  {d['name']}  [{exact}/{len(d['matches'])} exact]  genes: {genes or '-'}{pen}")
    lines.append("top genes: " + ", ".join(f"{g['gene']}({g['score']:.2f})" for g in res["genes"][:15]))
    for n in res["notes"]:
        lines.append(f"note: {n}")
    sources = [source_record("HPO annotations (phenotype.hpoa)", res["hpo_version"],
                             url="https://hpo.jax.org/data/annotations", note="local Resnik ranking")]
    return Outcome(res, sources=sources, text="\n".join(lines), query={"present": present, "excluded": excluded})


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("hpo", help="HPO terms: search, verify ids, offline disease/gene ranking")
    hs = p.add_subparsers(dest="action", metavar="<action>")

    q = hs.add_parser("fetch", help="download HPO release files (~80 MB) for offline use")
    q.add_argument("--force", action="store_true")
    q.set_defaults(func=_fetch)

    q = hs.add_parser("search", help="find HPO terms for a phrase (English labels and synonyms)")
    q.add_argument("text", nargs="+")
    q.add_argument("--limit", type=int, default=10)
    q.add_argument("--online", action="store_true", help="use the HPO web API even when local files exist")
    q.set_defaults(func=_search)

    q = hs.add_parser("term", help="verify HPO ids: name, obsolete/replaced, parents")
    q.add_argument("ids", nargs="+")
    q.add_argument("--online", action="store_true")
    q.set_defaults(func=_term)

    q = hs.add_parser("rank", help="rank diseases and genes for a phenotype profile (offline, Resnik BMA)")
    q.add_argument("hpo", nargs="*", help="present phenotypes HP:...")
    q.add_argument("--exclude", nargs="*", help="phenotypes known to be absent")
    q.add_argument("--from-case", action="store_true", help="use the phenotypes recorded in --case")
    q.add_argument("--top", type=int, default=20)
    q.add_argument("--db", default="OMIM,ORPHA", help="disease namespaces to rank (OMIM,ORPHA,DECIPHER)")
    q.set_defaults(func=_rank)
