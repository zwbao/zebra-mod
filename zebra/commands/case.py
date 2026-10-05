"""zebra case: the local case workspace (profile, evidence ledger)."""

from __future__ import annotations

import argparse
from typing import Any, Dict, List

from zebra import case as case_mod
from zebra.core import Outcome, UsageError


def _dir(args: argparse.Namespace) -> str:
    target = getattr(args, "dir", None) or args.case
    if not target:
        raise UsageError("which case? pass a directory or --case (or set ZEBRA_CASE)")
    return target


def _wrap(fn):
    def run(args: argparse.Namespace) -> Outcome:
        try:
            return fn(args)
        except case_mod.CaseError as err:
            raise UsageError(str(err)) from None

    return run


def _init(args: argparse.Namespace) -> Outcome:
    data = case_mod.init(_dir(args), title=args.title or "", role=args.role, language=args.lang)
    path = case_mod.case_file(_dir(args)).parent
    return Outcome({"path": str(path), "case": data}, text=f"case created at {path}\n  case.json, evidence/, records/, reports/")


def _show(args: argparse.Namespace) -> Outcome:
    data = case_mod.load(_dir(args))
    return Outcome(data)


def _summary(args: argparse.Namespace) -> Outcome:
    s = case_mod.summary(_dir(args))
    lines = [f"{s['title']}  ({s['role']}, {s['path']})"]
    present = [p for p in s["phenotypes"] if p["status"] == "present"]
    excluded = [p for p in s["phenotypes"] if p["status"] == "excluded"]
    lines.append(f"phenotypes: {len(present)} present, {len(excluded)} excluded")
    for p in present:
        lines.append(f"  + {p['id']} {p['label']}")
    for p in excluded:
        lines.append(f"  - {p['id']} {p['label']}")
    lines.append(f"variants: {len(s['variants'])}")
    for v in s["variants"]:
        lines.append(f"  {v['id']} {v['gene'] or ''} {v['label']} {v['zygosity'] or ''} [{v['classification'] or 'unclassified'}]")
    lines.append(f"hypotheses: {len(s['hypotheses'])}")
    for h in s["hypotheses"]:
        lines.append(f"  {h['id']} [{h['status']}] {h['disease']} (+{h['support']} / -{h['against']})")
    if s["therapy_leads"]:
        lines.append("therapy leads:")
        for t in s["therapy_leads"]:
            lines.append(f"  {t['id']} [{t['kind']}] {t['name']} {t['status'] or ''}")
    if s["questions"]:
        lines.append("open questions:")
        for q in s["questions"]:
            lines.append(f"  ? {q}")
    lines.append(f"evidence ledger: {s['evidence_count']} rows; protected identifiers: {s['identifiers']}")
    return Outcome(s, text="\n".join(lines))


def _add_hpo(args: argparse.Namespace) -> Outcome:
    added: List[Dict[str, Any]] = []
    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []
    for hpo_id in args.hpo:
        label = args.label if len(args.hpo) == 1 else None
        if not label:
            try:
                from zebra.sources import hpo as hpo_src

                term = hpo_src.term(hpo_id)
                label = term.result.get("name")
                sources.extend(term.sources)
                if term.result.get("obsolete"):
                    warnings.append(f"{hpo_id} is obsolete in HPO; replaced_by={term.result.get('replaced_by')}")
            except Exception as err:  # noqa: BLE001 - lookup is a convenience
                raise UsageError(f"could not verify {hpo_id} against HPO ({err}); pass --label after checking it yourself") from None
            if not label:
                raise UsageError(f"{hpo_id} was not found in HPO")
        added.append(case_mod.add_phenotype(_dir(args), hpo_id, label, status=args.status, onset=args.onset,
                                            source=args.source, note=args.note))
    text = "\n".join(f"{'+' if a['status'] == 'present' else '-'} {a['id']} {a['label']}" for a in added)
    return Outcome(added, sources=sources, warnings=warnings, text=text, query={"hpo": args.hpo})


def _add_variant(args: argparse.Namespace) -> Outcome:
    entry = case_mod.add_variant(
        _dir(args), gene=args.gene, hgvs_c=args.hgvs_c, hgvs_g=args.hgvs_g, hgvs_p=args.hgvs_p, vcf=args.vcf,
        assembly=args.assembly, zygosity=args.zygosity, inheritance=args.inheritance,
        classification_lab=args.lab_class, source=args.source, description=args.description,
    )
    return Outcome(entry, text=f"added {entry['id']}: {entry}")


def _add_hypothesis(args: argparse.Namespace) -> Outcome:
    ids = {}
    for pair in args.ids or []:
        if ":" not in pair:
            raise UsageError(f"--id takes PREFIX:VALUE (ORPHA:33069, OMIM:607208, MONDO:0100135), got {pair!r}")
        prefix, value = pair.split(":", 1)
        ids[prefix.upper()] = value
    entry = case_mod.add_hypothesis(_dir(args), args.disease, status=args.status, ids=ids, support=args.support,
                                    against=args.against, note=args.note)
    return Outcome(entry, text=f"{entry['id']} [{entry['status']}] {entry['disease']}")


def _add_lead(args: argparse.Namespace) -> Outcome:
    entry = case_mod.add_lead(_dir(args), args.name, args.kind, status=args.status, evidence=args.evidence, note=args.note)
    return Outcome(entry, text=f"{entry['id']} [{entry['kind']}] {entry['name']}")


def _add_question(args: argparse.Namespace) -> Outcome:
    qs = case_mod.add_question(_dir(args), args.text)
    return Outcome(qs, text="\n".join(f"? {q}" for q in qs))


def _identifiers(args: argparse.Namespace) -> Outcome:
    current = case_mod.load(_dir(args))["privacy"]["identifiers"]
    if args.add or args.clear:
        merged = [] if args.clear else list(current)
        merged += args.add or []
        current = case_mod.set_identifiers(_dir(args), merged)
    return Outcome({"count": len(current)}, text=f"{len(current)} protected identifiers (values are not printed)")


def _remove(args: argparse.Namespace) -> Outcome:
    ok = case_mod.remove(_dir(args), args.kind, args.item)
    return Outcome({"removed": ok}, text="removed" if ok else "nothing matched")


def _ledger(args: argparse.Namespace) -> Outcome:
    rows = case_mod.read_ledger(_dir(args))
    if args.eid:
        rows = [r for r in rows if r.get("eid") in set(args.eid)]
    text = "\n".join(f"{r['eid']}\t{r.get('db')}\t{r.get('record') or ''}\t{r.get('url') or ''}\t{r.get('retrieved_at') or ''}" for r in rows)
    return Outcome(rows, text=text or "(empty ledger)")


def _apply(args: argparse.Namespace) -> Outcome:
    """Several updates at once (what the mod's case_update tool sends)."""
    import json as _json

    try:
        ops = _json.loads(args.ops)
    except ValueError as err:
        raise UsageError(f"--ops is not JSON: {err}") from None
    if not isinstance(ops, dict):
        raise UsageError("--ops must be a JSON object")
    target = _dir(args)
    done: Dict[str, Any] = {}
    sources: List[Dict[str, Any]] = []
    warnings: List[str] = []
    errors: List[str] = []

    phenos = ops.get("phenotypes") or []
    if phenos:
        from zebra.sources import hpo as hpo_src

        added = []
        for p in phenos:
            hid = str(p.get("id", "")).strip()
            try:
                term = hpo_src.term(hid)
            except Exception as err:  # noqa: BLE001
                errors.append(f"{hid}: not verified against HPO ({err}); not recorded")
                continue
            if not term.result.get("name"):
                errors.append(f"{hid}: not found in HPO; not recorded")
                continue
            if term.result.get("obsolete"):
                warnings.append(f"{hid} is obsolete (replaced by {term.result.get('replaced_by')}); recorded as given")
            sources.extend(term.sources)
            try:
                added.append(case_mod.add_phenotype(target, hid, term.result["name"], status=p.get("status", "present"),
                                                    onset=p.get("onset"), source=p.get("source"), note=p.get("note")))
            except case_mod.CaseError as err:
                errors.append(f"{hid}: {err}")
        done["phenotypes"] = added
    for v in ops.get("variants") or []:
        try:
            done.setdefault("variants", []).append(case_mod.add_variant(target, **{k: v.get(k) for k in (
                "gene", "hgvs_c", "hgvs_g", "hgvs_p", "vcf", "assembly", "zygosity", "inheritance",
                "classification_lab", "source", "description")}))
        except case_mod.CaseError as err:
            errors.append(f"variant {v}: {err}")
    for h in ops.get("hypotheses") or []:
        ids = {}
        for pair in h.get("ids") or []:
            if ":" in str(pair):
                prefix, value = str(pair).split(":", 1)
                ids[prefix.upper()] = value
        try:
            done.setdefault("hypotheses", []).append(case_mod.add_hypothesis(
                target, h["disease"], status=h.get("status", "considered"), ids=ids,
                support=h.get("support"), against=h.get("against"), note=h.get("note")))
        except (case_mod.CaseError, KeyError) as err:
            errors.append(f"hypothesis {h}: {err}")
    for t in ops.get("leads") or []:
        try:
            done.setdefault("leads", []).append(case_mod.add_lead(target, t["name"], t["kind"], status=t.get("status"),
                                                                  evidence=t.get("evidence"), note=t.get("note")))
        except (case_mod.CaseError, KeyError) as err:
            errors.append(f"lead {t}: {err}")
    for a in ops.get("acmg") or []:
        try:
            from zebra import acmg as acmg_mod

            res = acmg_mod.classify(list(a.get("codes") or []))
            done.setdefault("acmg", []).append(case_mod.set_acmg(
                target, str(a["variant_id"]), res["classification"], res["points"],
                [c["label"] for c in res["codes"]], note=a.get("note")))
        except (case_mod.CaseError, KeyError, ValueError) as err:
            errors.append(f"acmg {a}: {err}")
    for q in ops.get("questions") or []:
        done["questions"] = case_mod.add_question(target, str(q))
    for r in ops.get("remove") or []:
        try:
            done.setdefault("removed", []).append({**r, "ok": case_mod.remove(target, r["kind"], r["id"])})
        except (case_mod.CaseError, KeyError) as err:
            errors.append(f"remove {r}: {err}")
    if errors:
        warnings.extend(errors)
    summary = case_mod.summary(target)
    result = {"applied": done, "errors": errors, "case": {k: summary[k] for k in ("title", "evidence_count")},
              "counts": {"phenotypes": len(summary["phenotypes"]), "variants": len(summary["variants"]),
                         "hypotheses": len(summary["hypotheses"]), "leads": len(summary["therapy_leads"])}}
    return Outcome(result, sources=sources, warnings=warnings, query={"ops": list(ops)})


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("case", help="local case workspace: profile, phenotypes, variants, hypotheses, evidence ledger")
    cs = p.add_subparsers(dest="action", metavar="<action>")

    def add(name: str, fn, help_: str, positional_dir: bool = True) -> argparse.ArgumentParser:
        q = cs.add_parser(name, help=help_)
        if positional_dir:
            q.add_argument("dir", nargs="?", help="case directory (default --case / $ZEBRA_CASE)")
        q.set_defaults(func=_wrap(fn), no_ledger=True)
        return q

    q = add("init", _init, "create a case folder")
    q.add_argument("--title")
    q.add_argument("--role", choices=case_mod.ROLES, default="family")
    q.add_argument("--lang", default="zh", help="report language: zh or en")
    add("show", _show, "print case.json")
    add("summary", _summary, "one-screen summary (what the case board shows)")

    q = add("add-hpo", _add_hpo, "add phenotypes by HPO id (label verified against HPO)", positional_dir=False)
    q.add_argument("hpo", nargs="+", help="HP:0001250 ...")
    q.add_argument("--label", help="only with one id, and only when HPO cannot be reached")
    q.add_argument("--status", choices=case_mod.PHENO_STATUS, default="present")
    q.add_argument("--onset", help="HPO onset term or free text (e.g. HP:0003593 or '6 months')")
    q.add_argument("--source", help="where it is documented (records/report.pdf p.2)")
    q.add_argument("--note")
    q.set_defaults(no_ledger=False)

    q = add("add-variant", _add_variant, "add a variant", positional_dir=False)
    q.add_argument("--gene")
    q.add_argument("--hgvs-c", dest="hgvs_c", help="NM_000492.4:c.1521_1523del")
    q.add_argument("--hgvs-g", dest="hgvs_g")
    q.add_argument("--hgvs-p", dest="hgvs_p")
    q.add_argument("--vcf", help="chrom-pos-ref-alt, e.g. 7-117559590-ATCT-A")
    q.add_argument("--assembly", choices=("GRCh38", "GRCh37"), default="GRCh38")
    q.add_argument("--zygosity", choices=("het", "hom", "hemi", "mosaic", "unknown"))
    q.add_argument("--inheritance", choices=("de_novo", "maternal", "paternal", "biparental", "unknown"))
    q.add_argument("--lab-class", dest="lab_class", help="classification on the lab report (P/LP/VUS/LB/B)")
    q.add_argument("--source")
    q.add_argument("--description", help="free text when there is no HGVS (e.g. 'exon 45-50 deletion')")

    q = add("add-hypothesis", _add_hypothesis, "add or update a diagnostic hypothesis", positional_dir=False)
    q.add_argument("disease")
    q.add_argument("--status", choices=case_mod.HYP_STATUS, default="considered")
    q.add_argument("--id", dest="ids", action="append", help="ORPHA:33069 / OMIM:607208 / MONDO:0100135 (repeatable)")
    q.add_argument("--support", nargs="*", help="evidence ids from the ledger (E3 E7)")
    q.add_argument("--against", nargs="*")
    q.add_argument("--note")

    q = add("add-lead", _add_lead, "add a therapy lead", positional_dir=False)
    q.add_argument("name")
    q.add_argument("--kind", choices=case_mod.LEAD_KINDS, required=True)
    q.add_argument("--status")
    q.add_argument("--evidence", nargs="*")
    q.add_argument("--note")

    q = add("add-question", _add_question, "add an open question for the care team", positional_dir=False)
    q.add_argument("text")

    q = add("identifiers", _identifiers, "strings the mod must never send out (names, birth dates, record numbers)", positional_dir=False)
    q.add_argument("--add", nargs="*")
    q.add_argument("--clear", action="store_true")

    q = add("remove", _remove, "remove a phenotype / variant / hypothesis / lead", positional_dir=False)
    q.add_argument("kind", choices=("phenotype", "variant", "hypothesis", "lead"))
    q.add_argument("item")

    q = add("apply", _apply, "apply several updates from JSON (phenotypes, variants, hypotheses, leads, questions, remove)")
    q.add_argument("--ops", required=True, help="JSON object")

    q = add("ledger", _ledger, "list evidence ledger rows", positional_dir=False)
    q.add_argument("eid", nargs="*", help="only these ids")
