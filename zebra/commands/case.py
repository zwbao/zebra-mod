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
    fields = {"gene": args.gene, "hgvs_c": args.hgvs_c, "hgvs_g": args.hgvs_g, "hgvs_p": args.hgvs_p}
    entry = case_mod.add_variant(
        _dir(args), gene=args.gene, hgvs_c=args.hgvs_c, hgvs_g=args.hgvs_g, hgvs_p=args.hgvs_p, vcf=args.vcf,
        assembly=args.assembly, zygosity=args.zygosity, inheritance=args.inheritance,
        classification_lab=args.lab_class, source=args.source, description=args.description,
        kind=args.kind, region=args.region, iscn=args.iscn, cnv_type=args.cnv_type, copy_number=args.copy_number,
        exons=args.exons, genes=args.genes, motif=args.motif, repeat_count=args.repeats, method=args.method,
    )
    return Outcome(entry, warnings=case_mod.variant_notes(fields), text=f"added {entry['id']}: {entry}")


def _add_hypothesis(args: argparse.Namespace) -> Outcome:
    target = _dir(args)
    known = case_mod.ledger_ids(target)
    ids = {}
    for pair in args.ids or []:
        prefix, value = case_mod.check_id(pair).split(":", 1)
        ids[prefix] = value
    support = case_mod.check_evidence_ids(args.support, "--support", known)
    against = case_mod.check_evidence_ids(args.against, "--against", known)
    entry = case_mod.add_hypothesis(target, args.disease, status=args.status, ids=ids, support=support,
                                    against=against, note=args.note)
    return Outcome(entry, text=f"{entry['id']} [{entry['status']}] {entry['disease']}")


def _add_lead(args: argparse.Namespace) -> Outcome:
    target = _dir(args)
    evidence = case_mod.check_evidence_ids(args.evidence, "--evidence", case_mod.ledger_ids(target))
    ids = {}
    for pair in args.ids or []:
        prefix, value = case_mod.check_id(pair).split(":", 1)
        ids[prefix] = value
    entry = case_mod.add_lead(target, args.name, args.kind, status=args.status, evidence=evidence,
                              note=args.note, ids=ids)
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
    text = "\n".join(f"{r.get('eid') or '(no id)'}\t{r.get('db')}\t{r.get('record') or ''}\t"
                     f"{r.get('url') or ''}\t{r.get('retrieved_at') or ''}" for r in rows)
    return Outcome(rows, text=text or "(empty ledger)")


OPS_KEYS = ("profile", "phenotypes", "variants", "hypotheses", "leads", "acmg", "questions", "remove")
# the fields each list item may carry; an unknown one is a typo that would
# otherwise be dropped in silence while the envelope reported success
ITEM_FIELDS = {
    "phenotypes": ("id", "status", "onset", "source", "note"),
    "hypotheses": ("disease", "status", "ids", "support", "against", "note"),
    "leads": ("name", "kind", "status", "evidence", "ids", "note"),
    "acmg": ("variant_id", "codes", "note"),
    "remove": ("kind", "id"),
}
# keys people and models reach for, and what they meant
OPS_ALIASES = {"phenotype": "phenotypes", "variant": "variants", "hypothesis": "hypotheses",
               "hypothese": "hypotheses", "lead": "leads", "therapy_leads": "leads", "therapy": "leads",
               "question": "questions", "removals": "remove", "delete": "remove", "acmg_codes": "acmg"}
_LIST_OPS = ("phenotypes", "variants", "hypotheses", "leads", "acmg", "remove")


def _ids_map(raw: Any, field: str) -> Dict[str, str]:
    """`["ORPHA:33069", …]` → {"ORPHA": "33069", …}, each shape-checked.

    The dict form `case show` prints is accepted too, so what is read out of a
    case can be fed straight back in.
    """
    out: Dict[str, str] = {}
    if isinstance(raw, dict):
        raw = [f"{k}:{v}" for k, v in raw.items()]
    for pair in case_mod.str_list(raw, field):
        prefix, value = case_mod.check_id(pair).split(":", 1)
        out[prefix] = value
    return out


def _check_ops(ops: Any, known_evidence: set) -> Dict[str, Any]:
    """Validate every op before anything is written. Raises UsageError; nothing is partially applied."""
    if not isinstance(ops, dict):
        raise UsageError("--ops must be a JSON object")
    unknown = [k for k in ops if k not in OPS_KEYS]
    if unknown:
        hints = [f"{k!r} (did you mean {OPS_ALIASES[k]!r}?)" if k in OPS_ALIASES else repr(k) for k in unknown]
        raise UsageError("unknown key(s) in --ops: " + ", ".join(hints)
                         + "; accepted keys are " + ", ".join(OPS_KEYS))
    clean: Dict[str, Any] = {}
    try:
        if ops.get("profile") is not None:
            if not isinstance(ops["profile"], dict):
                raise case_mod.CaseError(f"profile must be an object, not {type(ops['profile']).__name__}")
            bad = [k for k in ops["profile"] if k not in case_mod.PROFILE_FIELDS]
            if bad:
                raise case_mod.CaseError("unknown profile field(s) " + ", ".join(sorted(bad))
                                         + "; accepted: " + ", ".join(case_mod.PROFILE_FIELDS))
            clean["profile"] = dict(ops["profile"])
            for key_, value_ in clean["profile"].items():
                if key_ in ("title", "role", "language") and value_ is not None:
                    if not isinstance(value_, str) or not value_.strip():
                        raise case_mod.CaseError(f"profile.{key_} must be a non-empty string, got {value_!r}")
                    if key_ == "role" and value_ not in case_mod.ROLES:
                        raise case_mod.CaseError(f"role must be one of {', '.join(case_mod.ROLES)}")
                if key_ in ("sex", "age", "ancestry") and value_ is not None and not isinstance(value_, str):
                    raise case_mod.CaseError(f"profile.{key_} must be a string, got {type(value_).__name__}")
                if key_ == "consanguinity" and value_ is not None and not isinstance(value_, (bool, str)):
                    raise case_mod.CaseError("profile.consanguinity must be true, false or a short description")
        for key in _LIST_OPS:
            if ops.get(key) is None:
                continue
            if isinstance(ops[key], dict):
                raise case_mod.CaseError(f"{key} must be a list of objects; wrap the single object in a list")
            if not isinstance(ops[key], list):
                raise case_mod.CaseError(f"{key} must be a list of objects, not {type(ops[key]).__name__}")
            for item in ops[key]:
                if not isinstance(item, dict):
                    raise case_mod.CaseError(f"{key}: every item must be an object, got {item!r}")
                allowed = ITEM_FIELDS.get(key) or case_mod.VARIANT_FIELDS
                unknown_fields = [f for f in item if f not in allowed]
                if unknown_fields:
                    raise case_mod.CaseError(f"{key}: unknown field(s) " + ", ".join(sorted(unknown_fields))
                                             + "; accepted: " + ", ".join(allowed))
            clean[key] = list(ops[key])
        if ops.get("questions") is not None:
            clean["questions"] = case_mod.str_list(ops["questions"], "questions")

        for p in clean.get("phenotypes", []):
            hid = str(p.get("id", "")).strip()
            if not case_mod.HPO_RE.match(hid):
                raise case_mod.CaseError(f"phenotypes: {hid!r} is not an HPO id (HP:0001250); find it with hpo_search")
            if p.get("status") not in (None, *case_mod.PHENO_STATUS):
                raise case_mod.CaseError(f"phenotypes: status must be one of {', '.join(case_mod.PHENO_STATUS)}")
            for field in ("onset", "source", "note"):
                if p.get(field) is not None and not isinstance(p[field], str):
                    raise case_mod.CaseError(f"phenotypes[{hid}].{field} must be a string, "
                                             f"not {type(p[field]).__name__}")
        for v in clean.get("variants", []):
            bad = [k for k in v if k not in case_mod.VARIANT_FIELDS]
            if bad:
                raise case_mod.CaseError("unknown variant field(s) " + ", ".join(sorted(bad))
                                         + "; accepted: " + ", ".join(case_mod.VARIANT_FIELDS))
            v.update(case_mod.check_variant(**{k: v.get(k) for k in case_mod.VARIANT_FIELDS}))
            v["_notes"] = case_mod.variant_notes(v)
        for h in clean.get("hypotheses", []):
            if not isinstance(h.get("disease"), str) or not h["disease"].strip():
                raise case_mod.CaseError(f"hypotheses: disease must be a non-empty string, got {h.get('disease')!r}")
            if h.get("status") not in (None, *case_mod.HYP_STATUS):
                raise case_mod.CaseError(f"hypotheses: status must be one of {', '.join(case_mod.HYP_STATUS)}")
            if h.get("note") is not None and not isinstance(h["note"], str):
                raise case_mod.CaseError(f"hypotheses[{h['disease']}].note must be a string, "
                                         f"not {type(h['note']).__name__}")
            h["_ids"] = _ids_map(h.get("ids"), f"hypotheses[{h['disease']}].ids")
            for field in ("support", "against"):
                h[f"_{field}"] = case_mod.check_evidence_ids(h.get(field), f"hypotheses[{h['disease']}].{field}",
                                                             known_evidence)
        for t in clean.get("leads", []):
            if not isinstance(t.get("name"), str) or not t["name"].strip():
                raise case_mod.CaseError(f"leads: name must be a non-empty string, got {t.get('name')!r}")
            for field in ("status", "note"):
                if t.get(field) is not None and not isinstance(t[field], str):
                    raise case_mod.CaseError(f"leads[{t['name']}].{field} must be a string, "
                                             f"not {type(t[field]).__name__}")
            if t.get("kind") not in case_mod.LEAD_KINDS:
                raise case_mod.CaseError(f"leads[{t['name']}]: kind must be one of {', '.join(case_mod.LEAD_KINDS)}")
            t["_evidence"] = case_mod.check_evidence_ids(t.get("evidence"), f"leads[{t['name']}].evidence",
                                                         known_evidence)
            t["_ids"] = _ids_map(t.get("ids"), f"leads[{t['name']}].ids")
        for a in clean.get("acmg", []):
            if not str(a.get("variant_id") or "").strip():
                raise case_mod.CaseError("acmg: variant_id is required (v1, v2 … as case_status lists them)")
            a["_codes"] = case_mod.str_list(a.get("codes"), f"acmg[{a.get('variant_id')}].codes")
            if not a["_codes"]:
                raise case_mod.CaseError(f"acmg[{a.get('variant_id')}]: codes must list the ACMG codes you justified")
            from zebra import acmg as acmg_mod

            for code in a["_codes"]:
                try:
                    acmg_mod.parse_code(code)
                except (ValueError, KeyError) as err:
                    raise case_mod.CaseError(f"acmg[{a.get('variant_id')}]: {err}") from None
        for r in clean.get("remove", []):
            if r.get("kind") not in case_mod.REMOVE_KINDS:
                raise case_mod.CaseError("remove: kind must be one of " + ", ".join(case_mod.REMOVE_KINDS))
            if not str(r.get("id") or "").strip():
                raise case_mod.CaseError("remove: id is required (v1, h2, t1, or an HPO id)")
    except case_mod.CaseError as err:
        raise UsageError(f"{err} — nothing was written") from None
    return clean


def _apply(args: argparse.Namespace) -> Outcome:
    """Several updates at once (what the mod's case_update tool sends).

    Every op is validated first, the HPO labels are verified before the lock is
    taken, and then all of them are applied to one loaded case and written once.
    A bad op therefore never leaves the case half-updated.
    """
    import json as _json

    try:
        ops_raw = _json.loads(args.ops)
    except ValueError as err:
        raise UsageError(f"--ops is not JSON: {err}") from None
    target = _dir(args)
    case_mod.load(target)  # fail here, before any work, if there is no readable case
    # the positional directory IS the case: without this the mod's own call
    # (case apply <dir> --ops …, no --case) would verify HPO terms against the
    # live source and then throw the provenance away
    if not getattr(args, "case", None):
        args.case = target
    ops = _check_ops(ops_raw, case_mod.ledger_ids(target))

    sources: List[Dict[str, Any]] = []
    warnings: List[str] = []
    errors: List[str] = []

    # HPO verification is a network call: do it before taking the lock
    labels: Dict[str, str] = {}
    if ops.get("phenotypes"):
        from zebra.sources import hpo as hpo_src

        for p in ops["phenotypes"]:
            hid = str(p["id"]).strip()
            if hid in labels:
                continue
            try:
                term = hpo_src.term(hid)
            except Exception as err:  # noqa: BLE001 - any lookup failure means "not verified"
                errors.append(f"{hid}: not verified against HPO ({err}); not recorded")
                continue
            if not term.result.get("name"):
                errors.append(f"{hid}: not found in HPO; not recorded")
                continue
            if term.result.get("obsolete"):
                warnings.append(f"{hid} is obsolete (replaced by {term.result.get('replaced_by')}); recorded as given")
            sources.extend(term.sources)
            labels[hid] = term.result["name"]

    done: Dict[str, Any] = {}
    with case_mod.editing(target) as data:
        if "profile" in ops:
            done["profile"] = case_mod.apply_profile(data, **ops["profile"])
        if ops.get("phenotypes") is not None:
            added = []
            for p in ops["phenotypes"]:
                hid = str(p["id"]).strip()
                if hid not in labels:
                    continue
                # status=None keeps the status the term already has: re-recording a
                # term must not turn a clinician's "excluded" into "present"
                added.append(case_mod.apply_phenotype(data, hid, labels[hid], status=p.get("status"),
                                                      onset=p.get("onset"), source=p.get("source"),
                                                      note=p.get("note")))
            done["phenotypes"] = added
        for v in ops.get("variants", []):
            warnings.extend(v.get("_notes") or [])
            done.setdefault("variants", []).append(
                case_mod.apply_variant(data, **{k: v.get(k) for k in case_mod.VARIANT_FIELDS}))
        for h in ops.get("hypotheses", []):
            done.setdefault("hypotheses", []).append(case_mod.apply_hypothesis(
                data, h["disease"], status=h.get("status"), ids=h["_ids"],
                support=h["_support"], against=h["_against"], note=h.get("note")))
        for t in ops.get("leads", []):
            done.setdefault("leads", []).append(case_mod.apply_lead(
                data, t["name"], t["kind"], status=t.get("status"), evidence=t["_evidence"],
                note=t.get("note"), ids=t["_ids"]))
        for a in ops.get("acmg", []):
            from zebra import acmg as acmg_mod

            try:
                res = acmg_mod.classify(list(a["_codes"]))
            except (ValueError, KeyError) as err:
                errors.append(f"acmg {a.get('variant_id')}: {err}")
                continue
            try:
                done.setdefault("acmg", []).append(case_mod.apply_acmg(
                    data, str(a["variant_id"]), res["classification"], res["points"],
                    [c["label"] for c in res["codes"]], note=a.get("note")))
            except case_mod.CaseError as err:
                errors.append(f"acmg {a.get('variant_id')}: {err}")
        for q in ops.get("questions", []):
            done["questions"] = case_mod.apply_question(data, q)
        for r in ops.get("remove", []):
            done.setdefault("removed", []).append({"kind": r["kind"], "id": r["id"],
                                                   "ok": case_mod.apply_remove(data, r["kind"], r["id"])})
    if errors:
        warnings.extend(errors)
    summary = case_mod.summary(target)
    result = {"applied": done, "errors": errors, "case": {k: summary[k] for k in ("title", "evidence_count")},
              "counts": {"phenotypes": len(summary["phenotypes"]), "variants": len(summary["variants"]),
                         "hypotheses": len(summary["hypotheses"]), "leads": len(summary["therapy_leads"])}}
    return Outcome(result, sources=sources, warnings=warnings, query={"ops": sorted(ops)})


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
    q.add_argument("--kind", choices=case_mod.VARIANT_KINDS, default="small",
                   help="small (SNV/indel, default), cnv (CMA/CNV-seq), exon_cnv (MLPA exon del/dup), "
                        "copy_number (SMN1/SMN2), repeat_expansion")
    q.add_argument("--region", help="cnv: 15:23000000-28500000 (with --assembly)")
    q.add_argument("--iscn", help="cnv: the ISCN string as the report writes it")
    q.add_argument("--cnv-type", dest="cnv_type", choices=case_mod.CNV_TYPES, help="gain/loss as reported")
    q.add_argument("--copy-number", dest="copy_number", type=int, help="copy number as reported (CNV, SMN1/SMN2)")
    q.add_argument("--exons", help="exon_cnv: 45-50, or 7 (numbering is the report's transcript)")
    q.add_argument("--genes", nargs="*", help="genes the finding spans, as the report names them")
    q.add_argument("--motif", help="repeat_expansion: CGG, CAG …")
    q.add_argument("--repeats", help="repeat_expansion: the count or range as reported")
    q.add_argument("--method", help="how it was measured (CMA, CNV-seq, MLPA, ddPCR, repeat-primed PCR, …)")

    q = add("add-hypothesis", _add_hypothesis, "add or update a diagnostic hypothesis", positional_dir=False)
    q.add_argument("disease")
    q.add_argument("--status", choices=case_mod.HYP_STATUS, default=None,
                   help="leading/considered/excluded/confirmed; omitted keeps the status it already has")
    q.add_argument("--id", dest="ids", action="append", help="ORPHA:33069 / OMIM:607208 / MONDO:0100135 (repeatable)")
    q.add_argument("--support", nargs="*", help="evidence ids from the ledger (E3 E7)")
    q.add_argument("--against", nargs="*")
    q.add_argument("--note")

    q = add("add-lead", _add_lead, "add a therapy lead", positional_dir=False)
    q.add_argument("name")
    q.add_argument("--kind", choices=case_mod.LEAD_KINDS, required=True)
    q.add_argument("--status")
    q.add_argument("--evidence", nargs="*", help="evidence ids from the ledger (E3 E7)")
    q.add_argument("--id", dest="ids", action="append", help="NCT04006210 / PMID:28919360 (repeatable)")
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
    # the HPO verifications this performs are evidence: they must get E-ids, or a
    # hypothesis can never cite the phenotype it rests on
    q.set_defaults(no_ledger=False)

    q = add("ledger", _ledger, "list evidence ledger rows", positional_dir=False)
    q.add_argument("eid", nargs="*", help="only these ids")
