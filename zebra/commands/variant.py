"""zebra variant: one variant annotated (VEP, gnomAD, ClinVar, predictors, LitVar) with ACMG inputs."""

from __future__ import annotations

import argparse
from typing import Any, Dict, List

from zebra.core import Outcome, UsageError


def _f(x: Any, nd: int = 3) -> str:
    if x is None:
        return "-"
    if isinstance(x, float):
        if x != 0 and abs(x) < 0.001:
            return f"{x:.2e}"
        return f"{x:.{nd}g}" if abs(x) < 1 else f"{x:.{nd}f}".rstrip("0").rstrip(".")
    return str(x)


def _pct(x: Any) -> str:
    return f"{x:.0%}" if isinstance(x, (int, float)) else "?"


def render(r: Dict[str, Any]) -> str:
    lines: List[str] = []
    tx = r.get("transcript") or {}
    hg = r.get("hgvs") or {}
    head = " ".join(p for p in [r.get("gene") or "", hg.get("c") or r["input"], hg.get("p") or ""] if p)
    if tx:
        head += f"  [{tx.get('basis')} {tx.get('ensembl')}]"
    lines.append(head)
    vcf = r.get("vcf")
    loc = f"{r['assembly']} {vcf['id']}" if vcf else f"{r['assembly']} (VCF form unavailable)"
    oa = r.get("other_assembly")
    if oa:
        loc += f"  ({oa['assembly']} {oa['variant_id']}, gnomAD liftover)"
    if r.get("rsids"):
        loc += "  " + ",".join(r["rsids"][:3])
    if r.get("caid"):
        loc += f"  {r['caid']}"
    lines.append(loc)
    g38 = r.get("grch38")
    if g38:
        lines.append(f"GRCh38 view: {g38['vcf']['id']} ({g38['mapped_by']}); predictors and gnomAD v4 below are from it "
                     "where GRCh37 has none")
    if hg.get("note"):
        lines.append(f"note: {hg['note']}")
    c = r.get("consequence") or {}
    where = ", ".join(x for x in [f"exon {c['exon']}" if c.get("exon") else "", f"intron {c['intron']}" if c.get("intron") else ""] if x)
    lines.append(f"consequence: {', '.join(c.get('terms') or []) or c.get('most_severe') or '-'} ({c.get('impact') or '-'})"
                 + (f", {where}" if where else ""))
    p = r.get("predictors") or {}
    am = p.get("alphamissense") or {}
    sp = p.get("spliceai") or {}
    sp_txt = "-"
    if sp:
        sp_txt = f"{_f(sp.get('max'))} (AG {_f(sp.get('DS_AG'))} AL {_f(sp.get('DS_AL'))} DG {_f(sp.get('DS_DG'))} DL {_f(sp.get('DS_DL'))})"
    elif p.get("spliceai_gnomad_max") is not None:
        sp_txt = f"{_f(p['spliceai_gnomad_max'])} (gnomAD)"
    revel = _f(p.get("revel")) if p.get("revel") is not None else (f"{_f(p['revel_gnomad'])} (gnomAD)" if p.get("revel_gnomad") is not None else "-")
    lines.append(f"predictors: REVEL {revel} | AlphaMissense {_f(am.get('score'))} {am.get('class') or ''}".rstrip()
                 + f" | CADD {_f(p.get('cadd_phred'))} | SpliceAI {sp_txt}"
                 + (f" | SIFT {p['sift']['prediction']}" if p.get("sift") else "")
                 + (f" | PolyPhen {p['polyphen']['prediction']}" if p.get("polyphen") else ""))
    pop = r.get("population")
    if pop and pop.get("mitochondrial"):
        if pop.get("found"):
            lines.append(f"gnomAD mtDNA ({pop['dataset']}): homoplasmic {pop.get('ac_hom')}/{pop.get('an')} "
                         f"(AF_hom {_f(pop.get('af_hom'))}), heteroplasmic {pop.get('ac_het')}/{pop.get('an')} "
                         f"(AF_het {_f(pop.get('af_het'))}), max heteroplasmy {_f(pop.get('max_heteroplasmy'))}"
                         + (f"; filters {','.join(pop['filters'])}" if pop.get("filters") else ""))
            pr = pop.get("predictors") or {}
            if pr.get("mitotip") or pr.get("pon_mt_trna"):
                lines.append(f"  tRNA predictors (gnomAD): MitoTIP {pr.get('mitotip')} ({_f(pr.get('mitotip_score'))}), "
                             f"PON-mt-tRNA {pr.get('pon_mt_trna')}")
            lines.append(f"  {pop.get('definitions')}")
        else:
            lines.append(f"gnomAD mtDNA ({pop['dataset']}): {pop.get('note')}")
    elif pop and "dataset" in pop:
        cov = pop.get("coverage") or {}
        # F33: the covered fraction decides whether absence means anything; mean
        # depth is shown after it, never instead of it.
        cov_txt = ", ".join(
            f"{k} {_pct(v.get('over_20'))} of samples >=20x (median {_f(v.get('median'))}, mean {_f(v.get('mean'))})"
            for k, v in cov.items())
        if pop.get("found"):
            t = pop.get("total") or {}
            g = pop.get("grpmax") or {}
            faf = pop.get("faf95") or {}
            lines.append(f"gnomAD ({pop['dataset']}): AC {t.get('ac')}/{t.get('an')} AF {_f(t.get('af'))}, hom {t.get('hom')}"
                         + (f", hemi {t.get('hemi')}" if t.get("hemi") else "")
                         + f"; grpmax {_f(g.get('af'))} ({g.get('group') or '-'}, AN {g.get('an') or '-'}); faf95 {_f(faf.get('value'))} ({faf.get('group') or '-'})"
                         + (f"; filters {','.join(pop['filters'])}" if pop.get("filters") else ""))
            if g.get("basis"):
                lines.append(f"  grpmax basis: {g['basis']}")
            if faf.get("basis"):
                lines.append(f"  faf95 basis: {faf['basis']}")
        else:
            lines.append(f"gnomAD ({pop['dataset']}): absent; {'site covered' if pop.get('covered') else 'coverage low/unknown'} ({cov_txt or 'no coverage data'})")
            if pop.get("coverage_rule"):
                lines.append(f"  coverage rule: {pop['coverage_rule']}")
    elif pop and pop.get("fallback"):
        fb = pop["fallback"]
        lines.append(f"gnomAD (via VEP, fallback): exome AF {_f(fb.get('exome_af'))}, genome AF {_f(fb.get('genome_af'))}, grpmax {_f(fb.get('grpmax_af'))} ({fb.get('grpmax_group')})")
    else:
        lines.append("gnomAD: unavailable")
    g38pop = (g38 or {}).get("population")
    if g38pop and g38pop.get("dataset"):
        t = g38pop.get("total") or {}
        g = g38pop.get("grpmax") or {}
        faf = g38pop.get("faf95") or {}
        if g38pop.get("found"):
            lines.append(f"gnomAD v4 on GRCh38 ({g38pop['dataset']}): AC {t.get('ac')}/{t.get('an')} AF {_f(t.get('af'))}; "
                         f"grpmax {_f(g.get('af'))} ({g.get('group') or '-'}); faf95 {_f(faf.get('value'))}")
        else:
            cov38 = ", ".join(f"{k} {_pct(v.get('over_20'))} of samples >=20x"
                              for k, v in (g38pop.get("coverage") or {}).items())
            state = {True: "site covered", False: "coverage low"}.get(g38pop.get("covered"), "coverage unknown")
            lines.append(f"gnomAD v4 on GRCh38 ({g38pop['dataset']}): not observed; {state} ({cov38 or 'no coverage data'})")
    china = r.get("population_china")
    if china:
        ds = china.get("datasets") or []
        lines.append("Chinese cohorts: " + ("; ".join(
            f"{d.get('name')} AF {_f(d.get('af'))} ({d.get('ac')}/{d.get('an')})" for d in ds[:4])
            or f"no record in {', '.join(china.get('checked') or []) or 'the cohorts checked'}"))
    mt = r.get("mitochondrial") or {}
    if mt.get("heteroplasmy"):
        h = mt["heteroplasmy"]
        lines.append(f"heteroplasmy (given): {h['percent']}% — {h['class']}"
                     + (f"; gnomAD max heteroplasmy {_f(h.get('gnomad_max_heteroplasmy'))}" if h.get("gnomad_max_heteroplasmy") is not None else ""))
    mm = mt.get("mitomap")
    if mm:
        if mm.get("found"):
            lines.append(f"MITOMAP (via MITOMASTER): {mm.get('locus') or '-'}; disease: {mm.get('disease_reported') or 'none listed'}; "
                         f"GenBank {mm.get('genbank_sequences_with_variant')} sequences ({_f(mm.get('genbank_percent'))}%)")
            lines.append(f"  {mm.get('status_note')}")
        else:
            lines.append(f"MITOMAP (via MITOMASTER): {mm.get('note')}")
    cv = r.get("clinvar")
    if cv:
        stars = cv.get("stars")
        star_txt = ("★" * stars + "☆" * (4 - stars)) if isinstance(stars, int) else "?"
        conds = "; ".join(x["name"] for x in (cv.get("conditions") or [])[:3] if x.get("name"))
        lines.append(f"ClinVar: {cv.get('classification')} {star_txt} {cv.get('stars_text') or ''} ({cv.get('review_status')}); {cv.get('submissions_scv')} SCV; "
                     f"last evaluated {cv.get('last_evaluated') or '-'}; {cv.get('vcv')} — {conds or '-'}")
        lines.append(f"  {cv.get('url')}")
    else:
        lines.append(f"ClinVar: {r.get('clinvar_note') or 'unavailable'}")
    for o in (r.get("clinvar_others_at_locus") or [])[:3]:
        lines.append(f"  other at locus: {o['vcv']} {o['title']} — {o['classification']}")
    lit = r.get("literature")
    if lit and lit.get("litvar"):
        recs = lit["litvar"].get("records") or []
        excl = lit["litvar"].get("excluded") or []
        if recs:
            lines.append("literature (LitVar2): " + "; ".join(
                f"{x.get('litvar_id')} {x.get('pmid_count')} PMIDs ({x.get('allele')})" for x in recs[:3])
                + f"  (query {lit.get('query')})")
        else:
            lines.append(f"literature (LitVar2): no record for this allele ({lit.get('query')})")
        for x in excl[:3]:
            lines.append(f"  excluded: {x.get('litvar_id')} {x.get('pmid_count')} PMIDs — {x.get('reason')}")
    fs = (r.get("functional_scores") or {}).get("mavedb")
    if fs:
        sc = fs.get("scores") or []

        def interp(s: Dict[str, Any]) -> str:
            i = s.get("interpretation")
            if isinstance(i, dict) and i.get("label"):
                crit = " ".join(x for x in (i.get("acmg_criterion"), i.get("acmg_evidence_strength")) if x)
                return f" ({i['label']}" + (f"; calibration says {crit}" if crit else "") + ")"
            if isinstance(i, str) and i:
                return f" ({i})"
            return " (no functional class)" if s.get("interpretation_note") else ""

        searched = fs.get("datasets_searched")
        n = len(searched) if isinstance(searched, list) else len(fs.get("datasets") or [])
        lines.append("MaveDB: " + ("; ".join(f"{s.get('urn')} {s.get('hgvs')} score {_f(s.get('score'))}{interp(s)}"
                                             for s in sc[:3])
                                   or f"no score for this variant in the {n} score set(s) searched"))
    a = r.get("acmg_inputs") or {}
    keys = ("gnomad_ac", "grpmax_af", "grpmax_an", "faf95", "revel", "spliceai_max", "is_missense")
    if a.get("mitochondrial"):
        keys = ("mt_af_hom", "mt_af_het", "mt_max_heteroplasmy", "heteroplasmy")
    lines.append("acmg_inputs: " + ", ".join(f"{k}={_f(a.get(k))}" for k in keys)
                 + (f" (frequency: {a['frequency_source']})" if a.get("frequency_source") else ""))
    return "\n".join(lines)


def _variant(args: argparse.Namespace) -> Outcome:
    from zebra.sources import variant as variant_src

    text = (args.query or "").strip()
    if not text:
        raise UsageError("give a variant: NM_...:c. HGVS, an rsID, or chrom-pos-ref-alt with --assembly")
    try:
        out = variant_src.card(text, assembly=args.assembly, gene=args.gene, heteroplasmy=args.heteroplasmy)
    except ValueError as err:
        raise UsageError(str(err)) from None
    out.text = render(out.result)
    out.query = {"variant": text, "assembly": args.assembly, "gene": args.gene, "heteroplasmy": args.heteroplasmy}
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("variant", help="annotate one variant: VEP (MANE), gnomAD, ClinVar, REVEL/AlphaMissense/CADD/SpliceAI, LitVar")
    p.add_argument("query", help="NM_000492.4:c.1521_1523del | rs113993960 | 7-117559590-ATCT-A")
    p.add_argument("--assembly", choices=("GRCh38", "GRCh37"), default="GRCh38",
                   help="genome build of genomic coordinates (default GRCh38)")
    p.add_argument("--gene", default=None, help="expected gene symbol (checked against VEP)")
    from zebra.commands.acmg import _heteroplasmy_arg

    p.add_argument("--heteroplasmy", type=_heteroplasmy_arg, default=None,
                   help="mtDNA heteroplasmy level, e.g. 35 or 35%% (also accepted after the variant: 'm.3243A>G 35%%')")
    p.set_defaults(func=_variant)
