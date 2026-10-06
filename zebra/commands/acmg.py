"""zebra acmg: ACMG/AMP points and combining rules; code suggestions from data."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, List, Optional, Tuple

from zebra import acmg, stats
from zebra.core import Outcome, UsageError

# D-P0-1: the Whiffin et al. 2017 formula follows from the inheritance. Whiffin publishes a
# monoallelic and a biallelic formula; X-linked disorders fit neither exactly (a hemizygous male
# is affected by one allele, a carrier female usually is not), so for them the caller chooses.
WHIFFIN_MODE = {"AD": "monoallelic", "AR": "biallelic"}

# Walker et al. 2023 calibrated SpliceAI on raw scores with a maximum distance of 10,000 nt
# (+/-4,999 from the variant); the splice lookup in `suggest` is run at those settings.
WALKER_DISTANCE = 4999
# Per-request timeout for the lookup inside `acmg suggest`: the service is interactive and the
# suggestion must not wait minutes for a splice score it can fall back from (CP1-18).
LOOKUP_TIMEOUT = 30.0
SPLICE_LOOKUP_TERMS = frozenset(("splice_region_variant", "splice_donor_region_variant",
                                 "splice_donor_5th_base_variant", "splice_polypyrimidine_tract_variant",
                                 "intron_variant", "synonymous_variant"))
NO_TRANSCRIPT_TERMS = frozenset(("intergenic_variant", "upstream_gene_variant", "downstream_gene_variant"))


def _classify(args: argparse.Namespace) -> Outcome:
    if not args.codes:
        raise UsageError("give evidence codes, e.g. PVS1 PM2_Supporting PP3")
    try:
        res = acmg.classify(args.codes)
    except ValueError as err:
        raise UsageError(str(err)) from None
    lines = [" + ".join(c["label"] for c in res["codes"]),
             f"points: {res['points'] if res['points'] is not None else 'BA1 (stand-alone)'} → {res['classification']}",
             f"2015 combining rules → {res['classification_richards_2015']}"]
    if res.get("classification_richards_2015_svi_2020") != res["classification_richards_2015"]:
        lines.append(f"2015 rules + ClinGen SVI 2020 (PVS1 + one Supporting → LP) → "
                     f"{res['classification_richards_2015_svi_2020']}")
    if res.get("pm1_pp3_cap_applied"):
        lines.append("counted after the PM1 + PP3 cap: " + " + ".join(c["label"] for c in res["codes_after_cap"]))
    if res.get("capped_at_uncertain"):
        lines.append(f"capped: the points total alone would read {res['capped_at_uncertain']}; computational evidence "
                     f"cannot classify on its own")
    if res.get("note"):
        lines.append(f"note: {res['note']}")
    for w in res["warnings"]:
        lines.append(f"warning: {w}")
    lines.append("research-grade: confirm with an accredited laboratory before clinical use")
    return Outcome(res, text="\n".join(lines), query={"codes": args.codes})


def whiffin_mode(inheritance: Optional[str], explicit: Optional[str]) -> Tuple[str, str]:
    """(formula, why) for the maximum credible AF; refuses a formula that contradicts the inheritance."""
    implied = WHIFFIN_MODE.get(inheritance or "")
    if implied and explicit and explicit != implied:
        raise UsageError(
            f"--inheritance {inheritance} means the {implied} Whiffin et al. 2017 formula; --inheritance-mode "
            f"{explicit} contradicts it. Drop --inheritance-mode (it follows from --inheritance) or correct the "
            f"inheritance")
    if implied:
        return implied, f"{implied} formula, from --inheritance {inheritance}"
    if explicit:
        why = f"{explicit} formula, chosen with --inheritance-mode"
        if inheritance in ("XLR", "XLD"):
            why += (f"; Whiffin et al. 2017 publishes no X-linked formula, so for {inheritance} this is the "
                    f"caller's approximation (biallelic gives the higher ceiling, i.e. the harder BS1)")
        return explicit, why
    raise UsageError(
        f"--prevalence needs the inheritance to choose the Whiffin et al. 2017 formula: --inheritance AD "
        f"(monoallelic) or AR (biallelic). For {inheritance or 'unknown'} inheritance pass --inheritance-mode "
        f"monoallelic|biallelic explicitly (no X-linked formula is published; biallelic is the conservative one "
        f"for BS1)")


def _max_credible_af(args: argparse.Namespace) -> Dict[str, Any]:
    """The Whiffin 2017 maximum credible AF from the flags, or why it is absent."""
    inheritance = None if args.inheritance == "unknown" else args.inheritance
    if args.prevalence is None:
        # flags that only make sense with --prevalence are refused rather than silently dropped
        if args.allelic is not None:
            raise UsageError("--allelic needs --prevalence (the Whiffin et al. 2017 maximum credible AF uses both)")
        if args.inheritance_mode and WHIFFIN_MODE.get(inheritance or "") not in (None, args.inheritance_mode):
            whiffin_mode(inheritance, args.inheritance_mode)  # raises: contradicts --inheritance
        return {
            "max_credible_af": None,
            "why_absent": ("not given: pass --prevalence and --allelic (with --inheritance AD/AR, plus --genetic and "
                           "--penetrance if not 1) so the Whiffin et al. 2017 maximum credible AF can be computed. "
                           "Without it BS1 cannot be assessed."),
        }
    if args.allelic is None:
        raise UsageError("--prevalence needs --allelic (the largest share of cases one allele can explain)")
    mode, why = whiffin_mode(inheritance, args.inheritance_mode)
    try:
        res = stats.max_credible_af(args.prevalence, args.allelic, args.genetic, args.penetrance, mode)
    except ValueError as err:
        raise UsageError(str(err)) from None
    return {"max_credible_af": res["max_credible_af"], "model": res["model"], "inputs": res["inputs"],
            "formula_basis": why}


def splice_relevance(r: Dict[str, Any], data: Dict[str, Any]) -> Optional[str]:
    """Why the splice lookup should run for this variant, or None when it should not."""
    cons = data.get("consequence") or []
    terms = {cons} if isinstance(cons, str) else set(cons)
    if data.get("mitochondrial") or not r.get("transcript"):
        return None
    if terms & set(acmg.NULL_CONSEQUENCES):
        return None  # null / +/-1,2: PVS1 covers the mechanism; Walker 2023's thresholds do not apply there
    if terms and terms <= NO_TRANSCRIPT_TERMS:
        return None
    hit = sorted(terms & SPLICE_LOOKUP_TERMS)
    if hit:
        return f"splice-relevant class ({', '.join(hit)})"
    vcf = r.get("vcf") or {}
    indel = bool(vcf) and (len(vcf.get("ref") or "") != 1 or len(vcf.get("alt") or "") != 1)
    if indel and data.get("spliceai_max") is None:
        return "an indel with no precomputed SpliceAI score from VEP"
    if data.get("spliceai_max") is None:
        return "no precomputed SpliceAI score from VEP for this variant"
    return None


def _lookup_coordinates(r: Dict[str, Any]) -> Optional[Tuple[str, str]]:
    v38 = (r.get("grch38") or {}).get("vcf")
    if r.get("assembly") == "GRCh37" and v38:
        return v38["id"], "GRCh38"
    vcf = r.get("vcf")
    if vcf:
        return vcf["id"], r.get("assembly") or "GRCh38"
    return None


def splice_evidence(r: Dict[str, Any], data: Dict[str, Any], enabled: bool = True) -> Outcome:
    """CP1-1 / CP1-18: the splicing score that drives PP3/BP4/BP7, and which model and run it came from.

    For a splice-relevant variant the Broad SpliceAI-lookup is called at Walker et al. 2023's
    calibration settings (raw scores, +/-4,999 nt) for SpliceAI and Pangolin; SpliceAI's max delta on
    the gene's transcript drives the codes, Pangolin is shown beside it (it has no ClinGen
    calibration). When the lookup cannot be reached, VEP's precomputed SpliceAI is used and the
    output says so.
    """
    vep_value = data.get("spliceai_max")
    info: Dict[str, Any] = {"vep_precomputed_spliceai": vep_value}
    why = splice_relevance(r, data)
    info["splice_relevant"] = why
    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []
    fallback_note = ("SpliceAI precomputed scores served by Ensembl VEP (Illumina's precomputed set, computed at "
                     "SpliceAI's default window rather than the +/-4,999 nt Walker et al. 2023 calibrated)")
    if why is None:
        info["used"] = "VEP precomputed SpliceAI" if vep_value is not None else None
        return Outcome(info)
    coords = _lookup_coordinates(r)
    if not enabled or coords is None:
        reason = "the splice lookup was switched off (--no-splice-lookup)" if not enabled else "no VCF coordinates"
        info["used"] = "VEP precomputed SpliceAI" if vep_value is not None else None
        info["note"] = f"{reason}; " + (f"{fallback_note} used" if vep_value is not None else "splicing axis not assessed")
        data["spliceai_source"] = fallback_note if vep_value is not None else None
        return Outcome(info, warnings=[info["note"]])
    from zebra import s2f

    vid, asm = coords
    try:
        pred = s2f.predict(vid, assembly=asm, models=["spliceai", "pangolin"], distance=WALKER_DISTANCE, mask=0,
                           lookup_timeout=LOOKUP_TIMEOUT)
    except Exception as err:  # noqa: BLE001 - any failure of the lookup falls back to VEP (CP1-18)
        pred = None
        failure = f"{type(err).__name__}: {err}"
    else:
        failure = None
        sources += pred.sources
    rows = {row.get("model"): row for row in (pred.result.get("models") if pred else []) or []}
    sp = rows.get("spliceai") or {}
    pg = rows.get("pangolin") or {}
    gene = r.get("gene")
    card_tx = {str(x).split(".")[0] for x in ((r.get("transcript") or {}).get("ensembl"),
                                              (r.get("transcript") or {}).get("refseq"),
                                              ((r.get("grch38") or {}).get("transcript") or {}).get("ensembl")) if x}

    def same_gene(row: Dict[str, Any]) -> List[Dict[str, Any]]:
        txs = row.get("transcripts") or []
        return [t for t in txs if gene and str(t.get("gene") or "").upper() == str(gene).upper()]

    def gene_transcript(row: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """The card's own transcript if the lookup scored it, else the gene's highest-priority one.

        Never another gene's transcript: a score for a neighbouring or antisense gene (CFTR-AS2 at
        CFTR c.3718-2477C>T) says nothing about this gene's splicing.
        """
        same = same_gene(row)
        own = [t for t in same if {str(t.get("transcript") or "").split(".")[0],
                                   str(t.get("refseq") or "").split(".")[0]} & card_tx]
        return (own or same or [None])[0]

    def delta(t: Dict[str, Any]) -> Optional[float]:
        return (t.get("max_delta") or {}).get("delta")

    if sp.get("status") == "ran" and gene_transcript(sp) is None:
        failure = f"the lookup scored no transcript of {gene}"
    if sp.get("status") == "ran" and gene_transcript(sp) is not None:
        tx = gene_transcript(sp) or {}
        top = tx.get("max_delta") or {}
        value = top.get("delta")
        if value is not None:
            where = f"{top.get('score')} {value:.3f}" + (f" at position {top['position']}" if top.get("position") else "")
            src = (f"SpliceAI via the Broad SpliceAI-lookup, raw scores, +/-{WALKER_DISTANCE} nt (Walker et al. 2023's "
                   f"calibration settings), {asm} {vid}, {where} on {tx.get('refseq') or tx.get('transcript')} ({tx.get('gene')})")
            data["spliceai_max"] = float(value)
            data["spliceai_source"] = src
            info.update(used="SpliceAI-lookup", spliceai={"max_delta": value, "score": top.get("score"),
                                                           "position": top.get("position"),
                                                           "transcript": tx.get("transcript"), "refseq": tx.get("refseq"),
                                                           "gene": tx.get("gene")},
                        source=src)
            if vep_value is not None and abs(float(vep_value) - float(value)) >= 0.1:
                warnings.append(f"SpliceAI: the lookup at +/-{WALKER_DISTANCE} nt gives {value:.3f}, VEP's precomputed "
                                f"score is {vep_value}; the lookup (the calibrated settings) is used")
            # A benign splicing reading needs no transcript of the gene to predict an effect: when another
            # transcript of the gene, or VEP's own precomputed score, reaches the PP3 threshold, BP4/BP7
            # from splicing are withheld rather than chosen by transcript.
            others = [(t, delta(t)) for t in same_gene(sp) if t is not tx and delta(t) is not None]
            hi = [t for t, d in others if d >= acmg.SPLICEAI_PP3]
            if float(value) <= acmg.SPLICEAI_BP4 and (hi or (vep_value is not None and float(vep_value) >= acmg.SPLICEAI_PP3)):
                why = (f"another {gene} transcript ({hi[0].get('refseq') or hi[0].get('transcript')}) scores "
                       f"{delta(hi[0]):.3f}" if hi else f"VEP's precomputed SpliceAI is {vep_value}")
                data["splice_benign_blocked"] = (f"the lookup gives {value:.3f} on {tx.get('refseq') or tx.get('transcript')}, "
                                                 f"but {why} (>= {acmg.SPLICEAI_PP3}); a splice effect is not excluded")
                info["benign_withheld"] = data["splice_benign_blocked"]
    if info.get("used") != "SpliceAI-lookup":
        reason = failure or sp.get("reason") or "no SpliceAI score for the gene's transcript"
        info["lookup_failure"] = reason
        if vep_value is not None:
            data["spliceai_source"] = fallback_note
            info["used"] = "VEP precomputed SpliceAI"
            info["note"] = f"SpliceAI-lookup unavailable ({reason}); {fallback_note} used instead"
        else:
            data["spliceai_max"] = None
            info["used"] = None
            info["note"] = (f"SpliceAI-lookup unavailable ({reason}) and VEP has no precomputed score: the splicing "
                            "axis (PP3/BP4/BP7 from splicing) was not assessed — not a 'no effect'")
        warnings.append(info["note"])
    if pg.get("status") == "ran":
        tx = gene_transcript(pg) or {}
        top = tx.get("max_delta") or {}
        info["pangolin"] = {"score": top.get("score"), "delta": top.get("delta"), "position": top.get("position"),
                            "transcript": tx.get("refseq") or tx.get("transcript"),
                            "note": ("shown beside SpliceAI, not used for a code: Pangolin has no ClinGen-calibrated "
                                     "PP3/BP4 thresholds and is not independent of SpliceAI (Zeng & Li 2022)")}
    elif pred is not None:
        info["pangolin"] = {"status": pg.get("status"), "reason": pg.get("reason")}
    return Outcome(info, sources=sources, warnings=warnings)


def _fraction(raw: str) -> float:
    """argparse type: a finite number in (0, 1] (a frequency, prevalence or share)."""
    try:
        v = float(raw)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError(f"{raw!r} is not a number") from None
    if not (v == v and 0 < v <= 1):  # v == v rejects NaN; the bounds reject inf
        raise argparse.ArgumentTypeError(f"{raw!r} must be a number in (0, 1]")
    return v


def _heteroplasmy_arg(raw: str) -> float:
    t = str(raw).strip().rstrip("%").strip()
    try:
        v = float(t)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{raw!r} is not a percentage") from None
    if not (v == v and 0 < v <= 100):
        raise argparse.ArgumentTypeError("heteroplasmy is a percentage in (0, 100]")
    return v / 100.0


def _suggest(args: argparse.Namespace) -> Outcome:
    from zebra.sources import variant as variant_src

    maxaf = _max_credible_af(args)  # bad flags are refused before any request is made
    try:
        # the Chinese-cohort and MaveDB sections feed no code here: skipped, the card shows them
        card = variant_src.card(args.variant, assembly=args.assembly, heteroplasmy=args.heteroplasmy,
                                contracts=False)
    except UsageError:
        raise
    except ValueError as err:
        # a REF that does not match the build, an unreadable variant string: the
        # caller's input, not an internal fault
        raise UsageError(str(err)) from None
    r = card.result
    data: Dict[str, Any] = dict(r.get("acmg_inputs", {}))
    # The caller may know ClinVar's aggregate classification already; the variant card
    # carries it. BA1 and BS1 must never be offered against a P/LP assertion without saying so.
    clinvar = r.get("clinvar") or {}
    if isinstance(clinvar, dict) and not data.get("clinvar_classification"):
        assertion = clinvar.get("germline_classification") or clinvar.get("classification") or clinvar.get("clinical_significance")
        if assertion:
            data["clinvar_classification"] = assertion

    sources = list(card.sources)
    warnings = list(card.warnings)
    with ThreadPoolExecutor(max_workers=2) as pool:
        f_sp = pool.submit(splice_evidence, r, data, not args.no_splice_lookup)
        f_ji = pool.submit(variant_src.judgement_inputs, r)
        try:
            sp = f_sp.result()
        except Exception as err:  # noqa: BLE001 - never lose the suggestions over the splice lookup
            sp = Outcome({"used": None, "note": f"splice evidence failed: {type(err).__name__}: {err}"},
                         warnings=[f"splice evidence failed: {type(err).__name__}: {err}"])
        try:
            ji = f_ji.result()
        except Exception as err:  # noqa: BLE001
            ji = Outcome({}, warnings=[f"PVS1/PS1/PM5 inputs failed: {type(err).__name__}: {err}"])
    sources += sp.sources + ji.sources
    warnings += sp.warnings + ji.warnings

    computed: Optional[float] = maxaf["max_credible_af"]
    inheritance = None if args.inheritance == "unknown" else args.inheritance
    res_s = acmg.suggest(data, inheritance=inheritance, max_credible_af=computed, pm2_max_af=args.pm2_max_af)
    suggestions = res_s["suggested"]
    res = {
        "variant": r.get("variant"),
        "inputs": data,
        "max_credible_af": maxaf,
        "splicing": sp.result,
        "suggested": suggestions,
        "not_assessed": res_s["not_assessed"],
        "caveats": res_s["caveats"],
        "not_from_data": [
            "PVS1 (null variant in a gene where loss of function causes disease: ClinGen PVS1 decision tree, "
            "Abou Tayoun 2018 — its inputs are in pvs1_inputs)",
            "PS1/PM5 (same amino-acid change / another missense at the residue known pathogenic — the ClinVar "
            "records at the codon are in ps1_pm5_inputs)",
            "PS2/PM6 (de novo, confirmed or assumed)", "PS3/BS3 (functional studies, calibrated per Brnich 2019)",
            "PS4 (case-control enrichment)", "PM1 (hotspot / critical domain)", "PM3/BP2 (in trans / in cis)",
            "PM4/BP3 (in-frame length change)", "PP1/BS4 (segregation: zebra stats segregation)", "PP4 (phenotype specific to the gene)",
        ],
        "note": "Suggestions cover only codes that follow from numbers with ClinGen-calibrated thresholds; each still needs a look at coverage, transcript and gene context.",
    }
    for key in ("pvs1_inputs", "ps1_pm5_inputs"):
        if key in ji.result:
            res[key] = ji.result[key]
    if ji.result.get("note"):
        res["judgement_inputs_note"] = ji.result["note"]
    if r.get("mitochondrial"):
        res["mitochondrial"] = r["mitochondrial"]
    if r.get("builds"):
        res["builds"] = r["builds"]

    lines = [f"variant: {res['variant']}"]
    if r.get("builds"):
        b = r["builds"]
        lines.append(f"builds: GRCh37 {b['GRCh37']['id']} | GRCh38 {b['GRCh38']['id']} ({b['GRCh38']['mapped_by']})")
    for s in suggestions:
        lines.append(f"  {s['code']:<16} {s['basis']}\n                   rule: {s['rule']}")
        if s.get("conflict"):
            lines.append(f"                   CONFLICT: {s['conflict']}")
        for c in s.get("caveats") or []:
            lines.append(f"                   caveat: {c}")
    if not suggestions:
        lines.append("  no code follows from the numbers alone")
    if computed is not None:
        lines.append(f"maximum credible AF: {computed:.4g} ({maxaf['model']}; {maxaf['formula_basis']})")
    else:
        lines.append(f"maximum credible AF: {maxaf['why_absent']}")
    spl = sp.result or {}
    if spl.get("used") == "SpliceAI-lookup":
        lines.append(f"splicing score used: {spl.get('source')}")
    elif spl.get("note"):
        lines.append(f"splicing: {spl['note']}")
    if (spl.get("pangolin") or {}).get("delta") is not None:
        pgl = spl["pangolin"]
        lines.append(f"  Pangolin {pgl.get('score')} {pgl['delta']:+.3f} on {pgl.get('transcript')} (shown, not used for a code)")
    pv = res.get("pvs1_inputs")
    if pv:
        nmd = pv.get("nmd") or {}
        if nmd:
            lines.append(f"PVS1 inputs: stop at codon {nmd.get('ptc_codon')} in exon {nmd.get('ptc_exon')}, "
                         f"{nmd.get('position')}; NMD predicted: {nmd.get('nmd_predicted')}"
                         + (f"; {nmd['fraction_of_protein_after_ptc']:.0%} of the protein follows the stop"
                            if nmd.get("fraction_of_protein_after_ptc") is not None else ""))
        css = pv.get("canonical_splice_site")
        if css:
            lines.append(f"PVS1 inputs: canonical {css['site']} site of exon {css['exon_affected']} "
                         f"({css['exon_coding_length']} coding nt; skipping keeps frame: {css['skipping_keeps_frame']})")
        mech = pv.get("lof_mechanism") or {}
        lines.append(f"PVS1 inputs: LoF mechanism — ClinGen HI {mech.get('clingen_hi_score')} "
                     f"({mech.get('clingen_hi_description') or mech.get('clingen_note') or '-'}); "
                     f"LOEUF {mech.get('gnomad_loeuf')}, pLI {mech.get('gnomad_pli')}")
    cr = res.get("ps1_pm5_inputs")
    if cr:
        def fmt(rows: List[Dict[str, Any]]) -> str:
            return "; ".join(f"{x['vcv']} {x['title']} — {x['classification']}, {x.get('stars_text') or x.get('stars')}"
                             for x in rows[:4]) or "none"
        lines.append(f"PS1 inputs (same change {cr['change']}, ClinVar P/LP): {fmt(cr['ps1_inputs'])}")
        lines.append(f"PM5 inputs (other missense at {cr['residue']}, ClinVar P/LP): {fmt(cr['pm5_inputs'])}")
    for n in res["not_assessed"]:
        lines.append(f"not assessed: {n}")
    for c in res["caveats"]:
        lines.append(f"caveat: {c}")
    lines.append("judgement codes to consider: " + "; ".join(res["not_from_data"]))
    out = Outcome(res, sources=sources, warnings=warnings, text="\n".join(lines),
                  query={"variant": args.variant, "inheritance": args.inheritance,
                         "prevalence": args.prevalence, "allelic": args.allelic, "genetic": args.genetic,
                         "penetrance": args.penetrance, "inheritance_mode": args.inheritance_mode,
                         "pm2_max_af": args.pm2_max_af, "heteroplasmy": args.heteroplasmy,
                         "splice_lookup": not args.no_splice_lookup})
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("acmg", help="ACMG/AMP: points + combining rules (classify); data-driven code suggestions (suggest)")
    a = p.add_subparsers(dest="action", metavar="<action>")
    q = a.add_parser("classify", help="classify from evidence codes")
    q.add_argument("codes", nargs="*")
    q.set_defaults(func=_classify, no_ledger=True)
    q = a.add_parser("suggest", help="codes that follow from frequency and calibrated predictors")
    q.add_argument("variant")
    q.add_argument("--assembly", choices=("GRCh38", "GRCh37"), default="GRCh38")
    q.add_argument("--inheritance", choices=("AD", "AR", "XLD", "XLR", "unknown"), default="unknown")
    # BS1 needs the Whiffin et al. 2017 maximum credible AF; these are its inputs.
    q.add_argument("--prevalence", type=_fraction,
                   help="disease prevalence (affected fraction), e.g. 0.0004 for 1 in 2,500; needed for BS1")
    q.add_argument("--allelic", type=_fraction,
                   help="largest share of cases one allele can explain (maximum allelic contribution)")
    q.add_argument("--genetic", type=_fraction, default=1.0, help="share of cases due to this gene (default 1.0)")
    q.add_argument("--penetrance", type=_fraction, default=1.0, help="penetrance of the genotype (default 1.0)")
    q.add_argument("--inheritance-mode", choices=("monoallelic", "biallelic"), default=None,
                   help="Whiffin 2017 formula; follows from --inheritance (AD monoallelic, AR biallelic) and is "
                        "only needed for XLR/XLD/unknown; a value contradicting --inheritance is refused")
    q.add_argument("--pm2-max-af", type=_fraction, default=None,
                   help="gene-specific PM2 frequency ceiling from the gene's ClinGen VCEP specification (grpmax AF)")
    q.add_argument("--heteroplasmy", type=_heteroplasmy_arg, default=None,
                   help="mtDNA heteroplasmy level, e.g. 35 or 35%% (also accepted after the variant: 'm.3243A>G 35%%')")
    q.add_argument("--no-splice-lookup", action="store_true",
                   help="do not call the SpliceAI-lookup service; use VEP's precomputed SpliceAI only")
    q.set_defaults(func=_suggest)
