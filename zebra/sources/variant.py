"""Variant card: one variant, normalised and annotated from primary sources.

Ensembl VEP (consequence on the MANE Select transcript; REVEL, AlphaMissense,
CADD, SpliceAI on GRCh38) -> VCF-style left-normalised form checked against the
Ensembl reference -> gnomAD (frequencies, coverage), ClinVar (exact-allele
match) and LitVar (mention count) -> `acmg_inputs` for `zebra acmg suggest`.

Coordinates: VEP reports start/end on the forward strand but, for HGVS input
on a minus-strand transcript, the alleles on the transcript strand
(strand = -1). Alleles are reverse-complemented before the VCF form is built,
then anchored and left-aligned against the reference sequence, which is also
used to check REF.
"""

from __future__ import annotations

import re
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Dict, List, Optional, Tuple

from zebra.core import Outcome, attempt
from zebra.sources import clinvar, ensembl, gnomad

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
LEFT_PAD = 100
RIGHT_PAD = 20
CONTINENTAL = ("afr", "amr", "eas", "nfe", "sas")


class OutOfWindow(Exception):
    """A base outside the fetched reference window was needed."""


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


class Window:
    """Forward-strand reference sequence chrom:start..end (1-based, inclusive)."""

    def __init__(self, chrom: str, start: int, seq: str):
        self.chrom, self.start, self.seq = chrom, start, seq.upper()
        self.end = start + len(seq) - 1

    def base(self, p: int) -> str:
        if p < self.start or p > self.end:
            raise OutOfWindow(p)
        return self.seq[p - self.start]

    def slice(self, a: int, b: int) -> str:
        if a < self.start or b > self.end:
            raise OutOfWindow((a, b))
        return self.seq[a - self.start:b - self.start + 1]


def left_normalize(pos: int, ref: str, alt: str, base_at: Callable[[int], str]) -> Tuple[int, str, str]:
    """Parsimonious, left-aligned VCF representation (Tan et al. 2015)."""
    ref, alt = ref.upper(), alt.upper()
    if ref == alt:
        raise ValueError("REF equals ALT")
    for _ in range(100000):
        changed = False
        if not ref or not alt:
            b = base_at(pos - 1)
            ref, alt, pos = b + ref, b + alt, pos - 1
            changed = True
        if ref and alt and ref[-1] == alt[-1] and (len(ref) > 1 or len(alt) > 1):
            ref, alt = ref[:-1], alt[:-1]
            changed = True
        if not changed:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def vep_alleles(rec: Dict[str, Any]) -> Tuple[str, List[str]]:
    parts = str(rec.get("allele_string") or "").split("/")
    if len(parts) < 2:
        raise ValueError(f"VEP allele_string {rec.get('allele_string')!r} has no alternate allele")
    return parts[0], parts[1:]


def raw_vcf(rec: Dict[str, Any], alt_raw: str, base_at: Callable[[int], str]) -> Tuple[int, str, str, Optional[str]]:
    """VEP record + one alternate allele -> (pos, ref, alt) on the forward strand, before normalisation.

    The fourth value is a REF-check message when the forward-strand REF that VEP
    implies differs from the reference sequence at start..end, else None.
    """
    ref_raw, _ = vep_alleles(rec)
    strand = rec.get("strand", 1)
    ref = "" if ref_raw == "-" else ref_raw.upper()
    alt = "" if alt_raw == "-" else alt_raw.upper()
    if strand == -1:
        ref, alt = revcomp(ref), revcomp(alt)
    start, end = int(rec["start"]), int(rec["end"])
    if not ref:  # insertion between end and start (start = end + 1)
        b = base_at(end)
        return end, b, b + alt, None
    found = "".join(base_at(p) for p in range(start, end + 1))
    mismatch = None if found == ref else f"VEP implies REF {ref} at {start}-{end}, reference has {found}"
    if not alt:  # deletion of start..end
        b = base_at(start - 1)
        return start - 1, b + ref, b, mismatch
    return start, ref, alt, mismatch


def parse_spdi(spdi: str) -> Tuple[str, int, str, str]:
    acc, pos0, dele, ins = spdi.split(":")
    return acc, int(pos0), dele, ins


def spdi_to_vcf(pos0: int, dele: str, ins: str, base_at: Callable[[int], str]) -> Tuple[int, str, str]:
    if dele.isdigit():
        n = int(dele)
        dele = "".join(base_at(p) for p in range(pos0 + 1, pos0 + n + 1))
    if dele and ins:
        return pos0 + 1, dele.upper(), ins.upper()
    b = base_at(pos0)  # the base before the event (1-based pos0)
    return pos0, (b + dele).upper(), (b + ins).upper()


def _num(x: Any) -> Optional[float]:
    if isinstance(x, (int, float)):
        return float(x)
    if isinstance(x, str):
        for part in x.split(","):
            try:
                return float(part)
            except ValueError:
                continue
    return None


def pick_transcript(tcs: List[Dict[str, Any]], gene: Optional[str] = None) -> Optional[Dict[str, Any]]:
    pool = tcs
    if gene:
        same = [t for t in tcs if str(t.get("gene_symbol", "")).upper() == gene.upper()]
        if same:
            pool = same
    for pred in (lambda t: t.get("mane_select"), lambda t: t.get("canonical") and t.get("biotype") == "protein_coding",
                 lambda t: t.get("canonical"), lambda t: t.get("biotype") == "protein_coding", lambda t: True):
        for t in pool:
            if pred(t):
                return t
    return None


def predictors(tc: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    tc = tc or {}
    am = tc.get("alphamissense") if isinstance(tc.get("alphamissense"), dict) else {}
    sp = tc.get("spliceai") if isinstance(tc.get("spliceai"), dict) else None
    splice = None
    if sp:
        splice = {k: sp.get(k) for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")}
        splice["max"] = ensembl.spliceai_max(tc)
        splice["gene"] = sp.get("SYMBOL")
    out = {
        "revel": _num(tc.get("revel", tc.get("revel_score"))),
        "alphamissense": {"score": _num(am.get("am_pathogenicity")), "class": am.get("am_class")} if am else None,
        "cadd_phred": _num(tc.get("cadd_phred")),
        "spliceai": splice,
        "sift": {"prediction": tc.get("sift_prediction"), "score": _num(tc.get("sift_score"))} if tc.get("sift_prediction") else None,
        "polyphen": {"prediction": tc.get("polyphen_prediction"), "score": _num(tc.get("polyphen_score"))} if tc.get("polyphen_prediction") else None,
        "source": "Ensembl VEP plugins on the selected transcript",
    }
    return out


def vep_frequencies(rec: Dict[str, Any], alt_raw: str) -> Optional[Dict[str, Any]]:
    """gnomAD frequencies VEP attaches to the colocated variant (fallback when gnomAD is down)."""
    alt_key = alt_raw if alt_raw else "-"
    for cv in rec.get("colocated_variants") or []:
        freqs = cv.get("frequencies") or {}
        groups = freqs.get(alt_key) or (next(iter(freqs.values())) if len(freqs) == 1 else None)
        if not groups:
            continue
        cont = {}
        for k, v in groups.items():
            m = re.match(r"^gnomad([eg])_(\w+)$", k)
            if m and m.group(2) in CONTINENTAL and isinstance(v, (int, float)):
                cont[f"{m.group(1)}:{m.group(2)}"] = v
        best = max(cont.items(), key=lambda kv: kv[1]) if cont else None
        return {"rsid": cv.get("id"), "exome_af": groups.get("gnomade"), "genome_af": groups.get("gnomadg"),
                "grpmax_af": best[1] if best else None, "grpmax_group": best[0] if best else None,
                "note": "from VEP colocated_variants (gnomAD exome 'e:'/genome 'g:' groups; no allele numbers)"}
    return None


def _clinvar_candidates(rec: Dict[str, Any]) -> Tuple[List[str], Optional[str]]:
    vcvs, rsid = [], None
    for cv in rec.get("colocated_variants") or []:
        cid = str(cv.get("id") or "")
        if cid.startswith("rs") and rsid is None:
            rsid = cid
        for s in (cv.get("var_synonyms") or {}).get("ClinVar", []) or []:
            if str(s).startswith("VCV"):
                vcvs.append(str(s))
    return list(dict.fromkeys(vcvs))[:10], rsid


def match_clinvar(records: List[Dict[str, Any]], vcf: Optional[Tuple[str, int, str, str]], assembly: str,
                  window: Optional[Window], caid: Optional[str]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Split ClinVar records into those for exactly this allele and the others at the locus."""
    exact, others = [], []
    for r in records:
        if r.get("compound"):
            continue
        hit = False
        if caid and r.get("caid") and r["caid"] == caid:
            hit = True
        elif vcf and r.get("canonical_spdi"):
            try:
                _, pos0, dele, ins = parse_spdi(r["canonical_spdi"])
                if assembly == "GRCh37":
                    l38, l37 = r["locations"].get("GRCh38"), r["locations"].get("GRCh37")
                    if not (l38 and l37 and l38.get("start") and l37.get("start")):
                        raise ValueError("no GRCh37 location")
                    pos0 += l37["start"] - l38["start"]
                if window is not None:
                    got = left_normalize(*spdi_to_vcf(pos0, dele, ins, window.base), window.base)
                elif len(dele) == 1 and len(ins) == 1:
                    got = (pos0 + 1, dele.upper(), ins.upper())
                else:
                    raise ValueError("no reference window")
                hit = got == (vcf[1], vcf[2], vcf[3])
            except (ValueError, KeyError, OutOfWindow):
                hit = False
        (exact if hit else others).append(r)
    exact.sort(key=lambda r: -(r.get("stars") or 0))
    return exact, others


def _litvar(text: str, gene: Optional[str]) -> Outcome:
    from zebra.sources import litvar  # owned by another work package; may be absent

    return litvar.lookup(text, gene=gene)


def _litvar_compact(res: Any, rsids: List[str]) -> Optional[Dict[str, Any]]:
    """LitVar records for this variant: LitVar's first suggestion per spelling, and any with our rsID."""
    if not isinstance(res, dict):
        return None
    matches = res.get("matches")
    if not isinstance(matches, list):
        return None
    ours = [m for m in matches if isinstance(m, dict) and (m.get("top") or (m.get("rsid") and m.get("rsid") in rsids))]
    recs = [{"litvar_id": m.get("litvar_id"), "rsid": m.get("rsid"), "name": m.get("name"), "hgvs": m.get("hgvs"),
             "pmid_count": m.get("pmid_count")} for m in ours[:4]]
    counts = [r["pmid_count"] for r in recs if isinstance(r.get("pmid_count"), int)]
    return {"records": recs, "pmid_count_max": max(counts) if counts else None,
            "note": "LitVar keeps unlinked spellings as separate records; counts overlap, do not add them"}


def card(variant: str, assembly: str = "GRCh38", gene: Optional[str] = None) -> Outcome:
    text = variant.strip()
    if assembly not in ("GRCh38", "GRCh37"):
        raise ValueError("assembly must be GRCh38 or GRCh37")
    kind = ensembl.classify_input(text)
    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []

    v = ensembl.vep(text, assembly)  # the anchor of the card: failure propagates
    rec = v.result
    sources += v.sources
    warnings += v.warnings
    ref_raw, alts_raw = vep_alleles(rec)
    alt_raw = alts_raw[0]
    if len(alts_raw) > 1:
        warnings.append(f"{text} has {len(alts_raw)} alternate alleles ({'/'.join(alts_raw)}); annotated {alt_raw} — "
                        "give HGVS or chrom-pos-ref-alt to choose another")
    chrom = str(rec.get("seq_region_name"))

    # reference window: REF check, anchoring and left-alignment
    window: Optional[Window] = None
    vcf: Optional[Tuple[str, int, str, str]] = None
    start, end = int(rec["start"]), int(rec["end"])
    for pad in (LEFT_PAD, 2000):
        lo = max(1, min(start, end) - pad)
        hi = max(start, end) + RIGHT_PAD
        seq = attempt("Ensembl reference sequence", lambda: ensembl.sequence(chrom, lo, hi, assembly), warnings)
        if seq is None:
            break
        sources += seq.sources
        window = Window(chrom, lo, seq.result)
        try:
            pos0, ref0, alt0, mismatch = raw_vcf(rec, alt_raw, window.base)
            if mismatch:
                warnings.append(f"reference mismatch on chr{chrom} ({assembly}): {mismatch} — check the assembly and the HGVS")
            npos, nref, nalt = left_normalize(pos0, ref0, alt0, window.base)
            vcf = (chrom, npos, nref, nalt)
            break
        except OutOfWindow:
            vcf = None
            continue
    if vcf is None:
        if len(ref_raw) == 1 and len(alt_raw) == 1 and "-" not in (ref_raw, alt_raw):
            r, a = (revcomp(ref_raw), revcomp(alt_raw)) if rec.get("strand") == -1 else (ref_raw, alt_raw)
            vcf = (chrom, start, r.upper(), a.upper())
            warnings.append("REF not checked against the reference sequence (Ensembl sequence unavailable)")
        else:
            warnings.append("could not build the VCF form of this indel (reference sequence unavailable): "
                            "gnomAD lookup skipped; ClinVar allele match unverified")

    # transcript
    tcs_all = rec.get("transcript_consequences") or []
    tcs = [t for t in tcs_all if t.get("variant_allele") in (alt_raw, None)] or tcs_all
    tc = pick_transcript(tcs, gene)
    gene_symbol = (tc or {}).get("gene_symbol")
    if gene and gene_symbol and gene.upper() != str(gene_symbol).upper():
        genes_hit = sorted({t.get("gene_symbol") for t in tcs if t.get("gene_symbol")})
        warnings.append(f"gene mismatch: you gave {gene}; VEP places this variant in {', '.join(genes_hit) or '-'} "
                        f"(annotated on {gene_symbol}) — check the transcript and assembly")
    elif gene and not tc:
        warnings.append(f"gene {gene} given, but VEP reports no transcript consequence")

    hgvs: Dict[str, Any] = {}
    transcript = None
    consequence: Dict[str, Any] = {"most_severe": rec.get("most_severe_consequence")}
    if tc:
        mane = tc.get("mane_select")
        hgvsc = tc.get("hgvsc")
        cpart = hgvsc.split(":", 1)[1] if hgvsc and ":" in hgvsc else None
        ppart = tc["hgvsp"].split(":", 1)[1] if tc.get("hgvsp") and ":" in tc["hgvsp"] else None
        transcript = {"ensembl": hgvsc.split(":")[0] if hgvsc else tc.get("transcript_id"), "refseq": mane,
                      "mane_select": bool(mane or "MANE_Select" in (tc.get("mane") or [])),
                      "canonical": bool(tc.get("canonical")), "gene": gene_symbol, "hgnc_id": tc.get("hgnc_id"),
                      "basis": "MANE Select" if mane else ("Ensembl canonical" if tc.get("canonical") else "first protein-coding")}
        if transcript["hgnc_id"] is not None and not str(transcript["hgnc_id"]).startswith("HGNC:"):
            transcript["hgnc_id"] = f"HGNC:{transcript['hgnc_id']}"
        hgvs = {"c": f"{mane}:{cpart}" if mane and cpart else hgvsc, "c_ensembl": hgvsc,
                "p": ppart, "p_ensembl": tc.get("hgvsp")}
        if kind == "hgvs" and cpart:
            acc_in, _, change_in = text.partition(":")
            if mane and acc_in.split(".")[0].upper() == mane.split(".")[0].upper():
                notes = []
                if acc_in.upper() != mane.upper():
                    notes.append(f"input transcript version {acc_in}, MANE Select is {mane}")
                if change_in != cpart:
                    notes.append(f"VEP's normalised form is {mane}:{cpart} (input {text})")
                if notes:
                    hgvs["note"] = "; ".join(notes)
            elif mane and acc_in.upper().startswith(("NM_", "ENST")) and acc_in.split(".")[0] != (hgvsc or "").split(".")[0]:
                hgvs["note"] = f"input is on {acc_in}; MANE Select is {mane} ({hgvs['c']})"
        consequence.update({"terms": tc.get("consequence_terms") or [], "impact": tc.get("impact"),
                            "exon": tc.get("exon"), "intron": tc.get("intron"),
                            "protein_position": tc.get("protein_start"), "amino_acids": tc.get("amino_acids"),
                            "codons": tc.get("codons"), "biotype": tc.get("biotype")})
        if not mane:
            warnings.append(f"no MANE Select transcript in VEP output ({assembly}); annotated {transcript['basis']} "
                            f"{transcript['ensembl']}" + (" — GRCh37 VEP does not report MANE" if assembly == "GRCh37" else ""))
    else:
        consequence["terms"] = [rec.get("most_severe_consequence")] if rec.get("most_severe_consequence") else []
    preds = predictors(tc)
    if assembly == "GRCh37":
        warnings.append("GRCh37 VEP serves no AlphaMissense or SpliceAI scores; for splicing use `zebra s2f` "
                        "(or annotate the GRCh38 coordinates)")
    elif tc and preds["spliceai"] is None:
        warnings.append("SpliceAI: no precomputed score from VEP for this variant (VEP serves SNVs and short indels); see `zebra s2f`")
    if preds["spliceai"] and preds["spliceai"].get("gene") and gene_symbol and preds["spliceai"]["gene"] != gene_symbol:
        warnings.append(f"SpliceAI scores are for {preds['spliceai']['gene']}, not {gene_symbol}")

    vcvs, rsid = _clinvar_candidates(rec)
    rsids = [cv.get("id") for cv in rec.get("colocated_variants") or [] if str(cv.get("id", "")).startswith("rs")]

    # gnomAD, ClinVar, LitVar in parallel (different hosts)
    w_gn: List[str] = []
    w_cv: List[str] = []
    w_lv: List[str] = []

    def do_gnomad():
        if not vcf:
            return None
        return attempt("gnomAD", lambda: gnomad.variant(vcf[0], vcf[1], vcf[2], vcf[3], assembly), w_gn)

    def do_clinvar():
        if not (vcvs or rsid or vcf or kind == "hgvs"):
            w_cv.append("ClinVar not searched: no VCV id, rsID, HGVS or position to search with")
            return None

        def run():
            pos_q = None
            if vcf:
                p = vcf[1] + 1 if (len(vcf[2]) > len(vcf[3]) and vcf[2][0] == vcf[3][0]) else vcf[1]
                pos_q = p
            term = clinvar.term_for(vcv=vcvs, rsid=rsid, hgvs=text if kind == "hgvs" else None,
                                    chrom=chrom if pos_q else None, pos=pos_q, assembly=assembly)
            found = clinvar.search(term)
            out = clinvar.summaries(found.result["ids"])
            out.sources = found.sources + out.sources
            return out
        return attempt("ClinVar", run, w_cv)

    lit_text = hgvs.get("p") if (gene_symbol and hgvs.get("p")) else (rsid or hgvs.get("c") or text)

    def do_litvar():
        try:
            return _litvar(lit_text, gene_symbol)
        except Exception as err:  # noqa: BLE001 - sibling module may be missing or change shape
            w_lv.append(f"LitVar unavailable: {type(err).__name__}: {err}")
            return None

    with ThreadPoolExecutor(max_workers=3) as pool:
        f_gn, f_cv, f_lv = pool.submit(do_gnomad), pool.submit(do_clinvar), pool.submit(do_litvar)
        gn, cv, lv = f_gn.result(), f_cv.result(), f_lv.result()
    warnings += w_gn + w_cv + w_lv

    # population
    population: Optional[Dict[str, Any]] = None
    if gn is not None:
        sources += gn.sources
        warnings += gn.warnings
        population = gn.result
    fallback = None
    if population is None:
        fallback = vep_frequencies(rec, alt_raw)
        if fallback:
            warnings.append("gnomAD API unavailable: frequencies below come from VEP's colocated-variant data "
                            "(no allele numbers, no coverage; BA1/BS1 need the allele number — recheck in gnomAD)")
        else:
            warnings.append("no population frequency available (gnomAD unavailable and VEP reports none)")

    # ClinVar
    clin = None
    others: List[Dict[str, Any]] = []
    if cv is not None:
        sources += cv.sources
        warnings += cv.warnings
        caid = (population or {}).get("caid")
        exact, others = match_clinvar(cv.result, vcf, assembly, window, caid)
        if exact:
            clin = clinvar.compact(exact[0])
            if len(exact) > 1:
                warnings.append(f"ClinVar: {len(exact)} records for this allele ({', '.join(r['vcv'] for r in exact)}); showing the best reviewed")
        elif vcf is None and len([r for r in cv.result if not r.get('compound')]) == 1:
            only = [r for r in cv.result if not r.get("compound")][0]
            clin = clinvar.compact(only)
            warnings.append(f"ClinVar {only['vcv']}: allele match not verified (no VCF form)")

    literature = None
    if lv is not None:
        try:
            sources += list(getattr(lv, "sources", []) or [])
            warnings += list(getattr(lv, "warnings", []) or [])
            literature = {"query": f"{gene_symbol} {lit_text}" if gene_symbol and lit_text == hgvs.get("p") else lit_text,
                          "litvar": _litvar_compact(getattr(lv, "result", lv), rsids)}
        except Exception as err:  # noqa: BLE001
            warnings.append(f"LitVar result unreadable: {err}")

    # ACMG inputs
    terms = consequence.get("terms") or []
    revel = preds["revel"]
    sp_max = (preds["spliceai"] or {}).get("max")
    insil = (population or {}).get("in_silico") or {}
    if revel is None and insil.get("revel_max") not in (None, ""):
        revel = _num(insil.get("revel_max"))
        if revel is not None:
            preds["revel_gnomad"] = revel
    if sp_max is None and insil.get("spliceai_ds_max") not in (None, ""):
        sp_max = _num(insil.get("spliceai_ds_max"))
        if sp_max is not None:
            preds["spliceai_gnomad_max"] = sp_max
    grp = (population or {}).get("grpmax") or {}
    if population is not None:
        if population.get("found"):
            gac = (population.get("total") or {}).get("ac")
        else:
            gac = 0 if population.get("covered") else None
            if population.get("covered") is not True:
                warnings.append("absent from gnomAD but site coverage is low or unknown: PM2 not supported by absence alone")
    else:
        gac = None
    acmg_inputs = {
        "revel": revel,
        "spliceai_max": sp_max,
        "grpmax_af": grp.get("af") if population else (fallback or {}).get("grpmax_af"),
        "grpmax_an": grp.get("an") if population else None,
        "gnomad_ac": gac,
        "faf95": ((population or {}).get("faf95") or {}).get("value"),
        "consequence": list(terms),
        "is_missense": any("missense" in t for t in terms),
        "frequency_source": (population or {}).get("dataset") or ("VEP colocated gnomAD" if fallback else None),
    }

    other_assembly = (population or {}).get("liftover")
    label_parts = [hgvs.get("c") or text]
    if hgvs.get("p"):
        label_parts.append(f"({hgvs['p']})")
    if gene_symbol:
        label_parts.append(gene_symbol)
    if vcf:
        label_parts.append(f"| {assembly} {vcf[0]}-{vcf[1]}-{vcf[2]}-{vcf[3]}")
    result = {
        "variant": " ".join(label_parts),
        "input": text,
        "kind": kind,
        "assembly": assembly,
        "vcf": {"assembly": assembly, "chrom": vcf[0], "pos": vcf[1], "ref": vcf[2], "alt": vcf[3],
                "id": f"{vcf[0]}-{vcf[1]}-{vcf[2]}-{vcf[3]}"} if vcf else None,
        "other_assembly": other_assembly,
        "rsids": rsids,
        "caid": (population or {}).get("caid"),
        "gene": gene_symbol,
        "transcript": transcript,
        "hgvs": hgvs,
        "consequence": consequence,
        "predictors": preds,
        "population": population if population is not None else ({"fallback": fallback} if fallback else None),
        "clinvar": clin,
        "clinvar_others_at_locus": [{"vcv": r["vcv"], "title": r["title"], "classification": r["classification"],
                                     "stars": r["stars"], "url": r["url"]} for r in others[:5]],
        "literature": literature,
        "acmg_inputs": acmg_inputs,
    }
    if cv is not None and clin is None:
        result["clinvar_note"] = "no ClinVar record for this exact allele" + (
            f" ({len(others)} other record(s) at the locus)" if others else "")
    return Outcome(result, sources=sources, warnings=warnings)
