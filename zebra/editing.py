"""Base-editing feasibility screen: could an ABE or CBE put the patient's ALT back to REF?

A screen, not a design. It answers three questions from the reference sequence alone:

1. Is the correction one that a canonical base editor can make at all? Adenine base
   editors (ABE) deaminate adenine on the protospacer (non-target) strand -- the strand
   the guide RNA displaces -- so the edit reads A->G on that strand; cytosine base
   editors (CBE) read C->T there. A correction is therefore reachable only when the
   patient's allele is A (ABE) or C (CBE) on one of the two strands, with the reference
   allele being G or T on that same strand. Everything else (transversions, indels) is
   out of reach for ABE/CBE.
2. Is there a protospacer that puts that base inside the editing window? Positions are
   counted 1-20 from the PAM-distal end of the 20 nt protospacer. The window is not the
   same for the two chemistries, so each has its own default: CBE 4-8 (BE3/BE4; Komor
   et al. 2016, Nature 533:420) and ABE 4-7 (ABE7.10; Gaudelli et al. 2017, Nature
   551:464). ABE8e is wider and shifted, about 3-9 (Richter et al. 2020, Nat Biotechnol
   38:883), and is reported alongside as an alternative. `--window` overrides both.
3. What else would be edited? Every base of the same kind inside the window is a
   bystander and will be edited too.

What it does not answer: whether a guide works in a cell, how efficient it is, what the
off-target profile looks like, and -- the part that decides whether any of this is
possible for a person -- how the editor would reach the tissue that matters.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

from zebra.core import Outcome, UsageError
from zebra.s2f import normalise_variant
from zebra.sources import ensembl

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
PROTOSPACER_LEN = 20
# PAM, immediately 3' of the protospacer on the same strand. NGG is SpCas9; NG covers the
# engineered relaxed-PAM variants (SpCas9-NG, SpRY and friends), at lower efficiency.
PAMS = {
    "NGG": {"length": 3, "fixed": {1: "G", 2: "G"}, "nuclease": "SpCas9 (canonical NGG)"},
    "NG": {"length": 2, "fixed": {1: "G"}, "nuclease": "SpCas9-NG / SpRY-type relaxed PAM (lower efficiency)"},
}
EDITORS = {
    ("A", "G"): {
        "editor": "ABE",
        "name": "adenine base editor",
        "chemistry": "A->G on the protospacer (non-target) strand, by adenine deamination to inosine; the duplex ends as A*T -> G*C",
    },
    ("C", "T"): {
        "editor": "CBE",
        "name": "cytosine base editor",
        "chemistry": "C->T on the protospacer (non-target) strand, by cytosine deamination to uracil; the duplex ends as C*G -> T*A",
    },
}
# Editing windows are chemistry- and generation-specific; one window for both editors
# overstates feasibility for ABE at positions 8 and understates the bystander count for
# ABE8e at positions 3 and 9. Each entry carries the editor generation it belongs to.
EDITOR_WINDOWS = {
    "CBE": {
        "window": (4, 8),
        "generation": "BE3 / BE4",
        "citation": "Komor et al. 2016, Nature 533:420 (BE3): cytosine deamination window at protospacer positions 4-8",
    },
    "ABE": {
        "window": (4, 7),
        "generation": "ABE7.10",
        "citation": "Gaudelli et al. 2017, Nature 551:464 (ABE7.10): adenine deamination window at protospacer positions 4-7",
    },
}
# Alternative generations worth naming, because a guide that is bystander-free in the
# default window may not be in the editor that actually gets used.
EDITOR_WINDOW_ALTERNATIVES = {
    "ABE": [{"generation": "ABE8e", "window": (3, 9),
             "citation": "Richter et al. 2020, Nat Biotechnol 38:883 (ABE8e): wider and shifted window, about "
                         "protospacer positions 3-9 — rerun with --window 3-9 to see the bystanders it adds"}],
    "CBE": [{"generation": "BE4max / evoAPOBEC1", "window": (4, 8),
             "citation": "Koblan et al. 2018, Nat Biotechnol 36:843: BE4max keeps the BE3 position 4-8 window; "
                         "narrowed-window variants (YE1, eA3A) exist and edit fewer bystanders"}],
}
# Used only where no editor applies, so the output still carries a window for display.
DEFAULT_WINDOW = EDITOR_WINDOWS["CBE"]["window"]
FLANK = 40

LABELS = [
    "Feasibility screen, not a guide design: no efficiency, no off-target search, no validation.",
    "Delivery is the hard part. An editor that works in a dish still has to reach the tissue that drives the "
    "disease, in enough cells, at the right time; for most tissues that is unsolved.",
    "Bystander edits inside the window happen; a silent bystander is fine, a missense or splice one is not. "
    "Check every bystander on the transcript.",
    "Sequence here is the reference plus the patient's allele, not the patient's own read data: a second variant "
    "in cis inside the protospacer or PAM would change the answer.",
    "The editing window is chemistry- and generation-specific: CBE 4-8 (BE3/BE4; Komor et al. 2016, Nature "
    "533:420), ABE 4-7 (ABE7.10; Gaudelli et al. 2017, Nature 551:464), ABE8e about 3-9 (Richter et al. 2020, "
    "Nat Biotechnol 38:883). A guide with no bystander in one window may have one in another.",
]
# The two failure modes most likely to invalidate a "correctable" verdict, neither of
# which a protospacer-based bystander or off-target search can see.
CBE_BYPRODUCT_LABEL = (
    "CBE products are not only C->T. Cytosine base editors also make C->G and C->A at the target and at "
    "bystander positions at appreciable rates, and uracil excision can leave indels at the edited site "
    "(Komor et al. 2017, Sci Adv 3:eaao4774). \"CBE can revert this\" does not mean one clean product: the "
    "edited allele has to be sequenced, not assumed."
)
GUIDE_INDEPENDENT_LABEL = (
    "Guide-independent deamination is not visible to any protospacer-based off-target search. CBEs cause "
    "genome-wide Cas9-independent C->T mutations (Zuo et al. 2019, Science 364:289; Jin et al. 2019, Science "
    "364:292) and cytosine deaminase editors cause transcriptome-wide RNA edits (Grunewald et al. 2019, Nature "
    "569:433); ABEs do the same at lower rates. These need orthogonal assays (whole-genome or RNA sequencing of "
    "edited cells), not a guide search."
)
PRIME_EDITING_NOTE = (
    "Prime editing (a Cas9 nickase fused to a reverse transcriptase, with a pegRNA template) is not restricted to "
    "the four transitions and can write insertions and deletions; it is the route to look at when ABE/CBE cannot "
    "make the change."
)
CGBE_NOTE = (
    "C-to-G base editors (CGBE) and adenine transversion editors are reported for some C*G->G*C and A*T->C*G "
    "changes, but they are early-stage and context-dependent; this screen does not model them."
)


def complement(base: str) -> str:
    try:
        return COMPLEMENT[base.upper()]
    except KeyError:
        raise UsageError(f"{base!r} is not a DNA base (A/C/G/T/N)") from None


def revcomp(seq: str) -> str:
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


def pam_matches(seq: str, start: int, pam: str) -> Optional[str]:
    """The PAM read from `seq` at `start`, or None if it does not fit the pattern."""
    spec = PAMS[pam]
    end = start + spec["length"]
    if start < 0 or end > len(seq):
        return None
    text = seq[start:end]
    if "N" in text:
        return None
    for offset, base in spec["fixed"].items():
        if text[offset] != base:
            return None
    return text


def _guides(strand_seq: str, target_index: int, pam: str, window: Tuple[int, int],
            from_base: str, to_base: str, genomic_of) -> List[Dict[str, Any]]:
    """Every protospacer on this strand that puts `target_index` inside the editing window."""
    lo, hi = window
    out: List[Dict[str, Any]] = []
    # protospacer position p (1..20, from the PAM-distal end) of a base at index t for a
    # protospacer starting at index s is t - s + 1, so s = t - p + 1.
    for position in range(lo, hi + 1):
        start = target_index - position + 1
        if start < 0 or start + PROTOSPACER_LEN > len(strand_seq):
            continue
        protospacer = strand_seq[start:start + PROTOSPACER_LEN]
        if "N" in protospacer:
            continue
        pam_text = pam_matches(strand_seq, start + PROTOSPACER_LEN, pam)
        if pam_text is None:
            continue
        bystanders = []
        for p in range(lo, hi + 1):
            idx = start + p - 1
            if idx == target_index or idx >= len(strand_seq):
                continue
            if strand_seq[idx] == from_base:
                bystanders.append({
                    "protospacer_position": p,
                    "position": genomic_of(idx),
                    "edit_on_protospacer_strand": f"{from_base}>{to_base}",
                    "coding_effect": "check bystanders on the transcript",
                })
        same_base_in_protospacer = sum(1 for i, b in enumerate(protospacer) if b == from_base and start + i != target_index)
        out.append({
            "protospacer": protospacer,
            "protospacer_5to3": f"5'-{protospacer}-{pam_text.lower()}-3'",
            "pam": pam_text,
            "pam_pattern": pam,
            "nuclease": PAMS[pam]["nuclease"],
            "target_protospacer_position": position,
            "target_position": genomic_of(target_index),
            "bystanders_in_window": bystanders,
            "bystander_count": len(bystanders),
            "same_base_elsewhere_in_protospacer": same_base_in_protospacer,
            "pam_distal_position": genomic_of(start),
        })
    return out


def _annotate_bystanders(rows: Sequence[Dict[str, Any]], chrom: str, assembly: str,
                         plus_seq: str, window_start: int) -> Tuple[Dict[int, Dict[str, Any]], List[Dict[str, Any]], List[str]]:
    """VEP for each distinct bystander edit, expressed on the forward strand."""
    edits: Dict[int, Tuple[str, str]] = {}
    for guide in rows:
        for by in guide["bystanders_in_window"]:
            pos = by["position"]
            if pos is None or pos in edits:
                continue
            idx = pos - window_start
            if not 0 <= idx < len(plus_seq):
                continue
            ref = plus_seq[idx]
            # the protospacer-strand edit A>G is G>A on the forward strand when the
            # protospacer is the minus strand; derive it from the forward base
            on_proto = by["edit_on_protospacer_strand"].split(">")
            alt = on_proto[1] if ref == on_proto[0] else complement(on_proto[1])
            if ref == alt:
                continue
            edits[pos] = (ref, alt)
    if not edits:
        return {}, [], []
    variants = [(chrom, pos, ref, alt) for pos, (ref, alt) in sorted(edits.items())]
    out = ensembl.vep_batch(variants, assembly=assembly)
    annotated: Dict[int, Dict[str, Any]] = {}
    for rec in out.result:
        pos = rec.get("start")
        tc = ensembl.pick_transcript(rec) or {}
        annotated[int(pos)] = {
            "most_severe_consequence": rec.get("most_severe_consequence"),
            "gene": tc.get("gene_symbol"),
            "transcript": tc.get("transcript_id"),
            "hgvsc": tc.get("hgvsc"),
            "hgvsp": tc.get("hgvsp"),
            "silent": rec.get("most_severe_consequence") in ("synonymous_variant", "intron_variant",
                                                             "upstream_gene_variant", "downstream_gene_variant",
                                                             "intergenic_variant", "3_prime_UTR_variant",
                                                             "5_prime_UTR_variant"),
        }
    return annotated, out.sources, out.warnings


def _editor_for(ref: str, alt: str) -> Optional[str]:
    """Which editor could put `alt` back to `ref`, on either strand (exactly one can)."""
    for spec in (EDITORS.get((alt, ref)), EDITORS.get((complement(alt), complement(ref)))):
        if spec:
            return spec["editor"]
    return None


def screen(variant: str, assembly: str = "GRCh38", pam: str = "NGG",
           window: Optional[Tuple[int, int]] = None, flank: int = FLANK,
           annotate_bystanders: bool = False) -> Outcome:
    """Can a base editor revert this variant, and with which protospacers?

    `window` defaults to the editing window of the chemistry the correction needs --
    CBE 4-8, ABE 4-7 (see EDITOR_WINDOWS for the citations) -- rather than one window
    for both. An explicit window overrides it.
    """
    if assembly not in ("GRCh38", "GRCh37"):
        raise UsageError("assembly must be GRCh38 or GRCh37")
    pam = (pam or "NGG").upper()
    if pam not in PAMS:
        raise UsageError(f"--pam must be one of {', '.join(PAMS)}")
    explicit_window = window is not None
    if explicit_window:
        lo, hi = int(window[0]), int(window[1])
        if not 1 <= lo <= hi <= PROTOSPACER_LEN:
            raise UsageError(f"--window must be LO-HI with 1 <= LO <= HI <= {PROTOSPACER_LEN} "
                             f"(defaults: CBE 4-8, ABE 4-7)")
    else:
        lo, hi = DEFAULT_WINDOW

    norm = normalise_variant(variant, assembly=assembly)
    v = norm.result
    sources = list(norm.sources)
    warnings = list(norm.warnings)
    chrom, pos, ref, alt = v["chrom"], v["pos"], v["ref"], v["alt"]

    # Pick the window from the chemistry the correction needs, before anything is fetched,
    # so the sequence flank and the bystander search both use the right one.
    window_basis = "given with --window; it overrides the per-editor defaults"
    if not explicit_window and len(ref) == 1 and len(alt) == 1:
        editor_guess = _editor_for(ref.upper(), alt.upper())
        if editor_guess:
            spec = EDITOR_WINDOWS[editor_guess]
            lo, hi = spec["window"]
            window_basis = f"{editor_guess} default ({spec['generation']}): {spec['citation']}"
    flank = max(flank, hi + PROTOSPACER_LEN + PAMS[pam]["length"] + 2)

    base: Dict[str, Any] = {
        "variant": v,
        "assembly": assembly,
        "pam": pam,
        "editing_window": [lo, hi],
        "editing_window_basis": window_basis,
        "editor_windows": {k: {"window": list(spec["window"]), "generation": spec["generation"],
                               "citation": spec["citation"]} for k, spec in EDITOR_WINDOWS.items()},
        "window_convention": "protospacer positions 1-20 counted from the PAM-distal end",
        "labels": LABELS,
    }

    if len(ref) != 1 or len(alt) != 1:
        base.update(
            base_editable=False,
            reason=(f"{ref}>{alt} is an insertion/deletion. Canonical base editors change one base pair in place "
                    f"(A*T<->G*C, C*G<->T*A) and cannot restore lost or extra bases."),
            editor=None, guides=[], guide_count=0,
            alternatives=[PRIME_EDITING_NOTE],
        )
        return Outcome(base, sources=sources, warnings=warnings, text=render(base),
                       query={"variant": variant, "assembly": assembly, "pam": pam, "window": [lo, hi]})

    ref, alt = ref.upper(), alt.upper()
    if ref == alt:
        raise UsageError(f"REF and ALT are both {ref} at {chrom}:{pos}: that is not a variant")

    start, end = max(1, pos - flank), pos + flank
    seq_out = ensembl.sequence(chrom, start, end, assembly=assembly)
    sources.extend(seq_out.sources)
    plus_ref = seq_out.result.upper()
    idx = pos - start
    if idx >= len(plus_ref):
        raise UsageError(f"Ensembl returned {len(plus_ref)} bases for {chrom}:{start}-{end}; cannot place {pos}")
    observed = plus_ref[idx]
    if observed != ref:
        raise UsageError(
            f"REF mismatch: {variant} says {ref} at {chrom}:{pos} but the {assembly} reference has {observed}. "
            f"Check the assembly and the strand (transcript HGVS on a minus-strand gene is not the forward allele)."
        )
    if len(plus_ref) < end - start + 1:
        warnings.append(
            f"Ensembl returned {len(plus_ref)} of the {end - start + 1} bases asked for around {chrom}:{pos} "
            f"(contig edge): protospacers that would run past the end are not listed"
        )
    plus_patient = plus_ref[:idx] + alt + plus_ref[idx + 1:]
    minus_patient = revcomp(plus_patient)
    minus_index = len(plus_patient) - 1 - idx

    base["reference_window"] = {
        "chrom": chrom, "start": start, "end": end,
        "reference_plus": plus_ref,
        "patient_plus": plus_patient,
        "variant_index": idx,
        "ref_confirmed": True,
    }
    base["correction"] = {"from": alt, "to": ref, "on": "forward strand",
                          "description": f"restore {chrom}:{pos} {alt} (patient) to {ref} (reference)"}

    strands: List[Dict[str, Any]] = []
    guides: List[Dict[str, Any]] = []
    for strand, seq, target_index in (("+", plus_patient, idx), ("-", minus_patient, minus_index)):
        from_base = alt if strand == "+" else complement(alt)
        to_base = ref if strand == "+" else complement(ref)
        spec = EDITORS.get((from_base, to_base))
        row: Dict[str, Any] = {
            "strand": strand,
            "base_on_protospacer_strand": f"{from_base}>{to_base}",
            "editor": spec["editor"] if spec else None,
        }
        if not spec:
            row["reason"] = (f"on the {strand} strand the correction reads {from_base}>{to_base}, which is neither "
                             f"A>G (ABE) nor C>T (CBE)")
            row["guides"] = []
            row["guide_count"] = 0
            strands.append(row)
            continue
        if strand == "+":
            def genomic_of(i, _s=start):
                return _s + i
        else:
            def genomic_of(i, _s=start, _l=len(plus_patient)):
                return _s + (_l - 1 - i)
        found = _guides(seq, target_index, pam, (lo, hi), from_base, to_base, genomic_of)
        for g in found:
            g.update(strand=strand, editor=spec["editor"], editor_name=spec["name"], chemistry=spec["chemistry"],
                     edit_on_protospacer_strand=f"{from_base}>{to_base}")
        row.update(editor_name=spec["name"], chemistry=spec["chemistry"], guides=found, guide_count=len(found))
        strands.append(row)
        guides.extend(found)

    editors = sorted({r["editor"] for r in strands if r["editor"]})
    guides.sort(key=lambda g: (g["bystander_count"], abs(g["target_protospacer_position"] - (lo + hi) // 2),
                               g["strand"], g["target_protospacer_position"]))

    if annotate_bystanders and guides:
        annotated, src, warn = _annotate_bystanders(guides, chrom, assembly, plus_patient, start)
        sources.extend(src)
        warnings.extend(warn)
        for g in guides:
            for by in g["bystanders_in_window"]:
                hit = annotated.get(by["position"])
                if hit:
                    by["coding_effect"] = hit["most_severe_consequence"]
                    by["vep"] = hit

    # Chemistry-specific caveats: a CBE verdict implies products other than C->T, and
    # both chemistries deaminate without a guide, which no protospacer search can see.
    labels = list(LABELS)
    alternatives_note: List[Dict[str, Any]] = []
    if "CBE" in editors:
        labels.append(CBE_BYPRODUCT_LABEL)
    if editors:
        labels.append(GUIDE_INDEPENDENT_LABEL)
        for ed in editors:
            for alt_spec in EDITOR_WINDOW_ALTERNATIVES.get(ed, []):
                alternatives_note.append({"editor": ed, **alt_spec, "window": list(alt_spec["window"])})
    base.update(
        base_editable=bool(editors),
        editors=editors,
        editor=editors[0] if editors else None,
        strands=strands,
        guides=guides[:20],
        guide_count=len(guides),
        bystander_free_guides=sum(1 for g in guides if g["bystander_count"] == 0),
        bystanders_annotated=bool(annotate_bystanders),
        labels=labels,
        other_editor_generations=alternatives_note,
    )
    if alternatives_note and not explicit_window:
        base["window_note"] = (
            "bystander counts are for the window above. Other generations of the same chemistry have different "
            "windows: " + "; ".join(f"{a['generation']} {a['window'][0]}-{a['window'][1]}" for a in alternatives_note)
            + ". A guide listed with 0 bystanders here is not bystander-free for all of them."
        )
    if not editors:
        base["reason"] = (f"{alt}>{ref} is a transversion: it is A*T<->C*G or C*G<->G*C, and canonical ABE/CBE only "
                          f"make the transitions A*T<->G*C and C*G<->T*A")
        base["alternatives"] = [CGBE_NOTE, PRIME_EDITING_NOTE]
    elif not guides:
        base["reason"] = (f"an {editors[0]} could make this change, but no {pam} protospacer places the base in "
                          f"window {lo}-{hi}; widen --window, try --pam NG, or consider prime editing")
        base["alternatives"] = [PRIME_EDITING_NOTE]
    if not annotate_bystanders and any(g["bystander_count"] for g in guides):
        base["bystander_note"] = ("bystander coding effects were not annotated; rerun with --annotate-bystanders "
                                  "to put each one through VEP")
    return Outcome(base, sources=sources, warnings=warnings, text=render(base),
                   query={"variant": variant, "assembly": assembly, "pam": pam, "window": [lo, hi],
                          "annotate_bystanders": annotate_bystanders})


def render(result: Dict[str, Any]) -> str:
    v = result["variant"]
    lines = [f"{v['input']} -> {v['chrom']}-{v['pos']}-{v['ref']}-{v['alt']} ({result['assembly']}; {v['resolved_from']})"]
    if not result.get("base_editable"):
        lines.append(f"not base-editable with ABE/CBE: {result.get('reason')}")
        for alt in result.get("alternatives") or []:
            lines.append(f"  - {alt}")
        lines += ["", *(f"! {label}" for label in result["labels"])]
        return "\n".join(lines)
    corr = result["correction"]
    lo, hi = result["editing_window"]
    lines.append(f"correct {corr['from']} -> {corr['to']}: {', '.join(result['editors'])}; "
                 f"PAM {result['pam']}; window {lo}-{hi} ({result['window_convention']})")
    if result.get("editing_window_basis"):
        lines.append(f"  window: {result['editing_window_basis']}")
    for row in result["strands"]:
        if row["editor"]:
            lines.append(f"  {row['strand']} strand: {row['editor']} ({row['base_on_protospacer_strand']} on the "
                         f"protospacer strand) -> {row['guide_count']} protospacer(s)")
        else:
            lines.append(f"  {row['strand']} strand: {row['reason']}")
    lines.append(f"{result['guide_count']} protospacer(s), {result['bystander_free_guides']} with no bystander in the window")
    if result.get("reason"):
        lines.append(f"note: {result['reason']}")
    for g in result["guides"]:
        by = ", ".join(f"p{b['protospacer_position']}@{b['position']}"
                       + (f" [{b['coding_effect']}]" if b.get("coding_effect") and b["coding_effect"] != "check bystanders on the transcript" else "")
                       for b in g["bystanders_in_window"]) or "none in window"
        lines.append(f"  {g['editor']} {g['strand']} pos{g['target_protospacer_position']:>3}  {g['protospacer_5to3']}"
                     f"  target {g['target_position']}  bystanders: {by}")
    if result.get("bystander_note"):
        lines.append(f"note: {result['bystander_note']}")
    if result.get("window_note"):
        lines.append(f"note: {result['window_note']}")
    lines += ["", *(f"! {label}" for label in result["labels"])]
    return "\n".join(lines)
