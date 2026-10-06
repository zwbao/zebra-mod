"""Splice-switching antisense feasibility screen: is there anything on the pre-mRNA to aim at?

A research feasibility screen, not a design. It answers four questions, in order, and
stops at the first one that has no answer:

1. Does a splicing model predict an aberrant event at all? The screen reads the
   SpliceAI result for the variant (`zebra.s2f.predict`) and takes the gained acceptor
   and gained donor positions it reports. A gained acceptor AND a gained donor a
   plausible exon apart reconstruct a cryptic exon (pseudoexon); one of them alone is a
   cryptic acceptor or cryptic donor. Only losses, or nothing above the floor, means
   there is no aberrant site to block -- the screen says so and stops, because an
   antisense oligonucleotide blocks splice sites, it cannot create one that a variant
   destroyed.
2. Is the predicted site real in the patient's own sequence? Each site is checked
   against the Ensembl reference with the patient allele applied: a splice acceptor
   must have AG immediately 5' of the first exonic base and a donor GT immediately 3'
   of the last one, read on the transcript strand. A site that is canonical on the
   patient allele but not on the reference is the site the variant created, which is
   the most specific thing to aim at.
3. Which windows on the pre-mRNA cover those sites? For every length asked for, the
   screen enumerates the windows that fully cover the site (its dinucleotide included)
   and, inside a cryptic exon, tiles the exon body as well. For each window it reports
   the genomic coordinates of the target, the target as the pre-mRNA reads it, the
   reverse complement of that target, GC fraction, a self-complementarity (hairpin)
   flag, the longest homopolymer run, and whether the patient's variant lies inside it.
4. Is the target unique in the genome? With `--uniqueness` each candidate goes through
   one NCBI BLAST search (see `zebra.sources.ncbi_blast`). Without it, every candidate
   says uniqueness was not checked and the ranking leaves that component out: a
   sequence is never presented as unique because nobody looked.

What it does NOT answer, and what decides whether any of this could ever matter:
chemistry, backbone, modification pattern, length optimisation, delivery to the tissue,
dose, schedule, toxicity, immunogenicity, and whether the aberrant splicing is actually
happening in the patient. Nothing here is a drug, a dose or a recommendation.

The precondition the whole screen rests on: RNA evidence. A SpliceAI score is a
hypothesis about splicing. The aberrant transcript has to be shown in RNA from a tissue
that expresses the gene (`zebra expression <gene>` says which tissue that could be), or
in a minigene assay, before any of these coordinates mean anything.

Precedent and eligibility records come from Europe PMC, retrieved during the run and
listed as retrieved. Nothing about who is eligible for an individualized programme is
stated here from memory.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

from zebra.core import Outcome, UsageError
from zebra.sources import attempt, ensembl, europepmc

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

EVENTS = ("auto", "pseudoexon", "cryptic_acceptor", "cryptic_donor", "exon_skip")
DEFAULT_LENGTHS = (18, 25)
LENGTH_LIMITS = (12, 40)
# ClinGen SVI / Walker et al. 2023, AJHG 110:1046: a SpliceAI delta >= 0.2 supports a
# splice effect. Below 0.1 supports no effect. Between them is uninformative -- and that
# band is where some proven splice-disrupting variants sit, which is why the screen
# reports a site there instead of calling it absent (see DELTA_BANDS).
MIN_DELTA = 0.2
DELTA_FLOOR = 0.1
DELTA_BANDS = (
    "SpliceAI delta >= 0.2 supports a splice effect (PP3 at SUPPORTING weight only; Walker et al. 2023, AJHG "
    "110:1046, ClinGen SVI Splicing Subgroup).",
    "0.1 < delta < 0.2 is UNINFORMATIVE: it supports neither a splice effect nor its absence. A site in this band "
    "is reported here, flagged, and is NOT evidence that the event happens -- CFTR c.3718-2477C>T (3849+10kbC>T), "
    "a cystic-fibrosis-causing variant whose cryptic exon is established in RNA, scores 0.162 at this setting.",
    "delta <= 0.1 supports no splice effect; the screen stops there.",
)
MIN_PSEUDOEXON = 20
MAX_PSEUDOEXON = 1000
GC_BAND = (0.40, 0.60)
GC_HARD = (0.20, 0.80)
HAIRPIN_STEM = 5
HAIRPIN_MIN_LOOP = 3
HOMOPOLYMER_MAX = 3
DEFAULT_TOP = 12
SITE_STRIDE = 2
BODY_STRIDE = 10
# A ceiling on how many windows are built, so `--lengths 12-40 --body-stride 1` over a
# 1 kb cryptic exon cannot quietly turn into ~30,000 candidates in memory (measured:
# 26,738 windows, 1.8 s, 152 MB peak at a ceiling of 20,000; 8,000 keeps it near 60 MB).
# Only `--top` candidates are ever shown, so a lower ceiling costs nothing in coverage
# once the thinning is even. The number of windows cut is always reported.
MAX_WINDOWS = 8000
CANONICAL_ACCEPTOR = ("AG",)
# GT is the major class; GC-AG introns are a real minor class (~0.8% of human introns),
# so a GC donor is reported as canonical-minor rather than as "not a splice site".
CANONICAL_DONOR = ("GT", "GC")
MAJOR_DONOR = "GT"

REGION_ORDER_NOTE = (
    "Groups are ordered by how specific the target is to this variant, not by score: the site the patient's allele "
    "CREATED (canonical on the patient allele, not on the reference) comes first, then the other aberrant site, "
    "then the body of the cryptic exon. Within a group, candidates are ordered by the component scores shown."
)
RANK_COMPONENTS = {
    "site_centering": (2.0, "how close the splice site sits to the middle of the window; a site at the very edge "
                            "is covered in name only. 1.0 = centred, 0.0 = at the edge. Not defined for a window "
                            "inside the exon body, which covers no site."),
    "gc_in_band": (1.0, f"GC fraction inside {GC_BAND[0]:.0%}-{GC_BAND[1]:.0%} scores 1.0, falling linearly to 0 "
                        f"at {GC_HARD[0]:.0%} / {GC_HARD[1]:.0%}. A composition rule of thumb, not a melting "
                        "temperature."),
    "no_hairpin": (1.0, f"1.0 when the longest self-complementary stem is shorter than {HAIRPIN_STEM} bases "
                        f"(loop >= {HAIRPIN_MIN_LOOP}), else 0.0. Counted from the sequence, not from a "
                        "thermodynamic model."),
    "no_homopolymer": (0.5, f"1.0 when the longest single-base run is at most {HOMOPOLYMER_MAX}, else 0.0."),
    "uniqueness": (2.0, "1.0 when the target matches one GRCh38 primary-chromosome locus, 0.0 when it matches "
                        "more. Left out of the score entirely when uniqueness was not checked, and the maximum "
                        "is reduced to match."),
    "variant_in_target": (0.0, "WEIGHT 0, reported and not scored: a window covering the patient's own variant "
                               "would discriminate between the two alleles, which matters for a dominant or "
                               "gain-of-function mechanism and not for a recessive one. This screen does not "
                               "decide the mechanism, so it does not score the flag."),
}

PRECONDITIONS = (
    "RNA evidence of the aberrant splicing, in a tissue that expresses the gene, is a precondition for everything "
    "below -- not a confirmation step afterwards. A SpliceAI or Pangolin score is a hypothesis about splicing; the "
    "event itself has to be shown by RT-PCR or RNA-seq on patient RNA, or in a minigene assay. `zebra expression "
    "<gene>` says which tissue that could be. If the aberrant transcript is degraded by nonsense-mediated decay, "
    "the RNA test needs NMD inhibition or it can miss the event.",
    "Mechanism fit: an antisense oligonucleotide that binds a splice site BLOCKS it. That helps when the problem "
    "is a site being USED that should not be (a cryptic site, a cryptic exon, an exon that has to be skipped to "
    "restore the frame). It cannot restore a site the variant destroyed, and it cannot put back sequence.",
    "A reachable tissue: the tissue that drives the disease has to be one an oligonucleotide can reach at all. "
    "That is a question about the disease and the route, not about these sequences.",
    "Allele and phase: the screen applies the patient's single allele to the reference. A second variant in cis "
    "inside a target window, or the wrong assembly, changes every coordinate here.",
)
LIMITS = (
    "Feasibility screen, not a design: no chemistry, no backbone, no modification pattern, no length optimisation, "
    "no delivery route, no dose, no schedule, no toxicity or immunogenicity assessment. None of those are "
    "computable from sequence, and all of them decide whether anything here could matter.",
    "Sequence is the Ensembl reference with the patient's allele applied, not the patient's own reads.",
    "Uniqueness, where it was checked, is checked against the genome and not the transcriptome: a target unique in "
    "the genome can still be partly complementary to other transcripts, which is where off-target RNA effects "
    "come from.",
    "No splicing-regulatory-element prediction is done. Windows inside a cryptic exon are a tiling of its body; "
    "calling one an exonic splicing enhancer would require an ESE model, and none is used here.",
    "GC and hairpin flags are composition rules with their thresholds stated, not thermodynamics.",
    "The windows come from a predicted site. If the RNA shows a different boundary, these coordinates are wrong "
    "and the screen has to be rerun against the observed event.",
)
PRECEDENT_RULE = (
    "The records below were retrieved from Europe PMC during this run and are listed as retrieved, with their ids. "
    "Nothing about precedent or about who is eligible for an individualized antisense programme is stated here "
    "from memory; where a search returned nothing, the section is empty and says so. Read the records."
)
PRECEDENT_QUERIES = (
    {"key": "individualized_aso_precedent",
     "query": 'TITLE:"Patient-Customized Oligonucleotide Therapy for a Rare Genetic Disease"',
     "why": "the single-patient precedent for an individualized splice-switching antisense oligonucleotide "
            "(milasen, in CLN7/MFSD8 Batten disease)"},
    {"key": "n_lorem_eligibility",
     "query": '"n-Lorem" AND (eligibility OR eligible OR criteria OR "patient selection")',
     "why": "n-Lorem Foundation, which develops individualized antisense oligonucleotides free of charge; its own "
            "reports state whom it can take on"},
    {"key": "n1_collaborative",
     "query": '"N=1 Collaborative" AND (antisense OR oligonucleotide OR "nucleic acid")',
     "why": "the N=1 Collaborative, which publishes shared guidance for individualized nucleic acid therapies"},
)
PRECEDENT_LIMIT = 3
ABSTRACT_EXCERPT = 400


# ------------------------------------------------------------------ sequence primitives


def revcomp(seq: str) -> str:
    return "".join(COMPLEMENT.get(b, "N") for b in reversed((seq or "").upper()))


def gc_fraction(seq: str) -> Optional[float]:
    seq = (seq or "").upper()
    usable = [b for b in seq if b in "ACGT"]
    if not usable:
        return None
    return sum(1 for b in usable if b in "GC") / len(usable)


def hairpin(seq: str, min_loop: int = HAIRPIN_MIN_LOOP) -> Dict[str, Any]:
    """The longest self-complementary stem the sequence could fold back on.

    A stem is a pair of segments, `min_loop` or more bases apart, where one is the
    reverse complement of the other. This counts base pairs; it is not a free-energy
    calculation and says nothing about whether such a hairpin is stable.
    """
    s = (seq or "").upper()
    n = len(s)
    best = {"stem": 0, "loop": None, "five_prime_at": None, "three_prime_at": None, "stem_seq": None}
    for i in range(n):
        for j in range(i + min_loop, n):
            k = 0
            # s[i:i+k] pairs with s[j-k+1:j+1] reversed -> compare base by base
            while (i + k < j - k) and (COMPLEMENT.get(s[i + k]) == s[j - k]):
                k += 1
            loop = (j - k) - (i + k) + 1
            if k > best["stem"] and loop >= min_loop:
                best = {"stem": k, "loop": loop, "five_prime_at": i + 1, "three_prime_at": j + 1,
                        "stem_seq": s[i:i + k]}
    return best


def longest_run(seq: str) -> Dict[str, Any]:
    s = (seq or "").upper()
    best_base, best_len = None, 0
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j - i + 1 > best_len:
            best_base, best_len = s[i], j - i + 1
        i = j + 1
    return {"base": best_base, "length": best_len}


# ------------------------------------------------------------------ region bookkeeping


class Region:
    """A genomic window, plus the pre-mRNA (transcript-strand) reading of it.

    The pre-mRNA is what an antisense oligonucleotide pairs with, so every index in
    this screen is a pre-mRNA index and every reported coordinate is converted back
    to the forward genomic strand here.
    """

    def __init__(self, chrom: str, start: int, end: int, strand: int, forward: str):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = 1 if int(strand) >= 0 else -1
        self.forward = (forward or "").upper()
        self.premrna = self.forward if self.strand == 1 else revcomp(self.forward)
        # set by fetch_region: whether the patient's allele is actually in `forward`
        self.allele_applied = "not stated"

    def __len__(self) -> int:
        return len(self.premrna)

    def index(self, genomic_pos: int) -> int:
        """Pre-mRNA index of a forward-strand genomic position (can be out of range)."""
        return (genomic_pos - self.start) if self.strand == 1 else (self.end - genomic_pos)

    def genomic(self, index: int) -> int:
        return (self.start + index) if self.strand == 1 else (self.end - index)

    def holds(self, index: int) -> bool:
        return 0 <= index < len(self.premrna)

    def span(self, i0: int, i1: int) -> Tuple[int, int]:
        """Forward-strand (start, end) of the pre-mRNA slice [i0, i1]."""
        a, b = self.genomic(i0), self.genomic(i1)
        return (min(a, b), max(a, b))

    def slice(self, i0: int, i1: int) -> str:
        return self.premrna[i0:i1 + 1]


def fetch_region(chrom: str, lo: int, hi: int, strand: int, assembly: str,
                 variant_pos: int, ref: str, alt: str
                 ) -> Tuple[Region, Region, List[Dict[str, Any]], List[str]]:
    """The [lo, hi] window twice: with the patient's allele applied, and as the reference.

    Both come from one Ensembl request. The reference copy is what makes the
    "the variant CREATED this splice site" test possible: a dinucleotide that is
    canonical on the patient allele and not on the reference is the new site.
    """
    lo = max(1, int(lo))
    hi = int(hi)
    out = ensembl.sequence(chrom, lo, hi, assembly=assembly)
    reference = out.result.upper()
    warnings: List[str] = []
    expected = hi - lo + 1
    if len(reference) < expected:
        warnings.append(f"Ensembl returned {len(reference)} of the {expected} bases asked for at "
                        f"{chrom}:{lo}-{hi} (contig edge): windows that would run past the end are not listed")
        hi = lo + len(reference) - 1
    elif len(reference) > expected:
        raise UsageError(f"Ensembl returned {len(reference)} bases for {chrom}:{lo}-{hi}, which asked for "
                         f"{expected}: every coordinate derived from it would be shifted, so nothing is reported")
    forward = reference
    applied = "no: the variant lies outside the window fetched"
    idx = variant_pos - lo
    if 0 <= idx < len(reference):
        observed = reference[idx:idx + len(ref)]
        if len(observed) == len(ref) and observed != ref.upper():
            raise UsageError(
                f"REF mismatch: the variant says {ref} at {chrom}:{variant_pos} but the {assembly} reference has "
                f"{observed}. Check the assembly and the strand (transcript HGVS on a minus-strand gene is not the "
                f"forward allele)."
            )
        if len(observed) != len(ref):
            warnings.append(f"the {len(ref)}-base REF of {ref}>{alt} runs past the end of the sequence fetched at "
                            f"{chrom}:{lo}-{hi}, so it could not be checked against the reference")
        if len(ref) == 1 and len(alt) == 1 and len(observed) == 1:
            forward = reference[:idx] + alt.upper() + reference[idx + 1:]
            applied = f"yes: {ref}>{alt} at {chrom}:{variant_pos}"
        else:
            applied = (f"no: {ref}>{alt} is not a single-base substitution, so the sequence below is the plain "
                       f"reference")
            warnings.append(f"{ref}>{alt} is not a single-base substitution, so no allele was applied to the "
                            f"sequence: a window overlapping the variant shows reference bases, not the "
                            f"patient's, and the 'created by the variant' test cannot run. Read every coordinate "
                            f"here as approximate for this variant.")
    else:
        warnings.append(f"the variant at {chrom}:{variant_pos} lies outside the window fetched "
                        f"({chrom}:{lo}-{hi}); the sequence below is the plain reference")
    patient = Region(chrom, lo, hi, strand, forward)
    patient.allele_applied = applied
    return (patient, Region(chrom, lo, hi, strand, reference), list(out.sources), warnings)


# ------------------------------------------------------------------ the predicted event


def _strand_of(value: Any) -> Optional[int]:
    """+1 / -1, or None when the value is not a strand this screen recognises.

    Everything downstream -- the pre-mRNA, both dinucleotide checks, every reported
    coordinate -- depends on the strand, so guessing "+" for an absent or unreadable
    value would produce a complete, confident, wrong answer.
    """
    text = str(value if value is not None else "").strip()
    if text in ("-", "-1"):
        return -1
    if text in ("+", "1", "+1"):
        return 1
    return None


def _spliceai_row(models: Sequence[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    for row in models:
        if row.get("model") == "spliceai":
            return row
    return None


def _delta_value(delta: Dict[str, Any]) -> float:
    """The delta as a finite float; anything else is 0.0, which no threshold passes."""
    try:
        value = float(delta.get("delta"))
    except (TypeError, ValueError):
        return 0.0
    return value if value == value and abs(value) != float("inf") else 0.0


def _gains(transcript: Dict[str, Any]) -> Tuple[Dict[str, Dict[str, Any]], List[Dict[str, Any]]]:
    """(usable gains by kind, gains whose position cannot be used).

    A gain with no integer position cannot be turned into a target window, but it is
    still a gain: dropping it silently is how a 0.97 donor gain became "no aberrant
    splicing". `zebra.s2f` sets `position: None` whenever SpliceAI's offset is not an
    integer, so this is a state the upstream really produces.
    """
    out: Dict[str, Dict[str, Any]] = {}
    unusable: List[Dict[str, Any]] = []
    deltas = transcript.get("deltas")
    for delta in deltas if isinstance(deltas, list) else []:
        if not isinstance(delta, dict):
            continue
        kind = {"DS_AG": "acceptor", "DS_DG": "donor"}.get(delta.get("score"))
        if not kind:
            continue
        if isinstance(delta.get("position"), int):
            out[kind] = delta
        else:
            unusable.append({"kind": kind, "score": delta.get("score"), "delta": _delta_value(delta),
                             "position": delta.get("position"), "offset": delta.get("offset")})
    return out, unusable


def _losses(transcript: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    deltas = transcript.get("deltas")
    for delta in deltas if isinstance(deltas, list) else []:
        if not isinstance(delta, dict):
            continue
        kind = {"DS_AL": "acceptor", "DS_DL": "donor"}.get(delta.get("score"))
        if kind:
            out[kind] = dict(delta, delta=_delta_value(delta))
    return out


def _dinucleotide(region: Region, site_index: int, kind: str) -> Optional[str]:
    """The two intronic bases that make the site canonical, read on the transcript strand.

    Acceptor: the two bases immediately 5' of the first exonic base (AG).
    Donor: the two bases immediately 3' of the last exonic base (GT, or GC for the
    minor class).
    """
    if kind == "acceptor":
        i0, i1 = site_index - 2, site_index - 1
    else:
        i0, i1 = site_index + 1, site_index + 2
    if not (region.holds(i0) and region.holds(i1)):
        return None
    return region.premrna[i0:i1 + 1]


def _site_record(region: Region, reference_region: Region, kind: str, delta: Dict[str, Any],
                 min_delta: float, allele_applied: bool = True) -> Dict[str, Any]:
    position = int(delta["position"])
    index = region.index(position)
    patient = _dinucleotide(region, index, kind) if region.holds(index) else None
    ref_dinuc = _dinucleotide(reference_region, reference_region.index(position), kind) \
        if reference_region.holds(reference_region.index(position)) else None
    canonical = CANONICAL_ACCEPTOR if kind == "acceptor" else CANONICAL_DONOR
    is_canonical = patient in canonical if patient else False
    ref_canonical = ref_dinuc in canonical if ref_dinuc else False
    changed = bool(patient and ref_dinuc and patient != ref_dinuc)
    # "the variant created this site" has two ways of being true, and both are stated:
    # the reference dinucleotide was not a splice-site dinucleotide at all, or -- for a
    # donor -- it was the minor-class GC and the variant made it the major-class GT.
    # CFTR c.3718-2477C>T is the second kind: GC on the reference, GT on the patient
    # allele, and it is the donor the cryptic exon is reported to use.
    made_canonical = bool(allele_applied and is_canonical and changed and not ref_canonical)
    upgraded = bool(allele_applied and is_canonical and changed and kind == "donor"
                    and patient == MAJOR_DONOR and ref_dinuc != MAJOR_DONOR)
    if patient is None or ref_dinuc is None:
        basis = ("not determined: the two bases that would make this site canonical lie outside the sequence "
                 "that was fetched, so neither the reference nor the patient dinucleotide could be read")
    elif not allele_applied:
        basis = ("not determined: the patient's allele was not applied to the sequence (it is not a single-base "
                 "substitution), so the reference and 'patient' dinucleotides here are the same bases")
    elif made_canonical:
        basis = (f"the reference has {ref_dinuc} here, which is not a splice-site dinucleotide; the patient's "
                 f"allele makes it {patient}")
    elif upgraded:
        basis = (f"the reference has the minor-class {ref_dinuc} here; the patient's allele makes it the "
                 f"major-class {patient}")
    elif changed:
        basis = f"the variant changes this dinucleotide from {ref_dinuc} to {patient}, but not into a stronger class"
    else:
        basis = (f"the dinucleotide is {patient} on both the reference and the patient allele: the variant did not "
                 f"create this site, it is a site that the variant may cause to be USED")
    value = float(delta.get("delta") or 0.0)
    return {
        "kind": kind,
        "position": position,
        "premrna_index": index,
        "spliceai_score": delta.get("score"),
        "delta": value,
        "delta_band": ("supports a splice effect" if value >= min_delta else
                       "uninformative (0.1-0.2): neither supports nor excludes" if value > DELTA_FLOOR else
                       "supports no splice effect"),
        "passed_min_delta": value >= min_delta,
        "ref_probability": delta.get("ref_prob"),
        "alt_probability": delta.get("alt_prob"),
        "dinucleotide_patient_allele": patient,
        "dinucleotide_reference": ref_dinuc,
        "dinucleotide_expected": "AG immediately 5' of the first exonic base" if kind == "acceptor"
                                 else "GT (major class) or GC (minor class) immediately 3' of the last exonic base",
        "canonical_on_patient_allele": is_canonical,
        "canonical_on_reference": ref_canonical,
        "minor_class_donor": bool(kind == "donor" and patient == "GC"),
        "dinucleotide_changed_by_variant": changed,
        "created_by_variant": bool(made_canonical or upgraded),
        "created_by_variant_basis": basis,
    }


def detect_event(region: Region, reference_region: Region, transcript: Dict[str, Any], event: str,
                 min_delta: float) -> Dict[str, Any]:
    """Classify the aberrant event from the SpliceAI gains, checked against the sequence."""
    gains, unusable = _gains(transcript)
    allele_applied = str(getattr(region, "allele_applied", "")).startswith("yes")
    sites: Dict[str, Dict[str, Any]] = {}
    for kind, delta in gains.items():
        if _delta_value(delta) > DELTA_FLOOR:
            sites[kind] = _site_record(region, reference_region, kind, delta, min_delta,
                                       allele_applied=allele_applied)
    # the strongest gain is over EVERY gain the model reported, not only the ones above
    # the floor: "strongest gain 0.000" for a variant whose top gain is 0.09 is a lie
    strongest = max([_delta_value(d) for d in gains.values()] + [u["delta"] for u in unusable] + [0.0])
    out: Dict[str, Any] = {
        "requested": event,
        "gains_considered": sorted(sites.values(), key=lambda s: -s["delta"]),
        "gains_below_floor": sorted(
            [{"kind": k, "score": d.get("score"), "delta": _delta_value(d), "position": d.get("position")}
             for k, d in gains.items() if _delta_value(d) <= DELTA_FLOOR], key=lambda g: -g["delta"]),
        "gains_without_a_usable_position": unusable,
        "losses": [dict(d, kind=k) for k, d in sorted(_losses(transcript).items())],
        "min_delta": min_delta,
        "delta_floor": DELTA_FLOOR,
        "delta_bands": list(DELTA_BANDS),
        "strongest_gain": strongest,
        "patient_allele_applied": getattr(region, "allele_applied", "not stated"),
    }
    # A gain the model reported but whose position cannot be read is not an absence of a
    # gain. Refusing here is the only honest answer: there is nothing to aim at AND
    # something may well be happening.
    blocking = [u for u in unusable if u["delta"] > DELTA_FLOOR]
    # exon skipping does not rest on the predicted gains at all, so an unusable gain
    # position is not a reason to refuse it
    if blocking and not sites and event != "exon_skip":
        out.update(type=None, stop=True,
                   reason="SpliceAI reports " + "; ".join(
                       "a gained %s of %.3f" % (u["kind"], u["delta"]) for u in blocking)
                   + " for this variant, but with no usable position (offset %s), so no target window can be "
                     "placed. This is a prediction that could not be turned into coordinates, NOT an absence of "
                     "aberrant splicing: read the SpliceAI output directly."
                     % ", ".join(str(u.get("offset")) for u in blocking))
        return out
    if blocking:
        out["unusable_position_warning"] = (
            "a further gained site (%s) was reported above the floor but has no usable position, so no window "
            "covers it" % "; ".join("%s %.3f" % (u["kind"], u["delta"]) for u in blocking))
    acceptor, donor = sites.get("acceptor"), sites.get("donor")
    pseudoexon = None
    if acceptor and donor:
        size = (donor["premrna_index"] - acceptor["premrna_index"]) + 1
        if MIN_PSEUDOEXON <= size <= MAX_PSEUDOEXON:
            pseudoexon = {"acceptor": acceptor, "donor": donor, "size_bp": size,
                          "genomic": list(region.span(acceptor["premrna_index"], donor["premrna_index"])),
                          "in_frame": size % 3 == 0,
                          "frame_note": ("a length that is a multiple of 3 would be inserted in frame; the "
                                         "consequence still depends on whether it carries a stop codon"
                                         if size % 3 == 0 else
                                         "not a multiple of 3: inclusion shifts the frame downstream, which "
                                         "usually means a premature stop and nonsense-mediated decay")}
        elif size <= 0:
            out["pseudoexon_rejected"] = (
                f"the gained donor at {donor['position']} lies upstream of the gained acceptor at "
                f"{acceptor['position']} in transcript order, so the two cannot bound one exon; they are reported "
                f"as separate sites instead")
        else:
            out["pseudoexon_rejected"] = (
                f"a gained acceptor at {acceptor['position']} and a gained donor at {donor['position']} would make "
                f"an exon of {size} bp, outside the {MIN_PSEUDOEXON}-{MAX_PSEUDOEXON} bp this screen treats as a "
                f"cryptic exon; they are reported as separate sites instead")

    if event == "exon_skip":
        out.update(type="exon_skip", confidence="requested", sites=[],
                   why="exon skipping was asked for: the targets are the authentic splice sites of the exon, "
                       "taken from the transcript structure, not from the predicted gains")
        return out
    if event == "pseudoexon":
        if not pseudoexon:
            out.update(type=None, stop=True,
                       reason="a cryptic exon was asked for, but the prediction does not give a gained acceptor "
                              "AND a gained donor a plausible exon apart" +
                              (f" ({out['pseudoexon_rejected']})" if out.get("pseudoexon_rejected") else ""))
            return out
        out.update(type="pseudoexon", pseudoexon=pseudoexon, sites=[pseudoexon["acceptor"], pseudoexon["donor"]])
    elif event in ("cryptic_acceptor", "cryptic_donor"):
        kind = "acceptor" if event == "cryptic_acceptor" else "donor"
        site = sites.get(kind)
        if not site:
            out.update(type=None, stop=True,
                       reason=f"a gained {kind} was asked for, but SpliceAI reports no {kind} gain above "
                              f"{DELTA_FLOOR} for this variant on {transcript.get('refseq') or transcript.get('transcript')}")
            return out
        out.update(type=event, sites=[site])
    else:  # auto
        if not sites:
            out.update(type=None, stop=True,
                       reason=f"SpliceAI predicts no gained splice site above {DELTA_FLOOR} for this variant on "
                              f"{transcript.get('refseq') or transcript.get('transcript')} "
                              f"(strongest gain {strongest:.3f}). There is no aberrant site to block, so there is "
                              f"nothing for a splice-switching antisense oligonucleotide to aim at.")
            if strongest > 0:
                out["reason"] += (f" The strongest gain reported, {strongest:.3f}, is below the floor of "
                                  f"{DELTA_FLOOR} at which this screen will place a window.")
            return out
        if pseudoexon:
            out.update(type="pseudoexon", pseudoexon=pseudoexon, sites=[pseudoexon["acceptor"], pseudoexon["donor"]])
        else:
            site = max(sites.values(), key=lambda s: s["delta"])
            out.update(type=f"cryptic_{site['kind']}", sites=[site])

    chosen = out["sites"]
    failed = [s for s in chosen if not s["passed_min_delta"]]
    if not failed:
        out["confidence"] = ("above the ClinGen PP3 threshold of %s: every site this event rests on passes it"
                             % min_delta)
    elif len(failed) == len(chosen):
        out["confidence"] = (
            "BELOW the ClinGen PP3 threshold of %s: %s, in the uninformative band. The event is a hypothesis "
            "to test in RNA, not a finding."
            % (min_delta, "; ".join("the %s gain is %.3f" % (s["kind"], s["delta"]) for s in chosen)))
    else:
        out["confidence"] = (
            "PARTLY below the ClinGen PP3 threshold of %s: %s. A cryptic exon needs BOTH boundaries, so the "
            "event as a whole is a hypothesis to test in RNA, not a finding."
            % (min_delta, "; ".join("the %s gain is %.3f (%s)"
                                    % (s["kind"], s["delta"], "passes" if s["passed_min_delta"] else "below it")
                                    for s in chosen)))
    bad = [s for s in out["sites"] if not s["canonical_on_patient_allele"]]
    if bad:
        out["sequence_check"] = (
            "a predicted site is NOT canonical in the patient's own sequence: " +
            "; ".join(f"{s['kind']} at {s['position']} has {s['dinucleotide_patient_allele'] or '?'} where "
                      f"{s['dinucleotide_expected']} was expected" for s in bad) +
            ". Three things cause this: the site is genuinely non-canonical; the prediction's position convention "
            "does not match this transcript's strand; or the two bases that decide it lie outside the sequence "
            "that was fetched. Treat the windows below as unverified and check the position by hand."
        )
    else:
        out["sequence_check"] = ("every predicted site is canonical in the patient's own sequence "
                                 "(AG for an acceptor, GT/GC for a donor, read on the transcript strand)")
    return out


# ------------------------------------------------------------------ exon structure


def transcript_exons(transcript_id: str, assembly: str = "GRCh38") -> Outcome:
    """Exons of one Ensembl transcript, in transcript order."""
    from zebra.http import get_json
    from zebra.sources import record, validated_json

    tid = ensembl.strip_version(str(transcript_id or "").strip())
    if not tid:
        raise UsageError("give an Ensembl transcript id")
    url = f"{ensembl.host(assembly)}/lookup/id/{tid}"
    resp = get_json(url, source="Ensembl lookup", params={"expand": 1}, cache_ttl=30 * 86400, timeout=60)
    data = validated_json(resp, "Ensembl lookup", require="id")
    strand = 1 if int(data.get("strand") or 1) >= 0 else -1
    exons = []
    for e in data.get("Exon") or []:
        if isinstance(e, dict) and isinstance(e.get("start"), int):
            exons.append({"id": e.get("id"), "start": int(e["start"]), "end": int(e["end"])})
    exons.sort(key=lambda e: e["start"], reverse=strand < 0)
    for n, e in enumerate(exons, start=1):
        e["number"] = n
        e["length"] = e["end"] - e["start"] + 1
    return Outcome({"transcript": data.get("id"), "strand": strand, "chrom": data.get("seq_region_name"),
                    "exons": exons, "exon_count": len(exons)},
                   sources=[record("Ensembl lookup", tid, resp, note=f"{assembly}; {len(exons)} exons")])


def _exon_for(exons: Sequence[Dict[str, Any]], pos: int) -> Tuple[Optional[Dict[str, Any]], str]:
    for e in exons:
        if e["start"] <= pos <= e["end"]:
            return e, "the variant lies in this exon"
    if not exons:
        return None, "the transcript lookup returned no exons"
    nearest = min(exons, key=lambda e: min(abs(e["start"] - pos), abs(e["end"] - pos)))
    distance = min(abs(nearest["start"] - pos), abs(nearest["end"] - pos))
    return nearest, f"the variant is intronic; this is the nearest exon ({distance} bp away)"


# ------------------------------------------------------------------ candidate windows


def _site_windows(region: Region, site: Dict[str, Any], lengths: Sequence[int], stride: int) -> List[Dict[str, Any]]:
    """Windows that fully cover the site, its canonical dinucleotide included."""
    s = site["premrna_index"]
    out: List[Dict[str, Any]] = []
    for length in lengths:
        if site["kind"] == "acceptor":
            offsets = range(2, length)        # window start <= s-2 and window end >= s
        else:
            offsets = range(0, max(1, length - 2))  # window start <= s and window end >= s+2
        for k in offsets:
            if (k - (2 if site["kind"] == "acceptor" else 0)) % stride:
                continue
            i0 = s - k
            i1 = i0 + length - 1
            if not (region.holds(i0) and region.holds(i1)):
                continue
            out.append({"region": f"{site['kind']}_site", "site_kind": site["kind"], "site_position": site["position"],
                        "i0": i0, "i1": i1, "length": length, "site_premrna_index": s})
    return out


def _body_windows(region: Region, pseudoexon: Dict[str, Any], lengths: Sequence[int],
                  stride: int) -> List[Dict[str, Any]]:
    """Windows tiling the inside of the cryptic exon (a tiling, not an ESE prediction)."""
    a = pseudoexon["acceptor"]["premrna_index"]
    d = pseudoexon["donor"]["premrna_index"]
    out: List[Dict[str, Any]] = []
    for length in lengths:
        i0 = a
        while i0 + length - 1 <= d:
            i1 = i0 + length - 1
            # the exon's predicted boundaries can lie outside the sequence that actually
            # came back (a contig edge, or a donor beyond the fetched window); a window
            # there would slice to '' and be reported as a target
            if region.holds(i0) and region.holds(i1):
                out.append({"region": "pseudoexon_body", "site_kind": None, "site_position": None,
                            "i0": i0, "i1": i1, "length": length, "site_premrna_index": None})
            i0 += max(1, stride)
    return out


def _candidate(region: Region, window: Dict[str, Any], variant_pos: int,
               reference_region: Optional[Region] = None) -> Dict[str, Any]:
    i0, i1 = window["i0"], window["i1"]
    target = region.slice(i0, i1)
    start, end = region.span(i0, i1)
    aso = revcomp(target)
    hp = hairpin(aso)
    run = longest_run(aso)
    gc = gc_fraction(aso)
    variant_index = region.index(variant_pos)
    site_index = window["site_premrna_index"]
    out: Dict[str, Any] = {
        "region": window["region"],
        "target_genomic": {"chrom": region.chrom, "start": start, "end": end, "strand": "+"},
        "target_premrna_5to3": target,
        "antisense_5to3": aso,
        "length": window["length"],
        "gc_fraction": None if gc is None else round(gc, 3),
        "self_complementary_stem": hp["stem"],
        "hairpin_flag": hp["stem"] >= HAIRPIN_STEM,
        "hairpin_detail": hp,
        "longest_homopolymer": run["length"],
        "homopolymer_base": run["base"],
        "homopolymer_flag": run["length"] > HOMOPOLYMER_MAX,
        "contains_patient_variant": bool(i0 <= variant_index <= i1),
        "ambiguous_bases": sum(1 for b in target if b not in "ACGT") + max(0, window["length"] - len(target)),
        "target_complete": len(target) == window["length"],
    }
    if site_index is not None:
        offset_in_target = site_index - i0 + 1  # 1-based from the target's 5' end
        out.update(
            site_kind=window["site_kind"],
            site_position=window["site_position"],
            site_offset_in_target=offset_in_target,
            site_offset_in_antisense=window["length"] - (site_index - i0),
        )
    if out["contains_patient_variant"]:
        out["variant_offset_in_target"] = variant_index - i0 + 1
    if reference_region is not None:
        # The reference reading of the same window. Uniqueness is a property of the
        # GENOME, so that is what a genomic search has to be given: the patient's allele
        # is by definition not in the reference, and a search for it finds nothing (seen
        # live: the 25-mer carrying CFTR c.3718-2477C>T has no 16 nt exact seed on either
        # side of the substitution, so megablast returns no alignment at all).
        reference_target = reference_region.slice(i0, i1)
        out["reference_target_premrna_5to3"] = reference_target
        out["differs_from_reference"] = reference_target != target
    return out


# ------------------------------------------------------------------ ranking


def _linear_band(value: Optional[float], band: Tuple[float, float], hard: Tuple[float, float]) -> Optional[float]:
    if value is None:
        return None
    lo, hi = band
    hlo, hhi = hard
    if lo <= value <= hi:
        return 1.0
    if value < lo:
        return max(0.0, (value - hlo) / (lo - hlo)) if lo > hlo else 0.0
    return max(0.0, (hhi - value) / (hhi - hi)) if hhi > hi else 0.0


def _components(candidate: Dict[str, Any]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []

    def add(name: str, value: Optional[float], basis: str) -> None:
        # `basis` is this candidate's own reading; what the component MEANS is stated once
        # at result["ranking"]["components"][name], not repeated on every candidate
        weight, _why = RANK_COMPONENTS[name]
        rows.append({"name": name, "value": None if value is None else round(float(value), 3), "weight": weight,
                     "contribution": 0.0 if value is None else round(weight * float(value), 3),
                     "counted": value is not None, "basis": basis})

    length = candidate["length"]
    offset = candidate.get("site_offset_in_target")
    if offset is None:
        add("site_centering", None, "no splice site in this window (it tiles the exon body)")
    else:
        centre = (length + 1) / 2.0
        centering = 1.0 - (abs(offset - centre) / (centre - 1)) if centre > 1 else 1.0
        add("site_centering", max(0.0, min(1.0, centering)),
            f"the site sits at position {offset} of {length}; the middle is {centre:g}")
    gc = candidate.get("gc_fraction")
    add("gc_in_band", _linear_band(gc, GC_BAND, GC_HARD), f"GC fraction {gc}")
    add("no_hairpin", 0.0 if candidate["hairpin_flag"] else 1.0,
        f"longest self-complementary stem {candidate['self_complementary_stem']} "
        f"(flag at {HAIRPIN_STEM})")
    add("no_homopolymer", 0.0 if candidate["homopolymer_flag"] else 1.0,
        f"longest run {candidate['longest_homopolymer']}x{candidate['homopolymer_base']} "
        f"(flag above {HOMOPOLYMER_MAX})")
    uniq = candidate.get("uniqueness") or {}
    if uniq.get("checked") and uniq.get("unique") is not None:
        add("uniqueness", 1.0 if uniq["unique"] else 0.0,
            f"{uniq.get('locus_count')} GRCh38 primary-chromosome locus/loci")
    elif uniq.get("checked"):
        # the search ran but could not settle the question (a saturated subject list):
        # scoring that as 0 would penalise it as if it were known to be repeated
        add("uniqueness", None, "checked but not settled: %s"
            % uniq.get("unique_note", "the search could not determine uniqueness"))
    else:
        add("uniqueness", None, f"not checked: {uniq.get('reason', 'the uniqueness check was not run')}")
    add("variant_in_target", 1.0 if candidate["contains_patient_variant"] else 0.0,
        "the patient's variant is inside this window" if candidate["contains_patient_variant"]
        else "the patient's variant is outside this window")
    return rows


def _score(candidate: Dict[str, Any]) -> Dict[str, Any]:
    rows = _components(candidate)
    total = sum(r["contribution"] for r in rows)
    maximum = sum(r["weight"] for r in rows if r["counted"])
    left_out = [r["name"] for r in rows if not r["counted"] and r["weight"] > 0]
    return {
        "score": round(total, 3),
        "score_max": round(maximum, 3),
        "fraction": round(total / maximum, 3) if maximum else None,
        "components": rows,
        "components_not_counted": left_out,
    }


def _group_order(event: Dict[str, Any]) -> List[str]:
    sites = event.get("sites") or []
    created = [s for s in sites if s.get("created_by_variant")]
    others = [s for s in sites if not s.get("created_by_variant")]
    order = [f"{s['kind']}_site" for s in sorted(created, key=lambda s: -s["delta"])]
    order += [f"{s['kind']}_site" for s in sorted(others, key=lambda s: -s["delta"])]
    if event.get("type") == "pseudoexon":
        order.append("pseudoexon_body")
    seen, out = set(), []
    for name in order:
        if name not in seen:
            seen.add(name)
            out.append(name)
    return out


# ------------------------------------------------------------------ precedents


def precedents() -> Outcome:
    """Precedent and eligibility records, retrieved from Europe PMC now (never from memory)."""
    blocks: List[Dict[str, Any]] = []
    sources: List[Dict[str, Any]] = []
    warnings: List[str] = []
    for spec in PRECEDENT_QUERIES:
        out = attempt(f"Europe PMC ({spec['key']})",
                      lambda q=spec["query"]: europepmc.search(q, limit=PRECEDENT_LIMIT, result_type="core"),
                      warnings)
        block: Dict[str, Any] = {"key": spec["key"], "why_it_is_here": spec["why"], "query": spec["query"],
                                 "hits": []}
        if out is None:
            block["note"] = ("not retrieved in this session, so nothing is cited for it. Europe PMC was "
                             "unreachable or answered unusably; see the warnings.")
            blocks.append(block)
            continue
        sources.extend(out.sources)
        warnings.extend(out.warnings)
        for hit in out.result.get("hits") or []:
            abstract = hit.get("abstract")
            block["hits"].append({
                "pmid": hit.get("pmid"), "doi": hit.get("doi"), "title": hit.get("title"),
                "journal": hit.get("journal"), "year": hit.get("year"), "authors": hit.get("authors"),
                "url": hit.get("url"), "full_text_url": hit.get("fullTextUrl"),
                "abstract_excerpt": (abstract[:ABSTRACT_EXCERPT] + ("..." if len(abstract) > ABSTRACT_EXCERPT else ""))
                                    if abstract else None,
            })
        if not block["hits"]:
            block["note"] = ("the query returned no record in this session; nothing is cited for it")
        blocks.append(block)
    return Outcome({"rule": PRECEDENT_RULE, "records": blocks}, sources=sources, warnings=warnings)


# ------------------------------------------------------------------ the screen


def parse_lengths(raw: Optional[str]) -> Tuple[int, ...]:
    """`18-25`, `20`, or `18,20,25` -> the lengths to enumerate."""
    if raw is None or not str(raw).strip():
        return tuple(range(DEFAULT_LENGTHS[0], DEFAULT_LENGTHS[1] + 1))
    text = str(raw).strip()
    lo_limit, hi_limit = LENGTH_LIMITS
    m = re.match(r"^(\d{1,3})\s*[-:]\s*(\d{1,3})$", text)
    if m:
        lo, hi = int(m.group(1)), int(m.group(2))
        if not (lo_limit <= lo <= hi <= hi_limit):
            raise UsageError(f"--lengths must be LO-HI with {lo_limit} <= LO <= HI <= {hi_limit}")
        return tuple(range(lo, hi + 1))
    values = []
    for part in re.split(r"[,\s]+", text):
        if not part:
            continue
        if not part.isdigit():
            raise UsageError(f"--lengths takes LO-HI (e.g. 18-25) or a comma-separated list; cannot read {part!r}")
        value = int(part)
        if not lo_limit <= value <= hi_limit:
            raise UsageError(f"--lengths values must be between {lo_limit} and {hi_limit}")
        values.append(value)
    if not values:
        raise UsageError("--lengths is empty")
    return tuple(sorted(set(values)))


def screen(variant: str, assembly: str = "GRCh38", lengths: Optional[Sequence[int]] = None,
           event: str = "auto", min_delta: float = MIN_DELTA, top: int = DEFAULT_TOP,
           uniqueness: bool = False, distance: int = 500, stride: int = SITE_STRIDE,
           body_stride: int = BODY_STRIDE, timeout: float = 240.0,
           with_precedents: bool = True) -> Outcome:
    """Candidate antisense target windows for the aberrant splicing this variant may cause."""
    if assembly not in ("GRCh38", "GRCh37"):
        raise UsageError("assembly must be GRCh38 or GRCh37")
    if event not in EVENTS:
        raise UsageError(f"--event must be one of {', '.join(EVENTS)}")
    sizes = tuple(sorted(set(int(x) for x in (lengths or range(DEFAULT_LENGTHS[0], DEFAULT_LENGTHS[1] + 1)))))
    if not sizes or sizes[0] < LENGTH_LIMITS[0] or sizes[-1] > LENGTH_LIMITS[1]:
        raise UsageError(f"lengths must lie between {LENGTH_LIMITS[0]} and {LENGTH_LIMITS[1]}")
    stride = max(1, int(stride))
    body_stride = max(1, int(body_stride))
    try:
        top = int(top)
    except (TypeError, ValueError):
        raise UsageError("--top must be a whole number") from None
    if not 1 <= top <= 60:
        raise UsageError("--top must be between 1 and 60")
    try:
        min_delta = float(min_delta)
    except (TypeError, ValueError):
        raise UsageError("--min-delta must be a number") from None
    if not DELTA_FLOOR < min_delta <= 1.0:
        raise UsageError("--min-delta must be above the floor of %s and at most 1.0 (the ClinGen SVI threshold "
                         "is %s); below the floor the screen would place windows on sites that support no splice "
                         "effect at all" % (DELTA_FLOOR, MIN_DELTA))

    from zebra import s2f  # imported here so a change in s2f cannot break importing this module

    predicted = s2f.predict(variant, assembly=assembly, models=["spliceai", "pangolin"], distance=int(distance),
                            timeout=timeout)
    v = predicted.result["variant"]
    sources: List[Dict[str, Any]] = list(predicted.sources)
    warnings: List[str] = list(predicted.warnings)
    base: Dict[str, Any] = {
        "variant": v,
        "assembly": assembly,
        "lengths": list(sizes),
        "event_requested": event,
        "min_delta": min_delta,
        "preconditions": list(PRECONDITIONS),
        "limits": list(LIMITS),
        "what_this_is": ("a research feasibility screen for splice-switching antisense targets: predicted aberrant "
                         "splice sites, the windows on the pre-mRNA that cover them, and the reverse complement of "
                         "each. Not a design, not a drug, not advice."),
        "ranking": {"components": {k: {"weight": w, "what_it_means": why} for k, (w, why) in RANK_COMPONENTS.items()},
                    "group_order": REGION_ORDER_NOTE},
        "how_to_read_a_candidate": {
            "target_genomic": "forward-strand coordinates, 1-based inclusive",
            "target_premrna_5to3": "the target as the pre-mRNA reads it, on the transcript strand",
            "antisense_5to3": "the reverse complement of that target: a sequence, not a drug",
            "offsets": "1-based from the 5' end of the sequence named; the antisense strand runs the other way, "
                       "so site_offset_in_target and site_offset_in_antisense are mirror images",
            "rank.components": "each row is value x weight = contribution, with this candidate's own basis; what "
                               "the component measures is in ranking.components above",
            "rank.score": "the sum of the contributions and nothing else; score_max is the sum of the weights "
                          "that were counted, and components_not_counted names the ones left out of both",
        },
    }

    sai = _spliceai_row(predicted.result.get("models") or [])
    if not sai or sai.get("status") != "ran":
        base.update(event={"type": None, "stop": True,
                           "reason": f"SpliceAI did not run for this variant: "
                                     f"{(sai or {}).get('reason', 'the model was not in the result')}. "
                                     f"Without a splicing prediction there is no predicted event to screen against; "
                                     f"this is 'not measured', not 'no aberrant splicing'."},
                    candidates=[], candidate_count=0, stop=True)
        return Outcome(base, sources=sources, warnings=warnings, text=render(base),
                       query={"variant": variant, "assembly": assembly, "event": event})
    transcript = sai.get("transcript") if isinstance(sai.get("transcript"), dict) else {}
    base["transcript"] = {k: transcript.get(k) for k in
                          ("gene", "transcript", "refseq", "priority_label", "strand", "biotype")}
    strand = _strand_of(transcript.get("strand"))
    if strand is None:
        base.update(event={"type": None, "stop": True,
                           "reason": "the splicing result carries no readable strand for "
                                     "%s (it says %r). Every coordinate, the pre-mRNA itself and both "
                                     "dinucleotide checks depend on the strand, so nothing is reported rather "
                                     "than an answer built on a guessed '+'."
                                     % (transcript.get("refseq") or transcript.get("transcript") or "the "
                                        "transcript", transcript.get("strand"))},
                    candidates=[], candidate_count=0, stop=True)
        warnings.append("no readable transcript strand in the splicing result; the screen stopped")
        return Outcome(base, sources=sources, warnings=warnings, text=render(base),
                       query={"variant": variant, "assembly": assembly, "event": event})
    base["splice_prediction"] = {"model": "spliceai", "url": sai.get("url"), "distance": sai.get("distance"),
                                 "mask": sai.get("mask"), "deltas": transcript.get("deltas"),
                                 "reading": sai.get("reading")}
    pangolin = next((r for r in predicted.result.get("models") or [] if r.get("model") == "pangolin"), None)
    if pangolin and pangolin.get("status") == "ran":
        base["corroboration"] = {
            "model": "pangolin", "headline": pangolin.get("headline"), "url": pangolin.get("url"),
            "caveat": ("Pangolin is not independent of SpliceAI (Zeng & Li 2022, Genome Biol 23:103, state that "
                       "its architecture resembles SpliceAI's and both train on overlapping human splice-site "
                       "data), so agreement is weak corroboration, not a second axis of evidence."),
        }

    # Exon skipping aims at the exon's own splice sites, which can lie outside the
    # SpliceAI window, so the exon has to be known BEFORE the sequence window is chosen.
    exon_plan: Optional[Dict[str, Any]] = None
    if event == "exon_skip":
        exons_out = attempt("Ensembl transcript exons",
                            lambda: transcript_exons(str(transcript.get("transcript")), assembly), warnings)
        if exons_out is None:
            base.update(event={"type": None, "requested": event, "stop": True,
                               "reason": "the transcript's exon structure could not be retrieved, so the authentic "
                                         "splice sites of the exon to skip are unknown"},
                        candidates=[], candidate_count=0, stop=True)
            return Outcome(base, sources=sources, warnings=warnings, text=render(base),
                           query={"variant": variant, "assembly": assembly, "event": event})
        sources.extend(exons_out.sources)
        exon, why = _exon_for(exons_out.result["exons"], int(v["pos"]))
        if not exon:
            base.update(event={"type": None, "requested": event, "stop": True, "reason": why},
                        candidates=[], candidate_count=0, stop=True)
            return Outcome(base, sources=sources, warnings=warnings, text=render(base),
                           query={"variant": variant, "assembly": assembly, "event": event})
        exon_plan = {"exon": exon, "why": why, "exon_count": exons_out.result.get("exon_count")}

    # the window to fetch: everything a candidate could touch, plus the variant itself
    raw_deltas = transcript.get("deltas")
    positions = [int(d["position"]) for d in (raw_deltas if isinstance(raw_deltas, list) else [])
                 if isinstance(d, dict) and isinstance(d.get("position"), int)]
    positions.append(int(v["pos"]))
    if exon_plan:
        positions.extend([int(exon_plan["exon"]["start"]), int(exon_plan["exon"]["end"])])
    margin = sizes[-1] + 10
    lo, hi = min(positions) - margin, max(positions) + margin
    if hi - lo > 200_000:
        raise UsageError(f"the predicted positions span {hi - lo} bp around {v['chrom']}:{v['pos']}; "
                         f"lower --distance so the sequence window stays manageable")
    region, reference_region, seq_sources, seq_warnings = fetch_region(
        v["chrom"], lo, hi, strand, assembly, int(v["pos"]), v["ref"], v["alt"])
    sources.extend(seq_sources)
    warnings.extend(seq_warnings)
    base["sequence_window"] = {"chrom": region.chrom, "start": region.start, "end": region.end,
                               "transcript_strand": "+" if region.strand == 1 else "-",
                               "patient_allele_applied": region.allele_applied,
                               "length": len(region)}

    detected = attempt("reading the splicing prediction",
                       lambda: detect_event(region, reference_region, transcript, event, float(min_delta)),
                       warnings)
    if detected is None:
        base.update(event={"type": None, "stop": True,
                           "reason": "the splicing result could not be read into an event (its deltas are not in "
                                     "the shape this screen expects); see the warnings. Nothing is implied about "
                                     "whether this variant affects splicing."},
                    candidates=[], candidate_count=0, stop=True)
        return Outcome(base, sources=sources, warnings=warnings, text=render(base),
                       query={"variant": variant, "assembly": assembly, "event": event})
    base["event"] = detected

    windows: List[Dict[str, Any]] = []
    if detected.get("type") == "exon_skip":
        exon, why = exon_plan["exon"], exon_plan["why"]  # type: ignore[index]
        detected["exon"] = dict(exon, why_this_exon=why, of_exons=exon_plan.get("exon_count"))  # type: ignore[union-attr]
        terminal = ""
        if exon and exon.get("number") == 1:
            terminal = (" This is the FIRST exon: skipping it removes the transcription start and the start "
                        "codon, which frame arithmetic does not cover.")
        elif exon and exon_plan.get("exon_count") and exon["number"] == exon_plan["exon_count"]:
            terminal = (" This is the LAST exon: skipping it removes the stop codon and the 3' UTR, which frame "
                        "arithmetic does not cover.")
        detected["frame"] = None if not exon else {
            "exon_length": exon["length"], "in_frame": exon["length"] % 3 == 0,
            "note": ("removing this exon keeps the frame" if exon["length"] % 3 == 0 else
                     "removing this exon SHIFTS the frame: skipping it alone would not restore a reading frame")
                    + terminal,
        }
        sites = []
        if exon:
            # the authentic acceptor is the first exonic base in transcript order, the
            # authentic donor the last one
            first = exon["start"] if region.strand == 1 else exon["end"]
            last = exon["end"] if region.strand == 1 else exon["start"]
            for kind, pos in (("acceptor", first), ("donor", last)):
                index = region.index(pos)
                if not region.holds(index):
                    warnings.append(f"the authentic {kind} of exon {exon['number']} at {region.chrom}:{pos} lies "
                                    f"outside the sequence window fetched; no window is listed for it")
                    continue
                canonical = CANONICAL_ACCEPTOR if kind == "acceptor" else CANONICAL_DONOR
                patient = _dinucleotide(region, index, kind)
                on_ref = _dinucleotide(reference_region, index, kind)
                sites.append({
                    "kind": kind, "position": pos, "premrna_index": index, "spliceai_score": None, "delta": 0.0,
                    "delta_band": "not applicable: this is the exon's authentic site, not a predicted gain",
                    "passed_min_delta": False,
                    "dinucleotide_patient_allele": patient,
                    "dinucleotide_reference": on_ref,
                    "dinucleotide_expected": "AG 5' of the first exonic base" if kind == "acceptor"
                                             else "GT/GC 3' of the last exonic base",
                    "canonical_on_patient_allele": patient in canonical,
                    "canonical_on_reference": on_ref in canonical,
                    "minor_class_donor": bool(kind == "donor" and patient == "GC"),
                    "dinucleotide_changed_by_variant": bool(patient and on_ref and patient != on_ref),
                    "created_by_variant": False,
                    "created_by_variant_basis": "an authentic site of the transcript, not a site the variant made",
                    "authentic_site": True,
                })
        detected["sites"] = sites
        detected["confidence"] = ("not a prediction: these are the exon's authentic splice sites, taken from the "
                                  "Ensembl transcript structure. Whether skipping this exon helps is a question "
                                  "about the frame and the protein, not about the splicing model.")
        off = [s for s in sites if not s["canonical_on_patient_allele"]]
        detected["sequence_check"] = (
            "the exon's authentic sites are canonical in the sequence (AG for the acceptor, GT/GC for the donor)"
            if not off else
            "an authentic site is not canonical in the sequence: " +
            "; ".join(f"{s['kind']} at {s['position']} reads {s['dinucleotide_patient_allele'] or '?'}" for s in off)
            + ". The first exon has no acceptor and the last no donor, which is the usual reason; check which "
              "exon this is before reading the windows.")
        if not sites:
            base.update(candidates=[], candidate_count=0, stop=True)
            detected.update(stop=True, reason="no authentic splice site of the exon fell inside the sequence window")
            return Outcome(base, sources=sources, warnings=warnings, text=render(base),
                           query={"variant": variant, "assembly": assembly, "event": event})
        for site in sites:
            windows.extend(_site_windows(region, site, sizes, stride))
    elif detected.get("stop"):
        base.update(candidates=[], candidate_count=0, stop=True)
        return Outcome(base, sources=sources, warnings=warnings, text=render(base),
                       query={"variant": variant, "assembly": assembly, "event": event})
    else:
        for site in detected.get("sites") or []:
            windows.extend(_site_windows(region, site, sizes, stride))
        if detected.get("type") == "pseudoexon":
            windows.extend(_body_windows(region, detected["pseudoexon"], sizes, body_stride))

    if len(windows) > MAX_WINDOWS:
        # Truncating the tail would drop the exon body wholesale and, inside it, the
        # longest lengths: `windows` is site-then-body and length-ordered. Thin them out
        # round-robin over (region, length) instead, so what is kept still covers every
        # group and every length.
        buckets: Dict[Any, List[Dict[str, Any]]] = {}
        for w in windows:
            buckets.setdefault((w["region"], w["length"]), []).append(w)
        kept: List[Dict[str, Any]] = []
        order = list(buckets)
        while len(kept) < MAX_WINDOWS and any(buckets[k] for k in order):
            for key in order:
                if buckets[key] and len(kept) < MAX_WINDOWS:
                    kept.append(buckets[key].pop(0))
        warnings.append("%d windows were enumerated for these lengths and strides; %d were kept, thinned evenly "
                        "across every group and length rather than cut off at the end. Raise --stride / "
                        "--body-stride or narrow --lengths to cover the region evenly."
                        % (len(windows), len(kept)))
        base["windows_enumerated"] = len(windows)
        windows = kept
    # Drop windows that run into unknown bases or past the end of the sequence: a target
    # that is not the full length, or has N in it, is not a target. Partitioned in one
    # pass -- a `c not in dropped` filter over thousands of dicts is quadratic.
    candidates: List[Dict[str, Any]] = []
    short, ambiguous = 0, 0
    for w in windows:
        c = _candidate(region, w, int(v["pos"]), reference_region=reference_region)
        if not c["target_complete"]:
            short += 1
        elif c["ambiguous_bases"] or c["gc_fraction"] is None:
            ambiguous += 1
        else:
            candidates.append(c)
    if short or ambiguous:
        warnings.append("%d window(s) were dropped: %d ran past the end of the sequence available and %d "
                        "contained bases the reference does not resolve (N)"
                        % (short + ambiguous, short, ambiguous))
    # a site the sequence check failed must be visible on every window that covers it,
    # not only in the prose: a caller reading `--json | .candidates` sees these fields
    unverified = {s["position"] for s in (detected.get("sites") or [])
                  if not s.get("canonical_on_patient_allele")}
    for c in candidates:
        failed = c.get("site_position") in unverified if c.get("site_position") is not None else False
        c["site_canonical"] = None if c.get("site_position") is None else not failed
        if failed:
            c["sequence_check_failed"] = (
                "the site this window covers is not canonical in the patient's own sequence; the window is "
                "unverified (see event.sequence_check)")
    if not candidates:
        base.update(candidates=[], candidate_count=0, stop=True)
        base["event"]["stop"] = True
        base["event"].setdefault("reason", "no window of the requested lengths fits around the predicted site "
                                           "inside the sequence available")
        return Outcome(base, sources=sources, warnings=warnings, text=render(base),
                       query={"variant": variant, "assembly": assembly, "event": event})

    order = _group_order(detected) or sorted({c["region"] for c in candidates})
    rank_of = {name: i for i, name in enumerate(order)}
    for c in candidates:
        c["uniqueness"] = {"checked": False,
                           "reason": "not checked: rerun with --uniqueness to put the shortlist through one "
                                     "NCBI BLAST search against the human genome"}
        c["rank"] = _score(c)

    def sort_key(c):
        frac = c["rank"]["fraction"]
        return (rank_of.get(c["region"], len(order)), -(frac if frac is not None else 0.0),
                c["length"], c["target_genomic"]["start"])

    candidates.sort(key=sort_key)
    shortlist = candidates[:top]

    if uniqueness:
        from zebra.sources import ncbi_blast

        queries = {}
        names: Dict[int, str] = {}
        for i, c in enumerate(shortlist, start=1):
            name = "cand%02d" % i
            # search the REFERENCE reading of the window: uniqueness is about the genome,
            # and the patient's allele is not in it
            searched_for = c.get("reference_target_premrna_5to3") or c["target_premrna_5to3"]
            try:
                ncbi_blast.fasta({name: searched_for})
            except ValueError as err:
                c["uniqueness"] = {"checked": False,
                                   "reason": "this target cannot be sent to the search service: %s" % err}
                continue
            queries[name] = searched_for
            names[i] = name
        found = attempt("NCBI BLAST uniqueness", lambda: ncbi_blast.locate(queries), warnings) if queries else None
        if queries and found is None:
            for i, c in enumerate(shortlist, start=1):
                if i in names:
                    c["uniqueness"] = {"checked": False, "reason": "the NCBI BLAST uniqueness check failed; "
                                                                    "see the warnings"}
        elif found is not None:
            sources.extend(found.sources)
            warnings.extend(found.warnings)
            for i, c in enumerate(shortlist, start=1):
                if i in names:
                    row = dict(found.result.get(names[i],
                                                {"checked": False,
                                                 "reason": "BLAST returned no report for it"}))
                    row["searched_sequence"] = ("the reference reading of this window"
                                                if c.get("differs_from_reference")
                                                else "this window (reference and patient allele are the same here)")
                    if c.get("differs_from_reference"):
                        row["note"] = ("uniqueness is a property of the genome, so the reference reading of the "
                                       "window was searched; the antisense sequence above differs from it at the "
                                       "patient's variant, which makes the target MORE allele-specific, not less "
                                       "unique")
                    c["uniqueness"] = row
        for c in shortlist:
            c["rank"] = _score(c)
        shortlist.sort(key=sort_key)

    base.update(
        candidates=shortlist,
        candidate_count=len(candidates),
        candidates_shown=len(shortlist),
        groups=[{"region": name, "count": sum(1 for c in candidates if c["region"] == name)} for name in order],
        uniqueness_checked=bool(uniqueness),
        uniqueness_note=("uniqueness was checked for the shortlist only, with one NCBI BLAST search; the full "
                         f"candidate list ({len(candidates)}) was not checked"
                         if uniqueness else
                         "uniqueness was NOT checked for any candidate: nothing here is known to be unique in the "
                         "genome. Rerun with --uniqueness."),
        stop=False,
    )
    if with_precedents:
        prec = precedents()
        sources.extend(prec.sources)
        warnings.extend(prec.warnings)
        base["precedents"] = prec.result
    base["text_note"] = "every sequence here is a research hypothesis about a target, nothing more"
    return Outcome(base, sources=sources, warnings=warnings, text=render(base),
                   query={"variant": variant, "assembly": assembly, "event": event, "lengths": list(sizes),
                          "min_delta": min_delta, "top": top, "uniqueness": bool(uniqueness),
                          "distance": int(distance)})


# ------------------------------------------------------------------ rendering


def render(result: Dict[str, Any]) -> str:
    v = result["variant"]
    lines = [f"{v['input']} -> {v['chrom']}-{v['pos']}-{v['ref']}-{v['alt']} "
             f"({result['assembly']}; {v['resolved_from']})"]
    tx = result.get("transcript") or {}
    if tx.get("gene"):
        lines.append(f"{tx.get('gene')} {tx.get('refseq') or tx.get('transcript')} "
                     f"({tx.get('priority_label')}, strand {tx.get('strand')})")
    event = result.get("event") or {}
    if event.get("stop") or not event.get("type"):
        lines.append("")
        lines.append(f"no antisense target: {event.get('reason')}")
        for g in event.get("gains_considered") or []:
            lines.append(f"  gain considered: {g['spliceai_score']} {g['delta']:+.3f} at {g['position']} "
                         f"({g['delta_band']})")
        for loss in event.get("losses") or []:
            value = loss.get("delta")
            if isinstance(value, (int, float)) and value > 0:
                lines.append(f"  loss predicted: {loss.get('score')} {value:+.3f} at {loss.get('position')} "
                             f"- blocking a site cannot restore one the variant weakened; the route for a lost "
                             f"site is a different question (see --event exon_skip when removing the exon "
                             f"restores the frame)")
        lines += ["", *(f"! {p}" for p in result.get("preconditions") or [])]
        return "\n".join(lines)

    lines.append("")
    if event.get("type") == "pseudoexon":
        pe = event["pseudoexon"]
        lines.append(f"predicted event: cryptic exon of {pe['size_bp']} bp at {v['chrom']}:{pe['genomic'][0]}-"
                     f"{pe['genomic'][1]} ({'in frame' if pe['in_frame'] else 'frameshifting'})")
        lines.append(f"  {pe['frame_note']}")
    elif event.get("type") == "exon_skip":
        exon = event.get("exon") or {}
        lines.append(f"requested event: skip exon {exon.get('number')} "
                     f"({exon.get('start')}-{exon.get('end')}, {exon.get('length')} bp) — {exon.get('why_this_exon')}")
        if event.get("frame"):
            lines.append(f"  {event['frame']['note']}")
    else:
        lines.append(f"predicted event: {event.get('type')}")
    for site in event.get("sites") or []:
        tag = " (created by the patient's allele)" if site.get("created_by_variant") else ""
        delta = f"{site['spliceai_score']} {site['delta']:+.3f}" if site.get("spliceai_score") else "authentic site"
        lines.append(f"  {site['kind']:<8} {v['chrom']}:{site['position']}  {delta}  "
                     f"dinucleotide {site.get('dinucleotide_patient_allele')} on the patient allele / "
                     f"{site.get('dinucleotide_reference')} on the reference{tag}")
    lines.append(f"  confidence: {event.get('confidence')}")
    lines.append(f"  sequence check: {event.get('sequence_check')}")
    corr = result.get("corroboration")
    if corr:
        head = corr.get("headline") or {}
        lines.append(f"  pangolin: {head.get('score')} {head.get('value')} at {head.get('position')} "
                     f"(not an independent axis)")

    lines.append("")
    lines.append(f"{result['candidate_count']} candidate window(s); showing {result.get('candidates_shown')} "
                 f"(lengths {result['lengths'][0]}-{result['lengths'][-1]})")
    for g in result.get("groups") or []:
        lines.append(f"  group {g['region']:<16} {g['count']} window(s)")
    lines.append("")
    header = f"{'group':<16}{'len':>4} {'GC':>5} {'stem':>5} {'runs':>5} {'uniq':>6} {'score':>11}  target (genomic)"
    lines.append(header)
    for c in result.get("candidates") or []:
        uniq = c.get("uniqueness") or {}
        if not uniq.get("checked"):
            uniq_text = "-"                      # not checked
        elif uniq.get("unique") is None:
            uniq_text = "%s?" % uniq.get("locus_count")   # checked but not settled
        else:
            uniq_text = str(uniq.get("locus_count"))
        rank = c["rank"]
        tg = c["target_genomic"]
        if c.get("gc_fraction") is None:
            c = dict(c, gc_fraction=float("nan"))
        lines.append(f"{c['region']:<16}{c['length']:>4} {c['gc_fraction']:>5.2f} "
                     f"{c['self_complementary_stem']:>5} {c['longest_homopolymer']:>5} {uniq_text:>6} "
                     f"{rank['score']:>5.2f}/{rank['score_max']:<5.2f}  {tg['chrom']}:{tg['start']}-{tg['end']}")
        lines.append(f"{'':<16}  pre-mRNA 5'-{c['target_premrna_5to3']}-3'")
        lines.append(f"{'':<16}  antisense 5'-{c['antisense_5to3']}-3'"
                     + (f"   [covers the patient's variant at offset {c['variant_offset_in_target']} of the target]"
                        if c.get("contains_patient_variant") else ""))
    lines.append("")
    lines.append(f"uniqueness: {result.get('uniqueness_note')}")
    lines.append("ranking components (weight x value, nothing hidden):")
    for name, meta in (result.get("ranking") or {}).get("components", {}).items():
        lines.append(f"  {name:<20} weight {meta['weight']}  {meta['what_it_means']}")
    lines.append(f"  group order: {(result.get('ranking') or {}).get('group_order')}")
    prec = result.get("precedents")
    if prec:
        lines.append("")
        lines.append("precedent and eligibility records retrieved from Europe PMC in this run:")
        lines.append(f"  {prec['rule']}")
        for block in prec.get("records") or []:
            lines.append(f"  [{block['key']}] {block['why_it_is_here']}")
            if not block.get("hits"):
                lines.append(f"    {block.get('note')}")
            for hit in block.get("hits") or []:
                lines.append(f"    PMID {hit.get('pmid')} ({hit.get('year')}) {hit.get('journal')}: "
                             f"{hit.get('title')}")
                lines.append(f"      {hit.get('url')}")
    lines.append("")
    lines.append("preconditions:")
    for p in result.get("preconditions") or []:
        lines.append(f"  ! {p}")
    lines.append("limits:")
    for p in result.get("limits") or []:
        lines.append(f"  ! {p}")
    return "\n".join(lines)
