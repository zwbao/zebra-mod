"""Genomic QC for a family VCF: sex, relatedness, runs of homozygosity, Mendelian errors, UPD, mosaic de novo.

Everything here runs on this machine from the VCF alone: nothing is sent anywhere.
One streaming pass collects what every check needs; the checks then run on the
counts. Each method says where it comes from and what it assumes (`METHODS`).

Sites used (unless a check says otherwise): FILTER PASS or ".", biallelic SNVs,
each sample's call gated on depth (DP, else the AD sum) >= --min-dp and GQ >=
--min-gq when those fields are present. X/Y sex counts use every X/Y record,
as `zebra vcf triage` does.

Genotype codes below: 0 = hom-ref, 1 = het, 2 = hom-alt (a haploid ALT call is 2).

Mendelian / UPD table (trio child c, mother m, father f; autosomes):
    c  m  f   class
    1  0  0   de_novo_like   (an allele neither parent has)
    2  0  0   de_novo_like
    2  2  0   pat_absent     both alleles explained by the mother only: maternal UPD (hetero- or isodisomy),
    0  0  2   pat_absent     a paternal deletion, non-paternity or a swapped father sample
    2  1  0   pat_absent     (isodisomy-informative: the child is homozygous for one maternal allele)
    0  1  2   pat_absent
    0  2  0   mat_absent     the mirror images: paternal UPD, a maternal deletion, a swapped mother sample
    2  0  2   mat_absent
    2  0  1   mat_absent
    0  2  1   mat_absent
    1  2  2   de_novo_like   (a REF allele neither parent has: deletion or genotyping error)
    1  0  2 / 1 2 0 / 1 1 x / ...   consistent
Heterodisomy shows only at sites where the transmitting parent is homozygous
and the other parent opposite-homozygous (rows 3-4 / 7-8); isodisomy adds rows
5-6 / 9-10 and makes the child homozygous along the chromosome (a run of
homozygosity). So: excess pat_absent on one chromosome WITH a child ROH there =
maternal isodisomy; WITHOUT one = maternal heterodisomy (UPDio's logic, King
et al. 2014).
"""

from __future__ import annotations

import math
import os
import sys
import time
from array import array
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from zebra.core import Outcome, UsageError

AUTOSOMES = tuple(str(i) for i in range(1, 23))

METHODS: Dict[str, Dict[str, str]] = {
    "sex": {
        "cite": "zebra rule (cf. PLINK --check-sex)",
        "method": "X non-PAR heterozygous fraction, haploid X calls and Y non-PAR ALT calls (the rule zebra vcf "
                  "triage uses: zebra.vcf.sex_from_counts)",
        "source": "zebra's own rule; X/Y genotype sex checks as in PLINK --check-sex (Purcell et al. 2007, AJHG "
                  "81:559) without the F statistic, which needs population allele frequencies",
        "assumptions": "a germline sample; XXY, XO, mosaic aneuploidy and a caller that writes male X as diploid "
                       "het-free 1/1 are not distinguished from XY/XX",
    },
    "relatedness": {
        "cite": "KING-robust, Manichaikul 2010",
        "method": "KING-robust kinship (between-family estimator): phi = (N_AaAa - 2 N_AA,aa) / (2 N_Aa(i)) + 1/2 - "
                  "(N_Aa(i) + N_Aa(j)) / (4 N_Aa(i)), i the sample with fewer hets, over biallelic autosomal SNVs "
                  "called in both; degrees from KING's cut-offs >0.354 duplicate/MZ, [0.177,0.354] 1st, "
                  "[0.0884,0.177] 2nd, [0.0442,0.0884] 3rd",
        "source": "Manichaikul A et al. Robust relationship inference in genome-wide association studies. "
                  "Bioinformatics 2010;26:2867-2873; cut-offs from the KING manual "
                  "(kingrelatedness.com/manual.shtml, retrieved 2026-10-06)",
        "assumptions": "the estimator needs no allele frequencies and is robust to population structure; sites "
                       "where both samples are hom-ref contribute nothing, so a joint-called family VCF is enough. "
                       "Rare variants and genotyping errors pull phi slightly down, and phi is underestimated for "
                       "an inbred sample (Manichaikul 2010). Sites are taken up to a quota per chromosome (~18,000), "
                       "so a sorted genome is sampled across all autosomes rather than read from chr1 only. "
                       "Parent-offspring vs full siblings is told apart by opposite homozygotes per het (IBS0/het): "
                       "PO ~0 apart from errors, full sibs ~1/4 of an unrelated pair's value under HWE. The cut is "
                       "0.1 x the IBS0/het of the unrelated pairs in the same file when there are any, else 0.012, "
                       "with values up to twice the cut reported as unclear rather than as a problem "
                       "(zebra's thresholds, derived here, not KING's: KING's IBS0 is over all SNPs, which a VCF "
                       "does not hold). A verdict resting on fewer than 10 autosomes is marked inconclusive",
    },
    "roh": {
        "cite": "windowed, after PLINK --homozyg / AutoMap 2021",
        "method": "windowed homozygosity on the sample's own ALT-carrying biallelic SNVs (het vs hom-alt): a site "
                  "is in a run when some window of 20 consecutive such sites holds <= 1 het; runs are split at "
                  "gaps > 5 Mb and kept when >= 1 Mb and >= 20 sites; F_ROH = kept run length / autosomal length",
        "source": "in the spirit of PLINK --homozyg (Purcell et al. 2007) and AutoMap (Quinodoz M et al. Nat Commun "
                  "2021;12:518), not a re-implementation of either; no HMM (bcftools/RoH, Narasimhan V et al. "
                  "Bioinformatics 2016;32:1749) because that needs population allele frequencies",
        "assumptions": "exome caveats: informative sites cluster in genes and are absent from gene deserts, "
                       "centromeres and repeats, so run boundaries are only as precise as the nearest covered "
                       "exon, a run can bridge an uncovered stretch (up to 5 Mb) and small runs (< 1-2 Mb) are "
                       "unreliable; low-diversity regions look homozygous in everyone. Expected genome fraction in "
                       "runs for the offspring of: first cousins 1/16 (6.25%), second cousins 1/64 (1.6%), "
                       "uncle-niece or double first cousins 1/8 (12.5%) — averages with wide variance (Ceballos FC "
                       "et al. Nat Rev Genet 2018;19:220)",
    },
    "mendelian": {
        "cite": "UPDio logic, King 2014",
        "method": "per trio and chromosome, biallelic SNVs called in all three with an ALT in at least one "
                  "(informative sites): consistent / de_novo_like / pat_absent / mat_absent (table in zebra/qc.py); "
                  "a UPD flag needs >= 5 one-sided errors on one autosome, a Poisson tail p < 1e-5 against the "
                  "leave-one-chromosome-out background rate and a rate >= 5x that background; isodisomy when the "
                  "child's ROH covers >= 80% of the chromosome's informative span; 'segmental' (a deletion of the "
                  "other parent's allele, or segmental isodisomy) when >= 80% of the errors lie inside child ROH "
                  "that covers less; heterodisomy when the ROH covers < 10% and the errors lie outside it; mixed "
                  "otherwise. A whole-genome excess against one parent (one-sided binomial p < 1e-6 against the "
                  "mirror class, elevated on >= 10 chromosomes) is reported as parentage or a sample swap, never UPD; "
                  "a parent whose genotypes contradict the sex of the role disables the parent-of-origin call",
        "source": "the parent-of-origin logic of UPDio (King DA et al. Genome Res 2014;24:673-687); the thresholds "
                  "are zebra's",
        "assumptions": "autosomes only (X not assessed for UPD); a deletion on one parental chromosome gives the "
                       "same one-sided errors as UPD over its span, so a flag needs CNV review too; a whole-genome "
                       "excess against one parent is a sample swap or misattributed parentage, not UPD; segmental "
                       "or mosaic UPD (e.g. 11p15) can stay below a whole-chromosome test",
    },
    "mosaicism": {
        "cite": "binomial VAF, Acuna-Hidalgo 2015",
        "method": "for each de novo-like het (child 0/1, both parents 0/0 at depth), the child's ALT fraction "
                  "from AD with a 95% Wilson interval and a one-sided binomial test against 0.5; flagged as a "
                  "possible postzygotic mosaic when the interval's upper bound is below 0.40; a parent with >= 2 "
                  "ALT reads making >= 2% of its reads is flagged as possible parental mosaicism",
        "source": "allele-fraction testing as used by Acuna-Hidalgo R et al. Post-zygotic point mutations are an "
                  "underrecognized source of de novo genomic variation. AJHG 2015;97:67-74",
        "assumptions": "a diploid autosomal site, no copy-number change, no strand or mapping bias beyond the usual "
                       "slight reference bias (germline hets average ~0.45-0.5, hence 0.40, not 0.5); low depth "
                       "gives wide intervals and no flag; the de novo-like set of a raw VCF is dominated by "
                       "artefacts, so this is QC, not a candidate list (zebra vcf triage ranks candidates)",
    },
}

# imprinting disorders a UPD of that chromosome is known to cause (Eggermann T et al. Imprinting disorders:
# a group of congenital disorders with overlapping patterns of molecular changes affecting imprinted loci.
# Clin Epigenetics 2015;7:123)
IMPRINTED_UPD = {
    ("6", "paternal"): "transient neonatal diabetes mellitus (TNDM)",
    ("7", "maternal"): "Silver-Russell syndrome",
    ("11", "paternal"): "Beckwith-Wiedemann syndrome (usually segmental, mosaic 11p15)",
    ("11", "maternal"): "Silver-Russell syndrome (11p15, rare)",
    ("14", "maternal"): "Temple syndrome",
    ("14", "paternal"): "Kagami-Ogata syndrome",
    ("15", "maternal"): "Prader-Willi syndrome",
    ("15", "paternal"): "Angelman syndrome",
    ("20", "paternal"): "pseudohypoparathyroidism type 1B",
    ("20", "maternal"): "Mulchandani-Bhoj-Conlin syndrome (growth restriction)",
}

ROH_WINDOW = 20
ROH_MAX_HET = 1
ROH_MIN_BP = 1_000_000
ROH_MIN_SITES = 20
ROH_MAX_GAP = 5_000_000
PO_IBS0_PER_HET = 0.012  # below: parent-offspring-like; see METHODS["relatedness"]
KING_CUTS = ((0.354, "duplicate/MZ twin"), (0.177, "1st degree"), (0.0884, "2nd degree"), (0.0442, "3rd degree"))
UPD_MIN_ERRORS = 5
UPD_MAX_P = 1e-5
UPD_MIN_FOLD = 5.0
MOSAIC_UPPER = 0.40
KING_MAX_SITES = 400_000  # spread as a per-chromosome quota (a WGS has millions; kinship needs thousands)
KING_PER_CHROM = KING_MAX_SITES // 22
KING_MIN_CHROMS = 10  # fewer autosomes behind a kinship: the verdict is inconclusive (sibling variance)
MAX_SAMPLES = 12
MAX_ROH_SAMPLES = MAX_SAMPLES
ROH_MIN_INFORMATIVE = 300  # fewer ALT-carrying sites than this: ROH not assessed
QC_MAX_SECONDS = 900.0
MAX_POS = 3_000_000_000
PARENT_MOSAIC_MIN_FRACTION = 0.02

_G = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/0": 1, "1/1": 2, "1|1": 2, "0": 0, "1": 2}
_G_DIPLOID = {k: v for k, v in _G.items() if len(k) == 3}


# ------------------------------------------------------------------ PED

@dataclass
class Person:
    family: str
    iid: str
    father: Optional[str]
    mother: Optional[str]
    sex: Optional[str]  # "male" | "female" | None
    affected: Optional[bool]


class Pedigree:
    """A PED file: family, individual, father, mother, sex (1 male, 2 female, 0 unknown), phenotype (1/2/0/-9)."""

    def __init__(self, people: "OrderedDict[str, Person]", path: str = ""):
        self.people = people
        self.path = path

    def __contains__(self, iid: str) -> bool:
        return iid in self.people

    def get(self, iid: str) -> Optional[Person]:
        return self.people.get(iid)

    def children_of(self, iid: str) -> List[str]:
        return [p.iid for p in self.people.values() if iid in (p.father, p.mother)]

    def full_sibs(self, iid: str) -> List[str]:
        me = self.people.get(iid)
        if not me or not (me.father and me.mother):
            return []
        return [p.iid for p in self.people.values()
                if p.iid != iid and p.father == me.father and p.mother == me.mother]

    def trios(self, samples: Sequence[str]) -> List[Tuple[str, str, str]]:
        """(child, mother, father) for every child whose two parents are samples in the VCF."""
        have = set(samples)
        return [(p.iid, p.mother, p.father) for p in self.people.values()
                if p.iid in have and p.mother in have and p.father in have]

    def expected(self, a: str, b: str) -> Optional[str]:
        """'parent-offspring', 'full siblings', 'unrelated (parents of one child)' or None when the PED does not say."""
        pa, pb = self.people.get(a), self.people.get(b)
        if not pa or not pb:
            return None
        if a in (pb.father, pb.mother) or b in (pa.father, pa.mother):
            return "parent-offspring"
        if pa.father and pa.mother and pa.father == pb.father and pa.mother == pb.mother:
            return "full siblings"
        if (pa.father and pa.father == pb.father) or (pa.mother and pa.mother == pb.mother):
            return "half siblings"
        kids_a, kids_b = set(self.children_of(a)), set(self.children_of(b))
        if kids_a & kids_b:
            return "unrelated (parents of one child)"
        return None

    def as_dict(self) -> Dict[str, Any]:
        return {"path": self.path, "individuals": [
            {"family": p.family, "id": p.iid, "father": p.father, "mother": p.mother, "sex": p.sex,
             "affected": p.affected} for p in self.people.values()]}


def read_ped(path: str) -> Pedigree:
    p = os.path.expanduser(path)
    if not os.path.isfile(p):
        raise UsageError(f"PED file not found: {path}")
    people: "OrderedDict[str, Person]" = OrderedDict()
    with open(p, encoding="utf-8", errors="replace") as fh:
        for n, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = [c.strip() for c in line.split("\t")] if "\t" in line else line.split()
            if len(cols) < 6:
                raise UsageError(f"{path} line {n}: a PED line needs 6 columns (family, individual, father, mother, "
                                 f"sex, phenotype), got {len(cols)}")
            fam, iid, fat, mot, sex, phen = cols[:6]
            if iid in people:
                raise UsageError(f"{path} line {n}: individual {iid!r} appears twice")
            if iid in (fat, mot):
                raise UsageError(f"{path} line {n}: {iid!r} is named as their own parent")
            people[iid] = Person(fam, iid, None if fat in ("0", "-9", ".") else fat,
                                 None if mot in ("0", "-9", ".") else mot,
                                 {"1": "male", "2": "female"}.get(sex),
                                 {"2": True, "1": False}.get(phen))
    if not people:
        raise UsageError(f"{path}: no individuals in the PED file")
    for person in people.values():
        if person.father and person.father in people and people[person.father].sex == "female":
            raise UsageError(f"{path}: {person.father!r} is {person.iid!r}'s father but has sex 2 (female)")
        if person.mother and person.mother in people and people[person.mother].sex == "male":
            raise UsageError(f"{path}: {person.mother!r} is {person.iid!r}'s mother but has sex 1 (male)")
    return Pedigree(people, str(Path(p).resolve()))


# ------------------------------------------------------------ statistics

def wilson(k: int, n: int, z: float = 1.959964) -> Tuple[float, float]:
    """95% Wilson score interval for k successes in n trials."""
    if n <= 0:
        return 0.0, 1.0
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def binom_cdf(k: int, n: int, p: float) -> float:
    """P(X <= k) for X ~ Binomial(n, p), summed in log space."""
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    if p <= 0:
        return 1.0
    if p >= 1:
        return 0.0
    if n > 20_000:  # continuity-corrected normal approximation: exact sums over millions of reads hang
        z = (k + 0.5 - n * p) / math.sqrt(n * p * (1 - p))
        return 0.5 * math.erfc(-z / math.sqrt(2))
    lp, lq = math.log(p), math.log(1 - p)
    terms = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * lp + (n - i) * lq
             for i in range(0, k + 1)]
    top = max(terms)
    return min(1.0, math.exp(top) * sum(math.exp(t - top) for t in terms))


def poisson_sf(k: int, lam: float) -> float:
    """P(X >= k) for X ~ Poisson(lam)."""
    if k <= 0:
        return 1.0
    if lam <= 0:
        return 0.0
    if k <= lam:
        below = sum(math.exp(-lam + i * math.log(lam) - math.lgamma(i + 1)) for i in range(k))
        return max(0.0, 1.0 - below)
    total, i = 0.0, k
    while True:
        term = math.exp(-lam + i * math.log(lam) - math.lgamma(i + 1))
        total += term
        if term < 1e-300 or term < total * 1e-15 or i > k + 100_000:
            return min(1.0, total)
        i += 1


def mosaic_assessment(ref: Optional[int], alt: Optional[int]) -> Optional[Dict[str, Any]]:
    """ALT fraction, its 95% Wilson interval and the one-sided binomial p against a germline het (0.5)."""
    if ref is None or alt is None or ref < 0 or alt < 0 or ref + alt <= 0:
        return None
    n = ref + alt
    lo, hi = wilson(alt, n)
    return {"alt_reads": alt, "depth": n, "vaf": round(alt / n, 3), "ci95": [round(lo, 3), round(hi, 3)],
            "p_below_half": float(f"{binom_cdf(alt, n, 0.5):.3g}"), "possible_mosaic": hi < MOSAIC_UPPER}


# --------------------------------------------------------------------- ROH

def find_roh(positions: Sequence[int], hets: Sequence[int], window: int = ROH_WINDOW, max_het: int = ROH_MAX_HET,
             min_bp: int = ROH_MIN_BP, min_sites: int = ROH_MIN_SITES, max_gap: int = ROH_MAX_GAP
             ) -> List[Dict[str, Any]]:
    """Runs of homozygosity on one chromosome from sorted ALT-carrying sites (hets[i] 1 = het, 0 = hom-alt)."""
    n = len(positions)
    if n < max(window, min_sites):
        return []
    pre = [0] * (n + 1)
    for i in range(n):
        pre[i + 1] = pre[i] + (1 if hets[i] else 0)
    cover = [0] * (n + 1)
    for s in range(0, n - window + 1):
        if pre[s + window] - pre[s] <= max_het:
            cover[s] += 1
            cover[s + window] -= 1
    member = []
    run = 0
    for i in range(n):
        run += cover[i]
        member.append(run > 0)
    segs: List[Dict[str, Any]] = []
    i = 0
    while i < n:
        if not member[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and member[j + 1] and positions[j + 1] - positions[j] <= max_gap:
            j += 1
        a, b = i, j
        while a <= b and hets[a]:
            a += 1
        while b >= a and hets[b]:
            b -= 1
        if a <= b:
            span = positions[b] - positions[a] + 1
            sites = b - a + 1
            if span >= min_bp and sites >= min_sites:
                segs.append({"start": positions[a], "end": positions[b], "mb": round(span / 1e6, 2),
                             "sites": sites, "hets": int(pre[b + 1] - pre[a])})
        i = j + 1
    return segs


def autosomal_length(assembly: Optional[str]) -> int:
    from zebra.vcf import CHROM_LENGTHS

    table = CHROM_LENGTHS.get(assembly or "GRCh38") or CHROM_LENGTHS["GRCh38"]
    return sum(table[c] for c in AUTOSOMES)


def roh_summary(by_chrom: Dict[str, List[Dict[str, Any]]], spans: Dict[str, Tuple[int, int]],
                assembly: Optional[str]) -> Dict[str, Any]:
    segs = [dict(s, chrom=c) for c in AUTOSOMES for s in by_chrom.get(c, [])]
    total_bp = sum(s["end"] - s["start"] + 1 for s in segs)
    f_roh = total_bp / autosomal_length(assembly)
    per_chrom = {c: sum(s["end"] - s["start"] + 1 for s in by_chrom.get(c, [])) for c in AUTOSOMES}
    cover = {}
    for c, bp in per_chrom.items():
        lo_hi = spans.get(c)
        cover[c] = round(bp / max(1, lo_hi[1] - lo_hi[0] + 1), 3) if lo_hi else 0.0
    top_chrom = max(per_chrom, key=lambda c: per_chrom[c]) if segs else None
    pct = 100 * f_roh
    if pct < 0.5:
        reading = "no evidence of close parental relatedness"
    elif pct < 2:
        reading = "a little homozygosity: distant relatedness, a founder population, or noise in an exome"
    elif pct < 5:
        reading = "consistent with parents related about as second cousins or closer"
    elif pct < 10:
        reading = "consistent with parents related about as first cousins"
    else:
        reading = "consistent with parents related more closely than first cousins (uncle-niece, double first cousins)"
    upd_hint = None
    if top_chrom and per_chrom[top_chrom] >= 10_000_000 and per_chrom[top_chrom] >= 0.7 * total_bp \
            and cover.get(top_chrom, 0) >= 0.5 and (total_bp - per_chrom[top_chrom]) < 10_000_000:
        upd_hint = (f"chr{top_chrom} carries {per_chrom[top_chrom] / 1e6:.1f} of {total_bp / 1e6:.1f} Mb of run "
                    f"length ({cover[top_chrom]:.0%} of its covered span) and the rest of the genome little: that "
                    "pattern points to uniparental isodisomy of that chromosome rather than related parents — "
                    "test it with the parents' genotypes (zebra qc with a trio) or methylation/SNP array")
    segs.sort(key=lambda s: -(s["end"] - s["start"]))
    return {"segments": segs[:30], "n_segments": len(segs), "total_mb": round(total_bp / 1e6, 2),
            "f_roh": round(f_roh, 4), "reading": reading if not upd_hint else reading + "; but see upd_hint",
            "upd_hint": upd_hint,
            "chromosomes": {f"chr{c}": {"roh_mb": round(per_chrom[c] / 1e6, 2), "covered_span_fraction": cover[c]}
                            for c in AUTOSOMES if per_chrom[c]}}


# ---------------------------------------------------------------- kinship

def king_robust(het_het: int, opp_hom: int, het_a: int, het_b: int) -> Optional[float]:
    lo = min(het_a, het_b)
    if lo <= 0:
        return None
    return (het_het - 2 * opp_hom) / (2 * lo) + 0.5 - (het_a + het_b) / (4 * lo)


def degree(phi: Optional[float]) -> str:
    """KING's ranges (closed at the lower bound, as the manual writes them): >0.354, [0.177, 0.354], ..."""
    if phi is None:
        return "unknown"
    if phi > KING_CUTS[0][0]:
        return KING_CUTS[0][1]
    for cut, label in KING_CUTS[1:]:
        if phi >= cut:
            return label
    return "unrelated (beyond 3rd degree)"


def _pair_verdict(a: str, b: str, phi: Optional[float], ibs0: Optional[float], expected: Optional[str],
                  ped: Optional["Pedigree"], po_cut: float = PO_IBS0_PER_HET) -> Tuple[str, Optional[str]]:
    deg = degree(phi)
    kind = deg
    if deg == "1st degree" and ibs0 is not None:
        # between the cut and twice it the two cannot be told apart reliably (rare variants, genotyping errors)
        kind = ("1st degree, parent-offspring-like" if ibs0 < po_cut else
                "1st degree (parent-offspring vs full siblings unclear)" if ibs0 < 2 * po_cut else
                "1st degree, full-sibling-like")
    if expected is None:
        return kind, None
    if phi is None:
        return kind, (f"the PED expects {a}–{b} to be {expected}, but no kinship could be estimated (no shared "
                      "heterozygous SNVs): the relationship is unchecked")
    if expected == "parent-offspring":
        parent = None
        if ped is not None:
            pa, pb = ped.get(a), ped.get(b)
            parent = a if pb and a in (pb.father, pb.mother) else (b if pa and b in (pa.father, pa.mother) else None)
        role = "father" if ped is not None and parent and any(p.father == parent for p in ped.people.values()) \
            else ("mother" if parent else "parent")
        if deg == "duplicate/MZ twin":
            return kind, (f"{a} and {b} look like the same person (kinship {phi:.3f}): a duplicated or swapped "
                          "sample, not parent and child")
        if deg.startswith("unrelated"):
            why = ("non-paternity or a swapped/mislabelled father sample" if role == "father" else
                   "a swapped or mislabelled sample (or egg donation/surrogacy)" if role == "mother" else
                   "a swapped or mislabelled sample, or misattributed parentage")
            return kind, (f"the PED says {a}–{b} are parent and child but they look unrelated (kinship {phi:.3f}): "
                          f"{why}. Every inherited/de novo call that uses this {role} is unreliable")
        if deg in ("2nd degree", "3rd degree"):
            return kind, (f"the PED says {a}–{b} are parent and child but they look {deg} relatives (kinship "
                          f"{phi:.3f}): a sample from another relative (grandparent, aunt/uncle) or misattributed "
                          "parentage")
        if kind.endswith("full-sibling-like"):
            return kind, (f"{a}–{b} are first-degree relatives but share opposite homozygotes ({ibs0:.3f} per het, "
                          f"cut {po_cut:.4f}) like full siblings, not like parent and child: check the sample labels")
        return kind, None
    if expected == "full siblings":
        if deg.startswith("unrelated"):
            return kind, f"the PED says {a}–{b} are full siblings but they look unrelated (kinship {phi:.3f}): sample swap"
        if deg == "2nd degree":
            return kind, (f"the PED says {a}–{b} are full siblings but kinship {phi:.3f} is 2nd degree: half "
                          "siblings (a different father or mother) or a sample from another relative")
        if deg == "duplicate/MZ twin":
            return kind, f"{a}–{b} look identical (kinship {phi:.3f}): MZ twins or a duplicated sample"
        if kind.endswith("parent-offspring-like"):
            return kind, (f"{a}–{b} share almost no opposite homozygotes ({ibs0:.4f} per het, cut {po_cut:.4f}): "
                          "parent and child rather than siblings — check the labels")
        return kind, None
    if expected == "half siblings":
        if deg.startswith("unrelated"):
            return kind, f"the PED says {a}–{b} are half siblings but they look unrelated (kinship {phi:.3f})"
        return kind, None
    if expected.startswith("unrelated"):
        if not deg.startswith("unrelated"):
            return kind, (f"the parents {a} and {b} look related ({deg}, kinship {phi:.3f}): consanguinity — "
                          "expect runs of homozygosity in their children and weigh recessive causes")
    return kind, None


# --------------------------------------------------------------- the scan

def _me_class(c: int, m: int, f: int) -> str:
    pair = {0: (0, 0), 1: (0, 1), 2: (1, 1)}[c]
    am = {0: (0,), 1: (0, 1), 2: (1,)}[m]
    af = {0: (0,), 1: (0, 1), 2: (1,)}[f]
    x, y = pair
    if (x in am and y in af) or (y in am and x in af):
        return "ok"
    mat_uni = x in am and y in am
    pat_uni = x in af and y in af
    if mat_uni and not pat_uni:
        return "pat_absent"
    if pat_uni and not mat_uni:
        return "mat_absent"
    return "de_novo_like" if not (mat_uni or pat_uni) else "other"


ME_TABLE = {(c, m, f): _me_class(c, m, f) for c in (0, 1, 2) for m in (0, 1, 2) for f in (0, 1, 2)}


def _progress(progress: Optional[Callable[[str], None]], msg: str) -> None:
    if progress:
        progress(msg)


def stderr_progress(prefix: str) -> Callable[[str], None]:
    """A progress printer to stderr (stdout carries the JSON envelope); ZEBRA_PROGRESS=0 silences it."""
    def emit(msg: str) -> None:
        if os.environ.get("ZEBRA_PROGRESS", "1") == "0":
            return
        try:
            sys.stderr.write(f"[{prefix}] {msg}\n")
            sys.stderr.flush()
        except (OSError, ValueError):
            pass
    return emit


def run(path: str, ped: Optional[str] = None, proband: Optional[str] = None, mother: Optional[str] = None,
        father: Optional[str] = None, sex: Optional[str] = None, assembly: Optional[str] = None, min_dp: int = 10,
        min_gq: int = 20, max_samples: int = MAX_SAMPLES, max_seconds: float = QC_MAX_SECONDS,
        progress: Optional[Callable[[str], None]] = None) -> Outcome:
    """One pass over a (family) VCF: sex, KING-robust kinship, ROH, Mendelian errors and UPD, mosaic de novo."""
    from zebra import vcf as V

    t0 = time.monotonic()
    warnings: List[str] = []
    notes: List[str] = []
    if sex not in (None, "male", "female"):
        raise UsageError("--sex must be male or female")
    pedigree = read_ped(ped) if ped else None  # read before the VCF is opened: an error leaks no handle
    header, records = V.iter_records(path)
    try:
        samples_vcf = list(header.samples)
        if not samples_vcf:
            raise UsageError("the VCF has no sample columns: QC needs genotypes")
        for role, s in (("proband", proband), ("mother", mother), ("father", father)):
            if s is not None and s not in samples_vcf:
                raise UsageError(f"{role} {s!r} is not a sample in the VCF (samples: {', '.join(samples_vcf[:20])})")
        dup = sorted({s for s in samples_vcf if samples_vcf.count(s) > 1})
        if len({x for x in (proband, mother, father) if x}) < len([x for x in (proband, mother, father) if x]):
            raise UsageError("proband, mother and father must be different samples")
        # which samples, which trios: the proband's trio first, so a cap never drops it
        trios: List[Tuple[str, str, str]] = []
        if pedigree is not None:
            in_both = [s for s in pedigree.people if s in samples_vcf]
            if not in_both:
                raise UsageError(f"none of the PED individuals is a sample in the VCF (PED: "
                                 f"{', '.join(list(pedigree.people)[:10])}; VCF: {', '.join(samples_vcf[:10])})")
            absent = [s for s in pedigree.people if s not in samples_vcf]
            if absent:
                notes.append(f"PED individuals not in the VCF (not assessed): {', '.join(absent[:10])}")
            trios = pedigree.trios(samples_vcf)
            pool = in_both
            if proband and proband not in pedigree:
                warnings.append(f"--proband {proband} is not in the PED file: assessed from the VCF alone (no "
                                "expected relationships for it)")
                pool = [proband] + pool
        else:
            pool = [s for s in (proband, mother, father) if s] or list(samples_vcf)
        if proband and mother and father and (proband, mother, father) not in trios:
            trios.append((proband, mother, father))
        elif (mother or father) and proband and pedigree is None:
            notes.append("a duo: Mendelian errors and UPD need both parents")
        if sex and not proband:
            warnings.append("--sex states the proband's sex, but no --proband was given: it was not used")
        trios.sort(key=lambda t: t[0] != proband)
        chosen: List[str] = []
        for s in ([proband] if proband else []) + [x for t in trios for x in t] + list(pool):
            if s and s not in chosen:
                chosen.append(s)
        if len(chosen) > max_samples:
            warnings.append(f"{len(chosen)} samples: QC covers {max_samples} ({', '.join(chosen[:max_samples])}), the "
                            "proband's trio first; raise --max-samples to assess more")
            chosen = chosen[:max_samples]
            trios = [t for t in trios if all(s in chosen for s in t)]
        if [s for s in chosen if s in dup]:
            raise UsageError(f"sample name(s) {', '.join(s for s in chosen if s in dup)} appear more than once in the "
                             "VCF header: genotypes cannot be assigned to a person")
        build = V.guess_build(header)
        if assembly is None:
            assembly = build["guess"]
            if assembly is None:
                notes.append("genome build not stated in the header: GRCh38 pseudoautosomal bounds and lengths "
                             "assumed (pass --assembly)")
        elif build["guess"] and build["guess"] != assembly:
            raise UsageError(f"--assembly {assembly} contradicts the VCF header ({build['guess']})")
    except BaseException:
        records.close()
        raise
    asm = assembly or "GRCh38"
    idx = {s: samples_vcf.index(s) for s in chosen}
    roh_samples = [t[0] for t in trios] + ([proband] if proband else [])
    roh_samples = list(OrderedDict.fromkeys(roh_samples or chosen))[:MAX_ROH_SAMPLES]
    pairs = [(a, b) for i, a in enumerate(chosen) for b in chosen[i + 1:]]
    from zebra.vcf import GVCF_PLACEHOLDERS

    # accumulators
    sexc: Dict[str, Counter] = {s: Counter() for s in chosen}
    # per pair and chromosome: n, hethet, opposite hom, het a, het b (per chromosome, so a chromosome
    # found to be uniparental can be left out of the parent-child estimate)
    kin: Dict[Tuple[str, str], Dict[str, List[int]]] = {pr: defaultdict(lambda: [0, 0, 0, 0, 0]) for pr in pairs}
    king_by_chrom: Counter = Counter()
    roh_pos: Dict[str, Dict[str, array]] = {s: defaultdict(lambda: array("q")) for s in roh_samples}
    roh_het: Dict[str, Dict[str, bytearray]] = {s: defaultdict(bytearray) for s in roh_samples}
    me: Dict[Tuple[str, str, str], Dict[str, Counter]] = {t: defaultdict(Counter) for t in trios}
    me_pos: Dict[Tuple[str, str, str], Dict[Tuple[str, str], List[int]]] = {t: defaultdict(list) for t in trios}
    denovo: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = {t: [] for t in trios}
    denovo_n: Counter = Counter()
    per_sample: Dict[str, Counter] = {s: Counter() for s in chosen}
    n_rec = n_snv = no_depth = bad_pos = 0
    complete = True
    stopped_by = None
    fmt_cache: Dict[str, Tuple[int, int, int, int]] = {}
    last_note = time.monotonic()
    try:
        for rec in records:
            n_rec += 1
            if n_rec % 100_000 == 0:
                now = time.monotonic()
                if now - t0 > max_seconds:
                    complete, stopped_by = False, f"the {max_seconds:g} s time cap"
                    break
                if now - last_note > 5:
                    _progress(progress, f"{n_rec:,} records read, chr{rec.chrom}")
                    last_note = now
            V.sex_tally(rec, idx, sexc, asm)
            # a gVCF writes <NON_REF>/<*> after the real ALT: the site is still biallelic for this purpose
            if len(rec.alts) > 1 and all(a.upper() in GVCF_PLACEHOLDERS for a in rec.alts[1:]):
                alt0 = rec.alts[0]
            elif len(rec.alts) == 1:
                alt0 = rec.alts[0]
            else:
                continue
            if len(rec.ref) != 1 or len(alt0) != 1 or alt0 in ("*", ".") or not alt0.isalpha() \
                    or rec.filter not in ("PASS", ".") or not rec.fmt:
                continue
            chrom = rec.chrom
            if chrom not in AUTOSOMES:
                continue
            if rec.pos > MAX_POS:
                bad_pos += 1
                continue
            key = ":".join(rec.fmt)
            fi = fmt_cache.get(key)
            if fi is None:
                f = rec.fmt
                fi = (f.index("GT") if "GT" in f else -1, f.index("DP") if "DP" in f else -1,
                      f.index("GQ") if "GQ" in f else -1, f.index("AD") if "AD" in f else -1)
                fmt_cache[key] = fi
            gti, dpi, gqi, adi = fi
            if gti < 0:
                continue
            n_snv += 1
            g: Dict[str, Optional[int]] = {}
            ad: Dict[str, Tuple[Optional[int], Optional[int]]] = {}
            for s, i in idx.items():
                if i >= len(rec.sample_fields):
                    g[s] = None
                    continue
                parts = rec.sample_fields[i].split(":")
                code = _G_DIPLOID.get(parts[gti] if gti < len(parts) else ".")  # a haploid autosome call is not a genotype
                if code is None:
                    g[s] = None
                    continue
                ref_r = alt_r = None
                if 0 <= adi < len(parts) and parts[adi] not in ("", "."):
                    vals = parts[adi].split(",")
                    if len(vals) >= 2:
                        try:
                            ref_r, alt_r = int(vals[0]), int(vals[1])
                        except ValueError:
                            ref_r = alt_r = None
                        if ref_r is not None and (ref_r < 0 or alt_r < 0):  # type: ignore[operator]
                            ref_r = alt_r = None
                dp = None
                if 0 <= dpi < len(parts) and parts[dpi] not in ("", "."):
                    try:
                        dp = int(float(parts[dpi]))
                    except (ValueError, OverflowError):
                        dp = None
                if dp is None and ref_r is not None and alt_r is not None:
                    dp = ref_r + alt_r
                if dp is None:
                    no_depth += 1
                elif dp < min_dp:
                    g[s] = None
                    continue
                if 0 <= gqi < len(parts) and parts[gqi] not in ("", "."):
                    try:
                        gq = float(parts[gqi])
                        if gq != gq or gq < min_gq:  # NaN fails the gate too
                            g[s] = None
                            continue
                    except ValueError:
                        pass
                g[s] = code
                ad[s] = (ref_r, alt_r)
                per_sample[s]["called"] += 1
                if code == 1:
                    per_sample[s]["het"] += 1
                elif code == 2:
                    per_sample[s]["hom_alt"] += 1
            for s in roh_samples:
                code = g.get(s)
                if code in (1, 2):
                    roh_pos[s][chrom].append(rec.pos)
                    roh_het[s][chrom].append(1 if code == 1 else 0)
            if king_by_chrom[chrom] < KING_PER_CHROM:
                used = False
                for pr, per_chrom in kin.items():
                    ga, gb = g.get(pr[0]), g.get(pr[1])
                    if ga is None or gb is None:
                        continue
                    used = True
                    cnt = per_chrom[chrom]
                    cnt[0] += 1
                    if ga == 1:
                        cnt[3] += 1
                        if gb == 1:
                            cnt[1] += 1
                    if gb == 1:
                        cnt[4] += 1
                    if ga + gb == 2 and ga != 1:
                        cnt[2] += 1
                if used:
                    king_by_chrom[chrom] += 1
            for t in trios:
                c, m, f = g.get(t[0]), g.get(t[1]), g.get(t[2])
                if c is None or m is None or f is None or (c == m == f == 0):
                    continue  # uncalled, or no ALT in the trio: not informative
                cls = ME_TABLE[(c, m, f)]
                bucket = me[t][chrom]
                bucket["sites"] += 1
                bucket[cls] += 1
                if cls in ("pat_absent", "mat_absent"):
                    lst = me_pos[t][(chrom, cls)]
                    if len(lst) < 20_000:
                        lst.append(rec.pos)
                if cls == "de_novo_like" and c == 1 and m == 0 and f == 0:
                    denovo_n[t] += 1
                    if len(denovo[t]) < 2_000:
                        ref_r, alt_r = ad.get(t[0], (None, None))
                        mo, fa = ad.get(t[1], (None, None)), ad.get(t[2], (None, None))
                        denovo[t].append({"variant": f"{chrom}-{rec.pos}-{rec.ref}-{alt0}",
                                          "child": mosaic_assessment(ref_r, alt_r),
                                          "mother_alt_reads": mo[1], "mother_depth": _depth(mo),
                                          "father_alt_reads": fa[1], "father_depth": _depth(fa)})
    finally:
        records.close()
    seconds_scan = round(time.monotonic() - t0, 2)
    read_note = records.read_note()
    if read_note:
        bad_fraction = records.malformed / max(1, records.data_lines)
        if bad_fraction > 0.5 or (bad_fraction > 0.1 and records.malformed > 20):
            # the triage rule: no QC verdict ("no Mendelian errors", "no ROH") from a file that could not be read
            raise UsageError(f"{path}: {read_note}. That is {bad_fraction:.0%} of the data lines; QC refuses to "
                             "report on a file it cannot read (check it is tab-separated and complete)")
        warnings.append(read_note)
    if n_snv == 0:
        raise UsageError(f"{path}: no record is a usable biallelic autosomal SNV (FILTER PASS or '.', a GT field, "
                         f"{n_rec:,} record(s) read){'; the scan stopped by ' + str(stopped_by) if stopped_by else ''}: "
                         "QC has nothing to assess and reports no verdict (contig names zebra reads: 1-22 or "
                         "chr1-chr22)")
    if not complete:
        warnings.append(f"QC scan stopped by {stopped_by} after {n_rec:,} records: every count below is partial, and "
                        "a verdict resting on few chromosomes is marked inconclusive")
    if no_depth:
        notes.append(f"{no_depth:,} genotype(s) carried no DP/AD: depth was not checked for them")
    if bad_pos:
        warnings.append(f"{bad_pos} record(s) with a position beyond any human chromosome were skipped")
    _progress(progress, f"scan done: {n_rec:,} records, {n_snv:,} biallelic autosomal SNVs ({seconds_scan} s)")

    # ---- sex (a parent's role states a sex too: a swapped mother/father is caught here)
    meta = {"records_scanned": n_rec, "xy_records": None, "scan_complete": complete, "stopped_by": stopped_by}
    role_sex: Dict[str, str] = {}
    for c_, m_, f_ in trios:
        role_sex.setdefault(m_, "female")
        role_sex.setdefault(f_, "male")
    sex_out: Dict[str, Any] = {}
    for s in chosen:
        call = V.sex_from_counts(sexc[s], meta)
        stated, basis = None, None
        if pedigree is not None and pedigree.get(s) and pedigree.get(s).sex:  # type: ignore[union-attr]
            stated, basis = pedigree.get(s).sex, "PED"  # type: ignore[union-attr]
        if s == proband and sex:
            stated, basis = sex, "--sex"
        if stated is None and s in role_sex:
            stated, basis = role_sex[s], "named as a " + ("mother" if role_sex[s] == "female" else "father")
        agrees = None if (call["sex"] is None or stated is None) else call["sex"] == stated
        sex_out[s] = {"inferred": call["sex"], "stated": stated, "stated_by": basis, "agrees": agrees,
                      "basis": call["basis"]}
        if agrees is False:
            warnings.append(f"sex mismatch for {s}: stated {stated} ({basis}), genotypes look {call['sex']} "
                            f"({call['basis']}): a sample swap, swapped parent labels, or a PED/sex error")
    doubtful = {t for t in trios if any(sex_out.get(x, {}).get("agrees") is False and x in role_sex
                                        for x in (t[1], t[2]))}

    # ---- ROH
    roh_out: Dict[str, Any] = {}
    roh_segs: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
    roh_spans: Dict[str, Dict[str, Tuple[int, int]]] = {}
    for s in roh_samples:
        by_chrom: Dict[str, List[Dict[str, Any]]] = {}
        spans: Dict[str, Tuple[int, int]] = {}
        for c in AUTOSOMES:
            pos = roh_pos[s].get(c)
            if not pos:
                continue
            het = roh_het[s][c]
            if any(pos[i] > pos[i + 1] for i in range(len(pos) - 1)):
                order = sorted(range(len(pos)), key=lambda i: pos[i])
                pos = array("q", (pos[i] for i in order))
                het = bytearray(het[i] for i in order)
                roh_pos[s][c], roh_het[s][c] = pos, het
            spans[c] = (pos[0], pos[-1])
            by_chrom[c] = find_roh(pos, het)
        roh_segs[s] = by_chrom
        summ = roh_summary(by_chrom, spans, assembly)
        summ["informative_sites"] = sum(len(v) for v in roh_pos[s].values())
        if summ["informative_sites"] < ROH_MIN_INFORMATIVE:
            summ["assessed"] = False
            summ["reading"] = (f"not assessed: only {summ['informative_sites']} informative (ALT-carrying, "
                               f"quality-passing) SNVs, fewer than {ROH_MIN_INFORMATIVE}")
            summ["upd_hint"] = None
        roh_out[s] = summ
        roh_spans[s] = spans

    # ---- Mendelian errors and UPD
    me_out = []
    for t in trios:
        child, mom, dad = t
        per = me[t]
        tot: Counter = Counter()
        for c in AUTOSOMES:
            tot.update(per.get(c, Counter()))
        sites = tot["sites"]
        errors = tot["de_novo_like"] + tot["pat_absent"] + tot["mat_absent"] + tot["other"]
        rows = []
        flags = []
        for c in AUTOSOMES:
            b = per.get(c)
            if not b or not b["sites"]:
                continue
            rows.append({"chrom": c, "sites": b["sites"], "pat_absent": b["pat_absent"], "mat_absent": b["mat_absent"],
                         "de_novo_like": b["de_novo_like"],
                         "error_rate": round((b["sites"] - b["ok"]) / b["sites"], 4)})
        # a whole-genome excess against one parent is parentage or a swap, never UPD: decided first, by the
        # asymmetry between the two one-sided classes (robust to how many uninformative sites a cohort VCF adds)
        genome: List[str] = []
        swapped: List[str] = []
        for cls, mirror, who in (("pat_absent", "mat_absent", "father"), ("mat_absent", "pat_absent", "mother")):
            n_cls, n_mir = tot[cls], tot[mirror]
            if n_cls < 20 or not sites:
                continue
            p_asym = 1.0 - binom_cdf(n_cls - 1, n_cls + n_mir, 0.5)
            base = (n_mir + 1) / sites
            elevated = sum(1 for r in rows if r[cls] >= 3 and r[cls] / r["sites"] >= 5 * base)
            if p_asym < 1e-6 and elevated >= 10:
                swapped.append(cls)
                genome.append(f"{n_cls} of {sites} informative SNVs lack the {who}'s allele (against {n_mir} lacking "
                              f"the other parent's), elevated on {elevated} chromosomes: {who} "
                              f"{dad if who == 'father' else mom} is probably not the biological {who} of {child} "
                              "(sample swap or misattributed parentage); this is not UPD")
        if not swapped and sites:
            # both one-sided classes high and balanced: the asymmetry test cannot see it, the rates can
            both = [cls for cls in ("pat_absent", "mat_absent")
                    if tot[cls] >= 20 and tot[cls] / sites >= 0.01
                    and sum(1 for r in rows if r["sites"] >= 20 and r[cls] / r["sites"] >= 0.01) >= 10]
            if len(both) == 2:
                swapped = both
        if len(swapped) == 2:
            genome = [f"{child} lacks alleles of both its named parents across the genome ({tot['pat_absent']} and "
                      f"{tot['mat_absent']} one-sided errors in {sites} informative SNVs): the child's sample is "
                      "swapped or mislabelled (the parents may well be right); this is not UPD"]
        warnings.extend(genome)
        for c in AUTOSOMES:
            b = per.get(c)
            if not b or not b["sites"]:
                continue
            for cls, origin in (("pat_absent", "maternal"), ("mat_absent", "paternal")):
                if cls in swapped:
                    continue
                e = b[cls]
                bg_e, bg_n = tot[cls] - e, sites - b["sites"]
                rate = max(bg_e, 1) / max(bg_n, 1)
                lam = rate * b["sites"]
                p = poisson_sf(e, lam)
                if not (e >= UPD_MIN_ERRORS and p < UPD_MAX_P and e / b["sites"] >= UPD_MIN_FOLD * rate):
                    continue
                segs_c = roh_segs.get(child, {}).get(c, [])
                roh_bp = sum(sg["end"] - sg["start"] + 1 for sg in segs_c)
                pos = roh_pos.get(child, {}).get(c)
                span = (pos[-1] - pos[0] + 1) if pos else 0
                cover = roh_bp / span if span else 0.0
                errs = me_pos[t].get((c, cls), [])
                hit_segs = [sg for sg in segs_c if any(sg["start"] <= x <= sg["end"] for x in errs)]
                inside = sum(1 for x in errs if any(sg["start"] <= x <= sg["end"] for sg in segs_c))
                in_frac = inside / len(errs) if errs else 0.0
                region = f"chr{c}:{min(errs)}-{max(errs)}" if errs else None
                if cover >= 0.8:
                    kind = "isodisomy"
                elif in_frac >= 0.8 and hit_segs:
                    kind = "segmental (a deletion of the other parent's allele, or segmental isodisomy)"
                    region = f"chr{c}:{min(sg['start'] for sg in hit_segs)}-{max(sg['end'] for sg in hit_segs)}"
                elif cover < 0.1:
                    kind = "heterodisomy"
                else:
                    kind = "mixed (segmental isodisomy within heterodisomy, or partial UPD)"
                shown_origin = origin
                disorder = IMPRINTED_UPD.get((c, origin))
                if t in doubtful:
                    shown_origin = "unknown"
                    disorder = None
                flag = {"chrom": c, "origin": shown_origin, "type": kind, "one_sided_errors": e, "sites": b["sites"],
                        "background_rate": float(f"{rate:.3g}"), "poisson_p": float(f"{p:.3g}"),
                        "child_roh_fraction": round(cover, 3), "errors_inside_child_roh": f"{inside}/{len(errs)}",
                        "imprinting_disorder": disorder, "region": region}
                if t in doubtful:
                    flag["origin_withheld"] = ("a parent's genotypes contradict the sex of its role (see sex): the "
                                               "labels may be swapped, so the parent of origin is not stated")
                flags.append(flag)
                lacking = "father" if origin == "maternal" else "mother"
                if t in doubtful:
                    head = (f"possible UPD or deletion on chr{c} in {child} ({kind}): {e} of {b['sites']} informative "
                            "SNVs lack one named parent's allele; the parent of origin is withheld because the "
                            "parents' sexes contradict their labels")
                elif kind.startswith("segmental"):
                    head = (f"possible loss of the {lacking}'s allele over {region} in {child} ({e} one-sided errors, "
                            f"{in_frac:.0%} inside child ROH): a deletion on the {lacking}'s chromosome, or segmental "
                            f"{origin} isodisomy")
                else:
                    head = (f"possible {origin} UPD of chr{c} in {child} ({kind}): {e} of {b['sites']} informative "
                            f"SNVs lack the {lacking}'s allele (background {rate:.2%}, Poisson p {p:.2g}); child ROH "
                            f"covers {cover:.0%} of the chromosome")
                warnings.append(head + (f". UPD({c}){origin[:3]} causes {disorder}" if disorder else "")
                                + ". Confirm with CNV data, SNP array or methylation")
        # per-chromosome rows for the proband's trio and any trio with a finding (the envelope has a size budget)
        main_trio = t == trios[0]  # the proband's trio sorts first
        detail = main_trio or bool(flags) or bool(genome) or len(trios) <= 2
        me_out.append({
            "trio": {"child": child, "mother": mom, "father": dad}, "sites": sites,
            "sites_rule": "informative: called in all three with an ALT in at least one",
            "consistent": tot["ok"], "errors": errors,
            "error_rate": round(errors / sites, 4) if sites else None,
            "by_class": {k: tot[k] for k in ("de_novo_like", "pat_absent", "mat_absent", "other")},
            "per_chromosome": ([[r["chrom"], r["sites"], r["pat_absent"], r["mat_absent"], r["de_novo_like"]]
                                for r in rows] if detail else None),
            "per_chromosome_columns": ["chrom", "informative_sites", "pat_absent", "mat_absent", "de_novo_like"],
            "upd": flags, "genome_wide": genome or None,
        })
        if sites and errors / sites > 0.05 and not genome:
            warnings.append(f"trio {child}/{mom}/{dad}: Mendelian error rate {errors / sites:.1%} of informative SNVs "
                            "(> 5%): check sample identity, joint calling and the depth/GQ gates")

    # ---- relatedness (after UPD: a uniparental chromosome is left out of that child's parent estimate)
    upd_skip: Dict[frozenset, List[str]] = defaultdict(list)
    for m in me_out:
        t = m["trio"]
        for f in m["upd"]:
            if f["origin"] == "unknown":
                continue
            other = t["father"] if f["origin"] == "maternal" else t["mother"]
            upd_skip[frozenset((t["child"], other))].append(f["chrom"])
    raw = []
    for (a, b), per_chrom in kin.items():
        skip = upd_skip.get(frozenset((a, b)), [])
        tot5 = [0, 0, 0, 0, 0]
        used_chroms = 0
        for c, cnt in per_chrom.items():
            if c in skip:
                continue
            if cnt[0]:
                used_chroms += 1
            for i in range(5):
                tot5[i] += cnt[i]
        n, hh, opp, ha, hb = tot5
        phi = king_robust(hh, opp, ha, hb)
        ibs0 = (opp / min(ha, hb)) if min(ha, hb) > 0 else None
        raw.append((a, b, skip, used_chroms, n, hh, opp, phi, ibs0))
    # the parent-offspring / full-sibling cut, scaled to this file's unrelated pairs when there are any
    unrel = sorted(x[8] for x in raw if x[7] is not None and x[7] < 0.0442 and x[8] is not None and x[4] >= 1000)
    if unrel:
        po_cut = 0.1 * unrel[len(unrel) // 2]
        po_basis = f"0.1 x the median IBS0/het of {len(unrel)} unrelated pair(s) in this file ({unrel[len(unrel) // 2]:.4f})"
    else:
        po_cut = PO_IBS0_PER_HET
        po_basis = f"fixed {PO_IBS0_PER_HET} (no unrelated pair in this file to scale it on)"
    tp = pedigree if pedigree is not None else _trio_ped(trios)
    rel_pairs = []
    unlisted = 0
    for a, b, skip, used_chroms, n, hh, opp, phi, ibs0 in raw:
        expected = pedigree.expected(a, b) if pedigree is not None else None
        if expected is None and pedigree is None and trios:
            for c_, m_, f_ in trios:
                if {a, b} in ({c_, m_}, {c_, f_}):
                    expected = "parent-offspring"
                elif {a, b} == {m_, f_}:
                    expected = "unrelated (parents of one child)"
        kind, problem = _pair_verdict(a, b, phi, ibs0, expected, tp, po_cut)
        if problem and used_chroms < KING_MIN_CHROMS:
            problem = (f"inconclusive — kinship from only {used_chroms} autosome(s) (sibling sharing varies too much "
                       f"on so few): {problem}")
        row = {"a": a, "b": b, "kinship": None if phi is None else round(phi, 4), "relationship": kind,
               "ibs0_per_het": None if ibs0 is None else round(ibs0, 4), "snvs": n, "autosomes": used_chroms,
               "expected": expected, "problem": problem}
        if skip:
            row["excluded_chromosomes"] = [f"chr{c}" for c in skip]
            row["note"] = (f"chr{', chr'.join(skip)} left out: flagged as uniparental in the child, where no allele "
                           "of this parent is expected")
        if problem:
            warnings.append(problem)
        elif expected and n and n < 1000:
            warnings.append(f"{a}–{b}: kinship from only {n} shared SNVs — imprecise")
        # a cohort's unrelated pairs are counted, not listed (the envelope has a size budget)
        if len(raw) > 15 and expected is None and problem is None and degree(phi).startswith(("unrelated", "unknown")):
            unlisted += 1
            continue
        rel_pairs.append(row)
    # parentage problems found by kinship make that trio's de novo list unreliable
    bad_parent: Dict[Tuple[str, str, str], List[str]] = defaultdict(list)
    for row in rel_pairs:
        if row["problem"] and row["expected"] == "parent-offspring" and degree(row["kinship"]).startswith("unrelated"):
            for t in trios:
                if {row["a"], row["b"]} in ({t[0], t[1]}, {t[0], t[2]}):
                    bad_parent[t].append(t[1] if t[1] in (row["a"], row["b"]) else t[2])

    # ---- ROH readings, with any chromosome called uniparental isodisomy taken out of the relatedness reading
    iso: Dict[str, List[str]] = defaultdict(list)
    for m in me_out:
        for f in m["upd"]:
            if f["type"] == "isodisomy" or f["type"].startswith("mixed"):
                iso[m["trio"]["child"]].append(f["chrom"])
    for s, summ in roh_out.items():
        if summ.get("assessed") is False:
            continue
        if iso.get(s):
            keep = {c: v for c, v in roh_segs[s].items() if c not in iso[s]}
            excl = roh_summary(keep, {c: v for c, v in roh_spans[s].items() if c not in iso[s]}, assembly)
            summ["excluding_upd"] = {"chromosomes": [f"chr{c}" for c in iso[s]], "total_mb": excl["total_mb"],
                                     "f_roh": excl["f_roh"], "reading": excl["reading"]}
            summ["reading"] = (f"{excl['reading']} (F_ROH {excl['f_roh']:.2%} once chr{', chr'.join(iso[s])}, "
                               "called uniparental isodisomy, is left out)")
            summ["upd_hint"] = None
            if excl["f_roh"] >= 0.02:
                warnings.append(f"{s}: {excl['total_mb']} Mb in runs of homozygosity outside the uniparental "
                                f"chromosome(s) (F_ROH {excl['f_roh']:.1%}): {excl['reading']}")
            continue
        if summ["upd_hint"]:
            warnings.append(f"{s}: {summ['upd_hint']}")
        elif summ["f_roh"] >= 0.02:
            warnings.append(f"{s}: {summ['total_mb']} Mb in runs of homozygosity (F_ROH {summ['f_roh']:.1%}): "
                            f"{summ['reading']} — recessive candidates inside these runs gain weight")

    # ---- de novo mosaicism
    dn_out = []
    for t in trios:
        calls = denovo[t]
        for r in calls:
            r["flags"] = []
            if (r["child"] or {}).get("possible_mosaic"):
                ch = r["child"]
                r["flags"].append(f"ALT fraction {ch['vaf']} (95% CI {ch['ci95'][0]}–{ch['ci95'][1]}) below a germline "
                                  "het: possible postzygotic mosaic (or artefact)")
            for role, k, dk in (("mother", "mother_alt_reads", "mother_depth"), ("father", "father_alt_reads",
                                                                                 "father_depth")):
                alt_r, depth = r[k] or 0, r[dk] or 0
                if alt_r >= 2 and depth and alt_r / depth >= PARENT_MOSAIC_MIN_FRACTION:
                    r["flags"].append(f"{role} has {alt_r} of {depth} reads with the ALT: possible parental "
                                      "mosaicism (recurrence risk)")
        flagged = sum(1 for r in calls if (r["child"] or {}).get("possible_mosaic"))
        parental = sum(1 for r in calls if any("parental mosaicism" in f for f in r["flags"]))
        ordered = sorted(calls, key=lambda r: (not r["flags"], -(r["child"] or {}).get("depth", 0)))
        keep = 20 if (t == trios[0] or len(trios) <= 2) else 3
        entry = {"trio": {"child": t[0], "mother": t[1], "father": t[2]}, "de_novo_like_hets": denovo_n[t],
                 "listed": min(keep, len(ordered)), "possible_mosaic": flagged,
                 "possible_parental_mosaic": parental, "calls": ordered[:keep]}
        if bad_parent.get(t):
            entry["unreliable"] = (f"{', '.join(bad_parent[t])} looks unrelated to {t[0]} (kinship): most of these "
                                   "'de novo' calls are alleles from the biological parent")
            warnings.append(f"de novo-like calls in {t[0]} are unreliable: {entry['unreliable']}")
        dn_out.append(entry)
    result = {
        "vcf": str(Path(os.path.expanduser(path)).resolve()), "assembly": assembly, "samples": chosen,
        "pedigree": ({"path": pedigree.path, "columns": ["id", "father", "mother", "sex", "affected"],
                      "individuals": [[x.iid, x.father, x.mother, x.sex, x.affected] for x in pedigree.people.values()]}
                     if pedigree is not None else None),
        "trios": [{"child": c, "mother": m, "father": f} for c, m, f in trios],
        "scan": {"records": n_rec, "biallelic_autosomal_snvs": n_snv, "king_sites": sum(king_by_chrom.values()),
                 "king_autosomes": len(king_by_chrom), "complete": complete,
                 "stopped_by": stopped_by, "seconds": round(time.monotonic() - t0, 2),
                 "gates": {"min_dp": min_dp, "min_gq": min_gq, "filter": "PASS or ."}},
        "sex": sex_out, "relatedness": rel_pairs, "relatedness_unlisted_unrelated_pairs": unlisted,
        "po_fs_cut": {"ibs0_per_het": round(po_cut, 5), "basis": po_basis},
        "roh": roh_out, "mendelian": me_out, "de_novo": dn_out,
        "per_sample": {s: dict(c) for s, c in per_sample.items()},
        "methods": METHODS, "notes": notes,
    }
    from zebra.http import source_record

    src = source_record("local VCF", os.path.basename(path), url=None, note="QC computed on this machine; nothing sent")
    query = {"vcf": path, "ped": ped, "proband": proband, "mother": mother, "father": father, "sex": sex,
             "assembly": assembly, "min_dp": min_dp, "min_gq": min_gq, "max_samples": max_samples,
             "max_seconds": max_seconds}
    return Outcome(result, sources=[src], warnings=warnings, text=render(result, warnings), query=query)


def _depth(ad: Tuple[Optional[int], Optional[int]]) -> Optional[int]:
    return (ad[0] or 0) + (ad[1] or 0) if ad[0] is not None and ad[1] is not None else None


def _trio_ped(trios: Sequence[Tuple[str, str, str]]) -> Optional[Pedigree]:
    """A minimal pedigree from --proband/--mother/--father, so parent roles can be named."""
    if not trios:
        return None
    people: "OrderedDict[str, Person]" = OrderedDict()
    for c, m, f in trios:
        people.setdefault(m, Person("fam", m, None, None, "female", None))
        people.setdefault(f, Person("fam", f, None, None, "male", None))
        people[c] = Person("fam", c, f, m, None, None)
    return Pedigree(people)


def render(r: Dict[str, Any], warnings: Sequence[str] = ()) -> str:
    sc = r["scan"]
    lines = [f"zebra qc — {', '.join(r['samples'])} ({r['assembly'] or 'build not stated'}); "
             f"{sc['records']:,} records, {sc['biallelic_autosomal_snvs']:,} biallelic autosomal SNVs"
             + ("" if sc["complete"] else f" [scan stopped: {sc['stopped_by']}]") + f", {sc['seconds']} s"]
    lines.append("sex:")
    for s, v in r["sex"].items():
        st = f" (stated {v['stated']}{', AGREES' if v['agrees'] else (', MISMATCH' if v['agrees'] is False else '')})" \
            if v["stated"] else ""
        lines.append(f"  {s}: {v['inferred'] or 'undetermined'}{st}")
    if r["relatedness"]:
        lines.append("relatedness (KING-robust kinship; IBS0 per het):")
        for p in r["relatedness"]:
            k = "-" if p["kinship"] is None else f"{p['kinship']:.3f}"
            i0 = "-" if p["ibs0_per_het"] is None else f"{p['ibs0_per_het']:.4f}"
            lines.append(f"  {p['a']}–{p['b']}: {k} {p['relationship']} (IBS0/het {i0}, {p['snvs']:,} SNVs on "
                         f"{p['autosomes']} autosomes)"
                         + (f" [expected {p['expected']}]" if p["expected"] else "")
                         + (f"  ⚠ {p['problem']}" if p["problem"] else ""))
    for s, v in r["roh"].items():
        lines.append(f"ROH {s}: {v['n_segments']} run(s) ≥1 Mb, {v['total_mb']} Mb, F_ROH {v['f_roh']:.2%} — "
                     f"{v['reading']}")
        for seg in v["segments"][:6]:
            lines.append(f"  chr{seg['chrom']}:{seg['start']:,}-{seg['end']:,} {seg['mb']} Mb ({seg['sites']} sites, "
                         f"{seg['hets']} het)")
        if v["upd_hint"]:
            lines.append(f"  ⚠ {v['upd_hint']}")
    for m in r["mendelian"]:
        t = m["trio"]
        rate = "-" if m["error_rate"] is None else f"{m['error_rate']:.2%}"
        bc = m["by_class"]
        lines.append(f"Mendelian {t['child']} (mother {t['mother']}, father {t['father']}): {m['errors']} errors in "
                     f"{m['sites']:,} informative SNVs ({rate}): de novo-like {bc['de_novo_like']}, paternal allele "
                     f"absent {bc['pat_absent']}, maternal allele absent {bc['mat_absent']}")
        for f in m["upd"]:
            lines.append(f"  ⚠ UPD({f['chrom']}){f['origin'][:3]} {f['type']}: {f['one_sided_errors']}/{f['sites']} "
                         f"(p {f['poisson_p']:.2g}), child ROH {f['child_roh_fraction']:.0%}"
                         + (f" — {f['imprinting_disorder']}" if f["imprinting_disorder"] else ""))
        for g in m["genome_wide"] or []:
            lines.append(f"  ⚠ {g}")
    for d in r["de_novo"]:
        t = d["trio"]
        lines.append(f"de novo-like hets in {t['child']}: {d['de_novo_like_hets']} (possible mosaic "
                     f"{d['possible_mosaic']}, parental ALT reads {d['possible_parental_mosaic']})"
                     + (f"  ⚠ unreliable: {d['unreliable']}" if d.get("unreliable") else ""))
        for c in d["calls"][:8]:
            ch = c["child"] or {}
            vaf = f"VAF {ch.get('vaf')} CI {ch.get('ci95')}" if ch else "no AD"
            lines.append(f"  {c['variant']} {vaf}" + (f" — {'; '.join(c['flags'])}" if c["flags"] else ""))
    lines.extend(f"note: {n}" for n in r["notes"])
    lines.append("methods: " + "; ".join(f"{k}: {v['cite']}" for k, v in r["methods"].items()))
    return "\n".join(lines)
