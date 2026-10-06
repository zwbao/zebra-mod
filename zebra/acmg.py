"""ACMG/AMP variant classification arithmetic.

Two independent readings of the same evidence codes:
  * the naturally scaled point system (Tavtigian et al. 2020, Hum Mutat 41:1734):
    Supporting 1, Moderate 2, Strong 4, Very strong 8 (benign codes negative);
    P >= 10, LP 6..9, VUS 0..5, LB -1..-6, B <= -7;
  * the original combining rules (Richards et al. 2015, Genet Med 17:405).
They usually agree; when they differ the points reading is the one ClinGen's
Bayesian framework supports (Tavtigian et al. 2018), and the difference is said.

This module does arithmetic and structure only. Whether a code applies is a
judgement made from evidence (skills/zebra-variant holds that contract); the
numeric criteria that ClinGen has calibrated (PP3/BP4 for REVEL and SpliceAI,
BA1, BS1, PM2_Supporting) are offered by `suggest` as suggestions that carry
their thresholds and citations.

Every numeric threshold below carries its source. Where a threshold is zebra's
own conservative floor rather than a published number, the comment says so.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

POINTS = {"VeryStrong": 8, "Strong": 4, "Moderate": 2, "Supporting": 1}
STRENGTH_ORDER = {"Supporting": 1, "Moderate": 2, "Strong": 3, "VeryStrong": 4}
STRENGTH_ALIASES = {
    "verystrong": "VeryStrong", "very_strong": "VeryStrong", "vs": "VeryStrong", "vstr": "VeryStrong",
    "strong": "Strong", "str": "Strong", "s": "Strong",
    "moderate": "Moderate", "mod": "Moderate", "m": "Moderate",
    "supporting": "Supporting", "sup": "Supporting", "p": "Supporting",
    "standalone": "StandAlone", "stand_alone": "StandAlone", "sa": "StandAlone",
}
PATHOGENIC = {
    "PVS1": "VeryStrong",
    **{f"PS{i}": "Strong" for i in range(1, 5)},
    **{f"PM{i}": "Moderate" for i in range(1, 7)},
    **{f"PP{i}": "Supporting" for i in range(1, 6)},
}
BENIGN = {
    "BA1": "StandAlone",
    **{f"BS{i}": "Strong" for i in range(1, 5)},
    **{f"BP{i}": "Supporting" for i in range(1, 8)},
}
DEPRECATED = {
    "PP5": "ClinGen SVI recommends not using PP5 (Biesecker & Harrison 2018): weigh the underlying evidence instead.",
    "BP6": "ClinGen SVI recommends not using BP6 (Biesecker & Harrison 2018): weigh the underlying evidence instead.",
}
CODE_RE = re.compile(r"^(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])(?:[_\-\s]?([A-Za-z_]+))?$", re.I)

# Per-code strength ceilings, each with the source that sets it. A code absent from
# this table is left unrestricted, because ClinGen's published ceilings for those
# (PM3 and PS2/PM6, for instance) do reach VeryStrong.
STRENGTH_CEILINGS: Dict[str, Tuple[str, str]] = {
    # Richards et al. 2015 places PM2 at Moderate; ClinGen SVI (2020) recommends
    # Supporting by default. Moderate is therefore the ceiling, Supporting the default advice.
    "PM2": ("Moderate", "Richards et al. 2015, Genet Med 17:405 places PM2 at Moderate; ClinGen SVI (2020) recommends PM2_Supporting"),
    # Biesecker et al. 2024, AJHG 111:24 Table 3: co-segregation points top out at
    # 4.0 (= Strong on the Tavtigian scale) for the counts the guidance tabulates.
    "PP1": ("Strong", "Biesecker et al. 2024, AJHG 111:24 Table 3: co-segregation evidence tops out at 4 points (Strong)"),
    # Pejaver et al. 2022, AJHG 109:2163 Table 2: no missense predictor reached the
    # very-strong pathogenic tier, so PP3 above Strong has no calibration behind it.
    "PP3": ("Strong", "Pejaver et al. 2022, AJHG 109:2163 Table 2: no missense predictor reaches the very-strong pathogenic tier"),
    # Same table, benign side: REVEL <= 0.003 does reach very strong, so BP4 may.
    "BP4": ("VeryStrong", "Pejaver et al. 2022, AJHG 109:2163 Table 2: REVEL <= 0.003 reaches the very-strong benign tier"),
}

# Pejaver et al. 2022 (AJHG 109:2163), REVEL calibrated thresholds (Table 2).
REVEL_PP3 = [(0.932, "Strong"), (0.773, "Moderate"), (0.644, "Supporting")]
REVEL_BP4 = [(0.003, "VeryStrong"), (0.016, "Strong"), (0.183, "Moderate"), (0.290, "Supporting")]
# Cheng et al. 2023, Science 381:eadg7492: AlphaMissense class boundaries, "likely benign"
# below 0.34 and "likely pathogenic" above 0.564. Only the benign boundary is used here,
# and only to answer Walker 2023's question "has protein impact been excluded?" --
# AlphaMissense has no ClinGen PP3/BP4 strength calibration in this module.
ALPHAMISSENSE_BENIGN = 0.34
# Walker et al. 2023 (AJHG 110:1046), ClinGen SVI Splicing Subgroup.
SPLICEAI_PP3 = 0.2
SPLICEAI_BP4 = 0.1
# Walker 2023 leaves 0.1 < delta < 0.2 uninformative: neither PP3 nor BP4.
SPLICEAI_UNINFORMATIVE = (SPLICEAI_BP4, SPLICEAI_PP3)
# ClinGen SVI (Ghosh et al. 2018, Hum Mutat 39:1525): BA1 at AF > 5% in a continental
# group surveyed with at least 2,000 alleles.
BA1_AF = 0.05
BA1_MIN_AN = 2000
# zebra's own floor, not a published PM2 threshold: absence is only evidence if enough
# alleles were surveyed, and 2,000 is borrowed from Ghosh 2018's BA1 allele-number
# requirement for the same reason (a small group carries no information either way).
PM2_MIN_AN = 2000
# 0.07% is the recessive PM2 cut several ClinGen VCEPs use; it is a convention, not a
# calibrated threshold, hence the instruction to check it against the maximum credible AF.
RECESSIVE_PM2_AF = 0.0007

# Walker et al. 2023 (AJHG 110:1046), BP7 scope: synonymous variants outside the first base and the
# last 3 bases of the exon, and intronic variants "at or beyond positions +7/-21", each "having also
# met BP4". VEP's splice_region_variant marks the first and last 3 exonic bases, a superset of
# Walker's exonic exclusion, so it is used as the (conservative) synonymous gate.
BP7_DONOR_MIN = 7
BP7_ACCEPTOR_MIN = 21
SPLICE_REGION_TERMS = frozenset(("splice_region_variant", "splice_donor_region_variant",
                                 "splice_donor_5th_base_variant", "splice_polypyrimidine_tract_variant",
                                 "splice_donor_variant", "splice_acceptor_variant"))
EXONIC_TERMS = frozenset(("synonymous_variant", "missense_variant", "stop_gained", "frameshift_variant",
                          "inframe_insertion", "inframe_deletion", "coding_sequence_variant",
                          "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant",
                          "start_lost", "stop_lost", "stop_retained_variant", "protein_altering_variant"))
C_LOCATION_RE = re.compile(r"c\.([-*]?\d+(?:[+-]\d+)?(?:_[-*]?\d+(?:[+-]\d+)?)?)")
C_COORD_RE = re.compile(r"^([-*]?\d+)(?:([+-])(\d+))?$")
BP7_RULE = ("Walker et al. 2023 (AJHG 110:1046): \"BP7 is applicable after assignment of BP4 for no adverse "
            "splicing predictions\" -- synonymous variants outside the first base / last 3 bases of the exon, "
            "intronic variants at or beyond +7/-21. The paper cautions \"against inclusion of evolutionary "
            "conservation in assessment of silent or intronic variants for application of code BP7 without "
            "empirically derived justification\", so conservation is not checked.")
MITO_FREQUENCY_NOTE = (
    "BA1 / BS1 / PM2: not assessed for an mtDNA variant. Its population data are homoplasmic and heteroplasmic "
    "frequencies (gnomAD mtDNA, MITOMAP/MITOMASTER), not a nuclear allele count, and the nuclear cut-offs used "
    "here were not calibrated on them; the ClinGen/MSeqDR specification for mtDNA (McCormick et al. 2020, "
    "Hum Mutat 41:2028, PMID 32906214) sets its own frequency criteria. Read the homoplasmic/heteroplasmic "
    "frequencies on the variant card against that specification."
)

# The PVS1 classes (Abou Tayoun et al. 2018). stop_lost is not one of them: a protein extended at its
# C-terminus is PM4's question (protein length change), not loss of function.
NULL_CONSEQUENCES = ("stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant",
                     "start_lost", "transcript_ablation")
# Variant classes where splicing is the mechanism under consideration, so a splicing
# predictor can carry a benign code on its own (Walker et al. 2023).
SPLICE_MECHANISM_CONSEQUENCES = (
    "synonymous_variant", "intron_variant", "splice_region_variant",
    "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
    "splice_donor_5th_base_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "coding_sequence_variant",
)
# PP3 and BP4 are two readings of the same computational evidence; so is BP7.
COMPUTATIONAL_CODES = ("PP3", "BP4", "BP7")
PEJAVER_NO_SOLO_CLASS = (
    "Pejaver et al. 2022 (AJHG 109:2163): \"in silico predictors alone are not capable of classifying the "
    "pathogenicity of a variant\" and \"supporting evidence must then be combined with substantial other lines "
    "of evidence\". Computational evidence alone therefore cannot reach Benign or Likely benign here."
)
CLINVAR_PLP_RE = re.compile(r"\b(likely[_ ]pathogenic|pathogenic)\b", re.I)
_NOT_PATHOGENIC_RE = re.compile(r"\b(?:not|non)[-_ ]?pathogenic\b|probable[-_ ]non[-_ ]pathogenic", re.I)


def is_plp(classification: Optional[str]) -> bool:
    """A ClinVar classification that asserts Pathogenic / Likely pathogenic (not conflicting, not 'non-pathogenic')."""
    text = str(classification or "")
    if not text or "conflicting" in text.lower():
        return False
    return bool(CLINVAR_PLP_RE.search(_NOT_PATHOGENIC_RE.sub(" ", text)))
BA1_EXCEPTION_NOTE = (
    "BA1 exception list: Ghosh et al. 2018 (Hum Mutat 39:1525) recognised an initial list of nine high-frequency "
    "variants with evidence of pathogenicity, four of them classified Pathogenic. zebra does not bundle that list, "
    "so it has NOT been checked -- look the variant up before applying BA1 (low-penetrance alleles such as HFE "
    "p.Cys282Tyr and HBB/F5/SERPINA1 risk alleles are the recurring cases)."
)


@dataclass
class Code:
    code: str
    strength: str
    direction: str  # "pathogenic" | "benign"
    points: int
    modified: bool

    def label(self) -> str:
        return self.code if not self.modified else f"{self.code}_{self.strength}"


def parse_code(raw: str) -> Code:
    m = CODE_RE.match(raw.strip())
    if not m:
        raise ValueError(f"{raw!r} is not an ACMG/AMP code (e.g. PVS1, PM2_Supporting, PP3_Strong, BS1)")
    code = m.group(1).upper()
    default = PATHOGENIC.get(code) or BENIGN[code]
    strength = default
    if m.group(2):
        key = m.group(2).lower()
        if key not in STRENGTH_ALIASES:
            raise ValueError(f"{raw!r}: unknown strength {m.group(2)!r} (Supporting, Moderate, Strong, VeryStrong)")
        strength = STRENGTH_ALIASES[key]
    direction = "pathogenic" if code in PATHOGENIC else "benign"
    if code == "BA1":
        if strength != "StandAlone":
            raise ValueError("BA1 is stand-alone; it cannot be modified")
        return Code(code, "StandAlone", "benign", -8, False)
    if strength == "StandAlone":
        raise ValueError(f"{raw!r}: only BA1 is stand-alone")
    ceiling = STRENGTH_CEILINGS.get(code)
    if ceiling and STRENGTH_ORDER[strength] > STRENGTH_ORDER[ceiling[0]]:
        raise ValueError(
            f"{raw!r}: {code} cannot be applied above {ceiling[0]} -- {ceiling[1]}. "
            f"Use {code}_{ceiling[0]} or weaker."
        )
    pts = POINTS[strength]
    return Code(code, strength, direction, pts if direction == "pathogenic" else -pts, strength != default)


def points_class(total: int) -> str:
    if total >= 10:
        return "Pathogenic"
    if total >= 6:
        return "Likely pathogenic"
    if total >= 0:
        return "Uncertain significance"
    if total >= -6:
        return "Likely benign"
    return "Benign"


def richards_class(codes: Sequence[Code], svi_2020: bool = False) -> str:
    """The 2015 combining rules.

    Modified strengths are read as their new level on the pathogenic side
    (Supporting/Moderate/Strong/VeryStrong map onto pp/pm/ps/pvs). Richards 2015
    has no benign Moderate tier, so a benign code at Moderate is read as
    Supporting and one at Strong or VeryStrong as Strong -- the conservative
    mapping, since the rules were written for two benign tiers only.

    `svi_2020=True` adds the one combination ClinGen SVI's PM2 recommendation
    (v1.0, approved 4 Sep 2020) added to the 2015 rules: "the combination of
    one Very Strong criterion and one Supporting criterion reach a
    classification of Likely pathogenic" (PVS1 + PM2_Supporting, D-P2-6).
    """
    pvs = sum(1 for c in codes if c.direction == "pathogenic" and c.strength == "VeryStrong")
    ps = sum(1 for c in codes if c.direction == "pathogenic" and c.strength == "Strong")
    pm = sum(1 for c in codes if c.direction == "pathogenic" and c.strength == "Moderate")
    pp = sum(1 for c in codes if c.direction == "pathogenic" and c.strength == "Supporting")
    ba = any(c.code == "BA1" for c in codes)
    bs = sum(1 for c in codes if c.direction == "benign" and c.strength in ("Strong", "VeryStrong"))
    bp = sum(1 for c in codes if c.direction == "benign" and c.strength in ("Supporting", "Moderate"))

    pathogenic = (
        (pvs >= 1 and (ps >= 1 or pm >= 2 or (pm == 1 and pp == 1) or pp >= 2))
        or pvs >= 2
        or ps >= 2
        or (ps == 1 and (pm >= 3 or (pm == 2 and pp >= 2) or (pm == 1 and pp >= 4)))
    )
    likely_pathogenic = (
        (pvs >= 1 and pm == 1)
        or (ps == 1 and 1 <= pm <= 2)
        or (ps == 1 and pp >= 2)
        or pm >= 3
        or (pm == 2 and pp >= 2)
        or (pm == 1 and pp >= 4)
        or (svi_2020 and pvs >= 1 and pp >= 1)
    )
    benign = ba or bs >= 2
    likely_benign = (bs == 1 and bp >= 1) or bp >= 2

    any_p = pathogenic or likely_pathogenic
    any_b = benign or likely_benign
    if any_p and any_b:
        return "Uncertain significance"
    if pathogenic:
        return "Pathogenic"
    if likely_pathogenic:
        return "Likely pathogenic"
    if benign:
        return "Benign"
    if likely_benign:
        return "Likely benign"
    return "Uncertain significance"


def _tier_options(code: Code) -> List[Tuple[int, Optional[Code]]]:
    """The code at its own strength, at every weaker tier, and dropped (0 points)."""
    sign = 1 if code.direction == "pathogenic" else -1
    out: List[Tuple[int, Optional[Code]]] = [(abs(code.points), code)]
    for name, pts in sorted(POINTS.items(), key=lambda kv: -kv[1]):
        if pts < abs(code.points):
            out.append((pts, Code(code.code, name, code.direction, sign * pts, True)))
    out.append((0, None))
    return out


def cap_pm1_pp3(pm1: Code, pp3: Code, cap: int = POINTS["Strong"]) -> Tuple[Optional[Code], Optional[Code]]:
    """PM1 and PP3 weakened so that their sum stays within `cap` points (Pejaver et al. 2022).

    D-P0-3: the old rule weakened PP3 to exactly `points - excess` and deleted it
    when no tier had that many points, so PP3_Strong + PM1_Supporting (5 points,
    excess 1, target 3) lost PP3 altogether and the pair counted 1 point. Here
    every combination of tiers at or below each code's own strength is
    considered, the pair keeps the LARGEST total that does not exceed the cap,
    and among equal totals PM1 keeps its strength (PP3 is weakened first, as
    before). PP3_Strong + PM1_Supporting therefore becomes PP3_Strong alone
    (4 points); PP3_Strong + PM1 becomes PP3_Moderate + PM1 (4 points).
    """
    best: Optional[Tuple[Tuple[int, int, int], Optional[Code], Optional[Code]]] = None
    for p_pts, p_code in _tier_options(pp3):
        for m_pts, m_code in _tier_options(pm1):
            total = p_pts + m_pts
            if total > cap:
                continue
            key = (total, m_pts, p_pts)
            if best is None or key > best[0]:
                best = (key, m_code, p_code)
    assert best is not None  # (0, 0) always fits
    return best[1], best[2]


def classify(raw_codes: List[str]) -> Dict[str, Any]:
    codes: List[Code] = []
    warnings: List[str] = []
    seen: Dict[str, str] = {}
    for raw in raw_codes:
        c = parse_code(raw)
        if c.code in seen:
            warnings.append(f"{c.code} given twice ({seen[c.code]}, {raw}); each code counts once — kept the first")
            continue
        seen[c.code] = raw
        if c.code in DEPRECATED:
            warnings.append(DEPRECATED[c.code])
        codes.append(c)

    names = {c.code for c in codes}
    if "PVS1" in names and "PP3" in names:
        warnings.append("PVS1 with PP3: SVI guidance is not to add PP3 for the same predicted loss-of-function mechanism (double counting).")
    if "PM2" in names:
        pm2 = next(c for c in codes if c.code == "PM2")
        if pm2.strength != "Supporting":
            warnings.append("PM2 at Moderate: ClinGen SVI (2020) recommends PM2_Supporting by default.")
    # PP3 and BP4 (and BP7) are opposite readings of the same computational evidence:
    # Richards et al. 2015 defines BP4 as "multiple lines of computational evidence
    # suggest no impact", which cannot hold at the same time as PP3.
    for a, b in (("PP3", "BP4"), ("PP3", "BP7")):
        if a in names and b in names:
            warnings.append(
                f"{a} with {b}: these are opposite readings of the same computational evidence and cannot both "
                f"apply — resolve which axis (protein or splicing) the prediction speaks to and keep one."
            )
    if "BP4" in names and "BP7" in names:
        # Walker et al. 2023 (AJHG 110:1046): "BP7 is applicable after assignment of BP4 for no adverse
        # splicing predictions" -- for synonymous variants outside the first base and last 3 bases of the
        # exon, and intronic variants at or beyond +7/-21. Outside that scope the pair double counts.
        warnings.append(
            "BP4 with BP7: Walker et al. 2023 applies BP7 after BP4 only for a synonymous variant outside the "
            "first base / last 3 bases of the exon, or an intronic variant at or beyond +7/-21, both with no "
            "predicted splice effect; for any other variant the two double count one prediction — keep one."
        )
    if "BA1" in names and any(c.direction == "pathogenic" for c in codes):
        warnings.append("BA1 with pathogenic codes: BA1 is stand-alone; recheck the frequency data and the BA1 exception list (Ghosh et al. 2018) — zebra does not bundle that list.")

    # Pejaver et al. 2022: "we recommend that laboratories limit the sum of the evidence
    # strength of PP3 and PM1 to strong." Applied to the points total AND, below, to the
    # code list the 2015 combining rules read, so the two columns cannot disagree about it.
    capped: List[Code] = list(codes)
    cap_applied = False
    if "PM1" in names and "PP3" in names:
        pm1 = next(c for c in codes if c.code == "PM1")
        pp3 = next(c for c in codes if c.code == "PP3")
        if pm1.points + pp3.points > POINTS["Strong"]:
            cap_applied = True
            pm1_c, pp3_c = cap_pm1_pp3(pm1, pp3)
            kept = [c for c in (pm1_c, pp3_c) if c is not None]
            capped = [c for c in codes if c.code not in ("PM1", "PP3")] + kept
            became = " + ".join(c.label() for c in kept) or "nothing"
            dropped = [c.code for c, new in ((pm1, pm1_c), (pp3, pp3_c)) if new is None]
            warnings.append(
                f"PM1 + PP3 exceed Strong together ({pm1.label()} + {pp3.label()} = {pm1.points + pp3.points} points); "
                f"Pejaver et al. 2022 cap their sum at 4 points (Strong). Counted as {became} "
                f"({sum(c.points for c in kept)} points)"
                + (f"; {', '.join(dropped)} not counted" if dropped else "")
                + ". Capped in both readings; `codes_after_cap` lists what was counted."
            )

    total = sum(c.points for c in capped if c.code != "BA1")
    computational_only = bool(codes) and all(c.code in COMPUTATIONAL_CODES for c in codes)
    if "BA1" in names:
        point_class = "Benign"
        total_display = None
    else:
        point_class = points_class(total)
        total_display = total
    rules = richards_class(capped)
    rules_svi = richards_class(capped, svi_2020=True)

    # Pejaver et al. 2022: calibrated predictor strengths are contributions, not
    # classifications. Cap the benign side at Uncertain when nothing but computational
    # evidence is present. The pathogenic side needs no cap: PP3 is ceilinged at Strong
    # (4 points), which cannot reach Likely pathogenic (6) on its own.
    computational_cap = None
    if computational_only and point_class in ("Benign", "Likely benign"):
        computational_cap = point_class
        point_class = "Uncertain significance"
        warnings.append(
            f"computational evidence only ({', '.join(sorted(names))}): the points total ({total}) would read "
            f"{computational_cap}, capped to Uncertain significance. {PEJAVER_NO_SOLO_CLASS}"
        )
    if computational_only and rules in ("Benign", "Likely benign"):
        rules = "Uncertain significance"
    if computational_only and rules_svi in ("Benign", "Likely benign"):
        rules_svi = "Uncertain significance"

    def _rows(seq: Sequence[Code]) -> List[Dict[str, Any]]:
        return [{"code": c.code, "label": c.label(), "strength": c.strength, "direction": c.direction,
                 "points": c.points} for c in seq]

    out: Dict[str, Any] = {
        "codes": _rows(codes),
        # what the points total and both rule readings actually counted (differs from
        # `codes` only when the PM1 + PP3 cap weakened or dropped one of them)
        "codes_after_cap": _rows(capped),
        "points": total_display,
        "classification": point_class,
        "classification_richards_2015": rules,
        # the 2015 rules plus the one combination ClinGen SVI added in 2020
        # (PVS1 + one Supporting -> Likely pathogenic)
        "classification_richards_2015_svi_2020": rules_svi,
        "agree": point_class == rules,
        "scale": "P >= 10, LP 6..9, VUS 0..5, LB -1..-6, B <= -7 (Tavtigian 2020)",
        "computational_only": computational_only,
        "capped_at_uncertain": computational_cap,
        "pm1_pp3_cap_applied": cap_applied,
        "warnings": warnings,
    }
    if not out["agree"]:
        why = ""
        if len(codes) == 1:
            why = (f" The points reading here rests on the single criterion {codes[0].label()}; the 2015 rules "
                   f"require a combination, which is why they do not reach the same class. Do not report the points "
                   f"class on its own.")
        elif cap_applied:
            why = " The PM1 + PP3 cap (Pejaver 2022) is applied to both readings, so it is not the cause."
        if rules_svi != rules:
            why += (f" ClinGen SVI's PM2 recommendation (2020) added the combination one Very Strong + one "
                    f"Supporting -> Likely pathogenic to the 2015 rules, so read with that addition the rules give "
                    f"{rules_svi} (classification_richards_2015_svi_2020).")
        out["note"] = (
            f"Points give {point_class}; the 2015 combining rules give {rules}. ClinGen's Bayesian framework "
            f"(Tavtigian 2018/2020) supports the points reading; report both if the difference matters.{why}"
        )
    return out


def _band(value: float, bands: List[Tuple[float, str]], higher_is_stronger: bool) -> Optional[str]:
    for threshold, strength in bands:
        if (higher_is_stronger and value >= threshold) or (not higher_is_stronger and value <= threshold):
            return strength
    return None


def _code_label(base: str, strength: str) -> str:
    return base if strength == "Supporting" else f"{base}_{strength}"


def _intron_offset(hgvs_c: Optional[str]) -> Optional[Tuple[str, int]]:
    """The intronic offset nearest the exon: ('+', 15) for c.123+15A>G, ('-', 5) for c.124-30_124-5del.

    Every coordinate of a range is read (the second one carries no `c.`), so a deletion that reaches
    into the splice region is judged by its exon-side end. None when any coordinate is exonic (the
    variant is not purely intronic) or the c. cannot be read.
    """
    if not hgvs_c:
        return None
    m = C_LOCATION_RE.search(str(hgvs_c))
    if not m:
        return None
    best: Optional[Tuple[str, int]] = None
    for coord in m.group(1).split("_"):
        cm = C_COORD_RE.match(coord)
        if not cm:
            return None
        if not cm.group(2):
            return None  # an exonic coordinate
        off = (cm.group(2), int(cm.group(3)))
        if best is None or off[1] < best[1]:
            best = off
    return best


def _bp7_scope(consequence: Sequence[str], hgvs_c: Optional[str]) -> Tuple[bool, str]:
    """Whether Walker et al. 2023's BP7 scope holds, and why (or why not)."""
    terms = set(consequence)
    if "synonymous_variant" in terms:
        if terms & SPLICE_REGION_TERMS:
            return False, ("synonymous, but VEP places it in the splice region (the first or last 3 bases of the "
                           "exon); Walker et al. 2023: \"BP7 should not be applied for synonymous substitutions ... "
                           "located at the first base or the last 3 bases of the exon\"")
        return True, "synonymous, outside the splice region"
    if "intron_variant" in terms and not (terms & EXONIC_TERMS):
        off = _intron_offset(hgvs_c)
        if off is None:
            return False, ("intronic, but its distance from the exon could not be read from the c. HGVS, so "
                           "Walker et al. 2023's +7/-21 boundary could not be checked")
        sign, n = off
        if (sign == "+" and n >= BP7_DONOR_MIN) or (sign == "-" and n >= BP7_ACCEPTOR_MIN):
            return True, f"intronic at offset {sign}{n} from the exon, at or beyond +{BP7_DONOR_MIN}/-{BP7_ACCEPTOR_MIN}"
        return False, (f"intronic at offset {sign}{n} from the exon, inside the splice region Walker et al. 2023 designates "
                       f"(BP7 only at or beyond +{BP7_DONOR_MIN}/-{BP7_ACCEPTOR_MIN})")
    return False, "not a synonymous or intronic variant"


def suggest(data: Dict[str, Any], inheritance: Optional[str] = None,
            max_credible_af: Optional[float] = None, faf95: Optional[float] = None,
            pm2_max_af: Optional[float] = None) -> Dict[str, Any]:
    """Codes that follow from numbers alone, plus what could not be assessed and why.

    `data` keys used (all optional): revel, alphamissense, spliceai_max,
    spliceai_source (which model, run and score `spliceai_max` came from),
    grpmax_af (or popmax_af), grpmax_an (or gnomad_an), gnomad_ac, faf95,
    max_credible_af, consequence (VEP term list or str), is_missense, hgvs_c
    (for the intronic offset BP7 needs), site_covered / coverage_text (from the
    variant card's gnomAD coverage), mitochondrial (True for an mtDNA variant),
    clinvar_classification (the aggregate ClinVar germline classification string,
    passed in by the caller; this module never fetches it).

    `max_credible_af` and `faf95` given as arguments win over the same keys in
    `data`, so a caller that computes the Whiffin maximum credible AF (zebra
    stats maxaf) can reach BS1 without rewriting the card. `pm2_max_af` is a
    gene-specific PM2 frequency ceiling from a ClinGen VCEP specification,
    supplied by the caller; zebra has none of its own for dominant disorders.

    Returns {"suggested": [...], "not_assessed": [...], "caveats": [...]}.
    """
    out: List[Dict[str, Any]] = []
    not_assessed: List[str] = []
    caveats: List[str] = []

    revel = data.get("revel")
    alphamissense = data.get("alphamissense")
    spliceai = data.get("spliceai_max")
    splice_src = data.get("spliceai_source") or "SpliceAI"
    consequence = data.get("consequence") or []
    if isinstance(consequence, str):
        consequence = [consequence]
    is_missense = data.get("is_missense", any("missense" in c for c in consequence))
    null_like = any(c in NULL_CONSEQUENCES for c in consequence)
    inframe = any(c in ("inframe_deletion", "inframe_insertion") for c in consequence)
    splice_mechanism = any(c in SPLICE_MECHANISM_CONSEQUENCES for c in consequence)
    clinvar = str(data.get("clinvar_classification") or "")
    clinvar_plp = is_plp(clinvar)

    # ---------------------------------------------------------------- protein axis
    protein_pp3: Optional[Tuple[str, str, str]] = None   # (strength, basis, rule)
    protein_bp4: Optional[Tuple[str, str, str]] = None
    protein_excluded = False      # Walker 2023's "protein functional impact has been excluded"
    protein_assessed = False
    if revel is not None and is_missense:
        protein_assessed = True
        s = _band(float(revel), REVEL_PP3, True)
        if s:
            protein_pp3 = (s, f"REVEL {revel} (>= threshold for {s})",
                           "Pejaver et al. 2022: REVEL PP3 Supporting >=0.644, Moderate >=0.773, Strong >=0.932 "
                           "(no predictor reaches very strong on the pathogenic side)")
        b = _band(float(revel), REVEL_BP4, False)
        if b:
            protein_excluded = True
            protein_bp4 = (b, f"REVEL {revel} (<= threshold for {b})",
                           "Pejaver et al. 2022: REVEL BP4 Supporting <=0.290, Moderate <=0.183, Strong <=0.016, VeryStrong <=0.003")
    if alphamissense is not None and is_missense and revel is None:
        # Only when REVEL, the pre-chosen calibrated tool, is absent: reading AlphaMissense after an
        # indeterminate REVEL is the "scanning tools for the strongest evidence" Pejaver 2022 warns of.
        protein_assessed = True
        if float(alphamissense) < ALPHAMISSENSE_BENIGN:
            protein_excluded = True
            if protein_bp4 is None:
                # AlphaMissense has no ClinGen strength calibration here, so it is used only
                # to answer "is protein impact excluded?" -- never to set a BP4 strength.
                caveats.append(
                    f"AlphaMissense {alphamissense} < {ALPHAMISSENSE_BENIGN} (Cheng et al. 2023, Science 381:eadg7492, "
                    f"\"likely benign\" class) excludes protein impact, but carries no ClinGen BP4 strength "
                    f"calibration: it lets a splicing-based BP4 apply, it does not set one on its own."
                )

    # ---------------------------------------------------------------- splicing axis
    splice_pp3: Optional[Tuple[str, str, str]] = None
    splice_bp4_ok = False
    if spliceai is not None:
        sp = float(spliceai)
        if sp >= SPLICEAI_PP3:
            if null_like:
                caveats.append(
                    f"SpliceAI max delta {sp} ({splice_src}) >= {SPLICEAI_PP3}, but the variant class given is "
                    f"already a null / canonical-splice one: PVS1 covers that mechanism and PP3 is not added on top "
                    f"of it (Walker et al. 2023, echoing the SVI recommendation not to use PP3 with PVS1)."
                )
            else:
                # Walker et al. 2023 recommends applying PP3 from SpliceAI at supporting
                # weight only, despite likelihood ratios that would suggest moderate.
                splice_pp3 = ("Supporting", f"SpliceAI max delta {sp} >= {SPLICEAI_PP3} ({splice_src})",
                              "Walker et al. 2023 (ClinGen SVI Splicing): SpliceAI >=0.2 supports PP3 at Supporting "
                              "weight only; do not combine with PVS1 for the same splice effect")
        elif sp <= SPLICEAI_BP4:
            if data.get("splice_benign_blocked"):
                not_assessed.append(f"BP4/BP7 (splicing axis): withheld — {data['splice_benign_blocked']}")
            else:
                splice_bp4_ok = True
        else:
            caveats.append(
                f"SpliceAI max delta {sp} ({splice_src}) falls in the uninformative band {SPLICEAI_UNINFORMATIVE[0]} "
                f"< delta < {SPLICEAI_UNINFORMATIVE[1]} (Walker et al. 2023): it supports neither PP3 nor BP4."
            )

    # ---------------------------------------------------------------- PP3: one code
    # Pejaver et al. 2022 warns against "scanning multiple tools for the strongest
    # evidence", and PP3 counts once. Choose the axis that matches the variant's
    # mechanism -- protein for a missense, splicing otherwise -- and say so.
    if protein_pp3 and splice_pp3:
        caveats.append(
            "both the protein and the splicing predictor cross their PP3 thresholds; PP3 counts once. "
            "The missense predictor is used here (Pejaver et al. 2022: use one pre-chosen tool, do not scan "
            "tools for the strongest evidence). A splice mechanism for a missense variant is worth checking "
            "separately — it does not add a second PP3."
        )
    chosen_pp3 = protein_pp3 or splice_pp3
    if chosen_pp3:
        strength, basis, rule = chosen_pp3
        out.append({"code": _code_label("PP3", strength), "basis": basis, "rule": rule})

    # ---------------------------------------------------------------- BP4 (and BP7)
    splice_bp4_given = False
    if protein_bp4:
        strength, basis, rule = protein_bp4
        if spliceai is not None and float(spliceai) >= SPLICEAI_PP3:
            not_assessed.append(
                f"BP4 (protein axis, REVEL {revel}): withheld because SpliceAI {spliceai} ({splice_src}) >= "
                f"{SPLICEAI_PP3} predicts a splice effect. Richards et al. 2015 defines BP4 as no predicted impact on "
                f"the protein *or* splicing, so the two cannot be read as benign together."
            )
        elif spliceai is not None and float(spliceai) > SPLICEAI_BP4:
            not_assessed.append(
                f"BP4 (protein axis, REVEL {revel}): withheld because SpliceAI {spliceai} ({splice_src}) is in the "
                f"uninformative band {SPLICEAI_BP4} < delta < {SPLICEAI_PP3} (Walker et al. 2023): a splice effect is "
                f"neither predicted nor excluded, and BP4 needs no predicted impact on the protein or splicing."
            )
        elif data.get("splice_benign_blocked"):
            not_assessed.append(f"BP4 (protein axis, REVEL {revel}): withheld — {data['splice_benign_blocked']}")
        else:
            note = basis
            if spliceai is None:
                note += "; splice impact not assessed (no SpliceAI value) — BP4 also requires no predicted splice effect (Richards et al. 2015)"
            else:
                note += f"; no splice effect predicted (SpliceAI {spliceai}, {splice_src})"
            out.append({"code": _code_label("BP4", strength), "basis": note, "rule": rule})
    elif splice_bp4_ok and not null_like and not inframe:
        if is_missense and not protein_excluded:
            # Walker et al. 2023: the splicing benign code applies "for both intronic and
            # synonymous variants, and for missense variants if protein functional impact
            # has been excluded." A low SpliceAI delta on a missense says nothing about the
            # amino-acid change, which is the mechanism in question.
            not_assessed.append(
                f"BP4 (splicing axis, SpliceAI {spliceai} <= {SPLICEAI_BP4}, {splice_src}): not applicable to a "
                f"missense variant until protein impact is excluded (Walker et al. 2023: \"and for missense variants "
                f"if protein functional impact has been excluded\"). "
                + ("No REVEL or AlphaMissense value was supplied." if not protein_assessed
                   else "The missense predictor supplied does not fall in the benign range.")
                + " The low SpliceAI delta means no splicing effect is predicted; it is not an ACMG code."
            )
        elif is_missense or splice_mechanism or not consequence:
            basis = f"SpliceAI max delta {spliceai} <= {SPLICEAI_BP4} ({splice_src})"
            if is_missense:
                basis += "; protein impact excluded by the missense predictor supplied"
            rule = ("Walker et al. 2023: SpliceAI <=0.1 supports BP4 where splicing is the mechanism under "
                    "consideration (intronic, synonymous; missense only once protein impact is excluded). "
                    "Thresholds do not apply at canonical +/-1,+/-2 sites, where PVS1 covers the mechanism.")
            if not consequence:
                basis += "; variant class not supplied"
                caveats.append(
                    "no consequence term was supplied, so the BP4 scope condition (Walker et al. 2023) could not be "
                    "checked: confirm the variant is intronic/synonymous, or that protein impact is excluded."
                )
            out.append({"code": _code_label("BP4", "Supporting"), "basis": basis, "rule": rule})
            splice_bp4_given = True
        else:
            not_assessed.append(
                f"BP4 (splicing axis, SpliceAI {spliceai} <= {SPLICEAI_BP4}): the variant classes given "
                f"({', '.join(consequence)}) are not ones where splicing is the mechanism under consideration."
            )
    if splice_bp4_given:
        # Walker et al. 2023: "BP7 is applicable after assignment of BP4 for no adverse splicing
        # predictions", for synonymous variants outside the first base / last 3 bases of the exon
        # and for intronic variants at or beyond +7/-21. The same paper cautions "against inclusion
        # of evolutionary conservation in assessment of silent or intronic variants for application
        # of code BP7 without empirically derived justification", so conservation is not a gate.
        ok, why = _bp7_scope(consequence, data.get("hgvs_c"))
        if ok:
            out.append({"code": "BP7", "basis": f"{why}; BP4 assigned above from SpliceAI {spliceai} ({splice_src})",
                        "rule": BP7_RULE})
        else:
            not_assessed.append(f"BP7: not offered — {why}.")

    # ---------------------------------------------------------------- frequency axis
    if data.get("mitochondrial"):
        # E-5 / CP1-4: an mtDNA variant has homoplasmic and heteroplasmic frequencies, not a
        # nuclear allele count, and the nuclear BA1/BS1/PM2 cut-offs were never calibrated on them.
        not_assessed.append(MITO_FREQUENCY_NOTE)
        return {"suggested": out, "not_assessed": not_assessed, "caveats": caveats}

    af = _number(data.get("grpmax_af", data.get("popmax_af")))
    an = _number(data.get("grpmax_an", data.get("gnomad_an")))
    ac = _number(data.get("gnomad_ac"))
    if faf95 is None:
        faf95 = data.get("faf95")
    if max_credible_af is None:
        max_credible_af = data.get("max_credible_af")
    if pm2_max_af is None:
        pm2_max_af = data.get("pm2_max_af")

    ba1_basis: Optional[str] = None
    if faf95 is not None:
        # Whiffin et al. 2017 (Genet Med 19:1151) exists because the point estimate is the
        # wrong statistic: the filtering AF is the 95% lower bound, and the stand-alone
        # benign code is the one that most needs the statistical caution. So when a
        # filtering AF is available it DECIDES -- the point estimate is not a fallback that
        # can overrule it (a faf95 sixty-fold below the grpmax estimate used to be ignored).
        if float(faf95) > BA1_AF:
            ba1_basis = f"filtering AF (faf95) {float(faf95):.4g} > {BA1_AF}"
        elif af is not None and float(af) > BA1_AF:
            caveats.append(
                f"grpmax point estimate {float(af):.4g} is above {BA1_AF} but the filtering AF (faf95) is "
                f"{float(faf95):.4g}: BA1 is not offered. The filtering AF is the statistic ClinGen and Whiffin "
                f"et al. 2017 intend for benign frequency criteria; a high point estimate with a low faf95 means "
                f"the frequency rests on few alleles in a small group."
            )
    elif af is not None and float(af) > BA1_AF:
        if an is not None and an >= BA1_MIN_AN:
            ba1_basis = (f"grpmax AF {float(af):.4g} > {BA1_AF} at AN {int(an)} (point estimate; no filtering AF "
                         f"supplied — faf95 is the statistic ClinGen/Whiffin intend here)")
        else:
            not_assessed.append(
                f"BA1: grpmax AF {float(af):.4g} is above {BA1_AF}, but the allele number is "
                f"{'not reported' if an is None else int(an)} and Ghosh et al. 2018 requires >= {BA1_MIN_AN} alleles in "
                f"the group. Supply faf95 or a grpmax AN before applying BA1."
            )
    if ba1_basis:
        suggestion = {"code": "BA1", "basis": ba1_basis,
                      "rule": f"ClinGen SVI (Ghosh et al. 2018, Hum Mutat 39:1525): BA1 when AF > {BA1_AF} in a "
                              f"continental population with >= {BA1_MIN_AN} alleles",
                      "caveats": [BA1_EXCEPTION_NOTE]}
        if clinvar_plp:
            suggestion["conflict"] = (
                f"ClinVar reports this variant as {clinvar!r}. BA1 is stand-alone and would override every "
                f"pathogenic code. A high allele frequency with a Pathogenic/Likely pathogenic assertion is the "
                f"signature of a LOW-PENETRANCE or risk allele (the case Ghosh et al. 2018's exception list was "
                f"written for), not of a benign variant. DO NOT apply BA1 here without resolving the conflict."
            )
            suggestion["requires_review"] = True
            caveats.append(suggestion["conflict"])
        elif clinvar:
            suggestion["caveats"].append(f"ClinVar aggregate classification supplied: {clinvar!r}.")
        else:
            suggestion["caveats"].append(
                "no ClinVar classification was supplied; check ClinVar before applying BA1 — a P/LP assertion at "
                "this frequency means a low-penetrance allele, not a benign one."
            )
        out.append(suggestion)
        return {"suggested": out, "not_assessed": not_assessed, "caveats": caveats}

    # BS1 -- the maximum credible AF is a CEILING: above it the variant is too common to be a
    # fully penetrant cause. Below it, nothing follows (D-P0-2): it is not evidence of rarity.
    bs1_given = False
    if faf95 is not None and max_credible_af is not None:
        if float(faf95) > float(max_credible_af):
            bs1 = {"code": "BS1",
                   "basis": f"filtering AF (faf95) {float(faf95):.4g} > maximum credible AF {float(max_credible_af):.4g}",
                   "rule": "Whiffin et al. 2017 (Genet Med 19:1151): a variant more common than the disease "
                           "allows cannot be a fully penetrant cause",
                   "caveats": []}
            if clinvar_plp:
                # D-P0-1: the BA1 guard, applied to BS1 too. CFTR F508del (ClinVar Pathogenic,
                # 4 stars) reached BS1 when the maximum credible AF came from the wrong formula.
                bs1["conflict"] = (
                    f"ClinVar reports this variant as {clinvar!r}. A Pathogenic/Likely pathogenic assertion on a "
                    f"variant above the maximum credible AF means one of the inputs is wrong (the formula for the "
                    f"inheritance, the prevalence, the allelic or genetic heterogeneity, the penetrance) or the "
                    f"allele has reduced penetrance. DO NOT apply BS1 here without resolving the conflict."
                )
                bs1["requires_review"] = True
                caveats.append(bs1["conflict"])
            elif clinvar:
                bs1["caveats"].append(f"ClinVar aggregate classification supplied: {clinvar!r}.")
            else:
                bs1["caveats"].append("no ClinVar classification was supplied; check ClinVar before applying BS1.")
            out.append(bs1)
            bs1_given = True
        else:
            caveats.append(
                f"BS1 does not apply: filtering AF (faf95) {float(faf95):.4g} <= maximum credible AF "
                f"{float(max_credible_af):.4g}. That only says the variant is not too common for the disease; it is "
                f"not evidence that the variant is rare, and it is not PM2 (ClinGen SVI 2020: PM2 is \"absent from "
                f"controls, or at extremely low frequency if recessive\")."
            )
    else:
        missing = "no maximum credible AF supplied" if max_credible_af is None else "no filtering AF (faf95) supplied"
        not_assessed.append(
            f"BS1: {missing}, so the only calibrated benign frequency criterion could not be assessed. "
            f"Compute it with `zebra stats maxaf --prevalence ... --allelic ... --genetic ... --penetrance ...` "
            f"(or `zebra acmg suggest --prevalence ... --allelic ...`). The absence of BS1 here is not evidence "
            f"that the variant is rare enough."
        )

    # PM2 (ClinGen SVI recommendation v1.0, 2020): "Absent from controls, or at extremely low
    # frequency if recessive", at Supporting weight.
    if bs1_given:
        not_assessed.append("PM2: not assessed — BS1 is offered, and a variant too common for the disease is not rare.")
    elif ac == 0:
        if an == 0:
            not_assessed.append(
                "PM2: allele count is 0 but the allele number is 0 too — the site was not callable in gnomAD, "
                "which is not the same as the variant being absent. Check gnomAD coverage (fraction of "
                "samples at >= 20x) before applying PM2."
            )
        elif an is not None and an < PM2_MIN_AN:
            not_assessed.append(
                f"PM2: allele count is 0, but only {int(an)} alleles were surveyed (< {PM2_MIN_AN}); absence in so "
                f"small a sample carries little information. zebra's floor, not a published PM2 threshold — "
                f"it is Ghosh et al. 2018's BA1 allele-number requirement reused."
            )
        elif data.get("site_covered") is False:
            not_assessed.append(
                "PM2: absent from gnomAD, but the site is not adequately covered "
                f"({data.get('coverage_text') or 'coverage below the threshold'}); absence there is not evidence."
            )
        else:
            basis = f"absent from gnomAD (AC=0{'' if an is None else f', AN={int(an)}'})"
            if data.get("site_covered") is True:
                basis += f"; site covered ({data.get('coverage_text') or 'gnomAD coverage passes the threshold'})"
            elif an is None:
                basis += "; allele number not reported — confirm the site is covered (>= 20x in >= 80% of samples)"
            out.append({"code": "PM2_Supporting", "basis": basis,
                        "rule": "ClinGen SVI PM2 recommendation v1.0 (2020): PM2 at Supporting; absence is evidence "
                                "only at a site gnomAD could call"})
    elif pm2_max_af is not None and af is not None:
        # a gene-specific ceiling from the caller decides alone: the generic conventions below never
        # override it in either direction
        if float(af) <= float(pm2_max_af):
            out.append({"code": "PM2_Supporting",
                        "basis": f"grpmax AF {float(af):.4g} <= {float(pm2_max_af):.4g}, the gene-specific PM2 ceiling "
                                 f"supplied by the caller",
                        "rule": "ClinGen SVI PM2 recommendation v1.0 (2020) at Supporting; the threshold is the "
                                "caller's (a ClinGen VCEP specification for this gene), not zebra's"})
        else:
            not_assessed.append(f"PM2: not offered — grpmax AF {float(af):.4g} is above the gene-specific ceiling "
                                f"{float(pm2_max_af):.4g} supplied by the caller.")
    elif af is not None and inheritance == "AR" and float(af) < RECESSIVE_PM2_AF and an is not None and an < PM2_MIN_AN:
        not_assessed.append(f"PM2: grpmax AF {float(af):.4g} rests on only {int(an)} alleles (< {PM2_MIN_AN}); too few "
                            "to call the variant extremely rare.")
    elif af is not None and inheritance == "AR" and float(af) < RECESSIVE_PM2_AF:
        out.append({"code": "PM2_Supporting", "basis": f"grpmax AF {float(af):.4g} < {RECESSIVE_PM2_AF} (recessive)",
                    "rule": f"ClinGen SVI PM2 recommendation v1.0 (2020): \"at extremely low frequency if recessive\", "
                            f"at Supporting; {RECESSIVE_PM2_AF:.2%} is a convention several ClinGen VCEPs use for "
                            f"recessive genes, not a calibrated threshold — a gene-specific VCEP value wins "
                            f"(--pm2-max-af)"})
    elif ac is None and af is None and data.get("site_covered") is False:
        not_assessed.append(
            "PM2: absent from gnomAD, but the site is not adequately covered "
            f"({data.get('coverage_text') or 'coverage below the threshold'}); absence there is not evidence."
        )
    elif ac or af:
        freq = f"grpmax AF {float(af):.4g}" if af is not None else f"allele count {int(ac)} (no group AF supplied)"
        if inheritance == "AR" and af is not None:
            why = f"{freq} is not below the {RECESSIVE_PM2_AF:.2%} recessive convention"
        elif inheritance == "AR":
            why = f"the variant is present in gnomAD with {freq}, so the recessive convention cannot be checked"
        elif inheritance == "XLR":
            why = (f"the variant is present in gnomAD ({freq}); a hemizygous male is affected by one allele, so the "
                   f"autosomal-recessive {RECESSIVE_PM2_AF:.2%} convention is not applied to X-linked disorders — "
                   f"pass --pm2-max-af from the gene's ClinGen VCEP specification if there is one")
        elif inheritance in ("AD", "XLD"):
            why = (f"the variant is present in gnomAD ({freq}); for a dominant disorder ClinGen SVI 2020 gives PM2 "
                   f"for absence, and zebra has no gene-specific ceiling (pass --pm2-max-af from the gene's ClinGen "
                   f"VCEP specification if there is one)")
        else:
            why = (f"the variant is present in gnomAD ({freq}) and no inheritance was given, so neither the "
                   f"absence rule nor the recessive low-frequency rule could be applied")
        not_assessed.append(f"PM2: not offered — {why}. Being below the disease's maximum credible AF is not "
                            f"evidence of rarity.")
    else:
        not_assessed.append("PM2: not assessed — no gnomAD allele count or group frequency was available for this "
                            "variant (an unanswered lookup is not absence).")
    return {"suggested": out, "not_assessed": not_assessed, "caveats": caveats}


def _number(x: Any) -> Optional[float]:
    """A numeric input as a float; None for None. A non-numeric value is refused, not compared as a string."""
    if x is None:
        return None
    if isinstance(x, bool):
        raise ValueError(f"expected a number, got {x!r}")
    try:
        v = float(x)
    except (TypeError, ValueError):
        raise ValueError(f"expected a number, got {x!r}") from None
    if v != v:
        raise ValueError("expected a number, got NaN")
    return v


def suggest_from_data(data: Dict[str, Any], inheritance: Optional[str] = None,
                      max_credible_af: Optional[float] = None,
                      faf95: Optional[float] = None) -> List[Dict[str, Any]]:
    """`suggest`'s code list alone, for callers that only want the suggestions."""
    return suggest(data, inheritance=inheritance, max_credible_af=max_credible_af, faf95=faf95)["suggested"]


# ------------------------------------------------------------ PVS1 / PS1 / PM5 inputs (CP1-3)
# Structure only: these functions compute what the PVS1 decision tree and the PS1/PM5 criteria
# ask about. They never emit a code; the skill applies the decision tree.

NMD_RULE = (
    "Abou Tayoun et al. 2018 (Hum Mutat 39:1517), ClinGen SVI PVS1 decision tree: a premature termination "
    "codon is predicted to escape nonsense-mediated decay when it lies in the last exon or within the 3'-most "
    "50 nucleotides of the penultimate exon; elsewhere NMD is predicted. Computed here on the transcript's exon "
    "structure (Ensembl), with the first base of the new stop codon as the PTC position."
)
PVS1_MECHANISM_NOTE = (
    "PVS1 needs loss of function to be an established disease mechanism for this gene (Abou Tayoun et al. "
    "2018, first question of the decision tree). ClinGen haploinsufficiency score 3 means 'sufficient evidence'; "
    "LOEUF/pLI describe intolerance in the population and are supporting context, not proof of mechanism."
)
STOP_P_RE = re.compile(r"p\.\(?(?:[A-Z][a-z]{2}|[A-Z])(\d+)(?:Ter|\*|X)\)?$")
FS_P_RE = re.compile(r"p\.\(?(?:[A-Z][a-z]{2}|[A-Z])(\d+)(?:[A-Z][a-z]{2}|[A-Z])?fs(?:(?:Ter|\*|X)(\d+|\?))?\)?$")


def transcript_exons(exons: Sequence[Dict[str, Any]], strand: int) -> List[Dict[str, int]]:
    """Exons in transcript order (5'->3'), each with its genomic span, length and cDNA span."""
    order = sorted(({"start": int(e["start"]), "end": int(e["end"])} for e in exons),
                   key=lambda e: e["start"], reverse=(strand == -1))
    cum = 0
    out = []
    for i, e in enumerate(order, 1):
        length = e["end"] - e["start"] + 1
        out.append({"number": i, "start": e["start"], "end": e["end"], "length": length,
                    "cdna_start": cum + 1, "cdna_end": cum + length})
        cum += length
    return out


def cdna_position(exons: Sequence[Dict[str, int]], genomic: int, strand: int) -> Optional[int]:
    """The cDNA (transcript) coordinate of an exonic genomic position, or None if it is intronic."""
    for e in exons:
        if e["start"] <= genomic <= e["end"]:
            off = genomic - e["start"] if strand == 1 else e["end"] - genomic
            return e["cdna_start"] + off
    return None


def ptc_codon_from_hgvsp(hgvsp: Optional[str], consequence: Sequence[str]) -> Tuple[Optional[int], Optional[int], str]:
    """(first altered residue, codon of the premature stop, how it was read) from the p. HGVS."""
    if not hgvsp:
        return None, None, "no protein HGVS"
    p = hgvsp.split(":", 1)[-1]
    m = FS_P_RE.search(p)
    if m:
        first = int(m.group(1))
        ter = m.group(2)
        if ter and ter.isdigit() and int(ter) >= 1 and first >= 1:
            return first, first + int(ter) - 1, f"frameshift: new stop at residue {first} + {ter} - 1"
        return first, None, "frameshift with no stop codon in the new frame given (fsTer? or bare fs)"
    m = STOP_P_RE.search(p)
    if m and "stop_gained" in consequence and int(m.group(1)) >= 1:
        n = int(m.group(1))
        return n, n, "nonsense: the stop replaces residue " + str(n)
    return None, None, f"no premature stop read from {p}"


def nmd_inputs(exons: Sequence[Dict[str, int]], strand: int, cds_start_genomic: int, ptc_codon: int,
               protein_length: Optional[int]) -> Dict[str, Any]:
    """Where a premature stop at `ptc_codon` lies on the exon structure, and the NMD prediction."""
    if ptc_codon < 1:
        raise ValueError(f"premature stop codon {ptc_codon} is not a residue")
    cds_start = cdna_position(exons, cds_start_genomic, strand)
    if cds_start is None:
        raise ValueError("translation start is not inside an exon of this transcript")
    ptc = cds_start + (ptc_codon - 1) * 3
    n = len(exons)
    exon = next((e for e in exons if e["cdna_start"] <= ptc <= e["cdna_end"]), None)
    if exon is None:
        raise ValueError(f"premature stop at cDNA {ptc} lies beyond the transcript ({exons[-1]['cdna_end']} nt)")
    last_junction = exons[-1]["cdna_start"] - 1 if n > 1 else None
    if n == 1:
        nmd, where = False, "single-exon transcript: no exon junction downstream, NMD not expected"
    elif exon["number"] == n:
        nmd, where = False, "in the last exon"
    elif ptc > last_junction - 50:
        nmd, where = False, (f"{last_junction - ptc + 1} nt upstream of the last exon-exon junction, within the 50 nt "
                             f"window that escapes NMD (the 3'-most 50 nt of the penultimate exon)")
    else:
        nmd, where = True, f"{last_junction - ptc + 1} nt upstream of the last exon-exon junction"
    out: Dict[str, Any] = {
        "ptc_codon": ptc_codon, "ptc_cdna": ptc, "ptc_exon": f"{exon['number']}/{n}",
        "nmd_predicted": nmd, "position": where, "rule": NMD_RULE,
    }
    if protein_length:
        out["protein_length"] = protein_length
        out["fraction_of_protein_after_ptc"] = round(max(0, protein_length - ptc_codon + 1) / protein_length, 3)
    return out


def coding_overlap(exon: Dict[str, int], cds_cdna: Tuple[int, int]) -> int:
    lo, hi = max(exon["cdna_start"], cds_cdna[0]), min(exon["cdna_end"], cds_cdna[1])
    return max(0, hi - lo + 1)


P_TITLE_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|=|\*)\)")
P_CHANGE_RE = re.compile(r"^p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)?$")


def codon_records(records: Sequence[Dict[str, Any]], hgvsp: str, exclude_vcv: Optional[str] = None,
                  gene: Optional[str] = None, exclude_c: Optional[str] = None,
                  exclude_spdi: Optional[str] = None) -> Dict[str, Any]:
    """ClinVar records at the same codon as a missense change, split into PS1 and PM5 inputs.

    PS1: the same amino-acid change as an established pathogenic variant, from a different
    nucleotide change. PM5: a different missense change at a residue where another missense
    change is pathogenic (Richards et al. 2015). The records are matched on the protein change
    in their ClinVar title; the variant itself (`exclude_vcv`) is left out.
    """
    m = P_CHANGE_RE.match((hgvsp or "").split(":", 1)[-1])
    if not m or m.group(3) == "Ter" or m.group(3) == m.group(1):
        raise ValueError(f"not a missense protein change: {hgvsp!r}")
    ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
    own_c = (exclude_c or "").split(":", 1)[-1] or None
    ps1: List[Dict[str, Any]] = []
    pm5: List[Dict[str, Any]] = []
    other: List[Dict[str, Any]] = []
    for r in records:
        if exclude_vcv and r.get("vcv") == exclude_vcv:
            continue
        # PS1 needs a DIFFERENT nucleotide change: the variant's own record is never an input, even
        # when the card could not match it in ClinVar
        if exclude_spdi and str(r.get("canonical_spdi") or "").endswith(exclude_spdi):
            continue
        if own_c and re.search(r":" + re.escape(own_c) + r"(?:\s|$|\()", str(r.get("title") or "")):
            continue
        if r.get("compound"):
            continue
        title = str(r.get("title") or "")
        tm = P_TITLE_RE.search(title)
        if not tm or int(tm.group(2)) != pos or tm.group(1) != ref:
            continue
        if gene and f"({gene})" not in title:
            continue
        r_alt = tm.group(3)
        cls = str(r.get("classification") or "")
        plp = is_plp(cls)
        row = {"vcv": r.get("vcv"), "title": title, "classification": cls or None, "stars": r.get("stars"),
               "stars_text": r.get("stars_text"), "review_status": r.get("review_status"), "url": r.get("url")}
        if r_alt in ("Ter", "*", "="):
            continue  # not a missense change at this residue
        if plp and r_alt == alt:
            ps1.append(row)
        elif plp:
            pm5.append(row)
        else:
            other.append(row)
    key = lambda x: -(x.get("stars") or 0)  # noqa: E731
    return {
        "residue": f"p.{ref}{pos}", "change": f"p.{ref}{pos}{alt}",
        "ps1_inputs": sorted(ps1, key=key), "pm5_inputs": sorted(pm5, key=key),
        "other_missense_at_codon": [{"vcv": o["vcv"], "change": P_TITLE_RE.search(o["title"]).group(0)[1:-1],
                                     "classification": o["classification"], "stars": o["stars"]}
                                    for o in sorted(other, key=key)][:10],
        "rule": ("Inputs, not codes. PS1: same amino-acid change as an established Pathogenic variant from a "
                 "different nucleotide change; PM5: a novel missense change at a residue where a different missense "
                 "change is Pathogenic (Richards et al. 2015). Check that the reference variant's own classification "
                 "is solid (its review stars), that splicing does not explain it, and the gene's VCEP rules; "
                 "the skill decides."),
    }
