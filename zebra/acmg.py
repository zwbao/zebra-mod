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

NULL_CONSEQUENCES = ("stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant",
                     "start_lost", "stop_lost", "transcript_ablation")
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


def richards_class(codes: Sequence[Code]) -> str:
    """The 2015 combining rules.

    Modified strengths are read as their new level on the pathogenic side
    (Supporting/Moderate/Strong/VeryStrong map onto pp/pm/ps/pvs). Richards 2015
    has no benign Moderate tier, so a benign code at Moderate is read as
    Supporting and one at Strong or VeryStrong as Strong -- the conservative
    mapping, since the rules were written for two benign tiers only.
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


def _weaken(code: Code, by_points: int) -> Optional[Code]:
    """The same code with `by_points` taken off it, or None if nothing is left."""
    target = abs(code.points) - by_points
    for name, pts in sorted(POINTS.items(), key=lambda kv: kv[1]):
        if pts == target:
            sign = 1 if code.direction == "pathogenic" else -1
            return Code(code.code, name, code.direction, sign * pts, True)
    return None


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
    for a, b in (("PP3", "BP4"), ("PP3", "BP7"), ("BP4", "BP7")):
        if a in names and b in names:
            warnings.append(
                f"{a} with {b}: these are opposite (or duplicate) readings of the same computational evidence and "
                f"cannot both apply — resolve which axis (protein or splicing) the prediction speaks to and keep one."
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
        excess = pm1.points + pp3.points - POINTS["Strong"]
        if excess > 0:
            cap_applied = True
            warnings.append("PM1 + PP3 exceed Strong together; Pejaver et al. 2022 cap their sum at 4 points (Strong). Capped in both readings.")
            weaker = _weaken(pp3, excess)
            capped = [c for c in codes if c.code != "PP3"] + ([weaker] if weaker else [])

    total = sum(c.points for c in capped if c.code != "BA1")
    computational_only = bool(codes) and all(c.code in COMPUTATIONAL_CODES for c in codes)
    if "BA1" in names:
        point_class = "Benign"
        total_display = None
    else:
        point_class = points_class(total)
        total_display = total
    rules = richards_class(capped)

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

    out: Dict[str, Any] = {
        "codes": [{"code": c.code, "label": c.label(), "strength": c.strength, "direction": c.direction,
                   "points": c.points} for c in codes],
        "points": total_display,
        "classification": point_class,
        "classification_richards_2015": rules,
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


def suggest(data: Dict[str, Any], inheritance: Optional[str] = None,
            max_credible_af: Optional[float] = None, faf95: Optional[float] = None) -> Dict[str, Any]:
    """Codes that follow from numbers alone, plus what could not be assessed and why.

    `data` keys used (all optional): revel, alphamissense, spliceai_max,
    grpmax_af (or popmax_af), grpmax_an (or gnomad_an), gnomad_ac, faf95,
    max_credible_af, consequence (VEP term list or str), is_missense,
    clinvar_classification (the aggregate ClinVar germline classification string,
    passed in by the caller; this module never fetches it).

    `max_credible_af` and `faf95` given as arguments win over the same keys in
    `data`, so a caller that computes the Whiffin maximum credible AF (zebra
    stats maxaf) can reach BS1 without rewriting the card.

    Returns {"suggested": [...], "not_assessed": [...], "caveats": [...]}.
    """
    out: List[Dict[str, Any]] = []
    not_assessed: List[str] = []
    caveats: List[str] = []

    revel = data.get("revel")
    alphamissense = data.get("alphamissense")
    spliceai = data.get("spliceai_max")
    consequence = data.get("consequence") or []
    if isinstance(consequence, str):
        consequence = [consequence]
    is_missense = data.get("is_missense", any("missense" in c for c in consequence))
    null_like = any(c in NULL_CONSEQUENCES for c in consequence)
    inframe = any(c in ("inframe_deletion", "inframe_insertion") for c in consequence)
    splice_mechanism = any(c in SPLICE_MECHANISM_CONSEQUENCES for c in consequence)
    clinvar = str(data.get("clinvar_classification") or "")

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
    if alphamissense is not None and is_missense:
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
                    f"SpliceAI max delta {sp} >= {SPLICEAI_PP3}, but the variant class given is already a null / "
                    f"canonical-splice one: PVS1 covers that mechanism and PP3 is not added on top of it "
                    f"(Walker et al. 2023, echoing the SVI recommendation not to use PP3 with PVS1)."
                )
            else:
                # Walker et al. 2023 recommends applying PP3 from SpliceAI at supporting
                # weight only, despite likelihood ratios that would suggest moderate.
                splice_pp3 = ("Supporting", f"SpliceAI max delta {sp} >= {SPLICEAI_PP3}",
                              "Walker et al. 2023 (ClinGen SVI Splicing): SpliceAI >=0.2 supports PP3 at Supporting "
                              "weight only; do not combine with PVS1 for the same splice effect")
        elif sp <= SPLICEAI_BP4:
            splice_bp4_ok = True
            # Walker et al. 2023 permits a splicing predictor to carry BP7's splicing half, but
            # BP7 also needs the nucleotide not to be highly conserved (Richards et al. 2015),
            # which zebra does not check -- so BP7 is never emitted here.
            caveats.append(
                "BP7 is not offered by this module. It needs a synonymous or deep-intronic variant AND the "
                "nucleotide not to be highly conserved (Richards et al. 2015); zebra does not check conservation, "
                "so BP7 stays a judgement call."
            )
        else:
            caveats.append(
                f"SpliceAI max delta {sp} falls in the uninformative band {SPLICEAI_UNINFORMATIVE[0]} < delta < "
                f"{SPLICEAI_UNINFORMATIVE[1]} (Walker et al. 2023): it supports neither PP3 nor BP4."
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

    # ---------------------------------------------------------------- BP4
    if protein_bp4:
        strength, basis, rule = protein_bp4
        if spliceai is not None and float(spliceai) >= SPLICEAI_PP3:
            not_assessed.append(
                f"BP4 (protein axis, REVEL {revel}): withheld because SpliceAI {spliceai} >= {SPLICEAI_PP3} "
                f"predicts a splice effect. Richards et al. 2015 defines BP4 as no predicted impact on the "
                f"protein *or* splicing, so the two cannot be read as benign together."
            )
        else:
            note = basis
            if spliceai is None:
                note += "; splice impact not assessed (no SpliceAI value) — BP4 also requires no predicted splice effect (Richards et al. 2015)"
            out.append({"code": _code_label("BP4", strength), "basis": note, "rule": rule})
    elif splice_bp4_ok and not null_like and not inframe:
        if is_missense and not protein_excluded:
            # Walker et al. 2023: the splicing benign code applies "for both intronic and
            # synonymous variants, and for missense variants if protein functional impact
            # has been excluded." A low SpliceAI delta on a missense says nothing about the
            # amino-acid change, which is the mechanism in question.
            not_assessed.append(
                f"BP4 (splicing axis, SpliceAI {spliceai} <= {SPLICEAI_BP4}): not applicable to a missense variant "
                f"until protein impact is excluded (Walker et al. 2023: \"and for missense variants if protein "
                f"functional impact has been excluded\"). "
                + ("No REVEL or AlphaMissense value was supplied." if not protein_assessed
                   else "The missense predictor supplied does not fall in the benign range.")
                + " The low SpliceAI delta means no splicing effect is predicted; it is not an ACMG code."
            )
        elif is_missense or splice_mechanism or not consequence:
            basis = f"SpliceAI max delta {spliceai} <= {SPLICEAI_BP4}"
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
        else:
            not_assessed.append(
                f"BP4 (splicing axis, SpliceAI {spliceai} <= {SPLICEAI_BP4}): the variant classes given "
                f"({', '.join(consequence)}) are not ones where splicing is the mechanism under consideration."
            )

    # ---------------------------------------------------------------- frequency axis
    af = data.get("grpmax_af", data.get("popmax_af"))
    an = data.get("grpmax_an", data.get("gnomad_an"))
    ac = data.get("gnomad_ac")
    if faf95 is None:
        faf95 = data.get("faf95")
    if max_credible_af is None:
        max_credible_af = data.get("max_credible_af")

    ba1_basis: Optional[str] = None
    if faf95 is not None and float(faf95) > BA1_AF:
        # Whiffin et al. 2017 (Genet Med 19:1151) exists because the point estimate is the
        # wrong statistic: the filtering AF is the 95% lower bound, and the stand-alone
        # benign code is the one that most needs the statistical caution.
        ba1_basis = f"filtering AF (faf95) {float(faf95):.4g} > {BA1_AF}"
    elif af is not None and float(af) > BA1_AF:
        if an is not None and an >= BA1_MIN_AN:
            ba1_basis = (f"grpmax AF {float(af):.4g} > {BA1_AF} at AN {an} (point estimate; no filtering AF "
                         f"supplied — faf95 is the statistic ClinGen/Whiffin intend here)")
        else:
            not_assessed.append(
                f"BA1: grpmax AF {float(af):.4g} is above {BA1_AF}, but the allele number is "
                f"{'not reported' if an is None else an} and Ghosh et al. 2018 requires >= {BA1_MIN_AN} alleles in "
                f"the group. Supply faf95 or a grpmax AN before applying BA1."
            )
    if ba1_basis:
        suggestion = {"code": "BA1", "basis": ba1_basis,
                      "rule": f"ClinGen SVI (Ghosh et al. 2018, Hum Mutat 39:1525): BA1 when AF > {BA1_AF} in a "
                              f"continental population with >= {BA1_MIN_AN} alleles",
                      "caveats": [BA1_EXCEPTION_NOTE]}
        if CLINVAR_PLP_RE.search(clinvar):
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

    pm2_given = False
    if not ba1_basis:
        if faf95 is not None and max_credible_af is not None:
            if float(faf95) > float(max_credible_af):
                out.append({"code": "BS1",
                            "basis": f"filtering AF (faf95) {float(faf95):.4g} > maximum credible AF {float(max_credible_af):.4g}",
                            "rule": "Whiffin et al. 2017 (Genet Med 19:1151): a variant more common than the disease "
                                    "allows cannot be a fully penetrant cause"})
            else:
                # P1-4: the extremely-low-frequency arm of PM2 is not restricted to recessive
                # disorders. The right gate for a dominant disorder is the maximum credible AF
                # the same framework provides.
                out.append({"code": "PM2_Supporting",
                            "basis": f"filtering AF (faf95) {float(faf95):.4g} < maximum credible AF {float(max_credible_af):.4g}",
                            "rule": "ClinGen SVI (2020) PM2 at Supporting; the frequency gate is Whiffin et al. 2017's "
                                    "maximum credible AF, which applies to dominant and recessive disorders alike"})
                pm2_given = True
        else:
            missing = "no maximum credible AF supplied" if max_credible_af is None else "no filtering AF (faf95) supplied"
            not_assessed.append(
                f"BS1: {missing}, so the only calibrated benign frequency criterion could not be assessed. "
                f"Compute it with `zebra stats maxaf --prevalence ... --allelic ... --genetic ... --penetrance ...` "
                f"(or `zebra acmg suggest --prevalence ... --allelic ...`). The absence of BS1 here is not evidence "
                f"that the variant is rare enough."
            )

        if not pm2_given:
            if ac == 0:
                if an == 0:
                    not_assessed.append(
                        "PM2: allele count is 0 but the allele number is 0 too — the site was not callable in gnomAD, "
                        "which is not the same as the variant being absent. Check gnomAD coverage (fraction of "
                        "samples at >= 20x) before applying PM2."
                    )
                elif an is not None and an < PM2_MIN_AN:
                    not_assessed.append(
                        f"PM2: allele count is 0, but only {an} alleles were surveyed (< {PM2_MIN_AN}); absence in so "
                        f"small a sample carries little information. zebra's floor, not a published PM2 threshold — "
                        f"it is Ghosh et al. 2018's BA1 allele-number requirement reused."
                    )
                else:
                    basis = f"absent from gnomAD (AC=0{'' if an is None else f', AN={an}'})"
                    if an is None:
                        basis += "; allele number not reported — confirm the site is covered (>= 20x in >= 80% of samples)"
                    out.append({"code": "PM2_Supporting", "basis": basis,
                                "rule": "ClinGen SVI (2020): PM2 at Supporting; absence is evidence only at a site "
                                        "gnomAD could call"})
                    pm2_given = True
            elif af is not None and inheritance in ("AR", "XLR") and float(af) < RECESSIVE_PM2_AF:
                out.append({"code": "PM2_Supporting", "basis": f"grpmax AF {float(af):.4g} < {RECESSIVE_PM2_AF} (recessive)",
                            "rule": f"ClinGen SVI (2020) PM2_Supporting; {RECESSIVE_PM2_AF:.2%} is a convention several "
                                    f"ClinGen VCEPs use for recessive genes, not a calibrated threshold — confirm "
                                    f"against the maximum credible AF (zebra stats maxaf)"})
                pm2_given = True

    return {"suggested": out, "not_assessed": not_assessed, "caveats": caveats}


def suggest_from_data(data: Dict[str, Any], inheritance: Optional[str] = None,
                      max_credible_af: Optional[float] = None,
                      faf95: Optional[float] = None) -> List[Dict[str, Any]]:
    """`suggest`'s code list alone, for callers that only want the suggestions."""
    return suggest(data, inheritance=inheritance, max_credible_af=max_credible_af, faf95=faf95)["suggested"]
