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
BA1, PM2_Supporting) are offered by `suggest_from_data` as suggestions that
carry their thresholds and citations.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

POINTS = {"VeryStrong": 8, "Strong": 4, "Moderate": 2, "Supporting": 1}
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

# Pejaver et al. 2022 (AJHG 109:2163), REVEL calibrated thresholds.
REVEL_PP3 = [(0.932, "Strong"), (0.773, "Moderate"), (0.644, "Supporting")]
REVEL_BP4 = [(0.003, "VeryStrong"), (0.016, "Strong"), (0.183, "Moderate"), (0.290, "Supporting")]
# Walker et al. 2023 (AJHG 110:1046), ClinGen SVI Splicing Subgroup.
SPLICEAI_PP3 = 0.2
SPLICEAI_BP4 = 0.1
NULL_CONSEQUENCES = ("stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant",
                     "start_lost", "stop_lost", "transcript_ablation")


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


def richards_class(codes: List[Code]) -> str:
    """The 2015 combining rules, reading modified strengths as their new level."""
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
    if "PM1" in names and "PP3" in names:
        pm1 = next(c for c in codes if c.code == "PM1")
        pp3 = next(c for c in codes if c.code == "PP3")
        if pm1.points + pp3.points > 4:
            warnings.append("PM1 + PP3 exceed Strong together; Pejaver et al. 2022 cap their sum at 4 points (Strong). Capped.")
    if "BA1" in names and any(c.direction == "pathogenic" for c in codes):
        warnings.append("BA1 with pathogenic codes: BA1 is stand-alone; recheck the frequency data and any exception list.")

    total = sum(c.points for c in codes if c.code != "BA1")
    if "PM1" in names and "PP3" in names:
        pm1 = next(c for c in codes if c.code == "PM1")
        pp3 = next(c for c in codes if c.code == "PP3")
        excess = pm1.points + pp3.points - 4
        if excess > 0:
            total -= excess
    if "BA1" in names:
        point_class = "Benign"
        total_display = None
    else:
        point_class = points_class(total)
        total_display = total
    rules = richards_class(codes)
    out = {
        "codes": [{"code": c.code, "label": c.label(), "strength": c.strength, "direction": c.direction,
                   "points": c.points} for c in codes],
        "points": total_display,
        "classification": point_class,
        "classification_richards_2015": rules,
        "agree": point_class == rules,
        "scale": "P >= 10, LP 6..9, VUS 0..5, LB -1..-6, B <= -7 (Tavtigian 2020)",
        "warnings": warnings,
    }
    if not out["agree"]:
        out["note"] = (
            f"Points give {point_class}; the 2015 combining rules give {rules}. ClinGen's Bayesian framework "
            "(Tavtigian 2018/2020) supports the points reading; report both if the difference matters."
        )
    return out


def _band(value: float, bands: List[Tuple[float, str]], higher_is_stronger: bool) -> Optional[str]:
    for threshold, strength in bands:
        if (higher_is_stronger and value >= threshold) or (not higher_is_stronger and value <= threshold):
            return strength
    return None


def suggest_from_data(data: Dict[str, Any], inheritance: Optional[str] = None) -> List[Dict[str, Any]]:
    """Suggest the codes that follow from numbers alone, each with its rule.

    `data` keys used (all optional): revel, spliceai_max, grpmax_af (or popmax_af),
    grpmax_an, gnomad_ac, consequence (VEP term list or str), is_missense.
    """
    out: List[Dict[str, Any]] = []
    revel = data.get("revel")
    spliceai = data.get("spliceai_max")
    consequence = data.get("consequence") or []
    if isinstance(consequence, str):
        consequence = [consequence]
    is_missense = data.get("is_missense", any("missense" in c for c in consequence))

    if revel is not None and is_missense:
        s = _band(float(revel), REVEL_PP3, True)
        if s:
            out.append({"code": "PP3" if s == "Supporting" else f"PP3_{s}", "basis": f"REVEL {revel} (>= threshold for {s})",
                        "rule": "Pejaver et al. 2022: REVEL PP3 Supporting >=0.644, Moderate >=0.773, Strong >=0.932"})
        b = _band(float(revel), REVEL_BP4, False)
        if b:
            out.append({"code": "BP4" if b == "Supporting" else f"BP4_{b}", "basis": f"REVEL {revel} (<= threshold for {b})",
                        "rule": "Pejaver et al. 2022: REVEL BP4 Supporting <=0.290, Moderate <=0.183, Strong <=0.016, VeryStrong <=0.003"})
    null_like = any(c in NULL_CONSEQUENCES for c in consequence)
    inframe = any(c in ("inframe_deletion", "inframe_insertion") for c in consequence)
    if spliceai is not None:
        sp = float(spliceai)
        if sp >= SPLICEAI_PP3 and not null_like:
            out.append({"code": "PP3", "basis": f"SpliceAI max delta {sp} >= 0.2",
                        "rule": "Walker et al. 2023 (ClinGen SVI Splicing): SpliceAI >=0.2 supports PP3; do not combine with PVS1 for the same splice effect"})
        # Walker 2023: BP4/BP7 from SpliceAI only where splicing is the question — not for null variants or in-frame indels
        if sp <= SPLICEAI_BP4 and not null_like and not inframe and (revel is None or float(revel) <= 0.290 or not is_missense):
            out.append({"code": "BP4", "basis": f"SpliceAI max delta {sp} <= 0.1",
                        "rule": "Walker et al. 2023: SpliceAI <=0.1 supports BP4 (and BP7 for synonymous/deep intronic)"})

    af = data.get("grpmax_af", data.get("popmax_af"))
    an = data.get("grpmax_an")
    ac = data.get("gnomad_ac")
    faf95 = data.get("faf95")
    max_af = data.get("max_credible_af")
    if af is not None and float(af) > 0.05 and an is not None and an >= 2000:
        out.append({"code": "BA1", "basis": f"grpmax AF {float(af):.4g} > 0.05",
                    "rule": "ClinGen SVI (Ghosh et al. 2018): BA1 when AF > 0.05 in a continental population with >= 2000 alleles, outside the BA1 exception list"})
    elif faf95 is not None and max_af is not None and float(faf95) > float(max_af):
        out.append({"code": "BS1", "basis": f"filtering AF (faf95) {float(faf95):.4g} > maximum credible AF {float(max_af):.4g}",
                    "rule": "Whiffin et al. 2017: a variant more common than the disease allows cannot be fully penetrant and causal"})
    elif ac == 0:
        out.append({"code": "PM2_Supporting", "basis": "absent from gnomAD (AC=0)",
                    "rule": "ClinGen SVI (2020): PM2 at Supporting; confirm the site is covered in gnomAD"})
    elif af is not None and inheritance in ("AR", "XLR") and float(af) < 0.0007:
        out.append({"code": "PM2_Supporting", "basis": f"grpmax AF {float(af):.4g} < 0.0007 (recessive)",
                    "rule": "ClinGen SVI (2020) PM2_Supporting; 0.07% is a threshold several ClinGen VCEPs use for recessive genes — confirm against the maximum credible AF (zebra stats maxaf)"})
    return out
