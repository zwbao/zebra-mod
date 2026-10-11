"""Repeat sizes read by rule, not by recall: D4Z4 units from a 4q35 fragment size (FSHD1) and
DMPK CTG repeat classes (DM1).

Clinicians and models both convert these by memory and get them wrong (an EcoRI fragment of
17 kb read as 4 units instead of 3; a flank constant borrowed from the EcoRI/BlnI digest). The
conversion is arithmetic on published constants, so it belongs in code: the model reads the
result, never re-derives it.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from zebra.core import Outcome, UsageError
from zebra.http import source_record

D4Z4_UNIT_KB = 3.3
# kb of the fragment that is not D4Z4 repeat, per digest. EcoRI: (size - 6.8) / 3.3 (the formula
# of the Fujian FSHD1 cohort protocol, NCT04369209). EcoRI/BlnI fragments are ~3 kb shorter than
# EcoRI ones and EcoRI/HindIII ~2 kb shorter (Lemmers et al.), so their flanks are smaller by that.
D4Z4_FLANK_KB = {"ecori": 6.8, "ecori_blni": 3.8, "ecori_hindiii": 4.8}
DIGEST_ALIASES = {"ecori": "ecori", "e": "ecori", "ecori/blni": "ecori_blni", "e/b": "ecori_blni", "blni": "ecori_blni",
                  "ecori_blni": "ecori_blni", "ecori/hindiii": "ecori_hindiii", "hindiii": "ecori_hindiii",
                  "ecori_hindiii": "ecori_hindiii"}
SIZING_ERROR_KB = 1.0  # PFGE sizing of large fragments is rarely better than about ±1 kb

D4Z4_SOURCES = [
    ("FSHD1 registered cohort protocol (NCT04369209)", "D4Z4 units = (EcoRI fragment kb - 6.8) / 3.3",
     "https://cdn.clinicaltrials.gov/large-docs/09/NCT04369209/Prot_SAP_000.pdf"),
    ("Lemmers RJLF et al., digest offsets", "EcoRI/BlnI ~3 kb and EcoRI/HindIII ~2 kb shorter than EcoRI",
     "https://pmc.ncbi.nlm.nih.gov/articles/PMC1182148"),
    ("GeneReviews: Facioscapulohumeral Muscular Dystrophy (NBK1443)", "FSHD1: 1-10 units on a permissive 4qA allele",
     "https://www.ncbi.nlm.nih.gov/books/NBK1443/"),
]

DMPK_SOURCES = [
    ("GeneReviews: Myotonic Dystrophy Type 1 (NBK1165)", "DMPK CTG repeat classes",
     "https://www.ncbi.nlm.nih.gov/books/NBK1165/"),
]


def _d4z4_band(units: int) -> Tuple[str, str]:
    if units <= 0:
        return "below one unit", "a fragment this small is outside the sizes the formula was built on: check the digest and the report"
    if units <= 3:
        return "1-3 units (short)", ("FSHD1 range, short end: on average earlier onset and more severe disease; "
                                     "extramuscular involvement (retinal vasculopathy, hearing loss) is more frequent — "
                                     "the 'large deletion' group for which screening is advised")
    if units <= 6:
        return "4-6 units (intermediate)", "FSHD1 range: intermediate severity on average, wide variability"
    if units <= 10:
        return "7-10 units (long end of FSHD1)", ("FSHD1 range, long end: milder on average, reduced penetrance, "
                                                  "asymptomatic carriers common; methylation helps interpret")
    return "11 or more units", ("outside the FSHD1 range (normal 11-100); FSHD2 is still possible with D4Z4 "
                                "hypomethylation on a 4qA allele and an SMCHD1/DNMT3B/LRIF1 variant")


def d4z4(fragment_kb: float, digest: str = "ecori", haplotype: Optional[str] = None,
         methylation_pct: Optional[float] = None) -> Outcome:
    key = DIGEST_ALIASES.get(str(digest).strip().lower())
    if key is None:
        raise UsageError(f"digest must be one of EcoRI, EcoRI/BlnI, EcoRI/HindIII, got {digest!r}")
    if not 5 <= fragment_kb <= 400:
        raise UsageError(f"fragment_kb {fragment_kb} is outside 5-400 kb: check the value and unit")
    flank = D4Z4_FLANK_KB[key]
    exact = (fragment_kb - flank) / D4Z4_UNIT_KB
    units = max(0, round(exact))
    low = max(0, round((fragment_kb - SIZING_ERROR_KB - flank) / D4Z4_UNIT_KB))
    high = round((fragment_kb + SIZING_ERROR_KB - flank) / D4Z4_UNIT_KB)
    band, meaning = _d4z4_band(units)
    warnings: List[str] = []
    if key == "ecori":
        warnings.append("an EcoRI-only fragment does not tell a 4q from a 10q array; the report's BlnI digest or "
                        "4qA/4qB typing is what makes it 4q")
    if haplotype and haplotype.strip().upper().startswith("4QB"):
        warnings.append("a contraction on 4qB is not pathogenic for FSHD1: re-read which allele is contracted")
    if low != high:
        warnings.append(f"with ±{SIZING_ERROR_KB:g} kb of sizing error the estimate spans {low}-{high} units; "
                        "a laboratory-reported unit count takes precedence")
    result: Dict[str, Any] = {
        "locus": "4q35 D4Z4", "fragment_kb": fragment_kb, "digest": key, "flank_kb": flank,
        "units_exact": round(exact, 2), "units": units, "units_range": [low, high],
        "band": band, "meaning": meaning, "formula": f"({fragment_kb:g} - {flank:g}) / {D4Z4_UNIT_KB:g}",
    }
    if haplotype:
        result["haplotype"] = haplotype
    if methylation_pct is not None:
        if not 0 <= methylation_pct <= 100:
            raise UsageError("methylation_pct must be 0-100")
        result["methylation_pct"] = methylation_pct
        result["methylation_note"] = ("lower D4Z4 methylation goes with a shorter array and, on average, more severe "
                                      "disease; compare with the laboratory's own reference range")
    text = (f"D4Z4: {fragment_kb:g} kb ({key}) → {units} units (range {low}-{high}); {band}. {meaning}.")
    sources = [source_record(db, rec, url=url) for db, rec, url in D4Z4_SOURCES]
    return Outcome(result, sources=sources, warnings=warnings, text=text,
                   query={"locus": "d4z4", "fragment_kb": fragment_kb, "digest": key})


def _dmpk_band(n: int) -> Tuple[str, str]:
    if n < 5:
        return "below the usual range", "check the value"
    if n <= 34:
        return "normal (5-34)", "not DM1"
    if n <= 49:
        return "premutation (35-49)", "usually asymptomatic; can expand in the next generation"
    if n <= 150:
        return "mild / late-onset (50-~150)", "cataract, mild myotonia, often later onset"
    if n <= 1000:
        return "classic (~100-~1000)", "adult-onset classic DM1: weakness, myotonia, cardiac conduction disease, cataract"
    return "large (>1000)", "the range of childhood and congenital DM1"


def dmpk(ctg: int) -> Outcome:
    if ctg < 0 or ctg > 10000:
        raise UsageError("ctg repeat count must be 0-10000")
    band, meaning = _dmpk_band(ctg)
    warnings = ["classes overlap; blood repeat size underestimates muscle size (somatic instability) and predicts "
                "onset only loosely", "maternal transmission carries the risk of congenital DM1"]
    sources = [source_record(db, rec, url=url) for db, rec, url in DMPK_SOURCES]
    return Outcome({"locus": "DMPK CTG", "repeats": ctg, "band": band, "meaning": meaning}, sources=sources,
                   warnings=warnings, text=f"DMPK CTG {ctg}: {band}. {meaning}.", query={"locus": "dmpk", "ctg": ctg})
