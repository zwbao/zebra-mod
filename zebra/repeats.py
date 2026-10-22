"""Repeat sizes read by rule, not by recall: D4Z4 units from a 4q35 fragment size (FSHD1) and
DMPK CTG repeat classes (DM1).

Clinicians and models both convert these by memory and get them wrong (an EcoRI fragment of
17 kb read as 4 units instead of 3; a flank constant borrowed from the EcoRI/BlnI digest). The
conversion is arithmetic on published constants, so it belongs in code: the model reads the
result, never re-derives it.
"""

from __future__ import annotations

import math
import re
from typing import Any, Dict, List, Optional, Tuple

from zebra.core import Outcome, UsageError
from zebra.http import source_record

D4Z4_UNIT_KB = 3.3
# kb of the fragment that is not D4Z4 repeat, per digest. EcoRI: (size - 6.8) / 3.3 and EcoRI/HindIII:
# (size - 4.7) / 3.3 are both written in the Fujian FSHD1 cohort protocol (NCT04369209). That protocol
# gives no EcoRI/BlnI formula: 3.8 is derived, not quoted — the EcoRI flank less the ~3 kb by which
# EcoRI/BlnI fragments are shorter than EcoRI ones (Lemmers et al.).
D4Z4_FLANK_KB = {"ecori": 6.8, "ecori_blni": 3.8, "ecori_hindiii": 4.7}
DIGEST_ALIASES = {"ecori": "ecori", "e": "ecori", "ecori/blni": "ecori_blni", "e/b": "ecori_blni",
                  "ecori_blni": "ecori_blni", "ecori/hindiii": "ecori_hindiii", "e/h": "ecori_hindiii",
                  "ecori_hindiii": "ecori_hindiii"}
# a single-enzyme HindIII or BlnI digest is a different experiment with a different flank: refused by
# name, never quietly read as the double digest
SINGLE_DIGESTS = ("hindiii", "blni", "h", "b")
ACCEPTED_DIGESTS = "EcoRI, EcoRI/BlnI or EcoRI/HindIII"
SIZING_ERROR_KB = 1.0  # PFGE sizing of large fragments is rarely better than about ±1 kb

D4Z4_SOURCES = [
    ("FSHD1 registered cohort protocol (NCT04369209)", "D4Z4 units = (EcoRI fragment kb - 6.8) / 3.3",
     "https://cdn.clinicaltrials.gov/large-docs/09/NCT04369209/Prot_SAP_000.pdf"),
    ("FSHD1 registered cohort protocol (NCT04369209)", "D4Z4 units = (EcoRI/HindIII fragment kb - 4.7) / 3.3",
     "https://cdn.clinicaltrials.gov/large-docs/09/NCT04369209/Prot_SAP_000.pdf"),
    ("Lemmers RJLF et al., digest offsets", "EcoRI/BlnI fragments ~3 kb shorter than EcoRI (flank 3.8 derived)",
     "https://pmc.ncbi.nlm.nih.gov/articles/PMC1182148"),
    ("GeneReviews: Facioscapulohumeral Muscular Dystrophy (NBK1443)",
     "FSHD1: 1-10 units on a permissive 4qA allele; FSHD2: D4Z4 hypomethylation with SMCHD1, DNMT3B or LRIF1",
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


def _half_up(x: float) -> int:
    """Round .5 up. Python's round() rounds half to even, so 3.5 and 4.5 would both give 4."""
    return int(math.floor(x + 0.5))


def _haplotype(text: Optional[str]) -> Optional[str]:
    """"4qA", "4q A", "4qa", "A", "4qA161" → "4qA"; the same for B. None when it names neither."""
    if not text:
        return None
    s = re.sub(r"\s+", "", str(text)).upper()
    if s.startswith("4Q"):
        s = s[2:]
    if s[:1] in ("A", "B"):
        return "4q" + s[0]
    return None


def _locus(text: Optional[str]) -> Optional[str]:
    if not text:
        return None
    s = re.sub(r"\s+", "", str(text)).lower()
    for name in ("4q", "10q"):
        if s.startswith(name):
            return name
    raise UsageError(f"locus must be 4q or 10q as the laboratory states it, got {text!r}")


def _digest(digest: str) -> str:
    raw = re.sub(r"\s+", "", str(digest)).lower().replace("+", "/")
    if raw in SINGLE_DIGESTS:
        raise UsageError(f"a {digest.strip()}-only digest is a different experiment with a different flank, not "
                         "read here; accepted are EcoRI alone or the double digests EcoRI/BlnI and EcoRI/HindIII")
    key = DIGEST_ALIASES.get(raw)
    if key is None:
        raise UsageError(f"digest must be one of {ACCEPTED_DIGESTS}, got {digest!r}")
    return key


def d4z4(fragment_kb: float, digest: str = "ecori", haplotype: Optional[str] = None,
         methylation_pct: Optional[float] = None, locus: Optional[str] = None) -> Outcome:
    key = _digest(digest)
    if not 5 <= fragment_kb <= 400:
        raise UsageError(f"fragment_kb {fragment_kb} is outside 5-400 kb: check the value and unit")
    if methylation_pct is not None and not 0 <= methylation_pct <= 100:
        raise UsageError("methylation_pct must be 0-100")
    where = _locus(locus)
    hap = _haplotype(haplotype)
    flank = D4Z4_FLANK_KB[key]
    exact = (fragment_kb - flank) / D4Z4_UNIT_KB
    units = max(0, _half_up(exact))
    low = max(0, _half_up((fragment_kb - SIZING_ERROR_KB - flank) / D4Z4_UNIT_KB))
    high = _half_up((fragment_kb + SIZING_ERROR_KB - flank) / D4Z4_UNIT_KB)
    band, meaning = _d4z4_band(units)
    warnings: List[str] = []
    if where == "10q":
        warnings.append("the laboratory places this array on 10q: a 10q D4Z4 contraction is not FSHD1")
    elif where is None and key != "ecori_blni":
        # EcoRI and EcoRI/HindIII cut 4q and 10q arrays alike; BlnI cuts within 10q-type units only
        warnings.append(f"an {'EcoRI-only' if key == 'ecori' else 'EcoRI/HindIII'} fragment does not tell a 4q "
                        "from a 10q array; the report's EcoRI/BlnI digest or its statement that the allele is on "
                        "4q is what makes it 4q")
    if hap == "4qB":
        warnings.append("a contraction on 4qB is not pathogenic for FSHD1: re-read which allele is contracted")
    elif haplotype and hap is None:
        warnings.append(f"haplotype {haplotype!r} is not read as 4qA or 4qB: check the report")
    if low != high:
        warnings.append(f"with ±{SIZING_ERROR_KB:g} kb of sizing error the estimate spans {low}-{high} units; "
                        "a laboratory-reported unit count takes precedence")
        if _d4z4_band(low)[0] != _d4z4_band(high)[0]:
            warnings.append(f"that span crosses a band boundary ({_d4z4_band(low)[0]} / {_d4z4_band(high)[0]}): "
                            "the band is not settled by the fragment size alone")
    if key == "ecori" and 35 <= fragment_kb <= 41:
        warnings.append("many laboratories quote EcoRI 10-38 kb as about 1-10 units; this formula reads 38-40 kb "
                        "one unit lower (38 kb → 9), so at the upper boundary the laboratory's own count decides")
    warnings.append("not modelled here: mosaic contractions (a mosaic parent can carry a smaller fragment "
                    "alongside normal alleles) and the milder, less penetrant disease on average in women")
    result: Dict[str, Any] = {
        "locus": "4q35 D4Z4", "fragment_kb": fragment_kb, "digest": key, "flank_kb": flank,
        "units_exact": round(exact, 2), "units": units, "units_range": [low, high],
        "band": band, "meaning": meaning, "formula": f"({fragment_kb:g} - {flank:g}) / {D4Z4_UNIT_KB:g}",
        "rounding": "half up",
    }
    if where:
        result["locus_stated"] = where
    if haplotype:
        result["haplotype"] = haplotype
        result["haplotype_read"] = hap
    text = (f"D4Z4: {fragment_kb:g} kb ({key}) → {units} units (range {low}-{high}); {band}. {meaning}.")
    if methylation_pct is not None:
        result["methylation_pct"] = methylation_pct
        result["methylation_interpreted"] = False
        result["methylation_note"] = (
            "recorded as given, not interpreted by this tool: methylation cut-offs are specific to the laboratory "
            "and the method, so read the value against that laboratory's reference range. Low D4Z4 methylation "
            "with a 4qA allele that is not contracted is the pattern of FSHD2, which is caused by variants in "
            "SMCHD1, DNMT3B or LRIF1 (GeneReviews NBK1443).")
        text += f" Methylation {methylation_pct:g}% recorded, not interpreted here (cut-offs are lab-specific)."
        if units >= 11 and hap != "4qB":
            warnings.append("methylation was measured on an allele that is not contracted: if the laboratory calls "
                            "it low on a 4qA allele, FSHD2 is the question (SMCHD1, DNMT3B, LRIF1 testing)")
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
    if n < 100:
        return "mild (50-~150)", "cataract, mild myotonia, often later onset"
    if n <= 150:
        return "mild / classic overlap (~100-~150)", ("between the mild (50-~150) and classic (~100-~1000) classes, "
                                                       "which overlap here")
    if n <= 1000:
        return "classic (~100-~1000)", "adult-onset classic DM1: weakness, myotonia, cardiac conduction disease, cataract"
    return "congenital (>1000)", "the range GeneReviews (NBK1165) gives for congenital DM1"


def dmpk(ctg: int) -> Outcome:
    if ctg < 0 or ctg > 10000:
        raise UsageError("ctg repeat count must be 0-10000")
    band, meaning = _dmpk_band(ctg)
    warnings = ["classes overlap; blood repeat size underestimates muscle size (somatic instability) and predicts "
                "onset only loosely", "maternal transmission carries the risk of congenital DM1"]
    sources = [source_record(db, rec, url=url) for db, rec, url in DMPK_SOURCES]
    return Outcome({"locus": "DMPK CTG", "repeats": ctg, "band": band, "meaning": meaning}, sources=sources,
                   warnings=warnings, text=f"DMPK CTG {ctg}: {band}. {meaning}.", query={"locus": "dmpk", "ctg": ctg})
