"""The report forms rare-disease families actually hold, besides an SNV/indel list.

Four kinds of result can be entered and analysed here, each with its own result
type. None of them is classified by zebra: what is produced are the *inputs* a
clinician or a curator needs, every one of them traceable to Ensembl, ClinGen or
a retrieved guideline text.

- **cnv** — a CMA / CNV-seq / array result as coordinates or an ISCN string
  (`arr[GRCh38] 15q11.2q13.1(23123715_28193120)x1`). Reports the genes the
  interval spans (Ensembl overlap), every ClinGen-curated *region* it overlaps
  (the recurrent microdeletion/duplication regions and the population regions,
  with how much of each it covers), ClinGen gene-level dosage for every spanned
  gene, and the ACMG/ClinGen CNV evidence *inputs* (Riggs 2020 sections 1-5).
- **exon_cnv** — an MLPA-style exon deletion or duplication
  (`DMD exon 45-50 deletion`, `NM_004006.3:c.6439-?_7309+?del`). Resolves the
  exons on the reference transcript, counts the *coding* bases removed (exons
  clipped to the translation), states the frame consequence and — only for an
  out-of-frame deletion — which flanking exon's removal would restore the frame.
- **copy_number** — an SMN1/SMN2-style copy-number result. For SMN1/SMN2 the
  meaning of the combination (diagnosis, carrier status, SMN2 as the modifier)
  is given from the retrieved GeneReviews chapter, labelled research reference.
- **repeat_expansion** — a repeat-expansion result (`FMR1 CGG 230`). For FMR1,
  HTT, DMPK, FXN and C9orf72 the repeat categories are given from the retrieved
  GeneReviews chapters, with their citation, labelled research reference.

Nothing is inferred from the input: a form that does not carry coordinates (a
bare ISCN band, `del(15)(q11.2q13.1)`) is refused rather than guessed, and any
text the parser does not understand is refused rather than ignored.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

import json

from zebra.core import Outcome, UsageError
from zebra.http import SourceError, get_json
from zebra.sources import attempt, clingen, ensembl
from zebra.sources import record as source_record
from zebra.vcf import CHROM_LENGTHS, PAR_X  # Ensembl assembly lengths, GRC pseudoautosomal regions

ACMG_CNV_FRAMEWORK = ("ACMG/ClinGen technical standard for the interpretation and reporting of constitutional "
                      "copy-number variants (Riggs et al., Genet Med 2020;22:245-257)")
GENE_SYMBOL_LIMIT = 200
BENIGN_LIMIT = 15  # population regions listed (the ones containing the CNV first); the total is always given
SMALL_VARIANT_BP = 50  # below this, with exact breakpoints, it is a small variant, not an exon event
MT_LENGTH = 16569  # Ensembl /info/assembly/homo_sapiens/MT, GRCh38 and GRCh37 alike (retrieved 2026-10-06)
INPUT_LIMIT = 400

RESEARCH_LABEL = ("research reference: what the cited source states, retrieved on the date given. It is not a "
                  "diagnosis or a classification of this result; the laboratory's report and a genetics "
                  "professional decide.")

# --------------------------------------------------------------- sources
# Every threshold and reading below was taken from these texts, retrieved on
# 2026-10-06. NCBI Bookshelf answers scripted clients with a reCAPTCHA page, so the
# GeneReviews chapters were read in a browser and their page text saved; sha256 is
# of that saved text. tools/cnv/check_loci_sources.py re-checks every quoted
# sentence below against saved copies of these pages.
SOURCES: Dict[str, Dict[str, Any]] = {
    "GR_SMA": {"citation": "Prior TW, Leach ME, Finanger EL. Spinal Muscular Atrophy. GeneReviews (NBK1352)",
               "url": "https://www.ncbi.nlm.nih.gov/books/NBK1352/", "pmid": "20301526",
               "revision": "Last Revision: February 12, 2026", "retrieved_at": "2026-10-06T09:47:04Z",
               "sha256_of_page_text": "a20890fe57f39bd3f29c0ac99035eed7db5758a2ca6b9b30f0dbb8601b7ca3ad"},
    "GR_FMR1": {"citation": "Hunter JE, Berry-Kravis E, Hipp H, Todd PK. FMR1 Disorders. GeneReviews (NBK1384)",
                "url": "https://www.ncbi.nlm.nih.gov/books/NBK1384/", "pmid": "20301558",
                "revision": "Last Revision: May 16, 2024", "retrieved_at": "2026-10-06T09:47:11Z",
                "sha256_of_page_text": "f01337a19b4561e4ca6cb9cf21019138d0c532a594e475cfbcb358a64b58a778"},
    "GR_HD": {"citation": "Caldeira Bras I, Dawson J, Kay C, Caron NS, Hayden MR. Huntington Disease. "
                          "GeneReviews (NBK1305)",
              "url": "https://www.ncbi.nlm.nih.gov/books/NBK1305/", "pmid": "20301482",
              "revision": "Last Update: February 12, 2026", "retrieved_at": "2026-10-06T09:47:23Z",
              "sha256_of_page_text": "2e77cfaed6a7692260b08a1f78ed68127f3eef4dd16a633035ee74da3521fe96"},
    "GR_DM1": {"citation": "Bird TD. Myotonic Dystrophy Type 1. GeneReviews (NBK1165)",
               "url": "https://www.ncbi.nlm.nih.gov/books/NBK1165/", "pmid": "20301344",
               "revision": "Last Revision: November 14, 2024", "retrieved_at": "2026-10-06T09:47:29Z",
               "sha256_of_page_text": "3e29068eb42f4daad74ea4430bf4f5588da84ec96081cec087136774acf72d65"},
    "GR_FRDA": {"citation": "Bidichandani SI, Delatycki MB, Napierala M, Duquette A. Friedreich Ataxia. "
                            "GeneReviews (NBK1281)",
                "url": "https://www.ncbi.nlm.nih.gov/books/NBK1281/", "pmid": "20301458",
                "revision": "Last Revision: June 26, 2025", "retrieved_at": "2026-10-06T09:47:37Z",
                "sha256_of_page_text": "0f1c8bd1934ce57d266483faabc4ca2fb9e71d10bda1eeeb564143eefcf096cc"},
    "GR_C9": {"citation": "Gossye H, Engelborghs S, Van Broeckhoven C, van der Zee J. C9orf72 Frontotemporal "
                          "Dementia and/or Amyotrophic Lateral Sclerosis. GeneReviews (NBK268647)",
              "url": "https://www.ncbi.nlm.nih.gov/books/NBK268647/", "pmid": "25577942",
              "revision": "Last Update: December 17, 2020", "retrieved_at": "2026-10-06T09:47:43Z",
              "sha256_of_page_text": "36f8e7e189e54ef854f60c7472b98d6caa8e6baec76756556411dcfd9728f70c"},
    "EMQN_DM": {"citation": "Kamsteeg EJ et al. Best practice guidelines and recommendations on the molecular "
                            "diagnosis of myotonic dystrophy types 1 and 2. Eur J Hum Genet 2012 (PMC3499739)",
                "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3499739/fullTextXML", "pmid": "22643181",
                "retrieved_at": "2026-10-06", "sha256_of_body":
                "55291052facad3e4ef16012f659ee53b3e9182eb6a7ce150e690ff4d7b0beac9",
                "quotes": ["5–35 (normal range)", "36–50", "51–150", ">150"]},
    "SMN_EXONS": {"citation": "Akhkiamova M et al. Rare Variants of the SMN1 Gene Detected during Neonatal "
                              "Screening. Genes (Basel) 2024;15:956 (PMC11275604)",
                  "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11275604/fullTextXML",
                  "pmid": "39062735", "retrieved_at": "2026-10-06", "sha256_of_body":
                  "c0f1d870ef6238496cda54e2d8ae9580a2e1028db93fe2ea2138176c421211ae",
                  "quotes": ["Both genes include nine exons (1, 2a, 2b, 3–8)"]},
    "RIGGS_2020": {"citation": "Riggs ER, Andersen EF, Cherry AM, et al. Technical standards for the interpretation "
                               "and reporting of constitutional copy-number variants. Genet Med 2020;22:245-257",
                   "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7313390/", "pmid": "31690835",
                   "retrieved_at": "2026-10-06", "sha256_of_body":
                   "cc856b8e96075b61116219c1bb91f16b93573f095f06596899e9091266b6f324",
                   "quotes": ["3A. 0-24 genes 0", "3B. 25-34 genes 0.45", "3C. 35+ genes 0.90", "3A. 0-34 genes 0",
                              "3B. 35-49 genes 0.45", "3C. 50 or more genes 0.90",
                              "will indicate carrier status for autosomal recessive or X-linked disorders mapping within "
                              "the CNV interval",
                              "Given the significant reproductive risk to female carriers of X-linked conditions",
                              "females may manifest symptoms in many X-linked disorders"]},
}


def _cite(key: str) -> Dict[str, Any]:
    src = SOURCES[key]
    return {k: src[k] for k in ("citation", "url", "pmid", "revision", "retrieved_at") if src.get(k)}


# ACMG/ClinGen section 3 (gene number). Riggs 2020 Table 1 (copy-number loss) and
# Table 2 (copy-number gain), "Number of protein-coding RefSeq genes wholly or
# partially included": loss "3A. 0-24 genes 0", "3B. 25-34 genes 0.45",
# "3C. 35+ genes 0.90"; gain "3A. 0-34 genes 0", "3B. 35-49 genes 0.45",
# "3C. 50 or more genes 0.90".
SECTION3_BANDS: Dict[str, List[Tuple[str, int, Optional[int], float]]] = {
    "loss": [("3A", 0, 24, 0.0), ("3B", 25, 34, 0.45), ("3C", 35, None, 0.90)],
    "gain": [("3A", 0, 34, 0.0), ("3B", 35, 49, 0.45), ("3C", 50, None, 0.90)],
}
_BAND_TEXT = {"loss": "0-24, 25-34 and 35+ (Riggs 2020 Table 1, copy-number loss)",
              "gain": "0-34, 35-49 and 50+ (Riggs 2020 Table 2, copy-number gain)"}

# ------------------------------------------------------------- input forms

_NUM = r"\d[\d,]*"  # digits, optionally with thousands commas (checked in _int)
_BUILD_WORDS = r"GRCh3[78]|hg19|hg38|hg18|NCBI3[56]"
_ISCN_RE = re.compile(
    r"(?:arr|seq)\s*\[\s*(?P<build>" + _BUILD_WORDS + r")\s*\]\s*"
    r"(?P<band>[0-9XY]{1,2}[pq][\d.]*(?:[pq][\d.]*)?)?\s*"
    r"\(\s*(?P<start>" + _NUM + r")\s*[_-]\s*(?P<end>" + _NUM + r")\s*\)\s*x\s*(?P<cn>\d+)"
    r"(?:\s*~\s*(?P<cn2>\d+))?(?P<cn_rest>\s*-\s*\d+)?(?P<tail>.*)$",
    re.I)
# the separator never includes "." — "chr22:11.21-11.23" is a band pair, not coordinates
_REGION_RE = re.compile(
    r"^(?:chr(?P<chrom1>[0-9]{1,2}|X|Y|MT|M)[:\s]|(?P<chrom2>[0-9]{1,2}|X|Y|MT|M):)\s*(?P<start>" + _NUM + r")"
    r"\s*(?:-|--|_|\.\.)\s*"
    r"(?P<end>" + _NUM + r")(?P<rest>.*)$", re.I)
_EXON_RE = re.compile(
    r"^(?P<gene>[A-Za-z][A-Za-z0-9._-]{0,19})\s*(?::|\s)\s*(?:exons?|ex)\s*\.?\s*"
    r"(?P<first>\d{1,3}[ab]?)(?!\d)(?:\s*(?:[-_]|to)\s*(?P<last>\d{1,3}[ab]?)(?!\d))?\s*"
    r"(?P<zyg>homozygous|heterozygous|hemizygous|hom|het|hemi)?\s*"
    r"(?P<type>deletion|deleted|del|duplication|duplicated|dup|loss|gain)?\s*$", re.I)
_SMN_NO_EXON_RE = re.compile(
    r"^(?P<gene>SMN[12])\s*(?::|\s)\s*(?P<zyg>homozygous|heterozygous|hom|het)\s*"
    r"(?P<type>deletion|deleted|del|duplication|dup|loss|gain)\s*$", re.I)
_HGVS_EXON_RE = re.compile(
    r"^(?P<ref>[A-Za-z0-9_.]+)(?:\((?P<gene>[A-Za-z0-9_.-]+)\))?:c\.(?P<start>\d+)(?P<soff>[+-](?:\?|\d+))?_"
    r"(?P<end>\d+)(?P<eoff>[+-](?:\?|\d+))?(?P<type>del|dup)$", re.I)
# the HGVS uncertain-breakpoint form: c.(6438+1_6439-1)_(7309+1_7310-1)del
_HGVS_RANGE_RE = re.compile(
    r"^(?P<ref>[A-Za-z0-9_.]+)(?:\((?P<gene>[A-Za-z0-9_.-]+)\))?:c\.\(\s*(?P<a>\d+)\+1\s*_\s*(?P<b>\d+)-1\s*\)_"
    r"\(\s*(?P<c>\d+)\+1\s*_\s*(?P<d>\d+)-1\s*\)(?P<type>del|dup)$", re.I)
_COPY_RE = re.compile(
    r"^(?P<gene>[A-Za-z][A-Za-z0-9._-]{0,19})\s*(?::|\s)\s*"
    r"(?:(?:exon\s*(?P<exon>\d[\dab]*(?:\s*(?:,|and|&|-)\s*\d[\dab]*)*)\s*[:=]?\s*)?"
    r"(?:copy\s*number|copies|cn|x)\s*[:=]?\s*(?P<cn>\d+)"
    r"|(?:exon\s*(?P<exon2>\d[\dab]*)\s*[:=]?\s*)?(?P<cn2>\d+)\s*cop(?:y|ies))\s*$", re.I)
_COUNT = r"(?:[<>]=?|≥|≤)?\s*\d{1,5}"
_REPEAT_RE = re.compile(
    r"^(?P<gene>[A-Za-z][A-Za-z0-9._-]{0,19})\s*(?::|\s)\s*\(?(?P<motif>[ACGT]{2,10}|G4C2)\)?n?\s*"
    r"(?:repeats?\s*[:=]?\s*)?(?P<count>" + _COUNT + r"(?:\s*(?:[-_/]|,)\s*" + _COUNT + r")?)(?![\d,])\s*"
    r"(?:repeats?|units?)?\s*$", re.I)
_TRANSCRIPT_RE = re.compile(r"^(?:N[MR]_\d+(?:\.\d+)?|ENST\d+(?:\.\d+)?|LRG_\d+t\d+)$", re.I)
# only the report's own words; x-notation goes through _type_from_copies, which
# knows that one copy of X is normal in a male
_TYPE_WORDS = {"del": "loss", "deletion": "loss", "loss": "loss", "deleted": "loss",
               "dup": "gain", "duplication": "gain", "gain": "gain", "amplification": "gain",
               "duplicated": "gain", "triplication": "gain"}
_BUILDS = {"grch38": "GRCh38", "hg38": "GRCh38", "grch37": "GRCh37", "hg19": "GRCh37"}
_OLD_BUILDS = ("hg18", "ncbi36", "ncbi35")
# SMN1/SMN2 legacy clinical exon labels -> Ensembl ordinal exon on the 9-exon
# canonical transcript. GeneReviews NBK1352: "SMN1 and SMN2 each comprise nine
# exons"; PMC11275604: "Both genes include nine exons (1, 2a, 2b, 3–8)". Checked
# live 2026-10-06: ENST00000380707 has 9 exons and its 54 bp exon is ordinal 8.
SMN_LEGACY_TO_ORDINAL = {"1": 1, "2a": 2, "2b": 3, "3": 4, "4": 5, "5": 6, "6": 7, "7": 8, "8": 9}

FORMS = (
    'a CNV by coordinates:      "chr15:23123715-28193120 loss" (add --copies 1)',
    'a CNV as ISCN:             "arr[GRCh38] 15q11.2q13.1(23123715_28193120)x1"',
    'an exon deletion/dup:      "DMD exon 45-50 deletion"',
    'the same in HGVS:          "NM_004006.3:c.6439-?_7309+?del" or "NM_004006.3:c.(6438+1_6439-1)_(7309+1_7310-1)del"',
    'a copy-number result:      "SMN1 exon 7 copy number 0" or "SMN1 0 copies, SMN2 3 copies"',
    'a repeat expansion:        "FMR1 CGG 230" (two alleles: "FMR1 CGG 30/230")',
)


def _int(text: str) -> int:
    t = (text or "").strip()
    if "," in t and not re.match(r"^\d{1,3}(?:,\d{3})+$", t):
        raise UsageError(f"{text!r} is not a number zebra can read: thousands separators must group three digits "
                         "(18,648,855)")
    digits = t.replace(",", "")
    if not digits.isdigit():
        raise UsageError(f"{text!r} is not a coordinate")
    return int(digits)


def _normalise(text: str) -> str:
    """The typography reports and PDFs use, mapped to what the patterns read (B-P2-1)."""
    out = (text or "").replace("×", "x").replace("✕", "x")  # multiplication sign
    for dash in ("–", "—", "−", "‐", "‑"):
        out = out.replace(dash, "-")
    for wide, narrow in (("：", ":"), ("，", ","), ("；", ";"), ("（", "("), ("）", ")"),
                         (" ", " "), ("　", " ")):
        out = out.replace(wide, narrow)
    return " ".join(out.split())


def _sex_norm(sex: Optional[str]) -> Optional[str]:
    s = (sex or "").strip().lower()
    return {"male": "male", "m": "male", "female": "female", "f": "female"}.get(s)


def _in_par(chrom: Optional[str], start: Optional[int], end: Optional[int], assembly: str) -> bool:
    if chrom != "X" or start is None or end is None:
        return False
    return any(s <= start and end <= e for s, e in PAR_X.get(assembly, []))


def _normal_copies(chrom: Optional[str], sex: Optional[str], par: bool = False) -> Optional[int]:
    """The copy number that is normal here: 2 on an autosome, by sex on X and Y; None when it depends on sex."""
    if chrom in ("X", "Y") and not par:
        if sex is None:
            return None
        if chrom == "X":
            return 1 if sex == "male" else 2
        return 1 if sex == "male" else 0
    return 2


def _type_from_copies(chrom: Optional[str], copies: Optional[int], sex: Optional[str] = None,
                      par: bool = False) -> Tuple[str, Optional[str]]:
    """Loss or gain from a copy number — which on X and Y depends on the proband's sex."""
    if copies is None:
        return "unknown", None
    if chrom == "MT":
        return "unknown", ("a copy number on mtDNA is not a loss or gain by itself: mtDNA findings are read as "
                           "heteroplasmy levels; use the word the report prints")
    normal = _normal_copies(chrom, sex, par)
    if normal is None:
        if copies == 0:
            return "loss", None
        if copies >= 3:
            return "gain", None
        return "unknown", (f"{copies} cop{'y' if copies == 1 else 'ies'} of {chrom} is not a loss or a gain by "
                           f"itself: one copy is normal in a male and a loss in a female, two copies the other way "
                           "round. Use the word the report prints (loss/gain/deletion/duplication), or give the "
                           "proband's sex (--sex, or the case profile).")
    if chrom == "Y" and sex == "female" and not par:
        return "unknown", f"{copies} cop{'y' if copies == 1 else 'ies'} of Y in a proband recorded as female"
    if copies < normal:
        return "loss", None
    if copies > normal:
        return "gain", None
    where = "autosomal" if chrom not in ("X", "Y") or par else f"{chrom} in a {sex}"
    return "unknown", (f"{copies} copies is the normal {where} state: the report must say what it called abnormal")


def _possible_types(chrom: Optional[str], low: int, high: int, sex: Optional[str], par: bool) -> set:
    """Every reading (loss / gain / normal / spans) a copy number or range can have here, over the sexes still possible."""
    if chrom == "MT":
        return {"loss", "gain", "normal", "spans"}
    sexes = [sex] if (sex or chrom not in ("X", "Y") or par) else ["male", "female"]
    out = set()
    for one in sexes:
        normal = _normal_copies(chrom, one, par)
        if normal is None:
            out |= {"loss", "gain", "normal"}
        elif high < normal:
            out.add("loss")
        elif low > normal:
            out.add("gain")
        elif low == high == normal:
            out.add("normal")
        elif high <= normal:
            out.add("loss")  # a mosaic range reaching down from the normal state
        elif low >= normal:
            out.add("gain")
        else:
            out.add("spans")
    return out


def _type_from_range(chrom: Optional[str], low: int, high: int, sex: Optional[str],
                     par: bool = False) -> Tuple[str, Optional[str]]:
    """A mosaic copy-number range (x1~2, x2~3): a mosaic loss or gain relative to the normal state."""
    normal = _normal_copies(chrom, sex, par) if chrom != "MT" else None
    if chrom == "MT":
        return "unknown", ("a copy-number range on mtDNA is not a loss or gain by itself: mtDNA findings are read as "
                           "heteroplasmy levels; use the word the report prints")
    if normal is None:
        return "unknown", (f"a copy-number range of {low}~{high} on {chrom} cannot be read as a loss or a gain "
                           "without the proband's sex (--sex, or the case profile)")
    if high <= normal and low < normal:
        return "loss", None
    if low >= normal and high > normal:
        return "gain", None
    return "unknown", f"the range {low}~{high} spans the normal state ({normal} copies)"


# ------------------------------------------------- the words after the interval

_TAIL_PATTERNS: Tuple[Tuple[str, str, Any], ...] = (
    (r"\b(?:de\s*novo|dn)\b", "inheritance", "de_novo"),
    (r"\b(?:maternal(?:ly)?(?:\s+inherited)?|mat)\b", "inheritance", "maternal"),
    (r"\b(?:paternal(?:ly)?(?:\s+inherited)?|pat)\b", "inheritance", "paternal"),
    (r"\b(?:inh|inherited)\b", "inheritance_note", "inherited (the report does not say from which parent)"),
    (r"\blikely\s+pathogenic\b", "classification_lab", "likely pathogenic"),
    (r"\blikely\s+benign\b", "classification_lab", "likely benign"),
    (r"\b(?:vus|uncertain\s+significance|variant\s+of\s+uncertain\s+significance)\b", "classification_lab",
     "uncertain significance"),
    (r"\bpathogenic\b", "classification_lab", "pathogenic"),
    (r"\bbenign\b", "classification_lab", "benign"),
    (r"\b(?:homozygous|hom)\b", "zygosity", "homozygous"),
    (r"\b(?:heterozygous|het)\b", "zygosity", "heterozygous"),
    (r"\b(?:hemizygous|hemi)\b", "zygosity", "hemizygous"),
    (r"\binterstitial\b", "extent", "interstitial"),
    (r"\bterminal\b", "extent", "terminal"),
)
_FILLER = {"of", "the", "a", "an", "region", "cnv", "call", "on", "in", "with", "and", "copy", "number"}


def _read_tail(tail: str, raw: str) -> Dict[str, Any]:
    """Everything a report writes after the interval, read token by token; anything unread is refused (B-P1-4).

    Returns the fields found: assembly, type_word, copy_number, copy_range (low, high),
    mosaic, mosaic_fraction, inheritance, zygosity, classification_lab, extent.
    """
    got: Dict[str, Any] = {}
    low = " " + tail.lower() + " "
    if re.search(r"(?:\bchr|\b)(?:[0-9]{1,2}|x|y|mt?)\s*:\s*[\d,]{3,}\s*(?:-|_|\.\.)", low) or \
            re.search(r"\b(?:arr|seq)\s*\[", low):
        raise UsageError(f"{raw!r} holds more than one finding (a second interval follows the first). Give one "
                         "finding at a time, exactly as the report prints it, and run `zebra cnv` once for each.")
    builds = re.findall(r"\b(" + _BUILD_WORDS.lower() + r")\b", low)
    if builds:
        names = {b for b in builds}
        if any(b in _OLD_BUILDS for b in names):
            raise UsageError(f"{raw!r} is on {', '.join(sorted(names))} (NCBI36/hg18), a build zebra does not "
                             "read: lift the interval over to GRCh37 or GRCh38 first (the report's laboratory or "
                             "the UCSC liftOver tool can do this)")
        mapped = {_BUILDS[b] for b in names}
        if len(mapped) > 1:
            raise UsageError(f"{raw!r} names two builds ({', '.join(sorted(mapped))}): which one are the "
                             "coordinates on?")
        got["assembly"] = mapped.pop()
        low = re.sub(r"\b(?:" + _BUILD_WORDS.lower() + r")\b", " ", low)
    frac = re.search(r"(\d{1,3}(?:\.\d+)?)\s*%", low) or re.search(r"\[\s*(0?\.\d+)\s*\]", low)
    if frac:
        got["mosaic"] = True
        got["mosaic_fraction"] = frac.group(1) + ("%" if "%" in frac.group(0) else "")
        low = low.replace(frac.group(0), " ")
    if re.search(r"\bmos(?:aic(?:ism)?)?\b", low):
        got["mosaic"] = True
        low = re.sub(r"\bmos(?:aic(?:ism)?)?\b", " ", low)
    cns: List[Tuple[int, Optional[int]]] = []
    for pat in (r"\bx\s*(\d+)(?:\s*~\s*(\d+))?\b", r"\bcn\s*[:=]?\s*(\d+)\b()", r"\bcopy\s*number\s*[:=]?\s*(\d+)\b()",
                r"\b(\d+)\s*cop(?:y|ies)\b()"):
        for m in re.finditer(pat, low):
            cns.append((int(m.group(1)), int(m.group(2)) if m.group(2) else None))
        low = re.sub(pat, " ", low)
    if re.search(r"~", low):
        raise UsageError(f"{raw!r}: a '~' outside an x-notation range (x1~2) is not something zebra can read")
    if len(set(cns)) > 1:
        raise UsageError(f"{raw!r} gives more than one copy number ({', '.join(str(c[0]) for c in cns)})")
    if cns:
        cn, cn2 = cns[0]
        if cn2 is not None:
            got["copy_range"] = (min(cn, cn2), max(cn, cn2))
            got["mosaic"] = True
        else:
            got["copy_number"] = cn
    for word, mapped in sorted(_TYPE_WORDS.items(), key=lambda kv: -len(kv[0])):
        if re.search(rf"\b{word}\b", low):
            if got.get("type_word") and got["type_word"] != mapped:
                raise UsageError(f"{raw!r} says both loss and gain")
            got["type_word"] = mapped
            low = re.sub(rf"\b{word}\b", " ", low)
    for pat, key, value in _TAIL_PATTERNS:
        if re.search(pat, low):
            if got.get(key) and got[key] != value:
                raise UsageError(f"{raw!r} gives two different values for {key} ({got[key]}, {value})")
            got[key] = value
            low = re.sub(pat, " ", low)
    leftover = [w for w in re.split(r"[\s,;:()\[\]=/]+", low) if w and w not in _FILLER]
    if leftover:
        raise UsageError(
            f"{raw!r} carries {' '.join(leftover)!r}, which zebra cannot read, and it will not ignore it. After the "
            "interval zebra reads: loss/gain/deletion/duplication, a copy number (x1, cn=1, 1 copy, x1~2), a build "
            "(hg19/hg38/GRCh37/GRCh38), mosaic/mos/30%, dn/mat/pat, het/hom/hemi, the lab's classification. If it is "
            "a second finding, give one finding at a time.")
    return got


def _cnv_fields(raw: str, chrom: str, start: int, end: int, cn: Optional[int], tail: Dict[str, Any],
                iscn: Optional[str], band: Optional[str], assembly: Optional[str],
                sex: Optional[str], default_assembly: str = "GRCh38") -> Dict[str, Any]:
    build = assembly or tail.get("assembly") or default_assembly
    par = _in_par(chrom, start, end, build)
    copy_range = tail.get("copy_range")
    if copy_range and copy_range[0] == copy_range[1]:
        raise UsageError(f"{raw!r}: x{copy_range[0]}~{copy_range[1]} is not a range")
    frac = tail.get("mosaic_fraction")
    if frac:
        value = float(frac.rstrip("%")) / (100.0 if frac.endswith("%") else 1.0)
        if not 0 < value < 1:
            raise UsageError(f"{raw!r}: a mosaic fraction of {frac} is not a fraction of cells (it must lie strictly "
                             "between 0 and 100%)")
    zyg = tail.get("zygosity")
    if zyg and cn is not None and chrom not in ("X", "Y", "MT"):
        # on an autosome a homozygous loss is 0 copies and a heterozygous one is 1
        if (zyg == "homozygous" and cn == 1) or (zyg == "heterozygous" and cn == 0) or zyg == "hemizygous":
            raise UsageError(f"{raw!r}: '{zyg}' and copy number {cn} on chromosome {chrom} contradict each other")
    note = None
    if copy_range:
        kind, note = _type_from_range(chrom, copy_range[0], copy_range[1], sex, par)
    else:
        kind, note = _type_from_copies(chrom, cn, sex, par) if cn is not None else ("unknown", None)
    word = tail.get("type_word")
    if word and chrom != "MT" and (cn is not None or copy_range):
        lo, hi = (copy_range if copy_range else (cn, cn))
        if word not in _possible_types(chrom, lo, hi, sex, par):
            shown = cn if cn is not None else f"{lo}~{hi}"
            raise UsageError(f"the report says {word} and copy number {shown} on {chrom} cannot be a {word} "
                             f"({'in either sex' if chrom in ('X', 'Y') and not sex and not par else 'here'}): one of "
                             "the two is wrong, so zebra will not record either")
    if word:
        if kind in ("loss", "gain") and kind != word:
            raise UsageError(f"the report says {word} and the copy number says {kind} (copy number "
                             f"{cn if cn is not None else '~'.join(map(str, copy_range or ()))} on {chrom}): one of "
                             "the two is wrong, so zebra will not record either")
        kind, note = word, None  # the report's own word always wins over an ambiguous count
    return {"kind": "cnv", "input": raw, "iscn": iscn, "chrom": chrom, "band": band, "start": start, "end": end,
            "copy_number": cn, "copy_number_range": f"{copy_range[0]}~{copy_range[1]}" if copy_range else None,
            "mosaic": bool(tail.get("mosaic")), "mosaic_fraction": tail.get("mosaic_fraction"),
            "cnv_type": kind or "unknown", "cnv_type_note": note, "assembly": assembly or tail.get("assembly"),
            "inheritance": tail.get("inheritance"), "inheritance_note": tail.get("inheritance_note"),
            "zygosity": tail.get("zygosity"), "classification_lab": tail.get("classification_lab"),
            "extent": tail.get("extent"), "in_par": par}


def _smn_combo(raw: str) -> Optional[Dict[str, Any]]:
    """"SMN1 0 copies, SMN2 3 copies" (or with ';', 'and', '/'): one SMN1 result with its SMN2 modifier."""
    if not re.search(r"\bSMN1\b", raw, re.I) or not re.search(r"\bSMN2\b", raw, re.I):
        return None
    parts = [p.strip() for p in re.split(r"\s*(?:[,;/]|\band\b|\bwith\b)\s*", raw, flags=re.I) if p.strip()]
    if len(parts) < 2:
        return None
    found: Dict[str, Dict[str, Any]] = {}
    for part in parts:
        m = _COPY_RE.match(part)
        if not m:
            return None
        gene = m.group("gene").upper()
        if gene not in ("SMN1", "SMN2"):
            return None
        cn = int(m.group("cn") if m.group("cn") is not None else m.group("cn2"))
        exon = (m.group("exon") or m.group("exon2") or "").strip() or None
        if gene in found and found[gene]["copy_number"] != cn:
            raise UsageError(f"{raw!r} gives two copy numbers for {gene}")
        found[gene] = {"gene": gene, "copy_number": cn, "exon": exon}
    if "SMN1" not in found:
        return None
    smn1 = found["SMN1"]
    out = {"kind": "copy_number", "input": raw, "gene": "SMN1", "exon": smn1["exon"], "copy_number": smn1["copy_number"]}
    if "SMN2" in found:
        out["related_parsed"] = [found["SMN2"]]
    return out


def _hgvs_cds(start: int, soff: Optional[str], end: int, eoff: Optional[str]) -> Tuple[int, int]:
    """The first and last *coding* bases an exon-boundary HGVS deletion removes (B-P0-1).

    An intronic position is numbered from the nearest exon, so the sign of the
    offset says which side of the exon the breakpoint lies on: a start written
    `c.N+k` lies in the intron after N (the first lost coding base is N+1); an end
    written `c.M-k` lies in the intron before M (the last lost coding base is M-1).
    `c.N-k` / `c.M+k` already name the first/last lost coding base. The same holds
    for '?' offsets.
    """
    first = start + 1 if soff and soff.startswith("+") else start
    last = end - 1 if eoff and eoff.startswith("-") else end
    return first, last


_CHROM_TOKEN = re.compile(r"^(?:chr)?(?:[0-9]{1,2}|X|Y|MT?)$", re.I)


def _not_a_chromosome(token: str, raw: str) -> None:
    if _CHROM_TOKEN.match(token or ""):
        raise UsageError(f"{raw!r} names a chromosome where a gene belongs: a whole-chromosome copy number "
                         "(aneuploidy) is not read here. Give the interval the report prints "
                         "(\"chr21:1-46709983 x3\") or record the karyotype as the report writes it.")


def parse(text: str, sex: Optional[str] = None, assembly: str = "GRCh38") -> Dict[str, Any]:
    """Which result form this is, and the fields it carries. Nothing is guessed.

    `assembly` is the build to assume when the result does not name one (it decides
    where the X pseudoautosomal regions lie)."""
    raw = _normalise(text)
    if not raw:
        raise UsageError("give a result to record. Accepted forms:\n  " + "\n  ".join(FORMS))
    if len(raw) > INPUT_LIMIT:
        raise UsageError(f"the result is {len(raw)} characters long; give one finding as the report prints it "
                         f"(at most {INPUT_LIMIT} characters)")
    low = raw.lower()
    sex = _sex_norm(sex)

    if len(re.findall(r"\b(?:arr|seq)\s*\[", raw, re.I)) > 1:
        raise UsageError(f"{raw!r} holds more than one finding (several ISCN results). Give one finding at a time, "
                         "exactly as the report prints it, and run `zebra cnv` once for each.")
    combo = _smn_combo(raw)
    if combo:
        return combo

    iscn = _ISCN_RE.match(raw)
    if iscn:
        build_word = iscn.group("build").lower()
        if build_word in _OLD_BUILDS:
            raise UsageError(f"{raw!r} is on {iscn.group('build')} (NCBI36/hg18), a build zebra does not read: lift "
                             "the interval over to GRCh37 or GRCh38 first")
        if iscn.group("cn_rest"):
            raise UsageError(f"{raw!r} gives a copy number followed by {iscn.group('cn_rest').strip()!r}; ISCN "
                             "writes a mosaic copy-number range with '~' (x1~2). Give the result as the report "
                             "prints it.")
        tail = _read_tail(iscn.group("tail") or "", raw)
        if tail.get("assembly") and tail["assembly"] != _BUILDS[build_word]:
            raise UsageError(f"{raw!r} names {iscn.group('build')} in the brackets and {tail['assembly']} after it")
        if tail.get("copy_number") is not None or tail.get("copy_range"):
            raise UsageError(f"{raw!r} gives a second copy number after the ISCN x-notation")
        cn = int(iscn.group("cn"))
        if iscn.group("cn2") is not None:
            tail["copy_range"] = (min(cn, int(iscn.group("cn2"))), max(cn, int(iscn.group("cn2"))))
            tail["mosaic"] = True
        band = iscn.group("band")
        chrom = (band or "")[:2].upper().rstrip("PQ") or None
        fields = _cnv_fields(raw, chrom or "", _int(iscn.group("start")), _int(iscn.group("end")),
                             None if tail.get("copy_range") else cn, tail, raw, band,
                             _BUILDS[build_word], sex, assembly)
        fields["chrom"] = chrom
        return fields

    region = _REGION_RE.match(raw)
    if region and not _EXON_RE.match(raw):
        chrom = (region.group("chrom1") or region.group("chrom2")).upper()
        chrom = "MT" if chrom == "M" else chrom
        leftover = (region.group("rest") or "").strip()
        start, end = _int(region.group("start")), _int(region.group("end"))
        if re.match(r"^[-_.]?\d", leftover):
            raise UsageError(
                f"{raw!r} does not read as one interval: {leftover!r} is left over after "
                f"{chrom}:{start}-{end}. A cytogenetic band pair such as 22q11.21-q11.23 is not a base-pair "
                "interval; give the coordinates the report prints.")
        tail = _read_tail(leftover, raw)
        return _cnv_fields(raw, chrom, start, end, tail.get("copy_number"), tail, None, None, None, sex, assembly)

    hgvs = _HGVS_EXON_RE.match(raw)
    rng = _HGVS_RANGE_RE.match(raw)
    if hgvs or rng:
        if rng:
            a, b, c, d = (int(rng.group(k)) for k in "abcd")
            if not (b == a + 1 and d == c + 1):
                raise UsageError(f"{raw}: each bracket must name one intron — the last coding base before it and the "
                                 "first after it, e.g. (6438+1_6439-1). A bracket spanning exons leaves the deleted "
                                 "exons uncertain, so zebra gives no frame for it; read the exons from the report.")
            cds_start, cds_end = b, c
            ref, gene_in, kind = rng.group("ref"), rng.group("gene"), rng.group("type")
            notation = "uncertain-breakpoint ranges"
        else:
            start, end = int(hgvs.group("start")), int(hgvs.group("end"))
            soff, eoff = hgvs.group("soff"), hgvs.group("eoff")
            if end < start:
                raise UsageError(f"{raw}: the second c. position is before the first")
            if not soff and not eoff and end - start + 1 < SMALL_VARIANT_BP:
                raise UsageError(f"{raw} describes a {end - start + 1} bp change with exact breakpoints: that is a "
                                 "small variant, not an exon-level event. Use `zebra variant` for it (this command is "
                                 "for exon-boundary descriptions, whose breakpoints lie in the introns).")
            cds_start, cds_end = _hgvs_cds(start, soff, end, eoff)
            ref, gene_in, kind = hgvs.group("ref"), hgvs.group("gene"), hgvs.group("type")
            notation = ("unsequenced intronic breakpoints ('?')" if "?" in raw else
                        "intronic offsets" if (soff or eoff) else "coding positions")
        if cds_end < cds_start:
            raise UsageError(f"{raw} removes no coding base: both breakpoints lie in the same intron. Use "
                             "`zebra variant` (or `zebra s2f` for splicing) for an intronic change.")
        return {"kind": "exon_cnv", "input": raw, "gene": gene_in.upper() if gene_in else None, "transcript": ref,
                "first": None, "last": None, "cds_start": cds_start, "cds_end": cds_end,
                "cnv_type": _TYPE_WORDS[kind.lower()], "hgvs_notation": notation, "zygosity": None}

    exon = _EXON_RE.match(raw)
    smn_bare = _SMN_NO_EXON_RE.match(raw)
    if exon or smn_bare:
        m = exon or smn_bare
        token = m.group("gene")
        _not_a_chromosome(token, raw)
        gene = token.upper()
        word = (m.group("type") or "").lower()
        zyg = (m.group("zyg") or "").lower() or None
        zyg = {"hom": "homozygous", "het": "heterozygous", "hemi": "hemizygous"}.get(zyg or "", zyg)
        if gene in ("SMN1", "SMN2"):
            return _smn_exon_form(raw, gene, m.group("first") if exon else None,
                                  (m.group("last") if exon else None), _TYPE_WORDS.get(word, "unknown"), zyg)
        first_t, last_t = m.group("first"), m.group("last") or m.group("first")
        if re.search(r"[ab]$", first_t + last_t, re.I):
            raise UsageError(f"{raw}: exon labels such as 2a/2b are legacy numbering zebra knows only for SMN1 and "
                             "SMN2; give the exon numbers of the transcript the report names")
        first, last = int(first_t), int(last_t)
        if last < first:
            raise UsageError(f"{raw}: exon {last} is before exon {first}")
        if first < 1:
            raise UsageError(f"{raw}: exons are numbered from 1")
        is_transcript = bool(_TRANSCRIPT_RE.match(token))
        return {"kind": "exon_cnv", "input": raw, "gene": None if is_transcript else gene,
                "transcript": token if is_transcript else None,
                "first": first, "last": last, "cds_start": None, "cds_end": None,
                "cnv_type": _TYPE_WORDS.get(word, "unknown"), "zygosity": zyg}

    copies = _COPY_RE.match(raw)
    if copies and not _REPEAT_RE.match(raw):
        _not_a_chromosome(copies.group("gene"), raw)
        cn = int(copies.group("cn") if copies.group("cn") is not None else copies.group("cn2"))
        gene = copies.group("gene").upper()
        exon = (copies.group("exon") or copies.group("exon2") or "").strip() or None
        if gene in ("SMN1", "SMN2") and exon:
            for label in re.split(r"\s*(?:,|and|&|-)\s*", exon.lower()):
                if label and label not in SMN_LEGACY_TO_ORDINAL:
                    raise UsageError(f"{raw}: SMN1/SMN2 exons are numbered 1, 2a, 2b, 3-8 in the clinical numbering; "
                                     f"there is no exon {label}" + (" (exon 2 is split into 2a and 2b)"
                                                                   if label == "2" else ""))
        return {"kind": "copy_number", "input": raw, "gene": gene, "exon": exon, "copy_number": cn}

    repeat = _REPEAT_RE.match(raw)
    if repeat:
        _not_a_chromosome(repeat.group("gene"), raw)
        motif = repeat.group("motif").upper()
        motif = "GGGGCC" if motif == "G4C2" else motif
        alleles, count = _repeat_alleles(repeat.group("count"), raw)
        return {"kind": "repeat_expansion", "input": raw, "gene": repeat.group("gene").upper(),
                "motif": motif, "repeat_count": count, "alleles": alleles}

    if re.search(r"\b(?:arr|seq)\s*\[", raw, re.I):
        raise UsageError(f"{raw!r} has text before the ISCN result (or an ISCN form zebra cannot read): give the "
                         "ISCN string exactly as the report prints it, starting with arr[...] or seq[...]")
    if re.search(r"[,;]", raw) and (re.search(r"\bSMN[12]\b", raw, re.I) or _REGION_RE.match(raw.split(",")[0].strip())):
        raise UsageError(f"{raw!r} holds more than one finding, and zebra reads the parts together only for "
                         "SMN1 with SMN2 (\"SMN1 0 copies, SMN2 3 copies\"). Give one finding at a time.")
    if re.match(r"^(del|dup)\s*\(", low) or re.search(r"\b[0-9xy]{1,2}[pq][\d.]+\b", low) \
            or re.search(r"[:\s]\d{1,2}\.\d", low):
        raise UsageError(
            f"{raw!r} names a band but no coordinates, and zebra will not guess them: cytogenetic bands are not "
            "base-pair boundaries. Give the interval the report prints, e.g. "
            '"chr15:23123715-28193120 loss" or the full ISCN string with coordinates.')
    raise UsageError(f"cannot read {raw!r}. Accepted forms:\n  " + "\n  ".join(FORMS))


def _smn_exon_form(raw: str, gene: str, first: Optional[str], last: Optional[str], kind: str,
                   zyg: Optional[str]) -> Dict[str, Any]:
    """"SMN1 exon 7 deletion" is a dosage result in legacy exon numbering, not an exon-skipping question (B-P1-2)."""
    labels = []
    for label in (first, last):
        if label is None:
            continue
        label = label.lower()
        if label == "2":
            raise UsageError(f"{raw}: SMN1/SMN2 exon 2 is split into 2a and 2b in the clinical numbering; say which")
        if label not in SMN_LEGACY_TO_ORDINAL:
            raise UsageError(f"{raw}: SMN1/SMN2 exons are numbered 1, 2a, 2b, 3-8 in the clinical numbering; "
                             f"there is no exon {label}")
        if label not in labels:
            labels.append(label)
    copy_number: Optional[int] = None
    note = None
    if kind == "loss":
        if zyg == "homozygous":
            copy_number = 0
        elif zyg == "heterozygous":
            copy_number = 1
        else:
            note = ("the result says 'deletion' without saying whether one copy or both are missing: a homozygous "
                    "deletion is 0 copies (the molecular finding of SMA), a heterozygous one is 1 copy (a carrier). "
                    "Give the copy number the report prints (\"SMN1 exon 7 copy number 0\").")
    elif kind == "gain":
        note = "a duplication means more than two copies; give the copy number the report prints"
    else:
        note = "the result does not say deletion or duplication; give the copy number the report prints"
    return {"kind": "copy_number", "input": raw, "gene": gene, "exon": "-".join(labels) if labels else None,
            "copy_number": copy_number, "copy_number_note": note, "zygosity": zyg, "report_word": kind,
            "exon_numbering": "legacy clinical (1, 2a, 2b, 3-8)"}


def _repeat_alleles(text: str, raw: str) -> Tuple[List[Dict[str, Any]], Any]:
    """The repeat count as written: one number, a sizing range ('55-200'), two alleles ('30/230'), or a bound."""
    if re.search(r"\d,\d", text):
        raise UsageError(f"{raw!r}: '{text}' could be one number with a thousands comma or two alleles. Write 1000 "
                         "for one allele, or 30/230 (or '30, 230') for the two alleles of a female")
    t = text.replace("≥", ">=").replace("≤", "<=").replace(" ", "")
    two = re.match(r"^(\d+)[/,](\d+)$", t)
    if two:
        a, b = int(two.group(1)), int(two.group(2))
        return [{"allele": 1, "low": a, "high": a, "as_written": str(a)},
                {"allele": 2, "low": b, "high": b, "as_written": str(b)}], f"{a}/{b}"
    span = re.match(r"^(\d+)[-_](\d+)$", t)
    if span:
        a, b = int(span.group(1)), int(span.group(2))
        if b < a:
            raise UsageError(f"{raw}: the repeat range runs backwards ({a}-{b})")
        return [{"allele": 1, "low": a, "high": b, "as_written": f"{a}-{b}"}], f"{a}-{b}"
    bound = re.match(r"^(>=|<=|>|<)(\d+)$", t)
    if bound:
        op, n = bound.group(1), int(bound.group(2))
        low = n + 1 if op == ">" else n if op == ">=" else None
        high = n - 1 if op == "<" else n if op == "<=" else None
        return [{"allele": 1, "low": low, "high": high, "as_written": f"{op}{n}"}], f"{op}{n}"
    one = re.match(r"^(\d+)$", t)
    if one:
        n = int(one.group(1))
        return [{"allele": 1, "low": n, "high": n, "as_written": str(n)}], n
    raise UsageError(f"cannot read the repeat count in {raw!r}: give one number, a range (55-200), two alleles "
                     "(30/230) or a bound (>200)")


# ------------------------------------------------------------------- CNV

OVERLAP_WINDOW = 4_500_000  # Ensembl /overlap/region refuses more than 5 Mb per request
MAX_REGION_BP = 50_000_000  # beyond this the gene list is not walked one window at a time (summary only)


def _overlap_problem(text: str) -> Optional[str]:
    """Why an Ensembl /overlap body is unusable (an error object, HTML, a cut-off list), or None."""
    body = (text or "").lstrip()
    if not body.startswith("["):
        return "not a JSON list of genes (an error page or an error object)"
    try:
        data = json.loads(body)
    except ValueError:
        return "a truncated or malformed JSON list"
    return None if isinstance(data, list) else "not a JSON list"


def genes_in_region(chrom: str, start: int, end: int, assembly: str = "GRCh38") -> Outcome:
    """Every gene Ensembl places in [start, end] on `chrom` (in windows, because the endpoint caps at 5 Mb).

    When a window fails (the source, or the time budget), the genes of the windows
    already read are kept: `query["complete"]` is False and `query["covered_to"]`
    says how far the list reaches; the failure is a named warning."""
    if end - start + 1 > MAX_REGION_BP:
        raise UsageError(f"{chrom}:{start}-{end} is {(end - start + 1) / 1e6:.1f} Mb; zebra will not walk more than "
                         f"{MAX_REGION_BP / 1e6:.0f} Mb of Ensembl gene annotation one window at a time. Check the "
                         "coordinates (a digit too many is the usual cause), or query the genes you care about with "
                         "`zebra gene`.")
    by_id: Dict[str, Dict[str, Any]] = {}
    sources: List[Dict[str, Any]] = []
    warnings: List[str] = []
    window_start = start
    while window_start <= end:
        window_end = min(end, window_start + OVERLAP_WINDOW - 1)
        region = f"{chrom}:{window_start}-{window_end}"
        try:
            resp = get_json(f"{ensembl.host(assembly)}/overlap/region/human/{region}", source="Ensembl overlap",
                            params={"feature": "gene"}, cache_ttl=30 * 86400, timeout=120, validate=_overlap_problem)
            data = resp.json()
        except (SourceError, ValueError) as err:
            if not sources:
                raise
            warnings.append(f"Ensembl overlap stopped at {chrom}:{window_start} ({getattr(err, 'message', err)}): the "
                            f"gene list covers {chrom}:{start}-{window_start - 1} only, so gene counts are a lower "
                            "bound")
            genes = sorted(by_id.values(), key=lambda g: (g["start"] or 0))
            return Outcome(genes, sources=sources, warnings=warnings,
                           query={"complete": False, "covered_to": window_start - 1})
        if not isinstance(data, list):
            raise ValueError("Ensembl overlap returned no list")
        for g in data:
            if not isinstance(g, dict) or not g.get("id"):
                continue
            by_id[g["id"]] = {"symbol": g.get("external_name") or g.get("id"), "ensembl_id": g.get("id"),
                              "biotype": g.get("biotype"), "start": g.get("start"), "end": g.get("end"),
                              "strand": g.get("strand"), "named": bool(g.get("external_name"))}
        sources.append(source_record("Ensembl overlap", f"{assembly} {region}", resp, note="feature=gene"))
        window_start = window_end + 1
    genes = sorted(by_id.values(), key=lambda g: (g["start"] or 0))
    return Outcome(genes, sources=sources, query={"complete": True, "covered_to": end})


def _clingen_tables(assembly: str, warnings: List[str], sources: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Both ClinGen dosage tables for `assembly`, parsed once; a table that fails is one named warning."""
    out: Dict[str, Any] = {"genes": None, "regions": None, "gene_file_date": None, "region_file_date": None}
    for kind, db in (("gene", "ClinGen dosage sensitivity"), ("region", "ClinGen dosage sensitivity (regions)")):
        resp = attempt(f"ClinGen dosage map, {kind} curations ({assembly})",
                       lambda k=kind: clingen.fetch_table(k, assembly), warnings)
        if resp is None:
            continue
        parse = clingen.parse_gene_table if kind == "gene" else clingen.parse_region_table
        rows = attempt(f"ClinGen dosage map, {kind} curations ({assembly})", lambda: parse(resp.text), warnings)
        if rows is None:
            continue
        date = clingen.file_date(resp.text)
        out["genes" if kind == "gene" else "regions"] = rows
        out[f"{kind}_file_date"] = date
        sources.append(source_record(db, f"{kind} curation list, {assembly}" + (f" (file dated {date})" if date else ""),
                                     resp, note=f"bulk TSV, {len(rows)} {kind} curations, filtered by location"))
    return out


def _score_text(block: Dict[str, Any]) -> str:
    score = block.get("score")
    if score is None:
        return "-"
    desc = block.get("description") or clingen.DOSAGE_SCORES.get(score)
    return f"{score} ({desc})" if desc else score


def _relevant_scores(kind: str) -> List[str]:
    return ["haploinsufficiency"] if kind == "loss" else ["triplosensitivity"] if kind == "gain" else \
        ["haploinsufficiency", "triplosensitivity"]


def _region_entry(r: Dict[str, Any], kind: str, gene_pos: Dict[str, Tuple[str, int, int]],
                  chrom: str, start: int, end: int, genes_complete: bool = False,
                  copies: Optional[int] = None) -> Dict[str, Any]:
    ov = r["overlap"]
    named = []
    for g in r.get("named_genes") or []:
        pos = gene_pos.get(g.upper())
        inside = None
        if pos and pos[0] == chrom:
            inside = "whole gene" if start <= pos[1] and end >= pos[2] else (
                "partly" if pos[1] <= end and pos[2] >= start else "no")
        elif genes_complete:
            inside = "no"  # Ensembl's complete gene list for the interval does not contain it
        named.append({"gene": g, "in_cnv": inside})
    entry = {"isca_id": r["isca_id"], "name": r["name"], "cytoband": r["cytoband"], "location": r["location"],
             "haploinsufficiency": {"score": r["haploinsufficiency"]["score"],
                                    "description": r["haploinsufficiency"]["description"]},
             "triplosensitivity": {"score": r["triplosensitivity"]["score"],
                                   "description": r["triplosensitivity"]["description"]},
             "relation": ov["relation"], "region_fraction_covered": ov["item_fraction_covered"],
             "cnv_fraction_in_region": ov["cnv_fraction_in_item"], "overlap_bp": ov["overlap_bp"],
             "last_evaluated": r["last_evaluated"], "url": r["url"]}
    if "item_bp_outside_cnv" in ov:
        entry["region_bp_outside_cnv"] = ov["item_bp_outside_cnv"]
    if named:
        entry["named_genes"] = named
    for which in ("haploinsufficiency", "triplosensitivity"):
        if r[which]["score"] in ("1", "2", "3") and r[which].get("pmids"):
            entry[which]["pmids"] = r[which]["pmids"]
    entry["acmg_section_2"] = _region_row(entry, kind, copies)
    return entry


def _region_row(entry: Dict[str, Any], kind: str, copies: Optional[int] = None) -> Optional[str]:
    """Which section-2 row of Riggs 2020 this overlap corresponds to, when the standard's own wording settles it."""
    hi, ts = entry["haploinsufficiency"]["score"], entry["triplosensitivity"]["score"]
    whole = entry["relation"] in ("cnv_contains_it", "identical")
    inside = entry["relation"] in ("cnv_within_it", "identical")
    named_inside = [g["gene"] for g in entry.get("named_genes") or [] if g.get("in_cnv") == "whole gene"]
    partial_text = ("partial overlap of an established {} genomic region: the standard's 2B row is for a CNV that does "
                    "NOT contain the causative gene or critical region" +
                    (f"; the named gene(s) {', '.join(named_inside)} lie inside this CNV, so a curator decides "
                     "whether the critical region is covered (2A) or not (2B)" if named_inside else
                     " — see named_genes"))
    if kind == "loss":
        if hi == "3":
            return ("2A row (complete overlap of an established HI genomic region)" if whole else
                    "2A or 2B row: " + partial_text.format("HI"))
        if hi == "40":
            if copies == 0:
                return ("2F/2G rows describe population (heterozygous) deletions: they do not speak to a homozygous "
                        "loss (copy number 0) — curator judgement")
            return ("2F row (completely contained within an established benign CNV region)" if inside else
                    "2G row (overlaps an established benign CNV but includes additional genomic material)")
    if kind == "gain":
        if ts == "3":
            return ("2A row (the TS region is fully contained within the gain)" if whole else
                    "2A or 2B row: " + partial_text.format("TS"))
        if ts == "40":
            return ("2C-2G rows (established benign copy-number gain): compare the gene content of the two "
                    "intervals")
    return None


def _gene_overlap_side(g: Dict[str, Any], start: int, end: int, strand: Optional[int]) -> str:
    gs, ge = g["start"], g["end"]
    if start <= gs and end >= ge:
        return "whole gene"
    if start > gs and end < ge:
        return "intragenic (both breakpoints inside the gene)"
    covers_low = start <= gs
    if strand not in (1, -1):
        return "partial (strand not known here)"
    five_prime_low = strand == 1
    return "partial, 5' end" if covers_low == five_prime_low else "partial, 3' end"


def _gene_records(cg: List[Dict[str, Any]], genes: Optional[List[Dict[str, Any]]], chrom: str,
                  start: int, end: int) -> Tuple[List[Dict[str, Any]], int]:
    """ClinGen gene curations that overlap the CNV by location, plus spanned genes matched by symbol."""
    strand = {g["symbol"].upper(): g.get("strand") for g in genes or [] if g.get("named")}
    hits = {g["gene"].upper(): g for g in clingen.genes_overlapping(cg, chrom, start, end)}
    by_symbol = {g["gene"].upper(): g for g in cg}
    for g in genes or []:
        sym = g["symbol"].upper()
        if g.get("named") and g.get("biotype") == "protein_coding" and sym not in hits and sym in by_symbol:
            row = by_symbol[sym]
            if row["start"] is None or row["chrom"] != chrom:  # no usable ClinGen location: Ensembl's
                hits[sym] = dict(row, chrom=chrom, start=g["start"], end=g["end"])
    records = []
    for sym, row in sorted(hits.items(), key=lambda kv: (kv[1]["start"] or 0)):
        rec = {"gene": row["gene"], "location": row["location"],
               "overlap": _gene_overlap_side(row, start, end, strand.get(sym)),
               "haploinsufficiency": {"score": row["haploinsufficiency"]["score"]},
               "triplosensitivity": {"score": row["triplosensitivity"]["score"]},
               "last_evaluated": row["last_evaluated"]}
        for which in ("haploinsufficiency", "triplosensitivity"):
            if row[which]["score"] in ("1", "2", "3"):
                rec[which].update(description=row[which]["description"], disease=row[which]["disease"],
                                  pmids=row[which]["pmids"])
        records.append(rec)
    return records, len(hits)


def _dosage_split(rows: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]],
                                                       List[Dict[str, Any]]]:
    """ClinGen's established HI (score 3) and TS (score 3) genes, plus every gene with any score at all."""
    hi, ts, scored = [], [], []
    for row in rows:
        h = (row.get("haploinsufficiency") or {}).get("score")
        t = (row.get("triplosensitivity") or {}).get("score")
        if str(h) == "3":
            hi.append({"gene": row.get("gene"), "score": h, "overlap": row.get("overlap"),
                       "description": (row.get("haploinsufficiency") or {}).get("description")})
        if str(t) == "3":
            ts.append({"gene": row.get("gene"), "score": t, "overlap": row.get("overlap"),
                       "description": (row.get("triplosensitivity") or {}).get("description")})
        if (h not in (None, "", "0", "Not yet evaluated")) or (t not in (None, "", "0", "Not yet evaluated")):
            scored.append({"gene": row.get("gene"), "haploinsufficiency_score": h, "triplosensitivity_score": t})
    return hi, ts, scored


def _gene_rows_hint(hi: List[Dict[str, Any]], ts: List[Dict[str, Any]], kind: str) -> List[str]:
    rows = []
    if kind == "loss":
        for g in hi:
            if g["overlap"] == "whole gene":
                rows.append(f"{g['gene']}: 2A row (complete overlap of an established HI gene)")
            elif g["overlap"] and "5' end" in g["overlap"]:
                rows.append(f"{g['gene']}: 2C rows (partial overlap with the 5' end of an established HI gene): "
                            "whether coding sequence is involved decides 2C-1 vs 2C-2")
            elif g["overlap"] and "3' end" in g["overlap"]:
                rows.append(f"{g['gene']}: 2D rows (partial overlap with the 3' end of an established HI gene)")
            elif g["overlap"] and "intragenic" in g["overlap"]:
                rows.append(f"{g['gene']}: 2E row (both breakpoints within the same gene: PVS1 specifications)")
            else:
                rows.append(f"{g['gene']}: 2C-2E rows (partial overlap of an established HI gene; which end is "
                            "involved decides)")
    elif kind == "gain":
        for g in ts:
            if g["overlap"] == "whole gene":
                rows.append(f"{g['gene']}: 2A row (the TS gene is fully contained within the gain)")
        for g in hi:
            if g["overlap"] == "whole gene":
                rows.append(f"{g['gene']}: 2H row (HI gene fully contained within a gain: continue evaluation)")
            elif g["overlap"] and "intragenic" in g["overlap"]:
                rows.append(f"{g['gene']}: 2I row (both breakpoints within the same HI gene)")
            else:
                rows.append(f"{g['gene']}: 2J/2K rows (one breakpoint within an established HI gene)")
    return rows


def _section3(kind: str, count: Optional[int], complete: bool = True) -> Dict[str, Any]:
    kinds = [kind] if kind in ("loss", "gain") else ["loss", "gain"]
    out: Dict[str, Any] = {
        "protein_coding_genes": count,
        "counted_as": "Ensembl protein_coding genes wholly or partially inside the interval; the standard counts "
                      "protein-coding RefSeq genes, so the two can differ by a few near a band edge",
        "clingen_bands": "the standard's section 3 bands (3A/3B/3C) for " + (
            "a copy-number " + kind + " are " + _BAND_TEXT[kind] if kind in ("loss", "gain") else
            "a copy-number loss are " + _BAND_TEXT["loss"] + "; for a gain " + _BAND_TEXT["gain"] +
            " — the report does not say which applies"),
        "source": _cite("RIGGS_2020"),
    }
    if count is not None and not complete:
        out["protein_coding_genes"] = None
        out["protein_coding_genes_at_least"] = count
    if count is not None:
        found = {}
        for k in kinds:
            for band, lo, hi, points in SECTION3_BANDS[k]:
                # with a partial gene list only the open top band is certain
                if count >= lo and (hi is None or count <= hi) and (complete or hi is None):
                    found[k] = f"{band} ({lo}{'+' if hi is None else '-' + str(hi)} genes; {points:.2f} points in the table)"
        out["band_for_this_count"] = found
    return out


def cnv_card(parsed: Dict[str, Any], assembly: str = "GRCh38", inheritance: Optional[str] = None,
             sex: Optional[str] = None) -> Outcome:
    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []
    chrom, start, end = parsed.get("chrom"), parsed["start"], parsed["end"]
    if not chrom:
        raise UsageError("the ISCN string does not name a chromosome; give the interval as "
                         "chr15:23123715-28193120 instead")
    if end < start:
        raise UsageError(f"the interval ends before it starts ({start} > {end})")
    if start < 1:
        raise UsageError(f"{start} is not a 1-based coordinate")
    build = parsed.get("assembly") or assembly
    length_of_chrom = MT_LENGTH if chrom == "MT" else CHROM_LENGTHS.get(build, {}).get(chrom)
    if length_of_chrom is None:
        raise UsageError(f"{chrom!r} is not a chromosome zebra can place in {build} "
                         "(1-22, X, Y, MT); check the report's notation")
    if end > length_of_chrom:
        raise UsageError(f"{chrom}:{start}-{end} runs past the end of chromosome {chrom} in {build} "
                         f"({length_of_chrom:,} bp): check the build and the coordinates")
    if parsed.get("assembly") and parsed["assembly"] != assembly:
        warnings.append(f"the report states {parsed['assembly']} and --assembly says {assembly}: the "
                        f"{parsed['assembly']} coordinates in the report were used, because builds must never be "
                        "mixed")
        assembly = parsed["assembly"]
    if parsed.get("inheritance") and inheritance and parsed["inheritance"] != inheritance:
        raise UsageError(f"the report's tag says {parsed['inheritance']} and --inheritance says {inheritance}")
    inheritance = inheritance or parsed.get("inheritance")
    parsed["in_par"] = _in_par(chrom, start, end, assembly)
    length = end - start + 1
    kind = parsed.get("cnv_type") or "unknown"

    # the ClinGen tables first: two cached requests that carry the dosage answer, so a
    # slow gene walk under the time budget cannot cost them
    tables = _clingen_tables(assembly, warnings, sources) if chrom != "MT" else \
        {"genes": None, "regions": None, "gene_file_date": None, "region_file_date": None}

    # genes spanned (Ensembl): skipped beyond MAX_REGION_BP, where the ClinGen tables still answer
    genes: Optional[List[Dict[str, Any]]] = None
    genes_complete = False
    if length <= MAX_REGION_BP:
        try:
            got = genes_in_region(chrom, start, end, assembly)
        except (SourceError, ValueError) as err:
            warnings.append(f"Ensembl overlap unavailable ({getattr(err, 'message', err)}): the genes spanned are not "
                            "listed, so section 1 and the section-3 gene count are not filled; the ClinGen region and "
                            "gene curations below were still checked by position")
        else:
            sources.extend(got.sources)
            warnings.extend(got.warnings)
            genes = got.result
            genes_complete = bool(got.query.get("complete", True))
    else:
        warnings.append(f"the interval is {length / 1e6:.1f} Mb: the Ensembl gene list was not walked (zebra lists "
                        f"genes for intervals up to {MAX_REGION_BP / 1e6:.0f} Mb), so section 1 and the section-3 "
                        "gene count are not filled; the ClinGen region and gene curations below were still checked")
    coding = [g for g in genes or [] if g["biotype"] == "protein_coding"]
    # one symbol can carry two Ensembl gene ids (a readthrough locus): count it once
    named_coding = list(dict.fromkeys(g["symbol"] for g in coding if g["named"]))
    gene_status = "checked" if tables["genes"] is not None else ("not applicable" if chrom == "MT" else "unavailable")
    region_status = "checked" if tables["regions"] is not None else ("not applicable" if chrom == "MT" else
                                                                     "unavailable")
    rows: List[Dict[str, Any]] = []
    curated_genes = 0
    if tables["genes"] is not None:
        rows, curated_genes = _gene_records(tables["genes"], genes, chrom, start, end)
    checked = len(named_coding) if (genes_complete and tables["genes"] is not None) else None
    hi, ts, scored = _dosage_split(rows)
    gene_pos = {g["gene"].upper(): (g["chrom"], g["start"], g["end"]) for g in tables["genes"] or []
                if g["start"] is not None}
    for g in genes or []:
        if g.get("named") and g["symbol"].upper() not in gene_pos and g.get("start"):
            gene_pos[g["symbol"].upper()] = (chrom, g["start"], g["end"])

    curated, benign, other = [], [], 0
    no_location = 0
    if tables["regions"] is not None:
        no_location = sum(1 for r in tables["regions"] if r["start"] is None)
        # the cytoband still says which chromosome a location-less curation is on
        unplaced_here = [r for r in tables["regions"] if r["start"] is None and
                         re.match(rf"^{re.escape(chrom)}[pq]", r.get("cytoband") or "", re.I)]
        for r in clingen.regions_overlapping(tables["regions"], chrom, start, end):
            entry = _region_entry(r, kind, gene_pos, chrom, start, end, genes_complete, parsed.get("copy_number"))
            # read each region by the score that matches this CNV: HI for a loss, TS for a gain (both when unknown)
            scores = [entry[w]["score"] for w in _relevant_scores(kind)]
            if any(sc in ("1", "2", "3") for sc in scores):
                curated.append(entry)
            elif any(sc == "40" for sc in scores):
                benign.append(entry)
            else:
                other += 1
        if unplaced_here:
            warnings.append(f"{len(unplaced_here)} ClinGen region curation(s) on chromosome {chrom} have no {assembly} "
                            "location ('tbd') and could not be checked for overlap: "
                            + "; ".join(f"{r['isca_id']} {r['name']} ({r['cytoband']})" for r in unplaced_here[:5]))
    curated.sort(key=lambda e: (-max(_num(e["haploinsufficiency"]["score"]), _num(e["triplosensitivity"]["score"])),
                                -(e["region_fraction_covered"] or 0)))
    # population (score 40) regions: the fields a curator compares (section 2F/2G), nothing else
    benign = [{"isca_id": e["isca_id"], "name": e["name"], "location": e["location"], "relation": e["relation"],
               "region_fraction_covered": e["region_fraction_covered"],
               "haploinsufficiency_score": e["haploinsufficiency"]["score"],
               "triplosensitivity_score": e["triplosensitivity"]["score"], "acmg_section_2": e["acmg_section_2"]}
              for e in sorted(benign, key=lambda e: (e["relation"] not in ("cnv_within_it", "identical"),
                                                     -(e["region_fraction_covered"] or 0)))]
    benign_total = len(benign)
    benign = benign[:BENIGN_LIMIT]

    cn = parsed.get("copy_number")
    region_note = ("ClinGen's region curations (recurrent and population CNV regions) were checked against this "
                   "interval: each overlapping region is listed with how much of it the CNV covers. A curated region "
                   "is not the CNV's interpretation: whether its causative gene or critical region is inside is the "
                   "curator's question (see named_genes).")
    result: Dict[str, Any] = {
        "kind": "cnv",
        "input": parsed["input"],
        "iscn": parsed.get("iscn"),
        "assembly": assembly,
        "region": {"chrom": chrom, "start": start, "end": end, "length_bp": length,
                   "band_as_reported": parsed.get("band")},
        "cnv_type": kind,
        "copy_number": cn,
        "genes": ({"total": len(genes), "protein_coding": len(coding),
                   "protein_coding_symbols": named_coding[:GENE_SYMBOL_LIMIT],
                   "unnamed_protein_coding": len([g for g in coding if not g["named"]]),
                   "other_biotypes": len(genes) - len(coding),
                   "symbols_truncated": len(named_coding) > GENE_SYMBOL_LIMIT,
                   "complete": genes_complete} if genes is not None else
                  {"total": None, "protein_coding": None, "protein_coding_symbols": [],
                   "note": (f"not listed: the interval exceeds {MAX_REGION_BP / 1e6:.0f} Mb" if length > MAX_REGION_BP
                            else "not listed: Ensembl did not answer (see warnings)")}),
        "clingen_regions": {"status": region_status, "file_date": tables["region_file_date"],
                            "regions_in_file": len(tables["regions"]) if tables["regions"] is not None else None,
                            "regions_without_location": no_location,
                            "dosage_curated": curated, "benign_or_population": benign,
                            "benign_or_population_total": benign_total,
                            "other_overlapping_without_dosage_evidence": other, "note": region_note},
        "clingen_dosage": {"status": gene_status, "file_date": tables["gene_file_date"],
                           "genes_checked": checked, "genes_with_a_record": len(rows), "records": rows,
                           "note": "every spanned protein-coding gene (by symbol) and every ClinGen gene curation whose "
                                   "location overlaps the interval (by position) was checked: no cut-off. Scores: 3 "
                                   "sufficient, 2 some, 1 little evidence, 0 no evidence, 30 autosomal recessive gene, "
                                   "40 dosage sensitivity unlikely",
                           "record_url": clingen.DOSAGE_GENE_PAGE.format("<gene>")},
        "acmg_cnv_inputs": {
            "framework": ACMG_CNV_FRAMEWORK,
            "scope": kind if kind in ("loss", "gain") else
                     ("unknown: " + parsed["cnv_type_note"] if parsed.get("cnv_type_note")
                      else "unknown (the report does not say whether this is a loss or a gain)"),
            "section_1_variant_type": {
                "contains_protein_coding_genes": (bool(coding) if genes_complete else (True if coding else None)),
                "protein_coding_gene_count": len(coding) if genes_complete else None,
                "note": "section 1 asks only whether protein-coding or other important elements are contained",
            },
            "section_2_overlap_with_established_regions_or_genes": {
                "clingen_established_regions": [
                    {k: e[k] for k in ("isca_id", "name", "relation", "region_fraction_covered", "acmg_section_2")}
                    for e in curated if any(e[w]["score"] == "3" for w in _relevant_scores(kind))],
                "clingen_benign_regions": [
                    {k: e[k] for k in ("isca_id", "relation", "acmg_section_2")} for e in benign],
                "clingen_established_haploinsufficient_genes": hi,
                "clingen_established_triplosensitive_genes": ts,
                "clingen_genes_with_any_dosage_score": scored,
                "genes_checked_in_clingen": checked,
                "candidate_rows_for_genes": _gene_rows_hint(hi, ts, kind),
                "note": "the acmg_section_2 entries name the row of the standard whose wording the overlap matches; "
                        "they are inputs for a curator, who confirms them (and decides 2B/2C-2E) against the "
                        "critical region and the patient's phenotype",
            },
            "section_3_gene_number": _section3(kind, len(coding) if genes is not None else None,
                                               complete=genes_complete),
            "section_4_case_level_evidence": None,
            "section_5_inheritance_and_family_history": {"reported_inheritance": inheritance},
            "classification": None,
            "classification_note": "sections 4 and 5 need case-level and family data, and the section totals are a "
                                   "curator's judgement: zebra reports the inputs and does not score or classify "
                                   "this CNV.",
        },
        "caveats": [
            "an array/CNV-seq interval is a called boundary, not a breakpoint: the true endpoints lie between the "
            "last normal and the first abnormal probe",
            "a balanced rearrangement, a low-level mosaic and an intragenic event below the assay's resolution are "
            "all invisible here",
            "gene content is Ensembl's gene set for this build; a gene's clinical relevance is not implied by "
            "being spanned",
        ],
    }
    # what was not checked is null, never an empty list that reads as "none found"
    sec2 = result["acmg_cnv_inputs"]["section_2_overlap_with_established_regions_or_genes"]
    if region_status != "checked":
        result["clingen_regions"].update(dosage_curated=None, benign_or_population=None, benign_or_population_total=None,
                                         other_overlapping_without_dosage_evidence=None, regions_without_location=None)
        sec2.update(clingen_established_regions=None, clingen_benign_regions=None)
    if gene_status != "checked":
        result["clingen_dosage"].update(genes_with_a_record=None, records=None)
        sec2.update(clingen_established_haploinsufficient_genes=None, clingen_established_triplosensitive_genes=None,
                    clingen_genes_with_any_dosage_score=None, candidate_rows_for_genes=None)
    for key in ("copy_number_range", "mosaic_fraction", "zygosity", "classification_lab", "extent",
                "inheritance_note"):
        if parsed.get(key):
            result[key] = parsed[key]
    # what a one- or two-copy loss means for genes ClinGen scores 30 (autosomal recessive phenotype)
    if kind == "loss" and gene_status == "checked" and chrom not in ("Y", "MT"):
        ar = [r for r in rows if r["haploinsufficiency"]["score"] == "30"]
        if ar and cn == 0 and chrom != "X":
            result["recessive_genes_lost"] = {
                "copy_number": 0,
                "genes": [{"gene": r["gene"], "overlap": r["overlap"]} for r in ar],
                "reading": "copy number 0: both copies are lost, so genes ClinGen associates with an autosomal "
                           "recessive phenotype (score 30) are lost on both alleles where the whole gene is inside. "
                           "ClinGen's HI scores and the population (score 40) regions describe the loss of ONE copy and "
                           "do not speak to a homozygous loss.",
                "source": _cite("RIGGS_2020")}
        elif ar and cn in (1, None) and chrom != "X":
            result["recessive_genes_lost"] = {
                "copy_number": cn,
                "genes": [{"gene": r["gene"], "overlap": r["overlap"]} for r in ar],
                "reading": "a one-copy loss of a gene ClinGen associates with an autosomal recessive phenotype (score "
                           "30) indicates carrier status for that condition, not disease (Riggs 2020), unless a second "
                           "variant is found on the other allele",
                "source": _cite("RIGGS_2020")}
    if chrom == "X" and kind == "loss" and sex != "male" and (hi or any(
            e["haploinsufficiency"]["score"] == "3" for e in curated)):
        result["caveats"].append(
            ("" if sex == "female" else "if the proband is female: ") +
            "a loss of a dosage-sensitive X-linked gene in a female carries a significant reproductive risk (her sons) "
            "and females may manifest symptoms in many X-linked disorders (Riggs 2020)")
    if parsed.get("mosaic"):
        result["mosaic"] = True
        result["caveats"].append("mosaic: only a fraction of cells carry the change, so the dosage effect is partial "
                                 "and tissue-dependent; the ACMG/ClinGen points are written for constitutional, "
                                 "non-mosaic CNVs")
    if parsed.get("in_par"):
        result["caveats"].append("the interval lies in a pseudoautosomal region of X, where two copies are normal "
                                 "in both sexes")
    if chrom == "MT":
        result["caveats"].append("mtDNA: ClinGen dosage maps and the ACMG/ClinGen CNV standard cover nuclear "
                                 "constitutional CNVs; an mtDNA deletion is read by its heteroplasmy level and tissue")
    if sex:
        result["sex_used"] = sex
    if cn is None and not parsed.get("copy_number_range") and chrom != "MT":
        warnings.append("no copy number in the input: pass --copies, or write the ISCN x-notation, so a loss can be "
                        "told from a homozygous loss")
    if parsed.get("cnv_type_note"):
        warnings.append(parsed["cnv_type_note"])
    if gene_status == "unavailable":
        warnings.append("ClinGen gene-level dosage was NOT checked (the ClinGen file could not be read): the absence "
                        "of established HI/TS genes below means 'not checked', not 'none'")
    if region_status == "unavailable":
        warnings.append("ClinGen region curations were NOT checked (the ClinGen file could not be read): a recurrent "
                        "microdeletion/duplication region may be missed")
    text = _cnv_text(result, named_coding, curated, benign, other, hi, ts, scored, gene_status, region_status, kind)
    return Outcome(result, sources=sources, warnings=warnings, text=text,
                   query={"input": parsed["input"], "assembly": assembly, "inheritance": inheritance})


def _num(score: Optional[str]) -> float:
    try:
        value = float(score or "")
    except ValueError:
        return -1.0
    return value if value <= 3 else -1.0  # 30/40 are categories, not stronger evidence


def _cnv_text(result: Dict[str, Any], named_coding: List[str], curated: List[Dict[str, Any]],
              benign: List[Dict[str, Any]], other: int, hi: List[Dict[str, Any]], ts: List[Dict[str, Any]],
              scored: List[Dict[str, Any]], gene_status: str, region_status: str, kind: str) -> str:
    reg = result["region"]
    cn = result.get("copy_number")
    head = f"CNV {reg['chrom']}:{reg['start']}-{reg['end']} ({result['assembly']}, {reg['length_bp']:,} bp) {kind}"
    if cn is not None:
        head += f", copy number {cn}"
    if result.get("copy_number_range"):
        head += f", mosaic copy number {result['copy_number_range']}"
    elif result.get("mosaic"):
        head += ", mosaic" + (f" ({result['mosaic_fraction']})" if result.get("mosaic_fraction") else "")
    lines = [head]
    genes = result["genes"]
    if genes.get("total") is None:
        lines.append(f"genes spanned: not listed ({genes.get('note')})")
    else:
        lines.append(f"genes spanned: {genes['total']} total, {genes['protein_coding']} protein-coding"
                     + (": " + ", ".join(named_coding[:15]) + (" …" if len(named_coding) > 15 else "")
                        if named_coding else ""))
    regions = result["clingen_regions"]
    if region_status == "checked":
        lines.append(f"ClinGen curated regions ({regions['regions_in_file']} in the {result['assembly']} file"
                     + (f", dated {regions['file_date']}" if regions.get("file_date") else "") + "): "
                     + (f"{len(curated)} dosage-curated region(s) overlap" if curated else
                        "no curated (dosage-scored) region overlaps this interval"))
        for e in curated[:8]:
            pct = (e["region_fraction_covered"] or 0) * 100
            cover = ("the CNV contains the whole region" if e["relation"] in ("cnv_contains_it", "identical") else
                     f"the CNV lies inside the region ({pct:.1f}% of it)" if e["relation"] == "cnv_within_it" else
                     f"partial: the CNV covers {pct:.1f}% of the region")
            named = ", ".join(f"{g['gene']} {'inside' if g['in_cnv'] == 'whole gene' else g['in_cnv'] or '?'}"
                              for g in e.get("named_genes") or [])
            score_txt = (f"TS {_score_text(e['triplosensitivity'])}, HI {e['haploinsufficiency']['score'] or '-'}"
                         if kind == "gain" else
                         f"HI {_score_text(e['haploinsufficiency'])}, TS {e['triplosensitivity']['score'] or '-'}")
            lines.append(f"  {e['isca_id']} {e['name']} [{e['location']}]: {score_txt} — {cover}"
                         + (f"; named gene(s): {named}" if named else ""))
            if e.get("acmg_section_2"):
                lines.append(f"    ACMG section 2 input: {e['acmg_section_2']}")
        if len(curated) > 8:
            lines.append(f"  … {len(curated) - 8} more in clingen_regions.dosage_curated")
        if benign:
            inside = sum(1 for e in benign if e["relation"] in ("cnv_within_it", "identical"))
            lines.append(f"  population/benign regions (dosage sensitivity unlikely, score 40) overlapped: "
                         f"{regions['benign_or_population_total']} ({inside} contain the whole CNV"
                         + (": " + ", ".join(e["isca_id"] for e in benign
                                             if e["relation"] in ("cnv_within_it", "identical"))[:200] if inside else "")
                         + ")")
        if other:
            lines.append(f"  {other} other overlapping region(s) carry no dosage evidence for this CNV type")
    elif region_status == "unavailable":
        lines.append("ClinGen curated regions: NOT CHECKED (the ClinGen region file could not be read — see warnings)")
    else:
        lines.append("ClinGen curated regions: not applicable to mtDNA")
    if gene_status == "checked":
        n_checked = result['clingen_dosage']['genes_checked']
        lines.append(f"ClinGen gene dosage: {('all ' + str(n_checked) + ' named protein-coding genes spanned were checked, ' if n_checked else 'no protein-coding gene spanned, ') if n_checked is not None else ''}"
                     f"{result['clingen_dosage']['genes_with_a_record']} with a ClinGen curation; "
                     f"established HI: {', '.join(x['gene'] + ' (' + str(x['overlap']) + ')' for x in hi) or 'none'}; "
                     f"established TS: {', '.join(x['gene'] for x in ts) or 'none'}")
        if scored:
            lines.append("  any dosage score: " + ", ".join(
                f"{x['gene']} HI {x['haploinsufficiency_score'] or '-'}/TS {x['triplosensitivity_score'] or '-'}"
                for x in scored[:12]) + (" …" if len(scored) > 12 else "")
                + "  (3 sufficient, 2 some, 1 little evidence; 30 autosomal recessive gene; 40 dosage sensitivity "
                  "unlikely)")
    elif gene_status == "unavailable":
        lines.append("ClinGen gene dosage: NOT CHECKED (the ClinGen gene file could not be read — see warnings)")
    sec3 = result["acmg_cnv_inputs"]["section_3_gene_number"]
    if sec3.get("band_for_this_count"):
        count = sec3["protein_coding_genes"] if sec3.get("protein_coding_genes") is not None else \
            f"at least {sec3.get('protein_coding_genes_at_least')}"
        lines.append("ACMG section 3: " + "; ".join(f"{k} {v}" for k, v in sec3["band_for_this_count"].items())
                     + f" for {count} protein-coding genes (Riggs 2020; Ensembl count — the standard counts RefSeq "
                       "genes, which can differ by a few near a band edge)")
    rec = result.get("recessive_genes_lost")
    if rec:
        lines.append(f"recessive genes (ClinGen score 30) in the interval: "
                     f"{', '.join(g['gene'] + ' (' + str(g['overlap']) + ')' for g in rec['genes'][:12])}"
                     f"{' …' if len(rec['genes']) > 12 else ''} — {rec['reading']}")
    for caveat in result.get("caveats") or []:
        if caveat.startswith(("mosaic:", "mtDNA:", "the interval lies in a pseudoautosomal")) or "X-linked" in caveat:
            lines.append(f"caveat: {caveat}")
    lines.append("ACMG/ClinGen CNV inputs reported; no classification (sections 4-5 need case-level data)")
    return "\n".join(lines)


# -------------------------------------------------------------- exon CNV

def _lookup(symbol: str, assembly: str, expand: bool = False) -> Outcome:
    """Ensembl's record for `symbol`; a 400/404 (no such symbol) is the caller's input error, not a source failure."""
    try:
        return ensembl.lookup_symbol(symbol, assembly, expand=expand)
    except SourceError as err:
        if err.status in (400, 404):
            raise UsageError(f"Ensembl has no gene with the symbol {symbol!r} in {assembly}: check the spelling "
                             "(current HGNC symbol)") from None
        raise


def _transcript(gene: str, assembly: str, wanted: Optional[str]) -> Tuple[Dict[str, Any], Outcome, List[str]]:
    """The reference transcript whose exon numbering is used, with its exons."""
    notes: List[str] = []
    got = _lookup(gene, assembly, expand=True)
    data = got.result
    if not isinstance(data, dict) or not data.get("Transcript"):
        raise UsageError(f"Ensembl has no transcripts for {gene!r} in {assembly}; check the gene symbol")
    txs = [t for t in data["Transcript"] if isinstance(t, dict) and t.get("Exon")]
    chosen = None
    if wanted:
        base = wanted.split(".")[0].upper()
        chosen = next((t for t in txs if str(t.get("id", "")).split(".")[0].upper() == base), None)
        if chosen is None and wanted.upper().startswith("ENST"):
            notes.append(f"the report names {wanted}, which is not one of Ensembl's {gene} transcripts: the Ensembl "
                         "canonical transcript was used instead. Check the transcript and the gene the report names.")
        elif chosen is None:
            notes.append(f"the report names {wanted}, which is not an Ensembl transcript id: the exon numbering below "
                         "is the Ensembl canonical transcript's. Check that the two number exons the same way.")
    if chosen is None:
        chosen = next((t for t in txs if t.get("is_canonical")), None)
    if chosen is None:
        if not txs:
            raise UsageError(f"Ensembl lists no transcript with exons for {gene!r} in {assembly}")
        chosen = max(txs, key=lambda t: t.get("length") or 0)
        notes.append("no canonical transcript was marked; the longest one was used")
    return chosen, got, notes


def _exons_in_order(tx: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The transcript's exons in transcript order (5' to 3'), which is exon numbering order."""
    return sorted(tx["Exon"], key=lambda e: e["start"], reverse=(tx.get("strand") or 1) < 0)


def coding_exons(tx: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Each exon in transcript order with its coding part (clipped to the Translation) and c. positions (B-P1-1)."""
    tl = tx.get("Translation") or {}
    t_start, t_end = tl.get("start"), tl.get("end")
    strand = tx.get("strand") or 1
    start_codon = None if t_start is None else (t_start if strand > 0 else t_end)
    stop_codon = None if t_end is None else (t_end if strand > 0 else t_start)
    out: List[Dict[str, Any]] = []
    cum = 0
    for i, e in enumerate(_exons_in_order(tx), 1):
        length = e["end"] - e["start"] + 1
        rec: Dict[str, Any] = {"exon": i, "chrom": e.get("seq_region_name"), "start": e["start"], "end": e["end"],
                               "length_bp": length, "ensembl_id": e.get("id")}
        if t_start is None or t_end is None:
            rec["coding_bp"] = None
        else:
            lo, hi = max(e["start"], t_start), min(e["end"], t_end)
            rec["coding_bp"] = max(0, hi - lo + 1)
            rec["has_start_codon"] = e["start"] <= start_codon <= e["end"]
            rec["has_stop_codon"] = e["start"] <= stop_codon <= e["end"]
            if rec["coding_bp"]:
                rec["c_first"], rec["c_last"] = cum + 1, cum + rec["coding_bp"]
                cum += rec["coding_bp"]
        out.append(rec)
    return out


def _frame(picked: List[Dict[str, Any]], kind: str, basis: str, caveats: Optional[List[str]] = None,
           span: Optional[int] = None) -> Dict[str, Any]:
    """The frame consequence of removing (or duplicating) `picked` exons, counted on coding bases only."""
    extra = list(caveats or [])
    verb = "removed" if kind == "loss" else ("added" if kind == "gain" else "affected")
    base = {"basis": basis}
    if span is None:
        coding = [p.get("coding_bp") for p in picked]
        genomic = sum(p["length_bp"] for p in picked)
        if any(c is None for c in coding):
            return dict(base, bases=None, genomic_bases=genomic, modulo_3=None, consequence="not assessed",
                        reading="the transcript carries no translation (non-coding), so there is no reading frame",
                        caveats=extra)
        span = sum(coding)  # type: ignore[arg-type]
        base["genomic_bases"] = genomic
        utr = [p["exon"] for p in picked if p.get("coding_bp") == 0]
        partial = [p["exon"] for p in picked if p.get("coding_bp") and p["coding_bp"] < p["length_bp"]]
        if utr:
            extra.append(f"exon(s) {', '.join(map(str, utr))} are untranslated (UTR only): they add no coding bases")
        if partial:
            extra.append(f"exon(s) {', '.join(map(str, partial))} are partly untranslated: only their coding part "
                         "is counted")
        codons = [("start", p["exon"]) for p in picked if p.get("has_start_codon")] + \
                 [("stop", p["exon"]) for p in picked if p.get("has_stop_codon")]
        if codons:
            which = " and ".join(f"the {c} codon (exon {e})" for c, e in codons)
            return dict(base, bases=span, modulo_3=None, consequence="not assessed",
                        reading=f"the {'deletion' if kind == 'loss' else 'change'} includes {which}: whether a protein "
                                "is made at all, from which start, and how it ends is not a reading-frame question, "
                                "so no in-frame/out-of-frame verdict is given",
                        caveats=extra)
        if span == 0:
            return dict(base, bases=0, modulo_3=0, consequence="no coding sequence",
                        reading="only untranslated sequence is involved: the reading frame is not affected (the "
                                "untranslated region can still matter for expression)",
                        caveats=extra)
    rest = span % 3
    consequence = "in frame" if rest == 0 else "out of frame"
    tail = ("a multiple of 3: the reading frame downstream is preserved" if rest == 0 else
            f"not a multiple of 3 ({rest} over): the reading frame downstream shifts")
    return dict(base, bases=span, modulo_3=rest, consequence=consequence, reading=f"{span} coding bases {verb} is {tail}",
                caveats=extra + [
                    "frame arithmetic only: whether the protein that results is functional, and whether the transcript "
                    "escapes nonsense-mediated decay, is not predicted here",
                ])


def _restoration(cx: List[Dict[str, Any]], first: int, last: int, total: int) -> Dict[str, Any]:
    """Which flanking exon's additional removal would make an out-of-frame deletion a multiple of 3."""
    options = []
    for n in (first - 1, last + 1):
        if not 1 <= n <= len(cx):
            continue
        e = cx[n - 1]
        opt: Dict[str, Any] = {"exon": n, "length_bp": e["length_bp"], "coding_bp": e.get("coding_bp")}
        if e.get("has_start_codon") or e.get("has_stop_codon"):
            opt.update(restores_frame=False, total_bp=None,
                       note=f"carries the {'start' if e.get('has_start_codon') else 'stop'} codon: removing it is not "
                            "a frame repair")
        elif not e.get("coding_bp"):
            opt.update(restores_frame=False, total_bp=total, note="untranslated (UTR only): removing it changes no "
                                                                  "coding base")
        else:
            opt.update(total_bp=total + e["coding_bp"], restores_frame=(total + e["coding_bp"]) % 3 == 0)
        options.append(opt)
    return {
        "candidates": options,
        "note": "arithmetic only: these are the flanking exons whose additional removal (by an exon-skipping "
                "approach or a larger deletion) would make the lost coding span a multiple of 3. Whether an approved "
                "or investigational therapy exists for that exon, and whether this patient would be eligible, is a "
                "separate question (`zebra therapy`, `zebra trials`).",
    }


def _restoration_block(result: Dict[str, Any], cx: List[Dict[str, Any]], first: int, last: int,
                       kind: str) -> None:
    frame = result["frame"]
    if kind != "loss":
        result["frame_restoration"] = None
        result["frame_restoration_note"] = "computed for deletions only"
    elif frame.get("consequence") == "out of frame":
        result["frame_restoration"] = _restoration(cx, first, last, frame["bases"])
    elif frame.get("consequence") == "in frame":
        result["frame_restoration"] = None
        result["frame_restoration_note"] = ("already in frame: the deletion keeps the reading frame, so there is no "
                                            "frame to restore and no exon to skip for that purpose")
    elif frame.get("consequence") == "no coding sequence":
        result["frame_restoration"] = None
        result["frame_restoration_note"] = "not computed: no coding base is removed, so there is no frame to restore"
    else:
        result["frame_restoration"] = None
        result["frame_restoration_note"] = "not computed: no in-frame/out-of-frame verdict (see frame.reading)"


def exon_card(parsed: Dict[str, Any], assembly: str = "GRCh38", gene: Optional[str] = None) -> Outcome:
    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []
    kind = parsed.get("cnv_type") or "unknown"
    symbol = (parsed.get("gene") or gene or "").upper() or None
    if parsed.get("gene") and gene and parsed["gene"].upper() != gene.upper():
        raise UsageError(f"the input names {parsed['gene']} and --gene says {gene}")
    result: Dict[str, Any] = {"kind": "exon_cnv", "input": parsed["input"], "assembly": assembly,
                              "gene": symbol, "cnv_type": kind}
    if parsed.get("zygosity"):
        result["zygosity"] = parsed["zygosity"]

    if parsed.get("cds_start") is not None:
        _hgvs_exon(parsed, result, symbol, assembly, kind, warnings, sources)
        return Outcome(result, sources=sources, warnings=warnings, text=_exon_text(result),
                       query={"input": parsed["input"], "assembly": assembly, "gene": symbol})

    if not symbol:
        raise UsageError("which gene? write the symbol before the exons (\"DMD exon 45-50 deletion\"), or pass "
                         "--gene when the report names only a transcript")
    tx, got, notes = _transcript(symbol, assembly, parsed.get("transcript"))
    sources.extend(got.sources)
    warnings.extend(notes)
    cx = coding_exons(tx)
    first, last = parsed["first"], parsed["last"]
    if last > len(cx):
        raise UsageError(f"{symbol} transcript {tx['id']} has {len(cx)} exons; the report names exon {last}. "
                         "Check which transcript the report numbers against.")
    picked = cx[first - 1:last]
    result["transcript"] = _tx_block(tx, parsed.get("transcript"), cx)
    result["exons"] = {"first": first, "last": last, "count": len(picked),
                       "per_exon": [{"exon": p["exon"], "chrom": p["chrom"], "start": p["start"], "end": p["end"],
                                     "length_bp": p["length_bp"], "coding_bp": p.get("coding_bp"),
                                     "ensembl_id": p["ensembl_id"]} for p in picked]}
    result["coordinates"] = {"chrom": picked[0]["chrom"],
                             "start": min(p["start"] for p in picked), "end": max(p["end"] for p in picked),
                             "note": "the exon boundaries; the real breakpoints lie in the flanking introns"}
    lengths = [p.get("coding_bp") if p.get("coding_bp") is not None else p["length_bp"] for p in picked]
    result["frame"] = _frame(picked, kind, basis=f"coding bases of exons {first}-{last} on {tx.get('id')} "
                                                 f"({'+'.join(str(x) for x in lengths)})")
    _restoration_block(result, cx, first, last, kind)
    if result["frame"]["caveats"]:
        warnings.extend(c for c in result["frame"]["caveats"] if "untranslated" in c)
    warnings.append(f"exon numbers are Ensembl's ordinal numbering of {tx.get('id')} (5' to 3'); check that the "
                    "report numbers exons on the same transcript")
    return Outcome(result, sources=sources, warnings=warnings, text=_exon_text(result),
                   query={"input": parsed["input"], "assembly": assembly, "gene": symbol})


def _tx_block(tx: Dict[str, Any], as_reported: Optional[str], cx: List[Dict[str, Any]]) -> Dict[str, Any]:
    tl = tx.get("Translation") or {}
    return {"as_reported": as_reported, "resolved": tx.get("id"), "name": tx.get("display_name"),
            "biotype": tx.get("biotype"), "exon_total": len(cx), "strand": tx.get("strand"),
            "coding_length_bp": sum(p.get("coding_bp") or 0 for p in cx) if tl else None,
            "exon_numbering_source": "Ensembl canonical transcript" if tx.get("is_canonical")
            else "Ensembl transcript as named",
            "exon_numbering": "Ensembl ordinal, 5' to 3'; a laboratory may number exons on another transcript or "
                              "in a legacy scheme"}


def _hgvs_exon(parsed: Dict[str, Any], result: Dict[str, Any], symbol: Optional[str], assembly: str, kind: str,
               warnings: List[str], sources: List[Dict[str, Any]]) -> None:
    """The HGVS exon-boundary form: the c. span gives the coding bases; with a gene, the exons and restoration."""
    cds_start, cds_end = parsed["cds_start"], parsed["cds_end"]
    span = cds_end - cds_start + 1
    reported = parsed.get("transcript")
    result["transcript"] = {"as_reported": reported, "resolved": None,
                            "exon_numbering_source": "not needed for the frame: the c. positions give the span"}
    result["cds_span"] = {"from": cds_start, "to": cds_end, "length_bp": span}
    result["hgvs_reading"] = (f"coding bases c.{cds_start} to c.{cds_end} are {'lost' if kind == 'loss' else 'affected'}"
                              f" ({parsed.get('hgvs_notation')}; an intronic breakpoint written N+k lies after coding "
                              "base N, one written M-k before coding base M)")
    result["coordinates"] = None
    if "?" in parsed["input"]:
        warnings.append("the '?' in an exon-boundary HGVS description means the intronic breakpoints were not "
                        "sequenced: the span is the coding sequence lost, which is what the reading frame depends on")
    if not symbol:
        result["frame"] = _frame([], kind, basis=f"c.{cds_start}_{cds_end} spans {span} coding bases on {reported} as "
                                                 "the report writes it", span=span)
        result["frame_restoration"] = None
        result["frame_restoration_note"] = ("not computed: give the gene (NM_004006.3(DMD):c.… or --gene DMD) so the "
                                            "span can be placed on the transcript's exons")
        return
    tx, got, notes = _transcript(symbol, assembly, reported)
    sources.extend(got.sources)
    # for a RefSeq/LRG id the MANE warning below says it better; any other note is kept
    warnings.extend(n for n in notes if "not an Ensembl transcript id" not in n)
    cx = coding_exons(tx)
    result["transcript"] = _tx_block(tx, reported, cx)
    result["gene_location"] = {"seq_region_name": tx.get("seq_region_name"), "start": tx.get("start"),
                               "end": tx.get("end"), "strand": tx.get("strand"), "id": tx.get("id")}
    total_cds = result["transcript"].get("coding_length_bp")
    if not total_cds:
        result["frame"] = _frame([], kind, basis=f"c.{cds_start}_{cds_end} on {reported}", span=span)
        result["frame_restoration"] = None
        result["frame_restoration_note"] = f"not computed: {tx.get('id')} has no translation in Ensembl"
        return
    if cds_end > total_cds:
        raise UsageError(f"c.{cds_end} is beyond the coding sequence of {tx.get('id')} ({total_cds} coding bases): "
                         "check the transcript the report numbers against")
    if reported and not reported.upper().startswith("ENST"):
        warnings.append(f"the c. positions of {reported} were placed on Ensembl {tx.get('id')}: the two share CDS "
                        "numbering when the Ensembl transcript is the MANE Select match of the RefSeq one — check "
                        "that it is")
    first_ex = next((p for p in cx if p.get("c_first") is not None and p["c_first"] <= cds_start <= p["c_last"]), None)
    last_ex = next((p for p in cx if p.get("c_first") is not None and p["c_first"] <= cds_end <= p["c_last"]), None)
    whole = bool(first_ex and last_ex and first_ex["c_first"] == cds_start and last_ex["c_last"] == cds_end)
    if first_ex and last_ex:
        picked = cx[first_ex["exon"] - 1:last_ex["exon"]]
        result["exons"] = {"first": first_ex["exon"], "last": last_ex["exon"], "count": len(picked),
                           "whole_exons": whole,
                           "per_exon": [{"exon": p["exon"], "chrom": p["chrom"], "start": p["start"], "end": p["end"],
                                         "length_bp": p["length_bp"], "coding_bp": p.get("coding_bp"),
                                         "c_first": p.get("c_first"), "c_last": p.get("c_last")} for p in picked]}
        result["coordinates"] = {"chrom": picked[0]["chrom"], "start": min(p["start"] for p in picked),
                                 "end": max(p["end"] for p in picked),
                                 "note": "the exon boundaries; the real breakpoints lie in the flanking introns"}
        if whole:
            # the same arithmetic as the exon form, on the same exons: start/stop codons and UTRs are handled there
            result["frame"] = _frame(picked, kind, basis=f"c.{cds_start}_{cds_end} = exons "
                                                         f"{first_ex['exon']}-{last_ex['exon']} of {tx.get('id')}")
        elif cds_start <= 3 or cds_end >= total_cds - 2:
            result["frame"] = dict(_frame([], kind, basis=f"c.{cds_start}_{cds_end}", span=span),
                                   modulo_3=None, consequence="not assessed",
                                   reading="the span reaches the start or the stop codon: whether a protein is made, "
                                           "and how it ends, is not a reading-frame question; no in/out-of-frame "
                                           "verdict is given")
        else:
            result["frame"] = _frame([], kind, basis=f"c.{cds_start}_{cds_end} spans {span} coding bases on "
                                                     f"{tx.get('id')}", span=span)
        if whole:
            _restoration_block(result, cx, first_ex["exon"], last_ex["exon"], kind)
        else:
            result["frame_restoration"] = None
            result["frame_restoration_note"] = (f"not computed: c.{cds_start}_{cds_end} does not start and end at exon "
                                                f"boundaries of {tx.get('id')} (exon {first_ex['exon']} starts at "
                                                f"c.{first_ex['c_first']}, exon {last_ex['exon']} ends at "
                                                f"c.{last_ex['c_last']})")
            warnings.append(result["frame_restoration_note"])
    else:
        result["frame"] = _frame([], kind, basis=f"c.{cds_start}_{cds_end}", span=span)
        result["frame_restoration"] = None
        result["frame_restoration_note"] = "not computed: the span could not be placed on the transcript's exons"


def _exon_text(result: Dict[str, Any]) -> str:
    frame = result.get("frame") or {}
    exons = result.get("exons")
    span = "" if not exons else (f" exon {exons['first']}" if exons["first"] == exons["last"]
                                 else f" exons {exons['first']}-{exons['last']}")
    lines = [f"{result.get('gene') or result['input']} {result['cnv_type']}{span}"]
    tx = result.get("transcript") or {}
    if tx.get("resolved"):
        lines.append(f"transcript {tx['resolved']} ({tx.get('name')}), {tx.get('exon_total')} exons, "
                     f"{tx.get('exon_numbering_source')}; exon numbers are Ensembl ordinals")
        if tx.get("as_reported"):
            lines.append(f"  as reported: {tx['as_reported']}")
    elif tx.get("as_reported"):
        lines.append(f"transcript as reported: {tx['as_reported']}")
    if result.get("cds_span"):
        lines.append(f"coding span: c.{result['cds_span']['from']}-c.{result['cds_span']['to']} "
                     f"({result['cds_span']['length_bp']} bp)")
    coords = result.get("coordinates")
    if coords:
        lines.append(f"exon coordinates: {coords['chrom']}:{coords['start']}-{coords['end']} ({result['assembly']})")
    if frame.get("modulo_3") is None:
        lines.append(f"frame: {frame.get('consequence')} — {frame.get('reading')}")
    else:
        lines.append(f"frame: {frame.get('bases')} coding bp, mod 3 = {frame.get('modulo_3')} → {frame.get('consequence')}")
        lines.append(f"  {frame.get('reading')}")
    rest = result.get("frame_restoration")
    if rest:
        for c in rest["candidates"]:
            if c.get("note") and not c.get("restores_frame"):
                lines.append(f"  + exon {c['exon']}: {c['note']}")
            else:
                lines.append(f"  + exon {c['exon']} ({c['coding_bp']} coding bp) → {c['total_bp']} bp, "
                             f"{'restores the frame' if c['restores_frame'] else 'still out of frame'}")
        lines.append(f"  ({rest['note']})")
    elif result.get("frame_restoration_note"):
        lines.append(f"  restoration: {result['frame_restoration_note']}")
    for caveat in frame.get("caveats") or []:
        lines.append(f"  caveat: {caveat}")
    return "\n".join(lines)


# ------------------------------------------- copy number / repeat expansion

COPY_NUMBER_MECHANISM = {
    "SMN1": [
        "most SMN1 assays (MLPA, ddPCR, qPCR) read the copy number of exon 7, because SMN1 and SMN2 differ at only "
        "a few bases; the result is a count of SMN1 copies, not a sequence",
        "a copy-number count cannot tell two copies on one chromosome (2+0) from one on each (1+1), so a two-copy "
        "result (apparently not a carrier) does not exclude a 2+0 'silent carrier'",
        "SMN2 copy number is reported alongside as a modifier of the SMN1 finding; it is not itself the finding",
        "a copy-number assay does not see an intragenic SMN1 point variant: a 1-copy result with a matching "
        "phenotype is usually followed by sequencing of the remaining copy",
    ],
}
COPY_NUMBER_MECHANISM["SMN2"] = [
    "SMN2 copy number is read alongside SMN1: it modifies the phenotype once SMN1 is lost on both alleles, and is "
    "not a finding by itself",
]
COPY_NUMBER_GENERIC = [
    "a copy-number count is a dosage measurement: it says how many copies the assay found, not which sequence they "
    "carry",
    "the assay's own reference ranges and limits belong to the laboratory report",
]
REPEAT_MECHANISM = [
    "a repeat-expansion result is a length measurement of one repeat tract, not a sequence of it",
    "repeat-primed PCR and long-range PCR size the tract; sizing is approximate at large expansions, and a somatic "
    "mosaic range is reported as a range",
    "interruptions in the motif and the flanking haplotype can change what a given length means, and are not "
    "visible in a length alone",
]

# SMN2 copy number and SMA type with supportive care only: GeneReviews NBK1352
# Table 3, "Adapted from Calucho et al [2018]" (rows "1 96% 4% 0%", "2 79% 16% 5%",
# "3 15% 54% 31%", "≥4 1% 11% 88%").
SMN2_TABLE_QUOTES = ["1\t96%\t4%\t0%", "2\t79%\t16%\t5%", "3\t15%\t54%\t31%", "≥4 4\t1%\t11%\t88%",
                     "Adapted from Calucho et al [2018]", "Clinical phenotype with supportive care only"]
SMN2_PHENOTYPE_TABLE = {
    1: {"SMA I": "96%", "SMA II": "4%", "SMA III/IV": "0%"},
    2: {"SMA I": "79%", "SMA II": "16%", "SMA III/IV": "5%"},
    3: {"SMA I": "15%", "SMA II": "54%", "SMA III/IV": "31%"},
    4: {"SMA I": "1%", "SMA II": "11%", "SMA III/IV": "88%"},  # the table's row is ">=4"
}
SMA_QUOTES = {
    "diagnosis": "biallelic pathogenic (or likely pathogenic) variants in SMN1 identified by molecular genetic testing",
    "exon7_95": "Exon 7 of SMN1 is undetectable in more than 95% of individuals with SMA irrespective of the clinical "
                "subtype of SMA",
    "one_copy": "If one copy of SMN1 exon 7 is present in an individual suspected to have SMA, sequence analysis of "
                "SMN1 is performed.",
    "compound": "Detects the 2%-5% of individuals who are compound heterozygous for an intragenic pathogenic variant "
                "and an SMN1 deletion of at least exon 7",
    "two_plus_zero": "about 5%-8% of the population have two copies of SMN1 on a single chromosome and a deletion on "
                     "the other chromosome, known as a [2+0] SMN1 genotype",
    "residual_670": "If such an individual is found to have at least two SMN1 copies, the probability of being a "
                    "carrier is approximately 1/670",
    "african": "Black individuals of sub-Saharan African descent have a higher proportion of the [2+0] genotype and "
               "have a lower detection rate (70%) than other populations",
    "modifier": "Increases in SMN2 copy number often modify the phenotype.",
    "smn2_range": "The number of copies of SMN2 may range from zero to five.",
    "smn2_80": "The presence of two copies of SMN2 is approximately 80% predictive of the SMA I phenotype, whereas the "
               "presence of four or more copies of SMN2 is approximately 88% predictive of achieving the ability to "
               "ambulate with supportive care only (SMA III/IV)",
    "modifying_factors": "Modifying factors that are not fully understood are likely to contribute to the variability "
                         "in clinical severity",
    "five_copies": "Prior et al [2004] reported three asymptomatic, unrelated individuals homozygous for an SMN1 "
                   "deletion who had five copies of SMN2",
    "treat_234": "Targeted treatment is recommended for all individuals who have two, three, or four copies of SMN2, "
                 "regardless of whether symptoms are present [Glascock et al 2020];",
    "treat_1": "For individuals who have one copy of SMN2, targeted treatment is left to the discretion of the "
               "treating provider, taking into account the severity of symptoms, which may have been present "
               "prenatally or at birth;",
    "treat_5": "For individuals with five copies of SMN2, targeted treatment may be deferred until symptom onset, "
               "although careful monitoring for the development of symptoms by a neuromuscular expert is recommended.",
    "nine_exons": "SMN1 and SMN2 each comprise nine exons",
    "presymptomatic": "Treatment with an SMA-specific disease-modifying treatment is most efficacious when initiated "
                      "presymptomatically.",
    "no_family_history": "Carrier screening for persons not known to have a family history of SMA requires SMN1 "
                         "dosage analysis.",
    "two_plus_zero_false": "Individuals with a [2+0] SMN1 genotype will have a false negative carrier screening result",
    "parents_6pct": "Approximately 6% of parents of a child with SMA resulting from a homozygous SMN1 deletion have "
                    "normal results of SMN1 dosage testing",
    "nbs_context": "The decision of when to initiate targeted therapy after detection of an affected individual via "
                   "newborn screening relies on genotype and presence of symptoms [Glascock et al 2018]. After "
                   "confirmatory SMN1 genetic testing:",
    "c859": "has been identified as a disease modifier resulting in a milder disease",
    "de_novo": "about 2% of affected individuals have a de novo SMN1 pathogenic variant on one allele",
}


def _smn_exon_readable(exon: Optional[str]) -> bool:
    """The GeneReviews readings are for the SMN1 exon 7 count (what dosage assays measure)."""
    if not exon:
        return True
    labels = {x for x in re.split(r"\s*(?:-|,|and|&)\s*", exon.lower()) if x}
    return "7" in labels and labels <= {"7", "8"}


def smn_reference(smn1: Optional[int], smn2: Optional[int], exon: Optional[str] = None) -> Dict[str, Any]:
    """What an SMN1 copy number, with SMN2 alongside, means — from GeneReviews NBK1352 only (CP1-5)."""
    q = SMA_QUOTES
    block: Dict[str, Any] = {"label": RESEARCH_LABEL, "source": _cite("GR_SMA")}
    if not _smn_exon_readable(exon):
        block["smn1"] = {"copy_number": smn1, "reading": (
            f"this result counts SMN1 exon {exon}. GeneReviews reads the SMN1 exon 7 count (the dosage assays measure "
            "exon 7 first), so no diagnostic or carrier reading is given for this exon on its own: ask the laboratory "
            "for the exon 7 result.")}
        block["exon_note"] = f"the readings in GeneReviews are for exon 7; this result names exon {exon}"
        if smn2 is not None:
            block["smn2"] = {"copy_number": smn2, "reading": "SMN2 copy number is read once SMA is confirmed on the "
                                                             "SMN1 exon 7 result"}
        return block
    if smn1 is None:
        block["smn1"] = {"copy_number": None, "reading": "the result does not say how many SMN1 exon 7 copies were "
                                                         "found, so it cannot be read for diagnosis or carrier status"}
    elif smn1 == 0:
        block["smn1"] = {"copy_number": 0, "reading": (
            "no copy of SMN1 exon 7: loss of SMN1 on both alleles, which meets the molecular criterion GeneReviews "
            "gives for 5q spinal muscular atrophy (SMA): biallelic SMN1 pathogenic variants, with exon 7 of SMN1 "
            "undetectable in more than 95% of people with SMA. Whether the person has SMA is a clinical diagnosis "
            "made with this result; found before symptoms (newborn screening, a younger sibling) it is the genotype "
            "of SMA, and GeneReviews notes that disease-specific treatment is most efficacious when started before "
            "symptoms. SMN2 copy number modifies the phenotype (three people with five SMN2 copies have been "
            "reported without symptoms)."),
            "quotes": [q["diagnosis"], q["exon7_95"], q["presymptomatic"], q["five_copies"]]}
    elif smn1 == 1:
        block["smn1"] = {"copy_number": 1, "reading": (
            "one copy of SMN1 exon 7: one SMN1 allele lacks exon 7 — consistent with carrier status. If the person "
            "has signs of SMA, GeneReviews has the remaining copy sequenced: 2%-5% of affected people carry a "
            "deletion on one allele and an intragenic variant on the other, which a copy-number assay does not see. "
            "About 2% of affected people have a de novo SMN1 variant, so a parent can test negative."),
            "quotes": [q["one_copy"], q["compound"], q["de_novo"]]}
    else:
        block["smn1"] = {"copy_number": smn1, "reading": (
            f"{smn1} copies of SMN1 exon 7: in the general population most often at least one copy on each "
            "chromosome, but a dosage test cannot see a [2+0] arrangement (two copies on one chromosome, none on the "
            "other) or an intragenic variant, and a [2+0] carrier screens as a false negative. For a person NOT "
            "known to have a family history of SMA, GeneReviews gives a carrier probability of approximately 1/670 "
            "with at least two copies; about 5%-8% of the population have a [2+0] SMN1 genotype, and detection is "
            "lower (70%) in Black individuals of sub-Saharan African descent. In a family with SMA this figure does "
            "not apply: approximately 6% of parents of a child with SMA (homozygous SMN1 deletion) have normal SMN1 "
            "dosage results, so a parent's or relative's normal count does not exclude carrier status."),
            "quotes": [q["residual_670"], q["no_family_history"], q["two_plus_zero"], q["two_plus_zero_false"],
                       q["african"], q["parents_6pct"]]}
    if smn2 is not None:
        smn2_block: Dict[str, Any] = {"copy_number": smn2}
        if smn1 == 0:
            row_key = 4 if smn2 >= 4 else smn2
            row = SMN2_PHENOTYPE_TABLE.get(row_key)
            smn2_block["reading"] = (
                "SMN2 copy number modifies the phenotype once SMN1 is lost: more SMN2 copies are on average associated "
                "with a milder phenotype, but the ranges overlap and modifying factors that are not fully understood "
                "also contribute (one known modifier is the SMN2 variant c.859G>C, which a copy number does not "
                "show). It describes groups of people, not this person's course, and treatment changes the natural "
                "history.")
            smn2_block["quotes"] = [q["modifier"], q["smn2_80"], q["modifying_factors"], q["c859"]]
            if row:
                smn2_block["phenotype_distribution_supportive_care_only"] = dict(
                    row, smn2_copies=(">=4" if smn2 >= 4 else str(smn2)),
                    table="GeneReviews NBK1352 Table 3, adapted from Calucho et al 2018: the share of people with "
                          "this SMN2 copy number in each SMA type with supportive care only (no targeted treatment)",
                    not_a_prediction="a group distribution from untreated cohorts, never a prognosis for one person")
            else:
                smn2_block["phenotype_distribution_supportive_care_only"] = None
                if smn2 == 0:
                    smn2_block["table_note"] = "GeneReviews Table 3 has no row for 0 copies of SMN2"
            if smn2 >= 5:
                smn2_block["five_or_more"] = q["five_copies"]
            timing = q["treat_234"] if 2 <= smn2 <= 4 else q["treat_1"] if smn2 == 1 else \
                q["treat_5"] if smn2 >= 5 else None
            if timing:
                smn2_block["treatment_timing_as_genereviews_summarises_it"] = (
                    "for an infant detected by newborn screening, after confirmatory SMN1 testing (GeneReviews, "
                    "citing Glascock et al 2018 and 2020; the decision also rests on symptoms and the treating "
                    "team): " + timing)
                smn2_block["timing_context_quote"] = q["nbs_context"]
        else:
            smn2_block["reading"] = (
                "SMN2 copy number is read once SMA is confirmed (loss of SMN1 on both alleles); alongside this SMN1 "
                "result it carries no diagnostic or carrier meaning of its own.")
        block["smn2"] = smn2_block
    elif smn1 == 0:
        block["smn2"] = {"copy_number": None, "reading": "SMN2 copy number was not given: GeneReviews uses it to "
                                                         "modify the phenotype once SMN1 is lost (pass "
                                                         "--smn2-copies, or \"SMN1 0 copies, SMN2 3 copies\")"}
    return block


# Repeat categories per locus, from the retrieved GeneReviews chapter (and a second
# source where the two disagree). Every number and reading below was read in the
# sentence quoted next to it (`quote`, `risk_quotes`, `note_quotes`, `range_notes`);
# tools/cnv/check_loci_sources.py checks each quote against the saved page.
# `meaning_male` / `meaning_female` hold what the source says only for one sex.
REPEAT_LOCI: Dict[str, Dict[str, Any]] = {
    "FMR1": {
        "symbol": "FMR1", "motif": "CGG", "source": "GR_FMR1", "inheritance": "X-linked",
        "inheritance_quote": "FMR1 disorders are inherited in an X-linked manner.",
        "disease": "FMR1 disorders: fragile X syndrome (FXS), fragile X-associated tremor/ataxia syndrome (FXTAS), "
                   "fragile X-associated primary ovarian insufficiency (FXPOI)",
        "categories": [
            {"name": "normal", "min": 5, "max": 44, "quote": "Normal alleles. Approximately 5-44 repeats",
             "meaning": "a normal allele"},
            {"name": "intermediate (gray zone)", "min": 45, "max": 54,
             "quote": "Intermediate alleles (also termed \"gray zone\" or \"borderline\"). Approximately 45-54 repeats",
             "meaning": "does not cause FXS; about 14% of intermediate alleles are unstable and may expand into the "
                        "premutation range when transmitted by the mother; they are not known to expand to a full "
                        "mutation",
             "risk_quotes": ["Intermediate alleles do not cause FXS. However, about 14% of intermediate alleles are "
                             "unstable and may expand into the premutation range when transmitted by the mother",
                             "They are not known to expand to full mutations"]},
            {"name": "premutation", "min": 55, "max": 200, "quote": "Premutation alleles. Approximately 55-200 repeats",
             "meaning": "not associated with FXS, but an increased risk of FXTAS (late-onset tremor and ataxia) and, in "
                        "women, FXPOI",
             "meaning_male": "FXTAS is estimated at about 40% of males with a premutation who are older than 50; a "
                             "male passes his premutation to all of his daughters and to none of his sons, and through "
                             "him it does not become a full mutation",
             "meaning_female": "FXTAS penetrance over age 50 is lower in females (16.5%) than in males (45.5%); FXPOI "
                               "has been seen in 20% of women with a premutation (1% in the general population); a "
                               "woman with a premutation is at risk of having children with FXS — a risk that depends "
                               "on the repeat size and, for small premutations, on the AGG interruptions",
             "risk_quotes": ["The prevalence of FXTAS is estimated at approximately 40% overall for males with a "
                             "premutation who are older than age 50 years",
                             "The penetrance in individuals older than age 50 years is lower in females (16.5%) than in "
                             "males (45.5%)",
                             "The premutation is inherited by all of their daughters and none of their sons.",
                             "When premutations are transmitted by the father, small increases in trinucleotide repeat "
                             "number may occur but do not result in full mutations.",
                             "has been observed in 20% of women who carry a premutation allele compared to 1% in the "
                             "general population",
                             "women with alleles in this range are considered to be at risk of having children with "
                             "FXS, although this risk is heavily dependent on the number of AGG interspersions for "
                             "small premutation alleles"]},
            {"name": "full mutation", "min": 201, "max": None,
             "quote": "Full-mutation alleles. More than 200 CGG repeats",
             "meaning": "a size in the full-mutation range. Fragile X syndrome is established with a full-mutation size "
                        "and abnormal methylation (most alleles above 200 repeats are methylated): the laboratory's "
                        "methylation result is part of the answer. Methylation mosaicism and unmethylated full "
                        "mutations occur, rarely with normal intellect",
             "meaning_male": "in GeneReviews Table 3, males with a completely methylated full mutation: 100% have "
                             "intellectual disability",
             "meaning_female": "in GeneReviews Table 3, females with a full mutation: about 50% have intellectual "
                               "disability and about 50% normal intellect",
             "risk_quotes": ["with abnormal gene methylation for most alleles with >200 repeats",
                             "Completely methylated\t100% have ID.", "~50% w/ID, ~50% normal intellect",
                             "Rarely, individuals with methylation mosaicism or completely unmethylated full mutations "
                             "and normal intellect have been reported."]},
        ],
        "span_readings": {
            ("premutation", "full mutation"): {
                "reading": "a range spanning premutation and full-mutation sizes: laboratories may report the somatic "
                           "variation of a full mutation as a range of several hundred repeats, and GeneReviews lists "
                           "repeat-size mosaicism (premutation and full-mutation cell lines) separately — in males "
                           "nearly 100% have intellectual disability, possibly higher functioning than with a full "
                           "mutation alone; in females it is highly variable. Read it as a full-mutation range or size "
                           "mosaicism, not as a premutation, together with the laboratory's methylation result",
                "quotes": ["clinical laboratories may report this somatic variation as a range of several hundred "
                           "repeats", "Nearly 100% have ID; may be higher functioning"]},
        },
        "range_notes": [
            {"low": 201, "high": 230,
             "note": "201-230 sits at the upper edge of the premutation range, which is sometimes given as about 230 "
                     "(both 200 and 230 are Southern-blot estimates)",
             "quote": "The upper limit of the premutation range is sometimes noted as approximately 230."},
        ],
        "notes": [
            "the boundaries are approximate ('Approximately 5-44', '55-200') and the distinction between categories "
            "is not absolute",
            "laboratories usually state a precision of ±2-3 repeats; GeneReviews: \"it may be prudent to consider "
            "reported test results with 55 repeats as potential premutations\"",
            "a male has one FMR1 allele (X-linked); a female's report gives two",
        ],
        "note_quotes": ["Thus, it may be prudent to consider reported test results with 55 repeats as potential "
                        "premutations."],
        "precision": 3,
        "x_linked": True,
    },
    "HTT": {
        "symbol": "HTT", "motif": "CAG", "source": "GR_HD", "inheritance": "autosomal dominant",
        "inheritance_quote": "HD is inherited in an autosomal dominant manner.",
        "disease": "Huntington disease (HD)",
        "categories": [
            {"name": "normal", "min": 0, "max": 26, "quote": "Normal alleles. 26 or fewer CAG repeats",
             "meaning": "a normal allele"},
            {"name": "intermediate", "min": 27, "max": 35,
             "quote": "Intermediate alleles. 27 to 35 CAG repeats. An individual with an allele in this range is not "
                      "typically at risk of developing manifestations of HD but, because of germline instability in "
                      "the CAG tract, may be at risk of having a child with an allele in the HD-causing range.",
             "meaning": "not typically at risk of developing HD, but because of germline instability (expansion is "
                        "more likely on paternal transmission) may have a child with an HD-causing allele"},
            {"name": "reduced penetrance (HD-causing)", "min": 36, "max": 39,
             "quote": "Reduced-penetrance HD-causing alleles. 36 to 39 CAG repeats. An individual with an allele in "
                      "this range is at risk for HD but may not develop manifestations of HD.",
             "meaning": "at risk for HD but may not develop manifestations; people without symptoms who carry 36-39 "
                        "repeats are common. In the general population the lifetime penetrance of 36-38 CAG alleles "
                        "is estimated at 0.2%-2%, and it is significantly higher in families with HD. A child has a "
                        "50% chance of inheriting the allele",
             "risk_quotes": ["Asymptomatic individuals with CAG repeats in this range are common",
                             "0.2%-2% for 36-38 CAG alleles across a typical life span",
                             "Penetrance estimates for individuals with 36 to 39 CAG repeats are significantly higher in "
                             "individuals from clinically ascertained families with HD than in the general population"]},
            {"name": "full penetrance (HD-causing)", "min": 40, "max": None,
             "quote": "Full-penetrance HD-causing alleles. 40 or more CAG repeats. Alleles of this size are associated "
                      "with development of HD with increased certainty assuming a normal life span.",
             "meaning": "associated with developing HD with increased certainty assuming a normal life span; a child "
                        "has a 50% chance of inheriting the allele"},
        ],
        "range_notes": [
            {"low": 56, "high": None, "note": "GeneReviews: people with juvenile-onset HD usually have more than 55 CAG "
                                              "repeats (an association across groups, not a forecast of onset)",
             "quote": "Individuals with juvenile-onset HD usually have an HTT allele with CAG repeats greater than 55."},
        ],
        "notes": [
            "current diagnostic methods underestimate HTT CAG length by two repeats in people with loss of the CAA "
            "interruption (LOI variants): a result near a boundary (34-35, 38-39) may belong to the next category",
            "a predictive test in a person without symptoms is done within a genetic-counselling protocol",
        ],
        "note_quotes": ["Current diagnostic methods underestimate HTT CAG length by two CAG repeats in individuals "
                        "with loss of the CAA interruption in the CAG repeat"],
        "apparent_homozygosity": "Detection of an apparently homozygous repeat does not rule out the presence of an "
                                 "expanded CAG repeat; thus, testing by TP-PCR or Southern blot analysis is required to "
                                 "detect a repeat expansion.",
        "precision": 2,
    },
    "DMPK": {
        "symbol": "DMPK", "motif": "CTG", "source": "GR_DM1", "inheritance": "autosomal dominant",
        "inheritance_quote": "DM1 is inherited in an autosomal dominant manner.",
        "disease": "myotonic dystrophy type 1 (DM1)",
        "categories": [
            {"name": "normal", "min": 5, "max": 34, "quote": "Normal alleles. 5-34 CTG repeats",
             "meaning": "a normal allele"},
            {"name": "mutable normal (premutation)", "min": 35, "max": 49,
             "quote": "Mutable normal (premutation) alleles. 35-49 CTG repeats. Individuals with CTG expansions in the "
                      "premutation range have not been reported to have symptoms, but their children are at "
                      "increased risk of inheriting a larger repeat size and thus having symptoms.",
             "meaning": "not reported to cause symptoms, but children are at increased risk of inheriting a larger "
                        "repeat and having symptoms"},
            {"name": "full penetrance", "min": 50, "max": None,
             "quote": "Full-penetrance alleles. >50 CTG repeats. Full-penetrance alleles are associated with disease "
                      "manifestations.",
             "meaning": "associated with disease manifestations; small abnormal expansions (50-99) are often "
                        "associated with a mild or asymptomatic phenotype. GeneReviews Table 2 relates size to "
                        "phenotype only roughly (mild 50-~150, classic ~100-~1,000, congenital >1,000, with "
                        "considerable overlap) and repeat size should not be used to predict severity",
             "risk_quotes": ["50-~150", "CTG repeat size should not be used to predict disease severity",
                             "Small but abnormal repeats (50-99) are often associated with a mild or asymptomatic "
                             "phenotype"]},
        ],
        "notes": [
            "GeneReviews writes the full-penetrance range as '>50' and its Table 2 starts the mild phenotype at 50; "
            "50 is placed in full penetrance here",
            "a second source draws the lines one repeat higher: the EMQN best-practice guidelines (Kamsteeg 2012, "
            "Table 1) list 5-35 normal, 36-50 'may be unstable, no DM', 51-150 and >150 — a result of 35 or 50 is "
            "classified differently by the two",
        ],
        "boundary_conflicts": {35: "EMQN 2012 lists 35 as normal (5-35); GeneReviews as premutation (35-49)",
                               50: "EMQN 2012 lists 50 as 'may be unstable, no DM' (36-50); GeneReviews Table 2 "
                                   "includes 50 in the mild phenotype (50-~150)"},
        "second_source": "EMQN_DM",
        "precision": 0,
    },
    "FXN": {
        "symbol": "FXN", "motif": "GAA", "source": "GR_FRDA", "inheritance": "autosomal recessive",
        "inheritance_quote": "FRDA is inherited in an autosomal recessive manner.",
        "disease": "Friedreich ataxia (FRDA)",
        "categories": [
            {"name": "normal", "min": 5, "max": 33, "quote": "Normal. 5-33 GAA repeats.",
             "meaning": "a normal allele"},
            {"name": "intermediate (mutable normal)", "min": 34, "max": 65,
             "quote": "Intermediate (also considered mutable normal). 34-65 GAA repeats.",
             "meaning": "a mutable normal allele"},
            {"name": "borderline", "min": 44, "max": 65,
             "quote": "Borderline. 44-65 GAA repeats. The shortest repeat length associated with disease (i.e., the "
                      "exact demarcation between normal and full-penetrance alleles) has not been clearly determined",
             "meaning": "overlaps the intermediate range: the shortest disease-associated length is not clearly "
                        "determined (an affected individual with an allele of 56 repeats has been reported)"},
            {"name": "pathogenic (full penetrance)", "min": 66, "max": None,
             "quote": "Pathogenic (full penetrance). 66 to approximately 1,300 GAA repeats.",
             "meaning": "a pathogenic (full-penetrance) GAA expansion; what it means depends on the other allele (see "
                        "the reading of both alleles)"},
        ],
        "notes": ["Friedreich ataxia is autosomal recessive: the meaning lies in the two alleles together"],
        "note_quotes": ["Approximately 96% of individuals with FRDA have biallelic FXN GAA repeat expansions in "
                        "intron 1",
                        "Penetrance is complete in individuals with either biallelic FXN pathogenic GAA repeat sizes",
                        "Individuals who are compound heterozygous for a borderline allele and a full-penetrance allele "
                        "(range 66-1,300 repeats) may develop LOFA or VLOFA.",
                        "If an apparently heterozygous FXN expanded GAA repeat allele in the full-penetrance or "
                        "intermediate range is identified, sequence analysis of FXN is performed next"],
        "precision": 0,
        "biallelic": True,
    },
    "C9ORF72": {
        "symbol": "C9orf72", "motif": "GGGGCC", "source": "GR_C9", "inheritance": "autosomal dominant",
        "inheritance_quote": "C9orf72-FTD/ALS is inherited in an autosomal dominant manner.",
        "disease": "C9orf72 frontotemporal dementia and/or amyotrophic lateral sclerosis (C9orf72-FTD/ALS)",
        "categories": [
            {"name": "normal", "min": 2, "max": 24, "quote": "Normal. Range from 2 to 24 G4C2 repeats",
             "meaning": "a normal allele"},
            {"name": "uncertain significance (intermediate)", "min": 25, "max": 60,
             "quote": "Uncertain significance. Range from 25 to 60 G4C2 repeats",
             "meaning": "of uncertain significance: rare in the general population and typically not segregating with "
                        "C9orf72-FTD/ALS in families — but the shortest repeat reported to cosegregate with the "
                        "disorder in a family was 47 (in blood; much longer expansions were found in brain)",
             "risk_quotes": ["Repeats in this range are rare in the general population and typically do not segregate "
                             "in families with C9orf72-FTD/ALS.",
                             "reported to cosegregate with the disorder in a family with C9orf72-FTD/ALS was 47 G4C2 "
                             "repeats"]},
            {"name": "pathogenic", "min": 61, "max": None,
             "quote": "Pathogenic. Range from 61 to >4000 of G4C2 repeats",
             "meaning": "pathogenic, with age-dependent reduced penetrance; the expansion cannot predict the disease "
                        "course in any given individual",
             "risk_quotes": ["Pathogenic expansions >60 to hundreds or thousands of G4C2 repeats show age-dependent "
                             "reduced penetrance"]},
        ],
        "notes": ["sizing of large expansions is approximate and somatic mosaicism between tissues is common"],
        "note_quotes": ["(Note: The presence of a C9orf72 G4C2 repeat expansion cannot predict the disease course in "
                        "any given individual.)"],
        "apparent_homozygosity": "detection of apparent homozygosity for a normal G4C2 repeat does not rule out the "
                                 "presence of an expanded G4C2 repeat; thus, testing by RP-PCR or Southern blotting is "
                                 "required.",
        "precision": 0,
    },
}


def _motif_equivalent(a: str, b: str) -> bool:
    """Same repeat unit read on either strand or from another starting base (CTG = CAG on the other strand)."""
    a, b = a.upper(), b.upper()
    if len(a) != len(b):
        return False
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    rc = "".join(comp.get(ch, ch) for ch in reversed(b))
    rotations = {b[i:] + b[:i] for i in range(len(b))} | {rc[i:] + rc[:i] for i in range(len(rc))}
    return a in rotations


def _categories_for(locus: Dict[str, Any], low: Optional[int], high: Optional[int]) -> List[str]:
    """The categories a count (or a count range [low, high]) falls in, in table order."""
    names: List[str] = []
    lo = low if low is not None else 0
    for cat in locus["categories"]:
        cmin, cmax = cat["min"], cat["max"]
        if high is None:
            hit = cmax is None or cmax >= lo
        else:
            hit = (cmax is None or lo <= cmax) and high >= cmin
        if hit:
            names.append(cat["name"])
    return names


def _sex_meaning(cat: Dict[str, Any], sex: Optional[str]) -> List[str]:
    """A category's reading, with the sentences that hold for one sex only given for that sex."""
    out = [cat["meaning"]]
    male, female = cat.get("meaning_male"), cat.get("meaning_female")
    if sex == "male" and male:
        out.append(f"in a male: {male}")
    elif sex == "female" and female:
        out.append(f"in a female: {female}")
    else:
        if male:
            out.append(f"in a male: {male}")
        if female:
            out.append(f"in a female: {female}")
    return out


def _fxn_joint(placed: List[Dict[str, Any]]) -> str:
    """Friedreich ataxia is recessive: the reading belongs to the two alleles together (GeneReviews NBK1281)."""
    def tops(entry: Dict[str, Any]) -> set:
        return set(entry["categories"])

    path = "pathogenic (full penetrance)"
    if len(placed) == 1:
        entry = placed[0]
        head = ("one size is given: whether it is one allele (the other not stated) or both alleles of the same size "
                "is not said — ask the laboratory. ")
        if path in tops(entry):
            return head + ("If only one allele carries this expansion: one pathogenic GAA expansion; if there are "
                           "signs of FRDA, GeneReviews has sequence analysis of FXN, then deletion/duplication "
                           "analysis, look for a second variant (about 4% of people with FRDA carry one expansion and "
                           "one sequence variant or deletion). If both alleles are this size: biallelic pathogenic "
                           "expansions, with which penetrance is complete.")
        return head + "No full-penetrance expansion at this size."
    a, b = tops(placed[0]), tops(placed[1])
    if path in a and path in b:
        return ("both alleles are pathogenic GAA expansions (biallelic): penetrance is complete with biallelic "
                "pathogenic GAA repeats, the genotype of about 96% of people with FRDA; the age of onset varies widely "
                "(from under 5 to over 50 years) and is not predicted by this result")
    if path in a or path in b:
        other = b if path in a else a
        if "borderline" in other:
            return ("one full-penetrance expansion with a borderline allele: GeneReviews notes that compound "
                    "heterozygotes for a borderline and a full-penetrance allele may develop late-onset or very "
                    "late-onset FRDA (LOFA/VLOFA), and reduced penetrance is possible with borderline alleles")
        return ("one pathogenic expansion and one allele below the borderline range (heterozygous for the expansion). "
                "If there are signs of FRDA, GeneReviews has sequence analysis of FXN, then deletion/duplication "
                "analysis, look for a second variant (about 4% of people with FRDA carry one expansion and one "
                "sequence variant or deletion)")
    if "borderline" in a or "borderline" in b:
        return ("no full-penetrance expansion; borderline allele(s) (44-65), where the shortest disease-associated "
                "length is not clearly determined")
    return "no pathogenic GAA expansion on either allele"


def repeat_reference(gene: str, motif: str, alleles: List[Dict[str, Any]],
                     sex: Optional[str] = None) -> Tuple[Optional[Dict[str, Any]], List[str]]:
    """The retrieved repeat categories for `gene`, and where each reported allele falls (CP1-5)."""
    warnings: List[str] = []
    locus = REPEAT_LOCI.get(gene.upper())
    if locus is None:
        return None, warnings
    if not _motif_equivalent(motif, locus["motif"]):
        warnings.append(f"{locus['symbol']}'s repeat is {locus['motif']} (GeneReviews); the result names {motif}, so "
                        "no category is given — check the report")
        return {"label": RESEARCH_LABEL, "gene": locus["symbol"], "motif_expected": locus["motif"],
                "motif_reported": motif, "categories_for_this_result": None, "source": _cite(locus["source"])}, warnings
    cats = {c["name"]: c for c in locus["categories"]}
    placed = []
    for allele in alleles:
        low, high = allele.get("low"), allele.get("high")
        names = _categories_for(locus, low, high)
        entry: Dict[str, Any] = {"allele": allele["allele"], "as_written": allele["as_written"],
                                 "categories": names}
        span = locus.get("span_readings", {}).get(tuple(n for n in names if n != "normal"))
        if not names:
            entry["note"] = "outside the ranges the source lists: check the report"
        elif span:
            entry["meanings"] = [span["reading"]]
        else:
            entry["meanings"] = [m for n in names if n != "normal" or len(names) == 1 for m in _sex_meaning(cats[n], sex)]
            if len(names) > 1 and not (locus.get("biallelic") and set(names) <= {"intermediate (mutable normal)",
                                                                                 "borderline"}):
                entry["note"] = "the result spans more than one category: " + " to ".join(names)
        precision = locus.get("precision") or 0
        if precision and low is not None and high is not None:
            # a boundary lies between (min - 1) and min of every category after the first
            edges = sorted({c["min"] for c in locus["categories"][1:]})
            if any(edge - precision <= v <= edge - 1 + precision for edge in edges for v in (low, high)):
                entry["near_a_boundary"] = (f"within {precision} repeats of a category boundary — see the notes on "
                                            "sizing precision")
        for edge, why in (locus.get("boundary_conflicts") or {}).items():
            if low is not None and high is not None and low <= edge <= high:
                entry["sources_disagree"] = why
        for rn in locus.get("range_notes") or []:
            top = high if high is not None else 10 ** 9
            if (low or 0) <= (rn["high"] if rn["high"] is not None else 10 ** 9) and top >= rn["low"]:
                entry.setdefault("range_notes", []).append(rn["note"])
        placed.append(entry)
    ref: Dict[str, Any] = {
        "label": RESEARCH_LABEL,
        "gene": locus["symbol"], "motif": locus["motif"], "disease": locus["disease"],
        "inheritance": locus["inheritance"],
        "categories": [{"name": c["name"], "range": (f"<={c['max']}" if not c["min"] else f"{c['min']}-{c['max']}")
                        if c["max"] is not None else f">={c['min']}", "meaning": c["meaning"]}
                       for c in locus["categories"]],
        "categories_for_this_result": placed,
        "notes": list(locus["notes"]),
        "source": _cite(locus["source"]),
    }
    if locus.get("x_linked"):
        ref["sex_used"] = sex
        if sex is None:
            ref["notes"].append("the proband's sex was not given: the readings for males and for females are both "
                                "shown (--sex, or the case profile)")
        if sex == "male" and len(alleles) > 1:
            warnings.append(f"{locus['symbol']} is on X and a male has one allele: two sizes for a male suggest size "
                            "mosaicism or a sample issue — check the report")
    if locus.get("biallelic"):
        ref["reading_of_both_alleles"] = _fxn_joint(placed)
    homo = locus.get("apparent_homozygosity")
    if homo and placed and all(p["categories"] == ["normal"] for p in placed) and (
            len(alleles) == 1 or len({(a.get("low"), a.get("high")) for a in alleles}) == 1):
        ref["apparent_homozygosity"] = ("one normal size only (or two identical): an apparently homozygous normal "
                                        "result does not rule out a large expansion that PCR does not amplify — "
                                        "GeneReviews: " + homo)
    if locus.get("second_source"):
        ref["second_source"] = _cite(locus["second_source"])
    return ref, warnings


def _smn_exon_mapping(symbol: str, exon: Optional[str], assembly: str, warnings: List[str],
                      sources: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """Legacy SMN exon labels resolved onto Ensembl's canonical transcript (B-P1-2)."""
    labels = [x for x in re.split(r"\s*(?:-|,|and|&)\s*", (exon or "").lower()) if x]
    if not labels:
        return None
    got = attempt(f"Ensembl lookup {symbol} (transcripts)", lambda: _lookup(symbol, assembly, expand=True), warnings)
    out: Dict[str, Any] = {"numbering": "legacy clinical numbering 1, 2a, 2b, 3-8 (nine exons)",
                           "as_reported": exon, "source": [_cite("GR_SMA"), _cite("SMN_EXONS")]}
    if got is None:
        out["ensembl"] = None
        return out
    sources.extend(got.sources)
    txs = [t for t in (got.result or {}).get("Transcript") or [] if isinstance(t, dict) and t.get("is_canonical")]
    if not txs or len(txs[0].get("Exon") or []) != 9:
        warnings.append(f"Ensembl's canonical {symbol} transcript does not have the nine exons of the clinical "
                        "numbering: the legacy exon label was not placed on it")
        out["ensembl"] = None
        return out
    tx = txs[0]
    exons = _exons_in_order(tx)
    mapped = []
    for label in labels:
        n = SMN_LEGACY_TO_ORDINAL.get(label)
        if n is None:
            continue
        e = exons[n - 1]
        mapped.append({"legacy_exon": label, "ensembl_ordinal_exon": n, "transcript": tx.get("id"),
                       "chrom": e.get("seq_region_name"), "start": e["start"], "end": e["end"],
                       "length_bp": e["end"] - e["start"] + 1})
    out["ensembl"] = mapped
    if mapped:
        warnings.append(f"{symbol} exons are numbered clinically 1, 2a, 2b, 3-8: clinical exon "
                        + ", ".join(f"{m['legacy_exon']} is Ensembl exon {m['ensembl_ordinal_exon']} of {m['transcript']}"
                                    for m in mapped))
    return out


def finding_card(parsed: Dict[str, Any], assembly: str = "GRCh38", method: Optional[str] = None,
                 related: Optional[List[str]] = None, smn2_copies: Optional[int] = None,
                 sex: Optional[str] = None) -> Outcome:
    """A copy-number or repeat-expansion result: structure, mechanism, and for SMN/repeat loci the cited meaning."""
    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []
    symbol = parsed["gene"]
    got = attempt(f"Ensembl lookup {symbol}", lambda: _lookup(symbol, assembly), warnings)
    location = None
    if got is not None and isinstance(got.result, dict):
        sources.extend(got.sources)
        location = {k: got.result.get(k) for k in ("seq_region_name", "start", "end", "strand", "id", "biotype")}
    else:
        warnings.append(f"{symbol} could not be confirmed against Ensembl: the symbol is recorded as given")

    result: Dict[str, Any] = {"kind": parsed["kind"], "input": parsed["input"], "gene": symbol,
                              "assembly": assembly, "gene_location": location, "method": method,
                              "classification": None,
                              "classification_note": "zebra does not classify this result. Where a retrieved "
                                                     "reference gives the meaning of the count, it is shown as a "
                                                     "research reference with its source; the laboratory's own "
                                                     "interpretation governs."}
    if parsed["kind"] == "copy_number":
        result["copy_number"] = parsed.get("copy_number")
        result["exon"] = parsed.get("exon")
        if parsed.get("zygosity"):
            result["zygosity"] = parsed["zygosity"]
        if parsed.get("copy_number_note"):
            warnings.append(parsed["copy_number_note"])
        result["related_results"] = []
        for other in parsed.get("related_parsed") or []:
            result["related_results"].append(dict(other))
        for extra in related or []:
            try:
                other = parse(extra)
            except UsageError as err:
                raise UsageError(f"--related {extra!r}: {err}") from None
            if other["kind"] != "copy_number":
                raise UsageError(f"--related takes further copy-number results, e.g. \"SMN2 copy number 2\"; "
                                 f"{extra!r} reads as {other['kind']}")
            result["related_results"].append({"gene": other["gene"], "copy_number": other["copy_number"],
                                              "exon": other.get("exon")})
        if smn2_copies is not None:
            if smn2_copies < 0:
                raise UsageError("--smn2-copies cannot be negative")
            if symbol != "SMN1":
                raise UsageError(f"--smn2-copies goes with an SMN1 result; this result is for {symbol}")
            result["related_results"].append({"gene": "SMN2", "copy_number": smn2_copies, "exon": None})
        seen: Dict[str, Any] = {}
        for other in result["related_results"]:
            if other["gene"] == symbol:
                raise UsageError(f"{symbol} is given twice (the result and a related result)")
            if other["gene"] in seen and seen[other["gene"]] != other["copy_number"]:
                raise UsageError(f"{other['gene']} is given two different copy numbers "
                                 f"({seen[other['gene']]} and {other['copy_number']})")
            seen[other["gene"]] = other["copy_number"]
        result["related_results"] = [dict(r) for r in {r["gene"]: r for r in result["related_results"]}.values()]
        result["mechanism"] = COPY_NUMBER_MECHANISM.get(symbol, []) + COPY_NUMBER_GENERIC
        if symbol not in COPY_NUMBER_MECHANISM:
            warnings.append(f"zebra carries no gene-specific mechanism notes for {symbol}: only the general ones "
                            "about copy-number assays are shown")
        if symbol in ("SMN1", "SMN2"):
            mapping = _smn_exon_mapping(symbol, parsed.get("exon"), assembly, warnings, sources)
            if mapping:
                result["exon_numbering"] = mapping
        if symbol == "SMN1":
            smn2 = next((r["copy_number"] for r in result["related_results"] if r["gene"] == "SMN2"), None)
            result["interpretation"] = smn_reference(parsed.get("copy_number"), smn2, parsed.get("exon"))
        elif symbol == "SMN2":
            result["interpretation"] = {"label": RESEARCH_LABEL, "source": _cite("GR_SMA"),
                                        "reading": "SMN2 copy number alone does not diagnose SMA or carrier status: "
                                                   "it is read alongside SMN1 (give the SMN1 result, with this as "
                                                   "--smn2-copies)",
                                        "quotes": [SMA_QUOTES["modifier"], SMA_QUOTES["smn2_range"]]}
        cn_text = parsed.get("copy_number")
        text = [f"{symbol} copy number {cn_text if cn_text is not None else 'not stated'}"
                + (f" (exon {parsed['exon']})" if parsed.get("exon") else "")]
        for other in result["related_results"]:
            text.append(f"  with {other['gene']} copy number {other['copy_number']}")
    else:
        result["motif"] = parsed["motif"]
        result["repeat_count"] = parsed["repeat_count"]
        result["alleles"] = parsed.get("alleles")
        result["mechanism"] = REPEAT_MECHANISM
        ref, ref_warnings = repeat_reference(symbol, parsed["motif"], parsed.get("alleles") or [], sex=sex)
        warnings.extend(ref_warnings)
        if ref is not None:
            result["gene"] = ref["gene"]
            result["thresholds"] = ref
            result["thresholds_note"] = ("categories from the retrieved GeneReviews chapter, cited with its revision "
                                         "date; a research reference, not a classification — the laboratory's own "
                                         "ranges and interpretation govern")
        else:
            result["thresholds"] = None
            result["thresholds_note"] = (f"zebra carries no retrieved threshold table for {symbol}: the normal / "
                                         "intermediate / full-expansion boundaries are gene-specific. Take them from "
                                         "the laboratory report or this gene's GeneReviews chapter (`zebra disease`, "
                                         "`zebra gene`).")
        text = [f"{result['gene']} {parsed['motif']} repeat: {parsed['repeat_count']}"]
    if location:
        text.append(f"  {result['gene']} at {location['seq_region_name']}:{location['start']}-{location['end']} "
                    f"({assembly}, {location['id']})")
    text.extend(_reference_text(result))
    text.append("  " + "\n  ".join(result["mechanism"]))
    text.append("  no classification by zebra" + ("; the reference above is a research reference with its source"
                                                  if result.get("interpretation") or result.get("thresholds") else ""))
    return Outcome(result, sources=sources, warnings=warnings, text="\n".join(text),
                   query={"input": parsed["input"], "assembly": assembly, "method": method})


def _reference_text(result: Dict[str, Any]) -> List[str]:
    lines: List[str] = []
    interp = result.get("interpretation")
    if interp:
        src = interp["source"]
        lines.append(f"  meaning (research reference: {src['citation']}, {src.get('revision', '')}, retrieved "
                     f"{src['retrieved_at'][:10]}):")
        if interp.get("smn1"):
            lines.append(f"    SMN1: {interp['smn1']['reading']}")
        if interp.get("smn2"):
            s2 = interp["smn2"]
            lines.append(f"    SMN2: {s2['reading']}")
            dist = s2.get("phenotype_distribution_supportive_care_only")
            if dist:
                lines.append(f"    SMN2 copy number {dist['smn2_copies']} in untreated cohorts (Calucho 2018, GeneReviews "
                             f"Table 3): SMA I {dist['SMA I']}, SMA II {dist['SMA II']}, SMA III/IV "
                             f"{dist['SMA III/IV']} — {dist['not_a_prediction']}")
            if s2.get("treatment_timing_as_genereviews_summarises_it"):
                lines.append(f"    treatment timing (GeneReviews): {s2['treatment_timing_as_genereviews_summarises_it']}")
            if s2.get("table_note"):
                lines.append(f"    note: {s2['table_note']}")
            if s2.get("five_or_more"):
                lines.append(f"    note (GeneReviews): {s2['five_or_more']}")
        if interp.get("reading"):
            lines.append(f"    {interp['reading']}")
        if interp.get("exon_note"):
            lines.append(f"    note: {interp['exon_note']}")
    ref = result.get("thresholds")
    if ref and ref.get("categories_for_this_result") is not None:
        src = ref["source"]
        lines.append(f"  categories (research reference: {src['citation']}, {src.get('revision', '')}, retrieved "
                     f"{src['retrieved_at'][:10]}):")
        lines.append("    " + "; ".join(f"{c['name']} {c['range']}" for c in ref["categories"]))
        for allele in ref["categories_for_this_result"]:
            lines.append(f"    {allele['as_written']} → {', '.join(allele['categories']) or 'outside the listed ranges'}"
                         + (f" ({allele['note']})" if allele.get("note") else "")
                         + (f" [sources disagree: {allele['sources_disagree']}]" if allele.get("sources_disagree") else "")
                         + (f" [{allele['near_a_boundary']}]" if allele.get("near_a_boundary") else ""))
            for meaning in allele.get("meanings") or []:
                if meaning != "a normal allele":
                    lines.append(f"      {meaning}")
            for note in allele.get("range_notes") or []:
                lines.append(f"      note: {note}")
        if ref.get("reading_of_both_alleles"):
            lines.append(f"    both alleles: {ref['reading_of_both_alleles']}")
        if ref.get("apparent_homozygosity"):
            lines.append(f"    note: {ref['apparent_homozygosity']}")
        for note in ref.get("notes") or []:
            lines.append(f"    note: {note}")
    return lines


# ------------------------------------------------------------------ entry

def card(text: str, assembly: str = "GRCh38", gene: Optional[str] = None, copies: Optional[int] = None,
         inheritance: Optional[str] = None, method: Optional[str] = None,
         related: Optional[List[str]] = None, smn2_copies: Optional[int] = None,
         sex: Optional[str] = None) -> Outcome:
    """Read one non-SNV result and return what can be established about it."""
    sex = _sex_norm(sex)
    parsed = parse(text, sex=sex, assembly=assembly)
    if copies is not None and copies < 0:
        raise UsageError("--copies cannot be negative")
    kind = parsed["kind"]
    # an option that has no meaning for this form is refused, never silently dropped
    if related and kind != "copy_number":
        raise UsageError(f"--related takes further copy-number results of an SMN1-type result; {text!r} reads as "
                         f"{kind}")
    if gene and kind != "exon_cnv" and gene.upper() != str(parsed.get("gene") or "").upper():
        raise UsageError(f"--gene names the gene of an exon-level result; {text!r} reads as {kind}"
                         + (f" for {parsed['gene']}" if parsed.get("gene") else ""))
    if copies is not None and parsed["kind"] in ("cnv", "copy_number"):
        if parsed.get("copy_number_range"):
            raise UsageError(f"the input gives a mosaic copy-number range ({parsed['copy_number_range']}) and --copies "
                             f"says {copies}: give one of the two")
        if parsed.get("copy_number") is not None and parsed["copy_number"] != copies:
            raise UsageError(f"the input says copy number {parsed['copy_number']} and --copies says {copies}")
        parsed["copy_number"] = copies
        if parsed["kind"] == "cnv":
            par = parsed.get("in_par", False)
            chrom = parsed.get("chrom")
            expected, note = _type_from_copies(chrom, copies, sex, par)
            if parsed.get("cnv_type") in (None, "unknown"):
                parsed["cnv_type"], parsed["cnv_type_note"] = expected, note
            elif chrom != "MT" and parsed["cnv_type"] not in _possible_types(chrom, copies, copies, sex, par):
                raise UsageError(f"the report says {parsed['cnv_type']} and copy number {copies} on {chrom} cannot be "
                                 f"a {parsed['cnv_type']}: one of the two is wrong, so zebra will not record either")
            zyg = parsed.get("zygosity")
            if zyg and chrom not in ("X", "Y", "MT") and ((zyg == "homozygous" and copies == 1) or
                                                         (zyg == "heterozygous" and copies == 0)):
                raise UsageError(f"the report says {zyg} and --copies says {copies}: they contradict each other")
        elif parsed.get("report_word"):
            word = parsed["report_word"]
            if (word == "loss" and copies > 1) or (word == "gain" and copies < 3):
                raise UsageError(f"the report says {'deletion' if word == 'loss' else 'duplication'} and --copies "
                                 f"says {copies}: a deletion leaves 0 or 1 copies, a duplication 3 or more")
            parsed["copy_number_note"] = None
    if copies is not None and parsed["kind"] not in ("cnv", "copy_number"):
        raise UsageError(f"--copies applies to a CNV or a copy-number result; {text!r} reads as "
                         f"{parsed['kind']}, where a copy number has no meaning")
    if smn2_copies is not None and parsed["kind"] != "copy_number":
        raise UsageError(f"--smn2-copies goes with an SMN1 copy-number result; {text!r} reads as {parsed['kind']}")
    if parsed["kind"] == "copy_number" and smn2_copies is not None:
        given = next((r["copy_number"] for r in parsed.get("related_parsed") or [] if r["gene"] == "SMN2"), None)
        if given is not None:
            if given != smn2_copies:
                raise UsageError(f"the input says SMN2 copy number {given} and --smn2-copies says {smn2_copies}")
            smn2_copies = None  # already in the input
    if parsed["kind"] == "cnv":
        return cnv_card(parsed, assembly=assembly, inheritance=inheritance, sex=sex)
    if parsed["kind"] == "exon_cnv":
        return exon_card(parsed, assembly=assembly, gene=gene)
    return finding_card(parsed, assembly=assembly, method=method, related=related, smn2_copies=smn2_copies,
                        sex=sex)


# --------------------------------------------------------- case recording

def case_fields(result: Dict[str, Any], method: Optional[str] = None,
                inheritance: Optional[str] = None) -> Dict[str, Any]:
    """The finding as a case.json variant record (see `zebra.case.VARIANT_KINDS`)."""
    kind = result["kind"]
    out: Dict[str, Any] = {"kind": kind, "assembly": result.get("assembly"), "method": method,
                           "description": result["input"], "source": "zebra cnv",
                           "inheritance": inheritance, "zygosity": result.get("zygosity"),
                           "classification_lab": result.get("classification_lab")}
    notes: List[str] = []
    if kind == "cnv":
        r = result["region"]
        symbols = result["genes"].get("protein_coding_symbols") or []
        out.update({"region": f"{r['chrom']}:{r['start']}-{r['end']}", "iscn": result.get("iscn"),
                    "cnv_type": result.get("cnv_type"), "copy_number": result.get("copy_number"),
                    "genes": symbols[:50]})
        total = result["genes"].get("protein_coding")
        if total is not None and total > len(out["genes"]):
            notes.append(f"{total} protein-coding genes spanned; {len(out['genes'])} named ones are stored in genes")
        if result.get("copy_number_range"):
            notes.append(f"mosaic copy number {result['copy_number_range']}")
        elif result.get("mosaic"):
            notes.append("mosaic" + (f" ({result['mosaic_fraction']})" if result.get("mosaic_fraction") else ""))
        established = [e["isca_id"] + " " + e["name"] for e in
                       (result.get("clingen_regions") or {}).get("dosage_curated") or []
                       if e["relation"] in ("cnv_contains_it", "identical") and
                       "3" in (e["haploinsufficiency"]["score"], e["triplosensitivity"]["score"])]
        if established:
            notes.append("contains ClinGen region(s) " + "; ".join(established[:3]))
    elif kind == "exon_cnv":
        exons = result.get("exons") or {}
        out.update({"gene": result.get("gene"), "cnv_type": result.get("cnv_type"),
                    "exons": (f"{exons['first']}-{exons['last']}" if exons and exons["first"] != exons["last"]
                              else (str(exons["first"]) if exons else None))})
        if (result.get("transcript") or {}).get("as_reported") and ":c." in result["input"]:
            out["hgvs_c"] = result["input"]
    elif kind == "copy_number":
        out.update({"gene": result.get("gene"), "copy_number": result.get("copy_number"),
                    "exons": result.get("exon")})
        related = result.get("related_results") or []
        if related:
            # SMN2 copy number is the modifier the mechanism note names: it must
            # not survive only in the free-text description
            notes.append("; ".join(f"{r['gene']} copy number {r['copy_number']}"
                                   + (f" (exon {r['exon']})" if r.get("exon") else "") for r in related))
            out["genes"] = [result.get("gene")] + [r["gene"] for r in related]
    else:
        out.update({"gene": result.get("gene"), "motif": result.get("motif"),
                    "repeat_count": result.get("repeat_count")})
    if notes:
        out["note"] = "; ".join(notes)
    return {k: v for k, v in out.items() if v is not None}


def record_identity(fields: Dict[str, Any]) -> Tuple[Any, ...]:
    """What makes two recorded findings the same one, so a retried --record does not add it twice (B-P2-8)."""
    kind = fields.get("kind")
    if kind == "cnv":
        # the note carries the mosaic state; inheritance and the lab's classification are part of the finding too
        return (kind, fields.get("region"), fields.get("cnv_type"), fields.get("copy_number"), fields.get("assembly"),
                fields.get("note"), fields.get("inheritance"), fields.get("classification_lab"))
    if kind == "exon_cnv":
        return (kind, (fields.get("gene") or "").upper(), fields.get("exons"), fields.get("hgvs_c"),
                fields.get("cnv_type"))
    if kind == "copy_number":
        return (kind, (fields.get("gene") or "").upper(), fields.get("copy_number"), fields.get("exons"),
                fields.get("note"))
    count = str(fields.get("repeat_count"))
    if "/" in count:  # two alleles in either order are the same result
        count = "/".join(sorted(count.split("/"), key=lambda x: (len(x), x)))
    return (kind, (fields.get("gene") or "").upper(), (fields.get("motif") or "").upper(), count)
