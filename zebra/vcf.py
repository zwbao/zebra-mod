"""Local VCF reanalysis: read, inspect and triage a singleton, duo or trio VCF.

The VCF never leaves this machine. What goes out is only what annotation needs:
gene symbols to Ensembl lookup (to find gene regions) and candidate variants as
chrom-pos-ref-alt (no sample names, genotypes or depths) to Ensembl VEP.

Reading: plain text or gzip/bgzip (bgzip is multi-member gzip, which the
stdlib reads). Multi-allelic records are split per ALT; alleles are trimmed to
a minimal representation (no left-alignment: that needs the reference);
`chr` prefixes are dropped and chrM becomes MT (Ensembl naming).

Triage (see `triage`): quality -> inheritance class -> restriction to a gene
set or a capped, prioritised subset BEFORE any web call -> VEP annotation and
gnomAD frequency -> comp-het resolution -> score (phenotype fit, consequence /
predictors, inheritance fit) -> ranked TSV.
"""

from __future__ import annotations

import bisect
import gzip
import io
import os
import re
import zlib
from collections import Counter, OrderedDict, defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Sequence, Set, Tuple

from zebra.core import Outcome, UsageError, attempt
from zebra.http import post_json, source_record
from zebra.sources import ensembl

# Chromosome lengths from Ensembl /info/assembly/homo_sapiens (GRCh38.p14 and
# GRCh37.p13), retrieved 2026-10-05. MT is identical in both and is not used
# for the build guess.
CHROM_LENGTHS: Dict[str, Dict[str, int]] = {
    "GRCh38": {
        "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259, "6": 170805979,
        "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
        "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345, "17": 83257441, "18": 80373285,
        "19": 58617616, "20": 64444167, "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
    },
    "GRCh37": {
        "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260, "6": 171115067,
        "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
        "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753, "17": 81195210, "18": 78077248,
        "19": 59128983, "20": 63025520, "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
    },
}
# Pseudoautosomal regions on X (GRC definitions).
PAR_X = {
    "GRCh38": [(10001, 2781479), (155701383, 156030895)],
    "GRCh37": [(60001, 2699520), (154931044, 155260560)],
}
_BUILD_TAGS = (
    ("GRCh38", re.compile(r"grch38|hg38|hs38|\bb38\b|GCA_000001405\.(1[5-9]|2\d)", re.I)),
    ("GRCh37", re.compile(r"grch37|hg19|hs37|\bb37\b|g1k_v37|GCA_000001405\.1\b", re.I)),
)

# HPO inheritance-mode terms (names checked against the local hp.json, release 2026-09-01).
INHERITANCE_TERMS = {
    "HP:0000006": "AD", "HP:0012275": "AD", "HP:0012274": "AD",  # incl. AD with maternal/paternal imprinting
    "HP:0000007": "AR", "HP:0032113": "SD",  # semidominant
    "HP:0001417": "XL", "HP:0001419": "XLR", "HP:0001423": "XLD",
    "HP:0001427": "MT", "HP:0001450": "YL",
}

DOMINANT_CLASSES = ("de_novo", "possible_de_novo", "inherited_het", "het")
CLASS_ORDER = ("de_novo", "hom_recessive", "hom", "x_hemizygous", "comphet", "comphet_unphased", "possible_de_novo",
               "mitochondrial", "het", "inherited_het", "y_hemizygous")
# what to annotate first when the budget (--max-annotate) cannot cover everything
PRIORITY = ("de_novo", "hom_recessive", "hom", "x_hemizygous", "possible_de_novo", "mitochondrial", "inherited_het",
            "het", "y_hemizygous")
INH_BASE = {"de_novo": 1.0, "hom_recessive": 0.9, "comphet": 0.9, "x_hemizygous": 0.9, "hom": 0.75,
            "possible_de_novo": 0.6, "comphet_unphased": 0.6, "mitochondrial": 0.5, "het": 0.4,
            "inherited_het": 0.2, "y_hemizygous": 0.3}
FIT_FACTOR = {"fits": 1.0, "unknown": 0.8, "against": 0.5}
SEVERITY = {"HIGH": 1.0, "MODERATE": 0.6, "LOW": 0.25, "MODIFIER": 0.1}
IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
WEIGHTS = {"phenotype": 0.40, "variant": 0.35, "inheritance": 0.25}
PREDICTOR_CAP = 0.9  # in-silico support never outranks a HIGH-impact (LoF) consequence
GRPMAX_POPS = ("afr", "amr", "eas", "nfe", "sas")  # gnomAD grpmax groups; bottlenecked asj/fin/ami/mid/remaining excluded
REGION_PAD = 50
TOP_RESULT = 30
VEP_CHUNK = 200


# ----------------------------------------------------------------- reading

def norm_chrom(chrom: str) -> str:
    c = chrom.strip()
    if c[:3].lower() == "chr":
        c = c[3:]
    u = c.upper()
    if u in ("M", "MT"):
        return "MT"
    if u in ("X", "Y"):
        return u
    return c


def normalize_allele(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """Trim shared trailing then leading bases, keeping one anchor base (no left-alignment)."""
    ref, alt = ref.upper(), alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def is_symbolic(alt: str) -> bool:
    return alt in ("*", ".", "") or alt.startswith("<") or "[" in alt or "]" in alt


def compression(path: str) -> str:
    with open(path, "rb") as fh:
        head = fh.read(18)
    if head[:2] != b"\x1f\x8b":
        return "plain"
    # BGZF: FEXTRA flag set and a 'BC' subfield
    if len(head) >= 14 and head[3] & 4 and head[12:14] == b"BC":
        return "bgzip"
    return "gzip"


def open_text(path: str) -> io.TextIOBase:
    p = os.path.expanduser(path)
    if not os.path.isfile(p):
        raise UsageError(f"no such file: {path}")
    if compression(p) != "plain":
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


_BASES_RE = re.compile(r"^[ACGTNacgtn]+$")
_STRUCT_RE = re.compile(r'([A-Za-z_][A-Za-z0-9_.]*)=("(?:[^"\\]|\\.)*"|[^,]*)')


def _parse_struct(body: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for m in _STRUCT_RE.finditer(body):
        val = m.group(2)
        out.setdefault(m.group(1), val[1:-1] if val.startswith('"') and val.endswith('"') else val)
    return out


@dataclass
class Header:
    fileformat: str = ""
    samples: List[str] = field(default_factory=list)
    has_format: bool = False
    contigs: "OrderedDict[str, Optional[int]]" = field(default_factory=OrderedDict)
    contig_assembly: List[str] = field(default_factory=list)
    reference: Optional[str] = None
    info: Dict[str, Dict[str, str]] = field(default_factory=dict)
    formats: Dict[str, Dict[str, str]] = field(default_factory=dict)
    filters: Dict[str, str] = field(default_factory=dict)
    csq_fields: Optional[List[str]] = None
    ann_fields: Optional[List[str]] = None
    lines: int = 0
    has_chrom_line: bool = False

    def gnomad_info_keys(self) -> List[str]:
        """INFO keys holding gnomAD population *frequencies* (never the cohort's own INFO/AF).

        The key name is matched whole: `gnomAD_AF`, `gnomAD_AF_afr`, `gnomAD_faf95`
        are frequencies; `gnomAD_AC_afr`, `gnomAD_AN_afr`, `gnomAD_nhomalt_afr`
        are counts and must never be read as one. When the header declares a
        Type, an Integer field is refused as well.
        """
        keys = []
        for k in self.info:
            low = k.lower()
            if any(p in low for p in _BOTTLENECKED):
                continue
            declared = (self.info.get(k) or {}).get("Type")
            if declared and declared != "Float":
                continue  # Integer means a count, whatever the name says
            if low in _POPMAX_AF_KEYS:
                keys.append(k)
                continue
            if "gnomad" not in low:
                continue
            rest = _GNOMAD_PREFIX_RE.sub("", low, count=1).strip("_.")
            rest = _SUBSET_RE.sub("", rest, count=1).strip("_.")  # gnomad_exomes_AF, gnomAD_genomes_AF_nfe
            if not (_AF_TYPE_RE.match(rest) or _GROUP_AF_RE.match(rest)):
                continue
            keys.append(k)
        return keys

    def annotations(self) -> Dict[str, Any]:
        return {"vep_csq": self.csq_fields is not None, "snpeff_ann": self.ann_fields is not None,
                "gnomad_af_info_keys": self.gnomad_info_keys(),
                "csq_af_fields": [f for f in (self.csq_fields or []) if _is_af_field(f)]}


_BOTTLENECKED = ("asj", "fin", "ami", "mid", "oth", "remaining")
_CSQ_AF_RE = re.compile(r"^gnomad[eg]?_(?:(?:afr|amr|eas|nfe|sas)_)?af$", re.I)
# an AF-type field name once the gnomAD prefix is gone: AF, AF_<group>, faf95, faf99
_AF_TYPE_RE = re.compile(r"^(af|faf95|faf99)(_[a-z0-9_]+)?$")
_GNOMAD_PREFIX_RE = re.compile(r"gnomad[a-z0-9]*")
_POPMAX_AF_KEYS = ("af_popmax", "af_grpmax", "popmax_af", "grpmax_af")
_SUBSET_RE = re.compile(r"^(?:exomes?|genomes?|joint)(?=_|$)")
_GROUP_AF_RE = re.compile(r"^[a-z]{2,8}_(af|faf95|faf99)$")  # gnomAD_AFR_AF
# gnomAD never reports a frequency above 1; a larger value means the field is a count
MAX_PLAUSIBLE_AF = 1.0
# ClinGen SVI's BA1 threshold: above this, no allele is a fully penetrant cause,
# so a ClinVar assertion no longer exempts it from the frequency filter.
BA1_AF = 0.05
# a shared FORMAT/PS beyond this distance is not read-backed phasing (linked-read
# and statistical phasing write chromosome-wide blocks), so it does not settle phase
PHASE_TRUSTED_SPAN = 500_000


def _is_af_field(name: str) -> bool:
    """CSQ population-frequency fields used for prefiltering: gnomAD overall and non-bottlenecked groups."""
    return bool(_CSQ_AF_RE.match(name))


def read_header(fh: io.TextIOBase) -> Tuple[Header, Optional[str]]:
    """Read meta lines and the #CHROM line; return the header and the first data line (or None)."""
    h = Header()
    for line in fh:
        line = line.rstrip("\r\n")
        h.lines += 1
        if line.startswith("##"):
            key, _, val = line[2:].partition("=")
            if key == "fileformat":
                h.fileformat = val
            elif key == "reference":
                h.reference = val
            elif val.startswith("<") and val.endswith(">"):
                body = _parse_struct(val[1:-1])
                ident = body.get("ID", "")
                if key == "contig":
                    length = body.get("length")
                    h.contigs[ident] = int(length) if length and length.isdigit() else None
                    if body.get("assembly"):
                        h.contig_assembly.append(body["assembly"])
                elif key == "INFO":
                    h.info[ident] = body
                    desc = body.get("Description", "")
                    if ident == "CSQ" and "Format:" in desc:
                        h.csq_fields = [f.strip() for f in desc.split("Format:", 1)[1].strip().strip('"').split("|")]
                    elif ident == "ANN":
                        m = re.search(r"'([^']+)'", desc)
                        if m:
                            h.ann_fields = [f.strip() for f in m.group(1).split("|")]
                elif key == "FORMAT":
                    h.formats[ident] = body
                elif key == "FILTER":
                    h.filters[ident] = body.get("Description", "")
            continue
        if line.startswith("#CHROM"):
            h.has_chrom_line = True
            cols = line.lstrip("#").split("\t")
            if len(cols) < 8:
                raise UsageError("not a VCF: the #CHROM header line has fewer than 8 columns")
            h.has_format = len(cols) > 8 and cols[8] == "FORMAT"
            h.samples = cols[9:] if h.has_format else []
            continue
        if line.startswith("#") or not line.strip():
            continue
        return h, line
    return h, None


@dataclass
class Record:
    chrom_raw: str
    chrom: str
    pos: int
    vid: str
    ref: str
    alts: List[str]
    qual: str
    filter: str
    info_raw: str
    fmt: List[str]
    sample_fields: List[str]

    def info(self) -> Dict[str, str]:
        out: Dict[str, str] = {}
        if self.info_raw in (".", ""):
            return out
        for item in self.info_raw.split(";"):
            k, eq, v = item.partition("=")
            out[k] = v if eq else "1"
        return out


def parse_line(line: str, n_samples: int) -> Optional[Record]:
    """One data line, or None with a reason when it cannot be read (see `parse_line_why`)."""
    return parse_line_why(line, n_samples)[0]


def parse_line_why(line: str, n_samples: int) -> Tuple[Optional[Record], Optional[str]]:
    cols = line.rstrip("\r\n").split("\t")
    if len(cols) < 8:
        return None, (f"{len(cols)} column(s), at least 8 required"
                      + ("; the columns look space-separated, not tab-separated" if len(cols) == 1 and " " in line
                         else ""))
    try:
        pos = int(cols[1])
    except ValueError:
        return None, f"POS {cols[1]!r} is not an integer"
    if pos < 1:
        return None, f"POS {pos} is not a 1-based position"
    if not cols[3] or not cols[4]:
        return None, "empty REF or ALT"
    if not _BASES_RE.match(cols[3]):
        return None, f"REF {cols[3]!r} is not a DNA sequence"
    fmt = cols[8].split(":") if len(cols) > 8 else []
    return Record(cols[0], norm_chrom(cols[0]), pos, cols[2], cols[3].upper(), cols[4].split(","), cols[5], cols[6],
                  cols[7], fmt, cols[9:9 + n_samples]), None


class RecordStream:
    """Iterate data lines once; close() releases the file even if iteration never started or stopped early.

    Lines that cannot be read are counted in `malformed` (with the first few
    line numbers and reasons in `malformed_examples`) instead of disappearing,
    and a truncated or corrupt gzip member becomes a UsageError rather than an
    EOFError traceback.
    """

    def __init__(self, fh: io.TextIOBase, first: Optional[str], n_samples: int, path: str = "",
                 first_lineno: int = 1):
        self._fh, self._first, self._n = fh, first, n_samples
        self._path = path
        self._lineno = first_lineno
        self.malformed = 0
        self.malformed_examples: List[str] = []
        self.data_lines = 0

    def _take(self, line: str, lineno: int) -> Optional[Record]:
        self.data_lines += 1
        rec, why = parse_line_why(line, self._n)
        if rec is None:
            self.malformed += 1
            if len(self.malformed_examples) < 3:
                self.malformed_examples.append(f"line {lineno}: {why}")
        return rec

    def __iter__(self) -> Iterator[Record]:
        try:
            if self._first is not None:
                rec = self._take(self._first, self._lineno)
                if rec:
                    yield rec
            while True:
                self._lineno += 1
                try:
                    line = next(self._fh)
                except StopIteration:
                    break
                except (EOFError, OSError, zlib.error) as err:
                    raise UsageError(
                        f"{self._path or 'the VCF'} is truncated or corrupt: the compressed stream ended after "
                        f"{self.data_lines} data line(s) ({type(err).__name__}: {err}). "
                        "Download or re-write the file and try again."
                    ) from None
                if not line or line[0] == "#":
                    continue
                rec = self._take(line, self._lineno)
                if rec:
                    yield rec
        finally:
            self._fh.close()

    def read_note(self) -> Optional[str]:
        if not self.malformed:
            return None
        return (f"{self.malformed} of {self.data_lines} data line(s) could not be read as variants: "
                + "; ".join(self.malformed_examples)
                + ("; …" if self.malformed > len(self.malformed_examples) else ""))

    def close(self) -> None:
        self._fh.close()


def iter_records(path: str) -> Tuple[Header, RecordStream]:
    fh = open_text(path)
    try:
        header, first = read_header(fh)
    except UnicodeDecodeError as err:
        fh.close()
        raise UsageError(f"{path} is not a text or gzip VCF ({err})") from None
    except (EOFError, OSError, zlib.error) as err:
        fh.close()
        raise UsageError(f"{path} is truncated or corrupt: the compressed stream ended inside the header "
                         f"({type(err).__name__}: {err})") from None
    if not header.has_chrom_line:
        fh.close()
        raise UsageError(f"{path} is not a VCF: no #CHROM header line")
    return header, RecordStream(fh, first, len(header.samples), path=path, first_lineno=header.lines)


# ---------------------------------------------------------------- genotypes

def _int(s: Optional[str]) -> Optional[int]:
    if s is None or s in ("", "."):
        return None
    try:
        return int(s)
    except ValueError:
        try:
            return int(float(s))
        except ValueError:
            return None


def _float(s: Any) -> Optional[float]:
    if s is None or s in ("", "."):
        return None
    if isinstance(s, (int, float)):
        return float(s)
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


@dataclass
class Call:
    raw: str
    gt: Optional[Tuple[Optional[int], ...]]
    phased: bool
    dp: Optional[int]
    gq: Optional[float]
    ad: Optional[List[Optional[int]]]
    ps: Optional[int] = None  # FORMAT/PS: the phase set this call belongs to

    def haplotype(self, k: int) -> Optional[int]:
        """Which haplotype carries ALT k, when the call is phased and carries it exactly once."""
        if not self.phased or self.gt is None:
            return None
        where = [i for i, a in enumerate(self.gt) if a == k]
        return where[0] if len(where) == 1 else None

    def other_alt(self, k: int) -> bool:
        """True when the call carries a different ALT as well (GT 1/2): the two are in trans."""
        if self.gt is None:
            return False
        return any(a not in (None, 0, k) for a in self.gt)

    def depth(self) -> Optional[int]:
        if self.dp is not None:
            return self.dp
        if self.ad and any(a is not None for a in self.ad):
            return sum(a for a in self.ad if a is not None)
        return None

    def reads(self, k: int) -> Tuple[Optional[int], Optional[int]]:
        if not self.ad or len(self.ad) <= k:
            return None, None
        return self.ad[0], self.ad[k]

    def ab(self, k: int) -> Optional[float]:
        ref, alt = self.reads(k)
        if ref is None or alt is None or ref + alt == 0:
            return None
        return alt / (ref + alt)


def parse_gt(raw: str) -> Tuple[Optional[Tuple[Optional[int], ...]], bool]:
    raw = raw.strip()
    if raw in ("", "."):
        return None, False
    phased = "|" in raw
    alleles: List[Optional[int]] = []
    for a in re.split(r"[/|]", raw):
        index = None if a in (".", "") else _int(a)
        alleles.append(index if index is None or index >= 0 else None)
    return tuple(alleles), phased


def parse_call(fmt: Sequence[str], field_: str) -> Call:
    vals = field_.split(":")
    d = dict(zip(fmt, vals))
    raw = d.get("GT", ".")
    gt, phased = parse_gt(raw)
    ad = None
    if d.get("AD") not in (None, "", "."):
        ad = [_int(x) for x in d["AD"].split(",")]
    return Call(raw, gt, phased, _int(d.get("DP")), _float(d.get("GQ")), ad, _int(d.get("PS")))


def zygosity(call: Call, k: int) -> str:
    """Zygosity for ALT allele k (1-based): hom_ref, het, hom_alt, hemi, other (another ALT only), missing."""
    if call.gt is None or all(a is None for a in call.gt):
        return "missing"
    called = [a for a in call.gt if a is not None]
    n_k = called.count(k)
    ploidy = len(call.gt)
    if n_k == 0:
        if any(a != 0 for a in called):
            return "other"
        if len(called) < ploidy:
            return "missing"
        return "hom_ref"
    if ploidy == 1:
        return "hemi"
    if n_k == ploidy:
        return "hom_alt"
    return "het"


CARRIER = ("het", "hom_alt", "hemi")


def _gt_text(call: Optional[Call], k: int) -> str:
    if call is None:
        return ""
    parts = [call.raw or "."]
    dp = call.depth()
    if dp is not None:
        parts.append(f"DP{dp}")
    if call.gq is not None:
        parts.append(f"GQ{call.gq:g}")
    ref, alt = call.reads(k)
    if ref is not None and alt is not None:
        parts.append(f"AD{ref},{alt}")
    return " ".join(parts)


# ------------------------------------------------------------------ inspect

def guess_build(header: Header) -> Dict[str, Any]:
    votes: Counter = Counter()
    evidence: List[str] = []
    for raw, length in header.contigs.items():
        c = norm_chrom(raw)
        if length is None or c == "MT":
            continue
        for build, table in CHROM_LENGTHS.items():
            if table.get(c) == length:
                votes[build] += 1
    for build, n in votes.items():
        evidence.append(f"{n} ##contig lengths match {build}")
    tags = []
    if header.reference:
        tags.append(("##reference", header.reference))
    for a in sorted(set(header.contig_assembly)):
        tags.append(("##contig assembly", a))
    tag_votes: Counter = Counter()
    for where, text in tags:
        for build, rx in _BUILD_TAGS:
            if rx.search(text):
                tag_votes[build] += 1
                evidence.append(f"{where} {text!r} names {build}")
                break
    guess: Optional[str] = None
    conflict = False
    if votes:
        guess = votes.most_common(1)[0][0]
        conflict = len(votes) > 1 or any(b != guess for b in tag_votes)
    elif tag_votes:
        guess = tag_votes.most_common(1)[0][0]
        conflict = len(tag_votes) > 1
    if conflict:
        guess = None
    return {"guess": guess, "evidence": evidence, "conflict": conflict}


def _ctype(chrom: str, pos: int, assembly: Optional[str]) -> str:
    if chrom == "MT":
        return "mt"
    if chrom == "Y":
        return "y"
    if chrom == "X":
        for s, e in PAR_X.get(assembly or "", []):
            if s <= pos <= e:
                return "x_par"
        return "x_nonpar"
    return "auto"


def inspect(path: str, max_seconds: float = 20.0, max_records: int = 3_000_000, max_samples: int = 50) -> Outcome:
    """Samples, variant count, contig naming, build guess, genotype/annotation content, per-sample sex hint."""
    import time

    t0 = time.monotonic()
    comp = compression(os.path.expanduser(path)) if os.path.isfile(os.path.expanduser(path)) else None
    header, records = iter_records(path)
    build = guess_build(header)
    stats_samples = header.samples[:max_samples]
    idx = {s: i for i, s in enumerate(header.samples)}
    per = {s: Counter() for s in stats_samples}
    n = multi = 0
    filters: Counter = Counter()
    chroms: Counter = Counter()
    naming: Counter = Counter()
    complete = True
    fmt_seen: Set[str] = set()
    for rec in records:
        n += 1
        if len(rec.alts) > 1:
            multi += 1
        filters[rec.filter] += 1
        chroms[rec.chrom] += 1
        naming["chr" if rec.chrom_raw.lower().startswith("chr") else "plain"] += 1
        if rec.fmt:
            fmt_seen.update(rec.fmt)
        if rec.fmt and "GT" in rec.fmt:
            ctype = _ctype(rec.chrom, rec.pos, build["guess"])
            for s in stats_samples:
                i = idx[s]
                if i >= len(rec.sample_fields):
                    continue
                gt, _ = parse_gt(dict(zip(rec.fmt, rec.sample_fields[i].split(":"))).get("GT", "."))
                if not gt or all(a is None for a in gt):
                    per[s]["missing"] += 1
                    continue
                if ctype == "x_nonpar":
                    per[s]["x_called"] += 1
                    if len(gt) == 1:
                        per[s]["x_haploid_any"] += 1
                called = [a for a in gt if a is not None]
                if not any(a and a > 0 for a in called):
                    per[s]["ref"] += 1
                    continue
                per[s]["alt"] += 1
                if len(gt) == 1:
                    kind = "haploid"
                elif len(set(called)) == 1 and len(called) == len(gt):
                    kind = "hom"
                else:
                    kind = "het"
                per[s][kind] += 1
                if ctype == "x_nonpar":
                    per[s]["x_alt"] += 1
                    per[s][f"x_{kind}"] += 1
        if n >= max_records or (n % 20000 == 0 and time.monotonic() - t0 > max_seconds):
            complete = False
            break
    records.close()
    samples_out = {}
    for s in stats_samples:
        c = per[s]
        x_alt = c["x_alt"]
        hint = "unknown"
        if c["x_haploid_any"] and c["x_haploid_any"] >= 0.5 * c["x_called"]:
            hint = "male-like (haploid X genotypes)"
        elif x_alt >= 20:
            het_frac = c["x_het"] / x_alt
            hint = "female-like" if het_frac > 0.2 else ("male-like" if het_frac < 0.05 else "unclear")
        samples_out[s] = {"alt_calls": c["alt"], "het": c["het"], "hom": c["hom"], "haploid": c["haploid"],
                          "missing": c["missing"], "x_nonpar_alt": x_alt, "x_het": c["x_het"],
                          "sex_hint": hint}
    contig_naming = "none seen"
    if naming:
        contig_naming = "mixed" if len(naming) > 1 else ("chr-prefixed (chr1)" if "chr" in naming else "plain (1)")
    has_gt = header.has_format and "GT" in (set(header.formats) | fmt_seen)
    result = {
        "path": str(Path(os.path.expanduser(path)).resolve()),
        "compression": comp,
        "fileformat": header.fileformat,
        "samples": header.samples,
        "variants": n,
        "scan_complete": complete,
        "multiallelic_records": multi,
        "filters": dict(filters.most_common(10)),
        "chromosomes": dict(sorted(chroms.items(), key=lambda kv: _chrom_sort_key(kv[0]))[:40]),
        "contig_naming": contig_naming,
        "declared_contigs": len(header.contigs),
        "build": build,
        "reference": header.reference,
        "genotypes": has_gt,
        "format_fields": sorted(fmt_seen)[:20],
        "annotations": header.annotations(),
        "per_sample": samples_out,
        "malformed_lines": records.malformed,
        "seconds": round(time.monotonic() - t0, 2),
    }
    notes = []
    read_note = records.read_note()
    if read_note:
        notes.append(read_note)
    if not complete:
        notes.append(f"scan stopped after {n} records ({result['seconds']} s); the count is a lower bound")
    if build["guess"] is None:
        notes.append("genome build unknown: pass --assembly GRCh38|GRCh37 to triage (never mix builds)")
    if len(header.samples) > max_samples:
        notes.append(f"per-sample statistics for the first {max_samples} of {len(header.samples)} samples")
    result["notes"] = notes
    lines = [
        f"{result['path']} ({comp}, {header.fileformat or 'VCF'})",
        f"samples ({len(header.samples)}): {', '.join(header.samples[:20])}{' …' if len(header.samples) > 20 else ''}",
        f"variants: {n}{'' if complete else '+ (scan capped)'}; multi-allelic records: {multi}",
        f"build: {build['guess'] or 'unknown'}" + (f" ({'; '.join(build['evidence'][:3])})" if build["evidence"] else ""),
        f"contig naming: {contig_naming}; declared contigs: {len(header.contigs)}",
        f"genotypes: {'yes' if has_gt else 'no'} (FORMAT {':'.join(result['format_fields'][:8])})",
        "annotations: " + ", ".join(
            [x for x in ("VEP CSQ" if result["annotations"]["vep_csq"] else "",
                         "snpEff ANN" if result["annotations"]["snpeff_ann"] else "",
                         ("gnomAD AF " + ",".join(result["annotations"]["gnomad_af_info_keys"][:4])) if result["annotations"]["gnomad_af_info_keys"] else "") if x]
            or ["none (triage will need --genes or --hpo-genes, or it caps annotation at --max-annotate)"]),
        "filters: " + ", ".join(f"{k}={v}" for k, v in result["filters"].items()),
    ]
    for s, v in samples_out.items():
        lines.append(f"  {s}: {v['alt_calls']} ALT calls ({v['het']} het, {v['hom']} hom, {v['haploid']} haploid); "
                     f"X non-PAR ALT {v['x_nonpar_alt']} ({v['x_het']} het) → sex {v['sex_hint']}")
    lines.extend(f"note: {x}" for x in notes)
    src = source_record("local VCF", os.path.basename(path), url=None, note="read on this machine; never uploaded")
    return Outcome(result, sources=[src], text="\n".join(lines), query={"vcf": path})


def _chrom_sort_key(c: str) -> Tuple[int, Any]:
    if c.isdigit():
        return (0, int(c))
    order = {"X": 23, "Y": 24, "MT": 25}
    return (1, order.get(c, 99), c) if c in order else (2, c)  # type: ignore[return-value]


# ---------------------------------------------------------- INFO annotations

def _vep_allele(ref: str, alts: Sequence[str], k: int) -> str:
    """The ALT as VEP writes it in CSQ: the first base is stripped when every allele shares it ('-' if empty)."""
    alleles = [ref] + list(alts)
    if all(alleles) and len({a[0] for a in alleles}) == 1:
        return alts[k - 1][1:] or "-"
    return alts[k - 1]


def info_annotation(header: Header, rec: Record, k: int) -> Optional[Dict[str, Any]]:
    """What the VCF's own INFO says about ALT k: worst impact, gene, consequence, max gnomAD AF."""
    if not (header.csq_fields or header.ann_fields or header.gnomad_info_keys()):
        return None
    info = rec.info()
    out: Dict[str, Any] = {"impact": None, "gene": None, "consequence": None, "af": None, "from": [],
                           "implausible_af_fields": []}
    afs: List[float] = []

    def take(impact: Optional[str], gene: Optional[str], csq: Optional[str], src: str) -> None:
        if impact and IMPACT_RANK.get(impact, -1) > IMPACT_RANK.get(out["impact"] or "", -1):
            out["impact"], out["gene"], out["consequence"] = impact, gene or None, csq or None
        if src not in out["from"]:
            out["from"].append(src)

    if header.csq_fields and "CSQ" in info:
        f = {n: i for i, n in enumerate(header.csq_fields)}
        allele = _vep_allele(rec.ref, rec.alts, k)
        for entry in info["CSQ"].split(","):
            parts = entry.split("|")
            get = lambda name: parts[f[name]] if name in f and f[name] < len(parts) else ""  # noqa: E731
            if "ALLELE_NUM" in f:
                if get("ALLELE_NUM") not in ("", str(k)):
                    continue
            elif get("Allele") not in (allele, rec.alts[k - 1]):
                continue
            take(get("IMPACT") or None, get("SYMBOL"), get("Consequence"), "CSQ")
            for name in f:
                if _is_af_field(name):
                    for v in get(name).split("&"):
                        x = _float(v)
                        if x is None:
                            continue
                        if x > MAX_PLAUSIBLE_AF or x < 0:
                            if name not in out["implausible_af_fields"]:
                                out["implausible_af_fields"].append(name)
                            continue
                        afs.append(x)
    if header.ann_fields and "ANN" in info:
        f = {n: i for i, n in enumerate(header.ann_fields)}
        for entry in info["ANN"].split(","):
            parts = [p.strip() for p in entry.split("|")]
            get = lambda name: parts[f[name]] if name in f and f[name] < len(parts) else ""  # noqa: E731
            if get("Allele") not in (rec.alts[k - 1], ""):
                continue
            take(get("Annotation_Impact") or None, get("Gene_Name"), get("Annotation"), "ANN")
    for key in header.gnomad_info_keys():
        if key not in info:
            continue
        vals = info[key].split(",")
        if len(vals) == len(rec.alts):
            vals = [vals[k - 1]]
        for v in vals:
            x = _float(v)
            if x is None:
                continue
            if x > MAX_PLAUSIBLE_AF or x < 0:
                # an allele count or allele number in an "AF" slot: refuse it
                # rather than drop the variant as "too common"
                if key not in out["implausible_af_fields"]:
                    out["implausible_af_fields"].append(key)
                continue
            afs.append(x)
        if "INFO AF" not in out["from"]:
            out["from"].append("INFO AF")
    out["af"] = max(afs) if afs else None
    return out


# ----------------------------------------------------------------- gene sets

def read_gene_file(path: str) -> List[str]:
    p = os.path.expanduser(path)
    if not os.path.isfile(p):
        raise UsageError(f"gene list not found: {path}")
    genes: List[str] = []
    with open(p, encoding="utf-8") as fh:
        for line in fh:
            line = line.split("#", 1)[0]
            for tok in re.split(r"[\s,;]+", line):
                tok = tok.strip()
                if tok and tok not in genes:
                    genes.append(tok)
    if not genes:
        raise UsageError(f"no gene symbols in {path}")
    return genes


def gene_regions(symbols: Sequence[str], assembly: str = "GRCh38") -> Outcome:
    """Gene spans for symbols via one batched POST /lookup/symbol/homo_sapiens (≤1000 per request)."""
    base = ensembl.host(assembly)
    uniq = list(OrderedDict.fromkeys(s for s in symbols if s))
    found: Dict[str, Dict[str, Any]] = {}
    sources: List[Dict[str, Any]] = []
    warnings: List[str] = []
    for i in range(0, len(uniq), 1000):
        chunk = uniq[i:i + 1000]
        resp = post_json(f"{base}/lookup/symbol/homo_sapiens", {"symbols": chunk}, source="Ensembl lookup",
                         cache_ttl=30 * 86400, timeout=120)
        data = resp.json()
        if not isinstance(data, dict):
            raise ValueError("Ensembl lookup returned no object")
        by_upper = {k.upper(): v for k, v in data.items() if isinstance(v, dict)}
        for sym in chunk:
            rec = by_upper.get(sym.upper())
            if not rec or rec.get("start") is None:
                continue
            found[sym] = {"symbol": rec.get("display_name") or sym, "chrom": norm_chrom(str(rec.get("seq_region_name"))),
                          "start": int(rec["start"]), "end": int(rec["end"]), "id": rec.get("id"),
                          "biotype": rec.get("biotype")}
        preview = ",".join(chunk[:5]) + (f",… ({len(chunk)} symbols)" if len(chunk) > 5 else "")
        sources.append(source_record("Ensembl lookup", preview, resp, note=f"{assembly}; POST /lookup/symbol (batch)"))
    missing = [s for s in uniq if s not in found]
    odd = [s for s, r in found.items() if r["chrom"] not in CHROM_LENGTHS["GRCh38"] and r["chrom"] != "MT"]
    if missing:
        warnings.append(f"{len(missing)} gene symbol(s) not found in Ensembl {assembly} (not searched): "
                        + ", ".join(missing[:20]) + (" …" if len(missing) > 20 else ""))
    if odd:
        warnings.append("genes placed off the primary chromosomes (not searched): " + ", ".join(odd[:10]))
        for s in odd:
            found.pop(s, None)
    return Outcome({"regions": found, "missing": missing}, sources=sources, warnings=warnings)


class RegionIndex:
    """Interval lookup: which padded gene regions overlap [start, end] on a chromosome."""

    def __init__(self, regions: Dict[str, Dict[str, Any]], pad: int = REGION_PAD):
        by_chrom: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
        for r in regions.values():
            by_chrom[r["chrom"]].append((r["start"] - pad, r["end"] + pad, r["symbol"]))
        self.index: Dict[str, Tuple[List[int], List[int], List[Tuple[int, int, str]]]] = {}
        for c, ivs in by_chrom.items():
            ivs.sort()
            starts = [s for s, _, _ in ivs]
            maxend, m = [], 0
            for _, e, _ in ivs:
                m = max(m, e)
                maxend.append(m)
            self.index[c] = (starts, maxend, ivs)

    def hits(self, chrom: str, start: int, end: int) -> List[str]:
        got = self.index.get(chrom)
        if not got:
            return []
        starts, maxend, ivs = got
        i = bisect.bisect_right(starts, end) - 1
        out = []
        while i >= 0 and maxend[i] >= start:
            s, e, sym = ivs[i]
            if s <= end and e >= start:
                out.append(sym)
            i -= 1
        return out


# ---------------------------------------------------------- HPO / phenotype

def _phenotype(case_dir: Optional[str], hpo_terms: Optional[Sequence[str]], warnings: List[str],
               need: bool) -> Optional[Dict[str, Any]]:
    from zebra import hpo_local

    present = [t for t in (hpo_terms or []) if t]
    excluded: List[str] = []
    origin = "--hpo" if present else None
    if not present and case_dir:
        from zebra import case as case_mod

        try:
            data = case_mod.load(case_dir)
        except case_mod.CaseError as err:
            if need:
                raise UsageError(str(err)) from None
            warnings.append(f"case not readable, phenotype fit not scored: {err}")
            return None
        present = [p["id"] for p in data.get("phenotypes", []) if p.get("status", "present") == "present"]
        excluded = [p["id"] for p in data.get("phenotypes", []) if p.get("status") == "excluded"]
        origin = f"case {case_dir}"
    if not present:
        if need:
            raise UsageError("--hpo-genes needs phenotypes: record HPO terms in the case (zebra case add-hpo) or pass --hpo")
        return None
    try:
        idx = hpo_local.load()
    except hpo_local.HpoDataMissing as err:
        if need:
            raise UsageError(f"--hpo-genes needs the local HPO files: {err}") from None
        warnings.append(f"phenotype fit not scored: {err}")
        return None
    try:
        res = hpo_local.rank(idx, present, excluded, top=1000)
    except ValueError as err:
        if need:
            raise UsageError(f"--hpo-genes: {err}") from None
        warnings.append(f"phenotype fit not scored: {err}")
        return None
    maxp = res["max_possible"] or 0.0
    genes: Dict[str, Dict[str, Any]] = OrderedDict()
    for g in res["genes"]:
        rel = g["score"] / maxp if maxp else 0.0
        genes[g["gene"]] = {"score": round(max(0.0, min(1.0, rel)), 3), "via": g["via"], "via_name": g["via_name"]}
    return {"origin": origin, "present": present, "excluded": excluded, "hpo_version": res["hpo_version"],
            "genes": genes, "notes": res["notes"],
            "method": "gene score = best disease Resnik score / query maximum (zebra hpo rank); "
                      "0 for genes outside the 200 best-matching diseases"}


def gene_moi(genes: Set[str]) -> Dict[str, List[str]]:
    """Inheritance modes HPO annotates to each gene's diseases (genes_to_phenotype.txt)."""
    from zebra import hpo_local

    path = hpo_local.data_dir() / "genes_to_phenotype.txt"
    if not genes or not path.exists():
        return {}
    out: Dict[str, Set[str]] = defaultdict(set)
    with open(path, encoding="utf-8") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        try:
            gi, hi = cols.index("gene_symbol"), cols.index("hpo_id")
        except ValueError:
            return {}
        for line in fh:
            parts = line.split("\t", hi + 1)
            if len(parts) <= hi:
                continue
            code = INHERITANCE_TERMS.get(parts[hi])
            if code and parts[gi] in genes:
                out[parts[gi]].add(code)
    return {g: sorted(v) for g, v in out.items()}


def moi_fit(cls: str, moi: Optional[Sequence[str]], ctype: str) -> str:
    if not moi:
        return "unknown"
    m = set(moi)
    if cls in DOMINANT_CLASSES:
        # a heterozygote fits dominant modes; on X, plain "X-linked" counts only when not also annotated recessive
        if m & {"AD", "SD", "XLD"} or (ctype.startswith("x") and "XL" in m and "XLR" not in m):
            return "fits"
        return "against" if m & {"AR", "XLR", "MT", "YL"} else "unknown"
    if cls in ("hom_recessive", "hom", "comphet", "comphet_unphased"):
        ok = {"AR", "SD"} | ({"XL", "XLR", "XLD"} if ctype.startswith("x") else set())
        if m & ok:
            return "fits"
        return "against" if m & {"AD", "MT", "YL"} else "unknown"
    if cls == "x_hemizygous":
        return "fits" if m & {"XL", "XLR", "XLD"} else "against"
    if cls == "mitochondrial":
        return "fits" if "MT" in m else "against"
    if cls == "y_hemizygous":
        return "fits" if "YL" in m else "against"
    return "unknown"


# ----------------------------------------------------------- VEP annotation

def _num(x: Any) -> Optional[float]:
    if isinstance(x, (int, float)):
        return float(x)
    if isinstance(x, str):
        vals = [_float(v) for v in re.split(r"[,&]", x)]
        vals = [v for v in vals if v is not None]
        return max(vals) if vals else None
    return None


def af_summary(groups: Dict[str, Any]) -> Dict[str, Any]:
    """gnomAD frequencies as VEP serves them: overall exome/genome, a grpmax-like maximum, and the bound to filter on.

    `filter_af` keeps its meaning (the most conservative frequency seen) and is
    what the recessive threshold uses. `bound_af` is what a *dominant*
    threshold may be applied to:

    - a filtering AF (faf95/faf99, the 95% CI lower bound) when one is served:
      that is what Whiffin 2017 and gnomAD publish it for;
    - otherwise the overall exome/genome estimate, because a single
      observation in a small group (AN ≈ 5,000 → AF 2e-4) is not evidence that
      a variant is too common to be pathogenic. The per-group point estimate is
      still reported, and the caller flags the divergence.
    """
    exome = _num(groups.get("gnomade"))
    genome = _num(groups.get("gnomadg"))
    grp = []
    faf: List[Tuple[float, str]] = []
    for k, v in groups.items():
        value = _num(v)
        if value is None:
            continue
        if "faf95" in k or "faf99" in k:
            faf.append((value, k))
            continue
        if k.startswith(("gnomade_", "gnomadg_")) and k.split("_", 1)[1] in GRPMAX_POPS:
            grp.append((value, k))
    grpmax, grpmax_group = (max(grp) if grp else (None, None))
    faf95, faf95_group = (max(faf) if faf else (None, None))
    gn = [x for x in (exome, genome, grpmax) if x is not None]
    out: Dict[str, Any] = {"exome": exome, "genome": genome, "grpmax": grpmax, "grpmax_group": grpmax_group,
                           "faf95": faf95, "faf95_group": faf95_group,
                           "filter_af": None, "bound_af": None, "bound_basis": None, "source": None}
    if gn:
        out["filter_af"], out["source"] = max(gn), "gnomAD"
    else:
        bottlenecked = [(value, k) for k, v in groups.items()
                        for value in [_num(v)]
                        if value is not None and k.startswith(("gnomade_", "gnomadg_"))
                        and any(p in k for p in _BOTTLENECKED)]
        kg = [_num(groups.get(k)) for k in ("af", "afr", "amr", "eas", "eur", "sas")]
        kg = [x for x in kg if x is not None]
        if bottlenecked:
            # gnomAD's grpmax excludes these groups; a frequency filter must not,
            # or a variant seen only in one bottlenecked population reads as absent
            value, group = max(bottlenecked)
            out["filter_af"], out["source"] = value, f"gnomAD {group} (a bottlenecked group, outside grpmax)"
            out["bottlenecked_only"] = group
        elif kg:
            out["filter_af"], out["source"] = max(kg), "1000 Genomes (no gnomAD frequency)"
    if faf95 is not None:
        out["bound_af"], out["bound_basis"] = faf95, f"filtering AF {faf95_group} (95% CI lower bound)"
    else:
        overall = [x for x in (exome, genome) if x is not None]
        if overall:
            out["bound_af"] = max(overall)
            out["bound_basis"] = ("overall gnomAD exome/genome point estimate (no filtering AF served; a per-group "
                                  "point estimate has no allele-number control)")
        elif out["filter_af"] is not None:
            out["bound_af"], out["bound_basis"] = out["filter_af"], f"point estimate ({out['source']})"
    return out


def _allele_keys(alt: Optional[str]) -> Set[str]:
    """The spellings VEP uses for one allele: "" and "-" both mean "the allele is absent"."""
    if alt is None:
        return set()
    keys = {alt, alt.strip()}
    if alt.strip() in ("-", ""):
        keys |= {"-", ""}
    return keys


# ClinVar significances that must never be filtered away on frequency alone
CLINVAR_PLP = ("pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic",
               "pathogenic_low_penetrance", "likely_pathogenic_low_penetrance")


CLINVAR_BENIGN = ("benign", "likely_benign", "benign/likely_benign")


def clinvar_pathogenic(ann: Optional[Dict[str, Any]]) -> List[str]:
    """The P/LP assertions VEP's colocated ClinVar record carries for this allele.

    A record that also carries benign or likely benign (ClinVar's "conflicting
    classifications") is not an assertion anyone can lean on, so it returns
    nothing: see `clinvar_conflicting`.
    """
    sigs = [str(s).strip().lower().replace(" ", "_") for s in (ann or {}).get("clinvar") or []]
    if any(s in CLINVAR_BENIGN for s in sigs):
        return []
    return [str(s) for s in (ann or {}).get("clinvar") or []
            if str(s).strip().lower().replace(" ", "_") in CLINVAR_PLP]


def clinvar_conflicting(ann: Optional[Dict[str, Any]]) -> bool:
    """True when the ClinVar record carries both a pathogenic and a benign classification."""
    sigs = [str(s).strip().lower().replace(" ", "_") for s in (ann or {}).get("clinvar") or []]
    return any(s in CLINVAR_PLP for s in sigs) and any(s in CLINVAR_BENIGN for s in sigs)


def annotate_vep(rec: Dict[str, Any], prefer: Optional[Set[str]] = None) -> Dict[str, Any]:
    """Fields that matter from one VEP record (MANE > canonical > protein-coding transcript)."""
    tcs = rec.get("transcript_consequences") or []
    pool = [t for t in tcs if prefer and t.get("gene_symbol") in prefer] or tcs
    tc = ensembl.pick_transcript({"transcript_consequences": pool}) if pool else None
    gene = tc.get("gene_symbol") if tc else None
    same = [t for t in tcs if gene and t.get("gene_symbol") == gene] or ([tc] if tc else [])
    if tc:
        impact = tc.get("impact") or "MODIFIER"
        terms = tc.get("consequence_terms") or []
    else:
        ic = (rec.get("intergenic_consequences") or [{}])[0]
        impact = ic.get("impact") or "MODIFIER"
        terms = ic.get("consequence_terms") or [rec.get("most_severe_consequence") or "intergenic_variant"]
    notes: List[str] = []
    worst = max(same, key=lambda t: IMPACT_RANK.get(t.get("impact") or "", -1)) if same else None
    if worst is not None and tc is not None and IMPACT_RANK.get(worst.get("impact") or "", -1) > IMPACT_RANK.get(impact, -1):
        notes.append(f"more severe on {worst.get('transcript_id')}: {','.join(worst.get('consequence_terms') or [])} "
                     f"({worst.get('impact')})")

    def first(getter):
        vals = [getter(t) for t in ([tc] if tc else []) + same]
        vals = [v for v in vals if v is not None]
        return vals[0] if vals else None

    revel = first(lambda t: _num(t.get("revel")))
    am = first(lambda t: t.get("alphamissense") if isinstance(t.get("alphamissense"), dict) else None) or {}
    cadd = first(lambda t: _num(t.get("cadd_phred")))
    sai_vals = [ensembl.spliceai_max(t) for t in same]
    sai_vals = [v for v in sai_vals if v is not None]
    spliceai = max(sai_vals) if sai_vals else None
    allele_string = rec.get("allele_string") or ""
    vep_alt = allele_string.split("/")[1] if "/" in allele_string else None
    freq = ensembl.frequencies(rec, vep_alt) if vep_alt else {}
    if freq and freq.get("allele") != vep_alt:
        freq = {}
    af = af_summary(freq.get("groups") or {})
    rsid = freq.get("id") if freq else None
    clin: List[str] = []
    # VEP writes an *empty* allele in clin_sig_allele for a deletion
    # (":pathogenic;:risk_factor") while allele_string gives it as "-": both
    # spellings mean the same allele, so match them to each other.
    alt_keys = _allele_keys(vep_alt)
    for cv in rec.get("colocated_variants") or []:
        csa = cv.get("clin_sig_allele")
        if csa and vep_alt:
            for item in str(csa).split(";"):
                a, sep, sig = item.partition(":")
                if not sep:
                    continue
                if a.strip() in alt_keys and sig and sig not in clin:
                    clin.append(sig)
        elif cv.get("clin_sig") and cv.get("allele_string", "").endswith("/" + (vep_alt or "?")):
            clin.extend(s for s in cv["clin_sig"] if s not in clin)
        if not rsid and str(cv.get("id", "")).startswith("rs") and cv.get("allele_string", "").split("/")[-1:] == [vep_alt]:
            rsid = cv.get("id")
    return {
        "gene": gene, "gene_id": tc.get("gene_id") if tc else None, "transcript": tc.get("transcript_id") if tc else None,
        "mane": tc.get("mane_select") if tc else None, "hgvsc": tc.get("hgvsc") if tc else None,
        "hgvsp": tc.get("hgvsp") if tc else None, "consequence": ",".join(terms), "impact": impact,
        "revel": revel, "alphamissense": _num(am.get("am_pathogenicity")), "am_class": am.get("am_class"),
        "cadd": cadd, "spliceai_max": spliceai, "af": af, "rsid": rsid, "clinvar": clin,
        "annotation": "Ensembl VEP", "notes": notes,
    }


def _vep_all(keys: List[Tuple[str, int, str, str]], assembly: str, warnings: List[str],
             sources: List[Dict[str, Any]]) -> Dict[Tuple[str, int, str, str], Dict[str, Any]]:
    chunks = [keys[i:i + VEP_CHUNK] for i in range(0, len(keys), VEP_CHUNK)]
    out: Dict[Tuple[str, int, str, str], Dict[str, Any]] = {}
    if not chunks:
        return out

    def run(i_chunk: Tuple[int, List[Tuple[str, int, str, str]]]):
        i, chunk = i_chunk
        return attempt(f"Ensembl VEP batch {i + 1}/{len(chunks)} ({len(chunk)} variants)",
                       lambda: ensembl.vep_batch(chunk, assembly), warnings)

    with ThreadPoolExecutor(max_workers=min(3, len(chunks))) as pool:
        results = list(pool.map(run, enumerate(chunks)))
    lines = {ensembl.vcf_line(*k): k for k in keys}
    for got in results:
        if got is None:
            continue
        sources.extend(got.sources)
        for rec in got.result or []:
            k = lines.get(str(rec.get("input", "")).strip())
            if k is not None:
                out[k] = rec
    return out


# ------------------------------------------------------------------ scoring

def predictor_support(a: Dict[str, Any], assembly: str) -> Tuple[float, List[str]]:
    """In-silico support as a level in [0, 0.9]; thresholds: REVEL (Pejaver 2022), SpliceAI (ClinGen SVI 2023),
    AlphaMissense class, CADD PHRED (Pejaver 2022)."""
    levels: List[Tuple[float, str]] = []
    revel = a.get("revel")
    if revel is not None:
        if revel >= 0.932:
            levels.append((0.85, f"REVEL {revel:.3f} ≥0.932 (strong)"))
        elif revel >= 0.773:
            levels.append((0.75, f"REVEL {revel:.3f} ≥0.773 (moderate)"))
        elif revel >= 0.644:
            levels.append((0.65, f"REVEL {revel:.3f} ≥0.644 (supporting)"))
    sai = a.get("spliceai_max")
    if sai is not None:
        if sai >= 0.5:
            levels.append((0.85, f"SpliceAI {sai:.2f} ≥0.5"))
        elif sai >= 0.2:
            levels.append((0.65, f"SpliceAI {sai:.2f} ≥0.2"))
    if a.get("am_class") == "likely_pathogenic":
        levels.append((0.7, f"AlphaMissense likely_pathogenic ({a.get('alphamissense')})"))
    cadd = a.get("cadd")
    if cadd is not None:
        if cadd >= 28.1:
            levels.append((0.65, f"CADD {cadd:g} ≥28.1"))
        elif cadd >= 25.3:
            levels.append((0.6, f"CADD {cadd:g} ≥25.3"))
    if not levels:
        return 0.0, []
    score = max(l for l, _ in levels)
    kinds = {txt.split()[0] for _, txt in levels}
    if len(kinds) >= 2:
        score += 0.05
    return round(min(PREDICTOR_CAP, score), 3), [t for _, t in levels]


def _phase_relation(a: Dict[str, Any], b: Dict[str, Any]) -> Tuple[Optional[str], str]:
    """'cis', 'trans' or None for two candidates in one gene, from the proband's own call.

    Two ALT alleles called at one site (GT 1/2) are necessarily on different
    chromosomes. Within one phase set (FORMAT/PS), a phased GT says which
    haplotype each ALT sits on; the same haplotype is cis, the other is trans.
    """
    if a["variant"] == b["variant"]:
        return None, ""  # the same allele called twice is not a pair
    if a["site"] == b["site"]:
        if a.get("trans_at_site") and b.get("trans_at_site"):
            return "trans", (f"two ALT alleles called at {a['site']} (GT {a['phase'].get('gt') or '1/2'}): "
                             "necessarily in trans, no parental testing needed")
        return None, ""
    pa, pb = a.get("phase") or {}, b.get("phase") or {}
    if pa.get("ps") is not None and pa.get("ps") == pb.get("ps") and pa.get("hap") is not None \
            and pb.get("hap") is not None:
        span = abs(a["pos"] - b["pos"])
        if span > PHASE_TRUSTED_SPAN:
            # a phase set spanning this much is not read-backed: population or
            # panel phasing written back as PS is not evidence of phase here
            return None, (f"a phase set (PS {pa['ps']}) covers both this allele and {b['variant']}, "
                          f"{span:,} bp away: too far apart for read-backed phasing, so the phase is treated as "
                          "unknown")
        if pa["hap"] == pb["hap"]:
            return "cis", (f"same haplotype in phase set PS {pa['ps']} as {b['variant']} ({span:,} bp away): in cis, "
                           "not a compound heterozygote")
        return "trans", (f"opposite haplotype in phase set PS {pa['ps']} as {b['variant']} ({span:,} bp away): "
                         "read-backed phase, in trans")
    return None, ""


def _phase_pairs(vs: List[Dict[str, Any]]):
    """Split the pairs of a gene's candidates into phase-proven trans, proven cis and unknown."""
    trans, cis, unknown = [], [], []
    for i in range(len(vs)):
        for j in range(i + 1, len(vs)):
            rel, why = _phase_relation(vs[i], vs[j])
            if rel == "trans":
                trans.append((vs[i], vs[j], why))
            elif rel == "cis":
                cis.append((vs[i], vs[j], why))
            else:
                unknown.append((vs[i], vs[j], why))
    return trans, cis, unknown


def _pair_up(members: List[Dict[str, Any]], partners_of, pairing: str, why: str) -> None:
    """Record a compound-heterozygous pairing without destroying a stronger class label.

    A de novo variant that happens to share its gene with a second qualifying
    het is still a de novo: overwriting its class would cost it the strongest
    label in triage, flip its gene-MOI fit and drop its score. The pairing is
    recorded alongside the class instead, which is what exempts the partner
    from the dominant allele-frequency cut-off and what the scoring reads.
    """
    for v in members:
        partners = [p for p in partners_of(v) if p != v["variant"]]
        if not partners:
            continue
        v["partners"] = sorted(set(partners))
        v["comphet"] = pairing
        if not v.get("af_exempt"):
            v["af_exempt"] = "compound-heterozygous partner (recessive threshold)"
        if v["class"] in ("het", "inherited_het"):
            v["class"] = "comphet" if pairing == "trans" else "comphet_unphased"
        if why and why not in v["flags"]:
            v["flags"].append(why)


def _qualifies(a: Optional[Dict[str, Any]]) -> bool:
    """Counts toward a compound-heterozygous pair: protein-altering, or predicted to affect splicing/function."""
    if not a or a.get("impact") is None:
        return True  # not annotated: cannot be excluded
    if a.get("impact") in ("HIGH", "MODERATE"):
        return True
    return bool((a.get("spliceai_max") or 0) >= 0.2 or (a.get("revel") or 0) >= 0.644
                or a.get("am_class") == "likely_pathogenic")


# ------------------------------------------------------------------- triage

def _classify(ctype: str, sex: Optional[str], pz: str, mother: Optional[Dict[str, Any]],
              father: Optional[Dict[str, Any]]) -> Tuple[str, Optional[str], List[str]]:
    flags: List[str] = []

    def carries(p: Optional[Dict[str, Any]]) -> bool:
        return p is not None and p["z"] in CARRIER

    has_m, has_f = mother is not None, father is not None
    if ctype == "mt":
        return "mitochondrial", ("maternal" if carries(mother) else None), flags
    if ctype == "y":
        if sex == "female":
            flags.append("Y-chromosome call in a female proband: check sex / sample identity")
        return "y_hemizygous", ("paternal" if carries(father) else None), flags
    if ctype == "x_nonpar" and (sex == "male" or (sex is None and pz == "hemi")):
        if pz == "het":
            flags.append("heterozygous call on male X outside the PARs (mosaic, XXY, or artefact)")
        if sex is None:
            flags.append("haploid X call and --sex not given: treated as male")
        origin = None
        if has_m:
            if carries(mother):
                origin = "maternal"
            elif mother["adequate_ref"]:  # type: ignore[index]
                origin = "de_novo"
            else:
                flags.append("mother not adequately genotyped (" + "; ".join(mother["why"]) + ")")  # type: ignore[index]
        if carries(father):
            flags.append("father carries the allele, but X passes mother→son: check pedigree / sample labels")
        return "x_hemizygous", origin, flags
    z = pz
    if z == "hemi":
        flags.append("haploid call on X in a female proband: check the declared sex / caller ploidy"
                     if ctype == "x_nonpar" else f"haploid call on {ctype} chromosome")
        z = "hom_alt"
    if z == "hom_alt":
        x_unknown_sex = ctype == "x_nonpar" and sex is None
        if x_unknown_sex:
            # a diploid 1/1 call on male X is the caller's default, not a
            # homozygote: never read it as a Mendelian conflict
            flags.append("hom-alt call on X outside the PARs with the proband's sex unknown: in a male this is a "
                         "hemizygous call (callers emit diploid X by default). Pass --sex, or let zebra infer it "
                         "from the X/Y genotypes")
        for role, p in (("mother", mother), ("father", father)):
            if p is None:
                continue
            if p["adequate_ref"]:
                if x_unknown_sex:
                    continue  # the call may be a male hemizygote: say that instead of blaming a parent
                # for a female proband a hom-ref father IS a conflict (she has his X);
                # a male proband never reaches here, he is classed x_hemizygous above
                flags.append(f"Mendelian conflict: {role} is hom-ref (UPD, deletion in trans, or sample mix-up)")
            elif p["z"] == "hom_alt":
                flags.append(f"{role} is also hom-alt")
            elif p["z"] in ("missing", "other"):
                flags.append(f"{role} genotype {p['z']}: carrier status unknown (not 'not carried')")
            else:
                flags.append(f"{role} not adequately genotyped (" + "; ".join(p["why"]) + "): carrier status unknown")
        return ("hom_recessive" if (has_m or has_f) else "hom"), None, flags
    if not has_m and not has_f:
        return "het", None, flags
    mc, fc = carries(mother), carries(father)
    if has_m and has_f:
        if not mc and not fc:
            if mother["adequate_ref"] and father["adequate_ref"]:  # type: ignore[index]
                blind = [r for r, p in (("mother", mother), ("father", father)) if p.get("no_gq")]
                if blind:
                    flags.append(", ".join(blind) + " reported no genotype quality (GQ): this de novo call rests on "
                                 "depth and allele counts alone")
                return "de_novo", "de_novo", flags
            why = [f"{r}: {', '.join(p['why'])}" for r, p in (("mother", mother), ("father", father))
                   if not p["adequate_ref"]]  # type: ignore[index]
            flags.append("parents not adequately genotyped for a de novo call (" + "; ".join(why) + ")")
            return "possible_de_novo", None, flags
        if mc and fc:
            flags.append("both parents carry the allele")
            return "inherited_het", "both", flags
        return "inherited_het", ("maternal" if mc else "paternal"), flags
    p, role = (mother, "mother") if has_m else (father, "father")
    if carries(p):
        return "inherited_het", ("maternal" if role == "mother" else "paternal"), flags
    if not p["adequate_ref"]:  # type: ignore[index]
        # absence can only be read off an adequate hom-ref call
        flags.append(f"the {role}'s genotype here is not usable (" + "; ".join(p["why"])  # type: ignore[index]
                     + "): carrier status unknown, not 'not carried'")
        return "het", None, flags
    flags.append(f"not carried by the {role}; the other parent was not tested (de novo or inherited)")
    return "het", None, flags


def _parent(call: Optional[Call], k: int, min_dp: int, min_gq: int) -> Optional[Dict[str, Any]]:
    if call is None:
        return None
    z = zygosity(call, k)
    dp = call.depth()
    _, alt_reads = call.reads(k)
    why: List[str] = []
    if z == "missing":
        why.append(f"GT {call.raw or '.'}")
    elif z == "other":
        why.append("carries a different ALT at this site")
    if dp is None:
        why.append("no depth (DP/AD)")
    elif dp < min_dp:
        why.append(f"DP {dp} < {min_dp}")
    if call.gq is not None and call.gq < min_gq:
        why.append(f"GQ {call.gq:g} < {min_gq}")
    if z == "hom_ref" and alt_reads is not None and alt_reads >= 2:
        why.append(f"{alt_reads} ALT reads (parental mosaicism?)")
    return {"z": z, "adequate_ref": z == "hom_ref" and not why, "why": why,
            "no_gq": call.gq is None}


SEX_X_ALT_MIN = 20  # X non-PAR ALT calls needed before a het fraction decides
SEX_SCAN_MAX = 3_000_000  # records read before the pre-scan gives up
SEX_Y_ALT_MIN = 2  # Y non-PAR ALT calls needed before Y is read as evidence
SEX_Y_FRACTION = 0.25  # ... and the share of called Y sites they must make up (noise is sparser)


def infer_sex(path: str, proband: str, assembly: Optional[str]) -> Dict[str, Any]:
    """Infer the proband's sex from its own X non-PAR and Y genotypes.

    GATK HaplotypeCaller calls male X diploid unless the ploidy is set, so a
    male hemizygous variant arrives as `1/1`. Without this, such a call reads
    as a homozygote with a "Mendelian conflict" against both parents.

    Evidence used: haploid X calls, ALT calls on Y outside the PARs, and the
    heterozygous fraction of X non-PAR ALT calls. Returns the call and the
    counts it rests on; `sex` is None when the evidence is not enough.
    """
    header, records = iter_records(path)
    if proband not in header.samples:
        records.close()
        return {"sex": None, "basis": "the proband is not a sample in the VCF", "counts": {}}
    i = header.samples.index(proband)
    c: Counter = Counter()
    n = 0
    complete = True
    try:
        for rec in records:
            n += 1
            if n > SEX_SCAN_MAX:
                complete = False
                break
            if rec.chrom not in ("X", "Y") or not rec.fmt or "GT" not in rec.fmt or i >= len(rec.sample_fields):
                continue
            ctype = _ctype(rec.chrom, rec.pos, assembly)
            if ctype not in ("x_nonpar", "y"):
                continue
            gt, _ = parse_gt(dict(zip(rec.fmt, rec.sample_fields[i].split(":"))).get("GT", "."))
            if not gt or all(a is None for a in gt):
                continue
            called = [a for a in gt if a is not None]
            has_alt = any(a > 0 for a in called)
            if ctype == "y":
                c["y_called"] += 1
                if has_alt:
                    c["y_alt"] += 1
                continue
            c["x_called"] += 1
            if len(gt) == 1:
                c["x_haploid"] += 1
            if has_alt:
                c["x_alt"] += 1
                if len(set(called)) > 1:
                    c["x_het"] += 1
    finally:
        records.close()
    counts = dict(c)
    counts["records_scanned"] = n
    counts["scan_complete"] = complete
    x_alt, x_het, y_alt, y_called = c["x_alt"], c["x_het"], c["y_alt"], c["y_called"]
    het_frac = (x_het / x_alt) if x_alt else None
    y_frac = (y_alt / y_called) if y_called else None
    # the counts every answer rests on, so no basis can assert what was not measured
    seen = (f"X non-PAR: {x_alt} ALT call(s), {x_het} heterozygous, {c['x_haploid']} haploid; "
            f"Y non-PAR: {y_alt} ALT of {y_called} called"
            + ("" if complete else f"; scan stopped after {n} records"))
    if c["x_haploid"] and c["x_haploid"] >= 0.5 * max(1, c["x_called"]):
        return {"sex": "male", "basis": f"{c['x_haploid']}/{c['x_called']} X non-PAR calls are haploid ({seen})",
                "counts": counts}
    # Y ALT calls appear as noise in female samples, so they must be a real
    # share of the called Y sites, and must not contradict a female-like X
    if y_alt >= SEX_Y_ALT_MIN and (y_frac is None or y_frac >= SEX_Y_FRACTION) \
            and (het_frac is None or het_frac < 0.1):
        return {"sex": "male", "basis": f"{y_alt} of {y_called} called Y non-PAR sites carry an ALT allele ({seen})",
                "counts": counts}
    if x_alt >= SEX_X_ALT_MIN and het_frac is not None:
        if het_frac < 0.05 and (y_frac is None or y_frac >= SEX_Y_FRACTION or y_alt == 0):
            return {"sex": "male", "basis": f"only {x_het} of {x_alt} X non-PAR ALT calls are heterozygous ({seen})",
                    "counts": counts}
        if het_frac > 0.2 and (y_alt < SEX_Y_ALT_MIN or (y_frac is not None and y_frac < SEX_Y_FRACTION)):
            return {"sex": "female",
                    "basis": f"{x_het} of {x_alt} X non-PAR ALT calls are heterozygous and Y shows no consistent "
                             f"ALT calls ({seen})", "counts": counts}
    return {"sex": None, "basis": f"not enough consistent evidence ({seen})", "counts": counts}


def _parent_ref_ok(mother: Optional[Dict[str, Any]], father: Optional[Dict[str, Any]]) -> bool:
    """True when a tested parent who does not carry the allele has an adequate hom-ref call.

    "Not carried" can only be concluded from a real reference call with enough
    depth and quality; a `./.` or a DP-4 call means unknown, never absent.
    """
    tested = [p for p in (mother, father) if p is not None]
    if not tested:
        return False
    for p in tested:
        if p["z"] in CARRIER:
            continue
        if not p["adequate_ref"]:
            return False
    return True


def triage(path: str, proband: str, mother: Optional[str] = None, father: Optional[str] = None,
           sex: Optional[str] = None, max_af: float = 0.01, min_dp: int = 10, min_gq: int = 20,
           assembly: Optional[str] = None, genes: Optional[Sequence[str]] = None, hpo_genes: bool = False,
           case_dir: Optional[str] = None, max_annotate: int = 1500, out: Optional[str] = None,
           hpo_terms: Optional[Sequence[str]] = None, max_af_dominant: Optional[float] = None,
           today: Optional[str] = None) -> Outcome:
    """Filter, annotate and rank the proband's variants. See module docstring and the zebra-reanalysis skill."""
    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []
    notes: List[str] = []
    if sex not in (None, "male", "female"):
        raise UsageError("--sex must be male or female")
    if not 0 < max_af <= 1:
        raise UsageError("--max-af must be in (0, 1]")
    dom_af = min(max_af, 0.0001) if max_af_dominant is None else max_af_dominant
    dom_af_explicit = max_af_dominant is not None

    header, _records = iter_records(path)
    try:
        for _ in _records:  # touch the stream: a corrupt gzip must be named as such here
            break
    finally:
        _records.close()
    if not header.samples:
        raise UsageError("the VCF has no sample columns (no genotypes): triage needs at least the proband")
    for role, s in (("proband", proband), ("mother", mother), ("father", father)):
        if s is not None and s not in header.samples:
            raise UsageError(f"{role} {s!r} is not a sample in the VCF (samples: {', '.join(header.samples[:20])})")
    if len({x for x in (proband, mother, father) if x}) < len([x for x in (proband, mother, father) if x]):
        raise UsageError("proband, mother and father must be different samples")
    build = guess_build(header)
    if assembly is None:
        if build["guess"] is None:
            raise UsageError("genome build unknown from the VCF header"
                             + (" (conflicting evidence: " + "; ".join(build["evidence"]) + ")" if build["conflict"] else "")
                             + ": pass --assembly GRCh38 or GRCh37")
        assembly = build["guess"]
    elif assembly not in ensembl.HOSTS:
        raise UsageError("--assembly must be GRCh38 or GRCh37")
    elif build["guess"] and build["guess"] != assembly:
        raise UsageError(f"--assembly {assembly} contradicts the VCF header, which says {build['guess']} "
                         f"({'; '.join(build['evidence'][:2])}); builds must never be mixed")
    if not build["guess"]:
        warnings.append(f"build not stated in the VCF header; using --assembly {assembly} as given")
    mode = "trio" if mother and father else ("duo" if (mother or father) else "singleton")
    if mode != "trio":
        notes.append(f"{mode} analysis: de novo variants cannot be called without both parents")
    sex_inference: Optional[Dict[str, Any]] = None
    if sex is None:
        sex_inference = infer_sex(path, proband, assembly)
        if sex_inference["sex"]:
            sex = sex_inference["sex"]
            notes.append(f"--sex not given: inferred {sex} from the proband's own genotypes "
                         f"({sex_inference['basis']}). Pass --sex to override.")
            warnings.append(f"proband sex was not given and was inferred as {sex} from the X/Y genotypes "
                            f"({sex_inference['basis']}); X non-PAR calls are classed accordingly")
        else:
            notes.append("--sex not given and it could not be inferred from the X/Y genotypes "
                         f"({sex_inference['basis']}): X variants are treated as diploid unless the call is haploid")
            warnings.append("proband sex unknown: pass --sex male|female before trusting any X non-PAR call "
                            "(a male diploid X call reads as a homozygote)")
    if case_dir:
        case_dir = os.path.expanduser(case_dir)

    # ---- phenotype and gene set (before reading the VCF: restriction happens while streaming)
    pheno = _phenotype(case_dir, hpo_terms, warnings, need=hpo_genes)
    if isinstance(genes, str):
        genes = read_gene_file(genes)
    gene_set: List[str] = list(genes or [])
    hpo_gene_list: List[str] = []
    if hpo_genes and pheno:
        hpo_gene_list = list(pheno["genes"])[:200]
        gene_set += [g for g in hpo_gene_list if g not in gene_set]
    annotations = header.annotations()
    has_info_ann = bool(annotations["vep_csq"] or annotations["snpeff_ann"] or annotations["gnomad_af_info_keys"])
    regions: Optional[RegionIndex] = None
    region_info: Dict[str, Any] = {}
    if gene_set:
        got = gene_regions(gene_set, assembly)
        sources.extend(got.sources)
        warnings.extend(got.warnings)
        found = got.result["regions"]
        if not found:
            raise UsageError("none of the gene symbols could be placed on the genome (Ensembl lookup); nothing to search")
        regions = RegionIndex(found)
        region_info = {"genes_requested": len(OrderedDict.fromkeys(gene_set)), "genes_placed": len(found),
                       "from_file": len(genes or []), "from_hpo": len(hpo_gene_list), "pad_bp": REGION_PAD,
                       "missing": got.result["missing"][:20]}
    prefer: Optional[Set[str]] = None
    if gene_set:
        prefer = {r["symbol"] for r in found.values()}

    # ---- stream the VCF
    idx = {s: i for i, s in enumerate(header.samples)}
    pi, mi, fi = idx[proband], (idx[mother] if mother else None), (idx[father] if father else None)
    counts: "OrderedDict[str, int]" = OrderedDict((k, 0) for k in (
        "records", "alleles", "symbolic_skipped", "filter_pass", "proband_carries_alt", "quality_pass"))
    filtered_names: Counter = Counter()
    by_class_seen: Counter = Counter()
    stored: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    dropped_region = info_af_drop = info_modifier_drop = in_region = info_pass = 0
    no_depth = 0
    implausible_af: Set[str] = set()
    over_budget: Counter = Counter()
    no_gt = off_contig = 0
    chrom_lengths = CHROM_LENGTHS.get(assembly, {})
    de_novo_like = proband_hets_trio = 0
    _, records = iter_records(path)
    for rec in records:
        counts["records"] += 1
        if not rec.fmt or "GT" not in rec.fmt:
            no_gt += 1
            continue
        if chrom_lengths.get(rec.chrom) and rec.pos > chrom_lengths[rec.chrom]:
            off_contig += 1
        pcall = parse_call(rec.fmt, rec.sample_fields[pi]) if pi < len(rec.sample_fields) else None
        for k, alt_raw in enumerate(rec.alts, 1):
            if is_symbolic(alt_raw):
                counts["symbolic_skipped"] += 1
                continue
            counts["alleles"] += 1
            if rec.filter not in ("PASS", "."):
                filtered_names[rec.filter] += 1
                continue
            counts["filter_pass"] += 1
            if pcall is None:
                continue
            pz = zygosity(pcall, k)
            if pz not in CARRIER:
                continue
            counts["proband_carries_alt"] += 1
            dp = pcall.depth()
            if dp is None:
                no_depth += 1
            elif dp < min_dp:
                continue
            if pcall.gq is not None and pcall.gq < min_gq:
                continue
            counts["quality_pass"] += 1
            pos, ref, alt = normalize_allele(rec.pos, rec.ref, alt_raw)
            end = pos + max(len(ref), 1) - 1
            hits: List[str] = []
            if regions is not None:
                hits = regions.hits(rec.chrom, pos - 1, end + 1)
                if not hits:
                    dropped_region += 1
                    continue
                in_region += 1
            info_ann = info_annotation(header, rec, k) if has_info_ann else None
            if info_ann is not None:
                implausible_af.update(info_ann.get("implausible_af_fields") or [])
                if info_ann["af"] is not None and info_ann["af"] > max_af:
                    info_af_drop += 1
                    continue
                if info_ann["impact"] == "MODIFIER" and regions is None:
                    info_modifier_drop += 1
                    continue
                info_pass += 1
            mcall = parse_call(rec.fmt, rec.sample_fields[mi]) if mi is not None and mi < len(rec.sample_fields) else None
            fcall = parse_call(rec.fmt, rec.sample_fields[fi]) if fi is not None and fi < len(rec.sample_fields) else None
            mp, fp = _parent(mcall, k, min_dp, min_gq), _parent(fcall, k, min_dp, min_gq)
            if mother and mp is None:
                mp = {"z": "missing", "adequate_ref": False, "why": ["no sample column"]}
            if father and fp is None:
                fp = {"z": "missing", "adequate_ref": False, "why": ["no sample column"]}
            ctype = _ctype(rec.chrom, pos, assembly)
            cls, origin, flags = _classify(ctype, sex, pz, mp, fp)
            if mode == "trio" and pz == "het" and ctype == "auto":
                proband_hets_trio += 1
                if cls in ("de_novo",):
                    de_novo_like += 1
            if dp is None:
                flags.append("no depth in the proband call")
            ab = pcall.ab(k)
            if pz == "het" and ab is not None and ab < 0.2:
                flags.append(f"low allele balance {ab:.2f} (mosaic or artefact?)")
            by_class_seen[cls] += 1
            bucket = stored[cls]
            if len(bucket) >= max_annotate:
                # budget: never more than max_annotate per class are kept in
                # memory, with or without a gene set (an exome-wide gene list
                # over a WGS VCF would otherwise hold millions of dicts)
                over_budget[cls] += 1
                continue
            bucket.append({
                "chrom": rec.chrom, "pos": pos, "ref": ref, "alt": alt, "variant": f"{rec.chrom}-{pos}-{ref}-{alt}",
                "vcf_id": rec.vid if rec.vid not in (".", "") else None, "class": cls, "origin": origin,
                "flags": flags, "ctype": ctype, "ab": round(ab, 3) if ab is not None else None,
                "region_genes": hits, "info": info_ann,
                "site": f"{rec.chrom}:{pos}:{ref}",
                "trans_at_site": pcall.other_alt(k),
                "phase": {"phased": pcall.phased, "ps": pcall.ps, "hap": pcall.haplotype(k), "gt": pcall.raw},
                "parent_ref_ok": _parent_ref_ok(mp, fp),
                "gt": {"proband": _gt_text(pcall, k), "mother": _gt_text(mcall, k) if mother else None,
                       "father": _gt_text(fcall, k) if father else None},
            })
    read_note = records.read_note()
    if read_note:
        # malformed lines are never dropped silently: half a file must not vanish without a word
        bad_fraction = records.malformed / max(1, records.data_lines)
        if bad_fraction > 0.5 or (bad_fraction > 0.1 and records.malformed > 20):
            raise UsageError(f"{path}: {read_note}. That is {bad_fraction:.0%} of the data lines; triage refuses to "
                             "report a negative result from a file it cannot read. Check that it is tab-separated "
                             "and complete (zebra vcf inspect reports the same count without refusing).")
        warnings.append(read_note)
    counts["malformed_lines"] = records.malformed
    if no_gt:
        counts["records_without_gt"] = no_gt
        notes.append(f"{no_gt} record(s) carry no GT in FORMAT and were not genotyped (they are counted in "
                     "'records' but in no later step)")
    if off_contig:
        warnings.append(f"{off_contig} record(s) sit past the end of their chromosome in {assembly}: the file's "
                        "positions do not match the build it is being read as")
    if over_budget:
        counts["over_per_class_cap"] = dict(over_budget)
    if implausible_af:
        warnings.append("annotation field(s) " + ", ".join(sorted(implausible_af)) + " held values above 1 and were "
                        "not "
                        "read as frequencies (they look like allele counts); the frequency prefilter ignored them")
    if filtered_names:
        notes.append("removed by FILTER: " + ", ".join(f"{k} {v}" for k, v in filtered_names.most_common(5)))
    if no_depth:
        warnings.append(f"{no_depth} proband calls carry no DP/AD: depth not checked for them")
    if mode == "trio" and proband_hets_trio >= 50 and de_novo_like / proband_hets_trio > 0.05:
        warnings.append(f"{de_novo_like}/{proband_hets_trio} autosomal proband hets look de novo (>5%): "
                        "check parent labels, sample swaps and parental coverage before trusting de novo calls")

    # ---- restriction before any variant leaves the machine
    survivors_all = [v for cls in list(PRIORITY) + [c for c in stored if c not in PRIORITY] for v in stored.get(cls, [])]
    # One allele written on two lines (a merged call set, or two spellings of the
    # same deletion) must never become two candidates: it would read as a
    # compound heterozygote of an allele with itself.
    seen_allele: Dict[Tuple[str, int, str, str], Dict[str, Any]] = {}
    duplicates = 0
    deduped = []
    for v in survivors_all:
        key = (v["chrom"], v["pos"], v["ref"], v["alt"])
        first = seen_allele.get(key)
        if first is None:
            seen_allele[key] = v
            deduped.append(v)
            continue
        duplicates += 1
        first["flags"].append("the same allele is written on more than one line of the VCF "
                              f"(also as {v['chrom']}-{v['pos']}-{v['ref']}-{v['alt']} from another record); "
                              "counted once, with the strongest inheritance class")
        if first.get("origin") is None and v.get("origin"):
            first["origin"] = v["origin"]
    if duplicates:
        warnings.append(f"{duplicates} candidate line(s) repeated an allele already counted (normalised to the same "
                        "chrom-pos-ref-alt); each allele is counted once")
    survivors_all = deduped
    total_candidates = sum(by_class_seen.values()) - duplicates
    restriction: Dict[str, Any] = {}
    if regions is not None:
        counts["in_gene_regions"] = in_region
        restriction["mode"] = "gene regions"
        restriction.update(region_info)
        restriction["outside_regions"] = dropped_region
    if has_info_ann:
        counts["info_prefilter_pass"] = info_pass
        restriction["info_prefilter"] = {"dropped_af_gt_max": info_af_drop, "dropped_modifier_only": info_modifier_drop,
                                         "fields": annotations}
        restriction.setdefault("mode", "VCF INFO annotations")
    if regions is None and not has_info_ann:
        restriction["mode"] = "prioritised by inheritance class (no gene set, no INFO annotations)"
        warnings.append("no --genes/--hpo-genes and no annotations in the VCF: annotating at most "
                        f"{max_annotate} variants, de novo / homozygous / hemizygous first")
    survivors = survivors_all
    if len(survivors) > max_annotate or total_candidates > len(survivors_all):
        kept = survivors[:max_annotate]
        skipped = total_candidates - len(kept)
        skipped_by = Counter()
        kept_ids = {id(v) for v in kept}
        for cls in PRIORITY:
            n_seen = by_class_seen.get(cls, 0)
            n_kept = sum(1 for v in stored.get(cls, []) if id(v) in kept_ids)
            if n_seen - n_kept:
                skipped_by[cls] = n_seen - n_kept
        survivors = kept
        counts["within_budget"] = len(kept)
        if skipped:
            restriction["skipped_over_budget"] = dict(skipped_by)
            warnings.append(f"{skipped} candidate variants NOT annotated (over --max-annotate {max_annotate}): "
                            + ", ".join(f"{k} {v}" for k, v in skipped_by.items())
                            + "; compound-heterozygous pairs among them can be missed — "
                            + ("narrow the gene list or raise --max-annotate" if regions is not None
                               else "give --genes/--hpo-genes"))
    counts["sent_to_vep"] = len({(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in survivors})

    # ---- annotate
    keys = list(OrderedDict.fromkeys((v["chrom"], v["pos"], v["ref"], v["alt"]) for v in survivors))
    vep = _vep_all(keys, assembly, warnings, sources)
    counts["annotated"] = sum(1 for k in keys if k in vep)
    unannotated = 0
    for v in survivors:
        rec = vep.get((v["chrom"], v["pos"], v["ref"], v["alt"]))
        if rec is not None:
            v["ann"] = annotate_vep(rec, prefer=(set(v["region_genes"]) if v["region_genes"] else prefer))
        elif v.get("info") and (v["info"].get("impact") or v["info"].get("af") is not None):
            i = v["info"]
            v["ann"] = {"gene": i.get("gene"), "consequence": i.get("consequence"), "impact": i.get("impact"),
                        "af": {"filter_af": i.get("af"), "source": "VCF INFO"}, "annotation": "VCF INFO ("
                        + "+".join(i.get("from") or []) + "; VEP unavailable)", "notes": []}
        else:
            v["ann"] = None
            unannotated += 1
    if unannotated:
        warnings.append(f"{unannotated} variants could not be annotated (VEP unavailable for them); they are kept, "
                        "unscored, at the end of the list")

    # ---- frequency: recessive threshold for everyone first (comp-het partners must survive it)
    def filt_af(v: Dict[str, Any]) -> Optional[float]:
        a = v.get("ann") or {}
        return (a.get("af") or {}).get("filter_af")

    def bound_af(v: Dict[str, Any]) -> Optional[float]:
        """The frequency a dominant cut-off may be applied to (see `af_summary`)."""
        af = (v.get("ann") or {}).get("af") or {}
        return af.get("bound_af") if af.get("bound_af") is not None else af.get("filter_af")

    def label(v: Dict[str, Any]) -> str:
        gene = v.get("gene") or (v.get("ann") or {}).get("gene") or (v["region_genes"] or [None])[0]
        return f"{v['variant']}{' ' + gene if gene else ''}"

    kept1: List[Dict[str, Any]] = []
    af_dropped: List[Dict[str, Any]] = []
    af_kept_plp: List[Dict[str, Any]] = []
    for v in survivors:
        af = filt_af(v)
        plp = clinvar_pathogenic(v.get("ann"))
        if af is None or af <= max_af:
            kept1.append(v)
            continue
        if clinvar_conflicting(v.get("ann")):
            warnings.append(f"{label(v)} has conflicting ClinVar classifications "
                            f"({'/'.join(v['ann']['clinvar'])}) and gnomAD {af:.3g} above --max-af {max_af:g}: "
                            "removed on frequency, because a contested assertion is not a reason to keep a common "
                            "allele. Read the ClinVar submissions if this gene fits the phenotype.")
        elif plp and af > BA1_AF:
            warnings.append(f"{label(v)} carries a ClinVar {'/'.join(plp)} assertion but its gnomAD frequency is "
                            f"{af:.3g}, above the ClinGen BA1 threshold {BA1_AF:g}: removed on frequency, and the "
                            "assertion is worth reporting to ClinVar")
        if plp and af <= BA1_AF:
            # the commonest pathogenic recessive alleles (CFTR F508del grpmax
            # 0.015, GJB2 c.35delG, HFE C282Y) sit above a 1 % cut-off: a
            # frequency filter must not delete the textbook diagnosis. Above
            # the ClinGen BA1 threshold (5 %) the exemption stops: no allele
            # that common can be a fully penetrant cause, whatever a single
            # ClinVar submission asserts.
            v["af_exempt"] = "ClinVar " + "/".join(plp)
            v["flags"].append(f"gnomAD {af:.3g} is above --max-af {max_af:g}, but ClinVar reports "
                              f"{'/'.join(plp)} for this allele: kept for review rather than filtered "
                              "(check the assertion's review status yourself)")
            af_kept_plp.append(v)
            kept1.append(v)
            continue
        af_dropped.append(v)
    counts["af_pass_recessive"] = len(kept1)
    if af_dropped:
        warnings.append(f"{len(af_dropped)} annotated variant(s) removed by gnomAD AF > --max-af {max_af:g}: "
                        + ", ".join(f"{label(v)} {filt_af(v):.3g}" for v in af_dropped[:8])
                        + (" …" if len(af_dropped) > 8 else ""))
    if af_kept_plp:
        warnings.append(f"{len(af_kept_plp)} variant(s) above --max-af {max_af:g} were kept because ClinVar reports "
                        "them pathogenic/likely pathogenic: "
                        + ", ".join(f"{label(v)} {filt_af(v):.3g}" for v in af_kept_plp[:8]))

    # ---- compound heterozygotes, then the stricter dominant threshold
    by_gene: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for v in kept1:
        g = (v.get("ann") or {}).get("gene") or (v["region_genes"][0] if v["region_genes"] else None)
        v["gene"] = g
        if g and v["class"] in DOMINANT_CLASSES and _qualifies(v.get("ann")):
            by_gene[g].append(v)
    comphet_genes: List[str] = []
    for g, vs in by_gene.items():
        if len(vs) < 2:
            continue
        # ---- phase: a PS block or a 1/2 call settles cis/trans without the parents.
        # The parents' genotypes outrank it: when a trio puts two alleles in
        # trans and a phase block says cis, the two disagree, and the candidate
        # is kept with the disagreement named rather than quietly dropped.
        trans_pairs, cis_pairs, unknown_pairs = _phase_pairs(vs)
        trio_mat = [v for v in vs if v["origin"] == "maternal"] if mode == "trio" else []
        trio_pat = [v for v in vs if v["origin"] == "paternal"] if mode == "trio" else []
        if trio_mat and trio_pat and cis_pairs:
            keep = []
            for a, b, why in cis_pairs:
                if {a["origin"], b["origin"]} == {"maternal", "paternal"}:
                    for v, other in ((a, b), (b, a)):
                        v["flags"].append(
                            f"the phase set puts this allele on the same haplotype as {other['variant']}, but the "
                            f"parents' genotypes put them in trans ({a['origin']} / {b['origin']}): check the "
                            "phasing and the sample labels; the parental genotypes were used")
                else:
                    keep.append((a, b, why))
            cis_pairs = keep
        if trans_pairs and not (trio_mat and trio_pat):
            comphet_genes.append(g)
            partners: Dict[str, List[str]] = defaultdict(list)
            for a, b, why in trans_pairs:
                partners[a["variant"]].append(b["variant"])
                partners[b["variant"]].append(a["variant"])
                for v in (a, b):
                    if why not in v["flags"]:
                        v["flags"].append(why)
            _pair_up(vs, lambda v: partners.get(v["variant"], []), "trans", "")
            for a, b, why in cis_pairs:
                for v in (a, b):
                    if why not in v["flags"]:
                        v["flags"].append(why)
            continue
        if cis_pairs and not unknown_pairs and not (trio_mat and trio_pat) and not trans_pairs:
            for a, b, why in cis_pairs:
                for v in (a, b):
                    if why not in v["flags"]:
                        v["flags"].append(why)
            warnings.append(f"{g}: {len(vs)} qualifying het(s) are on one haplotype according to the phase set, so "
                            "they are NOT called compound heterozygous and are held to the dominant allele-frequency "
                            "cut-off. If the phase set came from population phasing rather than reads, re-run with "
                            "--max-af-dominant to see them.")
            continue  # every pair is on one haplotype: not compound heterozygous
        for a, b, why in cis_pairs:
            for v in (a, b):
                if why not in v["flags"]:
                    v["flags"].append(why)
        for a, b, why in trans_pairs:
            for v in (a, b):
                if why not in v["flags"]:
                    v["flags"].append(why)
        for a, b, why in unknown_pairs:
            if not why:
                continue  # no phase information at all needs no explanation
            for v in (a, b):
                if why not in v["flags"]:
                    v["flags"].append(why)
        if mode == "trio":
            mat, pat = trio_mat, trio_pat
            if mat and pat:
                comphet_genes.append(g)
                _pair_up(mat + pat,
                         lambda v: [w["variant"] for w in (pat if v in mat else mat)], "trans",
                         "in trans: one allele from each parent")
            rest = [v for v in vs if not v.get("comphet")]
            if any(v["class"] in ("de_novo", "possible_de_novo") or v["origin"] == "both" for v in rest):
                # a de novo plus an inherited allele in a recessive gene is a
                # real compound heterozygote until the phase is tested; holding
                # the inherited partner to the dominant cut-off deletes it
                _pair_up(rest, lambda v: [w["variant"] for w in vs], "unphased",
                         f"{len(vs)} qualifying hets in {g}, one de novo or of unknown parental origin: "
                         "phase unknown (possible comp-het); the recessive AF threshold is used, not the dominant one")
                if g not in comphet_genes:
                    comphet_genes.append(g)
            elif not (mat and pat):
                for v in rest:
                    v["flags"].append(f"{len(rest)} hets in {g} all from the {'mother' if mat else 'father'} "
                                      "(in cis): not comp-het")
            continue
        carried = [v for v in vs if v["class"] == "inherited_het"]
        # only a parent with an adequate hom-ref call establishes "not carried";
        # a missing or low-quality parental genotype leaves it unknown
        not_carried = [v for v in vs if v["class"] == "het" and v.get("parent_ref_ok")]
        unknown_parent = [v for v in vs if v["class"] == "het" and not v.get("parent_ref_ok")]
        if mode == "duo" and carried and not_carried:
            comphet_genes.append(g)
            role = "mother" if mother else "father"
            _pair_up(carried + not_carried,
                     lambda v: [w["variant"] for w in (not_carried if v in carried else carried)], "trans",
                     f"in trans inferred from one parent: one allele from the {role}, the other not "
                     "(from the untested parent or de novo)")
            for v in vs:
                if not v.get("comphet"):
                    v["flags"].append(f"another qualifying het in {g}")
            continue
        comphet_genes.append(g)
        _pair_up(vs, lambda v: [w["variant"] for w in vs], "unphased",
                 "phase unknown: test the parents to confirm the alleles are in trans")
        if mode == "duo" and unknown_parent:
            for v in unknown_parent:
                v["flags"].append("the tested parent's genotype here is missing or below the depth/quality "
                                  "thresholds: carrier status unknown, so trans was NOT inferred")
    counts["comphet_genes"] = len(comphet_genes)
    final: List[Dict[str, Any]] = []
    dom_dropped: List[Dict[str, Any]] = []
    for v in kept1:
        if v["class"] not in DOMINANT_CLASSES:
            final.append(v)
            continue
        bound = bound_af(v)
        if bound is not None and bound > dom_af and not v.get("af_exempt"):
            dom_dropped.append(v)
            continue
        af = (v.get("ann") or {}).get("af") or {}
        if af.get("grpmax") is not None and af["grpmax"] > dom_af and (bound is None or bound <= dom_af):
            v["flags"].append(
                f"gnomAD {af['grpmax_group']} point estimate {af['grpmax']:.3g} is above the dominant cut-off "
                f"{dom_af:g} while the bound used ({af.get('bound_basis') or 'overall estimate'}) is "
                f"{'absent' if bound is None else f'{bound:.3g}'}: a single observation in a small group is not "
                "evidence against pathogenicity (check AC/AN or faf95 in gnomAD)")
        final.append(v)
    counts["af_pass_dominant"] = len(final)
    if dom_dropped:
        warnings.append(f"{len(dom_dropped)} variant(s) in a dominant class removed by gnomAD AF > "
                        f"{dom_af:g} (the dominant cut-off): "
                        + ", ".join(f"{label(v)} {bound_af(v):.3g} [{v['class']}]" for v in dom_dropped[:8])
                        + (" …" if len(dom_dropped) > 8 else "")
                        + "; raise --max-af-dominant to keep them")

    # ---- scoring
    moi = gene_moi({v["gene"] for v in final if v.get("gene")})
    use_pheno = pheno is not None
    weights = dict(WEIGHTS) if use_pheno else {"variant": WEIGHTS["variant"] / (1 - WEIGHTS["phenotype"]),
                                               "inheritance": WEIGHTS["inheritance"] / (1 - WEIGHTS["phenotype"])}
    for v in final:
        a = v.get("ann")
        g = v.get("gene")
        v["moi"] = moi.get(g or "", [])
        # a variant paired with another in the same gene is scored as half of a
        # recessive pair even when it keeps a stronger class label (a de novo)
        pair = v.get("comphet")
        scored_as = ("comphet" if pair == "trans" else "comphet_unphased") if pair else v["class"]
        fit = moi_fit(scored_as, v["moi"], v["ctype"])
        base = max(INH_BASE.get(v["class"], 0.3), INH_BASE.get(scored_as, 0.3))
        if v["class"] == "x_hemizygous" and v.get("origin") == "de_novo":
            base = 1.0
        v["moi_fit"] = fit
        v["inheritance_score"] = round(base * FIT_FACTOR[fit], 3)
        if use_pheno:
            ph = pheno["genes"].get(g or "")  # type: ignore[index]
            v["phenotype_score"] = ph["score"] if ph else 0.0
            v["phenotype_via"] = f"{ph['via']} {ph['via_name']}" if ph else None
        else:
            v["phenotype_score"] = None
            v["phenotype_via"] = None
        if a is None or a.get("impact") is None:
            v["consequence_score"] = v["predictor_score"] = v["variant_score"] = None
            v["support"] = []
            v["score"] = None
            continue
        sev = SEVERITY.get(a.get("impact"), 0.1)
        pred, support = predictor_support(a, assembly)
        v["consequence_score"] = sev
        v["predictor_score"] = pred
        v["support"] = support
        v["variant_score"] = round(max(sev, pred), 3)
        total = weights["variant"] * v["variant_score"] + weights["inheritance"] * v["inheritance_score"]
        if use_pheno:
            total += weights["phenotype"] * v["phenotype_score"]
        v["score"] = round(total, 4)
    order = {c: i for i, c in enumerate(CLASS_ORDER)}
    final.sort(key=lambda v: (v["score"] is None, -(v["score"] or 0), order.get(v["class"], 99), -(v.get("variant_score") or 0)))
    counts["ranked"] = len(final)

    # ---- outputs
    tsv_path = out
    if not tsv_path and case_dir and os.path.isfile(os.path.join(case_dir, "case.json")):
        tsv_path = os.path.join(case_dir, "reports", f"triage-{today or date.today().isoformat()}.tsv")
    thresholds = {
        "max_af_recessive": max_af,
        "max_af_dominant": dom_af,
        "dominant_rule": ("as given (--max-af-dominant)" if dom_af_explicit else "min(--max-af, 0.0001)")
        + " for de_novo, possible_de_novo, inherited_het, het; --max-af for hom/comphet/hemizygous classes",
        "af_used": "max(gnomAD exome, gnomAD genome, grpmax-like over afr/amr/eas/nfe/sas) as served by Ensembl VEP;"
                   " absent = passes",
        "min_dp": min_dp, "min_gq": min_gq,
    }
    if tsv_path:
        _write_tsv(tsv_path, final, thresholds, assembly, proband, mode)
    else:
        warnings.append("no TSV written: pass --out or work in a case (zebra --case <dir>)")
    if use_pheno:
        sources.append(source_record("HPO annotations (local)", pheno["hpo_version"], url="https://hpo.jax.org/data/annotations",  # type: ignore[index]
                                     note="phenotype fit (Resnik) and gene inheritance modes from phenotype.hpoa / genes_to_phenotype.txt"))
    elif moi:
        sources.append(source_record("HPO annotations (local)", "genes_to_phenotype.txt", url="https://hpo.jax.org/data/annotations",
                                     note="gene inheritance modes"))
    sources.insert(0, source_record("local VCF", os.path.basename(path), url=None,
                                    note="read on this machine; only candidate chrom-pos-ref-alt sent to VEP"))
    result = {
        "vcf": str(Path(os.path.expanduser(path)).resolve()),
        "proband": proband, "mother": mother, "father": father, "sex": sex, "mode": mode, "assembly": assembly,
        "sex_inferred": sex_inference if (sex_inference and sex_inference.get("sex")) else None,
        "counts": dict(counts),
        "restriction": restriction,
        "thresholds": thresholds,
        "weights": {k: round(w, 3) for k, w in weights.items()},
        "scoring": ("score = Σ weight × component. variant = max(consequence severity HIGH 1/MODERATE 0.6/LOW 0.25/"
                    "MODIFIER 0.1, in-silico level ≤0.9); inheritance = class prior × gene-MOI fit (fits 1, unknown 0.8,"
                    " against 0.5); phenotype = gene's relative Resnik score for the case's HPO profile"),
        "phenotype": ({k: pheno[k] for k in ("origin", "present", "excluded", "hpo_version", "method", "notes")}
                      if pheno else None),
        "comphet_genes": comphet_genes,
        "candidates": [_row(v, i + 1) for i, v in enumerate(final[:TOP_RESULT])],
        "total_candidates": len(final),
        "tsv": str(Path(tsv_path).resolve()) if tsv_path else None,
        "notes": notes + ["SNV/indel VCF only: CNVs, repeat expansions, mtDNA (unless called), low-level mosaicism "
                          "and poorly covered regions are not assessed"],
    }
    text = _render(result, final)
    query = {"vcf": path, "proband": proband, "mother": mother, "father": father, "sex": sex, "max_af": max_af,
             "max_af_dominant": dom_af, "min_dp": min_dp, "min_gq": min_gq, "assembly": assembly,
             "genes": len(genes or []), "hpo_genes": hpo_genes, "max_annotate": max_annotate, "out": tsv_path}
    return Outcome(result, sources=sources, warnings=warnings, text=text, query=query)


def _row(v: Dict[str, Any], rank: int) -> Dict[str, Any]:
    a = v.get("ann") or {}
    af = a.get("af") or {}
    hgvsc = a.get("hgvsc")
    return {
        "rank": rank, "score": v.get("score"),
        "components": {"phenotype": v.get("phenotype_score"), "variant": v.get("variant_score"),
                       "consequence": v.get("consequence_score"), "predictors": v.get("predictor_score"),
                       "inheritance": v.get("inheritance_score")},
        "variant": v["variant"], "gene": v.get("gene"), "transcript": a.get("transcript"), "mane": a.get("mane"),
        "hgvsc": hgvsc.split(":", 1)[1] if hgvsc and ":" in hgvsc else hgvsc,
        "hgvsp": (a.get("hgvsp") or "").split(":", 1)[-1].replace("%3D", "=") or None,
        "consequence": a.get("consequence"), "impact": a.get("impact"),
        "class": v["class"], "origin": v.get("origin"), "moi": v.get("moi"), "moi_fit": v.get("moi_fit"),
        "partners": v.get("partners"),
        "gnomad": {"filter_af": af.get("filter_af"), "exome": af.get("exome"), "genome": af.get("genome"),
                   "grpmax": af.get("grpmax"), "grpmax_group": af.get("grpmax_group"), "source": af.get("source"),
                   "bound_af": af.get("bound_af"), "bound_basis": af.get("bound_basis"),
                   "faf95": af.get("faf95"), "faf95_group": af.get("faf95_group")},
        "af_exempt": v.get("af_exempt"), "comphet": v.get("comphet"),
        "phase": {k: x for k, x in (v.get("phase") or {}).items() if x not in (None, False)} or None,
        "predictors": {"revel": a.get("revel"), "alphamissense": a.get("alphamissense"), "am_class": a.get("am_class"),
                       "cadd": a.get("cadd"), "spliceai_max": a.get("spliceai_max")},
        "support": v.get("support"), "clinvar": a.get("clinvar") or None, "rsid": a.get("rsid") or v.get("vcf_id"),
        "phenotype_via": v.get("phenotype_via"), "genotypes": {k: g for k, g in v["gt"].items() if g},
        "flags": v["flags"] + (a.get("notes") or []), "annotation": a.get("annotation"),
    }


TSV_COLUMNS = ("rank", "score", "phenotype_score", "variant_score", "consequence_score", "predictor_score",
               "inheritance_score", "variant", "gene", "class", "origin", "moi", "moi_fit", "consequence", "impact",
               "transcript", "mane", "hgvsc", "hgvsp", "gnomad_filter_af", "gnomad_exome_af", "gnomad_genome_af", "grpmax_af",
               "grpmax_group", "revel", "alphamissense", "am_class", "cadd", "spliceai_max", "support", "clinvar",
               "rsid", "proband_gt", "mother_gt", "father_gt", "allele_balance", "partners", "phenotype_via", "flags",
               "annotation", "gnomad_bound_af", "gnomad_bound_basis", "gnomad_faf95", "af_exempt", "phase_set")


def _fmt(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.4g}"
    if isinstance(x, (list, tuple)):
        return ";".join(str(i) for i in x)
    return str(x).replace("\t", " ").replace("\n", " ")


def _write_tsv(path: str, final: List[Dict[str, Any]], thresholds: Dict[str, Any], assembly: str, proband: str,
               mode: str) -> None:
    p = Path(os.path.expanduser(path))
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(f"## zebra vcf triage {date.today().isoformat()} | proband {proband} | {mode} | {assembly}\n")
        fh.write(f"## AF thresholds: recessive classes <= {thresholds['max_af_recessive']}; dominant classes "
                 f"(de_novo, possible_de_novo, inherited_het, het) <= {thresholds['max_af_dominant']} "
                 f"[{thresholds['dominant_rule'].split(' for ')[0]}]; {thresholds['af_used']}\n")
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        for i, v in enumerate(final, 1):
            r = _row(v, i)
            g, pr = r["gnomad"], r["predictors"]
            vals = [r["rank"], r["score"], r["components"]["phenotype"], r["components"]["variant"],
                    r["components"]["consequence"], r["components"]["predictors"], r["components"]["inheritance"],
                    r["variant"], r["gene"], r["class"], r["origin"], r["moi"], r["moi_fit"], r["consequence"],
                    r["impact"], r["transcript"], r["mane"], r["hgvsc"], r["hgvsp"], g["filter_af"], g["exome"], g["genome"],
                    g["grpmax"], g["grpmax_group"], pr["revel"], pr["alphamissense"], pr["am_class"], pr["cadd"],
                    pr["spliceai_max"], r["support"], r["clinvar"], r["rsid"], v["gt"].get("proband"),
                    v["gt"].get("mother"), v["gt"].get("father"), v.get("ab"), r["partners"], r["phenotype_via"],
                    r["flags"], r["annotation"], g["bound_af"], g["bound_basis"], g["faf95"], r["af_exempt"],
                    (v.get("phase") or {}).get("ps")]
            fh.write("\t".join(_fmt(x) for x in vals) + "\n")


def _render(result: Dict[str, Any], final: List[Dict[str, Any]]) -> str:
    c = result["counts"]
    t = result["thresholds"]
    who = result["proband"] + (f" (mother {result['mother']}" if result["mother"] else "") + \
        (f", father {result['father']})" if result["father"] else (")" if result["mother"] else ""))
    lines = [f"zebra vcf triage — {result['mode']} {who}, {result['assembly']}, sex {result['sex'] or 'not given'}"]
    steps = [f"{c['records']} records", f"{c['alleles']} ALT alleles", f"{c['filter_pass']} FILTER PASS",
             f"{c['proband_carries_alt']} carried by proband", f"{c['quality_pass']} DP≥{t['min_dp']}/GQ≥{t['min_gq']}"]
    if "in_gene_regions" in c:
        r = result["restriction"]
        steps.append(f"{c['in_gene_regions']} in {r.get('genes_placed')} gene regions ±{r.get('pad_bp')} bp")
    if "info_prefilter_pass" in c:
        steps.append(f"{c['info_prefilter_pass']} after INFO prefilter")
    if "within_budget" in c:
        steps.append(f"{c['within_budget']} within --max-annotate")
    steps += [f"{c['sent_to_vep']} sent to VEP", f"{c['annotated']} annotated",
              f"{c['af_pass_recessive']} AF≤{t['max_af_recessive']:g}",
              f"{c['comphet_genes']} comp-het gene(s)", f"{c['af_pass_dominant']} after dominant AF≤{t['max_af_dominant']:g}"]
    lines.append("steps: " + " → ".join(steps))
    lines.append(f"AF rule: dominant classes ≤ {t['max_af_dominant']:g} ({t['dominant_rule'].split(' for ')[0]}), "
                 f"recessive/hemizygous ≤ {t['max_af_recessive']:g}; {t['af_used']}")
    w = result["weights"]
    lines.append("score = " + " + ".join(f"{v:.2f}×{k}" for k, v in w.items())
                 + ("" if result["phenotype"] else "  (no phenotype profile: phenotype fit not scored)"))
    if result["phenotype"]:
        p = result["phenotype"]
        lines.append(f"phenotype: {len(p['present'])} HPO terms from {p['origin']} (HPO {p['hpo_version']})")
    for r in result["candidates"]:
        comp = r["components"]
        af = r["gnomad"]["filter_af"]
        af_txt = "absent" if af is None else f"{af:.2g}"
        hg = " ".join(x for x in (r["hgvsc"], r["hgvsp"]) if x)
        cls = r["class"] + (f"/{r['origin']}" if r["origin"] and r["origin"] != r["class"] else "")
        comp_txt = " | ".join(f"{lab} {comp[k]:.2f}" for k, lab in (("phenotype", "pheno"), ("variant", "var"),
                                                                       ("inheritance", "inh")) if comp[k] is not None)
        score = f"{r['score']:.3f}" if r["score"] is not None else "  -  "
        lines.append(f"{r['rank']:>2}. {score} {r['gene'] or '-':<9} {r['variant']:<22} {hg[:48]:<48} "
                     f"{(r['consequence'] or '?').split(',')[0]} {r['impact'] or '?'}  {cls}  gnomAD {af_txt}  "
                     f"MOI {','.join(r['moi'] or []) or '?'} ({r['moi_fit']})  [{comp_txt}]")
        extra = list(r["support"] or [])
        if r["clinvar"]:
            extra.append("ClinVar " + "/".join(r["clinvar"]))
        if r["partners"]:
            extra.append("partner " + ",".join(r["partners"]))
        extra += r["flags"]
        if extra:
            lines.append("      " + "; ".join(extra))
    if result["total_candidates"] > len(result["candidates"]):
        lines.append(f"… {result['total_candidates'] - len(result['candidates'])} more in the TSV")
    if not result["candidates"]:
        lines.append("no candidate passed the filters (no candidate ≠ no genetic cause)")
    lines.append(f"TSV: {result['tsv'] or '(not written)'}")
    lines.extend(f"note: {n}" for n in result["notes"])
    return "\n".join(lines)
