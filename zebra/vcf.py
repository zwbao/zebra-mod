"""Local VCF reanalysis: read, inspect and triage a singleton, duo or trio VCF (PED-aware).

The VCF never leaves this machine. What goes out is only what annotation needs:
gene symbols to Ensembl lookup (to find gene regions) and candidate variants as
chrom-pos-ref-alt (no sample names, genotypes or depths) to Ensembl VEP; with
`--prefilter myvariant` every quality-passing allele's key also goes to
MyVariant.info (the whole-exome path), and with `--s2f-top N` up to N keys go
to the Broad SpliceAI-lookup. `result.sent_off_machine` counts each.

Reading: plain text or gzip/bgzip (bgzip is multi-member gzip, which the
stdlib reads). Multi-allelic records are split per ALT; alleles are trimmed to
a minimal representation (no left-alignment: that needs the reference);
`chr` prefixes are dropped and chrM becomes MT (Ensembl naming).

Triage (see `triage`): quality -> inheritance class -> restriction to a gene
set, or the MyVariant.info frequency prefilter (whole exome), or a capped,
prioritised subset -> VEP annotation and gnomAD frequency -> comp-het
resolution -> sibling segregation (PED) -> score (phenotype fit, consequence /
predictors, inheritance fit) -> SpliceAI/Pangolin on the best non-coding /
splice-region / synonymous candidates -> ranked TSV, with every stage counted
and timed (`result.counts`, `result.timings`).
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


GVCF_PLACEHOLDERS = ("<NON_REF>", "<*>", "<X>")  # gVCF "any other allele" placeholders, not variants


def is_structural(alt: str) -> bool:
    """A symbolic SV/CNV allele (<DEL>, <DUP:TANDEM>, <CNV>, a breakend) — not `*`, `.` or a gVCF <NON_REF>/<*>."""
    if alt.upper() in GVCF_PLACEHOLDERS:
        return False
    return alt.startswith("<") or "[" in alt or "]" in alt


def sv_type(alt: str, info: Dict[str, str]) -> str:
    """SVTYPE from INFO when present, else read off the allele token (<DUP:TANDEM> → DUP, a breakend → BND)."""
    declared = (info.get("SVTYPE") or "").strip()
    if declared:
        return declared.upper()
    if "[" in alt or "]" in alt:
        return "BND"
    token = alt.strip("<>").split(":")[0].upper()
    return token or "SV"


def sv_call(rec: "Record", k: int, sample_index: Optional[int]) -> Dict[str, Any]:
    """One symbolic allele as the `zebra cnv` command would take it, with the proband's genotype."""
    info = rec.info()
    alt = rec.alts[k - 1]
    kind = sv_type(alt, info)
    end = _int(info.get("END"))
    if end is None and _int(info.get("SVLEN")) is not None:
        end = rec.pos + abs(_int(info.get("SVLEN")) or 0)
    gt = None
    cn = _int(info.get("CN"))
    if sample_index is not None and sample_index < len(rec.sample_fields) and rec.fmt:
        fields = dict(zip(rec.fmt, rec.sample_fields[sample_index].split(":")))
        gt = fields.get("GT")
        if fields.get("CN") not in (None, "", "."):
            cn = _int(fields.get("CN"))
    direction = None
    if kind.startswith("DEL") or (cn is not None and cn < 2 and kind in ("CNV", "DEL")):
        direction = "loss"
    elif kind.startswith("DUP") or (cn is not None and cn > 2 and kind in ("CNV", "DUP")):
        direction = "gain"
    out: Dict[str, Any] = {"chrom": rec.chrom, "pos": rec.pos, "end": end, "svtype": kind, "alt": alt, "gt": gt,
                           "cn": cn, "filter": rec.filter}
    if end is not None and direction and kind != "BND":
        out["zebra_cnv"] = f'zebra cnv "chr{rec.chrom}:{rec.pos}-{end} {direction}"'
    return out


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
    sv_types: Counter = Counter()
    star = small = 0
    for rec in records:
        n += 1
        if len(rec.alts) > 1:
            multi += 1
        for alt in rec.alts:
            if is_structural(alt):
                sv_types[sv_type(alt, rec.info())] += 1
            elif alt == "*":
                star += 1
            elif not is_symbolic(alt):
                small += 1
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
        "symbolic_alleles": {"structural": sum(sv_types.values()), "types": dict(sv_types.most_common(10)),
                             "spanning_deletion_star": star, "snv_indel_alleles": small},
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
    if sv_types:
        n_sv = sum(sv_types.values())
        notes.append(f"{n_sv} symbolic CNV/SV allele(s) ({', '.join(f'{k} {v}' for k, v in sv_types.most_common(6))})"
                     + (": this is a CNV/SV VCF" if not small else "")
                     + " — `vcf triage` reads SNVs/indels only; give each CNV/SV call to `zebra cnv`")
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
        # the same scoring as `zebra phenotype rank` (excluded terms flagged, not scored: docs/BENCHMARK.md);
        # with the 0.1.0 penalty the true gene reached the top 10 in 21 of 60 benchmark cases, with this one 48
        from zebra.commands.phenotype import LOCAL_PARAMS

        res = hpo_local.rank(idx, present, excluded, top=1000, params=dict(LOCAL_PARAMS))
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


def _sig(s: Any) -> str:
    """A ClinVar significance in VEP's spelling: 'Pathogenic, low penetrance' -> 'pathogenic_low_penetrance'."""
    return re.sub(r"[\s,]+", "_", str(s).strip().lower())


def clinvar_pathogenic(ann: Optional[Dict[str, Any]]) -> List[str]:
    """The P/LP assertions VEP's colocated ClinVar record (or MyVariant's) carries for this allele.

    A record that also carries benign or likely benign, or that ClinVar itself
    calls conflicting ("Conflicting interpretations/classifications of
    pathogenicity"), is not an assertion anyone can lean on, so it returns
    nothing: see `clinvar_conflicting`.
    """
    sigs = [_sig(s) for s in (ann or {}).get("clinvar") or []]
    if any(s in CLINVAR_BENIGN or s.startswith("conflicting") for s in sigs):
        return []
    return [str(s) for s in (ann or {}).get("clinvar") or [] if _sig(s) in CLINVAR_PLP]


def clinvar_conflicting(ann: Optional[Dict[str, Any]]) -> bool:
    """True when the ClinVar record carries a pathogenic and a benign classification, or ClinVar says conflicting."""
    sigs = [_sig(s) for s in (ann or {}).get("clinvar") or []]
    return (any(s in CLINVAR_PLP for s in sigs) and any(s in CLINVAR_BENIGN for s in sigs)) \
        or any(s.startswith("conflicting") for s in sigs)


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

REF_CHECK_N = 50  # positions compared with the reference (one POST /sequence/region takes up to 50)


def ref_check(keys: Sequence[Tuple[str, int, str, str]], assembly: str, n: int = REF_CHECK_N) -> Outcome:
    """Compare the VCF's REF with the reference genome at up to `n` positions spread over `keys`."""
    pool = sorted({k for k in keys if 0 < len(k[2]) <= 20 and k[0] != "MT"}, key=lambda k: (_chrom_sort_key(k[0]), k[1]))
    if not pool:
        return Outcome({"checked": 0, "mismatch": 0, "examples": []})
    step = max(1, len(pool) // n)
    sample = pool[::step][:n]
    regions = [f"{c}:{p}..{p + len(r) - 1}:1" for c, p, r, _ in sample]
    resp = post_json(f"{ensembl.host(assembly)}/sequence/region/human", {"regions": regions},
                     source="Ensembl sequence", cache_ttl=90 * 86400, timeout=60)
    data = resp.json()
    if not isinstance(data, list):
        raise ValueError("Ensembl sequence/region returned no list")
    got = {str(d.get("query")): str(d.get("seq") or "").upper() for d in data if isinstance(d, dict)}
    checked = mismatch = 0
    examples: List[str] = []
    for (c, p, r, _), q in zip(sample, regions):
        seq = got.get(q)
        if not seq:
            continue
        checked += 1
        if seq != r.upper():
            mismatch += 1
            if len(examples) < 5:
                examples.append(f"{c}:{p} VCF {r} vs {assembly} {seq}")
    rec = source_record("Ensembl sequence", f"{checked} REF positions", resp,
                        note=f"{assembly}; REF check of positions only (no sample data)")
    return Outcome({"checked": checked, "mismatch": mismatch, "examples": examples, "asked": len(sample)},
                   sources=[rec])


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
    lookup = a.get("spliceai_lookup")
    if lookup is not None and (sai is None or lookup > sai):
        # the Broad lookup's raw score at ±500 nt (zebra s2f), used when it says more than VEP's precomputed one
        if lookup >= 0.5:
            levels.append((0.85, f"SpliceAI {lookup:.2f} ≥0.5 (Broad lookup, raw, ±500 nt)"))
        elif lookup >= 0.2:
            levels.append((0.65, f"SpliceAI {lookup:.2f} ≥0.2 (Broad lookup, raw, ±500 nt)"))
    elif sai is not None:
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
        carriers = 0
        for role, p in (("mother", mother), ("father", father)):
            if p is None:
                continue
            if p["adequate_ref"]:
                if x_unknown_sex:
                    continue  # the call may be a male hemizygote: say that instead of blaming a parent
                # for a female proband a hom-ref father IS a conflict (she has his X);
                # a male proband never reaches here, he is classed x_hemizygous above
                flags.append(f"Mendelian conflict: {role} is hom-ref (UPD, deletion in trans, or sample mix-up)")
            elif p["z"] == "hom_alt" and not (ctype == "x_nonpar" and role == "father"):
                carriers += 1
                flags.append(f"{role} is also hom-alt")
            elif p["z"] in ("het", "hemi", "hom_alt"):
                # the expected carrier parent of a recessive homozygote (B-P1-7); on X a
                # father's 1/1 is his hemizygous allele written diploid
                carriers += 1
            elif p["z"] in ("missing", "other"):
                flags.append(f"{role} genotype {p['z']}: carrier status unknown (not 'not carried')")
            else:
                flags.append(f"{role} not adequately genotyped (" + "; ".join(p["why"]) + "): carrier status unknown")
        # with the sex unknown a male's X call reads as 1/1 here, and a hemizygote has one parent, not two
        origin = "biparental" if (has_m and has_f and carriers == 2 and not x_unknown_sex) else None
        return ("hom_recessive" if (has_m or has_f) else "hom"), origin, flags
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
SEX_SCAN_MAX = 3_000_000  # X/Y records read before the scan gives up (autosomal records do not count: B-P2-5)
SEX_SCAN_SECONDS = 240.0  # wall-time cap on the sex scan, so a whole genome cannot stall a triage
SEX_Y_ALT_MIN = 2  # Y non-PAR ALT calls needed before Y is read as evidence
SEX_Y_FRACTION = 0.25  # ... and the share of called Y sites they must make up (noise is sparser)
_XY_NAMES = frozenset(("X", "Y", "chrX", "chrY", "x", "y", "CHRX", "CHRY", "chrx", "chry"))


def tabix_starts(path: str) -> Optional[Dict[str, int]]:
    """The first BGZF virtual offset of each sequence's records, read from `<vcf>.tbi`; None when unusable.

    Format: SAMtools tabix specification (TBI\\1 magic, n_ref, 6 int32 fields, l_nm and the names, then per
    reference its bins with their chunks and the linear index). The metadata pseudo-bin 37450 is skipped.
    Any parsing problem returns None and the caller reads the file from the start instead.
    """
    import struct

    tbi = os.path.expanduser(path) + ".tbi"
    if not os.path.isfile(tbi) or compression(os.path.expanduser(path)) != "bgzip":
        return None
    try:
        with gzip.open(tbi, "rb") as fh:
            data = fh.read(256 << 20)
        if data[:4] != b"TBI\x01":
            return None
        n_ref, _fmt, _cs, _cb, _ce, _meta, _skip, l_nm = struct.unpack_from("<8i", data, 4)
        if not 0 < n_ref < 100_000 or not 0 <= l_nm < len(data):
            return None
        names = data[36:36 + l_nm].split(b"\x00")[:n_ref]
        off = 36 + l_nm
        out: Dict[str, int] = {}
        for r in range(n_ref):
            (n_bin,) = struct.unpack_from("<i", data, off)
            off += 4
            best: Optional[int] = None
            for _ in range(n_bin):
                bin_id, n_chunk = struct.unpack_from("<Ii", data, off)
                off += 8
                for _ in range(n_chunk):
                    beg, _end = struct.unpack_from("<QQ", data, off)
                    off += 16
                    if bin_id != 37450 and (best is None or beg < best):
                        best = beg
            (n_intv,) = struct.unpack_from("<i", data, off)
            off += 4 + 8 * n_intv
            if best is not None and r < len(names):
                out[names[r].decode("utf-8", "replace")] = best
        return out
    except (OSError, EOFError, zlib.error, struct.error, ValueError):
        return None


def _lines_from(path: str, voffset: int) -> Iterator[str]:
    """Text lines of a BGZF file starting at a virtual offset (compressed block start << 16 | offset within it)."""
    raw = open(os.path.expanduser(path), "rb")
    try:
        raw.seek(voffset >> 16)
        gz = gzip.GzipFile(fileobj=raw, mode="rb")
        gz.read(voffset & 0xFFFF)
        for line in io.TextIOWrapper(gz, encoding="utf-8", errors="replace"):
            yield line
    finally:
        raw.close()


def sex_counts(path: str, samples: Sequence[str], assembly: Optional[str], max_xy: Optional[int] = None,
               max_seconds: Optional[float] = None) -> Tuple[Dict[str, Counter], Dict[str, Any]]:
    """X non-PAR / Y genotype counts for each named sample, in one pass over the file.

    Only X and Y records are parsed and only they count against `max_xy`
    (B-P2-5): in a sorted whole-genome VCF chrX follows 4-5 million autosomal
    records, and a cap that counted those stopped before X was reached. The
    other lines are skipped on their CHROM field without being split. With a
    usable tabix index (not older than the VCF, and pointing at lines of the
    right chromosome) the scan jumps straight to X and Y; any doubt about the
    index, or any read error through it, falls back to reading the file.
    """
    import time

    max_xy = SEX_SCAN_MAX if max_xy is None else max_xy
    max_seconds = SEX_SCAN_SECONDS if max_seconds is None else max_seconds
    t0 = time.monotonic()
    fh = open_text(path)
    try:
        try:
            header, first = read_header(fh)
        except (UnicodeDecodeError, EOFError, OSError, zlib.error) as err:
            raise UsageError(f"{path} is not a readable VCF ({type(err).__name__}: {err})") from None
        idx = {s: header.samples.index(s) for s in samples if s in header.samples}
        n_samples = len(header.samples)

        def count(lines: Iterator[str], meta: Dict[str, Any]) -> Dict[str, Counter]:
            out: Dict[str, Counter] = {s: Counter() for s in samples}
            for line in lines:
                if not line or line[0] == "#":
                    continue
                meta["records_scanned"] += 1
                tab = line.find("\t")
                if line[:tab] not in _XY_NAMES:
                    if meta["records_scanned"] % 200_000 == 0 and time.monotonic() - t0 > max_seconds:
                        meta.update(scan_complete=False, stopped_by=f"the {max_seconds:g} s time cap")
                        break
                    continue
                meta["xy_records"] += 1
                if meta["xy_records"] > max_xy:
                    meta.update(scan_complete=False, stopped_by=f"the cap of {max_xy} X/Y records")
                    break
                if meta["xy_records"] % 20_000 == 0 and time.monotonic() - t0 > max_seconds:
                    meta.update(scan_complete=False, stopped_by=f"the {max_seconds:g} s time cap")
                    break
                rec, _ = parse_line_why(line, n_samples)
                if rec is not None:
                    sex_tally(rec, idx, out, assembly)
            return out

        def fresh_meta() -> Dict[str, Any]:
            return {"records_scanned": 0, "xy_records": 0, "scan_complete": True, "stopped_by": None}

        xy_starts = _usable_xy_index(path)
        if xy_starts:
            def indexed() -> Iterator[str]:
                for name, off in sorted(xy_starts.items(), key=lambda kv: kv[1]):
                    for line in _lines_from(path, off):
                        if line[:line.find("\t")] != name:
                            break
                        yield line

            meta = fresh_meta()
            meta["indexed"] = True
            try:
                out = count(indexed(), meta)
                meta["seconds"] = round(time.monotonic() - t0, 2)
                return out, meta
            except (EOFError, OSError, zlib.error, ValueError):
                pass  # a stale or damaged index: read the file instead

        def linear() -> Iterator[str]:
            if first is not None:
                yield first
            while True:
                try:
                    line = next(fh)
                except StopIteration:
                    return
                except (EOFError, OSError, zlib.error) as err:
                    raise UsageError(f"{path} is truncated or corrupt ({type(err).__name__}: {err})") from None
                yield line

        meta = fresh_meta()
        out = count(linear(), meta)
    finally:
        fh.close()
    meta["seconds"] = round(time.monotonic() - t0, 2)
    return out, meta


def _usable_xy_index(path: str) -> Dict[str, int]:
    """X/Y start offsets from the tabix index, only when the index is not older than the VCF and each offset
    lands on a line of that chromosome; otherwise {} (read the file)."""
    p = os.path.expanduser(path)
    try:
        if os.path.getmtime(p + ".tbi") < os.path.getmtime(p):
            return {}
    except OSError:
        return {}
    starts = tabix_starts(path) or {}
    xy = {name: off for name, off in starts.items() if name in _XY_NAMES}
    for name, off in xy.items():
        gen = _lines_from(path, off)
        try:
            line = next(gen, "")
        except (EOFError, OSError, zlib.error, ValueError):
            return {}
        finally:
            gen.close()
        if line[:line.find("\t")] != name:
            return {}
    return xy


def sex_tally(rec: "Record", idx: Dict[str, int], out: Dict[str, Counter], assembly: Optional[str]) -> None:
    """Add one X/Y record to each sample's sex counts (any FILTER, any allele type, as `infer_sex` always did)."""
    if rec.chrom not in ("X", "Y") or not rec.fmt or "GT" not in rec.fmt:
        return
    ctype = _ctype(rec.chrom, rec.pos, assembly)
    if ctype not in ("x_nonpar", "y"):
        return
    gi = rec.fmt.index("GT")
    for s, i in idx.items():
        if i >= len(rec.sample_fields):
            continue
        parts = rec.sample_fields[i].split(":")
        gt, _ = parse_gt(parts[gi] if gi < len(parts) else ".")
        if not gt or all(a is None for a in gt):
            continue
        c = out[s]
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


def sex_from_counts(c: Counter, meta: Dict[str, Any]) -> Dict[str, Any]:
    """The sex call one sample's X/Y counts support, with the counts it rests on; None when they are not enough."""
    complete = meta.get("scan_complete", True)
    counts = dict(c)
    counts["records_scanned"] = meta.get("records_scanned", 0)
    counts["xy_records"] = meta.get("xy_records", 0)
    counts["scan_complete"] = complete
    x_alt, x_het, y_alt, y_called = c["x_alt"], c["x_het"], c["y_alt"], c["y_called"]
    het_frac = (x_het / x_alt) if x_alt else None
    y_frac = (y_alt / y_called) if y_called else None
    # the counts every answer rests on, so no basis can assert what was not measured
    seen = (f"X non-PAR: {x_alt} ALT call(s), {x_het} heterozygous, {c['x_haploid']} haploid; "
            f"Y non-PAR: {y_alt} ALT of {y_called} called"
            + ("" if complete else f"; scan stopped by {meta.get('stopped_by')} after "
                                   f"{meta.get('records_scanned')} records"))
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
    records.close()
    if proband not in header.samples:
        return {"sex": None, "basis": "the proband is not a sample in the VCF", "counts": {}}
    counts, meta = sex_counts(path, [proband], assembly)
    return sex_from_counts(counts[proband], meta)


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


def _resolve_ped(ped: str, samples: Sequence[str], proband: Optional[str], mother: Optional[str],
                 father: Optional[str], sex: Optional[str], notes: List[str]
                 ) -> Tuple[str, Optional[str], Optional[str], Optional[str], Dict[str, Any], List[Tuple[str, Optional[bool], Optional[str]]]]:
    """Proband, parents, sex and full siblings from a PED file; a flag given on the command line must agree with it."""
    from zebra import qc as qc_mod

    pedigree = qc_mod.read_ped(ped)
    if proband is None:
        affected = [p.iid for p in pedigree.people.values() if p.affected and p.iid in samples
                    and ((p.mother in samples) or (p.father in samples))]
        if len(affected) != 1:
            raise UsageError("--proband not given and the PED does not name exactly one affected individual with a "
                             f"parent in the VCF (affected with a parent here: {', '.join(affected) or 'none'}): "
                             "pass --proband")
        proband = affected[0]
    person = pedigree.get(proband)
    if person is None:
        raise UsageError(f"--proband {proband!r} is not in the PED file ({pedigree.path})")
    for role, given, stated in (("mother", mother, person.mother), ("father", father, person.father)):
        if given and stated and given != stated:
            raise UsageError(f"--{role} {given} contradicts the PED file, which names {stated} as {proband}'s {role}")
    if mother is None and person.mother:
        if person.mother in samples:
            mother = person.mother
        else:
            notes.append(f"the PED names {person.mother} as the mother, but she is not a sample in the VCF")
    if father is None and person.father:
        if person.father in samples:
            father = person.father
        else:
            notes.append(f"the PED names {person.father} as the father, but he is not a sample in the VCF")
    if sex and person.sex and sex != person.sex:
        raise UsageError(f"--sex {sex} contradicts the PED file, which gives {proband} sex {person.sex}")
    sex = sex or person.sex
    sibs = [(s, pedigree.get(s).affected, pedigree.get(s).sex)  # type: ignore[union-attr]
            for s in pedigree.full_sibs(proband) if s in samples]
    info = {"path": pedigree.path, "proband": proband, "mother": mother, "father": father, "sex": sex,
            "affected": person.affected, "siblings": [{"id": s, "affected": a, "sex": x} for s, a, x in sibs]}
    return proband, mother, father, sex, info, sibs


# consequences the splice models are asked about (non-coding / splice-region / synonymous); canonical
# donor/acceptor (HIGH, PVS1's domain) and missense (a protein question) are not
S2F_ELIGIBLE = frozenset((
    "splice_region_variant", "splice_donor_5th_base_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant", "synonymous_variant", "intron_variant",
    "non_coding_transcript_exon_variant", "start_retained_variant", "stop_retained_variant",
    "coding_sequence_variant"))
# where only a regulatory model (AlphaGenome, with the disease tissue) speaks
REGULATORY_TERMS = frozenset((
    "5_prime_UTR_variant", "3_prime_UTR_variant", "upstream_gene_variant", "downstream_gene_variant",
    "regulatory_region_variant", "TF_binding_site_variant", "intergenic_variant"))


def _s2f_rerank(final: List[Dict[str, Any]], assembly: str, top: int, max_seconds: float, warnings: List[str],
                sources: List[Dict[str, Any]], progress: Optional[Any]) -> Dict[str, Any]:
    """SpliceAI and Pangolin (zebra.s2f, Broad lookup) on the best-ranked non-coding/splice-region/synonymous
    candidates; the SpliceAI delta is put where the ranking reads it (`spliceai_lookup`)."""
    import time

    from zebra import s2f as s2f_mod
    from zebra.http import deadline_seconds

    t0 = time.monotonic()
    chosen = []
    lof = {"splice_donor_variant", "splice_acceptor_variant", "stop_gained", "frameshift_variant", "start_lost"}
    protein = {"missense_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant", "stop_lost"}
    for v in final:
        terms = set(((v.get("ann") or {}).get("consequence") or "").split(","))
        splice_region = {t for t in terms if t.startswith("splice_")} - lof
        # a LoF is PVS1's question and a plain missense a protein question; an exonic variant in the splice
        # region (the last bases of an exon) is asked about splicing as well
        if terms & S2F_ELIGIBLE and not (terms & lof) and (splice_region or not (terms & protein)):
            # the lookup service knows the primary chromosomes only; a call on GL000220.1, a decoy or an
            # HLA contig is left out here rather than failing the whole triage
            if ensembl.parse_vcf_like(f"{v['chrom']}-{v['pos']}-{v['ref']}-{v['alt']}") is None:
                continue
            chosen.append(v)
        if len(chosen) >= top:
            break
    ran: List[str] = []
    skipped: List[str] = []
    failed: List[str] = []
    asked: List[str] = []
    for n, v in enumerate(chosen):
        left = deadline_seconds()
        if (left is not None and left < 25) or time.monotonic() - t0 > max_seconds:
            skipped = [w["variant"] for w in chosen[n:]]
            why = "the call's time budget (ZEBRA_DEADLINE_MS)" if left is not None and left < 25 else \
                f"the {max_seconds:g} s cap on this stage"
            warnings.append(f"S2F: {len(skipped)} of {len(chosen)} selected candidate(s) not run within {why}: "
                            + ", ".join(skipped) + " — run zebra s2f predict on them")
            break
        if progress:
            progress(f"S2F {n + 1}/{len(chosen)}: {v['variant']}")
        key = f"{v['chrom']}-{v['pos']}-{v['ref']}-{v['alt']}"
        asked.append(v["variant"])
        try:
            out = attempt(f"S2F {key}", lambda: s2f_mod.predict(key, assembly=assembly,
                                                                models=["spliceai", "pangolin"], distance=500),
                          warnings)
        except UsageError as err:  # one unreadable variant costs its own row, never the triage
            warnings.append(f"S2F {key} not run: {err}")
            out = None
        if out is None:
            failed.append(v["variant"])
            continue
        sources.extend(out.sources)
        rows = {r.get("model"): r for r in (out.result or {}).get("models") or []}
        got: Dict[str, Any] = {}
        for model in ("spliceai", "pangolin"):
            r = rows.get(model) or {}
            if r.get("status") == "ran":
                h = r.get("headline") or {}
                got[model] = {"status": "ran", "score": h.get("score"), "value": h.get("value"),
                              "position": h.get("position"), "transcript": h.get("refseq") or h.get("transcript")}
            else:
                got[model] = {"status": r.get("status") or "not_run", "reason": r.get("reason")}
                warnings.append(f"S2F {key}: {model} {got[model]['status']} ({r.get('reason')})")
        got["source"] = "zebra s2f predict (Broad SpliceAI-lookup, raw scores, distance ±500 nt, GENCODE basic)"
        v["s2f"] = got
        val = (got.get("spliceai") or {}).get("value")
        if isinstance(val, (int, float)) and v.get("ann") is not None:
            v["ann"]["spliceai_lookup"] = float(val)
        ran.append(v["variant"])
    return {"selected": [v["variant"] for v in chosen], "ran": ran, "not_run": skipped, "failed": failed,
            "asked": asked, "seconds": round(time.monotonic() - t0, 2),
            "rule": f"the {top} best-ranked candidates whose consequence is splice-region, synonymous, intronic or "
                    "non-coding exon; SpliceAI's delta enters the variant score (Walker 2023 thresholds 0.2/0.5); "
                    "Pangolin is shown, not scored (no calibrated thresholds)"}


PREFILTER_TIERS = {0: "ClinVar P/LP or HIGH impact", 1: "MODERATE, splice-region, novel or phenotype gene",
                   2: "LOW (synonymous and other low-impact)", 3: "MODIFIER only (intronic, UTR, non-coding)",
                   4: "not looked up"}
_SPLICE_EFFECTS = ("splice_region_variant", "splice_donor", "splice_acceptor", "splice_donor_5th_base_variant",
                   "splice_donor_region_variant", "splice_polypyrimidine_tract_variant")


STRONG_CLASSES = ("de_novo", "possible_de_novo", "hom_recessive", "hom", "x_hemizygous", "mitochondrial",
                  "y_hemizygous")


def _myvariant_prefilter(cands: List[Dict[str, Any]], assembly: str, max_af: float, max_prefilter: int,
                         pheno_genes: Dict[str, Any], warnings: List[str], sources: List[Dict[str, Any]],
                         progress: Optional[Any], dom_af: float = 0.0001
                         ) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Drop alleles MyVariant.info shows to be common; order the rest for the VEP budget.

    Frequency: the recessive threshold (`--max-af`) on max(gnomAD exome, genome, grpmax over afr/amr/eas/nfe/
    sas) — the rule `af_summary` applies to VEP — with the ClinVar P/LP exemption up to BA1 (5%). The dominant
    threshold is applied later, on VEP's gnomAD, after compound heterozygotes are paired.
    Consequence (snpEff via MyVariant) only orders the VEP budget, it removes nothing: tier 0 ClinVar P/LP or
    HIGH; 1 MODERATE, splice-region, not in MyVariant (novel) or in a gene that fits the phenotype; 2 LOW;
    3 MODIFIER only. A variant MyVariant does not hold, or that a failed batch did not look up, passes.
    VEP order: tier 0; then the strong inheritance classes (de novo, homozygous, hemizygous) at any tier;
    then heterozygous classes tier 1, then 2-3; heterozygotes the dominant cut-off will remove (MyVariant
    frequency above it, and no second candidate in the gene to pair with) after those; unchecked last.
    """
    from zebra.sources import myvariant

    pool = cands[:max_prefilter]
    unqueried = cands[max_prefilter:]
    keys = [(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in pool]
    meta = attempt("MyVariant.info metadata", lambda: myvariant.metadata(assembly), warnings)
    if meta is not None:
        sources.extend(meta.sources)
    got = attempt("MyVariant.info", lambda: myvariant.batch(keys, assembly, progress=progress), warnings)
    if got is not None:
        sources.extend(got.sources)
        warnings.extend(got.warnings)
        records = got.result["records"]
        unanswered = set(got.result["unanswered"])
    else:
        records, unanswered = {}, set(keys)
    kept: List[Tuple[int, int, int, Dict[str, Any]]] = []
    dropped: List[Tuple[Dict[str, Any], float]] = []
    exempt = 0
    tiers: Counter = Counter()
    order = {c: i for i, c in enumerate(PRIORITY)}
    for n, v in enumerate(pool):
        key = (v["chrom"], v["pos"], v["ref"], v["alt"])
        info: Dict[str, Any] = {}
        if key in unanswered:
            tier, info["status"] = 4, "not looked up (batch failed or unanswered): not frequency-filtered"
        elif key not in records or records[key] is None:
            tier, info["status"] = 1, "not in MyVariant.info: no population frequency known (novel?)"
        else:
            rec = records[key]
            af = af_summary(rec["groups"])
            info.update(status="found", filter_af=af["filter_af"], af_source=af["source"], clinvar=rec["clinvar"],
                        snpeff_impact=rec["impact"], effects=rec["effects"][:4], genes=rec["genes"][:3])
            plp = clinvar_pathogenic(rec)
            if af["filter_af"] is not None and af["filter_af"] > max_af:
                if plp and af["filter_af"] <= BA1_AF and not clinvar_conflicting(rec):
                    exempt += 1
                    info["kept_because"] = f"ClinVar {'/'.join(plp)} (frequency {af['filter_af']:.3g} ≤ BA1 {BA1_AF:g})"
                else:
                    dropped.append((v, af["filter_af"]))
                    continue
            impact = rec["impact"]
            in_pheno = any((pheno_genes.get(g) or {}).get("score", 0) > 0 for g in rec["genes"])
            if plp or impact == "HIGH":
                tier = 0
            elif impact == "MODERATE" or impact is None or in_pheno \
                    or any(e.startswith(_SPLICE_EFFECTS) for e in rec["effects"]):
                tier = 1
            elif impact == "LOW":
                tier = 2
            else:
                tier = 3
        info["tier"] = tier
        v["prefilter"] = info
        tiers[PREFILTER_TIERS[tier]] += 1
        kept.append((tier, order.get(v["class"], 99), n, v))
    for v in unqueried:
        v["prefilter"] = {"tier": 4, "status": f"beyond --max-prefilter {max_prefilter}: not frequency-filtered"}
        kept.append((4, order.get(v["class"], 99), len(kept), v))
    # heterozygous candidates per snpEff gene: a dominant-class het above the dominant cut-off can still be
    # half of a compound heterozygote, so it is only sent late when its gene has no second candidate
    per_gene: Counter = Counter()
    for _, _, _, v in kept:
        if v["class"] in DOMINANT_CLASSES:
            for g in (v["prefilter"].get("genes") or [])[:1]:
                per_gene[g] += 1

    def group(t: Tuple[int, int, int, Dict[str, Any]]) -> int:
        tier, _, _, v = t
        strong = v["class"] in STRONG_CLASSES
        if tier == 0:
            return 0
        if tier == 4:
            return 5 if strong else 6
        if strong:
            return 1
        af = v["prefilter"].get("filter_af")
        genes = (v["prefilter"].get("genes") or [])[:1]
        doomed = af is not None and af > dom_af and not any(per_gene[g] >= 2 for g in genes)
        if doomed:
            v["prefilter"]["late"] = (f"frequency {af:.3g} is above the dominant cut-off {dom_af:g} and no second "
                                      "candidate in the gene to pair with")
            return 4
        return 2 if tier == 1 else 3

    kept.sort(key=lambda t: (group(t),) + t[:3])
    survivors = [t[3] for t in kept]
    res = got.result if got is not None else {"queried": 0, "batches": 0, "answered_batches": 0, "found": 0,
                                               "not_found": 0, "unanswered": keys, "cached_batches": 0}
    version = None
    if meta is not None:
        version = meta.result["versions"]
    summary = {
        "queried": res["queried"], "batches": res["batches"], "answered_batches": res["answered_batches"],
        "cached_batches": res["cached_batches"], "found": res["found"], "not_in_myvariant": res["not_found"],
        "not_looked_up": len(unanswered) + len(unqueried), "dropped_af_gt_max": len(dropped),
        "kept_clinvar_exempt": exempt, "passed": len(survivors), "tiers": dict(tiers),
        "beyond_max_prefilter": len(unqueried), "source_versions": version,
        "answered_keys": res.get("answered_keys", res["queried"] - len(unanswered)),
        # MyVariant's metadata names gnomAD 2.1.1 for both; on hg38 the genome AFs carry v3's ami/mid groups
        "gnomad": (f"gnomAD v2.1.1 exomes (hg38 liftover) and v3 genomes via MyVariant.info"
                   if assembly == "GRCh38" else "gnomAD v2.1.1 exomes and genomes via MyVariant.info")
                  + (f" (metadata: gnomad {(version or {}).get('gnomad')})" if (version or {}).get("gnomad") else ""),
        "rule": f"drop when max(gnomAD exome, genome, grpmax afr/amr/eas/nfe/sas) > --max-af {max_af:g}, unless "
                f"ClinVar P/LP and ≤ {BA1_AF:g}; not found = kept; consequence only orders the VEP budget",
        "sent": f"{res['queried']} variant keys (HGVS g., no sample data) to myvariant.info",
    }
    if dropped:
        warnings.append(f"MyVariant prefilter: {len(dropped)} allele(s) removed as common ({summary['gnomad']}, "
                        f"> --max-af {max_af:g}), e.g. "
                        + ", ".join(f"{v['variant']} {af:.3g}" for v, af in dropped[:5]))
    if len(unanswered) or unqueried:
        warnings.append(f"MyVariant prefilter: {len(unanswered) + len(unqueried)} candidate(s) were NOT frequency-"
                        f"checked ({len(unanswered)} in failed/unanswered batches, {len(unqueried)} beyond "
                        f"--max-prefilter {max_prefilter}); they are kept, last in the VEP queue. Re-run the same "
                        "command to resume: answered batches come from the on-disk cache")
    return survivors, summary


SIB_RECESSIVE = ("hom_recessive", "hom", "x_hemizygous")


def _sibling_flags(cls: str, sibs: List[Dict[str, Any]], ctype: str = "auto") -> List[str]:
    """Segregation in the proband's full siblings: an affected sib should share the genotype, an unaffected one not.

    On X outside the PARs a sister and a brother are read differently: a brother's 1/1 is his one X (hemizygous),
    and a heterozygous sister is a carrier, who may or may not be affected (X-inactivation), so neither her
    being affected nor unaffected argues against an X-linked cause.
    """
    flags: List[str] = []
    on_x = ctype == "x_nonpar"
    for sb in sibs:
        z, name, sex = sb["z"], sb["id"], sb.get("sex")
        usable = z in CARRIER or sb.get("adequate_ref")
        if not usable:
            continue
        if on_x and sex == "female" and z == "het" and cls in ("x_hemizygous", "hom_recessive", "hom"):
            flags.append(f"{'affected' if sb['affected'] else 'unaffected' if sb['affected'] is False else ''} "
                         f"sister {name} is a heterozygous carrier (on X a carrier sister may or may not be affected)"
                         .replace("  ", " ").strip())
            continue
        if cls in SIB_RECESSIVE:
            same = z in ("hom_alt", "hemi")
            if sb["affected"] is True and not same:
                flags.append(f"affected sibling {name} is {z.replace('_', '-')} here: does not share the genotype "
                             "(against this variant as the cause, unless a phenocopy)")
            elif sb["affected"] is False and same:
                flags.append(f"unaffected sibling {name} has the same genotype ({z.replace('_', '-')}): against a "
                             "fully penetrant cause")
            elif sb["affected"] is True and same:
                flags.append(f"shared with affected sibling {name}")
        elif cls in DOMINANT_CLASSES:
            carries = z in CARRIER
            if sb["affected"] is True and not carries:
                flags.append(f"affected sibling {name} does not carry it (adequate hom-ref call)")
            elif sb["affected"] is True and carries:
                flags.append(f"shared with affected sibling {name}"
                             + (" — a de novo seen in two sibs points to parental germline mosaicism"
                                if cls in ("de_novo", "possible_de_novo") else ""))
            elif sb["affected"] is False and carries:
                if on_x and sex == "male" and z in ("hemi", "hom_alt"):
                    flags.append(f"unaffected brother {name} is hemizygous for it: against an X-linked cause with "
                                 "full penetrance in males")
                else:
                    flags.append(f"unaffected sibling {name} also carries it (reduced penetrance, or not causal)")
    return flags


def triage(path: str, proband: Optional[str] = None, mother: Optional[str] = None, father: Optional[str] = None,
           sex: Optional[str] = None, max_af: float = 0.01, min_dp: int = 10, min_gq: int = 20,
           assembly: Optional[str] = None, genes: Optional[Sequence[str]] = None, hpo_genes: bool = False,
           case_dir: Optional[str] = None, max_annotate: int = 1500, out: Optional[str] = None,
           hpo_terms: Optional[Sequence[str]] = None, max_af_dominant: Optional[float] = None,
           today: Optional[str] = None, ped: Optional[str] = None, prefilter: str = "none",
           max_prefilter: int = 50_000, s2f_top: int = 0, s2f_seconds: float = 180.0,
           progress: Optional[Any] = None) -> Outcome:
    """Filter, annotate and rank the proband's variants. See module docstring and the zebra-reanalysis skill."""
    import time

    t_start = time.monotonic()
    timings: "OrderedDict[str, float]" = OrderedDict()
    t_mark = [t_start]

    def lap(stage: str) -> None:
        now = time.monotonic()
        timings[stage] = round(now - t_mark[0], 2)
        t_mark[0] = now
        if progress:
            progress(f"{stage}: {timings[stage]} s")

    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []
    notes: List[str] = []
    if prefilter not in ("none", "myvariant"):
        raise UsageError("--prefilter must be none or myvariant")
    if max_prefilter < 1:
        raise UsageError("--max-prefilter must be at least 1")
    if s2f_top < 0:
        raise UsageError("--s2f-top must be >= 0")
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
    ped_info: Optional[Dict[str, Any]] = None
    sibs: List[Tuple[str, Optional[bool], Optional[str]]] = []
    if ped:
        proband, mother, father, sex, ped_info, sibs = _resolve_ped(ped, header.samples, proband, mother, father,
                                                                    sex, notes)
    if not proband:
        raise UsageError("--proband is required (or a --ped file naming one affected individual)")
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
    sex_check: Optional[Dict[str, Any]] = None
    if sex is not None:
        # B-P2-4: a given --sex is still checked against the genotypes; a disagreement
        # usually means a sample-label or PED error, and every X call rests on it
        sex_check = infer_sex(path, proband, assembly)
        sex_check["given"] = sex
        sex_check["agrees"] = None if sex_check["sex"] is None else sex_check["sex"] == sex
        if sex_check["agrees"] is False:
            warnings.append(f"--sex {sex} contradicts the proband's own genotypes, which look {sex_check['sex']} "
                            f"({sex_check['basis']}): check the sample label, the PED file and the sex before "
                            f"trusting any X call; X non-PAR calls were classed as {sex} because you said so")
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
    sv_types: Counter = Counter()
    sv_calls: List[Dict[str, Any]] = []
    star_alleles = 0
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
    # in-memory bound per inheritance class while streaming (the VEP budget is applied later)
    stream_cap = max_prefilter if prefilter == "myvariant" else max_annotate
    weak_stored = 0
    sib_idx = [(s, idx[s], aff, sx) for s, aff, sx in sibs]
    from zebra import qc as qc_mod

    _, records = iter_records(path)
    for rec in records:
        counts["records"] += 1
        if progress and counts["records"] % 200_000 == 0:
            progress(f"read {counts['records']:,} records (chr{rec.chrom})")
        for k, alt_raw in enumerate(rec.alts, 1):
            if is_structural(alt_raw):
                # B-P1-8: a CNV/SV allele is not an SNV/indel; it is named, never dropped silently
                sv_types[sv_type(alt_raw, rec.info())] += 1
                if len(sv_calls) < 50:
                    sv_calls.append(sv_call(rec, k, pi))
        if not rec.fmt or "GT" not in rec.fmt:
            no_gt += 1
            continue
        if chrom_lengths.get(rec.chrom) and rec.pos > chrom_lengths[rec.chrom]:
            off_contig += 1
        pcall = parse_call(rec.fmt, rec.sample_fields[pi]) if pi < len(rec.sample_fields) else None
        for k, alt_raw in enumerate(rec.alts, 1):
            if is_symbolic(alt_raw):
                counts["symbolic_skipped"] += 1
                if alt_raw == "*":
                    star_alleles += 1
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
            mosaic = None
            if pz == "het" and cls in ("de_novo", "possible_de_novo"):
                ref_reads, alt_reads = pcall.reads(k)
                mosaic = qc_mod.mosaic_assessment(ref_reads, alt_reads)
                if mosaic and mosaic["possible_mosaic"]:
                    flags.append(f"ALT fraction {mosaic['vaf']} (95% CI {mosaic['ci95'][0]}–{mosaic['ci95'][1]}, "
                                 f"{mosaic['depth']} reads) is below a germline het: possible postzygotic mosaic "
                                 "(or an artefact) — confirm in a second tissue")
            if pz == "het" and ab is not None and ab < 0.2 and not (mosaic and mosaic["possible_mosaic"]):
                flags.append(f"low allele balance {ab:.2f} (mosaic or artefact?)")
            by_class_seen[cls] += 1
            bucket = stored[cls]
            if len(bucket) >= stream_cap or (prefilter == "myvariant" and cls not in STRONG_CLASSES
                                             and weak_stored >= max_prefilter):
                # (with --prefilter, heterozygous classes beyond what MyVariant will be asked about are
                # only counted: they would never be looked up, and they are what fills memory on a genome)
                # budget: never more than this many per class are kept in
                # memory, with or without a gene set (an exome-wide gene list
                # over a WGS VCF would otherwise hold millions of dicts)
                over_budget[cls] += 1
                continue
            if cls not in STRONG_CLASSES:
                weak_stored += 1
            sib_calls = []
            for s_name, s_i, s_aff, s_sex in sib_idx:
                scall = parse_call(rec.fmt, rec.sample_fields[s_i]) if s_i < len(rec.sample_fields) else None
                sp = _parent(scall, k, min_dp, min_gq) or {"z": "missing", "adequate_ref": False, "why": []}
                sib_calls.append({"id": s_name, "affected": s_aff, "sex": s_sex, "z": sp["z"],
                                  "adequate_ref": sp["adequate_ref"],
                                  "gt": _gt_text(scall, k)})
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
                "sibs": sib_calls, "mosaic": mosaic,
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
    n_sv = sum(sv_types.values())
    if n_sv:
        counts["structural_alleles"] = n_sv
        type_txt = ", ".join(f"{t} {n}" for t, n in sv_types.most_common())
        examples = "; ".join(c["zebra_cnv"] for c in sv_calls if c.get("zebra_cnv"))[:600]
        if counts["alleles"] == 0:
            raise UsageError(
                f"this is a CNV/SV VCF: all {n_sv} ALT allele(s) are symbolic ({type_txt}) and SNV/indel triage "
                "cannot read them, so it would report 'no candidate' for a file it never assessed. Give each call "
                "to `zebra cnv`" + (f", e.g. {examples}" if examples else " (chrom:start-end loss|gain)"))
        warnings.append(f"{n_sv} symbolic CNV/SV allele(s) ({type_txt}) were NOT assessed: this triage reads SNVs "
                        "and indels only. Give each call to `zebra cnv` (listed in result.structural_variants"
                        + (f", e.g. {examples.split('; ')[0]}" if examples else "") + ")")
    if star_alleles:
        counts["spanning_deletion_alleles"] = star_alleles
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
    # B-P2-6 check runs before anything else leaves the machine (MyVariant, VEP)
    ref_keys = [(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in survivors_all]
    ref_check_result: Optional[Dict[str, Any]] = None
    if not build["guess"] and ref_keys:
        # B-P2-6: VEP's region endpoint annotates whatever REF it is given, so a VCF read
        # against the wrong build is annotated as if right. With no build in the header the
        # given --assembly is checked against the reference bases before anything is scored.
        refc = attempt(f"Ensembl reference check ({assembly})", lambda: ref_check(ref_keys, assembly), warnings)
        if refc is None:
            warnings.append(f"the build could not be verified (no build in the header and the reference check "
                            f"failed): annotations assume --assembly {assembly}")
        else:
            sources.extend(refc.sources)
            rc = refc.result
            ref_check_result = rc
            if rc["checked"] >= 5 and rc["mismatch"] / rc["checked"] >= 0.5:
                raise UsageError(
                    f"the VCF's REF bases do not match {assembly} at {rc['mismatch']} of {rc['checked']} sampled "
                    f"positions (e.g. {'; '.join(rc['examples'][:3])}): the file is not on {assembly}. Pass the other "
                    "--assembly (the header names no build), or check the file")
            if rc["mismatch"]:
                warnings.append(f"{rc['mismatch']} of {rc['checked']} sampled REF bases differ from {assembly} "
                                f"({'; '.join(rc['examples'][:3])}): check the build and the file's normalisation")
            if rc.get("asked") and rc["checked"] < rc["asked"] / 2:
                warnings.append(f"Ensembl returned no reference base for {rc['asked'] - rc['checked']} of "
                                f"{rc['asked']} sampled positions (past the end of a {assembly} chromosome, or an "
                                "unknown contig): the build is not verified")
    lap("read VCF, sex check, classify")
    prefilter_info: Optional[Dict[str, Any]] = None
    if prefilter == "myvariant":
        survivors_all, prefilter_info = _myvariant_prefilter(
            survivors_all, assembly, max_af, max_prefilter, (pheno or {}).get("genes") or {}, warnings, sources,
            progress, dom_af=dom_af)
        restriction["myvariant_prefilter"] = prefilter_info
        restriction["mode"] = (restriction.get("mode", "whole file") + " + MyVariant.info frequency prefilter")
        counts["myvariant_queried"] = prefilter_info["queried"]
        counts["myvariant_pass"] = len(survivors_all)
        lap("MyVariant prefilter")
    elif regions is None and not has_info_ann:
        restriction["mode"] = "prioritised by inheritance class (no gene set, no INFO annotations)"
        warnings.append("no --genes/--hpo-genes and no annotations in the VCF: annotating at most "
                        f"{max_annotate} variants, de novo / homozygous / hemizygous first. For a whole-exome "
                        "reanalysis pass --prefilter myvariant (sends every quality-passing allele's chrom-pos-ref-"
                        "alt, no sample data, to myvariant.info to drop common alleles first)")
    survivors = survivors_all
    if prefilter == "myvariant":
        kept = survivors[:max_annotate]
        cut = survivors[max_annotate:]
        survivors = kept
        counts["within_budget"] = len(kept)
        never_read = sum(over_budget.values())
        if cut:
            by_tier = Counter(PREFILTER_TIERS[v["prefilter"]["tier"]] for v in cut)
            by_cls = Counter(v["class"] for v in cut)
            strong_cut = {k: n for k, n in by_cls.items() if k in STRONG_CLASSES}
            restriction["skipped_over_budget"] = dict(by_cls)  # by class, as without --prefilter
            restriction["skipped_over_budget_by_tier"] = dict(by_tier)
            warnings.append(f"{len(cut)} rare candidate(s) that passed the MyVariant prefilter were NOT sent to VEP "
                            f"(over --max-annotate {max_annotate}; the most plausible were sent first): "
                            + ", ".join(f"{k} {n}" for k, n in by_tier.most_common())
                            + (" — INCLUDING " + ", ".join(f"{k} {n}" for k, n in strong_cut.items())
                               if strong_cut else "")
                            + "; compound-heterozygous pairs among them can be missed; raise --max-annotate to "
                              "annotate them too")
        if never_read:
            warnings.append(f"{never_read} candidate(s) beyond the in-memory cap of {stream_cap} per inheritance "
                            "class were never looked at: " + ", ".join(f"{k} {n}" for k, n in over_budget.items())
                            + " (a genome-scale VCF: restrict with --genes/--hpo-genes or use an annotated VCF)")
    elif len(survivors) > max_annotate or total_candidates > len(survivors_all):
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
    # VEP's colocated record often cannot be matched to the allele (GRCh37 serves no clin_sig_allele, and an
    # indel's rs record is written in another representation); MyVariant's ClinVar is keyed on the exact
    # allele, so when VEP gives none it is used, and said so
    mv_clinvar = 0
    cv_version = (((prefilter_info or {}).get("source_versions") or {}).get("clinvar")) or "version not read"
    for v in survivors:
        a, pf = v.get("ann"), v.get("prefilter") or {}
        if a is not None and not a.get("clinvar") and pf.get("clinvar"):
            a["clinvar"] = list(pf["clinvar"])
            a.setdefault("notes", []).append(f"ClinVar from MyVariant.info (allele-matched by HGVS g.; snapshot "
                                             f"{cv_version}); VEP's colocated record did not name this allele")
            mv_clinvar += 1
    if mv_clinvar:
        notes.append(f"{mv_clinvar} variant(s) carry ClinVar significance from MyVariant.info because VEP's "
                     "colocated record could not be matched to the allele")
    lap("VEP annotation")

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

    # ---- siblings (from the PED): segregation flags on the final class / pairing
    if sibs:
        by_variant = {v["variant"]: v for v in final}
        for v in final:
            if v.get("comphet") and v.get("partners"):
                for sb in v.get("sibs") or []:
                    if not (sb["z"] in CARRIER or sb.get("adequate_ref")):
                        continue
                    # the sibling's call at each partner, matched by sample id
                    at_partners = [next((x for x in (by_variant.get(p_) or {}).get("sibs") or []
                                         if x["id"] == sb["id"]), None) for p_ in v["partners"]]
                    known = [x for x in at_partners if x and (x["z"] in CARRIER or x.get("adequate_ref"))]
                    if sb["z"] not in CARRIER:
                        both = False
                    elif any(x["z"] in CARRIER for x in known):
                        both = True
                    elif len(known) < len(at_partners):
                        v["flags"].append(f"sibling {sb['id']}: genotype at the partner allele unknown (missing or "
                                          "below the depth/quality thresholds): segregation not assessed")
                        continue
                    else:
                        both = False
                    if sb["affected"] is False and both:
                        v["flags"].append(f"unaffected sibling {sb['id']} carries this allele and its partner "
                                          "(if in trans, against a fully penetrant recessive cause)")
                    elif sb["affected"] is True and not both:
                        v["flags"].append(f"affected sibling {sb['id']} does not carry both alleles of the pair")
                    elif sb["affected"] is True and both:
                        v["flags"].append(f"affected sibling {sb['id']} carries both alleles of the pair")
            else:
                v["flags"].extend(_sibling_flags(v["class"], v.get("sibs") or [], v.get("ctype", "auto")))

    # ---- scoring
    moi = gene_moi({v["gene"] for v in final if v.get("gene")})
    use_pheno = pheno is not None
    weights = dict(WEIGHTS) if use_pheno else {"variant": WEIGHTS["variant"] / (1 - WEIGHTS["phenotype"]),
                                               "inheritance": WEIGHTS["inheritance"] / (1 - WEIGHTS["phenotype"])}

    def score(v: Dict[str, Any]) -> None:
        a = v.get("ann")
        if a is None or a.get("impact") is None:
            v["consequence_score"] = v["predictor_score"] = v["variant_score"] = None
            v["support"] = []
            v["score"] = None
            return
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

    for v in final:
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
        score(v)
    order = {c: i for i, c in enumerate(CLASS_ORDER)}

    def rank() -> None:
        # the variant key breaks ties, so the ranking does not depend on the order candidates were annotated in
        final.sort(key=lambda v: (v["score"] is None, -(v["score"] or 0), order.get(v["class"], 99),
                                  -(v.get("variant_score") or 0), _chrom_sort_key(v["chrom"]), v["pos"], v["alt"]))

    rank()
    lap("frequency, compound heterozygotes, scoring")
    s2f_info: Optional[Dict[str, Any]] = None
    if s2f_top > 0 and final:
        s2f_info = _s2f_rerank(final, assembly, s2f_top, s2f_seconds, warnings, sources, progress)
        for v in final:
            if v.get("s2f"):
                score(v)
        rank()
        lap("S2F splice models (top candidates)")
    for v in final[:TOP_RESULT]:
        terms = set(((v.get("ann") or {}).get("consequence") or "").split(","))
        reg = sorted(terms & REGULATORY_TERMS)
        if reg:
            v["flags"].append(f"non-coding regulatory candidate ({', '.join(reg)}): splice models do not cover this — "
                              "run s2f_predict with --models alphagenome and the disease tissue (--ontology, the "
                              "UBERON/CL term of the affected tissue); zebra does not choose a tissue")
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
                   " absent = passes"
                   + (f"; before VEP, the same rule at --max-af on {prefilter_info['gnomad']}" if prefilter_info else ""),
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
    # what actually left the machine: keys in requests that were answered (and attempted, when different)
    def _n(done: int, tried: int) -> str:
        return f"{done}" + (f" (of {tried} attempted)" if tried != done else "")

    sent = []
    if ref_check_result:
        sent.append(f"{ref_check_result.get('checked', 0)} positions (no alleles) to Ensembl /sequence/region")
    if prefilter_info:
        sent.append(f"{_n(prefilter_info['answered_keys'], prefilter_info['queried'])} quality-passing "
                    "chrom-pos-ref-alt to MyVariant.info")
    if counts.get("sent_to_vep"):
        sent.append(f"{_n(counts.get('annotated', 0), counts['sent_to_vep'])} candidate chrom-pos-ref-alt to "
                    "Ensembl VEP")
    if s2f_info and s2f_info["asked"]:
        sent.append(f"{_n(len(s2f_info['ran']), len(s2f_info['asked']))} to the Broad SpliceAI-lookup")
    sources.insert(0, source_record("local VCF", os.path.basename(path), url=None,
                                    note="read on this machine; sent: " + ("; ".join(sent) or "nothing")
                                         + " (variant keys only: no sample names, genotypes or depths)"))
    lap("outputs")
    timings["total"] = round(time.monotonic() - t_start, 2)
    result = {
        "vcf": str(Path(os.path.expanduser(path)).resolve()),
        "proband": proband, "mother": mother, "father": father, "sex": sex, "mode": mode, "assembly": assembly,
        "sex_inferred": sex_inference if (sex_inference and sex_inference.get("sex")) else None,
        "sex_check": ({k: sex_check[k] for k in ("given", "sex", "agrees", "basis")} if sex_check else None),
        "ref_check": ref_check_result,
        "pedigree": ped_info,
        "sent_off_machine": sent,
        "timings": dict(timings),
        "s2f": s2f_info,
        "counts": dict(counts),
        "restriction": restriction,
        "thresholds": thresholds,
        "weights": {k: round(w, 3) for k, w in weights.items()},
        "scoring": ("score = Σ weight × component. variant = max(consequence severity HIGH 1/MODERATE 0.6/LOW 0.25/"
                    "MODIFIER 0.1, in-silico level ≤0.9); inheritance = class prior × gene-MOI fit (fits 1, unknown 0.8,"
                    " against 0.5); phenotype = gene's relative Resnik score for the case's HPO profile"),
        "phenotype": ({k: pheno[k] for k in ("origin", "present", "excluded", "hpo_version", "method", "notes")}
                      if pheno else None),
        # capped: the JSON envelope has a size budget and the TSV lists every pair
        "comphet_genes": comphet_genes[:20],
        "comphet_gene_count": len(comphet_genes),
        "candidates": [_row(v, i + 1) for i, v in enumerate(final[:TOP_RESULT])],
        "total_candidates": len(final),
        "tsv": str(Path(tsv_path).resolve()) if tsv_path else None,
        "notes": notes + ["SNV/indel VCF only: CNVs, repeat expansions, mtDNA (unless called), low-level mosaicism "
                          "and poorly covered regions are not assessed"
                          + (f" ({n_sv} symbolic CNV/SV allele(s) in this file were skipped: see "
                             "structural_variants)" if n_sv else "")],
    }
    if n_sv:
        result["structural_variants"] = {"count": n_sv, "types": dict(sv_types), "calls": sv_calls[:20],
                                         "listed": min(20, len(sv_calls)),
                                         "next": "zebra cnv \"chrN:start-end loss|gain\" for each call"}
    text = _render(result, final)
    query = {"vcf": path, "proband": proband, "mother": mother, "father": father, "sex": sex, "max_af": max_af,
             "max_af_dominant": dom_af, "min_dp": min_dp, "min_gq": min_gq, "assembly": assembly,
             "genes": len(genes or []), "hpo_genes": hpo_genes, "max_annotate": max_annotate, "out": tsv_path,
             "ped": ped, "prefilter": prefilter, "max_prefilter": max_prefilter, "s2f_top": s2f_top}
    return Outcome(result, sources=sources, warnings=warnings, text=text, query=query)


def _row(v: Dict[str, Any], rank: int) -> Dict[str, Any]:
    a = v.get("ann") or {}
    af = a.get("af") or {}
    hgvsc = a.get("hgvsc")
    row = {
        "rank": rank, "score": v.get("score"),
        "components": {"phenotype": v.get("phenotype_score"), "variant": v.get("variant_score"),
                       "consequence": v.get("consequence_score"), "predictors": v.get("predictor_score"),
                       "inheritance": v.get("inheritance_score"),
                       "splice_s2f": a.get("spliceai_lookup")},
        "variant": v["variant"], "gene": v.get("gene"), "transcript": a.get("transcript"), "mane": a.get("mane"),
        "hgvsc": hgvsc.split(":", 1)[1] if hgvsc and ":" in hgvsc else hgvsc,
        "hgvsp": (a.get("hgvsp") or "").split(":", 1)[-1].replace("%3D", "=") or None,
        "consequence": a.get("consequence"), "impact": a.get("impact"),
        "class": v["class"], "origin": v.get("origin"), "moi": v.get("moi"), "moi_fit": v.get("moi_fit"),
        "partners": v.get("partners"),
        "gnomad": {k: x for k, x in (
            ("filter_af", af.get("filter_af")), ("exome", af.get("exome")), ("genome", af.get("genome")),
            ("grpmax", af.get("grpmax")), ("grpmax_group", af.get("grpmax_group")), ("source", af.get("source")),
            ("bound_af", af.get("bound_af")), ("bound_basis", af.get("bound_basis")), ("faf95", af.get("faf95")),
            ("faf95_group", af.get("faf95_group"))) if x is not None or k == "filter_af"},
        "af_exempt": v.get("af_exempt"), "comphet": v.get("comphet"),
        "phase": {k: x for k, x in (v.get("phase") or {}).items() if x not in (None, False)} or None,
        "predictors": {k: x for k, x in (("revel", a.get("revel")), ("alphamissense", a.get("alphamissense")),
                                         ("am_class", a.get("am_class")), ("cadd", a.get("cadd")),
                                         ("spliceai_max", a.get("spliceai_max"))) if x is not None},
        "support": v.get("support"), "clinvar": a.get("clinvar") or None, "rsid": a.get("rsid") or v.get("vcf_id"),
        "phenotype_via": v.get("phenotype_via"), "genotypes": {k: g for k, g in v["gt"].items() if g},
        "flags": v["flags"] + (a.get("notes") or []), "annotation": a.get("annotation"),
    }
    # the v0.2 fields are present only when they say something (the envelope has a size budget)
    extra = {
        "s2f": v.get("s2f"),
        "siblings": ({sb["id"]: f"{sb['gt'] or '.'} ({'affected' if sb['affected'] else 'unaffected' if sb['affected'] is False else 'status unknown'})"
                      for sb in v.get("sibs") or []} or None),
        "mosaic": v.get("mosaic") if (v.get("mosaic") or {}).get("possible_mosaic") else None,
        "prefilter": ({k: x for k, x in (v.get("prefilter") or {}).items() if k in ("tier", "status", "filter_af",
                                                                                  "kept_because", "late")} or None),
    }
    row.update({k: x for k, x in extra.items() if x is not None})
    return row


TSV_COLUMNS = ("rank", "score", "phenotype_score", "variant_score", "consequence_score", "predictor_score",
               "inheritance_score", "variant", "gene", "class", "origin", "moi", "moi_fit", "consequence", "impact",
               "transcript", "mane", "hgvsc", "hgvsp", "gnomad_filter_af", "gnomad_exome_af", "gnomad_genome_af", "grpmax_af",
               "grpmax_group", "revel", "alphamissense", "am_class", "cadd", "spliceai_max", "support", "clinvar",
               "rsid", "proband_gt", "mother_gt", "father_gt", "allele_balance", "partners", "phenotype_via", "flags",
               "annotation", "gnomad_bound_af", "gnomad_bound_basis", "gnomad_faf95", "af_exempt", "phase_set",
               "spliceai_lookup", "pangolin_lookup", "siblings", "prefilter_tier")


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
            g, pr = defaultdict(lambda: None, r["gnomad"]), defaultdict(lambda: None, r["predictors"])
            vals = [r["rank"], r["score"], r["components"]["phenotype"], r["components"]["variant"],
                    r["components"]["consequence"], r["components"]["predictors"], r["components"]["inheritance"],
                    r["variant"], r["gene"], r["class"], r["origin"], r["moi"], r["moi_fit"], r["consequence"],
                    r["impact"], r["transcript"], r["mane"], r["hgvsc"], r["hgvsp"], g["filter_af"], g["exome"], g["genome"],
                    g["grpmax"], g["grpmax_group"], pr["revel"], pr["alphamissense"], pr["am_class"], pr["cadd"],
                    pr["spliceai_max"], r["support"], r["clinvar"], r["rsid"], v["gt"].get("proband"),
                    v["gt"].get("mother"), v["gt"].get("father"), v.get("ab"), r["partners"], r["phenotype_via"],
                    r["flags"], r["annotation"], g["bound_af"], g["bound_basis"], g["faf95"], r["af_exempt"],
                    (v.get("phase") or {}).get("ps"),
                    ((v.get("s2f") or {}).get("spliceai") or {}).get("value"),
                    ((v.get("s2f") or {}).get("pangolin") or {}).get("value"),
                    [f"{k}={x}" for k, x in (r.get("siblings") or {}).items()],
                    (v.get("prefilter") or {}).get("tier")]
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
    if "myvariant_pass" in c:
        pf = result["restriction"].get("myvariant_prefilter") or {}
        steps.append(f"{c['myvariant_queried']} looked up in MyVariant.info → {c['myvariant_pass']} rare or unknown"
                     + (f" ({pf['gnomad'].split(' via ')[0]})" if pf.get("gnomad") else ""))
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
                                                                       ("inheritance", "inh"), ("splice_s2f", "s2f"))
                              if comp.get(k) is not None)
        score = f"{r['score']:.3f}" if r["score"] is not None else "  -  "
        lines.append(f"{r['rank']:>2}. {score} {r['gene'] or '-':<9} {r['variant']:<22} {hg[:48]:<48} "
                     f"{(r['consequence'] or '?').split(',')[0]} {r['impact'] or '?'}  {cls}  gnomAD {af_txt}  "
                     f"MOI {','.join(r['moi'] or []) or '?'} ({r['moi_fit']})  [{comp_txt}]")
        extra = list(r["support"] or [])
        if r.get("s2f"):
            sp, pg = r["s2f"].get("spliceai") or {}, r["s2f"].get("pangolin") or {}
            extra.append("s2f: " + "; ".join(
                f"{m} {x['score']} {x['value']:+.2f}" if x.get("status") == "ran" and isinstance(x.get("value"), (int, float))
                else f"{m} {x.get('status')}" for m, x in (("SpliceAI", sp), ("Pangolin", pg)) if x))
        if r.get("siblings"):
            extra.append("sibs " + ", ".join(f"{k} {x}" for k, x in r["siblings"].items()))
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
    sv = result.get("structural_variants")
    if sv:
        lines.append(f"CNV/SV: {sv['count']} symbolic allele(s) ({', '.join(f'{k} {v}' for k, v in sv['types'].items())}) "
                     "NOT assessed by this triage — give each to `zebra cnv`:")
        for c in sv["calls"][:5]:
            lines.append(f"      {c['chrom']}:{c['pos']}-{c.get('end') or '?'} {c['svtype']} GT {c.get('gt') or '?'}"
                         + (f"  → {c['zebra_cnv']}" if c.get("zebra_cnv") else ""))
    lines.append(f"TSV: {result['tsv'] or '(not written)'}")
    if result.get("timings"):
        lines.append("wall time: " + ", ".join(f"{k} {x:g} s" for k, x in result["timings"].items()))
    lines.extend(f"note: {n}" for n in result["notes"])
    return "\n".join(lines)
