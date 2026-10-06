"""Chinese-population allele frequencies for one variant, from the resources that can be queried by script.

Queried (each verified live on 2026-10-06 with rs671 12-111803962-G-A GRCh38):

* NyuWa NCVD -- 2,999 Chinese high-depth WGS (23 of 34 provincial divisions), GRCh38 only.
  `GET http://bigdata.ibp.ac.cn/NyuWa_variants/singlepage.php?variant=C-P-R-A` (no key; the
  database's own left-normalised chrom-pos-ref-alt id). The answer is an HTML page: the
  numbers sit in `<th>Allele Count: </th><th>1323</th>` rows and a variant it does not hold
  answers "There is no results for your searching request." with HTTP 200. robots.txt
  disallows `/*/searchact*` and `/*/browser*`, not `singlepage.php`. https does not answer
  (TLS reset), so this is plain http. GRCh37 input is lifted with Ensembl `/map` and checked
  against the GRCh38 reference before it is sent.
* WBBC (Westlake BioBank for Chinese) -- WGS of 4,480 Chinese individuals (the genotype counts
  it returns sum to 4,480), GRCh38 and GRCh37 both native, autosomes only (its search accepts
  numeric chromosomes only and its downloads are chr1-22).
  `POST https://wbbc.westlake.edu.cn/search_position.php` with JSON
  `{"position": "12:111803962", "grch": "grch38"}`: the endpoint the site's search page calls.
  It returns every allele at the position as a JSON list, and an EMPTY body when it holds
  none. Allele counts come from its `Genotype` field (`RR=2618|RA=1586|AA=276`) and are
  cross-checked against its `WBBC_AF`.
* 1000 Genomes phase 3 CHB / CHS / CDX (103 / 105 / 93 individuals) through Ensembl REST:
  `/overlap/region` finds the dbSNP record that is this allele, `/variation/<rs>?pops=1`
  gives allele counts per population. Small samples; labelled so.
* Taiwan Biobank WGS (TaiwanView) -- 1,492 individuals, GRCh38. `POST
  https://taiwanview.twbiobank.org.tw/ssss_query.php` (form, the page's own call). Searchable
  by rsID only: its region search answers HTTP 500 (tested 2026-10-06 with the site's own
  example region), so it is queried with the rsID the 1000 Genomes step found. It publishes
  AF rounded to 4 decimals and the number of samples (NS), not allele counts.

Not queried -- see NOT_QUERIED for the reason each is named in every answer: ChinaMAP, CMDB,
NyuWa Global 10K T2T, PGG.Han 2.0, PGG.SNV. gnomAD's East Asian group is not Chinese-specific
and is not counted here.

A dataset that answers but holds no record of the allele goes to `not_found_in` with a note,
never into `datasets` as af=0: these resources list observed alleles only, so a missing
record is not an observation of absence. The one exception is a 1000 Genomes population whose
genotypes at this dbSNP site Ensembl does report: there AC=0 out of AN called alleles is a
real count and is reported as such.

Before any dataset is asked, REF is checked against the Ensembl reference of the requested
build and the allele is left-aligned (the resources above store left-aligned alleles); a REF
that is not the reference is refused with `UsageError`, naming the other build when it matches
there.
"""

from __future__ import annotations

import json
import re
import urllib.parse
from typing import Any, Callable, Dict, List, Optional, Tuple

from zebra.core import Outcome, UsageError
from zebra.http import Response, SourceError, request
from zebra.sources import attempt, record, validated_json, validated_text

ENSEMBL = {"GRCh38": "https://rest.ensembl.org", "GRCh37": "https://grch37.rest.ensembl.org"}
OTHER = {"GRCh38": "GRCh37", "GRCh37": "GRCh38"}
NYUWA_URL = "http://bigdata.ibp.ac.cn/NyuWa_variants/singlepage.php"
WBBC_URL = "https://wbbc.westlake.edu.cn/search_position.php"
TWB_API = "https://taiwanview.twbiobank.org.tw/ssss_query.php"
TWB_PAGE = "https://taiwanview.twbiobank.org.tw/variantT.php"

STATIC_TTL = 30 * 86400  # the cohorts are frozen releases (NyuWa 2021, WBBC v20210103, TWB WGS)
ENSEMBL_TTL = 14 * 86400
SEQUENCE_TTL = 90 * 86400

NYUWA = "NyuWa"
WBBC = "WBBC"
TWB = "Taiwan Biobank WGS"
KG_POPS = (
    ("CHB", "Han Chinese in Beijing, China", 103),
    ("CHS", "Southern Han Chinese, China", 105),
    ("CDX", "Chinese Dai in Xishuangbanna, China", 93),
)
SAMPLES = {NYUWA: 2999, WBBC: 4480, TWB: 1492}
POPULATION = {
    NYUWA: "Chinese (NyuWa: 2,999 high-depth WGS individuals from 23 of 34 provincial divisions)",
    WBBC: "Chinese (Westlake BioBank for Chinese: WGS of 4,480 individuals)",
    TWB: "Taiwan (Taiwan Biobank WGS participants, predominantly Han Chinese ancestry; 1,492 individuals)",
}

# Resources tested and not queried: (name, short reason, detail). Named in every answer
# (one warning listing them all + one not_checked entry each).
NOT_QUERIED = (
    ("ChinaMAP", "login with captcha required",
     "every search and download redirects to a login form with a captcha (chinamapwgs.mbiobank.com "
     "/search/ -> /login/); no public API; its terms require an institutional account and forbid linking "
     "the data to third-party databases"),
    ("CMDB", "API token required",
     "its API (cmdb.bgi.com/api/v1.0/variant) needs a per-user token and answers 403 'token error' without "
     "one; variant pages redirect to login; GRCh37 only; terms forbid linking the data to third-party databases"),
    ("NyuWa Global 10K T2T", "T2T-CHM13 coordinates only",
     "T2T-CHM13 coordinates only and no verified GRCh38/GRCh37 to CHM13 liftover service; its Chinese subset "
     "is shown as a rounded AF without counts"),
    ("PGG.Han 2.0", "deferred until zebra.http paces www.biosino.org (robots.txt Crawl-delay: 10)",
     "scriptable (search page -> variant id -> JSON province table, 3 requests per variant) but "
     "www.biosino.org robots.txt sets Crawl-delay: 10 and zebra.http has no pacing for that host yet"),
    ("PGG.SNV", "GRCh37 only, ~800 KB per variant, Chinese groups mostly duplicate 1000 Genomes",
     "GRCh37 only, about 800 KB of JSON per variant; its Chinese groups are 1000 Genomes and HapMap (queried "
     "here through Ensembl) plus small AAGC/SGDP samples; www.pggsnv.org no longer serves the database "
     "(mirror pog.fudan.edu.cn)"),
)
NOTES = ["gnomAD's East Asian group (eas) is not Chinese-specific and is not counted here."]

NYUWA_NO_RESULT = "There is no results for your searching request."
CHROMS = frozenset([str(i) for i in range(1, 23)] + ["X", "Y", "MT"])
ALLELE_RE = re.compile(r"^[ACGT]+$")
SEQ_RE = re.compile(r"^[ACGTN]+$")
MAX_ALLELE = 1000  # longer alleles are structural variants: out of scope here, answered with named gaps
FLANK = 200  # reference bases fetched on each side for the REF check and left-alignment


class _OutOfWindow(Exception):
    """A reference base outside the fetched window was needed."""


class _Window:
    """Forward-strand reference chrom:start.., 1-based."""

    def __init__(self, chrom: str, start: int, seq: str):
        self.chrom, self.start, self.seq = chrom, start, seq.upper()
        self.end = start + len(seq) - 1

    def base(self, p: int) -> str:
        if p < self.start or p > self.end:
            raise _OutOfWindow(p)
        return self.seq[p - self.start]

    def slice(self, a: int, b: int) -> str:
        if a < self.start or b > self.end:
            raise _OutOfWindow((a, b))
        return self.seq[a - self.start:b - self.start + 1]


Var = Tuple[str, int, str, str]


def vid(v: Var) -> str:
    return f"{v[0]}-{v[1]}-{v[2]}-{v[3]}"


# ---------------------------------------------------------------- input

def _chrom(chrom: Any) -> str:
    if isinstance(chrom, bool) or not isinstance(chrom, (str, int)):
        raise UsageError(f"chromosome must be 1-22, X, Y or MT (got {chrom!r})")
    s = str(chrom).strip()
    if s[:3].lower() == "chr":
        s = s[3:]
    s = s.upper()
    if s == "M":
        s = "MT"
    if s not in CHROMS:
        raise UsageError(f"chromosome must be 1-22, X, Y or MT, with or without 'chr' (got {chrom!r})")
    return s


def _pos(pos: Any) -> int:
    if isinstance(pos, bool):
        raise UsageError(f"position must be a positive integer (got {pos!r})")
    if isinstance(pos, str) and pos.strip().isdigit():
        pos = int(pos.strip())
    if not isinstance(pos, int) or pos < 1 or pos > 250_000_000:
        raise UsageError(f"position must be a positive integer within a chromosome (got {pos!r})")
    return pos


def _allele(a: Any, what: str) -> str:
    if not isinstance(a, str):
        raise UsageError(f"{what} must be a string of A/C/G/T (got {a!r})")
    s = a.strip().upper()
    if not ALLELE_RE.match(s):
        raise UsageError(f"{what} must be A/C/G/T only (got {a!r}); symbolic, '-', '.', '*' and N alleles "
                         "cannot be looked up")
    return s


def _assembly(assembly: Any) -> str:
    for name in ENSEMBL:
        if isinstance(assembly, str) and assembly.strip().lower() == name.lower():
            return name
    raise UsageError(f"assembly must be GRCh38 or GRCh37 (got {assembly!r})")


def parse_input(chrom: Any, pos: Any, ref: Any, alt: Any, assembly: Any) -> Tuple[Var, str]:
    c, p, r, a, asm = _chrom(chrom), _pos(pos), _allele(ref, "ref"), _allele(alt, "alt"), _assembly(assembly)
    if r == a:
        raise UsageError(f"ref and alt are the same ({r})")
    return (c, p, r, a), asm


# ---------------------------------------------------------------- normalisation

def trim(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """Parsimony without reference sequence: drop shared trailing, then shared leading bases (one kept)."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def left_normalize(pos: int, ref: str, alt: str, base_at: Callable[[int], str]) -> Tuple[int, str, str]:
    """Parsimonious, left-aligned VCF representation (Tan et al. 2015)."""
    for _ in range(100000):
        changed = False
        if not ref or not alt:
            b = base_at(pos - 1)
            ref, alt, pos = b + ref, b + alt, pos - 1
            changed = True
        if ref and alt and ref[-1] == alt[-1] and (len(ref) > 1 or len(alt) > 1):
            ref, alt = ref[:-1], alt[:-1]
            changed = True
        if not changed:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def right_end(pos: int, ref: str, alt: str, base_at: Callable[[int], str]) -> int:
    """Last reference base a left-aligned indel can slide over (its 3'-shifted end); SNV/MNV: its own end."""
    if len(ref) > len(alt) and len(alt) == 1 and ref[0] == alt:
        d, p = ref[1:], pos
        while base_at(p + 1 + len(d)) == d[0]:
            d, p = d[1:] + d[0], p + 1
        return p + len(d)
    if len(alt) > len(ref) and len(ref) == 1 and alt[0] == ref:
        i, p = alt[1:], pos
        while base_at(p + 1) == i[0]:
            i, p = i[1:] + i[0], p + 1
        return p + 1
    return pos + len(ref) - 1


# ---------------------------------------------------------------- Ensembl helpers

def _seq_check(text: str) -> Optional[str]:
    t = (text or "").strip().upper()
    return None if t and SEQ_RE.match(t) else "not a reference sequence"


def _sequence(chrom: str, start: int, end: int, assembly: str) -> Tuple[str, Response]:
    region = f"{chrom}:{start}..{end}:1"
    resp = request(f"{ENSEMBL[assembly]}/sequence/region/human/{region}", source="Ensembl sequence",
                   accept="text/plain", cache_ttl=SEQUENCE_TTL, validate=_seq_check)
    seq = validated_text(resp, "Ensembl sequence").strip().upper()
    if not SEQ_RE.match(seq):
        raise SourceError("Ensembl sequence", resp.url, resp.status, "body is not a DNA sequence")
    return seq, resp


class _RefDiffers(Exception):
    def __init__(self, found: str):
        super().__init__(found)
        self.found = found


def _frame(v: Var, assembly: str, sources: List[Dict[str, Any]]) -> Tuple[Var, _Window, bool]:
    """REF checked against `assembly` and the allele left-aligned (False: the window was too short to finish).

    Raises _RefDiffers, SourceError, UsageError (position past the chromosome end).
    """
    chrom, pos, ref, alt = v
    start = max(1, pos - FLANK)
    try:
        seq, resp = _sequence(chrom, start, pos + len(ref) - 1 + FLANK, assembly)
    except SourceError as err:
        if err.status == 400 and "greater than" in (err.message or ""):
            raise UsageError(f"{chrom}:{pos} is beyond the end of chromosome {chrom} in {assembly}") from None
        raise
    sources.append(record("Ensembl sequence", f"{assembly} {chrom}:{start}-{start + len(seq) - 1}", resp,
                          note="REF check and left-alignment"))
    w = _Window(chrom, start, seq)
    try:
        found = w.slice(pos, pos + len(ref) - 1)
    except _OutOfWindow:
        raise UsageError(f"{chrom}:{pos}-{pos + len(ref) - 1} runs past the end of chromosome {chrom} "
                         f"in {assembly}") from None
    if found != ref:
        raise _RefDiffers(found)
    try:
        p, r, a = left_normalize(pos, ref, alt, w.base)
        aligned = True
    except _OutOfWindow:
        p, r, a = trim(pos, ref, alt)
        aligned = False
    return (chrom, p, r, a), w, aligned


def _other_build_bases(v: Var, assembly: str) -> Optional[str]:
    try:
        seq, _ = _sequence(v[0], v[1], v[1] + len(v[2]) - 1, OTHER[assembly])
        return seq
    except (SourceError, ValueError):
        return None


def _lift_37_to_38(v: Var, sources: List[Dict[str, Any]]) -> int:
    """GRCh38 position of a GRCh37 allele, through Ensembl /map; raises ValueError when not one clean block."""
    chrom, pos, ref, _ = v
    end = pos + len(ref) - 1
    region = f"{chrom}:{pos}..{end}:1"
    resp = request(f"{ENSEMBL['GRCh38']}/map/human/GRCh37/{region}/GRCh38", source="Ensembl assembly map",
                   cache_ttl=SEQUENCE_TTL)
    data = validated_json(resp, "Ensembl assembly map", require="mappings")
    sources.append(record("Ensembl assembly map", f"GRCh37 {region} -> GRCh38", resp))
    maps = data["mappings"]
    if len(maps) != 1:
        raise ValueError(f"GRCh37 {region} maps to {len(maps)} GRCh38 blocks; not lifted")
    o, m = maps[0]["original"], maps[0]["mapped"]
    if (str(m.get("seq_region_name")) != chrom or m.get("strand") != 1 or o.get("start") != pos
            or o.get("end") != end or m["end"] - m["start"] != end - pos):
        raise ValueError(f"GRCh37 {region} maps to {m.get('seq_region_name')}:{m.get('start')}-{m.get('end')} "
                         f"strand {m.get('strand')}; not a clean same-length block, not lifted")
    return int(m["start"])


# ---------------------------------------------------------------- result helpers

def _found(row: Dict[str, Any]) -> Dict[str, Any]:
    return {"status": "found", "name": row["name"], "row": row}


def _absent(name: str, note: str) -> Dict[str, Any]:
    return {"status": "absent", "name": name, "note": note}


def _absent_note(name: str, covers: str, extra: str = "") -> str:
    return (f"no record in {name}{extra}; {name} covers {covers} — not evidence the allele is absent")


# ---------------------------------------------------------------- NyuWa

def _nyuwa_check(text: str) -> Optional[str]:
    if NYUWA_NO_RESULT in (text or ""):
        return None
    if "Allele Count:" in text and "Variant:" in text:
        return None
    return "neither a NyuWa variant record nor its no-result message (login, error or maintenance page)"


def _th(text: str, label: str) -> str:
    m = re.search(rf"<th>\s*{re.escape(label)}:\s*</th>\s*<th>\s*([^<]*?)\s*</th>", text)
    if not m:
        raise ValueError(f"NyuWa page has no '{label}' row")
    return m.group(1)


def parse_nyuwa(text: str, variant_id: str) -> Optional[Dict[str, Any]]:
    """The NyuWa record on a singlepage.php answer, or None when it says it has no such variant."""
    if NYUWA_NO_RESULT in text:
        return None
    m = re.search(r"Variant:\s*([0-9XYMT]+-\d+-[ACGT]+-[ACGT]+)\s*\(([^)]*)\)", text)
    if not m:
        raise ValueError("NyuWa page has no 'Variant:' header")
    if m.group(1) != variant_id:
        raise ValueError(f"NyuWa answered for {m.group(1)}, not the {variant_id} that was asked")
    if m.group(2) != "GRCh38":
        raise ValueError(f"NyuWa page says assembly {m.group(2)!r}, expected GRCh38")
    ac, an = int(_th(text, "Allele Count")), int(_th(text, "Allele Number"))
    af_page = float(_th(text, "Allele Frequency"))
    hom = _th(text, "Number of Homozygotes")
    if an <= 0 or not 0 <= ac <= an:
        raise ValueError(f"NyuWa counts are inconsistent (AC {ac}, AN {an})")
    if abs(ac / an - af_page) > max(5e-5, 0.001 * af_page):
        raise ValueError(f"NyuWa AF {af_page} does not match AC/AN {ac}/{an}")
    rs = re.search(r"snp/\?term=(rs\d+)", text)
    return {"ac": ac, "an": an, "af_page": af_page, "homozygotes": int(hom) if hom.isdigit() else None,
            "rsid": rs.group(1) if rs else None}


def _nyuwa(v38: Var, lifted_from: Optional[str], sources: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    variant_id = vid(v38)
    url = f"{NYUWA_URL}?{urllib.parse.urlencode({'variant': variant_id})}"
    resp = request(url, source="NyuWa NCVD", accept="text/html", cache_ttl=STATIC_TTL, validate=_nyuwa_check)
    rec = parse_nyuwa(resp.text, variant_id)
    sources.append(record("NyuWa NCVD", variant_id, resp,
                          note="GRCh38; " + ("record found" if rec else "no record") +
                               (f"; lifted from {lifted_from}" if lifted_from else "")))
    if rec is None:
        return [_absent(NYUWA, _absent_note(NYUWA, "2,999 Chinese WGS individuals (GRCh38)",
                                            f" for {variant_id}"))]
    row = {"name": NYUWA, "population": POPULATION[NYUWA], "af": rec["ac"] / rec["an"], "ac": rec["ac"],
           "an": rec["an"], "url": url, "samples": SAMPLES[NYUWA], "homozygotes": rec["homozygotes"],
           "assembly": "GRCh38", "queried_as": variant_id, "rsid": rec["rsid"]}
    if lifted_from:
        row["lifted_from"] = lifted_from
        row["note"] = f"queried as GRCh38 {variant_id}, lifted from {lifted_from} through Ensembl /map " \
                      "and checked against the GRCh38 reference"
    return [_found(row)]


# ---------------------------------------------------------------- WBBC

def _json_list_or_empty(source: str) -> Callable[[str], Optional[str]]:
    def check(text: str) -> Optional[str]:
        t = (text or "").strip()
        if not t:
            return None  # this source's "no record" answer
        if t[:1] != "[":
            return f"{source} answered {t[:60]!r} where a JSON list was expected"
        try:
            json.loads(t)
        except ValueError as err:
            return f"{source} body is not JSON ({err})"
        return None
    return check


GENOTYPE_RE = re.compile(r"^RR=(\d+)\|RA=(\d+)\|AA=(\d+)$")


def _wbbc(v: Var, assembly: str, sources: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    chrom, pos, ref, alt = v
    body = {"position": f"{chrom}:{pos}", "grch": assembly.lower()}
    resp = request(WBBC_URL, source="WBBC", method="POST", body=body, accept="application/json",
                   cache_ttl=STATIC_TTL, validate=_json_list_or_empty("WBBC"))
    text = resp.text.strip()
    rows = json.loads(text) if text else []
    if not isinstance(rows, list) or not all(isinstance(r, dict) for r in rows):
        raise ValueError("WBBC answer is not a list of records")
    sources.append(record("WBBC", f"{assembly} {chrom}:{pos}", resp,
                          note=f"POST position={chrom}:{pos} grch={assembly.lower()}; {len(rows)} record(s) at the position"))
    hits = [r for r in rows if str(r.get("Chr")) == chrom and str(r.get("Position")) == str(pos)
            and str(r.get("Ref")).upper() == ref and str(r.get("Alt")).upper() == alt]
    covers = f"WGS of 4,480 Chinese individuals ({assembly}, autosomes)"
    if not hits:
        others = sorted({f"{r.get('Ref')}>{r.get('Alt')}" for r in rows})
        extra = f" for {ref}>{alt} at {chrom}:{pos}" + (f" (it holds {', '.join(others)} there)" if others else "")
        return [_absent(WBBC, _absent_note(WBBC, covers, extra))]
    r = hits[0]
    m = GENOTYPE_RE.match(str(r.get("Genotype") or ""))
    af_site = float(r["WBBC_AF"]) if r.get("WBBC_AF") not in (None, "") else None
    row: Dict[str, Any] = {"name": WBBC, "population": POPULATION[WBBC], "url": WBBC_URL,
                           "query": body, "assembly": assembly, "queried_as": vid(v), "rsid": r.get("ID") or None}
    if m:
        rr, ra, aa = (int(x) for x in m.groups())
        n = rr + ra + aa
        ac, an = ra + 2 * aa, 2 * n
        if an and (af_site is None or abs(ac / an - af_site) <= 1e-6):
            row.update(af=ac / an, ac=ac, an=an, samples=n, homozygotes=aa)
        elif af_site is None:
            raise ValueError(f"WBBC record has no usable genotype counts or AF ({r.get('Genotype')!r})")
        else:
            row.update(af=af_site, ac=None, an=None, samples=n,
                       note=f"genotype counts {r.get('Genotype')} do not reproduce WBBC_AF {af_site}; "
                            "WBBC_AF reported, counts withheld")
    elif af_site is not None:
        row.update(af=af_site, ac=None, an=None, samples=None,
                   note=f"no genotype counts in the record ({r.get('Genotype')!r}); WBBC_AF reported")
    else:
        raise ValueError(f"WBBC record has no usable genotype counts or AF ({r.get('Genotype')!r})")
    subs = []
    for key, label in (("North_AF", "North"), ("Central_AF", "Central"), ("South_AF", "South"),
                       ("Lingnan_AF", "Lingnan")):
        if r.get(key) not in (None, ""):
            subs.append({"population": label, "af": float(r[key])})
    if subs:
        row["subpopulations"] = subs
    if len(hits) > 1:
        row["note"] = (row.get("note", "") + f"; WBBC returned {len(hits)} records for this allele, "
                                              "the first is reported").lstrip("; ")
    return [_found(row)]


# ---------------------------------------------------------------- 1000 Genomes via Ensembl

def _ensembl_vcf(entry: Dict[str, Any], alt_str: str, w: Optional[_Window]) -> Optional[Tuple[int, str, str]]:
    """One Ensembl variation allele as a left-aligned VCF allele, or None when it cannot be read."""
    alleles = entry.get("alleles") or []
    ref_str = str(alleles[0]) if alleles else ""
    start, end = int(entry["start"]), int(entry["end"])
    if alt_str == ref_str or not all(ALLELE_RE.match(s) or s == "-" for s in (ref_str, alt_str)):
        return None
    if ref_str == "-" and alt_str == "-":
        return None
    if w is None:
        if "-" in (ref_str, alt_str) or len(ref_str) != len(alt_str):
            return None
        return start, ref_str, alt_str
    try:
        if ref_str == "-":  # insertion between `end` and `start` (start = end + 1)
            b = w.base(end)
            pos, ref, alt = end, b, b + alt_str
        elif alt_str == "-":
            if w.slice(start, end) != ref_str:
                return None
            b = w.base(start - 1)
            pos, ref, alt = start - 1, b + ref_str, b
        else:
            if w.slice(start, start + len(ref_str) - 1) != ref_str:
                return None
            pos, ref, alt = start, ref_str, alt_str
        return left_normalize(pos, ref, alt, w.base)
    except _OutOfWindow:
        return None


def _kg_matches(v: Var, assembly: str, w: Optional[_Window],
                sources: List[Dict[str, Any]]) -> List[Tuple[str, str]]:
    """(rsID, Ensembl allele string) of the dbSNP records that are this allele."""
    chrom, pos, ref, alt = v
    end = pos + len(ref) - 1
    if w is not None:
        try:
            end = max(end, right_end(pos, ref, alt, w.base))
        except _OutOfWindow:
            pass
    region = f"{chrom}:{pos}-{end + 1}"
    resp = request(f"{ENSEMBL[assembly]}/overlap/region/human/{region}?feature=variation",
                   source="Ensembl overlap", cache_ttl=ENSEMBL_TTL)
    data = validated_json(resp, "Ensembl overlap")
    if not isinstance(data, list):
        raise ValueError("Ensembl overlap answer is not a list")
    sources.append(record("Ensembl overlap", f"{assembly} {region} variation", resp,
                          note="dbSNP records at the site, to find the rsID that is this allele"))
    want = (pos, ref, alt)
    out: List[Tuple[str, str]] = []
    for e in data:
        rs = str(e.get("id") or "")
        if e.get("source") != "dbSNP" or not rs.startswith("rs") or str(e.get("seq_region_name")) != chrom:
            continue
        for a in (e.get("alleles") or [])[1:]:
            if _ensembl_vcf(e, str(a), w) == want and (rs, str(a)) not in out:
                out.append((rs, str(a)))
    return out


def _rsids(matches: List[Tuple[str, str]]) -> List[str]:
    out: List[str] = []
    for rs, _ in matches:
        if rs not in out:
            out.append(rs)
    return out


def _kg(v: Var, assembly: str, w: Optional[_Window], matches: List[Tuple[str, str]],
        sources: List[Dict[str, Any]], warnings: List[str]) -> List[Dict[str, Any]]:
    names = [f"1000 Genomes {code}" for code, _, _ in KG_POPS]
    rsids = _rsids(matches)
    if not matches:
        why = "" if w is not None or len(v[2]) == len(v[3]) else \
              " (indels cannot be matched without the reference sequence, which was unavailable)"
        return [_absent(n, _absent_note(n, f"{size} individuals (1000 Genomes phase 3)",
                                        f": no dbSNP record in Ensembl is this allele{why}"))
                for n, (_, _, size) in zip(names, KG_POPS)]
    chosen: Optional[Tuple[str, str, Dict[str, Any], Response]] = None
    with_data = []
    for rs, allele in matches[:3]:
        url = f"{ENSEMBL[assembly]}/variation/human/{rs}?pops=1"
        resp = request(url, source="Ensembl variation", cache_ttl=ENSEMBL_TTL)
        data = validated_json(resp, "Ensembl variation", require="populations")
        sources.append(record("Ensembl variation", f"{rs} ({assembly})", resp,
                              note="1000 Genomes phase 3 allele counts"))
        pops = data.get("populations") or []
        if any(str(p.get("population", "")).startswith("1000GENOMES:phase_3:") for p in pops):
            with_data.append(rs)
            if chosen is None:
                chosen = (rs, allele, data, resp)
    if len(with_data) > 1:
        warnings.append(f"1000 Genomes: {len(with_data)} dbSNP records carry this allele ({', '.join(with_data)}); "
                        f"counts are from {with_data[0]}")
    if chosen is None:
        return [_absent(n, _absent_note(n, f"{size} individuals (1000 Genomes phase 3)",
                                        f": {', '.join(rsids)} has no 1000 Genomes phase 3 genotypes in Ensembl"))
                for n, (_, _, size) in zip(names, KG_POPS)]
    rs, allele, data, resp = chosen
    pops = data.get("populations") or []
    url = f"{ENSEMBL[assembly]}/variation/human/{rs}?pops=1"
    # population rows name alleles in the record's own spelling; an allele outside that set
    # (or a record whose spelling differs from the overlap's) would read as AC=0, so refuse it
    known = set()
    for m in data.get("mappings") or []:
        if str(m.get("seq_region_name")) == v[0]:
            known.update(str(m.get("allele_string") or "").split("/"))
    kg_alleles = {str(p.get("allele")) for p in pops if str(p.get("population", "")).startswith("1000GENOMES:phase_3:")}
    if allele not in known or not kg_alleles <= known:
        raise ValueError(f"{rs}: allele {allele!r} or 1000 Genomes alleles {sorted(kg_alleles)} are not among the "
                         f"record's alleles {sorted(known)}; counts not read")
    out = []
    for name, (code, desc, size) in zip(names, KG_POPS):
        entries = [p for p in pops if p.get("population") == f"1000GENOMES:phase_3:{code}"]
        if not entries:
            out.append(_absent(name, _absent_note(name, f"{size} individuals (1000 Genomes phase 3)",
                                                  f": no {code} genotypes for {rs}")))
            continue
        row: Dict[str, Any] = {"name": name, "population": f"{code}: {desc} (1000 Genomes phase 3, {size} "
                                                          "individuals; small sample)",
                               "url": url, "rsid": rs, "assembly": assembly, "queried_as": vid(v),
                               "samples": size}
        counts = [p.get("allele_count") for p in entries]
        if all(isinstance(c, int) and not isinstance(c, bool) for c in counts) and sum(counts) > 0:
            an = sum(counts)
            ac = sum(p["allele_count"] for p in entries if p.get("allele") == allele)
            row.update(af=ac / an, ac=ac, an=an)
            if ac == 0:
                row["note"] = (f"site genotyped in {code} ({an} alleles called at {rs}); this allele was not "
                               "observed — a count of 0 in a small sample")
        else:
            af = sum(float(p.get("frequency") or 0) for p in entries if p.get("allele") == allele)
            row.update(af=af, ac=None, an=None, note="Ensembl gives frequencies without allele counts here")
        out.append(_found(row))
    return out


# ---------------------------------------------------------------- Taiwan Biobank (TaiwanView)

def _twb(v38: Var, rsids: List[str], w38: Optional[_Window], sources: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    chrom, pos, ref, alt = v38
    seen: List[str] = []
    page = ""
    for rs in rsids[:3]:
        form = {"type": "variant", "database": "WGS", "chr": "", "chrFrom": "", "chrTo": "", "searchBy": "R",
                "geneOrVariant": rs, "page": "1", "pathogenic": "false", "likely_pathogenic": "false",
                "drug_response": "false", "pLoF": "false", "MI": "false", "synonymous_variant": "false",
                "other": "false", "CO": "", "fromV": "", "toV": ""}
        resp = request(TWB_API, source="Taiwan Biobank (TaiwanView)", method="POST",
                       body=urllib.parse.urlencode(form),
                       headers={"Content-Type": "application/x-www-form-urlencoded"},
                       accept="application/json", cache_ttl=STATIC_TTL,
                       validate=_json_list_or_empty("TaiwanView"))
        text = resp.text.strip()
        rows = json.loads(text) if text else []
        if not isinstance(rows, list) or not all(isinstance(r, dict) for r in rows):
            raise ValueError("TaiwanView answer is not a list of records")
        page = f"{TWB_PAGE}?{urllib.parse.urlencode({'database': 'WGS', 'geneOrVariant': rs, 'searchBy': 'R'})}"
        sources.append(record("Taiwan Biobank (TaiwanView)", f"WGS {rs}", resp, url=page,
                              note=f"POST ssss_query.php database=WGS searchBy=R; {len(rows)} record(s)"))
        seen.append(rs)
        for r in rows:
            if str(r.get("chr")) != f"chr{chrom}" or str(r.get("POS")) != str(pos):
                continue
            r_ref, r_alt = str(r.get("REF", "")).upper(), str(r.get("ALT", "")).upper()
            if not (ALLELE_RE.match(r_ref) and ALLELE_RE.match(r_alt)):
                continue
            cand = (int(r["POS"]), r_ref, r_alt)
            if w38 is not None:
                try:
                    cand = left_normalize(cand[0], r_ref, r_alt, w38.base)
                except _OutOfWindow:
                    pass
            if cand != (pos, ref, alt):
                continue
            ns = int(r["NS"])
            row = {"name": TWB, "population": POPULATION[TWB], "af": float(r["AF"]), "ac": None,
                   "an": 2 * ns if chrom not in ("X", "Y", "MT") else None, "url": page, "samples": ns,
                   "rsid": rs, "assembly": "GRCh38", "queried_as": vid(v38),
                   "note": "TaiwanView publishes AF rounded to 4 decimals and NS (samples with data), "
                           "not allele counts"}
            return [_found(row)]
    return [_absent(TWB, _absent_note(TWB, "1,492 WGS individuals (GRCh38)", f" under {', '.join(seen)}"))]


# ---------------------------------------------------------------- lookup

def lookup(chrom: Any, pos: Any, ref: Any, alt: Any, assembly: str = "GRCh38") -> Outcome:
    """Chinese-population frequencies of one allele. Bad input raises UsageError; a failing resource is a warning."""
    given, asm = parse_input(chrom, pos, ref, alt, assembly)
    sources: List[Dict[str, Any]] = []
    warnings: List[str] = []
    results: List[Dict[str, Any]] = []
    not_checked: List[Dict[str, str]] = []

    def skip(name: str, reason: str) -> None:
        not_checked.append({"name": name, "reason": reason})
        warnings.append(f"{name} not queried: {reason}")

    if max(len(given[2]), len(given[3])) > MAX_ALLELE:
        # a valid allele, but a structural variant: not an input error (a variant card must not fail on it)
        for name in [NYUWA, WBBC] + [f"1000 Genomes {c}" for c, _, _ in KG_POPS] + [TWB]:
            skip(name, f"alleles longer than {MAX_ALLELE} bp are structural variants, outside these short-variant "
                       "call sets; not looked up")
        for name, _, reason in NOT_QUERIED:
            not_checked.append({"name": name, "reason": reason})
        return Outcome({"variant": vid(given), "assembly": asm, "datasets": [], "checked": [], "not_found_in": [],
                        "not_checked": not_checked, "notes": list(NOTES)}, sources=sources, warnings=warnings,
                       query={"chrom": given[0], "pos": given[1], "ref": given[2], "alt": given[3], "assembly": asm})

    # 1. REF check + left-alignment in the requested build
    v: Var = given
    w: Optional[_Window] = None
    unaligned = (f"the allele could not be left-aligned within {FLANK} bp (a long repeat or a chromosome end); "
                 "datasets were asked for {} as given, so a 'no record' there may only be a different spelling "
                 "of the same allele")
    try:
        v, w, aligned = _frame(given, asm, sources)
        if not aligned:
            warnings.append(unaligned.format(vid(v)))
    except _RefDiffers as diff:
        other = _other_build_bases(given, asm)
        msg = f"REF {given[2]} does not match the {asm} reference at {given[0]}:{given[1]}, which has {diff.found}."
        if other is not None and other == given[2] and len(given[2]) >= 3:
            msg += f" {OTHER[asm]} has {other} there, which matches: the variant is probably given in {OTHER[asm]}."
        elif other is not None and other == given[2]:
            # one or two bases also match by chance (1 in 4 for a single base): never a build verdict
            msg += (f" {OTHER[asm]} has {other} at the same number, which a short REF matches by chance one time "
                    "in four or more: check the build your report uses before asking again.")
        elif other is not None:
            msg += f" {OTHER[asm]} has {other} there, which does not match either: check position, alleles and build."
        raise UsageError(msg + " Not looked up: a REF that is not the reference is an allele that does not "
                               "exist, and its absence from these datasets would not be evidence.") from None
    except SourceError as err:
        v = (given[0],) + trim(*given[1:])  # type: ignore[assignment]
        warnings.append(f"Ensembl reference sequence unavailable ({err.message}): REF {given[2]} was not checked "
                        f"against {asm} and the allele was not left-aligned; datasets were asked for {vid(v)}")

    # 2. the GRCh38 frame NyuWa and Taiwan Biobank need
    v38: Optional[Var] = v if asm == "GRCh38" else None
    w38: Optional[_Window] = w if asm == "GRCh38" else None
    lifted_from: Optional[str] = None
    no38 = ""
    if asm == "GRCh37":
        if v[0] == "MT":
            v38, w38 = v, None  # rCRS in both builds
        else:
            try:
                p38 = _lift_37_to_38(v, sources)
                v38, w38, aligned38 = _frame((v[0], p38, v[2], v[3]), "GRCh38", sources)
                lifted_from = f"GRCh37 {vid(v)}"
                if not aligned38:
                    warnings.append("GRCh38: " + unaligned.format(vid(v38)))
            except _RefDiffers as diff:
                no38 = (f"GRCh37 {vid(v)} lifts to a GRCh38 position whose reference is {diff.found}, not {v[2]}; "
                        "the allele differs between builds, so GRCh38-only datasets were not asked")
            except (SourceError, ValueError, KeyError, TypeError) as err:
                no38 = f"could not lift GRCh37 {vid(v)} to GRCh38 ({getattr(err, 'message', err)}); " \
                       "GRCh38-only datasets were not asked"

    # 3. the datasets
    if v38 is None:
        skip(NYUWA, f"GRCh38 only; {no38}")
    elif v38[0] == "MT":
        skip(NYUWA, "NyuWa holds chr1-22, X and Y; no mitochondrial variants")
    else:
        got = attempt("NyuWa (NCVD)", lambda: _nyuwa(v38, lifted_from, sources), warnings)  # type: ignore[arg-type]
        if got is None:
            not_checked.append({"name": NYUWA, "reason": warnings[-1]})
        else:
            results += got

    if v[0] in ("X", "Y", "MT"):
        skip(WBBC, "WBBC holds autosomes only (its search takes numeric chromosomes; downloads are chr1-22)")
    else:
        got = attempt("WBBC", lambda: _wbbc(v, asm, sources), warnings)
        if got is None:
            not_checked.append({"name": WBBC, "reason": warnings[-1]})
        else:
            results += got

    # dbSNP record(s) that are this allele: 1000 Genomes counts hang off them, TaiwanView searches by them
    matches = attempt("dbSNP rsID lookup (Ensembl overlap)", lambda: _kg_matches(v, asm, w, sources), warnings)
    kg = None
    if matches is not None:
        kg = attempt("1000 Genomes (Ensembl)", lambda: _kg(v, asm, w, matches, sources, warnings), warnings)  # type: ignore[arg-type]
    if kg is None:
        for code, _, _ in KG_POPS:
            not_checked.append({"name": f"1000 Genomes {code}", "reason": warnings[-1]})
    else:
        results += kg
    rsids = _rsids(matches) if matches is not None else []

    if v38 is None:
        skip(TWB, f"GRCh38 only; {no38}")
    elif matches is None:
        skip(TWB, "searchable by rsID only and the rsID lookup through Ensembl failed")
    elif not rsids:
        skip(TWB, "searchable by rsID only (its region search answers HTTP 500) and no dbSNP rsID is this allele")
    else:
        got = attempt("Taiwan Biobank (TaiwanView)", lambda: _twb(v38, rsids, w38, sources), warnings)  # type: ignore[arg-type]
        if got is None:
            not_checked.append({"name": TWB, "reason": warnings[-1]})
        else:
            results += got

    for name, _, reason in NOT_QUERIED:
        not_checked.append({"name": name, "reason": reason})
    warnings.append("Not queried (no scriptable public interface, or deferred): "
                    + "; ".join(f"{n} ({short})" for n, short, _ in NOT_QUERIED)
                    + " — reasons in result.not_checked")

    result: Dict[str, Any] = {
        "variant": vid(v),
        "assembly": asm,
        "datasets": [r["row"] for r in results if r["status"] == "found"],
        "checked": [r["name"] for r in results],
        "not_found_in": [{"name": r["name"], "note": r["note"]} for r in results if r["status"] == "absent"],
        "not_checked": not_checked,
        "notes": list(NOTES),
    }
    if len(v[2]) == len(v[3]) > 1:
        result["notes"].append("multi-nucleotide variant: these call sets usually store each base change as a "
                               "separate SNV, so a 'no record' here does not cover the individual changes")
    if vid(v) != vid(given):
        result["input"] = vid(given)
        result["notes"].append(f"input {vid(given)} was left-aligned to {vid(v)} before lookup")
    if v38 is not None and asm == "GRCh37":
        result["grch38"] = vid(v38)
    if NYUWA in result["checked"]:
        # privacy: say it where the reader sees it, not only in this module's docstring
        result["notes"].append("NyuWa answers over plain http only: the variant's coordinates "
                               f"({vid(v38) if v38 else vid(v)}) were sent unencrypted to bigdata.ibp.ac.cn")
    counts = {r["row"]["name"]: r["row"].get("an") for r in results if r["status"] == "found"}
    if v[0] == "X" and counts:
        result["notes"].append("chrX allele numbers are counted differently by each dataset (NyuWa counts two "
                               "alleles per person; 1000 Genomes counts one for males): compare AF, not AN")
    return Outcome(result, sources=sources, warnings=warnings,
                   query={"chrom": given[0], "pos": given[1], "ref": given[2], "alt": given[3], "assembly": asm})
