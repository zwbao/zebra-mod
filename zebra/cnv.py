"""The report forms rare-disease families actually hold, besides an SNV/indel list.

Four kinds of result can be entered and analysed here, each with its own result
type. None of them is classified by zebra: what is produced are the *inputs* a
clinician or a curator needs, every one of them traceable to Ensembl or ClinGen.

- **cnv** — a CMA / CNV-seq / array result as coordinates or an ISCN string
  (`arr[GRCh38] 15q11.2q13.1(23123715_28193120)x1`). Reports the genes the
  interval spans (Ensembl overlap), ClinGen dosage sensitivity for them, and
  the ACMG/ClinGen CNV evidence *inputs* (Riggs 2020 sections 1-5).
- **exon_cnv** — an MLPA-style exon deletion or duplication
  (`DMD exon 45-50 deletion`, `NM_004006.3:c.6439-?_7309+?del`). Resolves the
  exon coordinates on the reference transcript and states the frame
  consequence, plus which flanking exon's removal would restore the frame
  (arithmetic only).
- **copy_number** — an SMN1/SMN2-style copy-number result. Structure and
  mechanism only; no prediction, no thresholds.
- **repeat_expansion** — a repeat-expansion result (`FMR1 CGG 230`). Structure
  and mechanism only; the interpretation thresholds are gene-specific and are
  not carried here.

Nothing is inferred from the input: a form that does not carry coordinates (a
bare ISCN band, `del(15)(q11.2q13.1)`) is refused rather than guessed.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

from zebra.core import Outcome, UsageError, attempt
from zebra.http import get_json, request, source_record
from zebra.sources import ensembl
from zebra.vcf import CHROM_LENGTHS  # Ensembl assembly lengths, checked for both builds

ACMG_CNV_FRAMEWORK = ("ACMG/ClinGen technical standard for the interpretation and reporting of constitutional "
                      "copy-number variants (Riggs et al., Genet Med 2020;22:245-257)")
CLINGEN_DOSAGE_MAX_GENES = 40  # how many spanned genes are looked up in ClinGen before the list is cut
GENE_SYMBOL_LIMIT = 200
SMALL_VARIANT_BP = 50  # below this, with exact breakpoints, it is a small variant, not an exon event

_NUM = r"[\d,_]+"
_ISCN_RE = re.compile(
    r"(?:arr|seq)\s*\[\s*(?P<build>GRCh3[78]|hg19|hg38)\s*\]\s*"
    r"(?P<band>[0-9XY]{1,2}[pq][\d.]*(?:[pq][\d.]*)?)?\s*"
    r"\(\s*(?P<start>" + _NUM + r")\s*[_-]\s*(?P<end>" + _NUM + r")\s*\)\s*x\s*(?P<cn>\d+)"
    r"(?P<cn_rest>\s*[~-]\s*\d+)?(?P<tail>.*)$",
    re.I)
# the separator never includes "." — "chr22:11.21-11.23" is a band pair, not coordinates
_REGION_RE = re.compile(
    r"^(?:chr)?(?P<chrom>[0-9]{1,2}|X|Y|MT)[:\s]\s*(?P<start>" + _NUM + r")\s*(?:-|--|_|\.\.)\s*"
    r"(?P<end>" + _NUM + r")(?P<rest>.*)$", re.I)
_EXON_RE = re.compile(
    r"^(?P<gene>[A-Za-z][A-Za-z0-9._-]{0,19})\s*(?::|\s)\s*(?:exons?|ex)\s*\.?\s*"
    r"(?P<first>\d{1,3})\s*(?:[-_]|\s*to\s*|\s*–\s*)?\s*(?P<last>\d{1,3})?\s*"
    r"(?P<type>deletion|del|duplication|dup|loss|gain)?\s*$", re.I)
_HGVS_EXON_RE = re.compile(
    r"^(?P<ref>[A-Za-z0-9_.]+)(?:\([A-Za-z0-9_.-]+\))?:c\.(?P<start>\d+)(?:[+-]\?|[+-]\d+)?_"
    r"(?P<end>\d+)(?:[+-]\?|[+-]\d+)?(?P<type>del|dup)(?:[A-Za-z]*)$", re.I)
_COPY_RE = re.compile(
    r"^(?P<gene>[A-Za-z][A-Za-z0-9._-]{0,19})\s*(?::|\s)\s*"
    r"(?:(?:exon\s*(?P<exon>[\d,and ]+?)\s*)?(?:copy\s*number|copies|cn|x)\s*[:=]?\s*(?P<cn>\d+)"
    r"|(?P<cn2>\d+)\s*cop(?:y|ies))\s*$", re.I)
_REPEAT_RE = re.compile(
    r"^(?P<gene>[A-Za-z][A-Za-z0-9._-]{0,19})\s*(?::|\s)\s*\(?(?P<motif>[ACGTacgt]{2,10})\)?n?\s*"
    r"(?:repeats?\s*[:=]?\s*)?(?P<count>\d{1,5})(?:\s*[-–_]\s*(?P<count2>\d{1,5}))?\s*"
    r"(?:repeats?|units?)?\s*$", re.I)
# only the report's own words; x-notation goes through _type_from_copies, which
# knows that one copy of X is normal in a male
_TRANSCRIPT_RE = re.compile(r"^(?:N[MR]_\d+(?:\.\d+)?|ENST\d+(?:\.\d+)?|LRG_\d+t\d+)$", re.I)
_TYPE_WORDS = {"del": "loss", "deletion": "loss", "loss": "loss", "deleted": "loss",
               "dup": "gain", "duplication": "gain", "gain": "gain", "amplification": "gain"}
_BUILDS = {"grch38": "GRCh38", "hg38": "GRCh38", "grch37": "GRCh37", "hg19": "GRCh37"}

FORMS = (
    'a CNV by coordinates:      "chr15:23123715-28193120 loss" (add --copies 1)',
    'a CNV as ISCN:             "arr[GRCh38] 15q11.2q13.1(23123715_28193120)x1"',
    'an exon deletion/dup:      "DMD exon 45-50 deletion"',
    'the same in HGVS:          "NM_004006.3:c.6439-?_7309+?del"',
    'a copy-number result:      "SMN1 exon 7 copy number 0"',
    'a repeat expansion:        "FMR1 CGG 230"',
)


def _int(text: str) -> int:
    return int(re.sub(r"[,_\s]", "", text))


def _type_from_copies(chrom: Optional[str], copies: Optional[int]) -> Tuple[str, Optional[str]]:
    """Loss or gain from a copy number — which on X and Y depends on the proband's sex."""
    if copies is None:
        return "unknown", None
    if chrom in ("X", "Y"):
        if copies == 0:
            return "loss", None
        if copies >= 3:
            return "gain", None
        return "unknown", (f"{copies} cop{'y' if copies == 1 else 'ies'} of {chrom} is not a loss or a gain by "
                           f"itself: one copy is normal in a male and a loss in a female, two copies the other way "
                           "round. Use the word the report prints (loss/gain/deletion/duplication), or give the "
                           "proband's sex to a clinician reading this.")
    if copies < 2:
        return "loss", None
    if copies > 2:
        return "gain", None
    return "unknown", "two copies is the normal autosomal state: the report must say what it called abnormal"


def parse(text: str) -> Dict[str, Any]:
    """Which result form this is, and the fields it carries. Nothing is guessed."""
    raw = (text or "").strip()
    if not raw:
        raise UsageError("give a result to record. Accepted forms:\n  " + "\n  ".join(FORMS))
    low = raw.lower()

    iscn = _ISCN_RE.match(raw.strip())
    if iscn:
        if iscn.group("cn_rest"):
            raise UsageError(f"{raw!r} gives a copy-number range ({iscn.group('cn')}{iscn.group('cn_rest').strip()}), "
                             "which an array reports for a mosaic result. zebra will not reduce it to one number: "
                             "record it with `zebra case add-variant --kind cnv --iscn` and the range as the report "
                             "writes it.")
        tail = (iscn.group("tail") or "").strip()
        if tail:
            raise UsageError(f"{raw!r} carries {tail!r} after the copy number; that often marks mosaicism (mos) or a "
                             "second finding. Give one finding at a time, exactly as the report prints it.")
        cn = int(iscn.group("cn"))
        chrom = (iscn.group("band") or "")[:2].rstrip("pq").upper() or None
        kind, note = _type_from_copies(chrom, cn)
        return {"kind": "cnv", "input": raw, "iscn": raw, "chrom": chrom,
                "band": iscn.group("band"), "start": _int(iscn.group("start")), "end": _int(iscn.group("end")),
                "copy_number": cn, "cnv_type": kind, "cnv_type_note": note,
                "assembly": _BUILDS[iscn.group("build").lower()]}

    region = _REGION_RE.match(raw)
    if region and not _EXON_RE.match(raw):
        rest = (region.group("rest") or "").strip().lower()
        cn: Optional[int] = None
        cn_match = re.search(r"\bx\s*(\d+)\b|\bcn\s*[:=]?\s*(\d+)\b|\b(\d+)\s*cop(?:y|ies)\b", rest)
        if cn_match:
            cn = int(next(g for g in cn_match.groups() if g))
        chrom = region.group("chrom").upper()
        leftover = (region.group("rest") or "").strip()
        if re.match(r"^[-_.]?\d", leftover):
            raise UsageError(
                f"{raw!r} does not read as one interval: {leftover!r} is left over after "
                f"{chrom}:{_int(region.group('start'))}-{_int(region.group('end'))}. A cytogenetic band pair such as "
                "22q11.21-q11.23 is not a base-pair interval; give the coordinates the report prints.")
        kind = None
        for word, mapped in _TYPE_WORDS.items():
            if re.search(rf"\b{word}\b", rest):
                kind = mapped  # the report's own word always wins
                break
        note = None
        if kind is None:
            kind, note = _type_from_copies(chrom, cn)
        return {"kind": "cnv", "input": raw, "iscn": None, "chrom": chrom,
                "band": None, "start": _int(region.group("start")), "end": _int(region.group("end")),
                "copy_number": cn, "cnv_type": kind or "unknown", "cnv_type_note": note, "assembly": None}

    hgvs = _HGVS_EXON_RE.match(raw)
    if hgvs:
        start, end = int(hgvs.group("start")), int(hgvs.group("end"))
        if end < start:
            raise UsageError(f"{raw}: the second c. position is before the first")
        if "?" not in raw and end - start + 1 < SMALL_VARIANT_BP:
            raise UsageError(f"{raw} describes a {end - start + 1} bp change with exact breakpoints: that is a small "
                             "variant, not an exon-level event. Use `zebra variant` for it (this command is for "
                             "exon-boundary descriptions, which carry '?' for the unsequenced intronic breakpoints).")
        return {"kind": "exon_cnv", "input": raw, "gene": None, "transcript": hgvs.group("ref"),
                "first": None, "last": None, "cds_start": start, "cds_end": end,
                "cnv_type": _TYPE_WORDS[hgvs.group("type").lower()]}

    exon = _EXON_RE.match(raw)
    if exon:
        first = int(exon.group("first"))
        last = int(exon.group("last")) if exon.group("last") else first
        if last < first:
            raise UsageError(f"{raw}: exon {last} is before exon {first}")
        if first < 1:
            raise UsageError(f"{raw}: exons are numbered from 1")
        word = (exon.group("type") or "").lower()
        token = exon.group("gene")
        is_transcript = bool(_TRANSCRIPT_RE.match(token))
        return {"kind": "exon_cnv", "input": raw,
                "gene": None if is_transcript else token.upper(),
                "transcript": token if is_transcript else None,
                "first": first, "last": last, "cds_start": None, "cds_end": None,
                "cnv_type": _TYPE_WORDS.get(word, "unknown")}

    copies = _COPY_RE.match(raw)
    if copies and not _REPEAT_RE.match(raw):
        cn = int(copies.group("cn") or copies.group("cn2"))
        return {"kind": "copy_number", "input": raw, "gene": copies.group("gene").upper(),
                "exon": (copies.group("exon") or "").strip() or None, "copy_number": cn}

    repeat = _REPEAT_RE.match(raw)
    if repeat:
        count: Any = int(repeat.group("count"))
        if repeat.group("count2"):
            count = f"{repeat.group('count')}-{repeat.group('count2')}"
        return {"kind": "repeat_expansion", "input": raw, "gene": repeat.group("gene").upper(),
                "motif": repeat.group("motif").upper(), "repeat_count": count}

    if re.match(r"^(del|dup)\s*\(", low) or re.search(r"\b[0-9XY]{1,2}[pq][\d.]+\b", low) \
            or re.search(r"[:\s]\d{1,2}\.\d", low):
        raise UsageError(
            f"{raw!r} names a band but no coordinates, and zebra will not guess them: cytogenetic bands are not "
            "base-pair boundaries. Give the interval the report prints, e.g. "
            '"chr15:23123715-28193120 loss" or the full ISCN string with coordinates.')
    raise UsageError(f"cannot read {raw!r}. Accepted forms:\n  " + "\n  ".join(FORMS))


# ------------------------------------------------------------------- CNV

OVERLAP_WINDOW = 4_500_000  # Ensembl /overlap/region refuses more than 5 Mb per request
MAX_REGION_BP = 50_000_000  # a CNV larger than this is a typo, not a finding to enumerate gene by gene


def genes_in_region(chrom: str, start: int, end: int, assembly: str = "GRCh38") -> Outcome:
    """Every gene Ensembl places in [start, end] on `chrom` (in windows, because the endpoint caps at 5 Mb)."""
    if end - start + 1 > MAX_REGION_BP:
        raise UsageError(f"{chrom}:{start}-{end} is {(end - start + 1) / 1e6:.1f} Mb; zebra will not walk more than "
                         f"{MAX_REGION_BP / 1e6:.0f} Mb of Ensembl gene annotation one window at a time. Check the "
                         "coordinates (a digit too many is the usual cause), or query the genes you care about with "
                         "`zebra gene`.")
    by_id: Dict[str, Dict[str, Any]] = {}
    sources: List[Dict[str, Any]] = []
    window_start = start
    while window_start <= end:
        window_end = min(end, window_start + OVERLAP_WINDOW - 1)
        region = f"{chrom}:{window_start}-{window_end}"
        resp = get_json(f"{ensembl.host(assembly)}/overlap/region/human/{region}", source="Ensembl overlap",
                        params={"feature": "gene"}, cache_ttl=30 * 86400, timeout=120)
        data = resp.json()
        if not isinstance(data, list):
            raise ValueError("Ensembl overlap returned no list")
        for g in data:
            if not g.get("id"):
                continue
            by_id[g["id"]] = {"symbol": g.get("external_name") or g.get("id"), "ensembl_id": g.get("id"),
                              "biotype": g.get("biotype"), "start": g.get("start"), "end": g.get("end"),
                              "strand": g.get("strand"), "named": bool(g.get("external_name"))}
        sources.append(source_record("Ensembl overlap", f"{assembly} {region}", resp, note="feature=gene"))
        window_start = window_end + 1
    genes = sorted(by_id.values(), key=lambda g: (g["start"] or 0))
    return Outcome(genes, sources=sources)


def _dosage_rows(symbols: List[str], assembly: str, warnings: List[str],
                 sources: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """ClinGen dosage records for the genes a CNV spans.

    ClinGen publishes one bulk TSV, so it is fetched once and parsed for every
    gene rather than requested per gene (which would re-download it whenever
    the cache is off).
    """
    from zebra.sources import clingen

    rows: List[Dict[str, Any]] = []
    if not symbols:
        return rows
    url = clingen.DOSAGE_TSV.get(assembly)
    bulk = attempt(
        f"ClinGen dosage map ({assembly})",
        lambda: request(url, source="ClinGen dosage", accept="text/tab-separated-values,text/plain,*/*",
                        cache_ttl=clingen.CACHE_TTL, timeout=120),
        warnings) if url else None
    if bulk is not None:
        sources.append(source_record("ClinGen dosage sensitivity", f"{len(symbols)} spanned genes", bulk,
                                     note="bulk TSV, filtered per gene"))
        for symbol in symbols:
            row = attempt(f"ClinGen dosage {symbol}",
                          lambda s=symbol: clingen.parse_dosage_tsv(bulk.text, s), warnings)
            if row:
                row["url"] = f"https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen/clingen_gene.cgi?sym={symbol}"
                rows.append(row)
        return rows
    # the bulk file was unreachable: fall back to the per-gene helper, which
    # may still be served from the cache
    seen_source = False
    for symbol in symbols:
        got = attempt(f"ClinGen dosage {symbol}", lambda s=symbol: clingen.dosage(s, assembly), warnings)
        if got is None:
            continue
        if not seen_source:
            sources.extend(got.sources)
            seen_source = True
        if got.result:
            rows.append(got.result)
    return rows


def _dosage_split(rows: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]],
                                                       List[Dict[str, Any]]]:
    """ClinGen's established HI (score 3) and TS (score 3) genes, plus every gene with any score at all."""
    hi, ts, scored = [], [], []
    for row in rows:
        h = (row.get("haploinsufficiency") or {}).get("score")
        t = (row.get("triplosensitivity") or {}).get("score")
        if str(h) == "3":
            hi.append({"gene": row.get("gene"), "score": h,
                       "description": (row.get("haploinsufficiency") or {}).get("description")})
        if str(t) == "3":
            ts.append({"gene": row.get("gene"), "score": t,
                       "description": (row.get("triplosensitivity") or {}).get("description")})
        if (h not in (None, "", "0")) or (t not in (None, "", "0")):
            scored.append({"gene": row.get("gene"), "haploinsufficiency_score": h, "triplosensitivity_score": t})
    return hi, ts, scored


def cnv_card(parsed: Dict[str, Any], assembly: str = "GRCh38", inheritance: Optional[str] = None) -> Outcome:
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
    length_of_chrom = CHROM_LENGTHS.get(build, {}).get(chrom)
    if length_of_chrom is None:
        raise UsageError(f"{chrom!r} is not a chromosome zebra can place in {build} "
                         "(1-22, X, Y); check the report's notation")
    if end > length_of_chrom:
        raise UsageError(f"{chrom}:{start}-{end} runs past the end of chromosome {chrom} in {build} "
                         f"({length_of_chrom:,} bp): check the build and the coordinates")
    if parsed.get("assembly") and parsed["assembly"] != assembly:
        warnings.append(f"the report states {parsed['assembly']} and --assembly says {assembly}: the "
                        f"{parsed['assembly']} coordinates in the report were used, because builds must never be "
                        "mixed")
        assembly = parsed["assembly"]

    got = genes_in_region(chrom, start, end, assembly)
    sources.extend(got.sources)
    genes: List[Dict[str, Any]] = got.result
    coding = [g for g in genes if g["biotype"] == "protein_coding"]
    # one symbol can carry two Ensembl gene ids (a readthrough locus): look it up once
    named_coding = list(dict.fromkeys(g["symbol"] for g in coding if g["named"]))
    dosage_targets = named_coding[:CLINGEN_DOSAGE_MAX_GENES]
    if len(named_coding) > CLINGEN_DOSAGE_MAX_GENES:
        warnings.append(f"{len(named_coding)} protein-coding genes are spanned; ClinGen dosage was looked up for the "
                        f"first {CLINGEN_DOSAGE_MAX_GENES} by position. Query the rest with `zebra gene <symbol>`.")
    rows = _dosage_rows(dosage_targets, assembly, warnings, sources)
    hi, ts, scored = _dosage_split(rows)

    cn = parsed.get("copy_number")
    kind = parsed.get("cnv_type") or "unknown"
    length = end - start + 1
    result: Dict[str, Any] = {
        "kind": "cnv",
        "input": parsed["input"],
        "iscn": parsed.get("iscn"),
        "assembly": assembly,
        "region": {"chrom": chrom, "start": start, "end": end, "length_bp": length,
                   "band_as_reported": parsed.get("band")},
        "cnv_type": kind,
        "copy_number": cn,
        "genes": {"total": len(genes), "protein_coding": len(coding),
                  "protein_coding_symbols": named_coding[:GENE_SYMBOL_LIMIT],
                  "unnamed_protein_coding": len([g for g in coding if not g["named"]]),
                  "other_biotypes": len(genes) - len(coding),
                  "symbols_truncated": len(named_coding) > GENE_SYMBOL_LIMIT},
        "clingen_dosage": {"genes_checked": len(dosage_targets), "genes_with_a_record": len(rows), "records": rows},
        "acmg_cnv_inputs": {
            "framework": ACMG_CNV_FRAMEWORK,
            "scope": kind if kind in ("loss", "gain") else
                     ("unknown: " + parsed["cnv_type_note"] if parsed.get("cnv_type_note")
                      else "unknown (the report does not say whether this is a loss or a gain)"),
            "section_1_variant_type": {
                "contains_protein_coding_genes": bool(coding),
                "protein_coding_gene_count": len(coding),
                "note": "section 1 asks only whether protein-coding or other important elements are contained",
            },
            "section_2_overlap_with_established_regions_or_genes": {
                "clingen_established_haploinsufficient_genes": hi,
                "clingen_established_triplosensitive_genes": ts,
                "clingen_genes_with_any_dosage_score": scored,
                "genes_checked_in_clingen": len(dosage_targets),
                "note": "zebra checked the ClinGen dosage-sensitivity map gene by gene; established recurrent "
                        "*regions* (and the breakpoint rules for partial overlap) must be read from the ClinGen "
                        "dosage map itself",
            },
            "section_3_gene_number": {
                "protein_coding_genes": len(coding),
                "clingen_bands": "the standard's section 3 bands (3A/3B/3C) are <25, 25-34 and >=35 "
                                 "protein-coding genes; see the framework citation above for the wording",
            },
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
    if cn is None:
        warnings.append("no copy number in the input: pass --copies, or write the ISCN x-notation, so a loss can be "
                        "told from a homozygous loss")
    if parsed.get("cnv_type_note"):
        warnings.append(parsed["cnv_type_note"])
    lines = [f"CNV {chrom}:{start}-{end} ({assembly}, {length:,} bp) {kind}"
             + (f", copy number {cn}" if cn is not None else "")]
    lines.append(f"genes spanned: {len(genes)} total, {len(coding)} protein-coding"
                 + (": " + ", ".join(named_coding[:15]) + (" …" if len(named_coding) > 15 else "")
                    if named_coding else ""))
    lines.append(f"ClinGen dosage: {len(rows)} of {len(dosage_targets)} checked genes have a record; "
                 f"established HI: {', '.join(x['gene'] for x in hi) or 'none'}; "
                 f"established TS: {', '.join(x['gene'] for x in ts) or 'none'}")
    if scored:
        lines.append("  any dosage score: " + ", ".join(
            f"{x['gene']} HI {x['haploinsufficiency_score'] or '-'}/TS {x['triplosensitivity_score'] or '-'}"
            for x in scored[:10]))
    lines.append("ACMG/ClinGen CNV inputs reported; no classification (sections 4-5 need case-level data)")
    return Outcome(result, sources=sources, warnings=warnings, text="\n".join(lines),
                   query={"input": parsed["input"], "assembly": assembly, "inheritance": inheritance})


# -------------------------------------------------------------- exon CNV

def _transcript(gene: str, assembly: str, wanted: Optional[str]) -> Tuple[Dict[str, Any], Outcome, List[str]]:
    """The reference transcript whose exon numbering is used, with its exons."""
    notes: List[str] = []
    got = ensembl.lookup_symbol(gene, assembly, expand=True)
    data = got.result
    if not isinstance(data, dict) or not data.get("Transcript"):
        raise UsageError(f"Ensembl has no transcripts for {gene!r} in {assembly}; check the gene symbol")
    txs = [t for t in data["Transcript"] if t.get("Exon")]
    chosen = None
    if wanted:
        base = wanted.split(".")[0].upper()
        chosen = next((t for t in txs if str(t.get("id", "")).split(".")[0].upper() == base), None)
        if chosen is None:
            notes.append(f"the report names {wanted}, which is not an Ensembl transcript id: the exon numbering below "
                         "is the Ensembl canonical transcript's. Check that the two number exons the same way.")
    if chosen is None:
        chosen = next((t for t in txs if t.get("is_canonical")), None)
    if chosen is None:
        chosen = max(txs, key=lambda t: t.get("length") or 0)
        notes.append("no canonical transcript was marked; the longest one was used")
    return chosen, got, notes


def _exons_in_order(tx: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The transcript's exons in transcript order (5' to 3'), which is exon numbering order."""
    return sorted(tx["Exon"], key=lambda e: e["start"], reverse=(tx.get("strand") or 1) < 0)


def exon_card(parsed: Dict[str, Any], assembly: str = "GRCh38", gene: Optional[str] = None) -> Outcome:
    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []
    kind = parsed.get("cnv_type") or "unknown"
    symbol = parsed.get("gene") or gene
    result: Dict[str, Any] = {"kind": "exon_cnv", "input": parsed["input"], "assembly": assembly,
                              "gene": symbol, "cnv_type": kind}

    if parsed.get("cds_start") is not None:
        # the HGVS exon-boundary form carries the deleted CDS span itself
        span = parsed["cds_end"] - parsed["cds_start"] + 1
        result["transcript"] = {"as_reported": parsed.get("transcript"), "resolved": None,
                               "exon_numbering_source": "not needed: the c. positions give the span"}
        result["cds_span"] = {"from": parsed["cds_start"], "to": parsed["cds_end"], "length_bp": span}
        result["frame"] = _frame(span, kind, basis=f"c.{parsed['cds_start']}_{parsed['cds_end']} spans {span} coding "
                                                   f"bases on {parsed.get('transcript')} as the report writes it")
        result["coordinates"] = None
        warnings.append("the '?' in an exon-boundary HGVS description means the intronic breakpoints were not "
                        "sequenced: the span is the coding sequence lost, which is what the reading frame depends on")
        if symbol:
            got = ensembl.lookup_symbol(symbol, assembly)
            sources.extend(got.sources)
            result["gene_location"] = {k: got.result.get(k) for k in ("seq_region_name", "start", "end", "strand", "id")}
        return Outcome(result, sources=sources, warnings=warnings, text=_exon_text(result),
                       query={"input": parsed["input"], "assembly": assembly})

    if not symbol:
        raise UsageError("which gene? write the symbol before the exons (\"DMD exon 45-50 deletion\"), or pass "
                         "--gene when the report names only a transcript")
    tx, got, notes = _transcript(symbol, assembly, parsed.get("transcript"))
    sources.extend(got.sources)
    warnings.extend(notes)
    exons = _exons_in_order(tx)
    first, last = parsed["first"], parsed["last"]
    if last > len(exons):
        raise UsageError(f"{symbol} transcript {tx['id']} has {len(exons)} exons; the report names exon {last}. "
                         "Check which transcript the report numbers against.")
    picked = exons[first - 1:last]
    lengths = [e["end"] - e["start"] + 1 for e in picked]
    total = sum(lengths)
    result["transcript"] = {"as_reported": parsed.get("transcript"), "resolved": tx.get("id"),
                            "name": tx.get("display_name"), "biotype": tx.get("biotype"),
                            "exon_total": len(exons), "strand": tx.get("strand"),
                            "exon_numbering_source": "Ensembl canonical transcript"
                            if tx.get("is_canonical") else "Ensembl transcript as named"}
    result["exons"] = {"first": first, "last": last, "count": len(picked),
                       "per_exon": [{"exon": first + i, "chrom": e.get("seq_region_name"), "start": e["start"],
                                     "end": e["end"], "length_bp": lengths[i], "ensembl_id": e.get("id")}
                                    for i, e in enumerate(picked)]}
    result["coordinates"] = {"chrom": picked[0].get("seq_region_name"),
                             "start": min(e["start"] for e in picked), "end": max(e["end"] for e in picked),
                             "note": "the exon boundaries; the real breakpoints lie in the flanking introns"}
    edge = []
    if first == 1:
        edge.append("exon 1 carries the 5' UTR, so its length is not all coding sequence")
    if last == len(exons):
        edge.append("the last exon carries the 3' UTR, so its length is not all coding sequence")
    result["frame"] = _frame(total, kind, basis=f"sum of exon lengths {first}-{last} on {tx.get('id')} "
                                                f"({'+'.join(str(x) for x in lengths)} = {total} bp)",
                             caveats=edge)
    result["frame_restoration"] = _restoration(exons, first, last, total) if kind == "loss" else None
    if edge:
        warnings.extend(edge)
    return Outcome(result, sources=sources, warnings=warnings, text=_exon_text(result),
                   query={"input": parsed["input"], "assembly": assembly, "gene": symbol})


def _frame(span: int, kind: str, basis: str, caveats: Optional[List[str]] = None) -> Dict[str, Any]:
    rest = span % 3
    consequence = "in frame" if rest == 0 else "out of frame"
    verb = "removed" if kind == "loss" else ("added" if kind == "gain" else "affected")
    tail = ("a multiple of 3: the reading frame downstream is preserved" if rest == 0 else
            f"not a multiple of 3 ({rest} over): the reading frame downstream shifts")
    return {
        "bases": span, "modulo_3": rest, "consequence": consequence,
        "reading": f"{span} bases {verb} is {tail}",
        "basis": basis,
        "caveats": (caveats or []) + [
            "frame arithmetic only: whether the protein that results is functional, and whether the transcript "
            "escapes nonsense-mediated decay, is not predicted here",
            "the arithmetic assumes the whole span is coding sequence of this transcript",
        ],
    }


def _restoration(exons: List[Dict[str, Any]], first: int, last: int, total: int) -> Dict[str, Any]:
    """Which flanking exon's additional removal would make the deletion a multiple of 3."""
    options = []
    for n in (first - 1, last + 1):
        if 1 <= n <= len(exons):
            e = exons[n - 1]
            length = e["end"] - e["start"] + 1
            options.append({"exon": n, "length_bp": length, "total_bp": total + length,
                            "restores_frame": (total + length) % 3 == 0})
    return {
        "candidates": options,
        "note": "arithmetic only: these are the flanking exons whose additional removal (by an exon-skipping "
                "approach or a larger deletion) would make the lost span a multiple of 3. Whether an approved or "
                "investigational therapy exists for that exon, and whether this patient would be eligible, is a "
                "separate question (`zebra therapy`, `zebra trials`).",
    }


def _exon_text(result: Dict[str, Any]) -> str:
    frame = result.get("frame") or {}
    lines = [f"{result.get('gene') or result['input']} {result['cnv_type']}"
             + (f" exon {result['exons']['first']}-{result['exons']['last']}" if result.get("exons") else "")]
    tx = result.get("transcript") or {}
    if tx.get("resolved"):
        lines.append(f"transcript {tx['resolved']} ({tx.get('name')}), {tx.get('exon_total')} exons, "
                     f"{tx.get('exon_numbering_source')}")
    elif tx.get("as_reported"):
        lines.append(f"transcript as reported: {tx['as_reported']}")
    coords = result.get("coordinates")
    if coords:
        lines.append(f"exon coordinates: {coords['chrom']}:{coords['start']}-{coords['end']} ({result['assembly']})")
    lines.append(f"frame: {frame.get('bases')} bp, mod 3 = {frame.get('modulo_3')} → {frame.get('consequence')}")
    lines.append(f"  {frame.get('reading')}")
    rest = result.get("frame_restoration")
    if rest:
        for c in rest["candidates"]:
            lines.append(f"  + exon {c['exon']} ({c['length_bp']} bp) → {c['total_bp']} bp, "
                         f"{'restores the frame' if c['restores_frame'] else 'still out of frame'}")
    return "\n".join(lines)


# ------------------------------------------- copy number / repeat expansion

COPY_NUMBER_MECHANISM = {
    "SMN1": [
        "most SMN1 assays (MLPA, ddPCR, qPCR) read the copy number of exon 7, because SMN1 and SMN2 differ at only "
        "a few bases; the result is a count of SMN1 copies, not a sequence",
        "a copy-number count cannot tell two copies on one chromosome (2+0) from one on each (1+1), so a carrier "
        "result does not exclude the 2+0 'silent carrier' arrangement",
        "SMN2 copy number is reported alongside as a modifier of the SMN1 finding; it is not itself the finding",
        "a copy-number assay does not see an intragenic SMN1 point variant: a 1-copy result with a matching "
        "phenotype is usually followed by sequencing of the remaining copy",
    ],
}
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


def finding_card(parsed: Dict[str, Any], assembly: str = "GRCh38", method: Optional[str] = None,
                 related: Optional[List[str]] = None) -> Outcome:
    """A copy-number or repeat-expansion result, recorded as structure plus mechanism. Nothing is predicted."""
    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []
    symbol = parsed["gene"]
    got = attempt(f"Ensembl lookup {symbol}", lambda: ensembl.lookup_symbol(symbol, assembly), warnings)
    location = None
    if got is not None:
        sources.extend(got.sources)
        location = {k: got.result.get(k) for k in ("seq_region_name", "start", "end", "strand", "id", "biotype")}
    else:
        warnings.append(f"{symbol} could not be confirmed against Ensembl: the symbol is recorded as given")

    result: Dict[str, Any] = {"kind": parsed["kind"], "input": parsed["input"], "gene": symbol,
                              "assembly": assembly, "gene_location": location, "method": method,
                              "classification": None,
                              "classification_note": "zebra records this result and what it means for the "
                                                     "mechanism; it does not classify it or compare it against a "
                                                     "threshold. Read the laboratory's own interpretation."}
    if parsed["kind"] == "copy_number":
        result["copy_number"] = parsed["copy_number"]
        result["exon"] = parsed.get("exon")
        result["related_results"] = []
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
        result["mechanism"] = COPY_NUMBER_MECHANISM.get(symbol, []) + COPY_NUMBER_GENERIC
        if symbol not in COPY_NUMBER_MECHANISM:
            warnings.append(f"zebra carries no gene-specific mechanism notes for {symbol}: only the general ones "
                            "about copy-number assays are shown")
        text = [f"{symbol} copy number {parsed['copy_number']}"
                + (f" (exon {parsed['exon']})" if parsed.get("exon") else "")]
        for other in result["related_results"]:
            text.append(f"  with {other['gene']} copy number {other['copy_number']}")
    else:
        result["motif"] = parsed["motif"]
        result["repeat_count"] = parsed["repeat_count"]
        result["mechanism"] = REPEAT_MECHANISM
        result["thresholds"] = None
        result["thresholds_note"] = ("the normal / intermediate / full-expansion boundaries are gene-specific and "
                                     "are not carried in zebra: take them from the laboratory report or from this "
                                     "gene's GeneReviews chapter (`zebra disease`, `zebra gene`).")
        text = [f"{symbol} {parsed['motif']} repeat: {parsed['repeat_count']}"]
    if location:
        text.append(f"  {symbol} at {location['seq_region_name']}:{location['start']}-{location['end']} "
                    f"({assembly}, {location['id']})")
    text.append("  " + "\n  ".join(result["mechanism"]))
    text.append("  no classification: zebra does not compare this result against a threshold")
    return Outcome(result, sources=sources, warnings=warnings, text="\n".join(text),
                   query={"input": parsed["input"], "assembly": assembly, "method": method})


# ------------------------------------------------------------------ entry

def card(text: str, assembly: str = "GRCh38", gene: Optional[str] = None, copies: Optional[int] = None,
         inheritance: Optional[str] = None, method: Optional[str] = None,
         related: Optional[List[str]] = None) -> Outcome:
    """Read one non-SNV result and return what can be established about it."""
    parsed = parse(text)
    if copies is not None and parsed["kind"] in ("cnv", "copy_number"):
        if parsed.get("copy_number") is not None and parsed["copy_number"] != copies:
            raise UsageError(f"the input says copy number {parsed['copy_number']} and --copies says {copies}")
        parsed["copy_number"] = copies
        if parsed["kind"] == "cnv" and parsed.get("cnv_type") in (None, "unknown"):
            parsed["cnv_type"], parsed["cnv_type_note"] = _type_from_copies(parsed.get("chrom"), copies)
    if copies is not None and parsed["kind"] not in ("cnv", "copy_number"):
        raise UsageError(f"--copies applies to a CNV or a copy-number result; {text!r} reads as "
                         f"{parsed['kind']}, where a copy number has no meaning")
    if parsed["kind"] == "cnv":
        cn, kind = parsed.get("copy_number"), parsed.get("cnv_type")
        if cn is not None and kind in ("loss", "gain"):
            expected, _ = _type_from_copies(parsed.get("chrom"), cn)
            if expected in ("loss", "gain") and expected != kind:
                raise UsageError(f"the report says {kind} and the copy number says {expected} (copy number {cn} on "
                                 f"{parsed.get('chrom')}): one of the two is wrong, so zebra will not record either")
        return cnv_card(parsed, assembly=assembly, inheritance=inheritance)
    if parsed["kind"] == "exon_cnv":
        return exon_card(parsed, assembly=assembly, gene=gene)
    return finding_card(parsed, assembly=assembly, method=method, related=related)


# --------------------------------------------------------- case recording

def case_fields(result: Dict[str, Any], method: Optional[str] = None) -> Dict[str, Any]:
    """The finding as a case.json variant record (see `zebra.case.VARIANT_KINDS`)."""
    kind = result["kind"]
    out: Dict[str, Any] = {"kind": kind, "assembly": result.get("assembly"), "method": method,
                           "description": result["input"], "source": "zebra cnv"}
    if kind == "cnv":
        r = result["region"]
        out.update({"region": f"{r['chrom']}:{r['start']}-{r['end']}", "iscn": result.get("iscn"),
                    "cnv_type": result.get("cnv_type"), "copy_number": result.get("copy_number"),
                    "genes": result["genes"]["protein_coding_symbols"][:50]})
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
            out["note"] = "; ".join(f"{r['gene']} copy number {r['copy_number']}"
                                    + (f" (exon {r['exon']})" if r.get("exon") else "") for r in related)
            out["genes"] = [result.get("gene")] + [r["gene"] for r in related]
    else:
        out.update({"gene": result.get("gene"), "motif": result.get("motif"),
                    "repeat_count": result.get("repeat_count")})
    return {k: v for k, v in out.items() if v is not None}
