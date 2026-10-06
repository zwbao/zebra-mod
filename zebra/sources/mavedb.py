"""MaveDB API: multiplexed assay of variant effect (MAVE) score sets for a gene, and one variant's score.

API: https://api.mavedb.org/api/v1 (no key; OpenAPI at https://api.mavedb.org/openapi.json).

  POST /score-sets/search            {"targets": ["<gene>"], "published": true}
  GET  /score-sets/{urn}/scores      the score table as CSV (start/limit paginate it)
  GET  /score-calibrations/score-set/{urn}/primary
                                     the score set's own primary calibration, when it has one

Two things make a MAVE score dangerous to look up naively, and both are refused here
rather than guessed:

1. **Whose numbering is it?** A score set's target is either an ACCESSION
   (`targetAccession`, e.g. NM_007294.3) or a bare SEQUENCE (`targetSequence`) -- a
   construct, whose HGVS is numbered from position 1 of that construct, not of the
   gene. Of the seven published BRCA1 score sets, one is accession-numbered and the
   rest are sequence targets, two of them 191 aa and 287 aa BRCT fragments. Matching a
   patient's `p.Glu10Lys` against a 191 aa construct's `p.Glu10Lys` returns the score
   of BRCA1 residue 1682 under the label of residue 10. So a sequence-target score set
   is listed but NOT searched by HGVS unless the caller explicitly asks, and any hit
   from one is reported with `numbering: "construct"` and no functional class.
2. **A raw score is not a functional class**, and nothing here derives one. A score is
   labelled only from the score set's own published calibration, with that
   calibration's label, range and ACMG criterion as the dataset stated them (verified
   live 2026-10-06 on urn:mavedb:00000097-0-2, whose "Fayer calibration" gives
   Functional -> BS3 STRONG and Non-functional -> PS3 STRONG, with OddsPath 0.02 and
   52.4). Turning a calibrated assay into PS3 or BS3 for a real classification is the
   ClinGen SVI's published framework (Brnich et al. 2019, Genome Medicine); the
   calibration is an input to that process, not a substitute for it.

A calibration band with no finite bound, a band that overlaps another band with a
different verdict, a non-finite score, and a calibration marked research-use-only all
yield NO criterion, each with its own stated reason.
"""

from __future__ import annotations

import csv
import io
import math
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

from zebra.core import Outcome
from zebra.http import SourceError, get_json, post_json, request
from zebra.sources import attempt, record, validated_json, validated_text

BASE = "https://api.mavedb.org/api/v1"
# Rows read per score set before the scan stops and says so. The published sets run
# from a few hundred to tens of thousands of variants.
MAX_ROWS = 30000
PAGE_ROWS = 10000
MAX_DATASETS_SCORED = 6

BRNICH_NOTE = (
    "A raw score is not a functional class, and nothing here assigns one. Where a score set publishes its own "
    "calibration, that calibration's label, score range and ACMG criterion are reported as the dataset stated "
    "them. Applying PS3 or BS3 to a variant from a functional assay is the ClinGen Sequence Variant "
    "Interpretation working group's published framework (Brnich et al. 2019, Genome Medicine) -- the assay has to "
    "be shown to model the disease mechanism and be calibrated against known controls first, and the evidence "
    "strength comes out of that calibration. Read the framework and the score set's own paper before using a "
    "number here as evidence."
)
NUMBERING_NOTE = (
    "A score set whose target is a bare sequence numbers its variants from position 1 of THAT construct, which is "
    "not the gene's numbering unless the construct happens to be the whole gene. Such a score set is listed here "
    "but is not searched by HGVS by default, because a match would silently return a different residue's score."
)
NO_CALIBRATION_NOTE = (
    "this score set publishes no primary calibration in MaveDB, so the score carries no functional class here; "
    "read the score set's own paper for how its authors drew the line"
)
NO_SCORE_NOTE = (
    "this row carries no usable score (the cell was empty, NA, or not a finite number), so there is nothing to "
    "place in the calibration"
)
RUO_NOTE = (
    "this calibration is marked research-use-only by its depositors, so its ACMG criterion is reported as the "
    "dataset's own label and is NOT offered as usable evidence"
)
# MaveDB answers a missing calibration, a missing score set and a wrong route all with
# HTTP 404; only the first means "this score set has no calibration".
NO_CALIBRATION_404 = re.compile(r"no primary score calibration", re.I)

# One-letter to three-letter amino acid codes, so p.R1699W finds MaveDB's p.Arg1699Trp.
AA1_TO_3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
    "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
    "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "U": "Sec", "O": "Pyl", "B": "Asx", "Z": "Glx",
    # '*' and the legacy 'X' both spell a stop in a one-letter protein change
    "*": "Ter", "X": "Ter",
}
AA3 = {v.upper(): v for v in AA1_TO_3.values()}
AA3["TER"] = "Ter"
_ACCESSION_SPLIT = re.compile(r"^([^:]+):(?=[cgnmrp]\.)", re.I)
_P_CHANGE = re.compile(r"^p\.\(?([A-Za-z]{3}|[A-Za-z*])(\d+)([A-Za-z]{3}|[A-Za-z*=])?\)?$", re.I)
_PREFIX = re.compile(r"^([cgnmrp]\.)(.*)$", re.I)
_BARE_PROTEIN = re.compile(r"^([A-Za-z*])(\d+)([A-Za-z*])$")


def strip_version(accession: Optional[str]) -> Optional[str]:
    return accession.split(".")[0] if accession else accession


def split_hgvs(hgvs: Optional[str]) -> Tuple[Optional[str], Optional[str]]:
    """`NM_007294.3:c.5565A>T` -> ("NM_007294.3", "c.5565A>T"); a bare change -> (None, change)."""
    text = (hgvs or "").strip().replace(" ", "")
    if not text:
        return None, None
    m = _ACCESSION_SPLIT.match(text)
    if m:
        return m.group(1), text[m.end():]
    return None, text


def normalise_change(change: Optional[str]) -> Optional[str]:
    """The comparable form of an HGVS change: no accession, no parentheses, one case.

    Protein changes are expanded from one-letter to three-letter code, because MaveDB
    writes `p.Arg1699Trp` where a report writes `p.R1699W`. `p.=` with no position is
    refused: it is "the protein is unchanged", which matches thousands of rows and
    identifies no variant.
    """
    if not change:
        return None
    text = change.strip().replace(" ", "")
    if not text:
        return None
    bare = _BARE_PROTEIN.match(text)
    if bare:  # "R1699W" with the p. left off
        text = "p." + text
    if text.lower().startswith("p."):
        body = text[2:]
        if body.strip("()") in ("=", ""):
            return None  # p.= / p.(=): no position, so it identifies nothing
        m = _P_CHANGE.match(text)
        if not m:
            return text.upper()
        ref, number, alt = m.group(1), m.group(2), m.group(3)

        def three(token: Optional[str]) -> Optional[str]:
            if token is None or token == "=":
                return token
            if len(token) == 1:
                return AA1_TO_3.get(token.upper())
            return AA3.get(token.upper())

        ref3, alt3 = three(ref), three(alt)
        if ref3 is None or (alt is not None and alt3 is None):
            return text.upper()
        return ("p." + ref3 + number + (alt3 or "")).upper()
    m = _PREFIX.match(text)
    if m:
        return (m.group(1).lower() + m.group(2).upper()).upper()
    return text.upper()


def normalise_hgvs(hgvs: Optional[str]) -> Optional[str]:
    """Back-compatible helper: the normalised change part of an HGVS string."""
    _, change = split_hgvs(hgvs)
    norm = normalise_change(change)
    if norm is None:
        return None
    # return it in MaveDB's own casing for protein changes, which is what callers display
    m = _P_CHANGE.match(norm)
    if norm.upper().startswith("P.") and m:
        ref, number, alt = m.group(1), m.group(2), m.group(3)
        ref3 = AA3.get(ref.upper()) or ref.capitalize()
        alt3 = AA3.get((alt or "").upper()) or ((alt or "").capitalize() or None)
        return "p." + ref3 + number + (alt3 or "")
    pre = _PREFIX.match(norm)
    if pre:   # a coding/genomic change: show the prefix as HGVS writes it
        return pre.group(1).lower() + pre.group(2)
    return norm


def match_hgvs(cell: Optional[str], query: Optional[str]) -> Optional[str]:
    """How `cell` matches `query`: "accession" , "change_only", or None for no match.

    "change_only" means one of the two carried no accession, so the match rests on the
    change alone. For a `c.` change that is a real risk -- coding numbering is
    transcript-specific -- so the caller labels such a hit instead of presenting it as
    the same variant.
    """
    cell_acc, cell_change = split_hgvs(cell)
    query_acc, query_change = split_hgvs(query)
    a, b = normalise_change(cell_change), normalise_change(query_change)
    if a is None or b is None or a != b:
        return None
    if cell_acc and query_acc:
        return "accession" if strip_version(cell_acc) == strip_version(query_acc) else None
    return "change_only"


def _target_rows(score_set: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Each target gene of a score set, with how its variants are numbered."""
    out = []
    for t in score_set.get("targetGenes") or []:
        if not isinstance(t, dict) or not t.get("name"):
            continue
        accession = (t.get("targetAccession") or {}) if isinstance(t.get("targetAccession"), dict) else {}
        sequence = (t.get("targetSequence") or {}) if isinstance(t.get("targetSequence"), dict) else {}
        acc = accession.get("accession")
        out.append({
            "name": str(t["name"]),
            "category": t.get("category"),
            "accession": acc,
            "assembly": accession.get("assembly"),
            "numbering": "accession" if acc else "construct",
            "sequence_type": sequence.get("sequenceType"),
            "sequence_length": len(sequence.get("sequence") or "") or None,
        })
    return out


def _publication(score_set: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    for p in score_set.get("primaryPublicationIdentifiers") or score_set.get("publicationIdentifiers") or []:
        if isinstance(p, dict):
            return {"pmid": p.get("pubmedId"), "doi": p.get("doi"), "title": p.get("title"),
                    "journal": p.get("publicationJournal"), "year": p.get("publicationYear"),
                    "url": p.get("url") or p.get("referenceHtml")}
    return None


def score_set_url(urn: Optional[str]) -> str:
    return "https://mavedb.org/score-sets/%s" % (urn or "")


def search_score_sets(gene: str, published: bool = True) -> Outcome:
    """Published MaveDB score sets whose target gene is `gene`."""
    gene = (gene or "").strip()
    if not gene:
        raise ValueError("give a gene symbol")
    payload: Dict[str, Any] = {"targets": [gene]}
    if published:
        payload["published"] = True
    resp = post_json("%s/score-sets/search" % BASE, payload, source="MaveDB score-sets/search",
                     cache_ttl=7 * 86400, timeout=90)
    data = validated_json(resp, "MaveDB score-sets/search")
    rows = data.get("scoreSets") if isinstance(data, dict) else data
    if not isinstance(rows, list):
        raise ValueError("MaveDB score-sets/search: no scoreSets list in the response")
    declared = data.get("numScoreSets") if isinstance(data, dict) else None
    wanted = gene.upper()
    datasets: List[Dict[str, Any]] = []
    warnings: List[str] = []
    for s in rows:
        if not isinstance(s, dict):
            continue
        targets = _target_rows(s)
        mine = [t for t in targets if t["name"].upper() == wanted]
        if not mine:
            continue  # the search can return a set whose OTHER target matched
        target = mine[0]
        datasets.append({
            "urn": s.get("urn"),
            "title": s.get("title"),
            "target": target["name"],
            "targets": [t["name"] for t in targets],
            "numbering": target["numbering"],
            "target_accession": target["accession"],
            "target_sequence_type": target["sequence_type"],
            "target_sequence_length": target["sequence_length"],
            "searchable_by_hgvs": target["numbering"] == "accession",
            "num_variants": s.get("numVariants"),
            "short_description": s.get("shortDescription"),
            "published_date": s.get("publishedDate"),
            "publication": _publication(s),
            "url": score_set_url(s.get("urn")),
        })
        if target["numbering"] != "accession":
            datasets[-1]["numbering_note"] = (
                "the target is a %s sequence of %s, not an accession, so this set's HGVS is numbered from "
                "position 1 of that construct" % (target["sequence_type"] or "bare",
                                                  target["sequence_length"] or "unknown length"))
    datasets.sort(key=lambda d: (-(d.get("num_variants") or 0), str(d.get("urn"))))
    if isinstance(declared, int) and declared != len(rows):
        warnings.append("MaveDB said it had %d score set(s) for %s and returned %d; the listing may be paginated"
                        % (declared, gene.upper(), len(rows)))
    note = "%d of %d returned score sets target %s" % (len(datasets), len(rows), gene.upper())
    return Outcome(datasets, sources=[record("MaveDB score-sets/search", gene, resp, note=note)],
                   warnings=warnings)


def calibration(urn: str) -> Outcome:
    """The score set's primary calibration, or None when it publishes none."""
    url = "%s/score-calibrations/score-set/%s/primary" % (BASE, urn)
    resp = get_json(url, source="MaveDB score calibration", cache_ttl=7 * 86400, timeout=60,
                    ok_statuses=(200, 404), not_found_ttl=86400)
    src = record("MaveDB score calibration", urn, resp)
    if resp.status == 404:
        detail = ""
        try:
            body = resp.json()
            detail = str(body.get("detail") or "") if isinstance(body, dict) else ""
        except ValueError:
            detail = (resp.text or "")[:200]
        if NO_CALIBRATION_404.search(detail) or not detail:
            return Outcome(None, sources=[src])
        # a 404 about the score set or the route is a lookup failure, not "no calibration"
        raise SourceError("MaveDB score calibration", url, 404, detail or "404 with no detail")
    if not (resp.text or "").strip():
        return Outcome(None, sources=[src])
    data = validated_json(resp, "MaveDB score calibration")
    if isinstance(data, list):
        data = data[0] if data else None
    if not isinstance(data, dict):
        return Outcome(None, sources=[src])
    classes = []
    for c in data.get("functionalClassifications") or []:
        if not isinstance(c, dict):
            continue
        acmg = c.get("acmgClassification") if isinstance(c.get("acmgClassification"), dict) else {}
        raw = c.get("range")
        bounds = list(raw) if isinstance(raw, (list, tuple)) else []
        lo = bounds[0] if len(bounds) > 0 else None
        hi = bounds[1] if len(bounds) > 1 else None
        lo = lo if isinstance(lo, (int, float)) and math.isfinite(lo) else None
        hi = hi if isinstance(hi, (int, float)) and math.isfinite(hi) else None
        classes.append({
            "label": c.get("label"),
            "functional_classification": c.get("functionalClassification"),
            "range": [lo, hi],
            "range_as_given": bounds,
            "is_score_band": lo is not None or hi is not None,
            "class": c.get("class"),
            "inclusive_lower_bound": bool(c.get("inclusiveLowerBound")),
            "inclusive_upper_bound": bool(c.get("inclusiveUpperBound")),
            "acmg_criterion": acmg.get("criterion"),
            "acmg_evidence_strength": acmg.get("evidenceStrength"),
            "odds_path": c.get("oddspathsRatio"),
            "positive_likelihood_ratio": c.get("positiveLikelihoodRatio"),
            "variant_count": c.get("variantCount"),
            "description": c.get("description"),
        })
    return Outcome({
        "title": data.get("title"),
        "research_use_only": bool(data.get("researchUseOnly")),
        "baseline_score": data.get("baselineScore"),
        "baseline_score_description": data.get("baselineScoreDescription"),
        "notes": data.get("notes"),
        "classifications": classes,
        "score_bands": sum(1 for c in classes if c["is_score_band"]),
        "url": score_set_url(urn),
    }, sources=[src])


def _in_range(score: float, cls: Dict[str, Any]) -> bool:
    lo, hi = cls["range"][0], cls["range"][1]
    if lo is None and hi is None:
        return False  # not a score band at all; it cannot contain a score
    if lo is not None:
        if cls.get("inclusive_lower_bound"):
            if score < lo:
                return False
        elif score <= lo:
            return False
    if hi is not None:
        if cls.get("inclusive_upper_bound"):
            if score > hi:
                return False
        elif score >= hi:
            return False
    return True


def interpret(score: Optional[float], calib: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """The score set's own class for this score, or None when none can be read.

    None means "no class here", and the caller says why. A class is withheld -- with
    its reason in `note` -- when the score is not finite, when the calibration has no
    score band, when the score falls in none of them, and when two bands that contain
    it disagree.
    """
    if not calib:
        return None
    base = {"source": "MaveDB primary score calibration", "calibration_title": calib.get("title"),
            "research_use_only": calib.get("research_use_only"),
            "baseline_score": calib.get("baseline_score"),
            "baseline_score_description": calib.get("baseline_score_description"),
            "label": None, "functional_classification": None}
    if score is None or not isinstance(score, (int, float)) or not math.isfinite(score):
        return dict(base, note="the score is not a finite number, so it cannot be placed in this calibration")
    bands = [c for c in calib.get("classifications") or [] if c.get("is_score_band")]
    if not bands:
        categorical = [c for c in calib.get("classifications") or [] if c.get("class")]
        return dict(base, note=(
            "this calibration declares no numeric score band"
            + (" (its classifications are categorical, keyed on `class`), so a score cannot be placed in it"
               if categorical else ", so a score cannot be placed in it")))
    hits = [c for c in bands if _in_range(score, c)]
    if not hits:
        return dict(base, note="the score falls in no band of this score set's calibration")
    verdicts = {(c.get("acmg_criterion"), c.get("functional_classification")) for c in hits}
    if len(verdicts) > 1:
        return dict(base, note=(
            "the score falls in %d bands of this calibration that do not agree (%s); no class is reported"
            % (len(hits), "; ".join(str(c.get("label")) for c in hits))))
    cls = hits[0]
    out = dict(base)
    out.update({k: cls.get(k) for k in ("label", "functional_classification", "range", "acmg_criterion",
                                        "acmg_evidence_strength", "odds_path", "positive_likelihood_ratio",
                                        "variant_count")})
    if calib.get("research_use_only"):
        out["acmg_criterion_withheld"] = out.pop("acmg_criterion", None)
        out["acmg_criterion"] = None
        out["note"] = RUO_NOTE
    return out


def _header_ok(text: str) -> Optional[str]:
    """Reject a body that is not the scores CSV (MaveDB answers errors with JSON and a 200)."""
    first = (text or "").lstrip().splitlines()[0] if (text or "").strip() else ""
    if not first.startswith("accession"):
        return "not the scores CSV: the first line is %r" % first[:80]
    columns = [c.strip() for c in first.split(",")]
    if "score" not in columns:
        return "the scores CSV has no `score` column (columns: %s)" % ", ".join(columns[:8])
    return None


def _scores_page(urn: str, start: int, limit: int):
    return request("%s/score-sets/%s/scores" % (BASE, urn), source="MaveDB scores",
                   params={"start": start, "limit": limit}, accept="text/csv",
                   cache_ttl=7 * 86400, timeout=180, validate=_header_ok)


def _number(cell: Any) -> Optional[float]:
    text = str(cell if cell is not None else "").strip()
    if not text or text.upper() in ("NA", "NAN", "NULL", "NONE", ""):
        return None
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def variant_rows(urn: str, queries: Sequence[str], max_rows: int = MAX_ROWS) -> Outcome:
    """Rows of this score set whose hgvs_nt / hgvs_splice / hgvs_pro matches one of `queries`."""
    wanted = [q for q in queries if q and normalise_change(split_hgvs(q)[1])]
    found: List[Dict[str, Any]] = []
    sources: List[Dict[str, Any]] = []
    warnings: List[str] = []
    read = 0
    truncated = False
    if not wanted:
        return Outcome(found, sources=sources, warnings=warnings)
    while read < max_rows:
        page = min(PAGE_ROWS, max_rows - read)
        resp = _scores_page(urn, read, page)
        text = validated_text(resp, "MaveDB scores")
        reader = csv.DictReader(io.StringIO(text))
        rows = list(reader)
        sources.append(record("MaveDB scores", "%s rows %d..%d" % (urn, read, read + len(rows)), resp))
        if not rows:
            break
        for row in rows:
            read += 1
            for column in ("hgvs_nt", "hgvs_splice", "hgvs_pro"):
                cell = (row.get(column) or "").strip()
                if not cell or cell.upper() == "NA":
                    continue
                kind = next((match_hgvs(cell, q) for q in wanted if match_hgvs(cell, q)), None)
                if kind:
                    found.append({"accession": row.get("accession"), "matched_column": column, "hgvs": cell,
                                  "match_kind": kind, "score": _number(row.get("score")),
                                  "score_cell": row.get("score"),
                                  "hgvs_nt": row.get("hgvs_nt"), "hgvs_pro": row.get("hgvs_pro"),
                                  "hgvs_splice": row.get("hgvs_splice")})
                    break
        if len(rows) < page:
            break
        if read >= max_rows:
            truncated = True
    if truncated:
        warnings.append("MaveDB %s: stopped after reading %d rows; a match beyond that was not looked for"
                        % (urn, read))
    return Outcome(found, sources=sources, warnings=warnings,
                   query={"urn": urn, "rows_read": read, "truncated": truncated})


def lookup(gene: str, hgvs_p: Optional[str] = None, hgvs_c: Optional[str] = None,
           max_datasets: int = MAX_DATASETS_SCORED, include_construct_numbered: bool = False) -> Outcome:
    """MaveDB score sets for `gene`, and this variant's score in them when it has one.

    With no HGVS, only the score sets are listed. With `hgvs_p` and/or `hgvs_c`, each
    ACCESSION-numbered score set's table is searched for that variant and the score is
    reported with the score set's OWN calibration, or with no class at all when it has
    none. Sequence-target (construct-numbered) score sets are searched only when
    `include_construct_numbered` is set, and their hits are never given a class.
    """
    sets_out = search_score_sets(gene)
    datasets = sets_out.result
    sources = list(sets_out.sources)
    warnings = list(sets_out.warnings)
    queries = [q for q in (hgvs_p, hgvs_c) if q]
    comparable = [normalise_hgvs(q) for q in queries]
    comparable = [c for c in comparable if c]
    result: Dict[str, Any] = {
        "gene": gene.upper(),
        "query_hgvs": {"hgvs_p": hgvs_p, "hgvs_c": hgvs_c, "compared_as": comparable},
        "datasets": datasets,
        "dataset_count": len(datasets),
        "scores": [],
        "scores_searched": False,
        "numbering_rule": NUMBERING_NOTE,
        "how_to_use": BRNICH_NOTE,
    }
    if not datasets:
        result["note"] = ("MaveDB has no published score set targeting %s. That is an absence of a multiplexed "
                          "assay in MaveDB, not evidence about any variant." % gene.upper())
        return Outcome(result, sources=sources, warnings=warnings,
                       query={"gene": gene, "hgvs_p": hgvs_p, "hgvs_c": hgvs_c})
    if queries and not comparable:
        result["note"] = ("%s identifies no single variant (an HGVS such as `p.=` with no position matches any "
                          "number of rows), so no score was looked up" % " / ".join(queries))
        return Outcome(result, sources=sources, warnings=warnings,
                       query={"gene": gene, "hgvs_p": hgvs_p, "hgvs_c": hgvs_c})
    if not comparable:
        result["note"] = ("no variant given, so only the score sets are listed; pass the protein or coding HGVS to "
                          "look the variant up in them")
        return Outcome(result, sources=sources, warnings=warnings,
                       query={"gene": gene, "hgvs_p": hgvs_p, "hgvs_c": hgvs_c})

    eligible = [d for d in datasets if d["searchable_by_hgvs"] or include_construct_numbered]
    skipped = [d for d in datasets if d not in eligible]
    searched = eligible[:max(1, int(max_datasets))]
    result["scores_searched"] = bool(searched)
    result["datasets_searched"] = [d["urn"] for d in searched]
    result["datasets_not_searched"] = [
        {"urn": d["urn"], "reason": "construct-numbered target (%s); its HGVS is not the gene's numbering. "
                                    "Pass include_construct_numbered to search it anyway, knowing a hit names the "
                                    "construct's residue, not the gene's." % d.get("numbering_note", "no accession")}
        for d in skipped]
    if len(eligible) > len(searched):
        warnings.append("MaveDB: %d searchable score set(s) target %s; the variant was looked for in the %d largest "
                        "(%s)" % (len(eligible), gene.upper(), len(searched),
                                  ", ".join(str(d["urn"]) for d in searched)))
    if skipped:
        warnings.append("MaveDB: %d score set(s) for %s are construct-numbered and were not searched by HGVS (%s)"
                        % (len(skipped), gene.upper(), ", ".join(str(d["urn"]) for d in skipped)))

    scores: List[Dict[str, Any]] = []
    incomplete: List[str] = []
    for d in searched:
        urn = str(d["urn"])
        rows_out = attempt("MaveDB scores %s" % urn, lambda u=urn: variant_rows(u, queries), warnings)
        if rows_out is None:
            incomplete.append(urn)
            continue
        sources.extend(rows_out.sources)
        warnings.extend(rows_out.warnings)
        if (rows_out.query or {}).get("truncated"):
            incomplete.append(urn)
        if not rows_out.result:
            continue
        calib_out = attempt("MaveDB calibration %s" % urn, lambda u=urn: calibration(u), warnings)
        calib = None
        if calib_out is None:
            incomplete.append(urn)
        else:
            sources.extend(calib_out.sources)
            calib = calib_out.result
        for row in rows_out.result:
            entry: Dict[str, Any] = {
                "urn": urn, "title": d.get("title"), "hgvs": row["hgvs"], "score": row["score"],
                "matched_column": row["matched_column"], "variant_accession": row.get("accession"),
                "numbering": d["numbering"], "target_accession": d.get("target_accession"),
                "url": ("https://mavedb.org/variants/%s" % row["accession"]) if row.get("accession") else d["url"],
                "interpretation": None,
            }
            if row["match_kind"] == "change_only":
                entry["accession_checked"] = False
                entry["match_note"] = (
                    "matched on the change alone: one side carried no transcript or protein accession. Coding and "
                    "protein numbering are accession-specific, so confirm this is the same variant before using it."
                    + (" The score set's target is %s." % d["target_accession"] if d.get("target_accession") else ""))
                warnings.append("MaveDB %s: %s matched %s without an accession on both sides"
                                % (urn, row["hgvs"], " / ".join(queries)))
            else:
                entry["accession_checked"] = True
            if d["numbering"] != "accession":
                entry["interpretation_note"] = (
                    "no functional class: this score set's target is a construct, so position %s is the "
                    "construct's residue, not the gene's" % row["hgvs"])
            elif row["score"] is None:
                entry["interpretation_note"] = NO_SCORE_NOTE
            elif calib is None:
                entry["interpretation_note"] = NO_CALIBRATION_NOTE
            else:
                interp = interpret(row["score"], calib)
                entry["interpretation"] = interp
                if interp is not None and not interp.get("label"):
                    entry["interpretation_note"] = "no functional class: %s" % interp.get("note")
            scores.append(entry)
    result["scores"] = scores
    result["score_count"] = len(scores)
    result["datasets_incomplete"] = sorted(set(incomplete))
    if not scores:
        complete = [d["urn"] for d in searched if d["urn"] not in incomplete]
        if complete:
            result["note"] = ("the variant was not found in the %d score set(s) searched in full (%s). A MAVE "
                              "usually covers only part of a gene, so this means 'not assayed there', not 'no "
                              "effect'." % (len(complete), ", ".join(str(u) for u in complete)))
        else:
            result["note"] = ("no score set was searched in full, so nothing can be said about this variant here; "
                              "see the warnings")
        if incomplete:
            result["note"] += (" %d score set(s) were searched only partly or failed (%s) and are not covered by "
                               "that statement." % (len(set(incomplete)), ", ".join(sorted(set(incomplete)))))
    return Outcome(result, sources=sources, warnings=warnings,
                   query={"gene": gene, "hgvs_p": hgvs_p, "hgvs_c": hgvs_c,
                          "include_construct_numbered": bool(include_construct_numbered)})


def render(result: Dict[str, Any]) -> str:
    lines = ["MaveDB %s: %d published score set(s) targeting this gene"
             % (result["gene"], result["dataset_count"])]
    for d in result.get("datasets") or []:
        pub = d.get("publication") or {}
        tail = "  PMID %s" % pub["pmid"] if pub.get("pmid") else ""
        mark = "" if d.get("searchable_by_hgvs") else "  [construct-numbered: not searched by HGVS]"
        lines.append("  %-26s %6s variants  %s%s%s"
                     % (str(d.get("urn") or "-"), str(d.get("num_variants") or "-"),
                        str(d.get("title"))[:60], tail, mark))
    if result.get("scores"):
        lines.append("")
        lines.append("this variant:")
        for s in result["scores"]:
            interp = s.get("interpretation") or {}
            if interp.get("label"):
                verdict = "%s (%s) per %s" % (interp["label"], interp.get("functional_classification"),
                                              interp.get("calibration_title"))
                if interp.get("acmg_criterion"):
                    verdict += " -> %s %s" % (interp["acmg_criterion"], interp.get("acmg_evidence_strength"))
                elif interp.get("acmg_criterion_withheld"):
                    verdict += (" -> the dataset labels it %s, withheld: %s"
                                % (interp["acmg_criterion_withheld"], interp.get("note")))
            else:
                verdict = "no functional class: %s" % (s.get("interpretation_note") or interp.get("note")
                                                       or "no calibration")
            score = s.get("score")
            lines.append("  %s %s  score=%s  %s" % (s["urn"], s["hgvs"],
                                                    "-" if score is None else round(score, 4), verdict))
            if s.get("match_note"):
                lines.append("      ! %s" % s["match_note"])
    elif result.get("note"):
        lines.append("")
        lines.append(result["note"])
    for item in result.get("datasets_not_searched") or []:
        lines.append("  not searched: %s - %s" % (item["urn"], item["reason"]))
    lines.append("")
    lines.append("! " + result["numbering_rule"])
    lines.append("! " + result["how_to_use"])
    return "\n".join(lines)
