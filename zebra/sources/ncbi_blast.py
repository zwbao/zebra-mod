"""NCBI BLAST URL API: how many places in the human genome a short sequence matches.

Service: https://blast.ncbi.nlm.nih.gov/Blast.cgi (no key). It is asynchronous --
`CMD=Put` returns a request id (RID), `CMD=Get&FORMAT_OBJECT=SearchInfo` reports
`Status=WAITING|READY|FAILED|UNKNOWN`, and only then does the report come back.
Asking for the report early returns an HTML waiting page with HTTP 200, which is why
the status is always checked first.

Why this service and not UCSC BLAT: every UCSC BLAT host (genome, genome-euro,
genome-asia) answers `hgBlat?...&output=json` with a Cloudflare Turnstile challenge
page instead of a result, so it cannot be scripted (checked live 2026-10-06 on all
three). NCBI's URL API has no such gate.

Report format: `FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular`, one line per alignment.
`FORMAT_TYPE=Tabular` on its own is accepted and returns only a 62-byte QBlastInfo
stub, and `FORMAT_TYPE=CSV` returns an empty body (both measured 2026-10-06), so
neither is used. The JSON report is not used either: it is unbounded in practice -- a
single 35 nt repeat-derived query with the low-complexity filter off produced a JSON
body still growing past 120 MB when it was cut off, which `zebra.http` would have read
entirely into memory. The tabular body for the same search was 21 MB, and most of that
was the one repeat-derived query, so `HITLIST_SIZE` is kept small to bound it. The
low-complexity (DUST) filter is left ON for the same reason; masking can no longer fake
an answer, because a query that comes back with no locus at all -- not even its own --
is reported as VOID rather than as a count (see `_void_if_self_missing`).

The tabular body arrives wrapped in the service's `<p><!-- QBlastInfo ... --><p><PRE>`
preamble, carries `# Query:`, `# Fields:` and `# N hits found` comment lines, and names
each query by its FASTA title in the first column. Verified live on 2026-10-06: a 25 nt
query taken from GRCh38 7:117,639,947-117,639,971 came back with a 25/25 identity
alignment at exactly `NC_000007.14 117639947 117639971`, alongside a T2T-CHM13 contig
(`NC_060931.1`) that this module excludes.

NCBI's usage policy for this API: one search submitted at a time and no more than one
every 10 seconds, no polling of a single RID more often than once a minute, and
batches of more than ~100 searches outside US business hours. Every call here sends
ALL the sequences of one screen as a single multi-FASTA search, so one screen costs
one search; the intervals are enforced below.

Uniqueness is read conservatively: only a full-length alignment (at most one base
short of the query), with at most one mismatch, no gap, and on a GRCh38 primary
chromosome counts as a locus. Hits on alternate haplotypes, patch scaffolds and the
T2T-CHM13 assembly -- the search database holds all of them -- are counted separately
and never make a sequence look non-unique.
"""

from __future__ import annotations

import os
import re
import time
import urllib.parse
from typing import Any, Dict, List, Mapping, Optional, Tuple

from zebra.core import Outcome
from zebra.http import SourceError, deadline_seconds, request
from zebra.sources import record

BASE = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
# The Genome Data Hub set for Homo sapiens: GRCh38 primary chromosomes plus alternate
# haplotypes, patch scaffolds and T2T-CHM13. The GRCh38 primary chromosomes are the
# RefSeq accessions NC_000001..NC_000024; chrM is NC_012920.
DATABASE = "GPIPE/9606/current/all_top_level"
PRIMARY_RE = re.compile(r"^NC_0000(0[1-9]|1[0-9]|2[0-4])(\.\d+)?$")
CHRM_RE = re.compile(r"^NC_012920(\.\d+)?$")
CHROM_OF = {"NC_0000%02d" % i: (str(i) if i <= 22 else {23: "X", 24: "Y"}[i]) for i in range(1, 25)}

RID_RE = re.compile(r"RID\s*=\s*([A-Z0-9]+)", re.I)
STATUS_RE = re.compile(r"Status\s*=\s*([A-Z]+)", re.I)
QUERY_LINE_RE = re.compile(r"^#\s*Query:\s*(.*)$")
FIELDS_LINE_RE = re.compile(r"^#\s*Fields:\s*(.*)$")

MAX_QUERIES = 50          # one multi-FASTA search; more than this is a different job
MAX_QUERY_LEN = 200       # this is a short-oligonucleotide check, not a sequence search
# Bounds the subject count, and with it the report size. 20 is ample for a uniqueness
# question: a target with full-length matches at more than 20 places is not unique by
# any reading, and the target's own locus is a perfect match so it is never ranked out.
HITLIST_SIZE = 20
MAX_REPORT_LINES = 400000  # a bound on parsing work even if the service ignores HITLIST_SIZE

def _interval(name: str, default: float, floor: float) -> float:
    """An env override, clamped to the floor NCBI's usage policy sets.

    The intervals are a condition of using this service, so an environment variable may
    lengthen them and may not shorten them below the published limit.
    """
    try:
        value = float(os.environ.get(name, default))
    except (TypeError, ValueError):
        return default
    return max(floor, value)


SUBMIT_INTERVAL = _interval("ZEBRA_BLAST_SUBMIT_INTERVAL", 10.0, 10.0)
FIRST_POLL_AFTER = _interval("ZEBRA_BLAST_FIRST_POLL", 20.0, 5.0)
POLL_INTERVAL = _interval("ZEBRA_BLAST_POLL_INTERVAL", 60.0, 60.0)
MAX_WAIT = _interval("ZEBRA_BLAST_MAX_WAIT", 240.0, 0.0)
# A submitted search is cached for a few hours, so a rerun of the same screen reuses the
# RID instead of queueing a second search. The queue's latency is not predictable --
# the same four-sequence search came back in under a minute once and was still running
# after thirteen minutes two hours later (both measured 2026-10-06, US business hours) --
# so a run that gives up waiting says which RID to collect, and the rerun collects it.
# NCBI keeps a RID's result for about a day.
SUBMIT_TTL = 6 * 3600

TERMS = ("NCBI BLAST URL API: one search at a time, at most one every 10 s, one poll per RID per minute; "
         "large batches outside US business hours. https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html")
UNIQUENESS_CRITERION = (
    "a locus is an alignment covering the whole query (at most 1 base short), with at most 1 mismatch and no gap, "
    "on a GRCh38 primary chromosome (RefSeq NC_000001-NC_000024, NC_012920 for chrM). Alternate haplotypes, patch "
    "scaffolds and T2T-CHM13 contigs in the same database are counted separately and do not count against "
    "uniqueness. The search runs with NCBI's low-complexity filter on, so a repeat-derived target can come back "
    "with no alignment at all; that is reported as a void check, never as a count."
)

_last_submit = 0.0


def _pace_submit() -> None:
    """Wait out the 10 s NCBI asks for between submissions, unless there is no time left."""
    global _last_submit
    wait = _last_submit + SUBMIT_INTERVAL - time.monotonic()
    if wait <= 0:
        return
    left = deadline_seconds()
    if left is not None and left - wait <= 0.5:
        raise SourceError("NCBI BLAST (submit)", BASE, None,
                          "the %g s courtesy interval between searches would use the whole remaining time "
                          "budget; no search was submitted" % SUBMIT_INTERVAL)
    time.sleep(wait)


def fasta(queries: Mapping[str, str]) -> str:
    """A multi-FASTA body; the names become the BLAST query titles."""
    if len(queries) > MAX_QUERIES:
        raise ValueError(f"{len(queries)} sequences in one search; this module sends at most {MAX_QUERIES}")
    out = []
    for name, seq in queries.items():
        if not re.match(r"^[A-Za-z0-9_.:-]{1,60}$", str(name)):
            raise ValueError(f"{name!r} is not usable as a FASTA title (letters, digits, _ . : - only)")
        clean = re.sub(r"[^ACGTNacgtn]", "", seq or "").upper()
        if not clean:
            raise ValueError(f"{name}: no DNA bases to search")
        if len(clean) > MAX_QUERY_LEN:
            raise ValueError(f"{name}: {len(clean)} bases; this check takes sequences up to {MAX_QUERY_LEN}")
        out.append(">%s\n%s" % (name, clean))
    if not out:
        raise ValueError("no sequences to search")
    return "\n".join(out) + "\n"


def _form(params: Dict[str, Any]) -> str:
    return urllib.parse.urlencode({k: v for k, v in params.items() if v is not None})


def submit(queries: Mapping[str, str], *, refresh: bool = False) -> Tuple[str, Dict[str, Any]]:
    """Queue one multi-FASTA megablast search; returns (RID, source record)."""
    body = _form({
        "CMD": "Put",
        "PROGRAM": "blastn",
        "MEGABLAST": "on",
        "DATABASE": DATABASE,
        "QUERY": fasta(queries),
        "HITLIST_SIZE": HITLIST_SIZE,
        "WORD_SIZE": 16,
        # an 18-25 nt exact match scores ~36-50 bits, which is not significant at the
        # default E=10 in a 3 Gb database: a tight expect value would silently drop
        # exactly the alignments a uniqueness check is looking for
        "EXPECT": 1000,
        "FORMAT_TYPE": "Tabular",
        "tool": "zebra-mod",
    })
    _pace_submit()
    resp = request(BASE, source="NCBI BLAST (submit)", method="POST", body=body,
                   headers={"Content-Type": "application/x-www-form-urlencoded"},
                   accept="text/html", cache_ttl=SUBMIT_TTL, timeout=120, retries=1, refresh=refresh,
                   validate=lambda t: None if RID_RE.search(t or "") else
                   "the submit response carries no RID (NCBI refused or rate-limited the search)")
    global _last_submit
    if not resp.cached:
        _last_submit = time.monotonic()
    m = RID_RE.search(resp.text)
    if not m:  # validate() already guarantees this; kept so a contract change cannot pass silently
        raise SourceError("NCBI BLAST (submit)", BASE, resp.status, "no RID in the submit response")
    rid = m.group(1)
    return rid, record("NCBI BLAST", rid, resp, url="%s?CMD=Get&RID=%s" % (BASE, rid), note=TERMS)


def status(rid: str) -> str:
    """WAITING, READY, FAILED or UNKNOWN (an expired or unknown RID)."""
    resp = request("%s?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=%s" % (BASE, urllib.parse.quote(rid)),
                   source="NCBI BLAST (status)", accept="text/html", cache_ttl=0, timeout=60, retries=1,
                   validate=lambda t: None if STATUS_RE.search(t or "") else "no QBlastInfo Status in the response")
    m = STATUS_RE.search(resp.text)
    return m.group(1).upper() if m else "UNKNOWN"


def _is_tabular(text: Optional[str]) -> Optional[str]:
    """None when this body is a tabular report; otherwise why it is not.

    The service answers a report request that it cannot serve with a short plain
    message ("SYSTEM CAN\'T PROCESS YOUR REQUEST", seen after repeated large
    downloads of one RID) and a report request made too early with a QBlastInfo stub.
    Neither contains a field list or a hit count, and neither may be parsed as "no
    alignments" -- which is why this is a `validate` callback: a body it rejects is
    never cached and never returned.
    """
    body = text or ""
    if "hits found" in body or "# Fields:" in body:
        return None
    return ("not a BLAST tabular report: no field list and no hit count "
            "(the search may be unfinished, or the service declined to serve it)")


def report(rid: str) -> Tuple[str, Dict[str, Any]]:
    """The finished tabular report. Only call this once the status is READY."""
    url = "%s?CMD=Get&RID=%s&FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular" % (BASE, urllib.parse.quote(rid))
    resp = request(url, source="NCBI BLAST (report)", accept="text/plain", cache_ttl=SUBMIT_TTL,
                   timeout=180, retries=1, validate=_is_tabular)
    return resp.text, record("NCBI BLAST", rid, resp, note=TERMS)


# The default tabular field order, used only when the report carries no `# Fields:` line.
DEFAULT_FIELDS = ("query acc.ver", "subject acc.ver", "% identity", "alignment length", "mismatches",
                  "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")
WANTED = {
    "subject": ("subject acc.ver", "subject id", "subject accession.version", "subject accession"),
    "identity": ("% identity",),
    "align_len": ("alignment length",),
    "mismatches": ("mismatches",),
    "gaps": ("gap opens", "gaps"),
    "s_start": ("s. start",),
    "s_end": ("s. end",),
}


def parse_tabular(text: str) -> Tuple[Dict[str, List[Dict[str, Any]]], List[str], bool]:
    """Tabular report -> {query title: [alignment rows]}, the titles in order, truncated?

    The column order is read from the report's own `# Fields:` line rather than
    assumed, so a change in the service's default field set cannot silently shift
    every number by one column.
    """
    blocks: Dict[str, List[Dict[str, Any]]] = {}
    order: List[str] = []
    fields = list(DEFAULT_FIELDS)
    current: Optional[str] = None
    lines = 0
    truncated = False
    for raw in (text or "").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            m = QUERY_LINE_RE.match(line)
            if m:
                current = m.group(1).strip().split()[0] if m.group(1).strip() else None
                if current is not None and current not in blocks:
                    blocks[current] = []
                    order.append(current)
                continue
            m = FIELDS_LINE_RE.match(line)
            if m:
                fields = [f.strip() for f in m.group(1).split(",")]
            continue
        lines += 1
        if lines > MAX_REPORT_LINES:
            truncated = True
            break
        if current is None:
            continue
        cells = line.split("\t")
        if len(cells) < len(fields):
            continue
        row: Dict[str, Any] = {}
        ok = True
        for key, names in WANTED.items():
            index = next((fields.index(n) for n in names if n in fields), None)
            if index is None or index >= len(cells):
                ok = False
                break
            row[key] = cells[index]
        if not ok:
            continue
        try:
            row["identity"] = float(row["identity"])
            row["align_len"] = int(row["align_len"])
            row["mismatches"] = int(row["mismatches"])
            row["gaps"] = int(row["gaps"])
            row["s_start"] = int(row["s_start"])
            row["s_end"] = int(row["s_end"])
        except (TypeError, ValueError):
            continue
        blocks[current].append(row)
    return blocks, order, truncated


def _chrom_label(accession: str) -> str:
    bare = accession.split(".")[0]
    if CHRM_RE.match(accession):
        return "MT"
    return CHROM_OF.get(bare, accession)


def _loci(rows: List[Dict[str, Any]], query_len: int) -> Dict[str, Any]:
    """Full-length near-exact loci for one query, split by primary / other contigs."""
    primary: List[Dict[str, Any]] = []
    other = 0
    for row in rows:
        accession = str(row.get("subject") or "")
        # some field sets spell the subject as gi|...|ref|NC_000007.14|
        if "|" in accession:
            parts = [p for p in accession.split("|") if p]
            accession = parts[-1] if parts else accession
        if row["align_len"] < query_len - 1 or row["mismatches"] > 1 or row["gaps"]:
            continue
        if not (PRIMARY_RE.match(accession) or CHRM_RE.match(accession)):
            other += 1
            continue
        start, end = row["s_start"], row["s_end"]
        primary.append({
            "accession": accession,
            "chrom": _chrom_label(accession),
            "start": min(start, end),
            "end": max(start, end),
            "strand": "-" if start > end else "+",
            "identity_pct": row["identity"],
            "align_len": row["align_len"],
            "mismatches": row["mismatches"],
        })
    # The same locus can come back as two alignments. Collapse by OVERLAP: two
    # alignments are one locus when they are on the same accession and within a query
    # length of each other. A fixed-width bucket is not a distance test -- it merges two
    # genuinely distinct loci that fall in one bucket (a false "unique") and splits one
    # locus that straddles a boundary (a false "not unique").
    unique_loci: List[Dict[str, Any]] = []
    for locus in sorted(primary, key=lambda l: (str(l["accession"]), l["start"])):
        previous = unique_loci[-1] if unique_loci else None
        if (previous is not None and previous["accession"] == locus["accession"]
                and locus["start"] - previous["start"] <= max(1, query_len)):
            previous["end"] = max(previous["end"], locus["end"])
            previous["alignments"] = previous.get("alignments", 1) + 1
            continue
        unique_loci.append(dict(locus, alignments=1))
    out = {
        "query_len": query_len,
        "loci": unique_loci[:20],
        "locus_count": len(unique_loci),
        "unique": len(unique_loci) == 1,
        "non_primary_alignments": other,
        "criterion": UNIQUENESS_CRITERION,
    }
    subjects = len({str(row.get("subject") or "") for row in rows})
    if subjects >= HITLIST_SIZE or len(unique_loci) >= HITLIST_SIZE:
        # The hit list caps SUBJECTS, so a saturated subject list means further loci may
        # exist that were never reported -- including primary-chromosome ones crowded out
        # by alternate haplotypes and patches. "unique" cannot be asserted from that.
        out["locus_count_is_floor"] = (
            "the search returns at most %d subject sequences and %d came back, so the count above is a floor: "
            "further loci may exist that were not reported." % (HITLIST_SIZE, subjects))
        out["unique"] = None
        out["unique_note"] = ("not determined: the subject list was saturated, so a target that looks unique here "
                              "may not be")
    return out


def _void_if_self_missing(row: Dict[str, Any], warnings: List[str], name: str) -> Dict[str, Any]:
    """A sequence taken from the genome must match the genome at least once.

    Zero loci therefore means the search did not see the target's own locus -- the
    low-complexity mask, a reporting threshold, or a truncated hit list -- and the
    answer is void. Reporting it as "0 loci" would read as a result; it is not one.
    """
    if row.get("checked") and row.get("locus_count") == 0:
        warnings.append("%s: BLAST reported no full-length alignment anywhere, not even the target's own locus; "
                        "the uniqueness check is void for it" % name)
        return {"checked": False, "rid": row.get("rid"), "url": row.get("url"),
                "reason": "BLAST returned no full-length alignment for this target, not even its own locus in the "
                          "genome. A sequence read off the genome must match it at least once, so this answer is "
                          "void rather than a count of zero (a repeat-derived or low-complexity target is the "
                          "usual cause).",
                "non_primary_alignments": row.get("non_primary_alignments"),
                "criterion": row.get("criterion")}
    return row


def locate(queries: Mapping[str, str], max_wait: Optional[float] = None) -> Outcome:
    """Where each sequence matches the human genome: one BLAST search for all of them.

    Returns `result = {name: {"checked": True, "loci": [...], "locus_count": n,
    "unique": bool, ...}}`, or `{"checked": False, "reason": ...}` for a name whose
    answer could not be obtained or is void -- never a silent "unique".
    """
    budget = MAX_WAIT if max_wait is None else float(max_wait)
    lengths = {name: len(re.sub(r"[^ACGTNacgtn]", "", seq or "")) for name, seq in queries.items()}
    names = list(queries)
    sources: List[Dict[str, Any]] = []
    warnings: List[str] = []
    rid, src = submit(queries)
    sources.append(src)
    deadline = time.monotonic() + budget
    state = "WAITING"
    reason: Optional[str] = None
    first = True
    resubmitted = False
    while True:
        left_budget = deadline - time.monotonic()
        left_deadline = deadline_seconds()
        pause = FIRST_POLL_AFTER if first else POLL_INTERVAL
        if left_budget < pause or (left_deadline is not None and left_deadline - pause <= 5.0):
            reason = ("the BLAST search %s was still running after %.0fs and the time budget for this call ran "
                      "out. Rerun the same command to collect it: the submitted search is cached for %dh, so the "
                      "rerun picks up this RID instead of queueing another. It is also readable at "
                      "%s?CMD=Get&RID=%s for about a day."
                      % (rid, budget - max(0.0, left_budget), SUBMIT_TTL // 3600, BASE, rid))
            break
        time.sleep(pause)
        first = False
        state = status(rid)
        if state == "READY":
            break
        if state == "FAILED":
            reason = "NCBI reported the BLAST search %s as FAILED" % rid
            break
        if state == "UNKNOWN":
            if not resubmitted:
                # the RID came out of this module's own 6 h submit cache and has since
                # expired upstream; telling the caller to "rerun" would hand them the
                # same cached id for ever, so re-submit past the cache once
                resubmitted = True
                warnings.append("the BLAST request id %s had expired; a new search was submitted" % rid)
                rid, src = submit(queries, refresh=True)
                sources.append(src)
                first = True
                continue
            reason = ("NCBI does not know the BLAST request id %s even after re-submitting; the service is not "
                      "answering usably" % rid)
            break
    result: Dict[str, Any] = {}
    if state != "READY":
        for name in names:
            result[name] = {"checked": False, "reason": reason or "BLAST status %s" % state, "rid": rid}
        warnings.append("genomic uniqueness not checked: %s" % (reason or state))
        return Outcome(result, sources=sources, warnings=warnings)
    text, rep_src = report(rid)
    sources.append(rep_src)
    blocks, order, truncated = parse_tabular(text)
    if truncated:
        warnings.append("the BLAST report was longer than %d alignment lines and was not read past that; "
                        "uniqueness counts below may be low" % MAX_REPORT_LINES)
    for i, name in enumerate(names):
        rows = blocks.get(name)
        if rows is None and i < len(order):   # the title did not survive the round trip
            rows = blocks.get(order[i])
        if rows is None:
            result[name] = {"checked": False, "rid": rid, "url": "%s?CMD=Get&RID=%s" % (BASE, rid),
                            "reason": "the BLAST report carries no block for this sequence"}
            warnings.append("%s: the BLAST report carries no block for this sequence" % name)
            continue
        row = _loci(rows, lengths.get(name, 0))
        row.update(checked=True, rid=rid, url="%s?CMD=Get&RID=%s" % (BASE, rid), query_title=name)
        result[name] = _void_if_self_missing(row, warnings, name)
    if len(order) != len(names):
        warnings.append("BLAST reported %d query block(s) for %d sequence(s)" % (len(order), len(names)))
    return Outcome(result, sources=sources, warnings=warnings)
