"""MyVariant.info (BioThings): population frequency, predictors and ClinVar for many variants at once.

    POST https://myvariant.info/v1/variant
    {"ids": [<= 1000 genomic HGVS ids], "fields": "...", "assembly": "hg38" | "hg19"}

Used by `zebra vcf triage --prefilter myvariant` to drop common alleles before
anything goes to Ensembl VEP, so a whole exome can be triaged instead of a
capped subset. What is sent is the variant key only (chrom, position, REF,
ALT spelled as an HGVS g. id): no sample names, genotypes or depths.

Verified live on 2026-10-06 (field paths, single-element collapsing, `notfound`):
- ids: `chr7:g.117559479G>A` (SNV), `chr7:g.117559591_117559593del` (the
  deleted bases, anchor dropped: VCF 117559590 ATCT>A), `chr2:g.176093093_176093094insA`
  (VCF 176093093 C>CA), `chr1:g.12282655del` (one base); delins as
  `chrN:g.S_EdelinsALT`. The docs cap a request at 1000 ids ("the rest will be
  omitted"), so every answer is matched back by its `query` field and a
  missing one is reported, never assumed absent.
- `gnomad_exome.af.*` is gnomAD **v2.1.1** exomes (on hg38 a liftover), and
  `gnomad_genome.af.*` is gnomAD v3 genomes on hg38 (it has ami/mid groups),
  v2.1.1 genomes on hg19. Neither carries a grpmax/popmax field: the caller
  computes a grpmax-like maximum over afr/amr/eas/nfe/sas exactly as for VEP.
  The versions are read from /v1/metadata at run time, not assumed.
- `dbnsfp.revel.score`, `dbnsfp.alphamissense.score|pred` are lists (one per
  transcript); `dbnsfp.cadd.phred` on hg38, top-level `cadd.phred` on hg19.
  BioThings collapses a one-element list to its element (`snpeff.ann` and
  `clinvar.rcv` arrive as a dict when there is one), so every list field is
  normalised here.
- A variant MyVariant does not hold comes back `{"query": id, "notfound": true}`:
  that means "not in MyVariant" (no frequency known), never "absent from gnomAD".
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from zebra.core import Outcome
from zebra.http import post_json, get_json
from zebra.sources import attempt, record, validated_json

BASE = "https://myvariant.info/v1"
MAX_IDS = 1000  # the service's documented per-request maximum
CACHE_TTL = 30 * 86400
ASSEMBLY = {"GRCh38": "hg38", "GRCh37": "hg19"}
POPS = ("afr", "amr", "eas", "nfe", "sas", "asj", "fin", "oth", "ami", "mid")
FIELDS = ",".join(
    [f"gnomad_exome.af.{k}" for k in ("af",) + tuple(f"af_{p}" for p in POPS)]
    + [f"gnomad_genome.af.{k}" for k in ("af",) + tuple(f"af_{p}" for p in POPS)]
    + ["gnomad_exome.an.an", "gnomad_genome.an.an",
       "clinvar.rcv.clinical_significance", "clinvar.rcv.review_status", "clinvar.variant_id",
       "dbnsfp.revel.score", "dbnsfp.alphamissense.score", "dbnsfp.alphamissense.pred", "dbnsfp.cadd.phred",
       "cadd.phred", "snpeff.ann.putative_impact", "snpeff.ann.effect", "snpeff.ann.gene_id", "dbsnp.rsid"]
)
IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
TERMS = ("MyVariant.info (BioThings, Su lab, Scripps): aggregates gnomAD (ODbL), dbNSFP, ClinVar, dbSNP and "
         "snpEff; free, no key; batch queries of up to 1000 ids")


def hgvs_id(chrom: str, pos: int, ref: str, alt: str) -> Optional[str]:
    """MyVariant's genomic HGVS id for a VCF allele already trimmed to a minimal form (one anchor base kept).

    SNV `chrN:g.POSR>A`; deletion `chrN:g.S_Edel` (S..E the deleted bases);
    insertion `chrN:g.P_P+1insSEQ`; anything else `chrN:g.S_EdelinsSEQ`.
    None for alleles that have no such spelling (symbolic, empty).
    """
    ref, alt = ref.upper(), alt.upper()
    if not ref or not alt or not ref.isalpha() or not alt.isalpha() or ref == alt:
        return None
    c = "MT" if chrom.upper() in ("M", "MT") else chrom
    head = f"chr{c}:g."
    if len(ref) == 1 and len(alt) == 1:
        return f"{head}{pos}{ref}>{alt}"
    if len(alt) == 1 and ref[0] == alt[0]:  # deletion after an anchor base
        s, e = pos + 1, pos + len(ref) - 1
        return f"{head}{s}del" if s == e else f"{head}{s}_{e}del"
    if len(ref) == 1 and alt[0] == ref[0]:  # insertion after an anchor base
        return f"{head}{pos}_{pos + 1}ins{alt[1:]}"
    e = pos + len(ref) - 1
    return f"{head}{pos}delins{alt}" if e == pos else f"{head}{pos}_{e}delins{alt}"


def _as_list(x: Any) -> List[Any]:
    if x is None:
        return []
    return list(x) if isinstance(x, list) else [x]


def _max_num(x: Any) -> Optional[float]:
    vals = []
    for v in _as_list(x):
        try:
            vals.append(float(v))
        except (TypeError, ValueError):
            continue
    return max(vals) if vals else None


def parse(hit: Dict[str, Any]) -> Dict[str, Any]:
    """The fields triage reads from one MyVariant hit; frequency groups keyed the way Ensembl VEP keys them.

    `groups` uses VEP's names (`gnomade`, `gnomade_afr`, `gnomadg`, …) so
    `zebra.vcf.af_summary` reads both sources with one set of rules.
    """
    groups: Dict[str, float] = {}
    an: Dict[str, Optional[float]] = {}
    for src, prefix in (("gnomad_exome", "gnomade"), ("gnomad_genome", "gnomadg")):
        block = hit.get(src)
        if not isinstance(block, dict):
            continue
        af = block.get("af") if isinstance(block.get("af"), dict) else {}
        for key, val in af.items():
            num = _max_num(val)
            if num is None:
                continue
            if key == "af":
                groups[prefix] = num
            elif key.startswith("af_") and key[3:] in POPS:
                groups[f"{prefix}_{key[3:]}"] = num
        an_block = block.get("an") if isinstance(block.get("an"), dict) else {}
        an[prefix] = _max_num(an_block.get("an"))
    clin: List[str] = []
    review: List[str] = []
    cv = hit.get("clinvar") if isinstance(hit.get("clinvar"), dict) else {}
    for rcv in _as_list(cv.get("rcv")):
        if not isinstance(rcv, dict):
            continue
        for sig in _as_list(rcv.get("clinical_significance")):
            s = str(sig).strip()
            if s and s not in clin:
                clin.append(s)
        rs = rcv.get("review_status")
        if rs and str(rs) not in review:
            review.append(str(rs))
    db = hit.get("dbnsfp") if isinstance(hit.get("dbnsfp"), dict) else {}
    am = db.get("alphamissense") if isinstance(db.get("alphamissense"), dict) else {}
    am_pred = [str(p) for p in _as_list(am.get("pred"))]
    cadd = _max_num((db.get("cadd") or {}).get("phred") if isinstance(db.get("cadd"), dict) else None)
    if cadd is None and isinstance(hit.get("cadd"), dict):
        cadd = _max_num(hit["cadd"].get("phred"))
    impact: Optional[str] = None
    effects: List[str] = []
    genes: List[str] = []
    sn = hit.get("snpeff") if isinstance(hit.get("snpeff"), dict) else {}
    for ann in _as_list(sn.get("ann")):
        if not isinstance(ann, dict):
            continue
        imp = str(ann.get("putative_impact") or "").upper() or None
        if imp in IMPACT_RANK and IMPACT_RANK[imp] > IMPACT_RANK.get(impact or "", -1):
            impact = imp
        for eff in str(ann.get("effect") or "").split("&"):
            if eff and eff not in effects:
                effects.append(eff)
        g = ann.get("gene_id")
        if g and str(g) not in genes:
            genes.append(str(g))
    rsid = (hit.get("dbsnp") or {}).get("rsid") if isinstance(hit.get("dbsnp"), dict) else None
    return {
        "found": True, "id": hit.get("_id"), "groups": groups, "an": an,
        "clinvar": clin, "clinvar_review": review, "clinvar_variation": cv.get("variant_id"),
        "revel": _max_num(db.get("revel", {}).get("score") if isinstance(db.get("revel"), dict) else None),
        "alphamissense": _max_num(am.get("score")),
        "am_class": ("likely_pathogenic" if "P" in am_pred else ("likely_benign" if am_pred and set(am_pred) <= {"B"}
                                                                 else ("ambiguous" if am_pred else None))),
        "cadd": cadd, "impact": impact, "effects": effects[:8], "genes": genes[:5], "rsid": rsid,
    }


def metadata(assembly: str = "GRCh38") -> Outcome:
    """Source versions MyVariant serves for this build (gnomAD, dbNSFP, ClinVar, build date)."""
    hg = ASSEMBLY.get(assembly)
    if hg is None:
        raise ValueError("assembly must be GRCh38 or GRCh37")
    resp = get_json(f"{BASE}/metadata", source="MyVariant.info", params={"assembly": hg}, cache_ttl=7 * 86400,
                    timeout=60)
    data = validated_json(resp, "MyVariant.info", require="src")
    src = data.get("src") or {}
    versions = {k: (src.get(k) or {}).get("version") for k in ("gnomad", "dbnsfp", "clinvar", "dbsnp", "snpeff")}
    out = {"assembly": hg, "build_date": data.get("build_date"), "build_version": data.get("build_version"),
           "versions": versions}
    return Outcome(out, sources=[record("MyVariant.info", f"metadata {hg}", resp,
                                        note="source versions served by MyVariant.info")])


def _post(ids: Sequence[str], hg: str, ttl: float = CACHE_TTL):
    payload = {"ids": list(ids), "fields": FIELDS, "assembly": hg}

    def go(t: float):
        return post_json(f"{BASE}/variant", payload, source="MyVariant.info", cache_ttl=t, timeout=120)

    resp = go(ttl)
    data = validated_json(resp, "MyVariant.info", refetch=(lambda: go(0)) if ttl else None)
    if not isinstance(data, list):
        raise ValueError(f"MyVariant.info answered a {type(data).__name__}, not a list of hits")
    return resp, data


def _post_complete(ids: Sequence[str], hg: str):
    """One batch; ids the answer left out are asked once more without the cache.

    A partial answer is otherwise stored for 30 days and every re-run would
    replay it, so "re-run to resume" would never fill the gap.
    """
    resp, hits = _post(ids, hg)
    answered = {str(h.get("query")) for h in hits if isinstance(h, dict)}
    missing = [i for i in ids if i not in answered]
    if missing:
        _, more = _post(missing, hg, ttl=0)
        hits = list(hits) + [h for h in more if isinstance(h, dict)]
    return resp, hits


def batch(variants: Sequence[Tuple[str, int, str, str]], assembly: str = "GRCh38", chunk: int = MAX_IDS,
          workers: int = 3, progress: Optional[Callable[[str], None]] = None) -> Outcome:
    """Look up many VCF-style alleles, `chunk` (<= 1000) per POST, a few POSTs at a time.

    The keys are sorted before chunking, so the same input always makes the
    same request bodies: a run that stopped part-way (deadline, network) is
    resumed by running it again, every answered batch coming from the on-disk
    cache (30 days). `result["records"]` maps each key to the parsed hit, or to
    None when MyVariant does not hold it; keys of batches that failed are
    listed in `result["unanswered"]` (never treated as not found).
    """
    hg = ASSEMBLY.get(assembly)
    if hg is None:
        raise ValueError("assembly must be GRCh38 or GRCh37")
    chunk = max(1, min(MAX_IDS, int(chunk)))
    keys = sorted({(str(c), int(p), str(r), str(a)) for c, p, r, a in variants})
    ids: Dict[Tuple[str, int, str, str], Optional[str]] = {k: hgvs_id(*k) for k in keys}
    queryable = [k for k in keys if ids[k]]
    chunks = [queryable[i:i + chunk] for i in range(0, len(queryable), chunk)]
    warnings: List[str] = []
    sources: List[Dict[str, Any]] = []
    records: Dict[Tuple[str, int, str, str], Optional[Dict[str, Any]]] = {}
    unanswered: List[Tuple[str, int, str, str]] = [k for k in keys if not ids[k]]
    done = [0]

    def run(i_part: Tuple[int, List[Tuple[str, int, str, str]]]):
        i, part = i_part
        got = attempt(f"MyVariant.info batch {i + 1}/{len(chunks)} ({len(part)} variants)",
                      lambda: _post_complete([ids[k] for k in part], hg), warnings)
        done[0] += 1
        if progress:
            progress(f"MyVariant.info batch {done[0]}/{len(chunks)} "
                     f"{'answered' + (' (cache)' if got and got[0].cached else '') if got else 'FAILED'}")
        return got

    results: List[Any] = []
    if chunks:
        with ThreadPoolExecutor(max_workers=max(1, min(workers, len(chunks)))) as pool:
            results = list(pool.map(run, enumerate(chunks)))
    answered = cached = 0
    answered_resps: List[Any] = []
    for part, got in zip(chunks, results):
        if got is None:
            unanswered.extend(part)
            continue
        resp, hits = got
        answered += 1
        cached += 1 if resp.cached else 0
        by_query: Dict[str, Dict[str, Any]] = {}
        for h in hits:
            if isinstance(h, dict) and h.get("query") is not None:
                q = str(h["query"])
                # one id can match more than one document; the first that is not `notfound` wins
                if q not in by_query or by_query[q].get("notfound"):
                    by_query[q] = h
        missing = 0
        for k in part:
            h = by_query.get(str(ids[k]))
            if h is None:
                missing += 1
                unanswered.append(k)
            elif h.get("notfound"):
                records[k] = None
            else:
                records[k] = parse(h)
        if missing:
            warnings.append(f"MyVariant.info answered a batch without {missing} of its {len(part)} ids: those were "
                            "not looked up (kept, unfiltered)")
        answered_resps.append(resp)
    if answered_resps:
        # one provenance record for the whole run (an exome is ~20 requests to one endpoint): the
        # request count, the id range and the oldest retrieval time, which is what dates the answers
        oldest = min(answered_resps, key=lambda r: r.retrieved_at)
        rec = record("MyVariant.info", f"{len(queryable)} variants in {answered}/{len(chunks)} POST(s), "
                                       f"{ids[queryable[0]]} … {ids[queryable[-1]]}", oldest,
                     url=f"{BASE}/variant",
                     note=f"{hg}; POST /v1/variant (≤{chunk} ids each); fields gnomAD exome/genome AF, dbNSFP REVEL/"
                          "AlphaMissense/CADD, ClinVar RCV, snpEff impact; variant keys only, no sample data; "
                          f"{cached} of {answered} answered from the local cache")
        sources.append(rec)
    result = {"records": records, "unanswered": unanswered, "batches": len(chunks), "answered_batches": answered,
              "answered_keys": len(records),
              "cached_batches": cached, "queried": len(queryable), "found": sum(1 for v in records.values() if v),
              "not_found": sum(1 for v in records.values() if v is None), "assembly": hg}
    return Outcome(result, sources=sources, warnings=warnings)


def describe(result: Dict[str, Any]) -> str:
    """One line for a run's counts (used in notes and progress)."""
    return (f"{result['queried']} variant(s) in {result['batches']} batch(es): {result['found']} found, "
            f"{result['not_found']} not in MyVariant, {len(result['unanswered'])} not looked up "
            f"({result['answered_batches']}/{result['batches']} batches answered, {result['cached_batches']} from cache)")


