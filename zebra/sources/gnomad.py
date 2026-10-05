"""gnomAD GraphQL: variant frequencies (with site coverage) and gene constraint.

API: POST https://gnomad.broadinstitute.org/api (no key). gnomAD asks for
light use (about 10 requests a minute; no bulk scraping), so each variant is
one request (frequencies, coverage at the site and liftover together) and
answers are cached for 30 days.

Datasets: GRCh38 -> gnomad_r4 (v4.1, exomes + genomes + joint);
GRCh37 -> gnomad_r2_1 (v2.1.1).
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from zebra.core import Outcome
from zebra.http import SourceError, post_json, source_record

API = "https://gnomad.broadinstitute.org/api"
BROWSER = "https://gnomad.broadinstitute.org"
DATASETS = {"GRCh38": "gnomad_r4", "GRCh37": "gnomad_r2_1"}
CACHE_TTL = 30 * 86400
# grpmax / popmax is taken over non-bottlenecked continental groups only, as
# gnomAD does: v4 excludes ami, asj, fin, mid and "remaining"; v2 excludes asj,
# fin and oth.
GRPMAX_GROUPS = ("afr", "amr", "eas", "nfe", "sas")
GROUP_NAMES = {
    "afr": "African/African American", "amr": "Admixed American", "asj": "Ashkenazi Jewish",
    "eas": "East Asian", "fin": "Finnish", "nfe": "Non-Finnish European", "sas": "South Asian",
    "mid": "Middle Eastern", "ami": "Amish", "remaining": "Remaining", "oth": "Other",
}
COVERED_MEAN = 20.0  # a site counts as covered when exome or genome mean depth >= 20

_SEQ_FIELDS = """ac an af homozygote_count hemizygote_count filters flags
      faf95 { popmax popmax_population }
      populations { id ac an homozygote_count hemizygote_count }"""

VARIANT_QUERY = """
query ZebraVariant($id: String!, $ds: DatasetId!, $chrom: String!, $pos: Int!, $rg: ReferenceGenomeId!) {
  variant(variantId: $id, dataset: $ds) {
    variant_id reference_genome rsids caid flags
    exome { %(exome)s }
    genome { %(genome)s }
    %(joint)s
    coverage { exome { mean median over_20 } genome { mean median over_20 } }
    in_silico_predictors { id value flags }
  }
  region(chrom: $chrom, start: $pos, stop: $pos, reference_genome: $rg) {
    coverage(dataset: $ds) { exome { pos mean median over_20 } genome { pos mean median over_20 } }
  }
  %(liftover)s
}
"""

CONSTRAINT_QUERY = """
query ZebraConstraint(%(var)s: String!, $rg: ReferenceGenomeId!) {
  gene(%(arg)s: %(var)s, reference_genome: $rg) {
    gene_id symbol hgnc_id chrom start stop
    gnomad_constraint { pLI oe_lof oe_lof_lower oe_lof_upper mis_z oe_mis oe_mis_lower oe_mis_upper syn_z oe_syn
                        exp_lof obs_lof exp_mis obs_mis flags }
  }
}
"""


def dataset_for(assembly: str) -> str:
    if assembly not in DATASETS:
        raise ValueError("assembly must be GRCh38 or GRCh37")
    return DATASETS[assembly]


def variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    chrom = str(chrom).upper().replace("CHR", "")
    return f"{'MT' if chrom == 'M' else chrom}-{int(pos)}-{ref.upper()}-{alt.upper()}"


def build_variant_query(dataset: str) -> str:
    v4 = dataset == "gnomad_r4"
    seq = _SEQ_FIELDS
    joint = ("joint { ac an homozygote_count hemizygote_count filters faf95 { popmax popmax_population } "
             "populations { id ac an homozygote_count hemizygote_count } }") if v4 else ""
    if v4:
        lift = "liftover(liftover_variant_id: $id, reference_genome: GRCh38) { source { variant_id reference_genome } }"
    else:
        lift = "liftover(source_variant_id: $id, reference_genome: GRCh37) { liftover { variant_id reference_genome } }"
    return VARIANT_QUERY % {"exome": seq, "genome": seq, "joint": joint, "liftover": lift}


EXPECTED_ERRORS = ("Variant not found", "Gene not found")


def _graphql(query: str, variables: Dict[str, Any], label: str) -> Tuple[Any, Dict[str, Any]]:
    """POST a query; HTTP 200 can still carry `errors`. An unexpected error read from the cache is
    refetched once, so a transient server error is not served from the cache for 30 days."""
    payload = {"query": query, "variables": variables}
    resp = post_json(API, payload, source="gnomAD", cache_ttl=CACHE_TTL, timeout=60)
    body = resp.json()
    errs = [e.get("message", "") for e in body.get("errors") or []]
    if resp.cached and any(e not in EXPECTED_ERRORS for e in errs):
        resp = post_json(API, payload, source="gnomAD", cache_ttl=0, timeout=60)
        body = resp.json()
    return resp, body


def variant(chrom: str, pos: int, ref: str, alt: str, assembly: str = "GRCh38") -> Outcome:
    """Frequencies of one VCF-style (left-normalised) variant, with coverage at the site."""
    ds = dataset_for(assembly)
    vid = variant_id(chrom, pos, ref, alt)
    variables = {"id": vid, "ds": ds, "chrom": vid.split("-")[0], "pos": int(pos), "rg": assembly}
    resp, body = _graphql(build_variant_query(ds), variables, vid)
    errors = [e.get("message", "") for e in body.get("errors") or []]
    other = [e for e in errors if e != "Variant not found"]
    data = body.get("data") or {}
    if other and not data.get("variant") and not data.get("region"):
        raise SourceError("gnomAD", resp.url, resp.status, "; ".join(other))
    result = parse_variant(data, ds, vid)
    warnings: List[str] = []
    if other:
        warnings.append(f"gnomAD reported: {'; '.join(other)}")
    return Outcome(result, sources=[source_record("gnomAD", f"{vid} ({ds})", resp, url=result["url"])],
                   warnings=warnings)


def _seq(d: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not d:
        return None
    faf = d.get("faf95") or {}
    an = d.get("an")
    ac = d.get("ac")
    af = d.get("af")
    if af is None and ac is not None and an:
        af = ac / an
    return {"ac": ac, "an": an, "af": af, "hom": d.get("homozygote_count"), "hemi": d.get("hemizygote_count"),
            "filters": d.get("filters") or [], "faf95": faf.get("popmax"), "faf95_group": faf.get("popmax_population")}


def _top_groups(pops: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """Genetic-ancestry groups only: drop sex splits, sub-groups (nfe_swe), HGDP/1KG labels."""
    out: Dict[str, Dict[str, Any]] = {}
    for p in pops or []:
        gid = p.get("id") or ""
        if not gid or "_" in gid or ":" in gid or gid in ("XX", "XY"):
            continue
        if gid in out:  # joint lists some rows twice
            continue
        out[gid] = {"ac": p.get("ac") or 0, "an": p.get("an") or 0, "hom": p.get("homozygote_count") or 0,
                    "hemi": p.get("hemizygote_count")}
    return out


def _sum_groups(a: Dict[str, Dict[str, Any]], b: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for gid in set(a) | set(b):
        x, y = a.get(gid, {}), b.get(gid, {})
        out[gid] = {"ac": (x.get("ac") or 0) + (y.get("ac") or 0), "an": (x.get("an") or 0) + (y.get("an") or 0),
                    "hom": (x.get("hom") or 0) + (y.get("hom") or 0),
                    "hemi": (x.get("hemi") or 0) + (y.get("hemi") or 0) if (x.get("hemi") is not None or y.get("hemi") is not None) else None}
    return out


def grpmax(groups: Dict[str, Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    best: Optional[Dict[str, Any]] = None
    for gid in GRPMAX_GROUPS:
        g = groups.get(gid)
        if not g or not g.get("an"):
            continue
        af = g["ac"] / g["an"]
        if best is None or af > best["af"]:
            best = {"group": gid, "name": GROUP_NAMES.get(gid), "af": af, "ac": g["ac"], "an": g["an"]}
    return best


def _site_coverage(region: Optional[Dict[str, Any]], pos: int) -> Dict[str, Any]:
    cov = (region or {}).get("coverage") or {}
    out: Dict[str, Any] = {}
    for kind in ("exome", "genome"):
        bins = cov.get(kind) or []
        hit = None
        for b in bins:
            if b.get("pos") == pos:
                hit = b
                break
        if hit is None and bins:  # coverage is binned; take the nearest bin
            hit = min(bins, key=lambda b: abs((b.get("pos") or 0) - pos))
        if hit is not None:
            out[kind] = {"mean": _r(hit.get("mean")), "median": hit.get("median"), "over_20": _r(hit.get("over_20")),
                         "pos": hit.get("pos")}
    return out


def _r(x: Any, n: int = 3) -> Any:
    return round(x, n) if isinstance(x, float) else x


def covered(coverage: Dict[str, Any]) -> Optional[bool]:
    means = [(coverage.get(k) or {}).get("mean") for k in ("exome", "genome")]
    means = [m for m in means if isinstance(m, (int, float))]
    if not means:
        return None
    return any(m >= COVERED_MEAN for m in means)


def parse_variant(data: Dict[str, Any], dataset: str, vid: str) -> Dict[str, Any]:
    v = data.get("variant")
    pos = int(vid.split("-")[1])
    url = f"{BROWSER}/variant/{vid}?dataset={dataset}"
    lift = data.get("liftover") or []
    other = None
    for row in lift:
        node = row.get("source") or row.get("liftover") or {}
        if node.get("variant_id"):
            other = {"assembly": node.get("reference_genome"), "variant_id": node.get("variant_id"),
                     "source": "gnomAD liftover"}
            break
    if not v:
        cov = _site_coverage(data.get("region"), pos)
        is_cov = covered(cov)
        return {
            "dataset": dataset, "variant_id": vid, "url": url, "found": False, "absent": True,
            "coverage": cov, "covered": is_cov,
            "coverage_rule": f"covered = exome or genome mean depth >= {COVERED_MEAN:g} at the site",
            "total": {"ac": 0 if is_cov else None, "an": None, "af": 0.0 if is_cov else None},
            "grpmax": None, "faf95": None, "liftover": other,
            "note": ("absent from gnomAD at a covered site" if is_cov else
                     "absent from gnomAD, but coverage at the site is low or unknown: absence is weak evidence"),
        }
    exome, genome = _seq(v.get("exome")), _seq(v.get("genome"))
    joint_raw = v.get("joint")
    joint = _seq(joint_raw) if joint_raw else None
    ex_groups = _top_groups((v.get("exome") or {}).get("populations") or [])
    ge_groups = _top_groups((v.get("genome") or {}).get("populations") or [])
    if joint_raw:
        groups = _top_groups(joint_raw.get("populations") or [])
        basis = "joint exome+genome (gnomAD v4)"
    else:
        groups = _sum_groups(ex_groups, ge_groups)
        basis = "exome+genome summed per group"
    if joint:
        total = {"ac": joint["ac"], "an": joint["an"], "af": joint["af"], "hom": joint["hom"], "hemi": joint["hemi"]}
    else:
        ac = sum((x or {}).get("ac") or 0 for x in (exome, genome))
        an = sum((x or {}).get("an") or 0 for x in (exome, genome))
        total = {"ac": ac, "an": an, "af": (ac / an) if an else None,
                 "hom": sum((x or {}).get("hom") or 0 for x in (exome, genome)),
                 "hemi": sum((x or {}).get("hemi") or 0 for x in (exome, genome))}
    gmax = grpmax(groups)
    if gmax:
        gmax["basis"] = f"{basis}; groups {'/'.join(GRPMAX_GROUPS)} (bottlenecked groups excluded)"
    if joint and joint.get("faf95") is not None:
        faf = {"value": joint["faf95"], "group": joint["faf95_group"], "basis": "joint faf95 (grpmax filtering AF)"}
    else:
        cands = [(x["faf95"], x["faf95_group"], k) for k, x in (("exome", exome), ("genome", genome))
                 if x and x.get("faf95") is not None]
        faf = None
        if cands:
            val, grp, kind = max(cands)
            faf = {"value": val, "group": grp, "basis": f"higher of exome/genome faf95 popmax ({kind})"}
    populations = []
    for gid, g in sorted(groups.items(), key=lambda kv: -(kv[1]["ac"] / kv[1]["an"] if kv[1]["an"] else 0)):
        populations.append({"id": gid, "name": GROUP_NAMES.get(gid), "ac": g["ac"], "an": g["an"],
                            "af": (g["ac"] / g["an"]) if g["an"] else None, "hom": g["hom"]})
    cov = {}
    for kind in ("exome", "genome"):
        c = (v.get("coverage") or {}).get(kind)
        if c:
            cov[kind] = {"mean": _r(c.get("mean")), "median": c.get("median"), "over_20": _r(c.get("over_20"))}
    if not cov:
        cov = _site_coverage(data.get("region"), pos)
    predictors = {p["id"]: p.get("value") for p in v.get("in_silico_predictors") or [] if p.get("id")}
    filters = sorted(set((exome or {}).get("filters", []) + (genome or {}).get("filters", []) + ((joint or {}).get("filters") or [])))
    return {
        "dataset": dataset, "variant_id": v.get("variant_id") or vid, "url": url, "found": True, "absent": False,
        "rsids": v.get("rsids") or [], "caid": v.get("caid"), "flags": v.get("flags") or [], "filters": filters,
        "exome": exome, "genome": genome, "joint": joint, "total": total,
        "grpmax": gmax, "faf95": faf, "populations": populations[:12],
        "coverage": cov, "covered": covered(cov), "in_silico": predictors or None, "liftover": other,
    }


def gene_constraint(symbol: str, assembly: str = "GRCh38", gene_id: Optional[str] = None) -> Outcome:
    """pLI, LOEUF (oe_lof_upper), missense Z ... GRCh38 = gnomAD v4 constraint, GRCh37 = v2.1.1.

    `gene_id` (Ensembl ENSG) is used when given: gnomAD's symbols follow GENCODE and can lag HGNC.
    """
    dataset_for(assembly)
    if gene_id:
        query = CONSTRAINT_QUERY % {"var": "$gid", "arg": "gene_id"}
        variables = {"gid": gene_id.split(".")[0], "rg": assembly}
    else:
        query = CONSTRAINT_QUERY % {"var": "$symbol", "arg": "gene_symbol"}
        variables = {"symbol": symbol, "rg": assembly}
    resp, body = _graphql(query, variables, symbol)
    errors = [e.get("message", "") for e in body.get("errors") or []]
    gene = (body.get("data") or {}).get("gene")
    if not gene:
        raise SourceError("gnomAD", resp.url, resp.status, "; ".join(errors) or f"no gene {symbol}")
    result = parse_constraint(gene, assembly)
    warnings = []
    if result.get("flags"):
        warnings.append(f"gnomAD constraint flags for {symbol}: {', '.join(result['flags'])} — read the affected metrics "
                        "with caution; a null value (e.g. no_exp_lof) is missing data, not tolerance")
    return Outcome(result, sources=[source_record("gnomAD constraint", f"{symbol} ({'v4' if assembly == 'GRCh38' else 'v2.1.1'})",
                                                  resp, url=result["url"])], warnings=warnings)


def parse_constraint(gene: Dict[str, Any], assembly: str) -> Dict[str, Any]:
    c = gene.get("gnomad_constraint") or {}
    ds = dataset_for(assembly)
    return {
        "gene_id": gene.get("gene_id"), "symbol": gene.get("symbol"), "hgnc_id": gene.get("hgnc_id"),
        "version": "gnomAD v4.1" if assembly == "GRCh38" else "gnomAD v2.1.1",
        "pLI": _r(c.get("pLI"), 4), "loeuf": _r(c.get("oe_lof_upper")), "oe_lof": _r(c.get("oe_lof")),
        "oe_lof_ci": [_r(c.get("oe_lof_lower")), _r(c.get("oe_lof_upper"))],
        "obs_lof": c.get("obs_lof"), "exp_lof": _r(c.get("exp_lof"), 1),
        "mis_z": _r(c.get("mis_z"), 2), "oe_mis": _r(c.get("oe_mis")), "oe_mis_ci": [_r(c.get("oe_mis_lower")), _r(c.get("oe_mis_upper"))],
        "syn_z": _r(c.get("syn_z"), 2), "oe_syn": _r(c.get("oe_syn")),
        "flags": c.get("flags") or [],
        "url": f"{BROWSER}/gene/{gene.get('gene_id')}?dataset={ds}",
        "reading": ("lower LOEUF / higher pLI = fewer loss-of-function variants observed than expected; "
                    "read with the gene's established disease mechanism, not alone"),
    }
