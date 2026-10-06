"""GTEx Portal API v2: median gene expression per tissue, with the tissue's ontology id.

API: https://gtexportal.org/api/v2 (no key). Two calls per gene:

  /reference/gene?geneId=<symbol>&gencodeVersion=<v>&genomeBuild=GRCh38/hg38
      resolves a symbol to the GENCODE id of the release the dataset was built on.
  /expression/medianGeneExpression?gencodeId=<id>&datasetId=<dataset>
      median TPM per tissue, each row carrying GTEx's own `ontologyId`.

The GENCODE id must match the dataset's GENCODE release or the expression call
answers an empty list with HTTP 200 (verified live 2026-10-06: CFTR is
ENSG00000001626.14 in GENCODE v26 / `gtex_v8` and ENSG00000001626.17 in v39 /
`gtex_v10`; asking for `.16`, or omitting `datasetId`, returns no rows at all).
That is why this module resolves the id per dataset instead of taking one from
Ensembl, and why `datasetId` is always sent.

Two of the 54 GTEx v8 sites are cell lines and carry EFO rather than UBERON ids
(cultured fibroblasts EFO:0002009, EBV-transformed lymphocytes EFO:0000572).
The contract field is still called `uberon` and carries the id GTEx gave, so it
can be handed to an ontology-keyed model as it stands; `ontology_is_uberon`
says which prefix it actually is.

What this is for: deciding which tissue an RNA test could be run on, and which
ontology term a tissue-aware model should be asked about. A median TPM is a
bulk-tissue summary from adult post-mortem donors -- not a detection limit, not
cell-type resolved, and systematically low for a transcript that is degraded by
nonsense-mediated decay.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

from zebra.core import Outcome, UsageError
from zebra.http import SourceError, get_json
from zebra.sources import attempt, record, validated_json

BASE = "https://gtexportal.org/api/v2"
# datasetId -> the GENCODE release that dataset's gene ids come from (verified live
# against /reference/gene for both releases on 2026-10-06).
DATASETS = {"gtex_v8": "v26", "gtex_v10": "v39"}
DEFAULT_DATASET = "gtex_v8"
GENOME_BUILD = "GRCh38/hg38"

# The tissues that decide which RNA test is practical, in the order a clinical lab
# would consider them: blood needs only a tube, fibroblasts need a skin biopsy,
# muscle needs a muscle biopsy.
RNA_TEST_TISSUES = (
    ("Whole_Blood", "blood: a PAXgene/EDTA tube, no biopsy"),
    ("Cells_EBV-transformed_lymphocytes", "lymphoblastoid line from blood (GTEx measures a cell line, not fresh blood)"),
    ("Cells_Cultured_fibroblasts", "cultured fibroblasts: needs a skin biopsy, then weeks of culture"),
    ("Skin_Sun_Exposed_Lower_leg", "skin itself, the tissue a fibroblast culture is started from"),
    ("Skin_Not_Sun_Exposed_Suprapubic", "skin itself, the tissue a fibroblast culture is started from"),
    ("Muscle_Skeletal", "skeletal muscle: needs a muscle biopsy"),
)
# A stated heuristic for reading a median TPM, not a validated assay limit: where the
# detection floor of a given RT-PCR or RNA-seq run sits depends on the run.
EXPRESSION_BANDS = (
    (10.0, "expressed", "comfortably detectable in bulk RNA from this tissue"),
    (1.0, "expressed_low", "detectable by RT-PCR in most hands; RNA-seq coverage may be thin"),
    (0.1, "trace", "at or below what bulk RNA-seq resolves; RT-PCR may still work, but expect to optimise"),
    (0.0, "not_detected", "GTEx sees essentially nothing here; do not plan an RNA test on this tissue"),
)
CAVEATS = (
    "Median TPM over GTEx donors for a bulk tissue: it is not cell-type resolved, so a gene expressed only in a "
    "minority cell type of that tissue reads low.",
    "GTEx donors are adult post-mortem. There is no fetal tissue and no disease tissue here.",
    "Cultured fibroblasts and EBV-transformed lymphocytes are cell lines GTEx profiled, not fresh patient cells; "
    "expression in culture can differ from the tissue they came from.",
    "An allele whose transcript is degraded by nonsense-mediated decay reads LOWER than the gene's expression "
    "suggests. An RNA test for an NMD-targeted event needs NMD inhibition (e.g. cycloheximide or puromycin on "
    "cultured cells) or it can miss the aberrant transcript entirely.",
    "The expressed/trace/not_detected bands here are a stated reading heuristic, not a detection limit: the floor "
    "of a particular RT-PCR or RNA-seq run is a property of that run.",
)


def _band(tpm: Optional[float]) -> Dict[str, Optional[str]]:
    if tpm is None:
        return {"band": None, "meaning": None}
    for floor, name, meaning in EXPRESSION_BANDS:
        if tpm >= floor:
            return {"band": name, "meaning": meaning}
    return {"band": "not_detected", "meaning": EXPRESSION_BANDS[-1][2]}


def resolve_gene(symbol: str, dataset: str = DEFAULT_DATASET) -> Outcome:
    """Symbol -> the GENCODE id of the release `dataset` was built on."""
    symbol = (symbol or "").strip()
    if not symbol:
        raise ValueError("give a gene symbol")
    if dataset not in DATASETS:
        raise ValueError(f"dataset must be one of {', '.join(sorted(DATASETS))}")
    version = DATASETS[dataset]
    params = {"geneId": symbol, "gencodeVersion": version, "genomeBuild": GENOME_BUILD, "itemsPerPage": 50}
    resp = get_json(f"{BASE}/reference/gene", source="GTEx reference/gene", params=params,
                    cache_ttl=30 * 86400, timeout=60)
    data = validated_json(resp, "GTEx reference/gene", require="data",
                          refetch=lambda: get_json(f"{BASE}/reference/gene", source="GTEx reference/gene",
                                                   params=params, cache_ttl=30 * 86400, timeout=60, refresh=True))
    rows = [r for r in (data.get("data") or []) if isinstance(r, dict)]
    src = record("GTEx reference/gene", f"{symbol} {version}", resp,
                 note=f"GENCODE {version} ({dataset}), {GENOME_BUILD}")
    if not rows:
        return Outcome({"gene": symbol, "gencode_id": None, "dataset": dataset, "gencode_version": version,
                        "found": False},
                       sources=[src],
                       warnings=[f"GTEx ({dataset}, GENCODE {version}) has no gene record for {symbol!r}: it may be "
                                 f"an alias, a withdrawn symbol, or not in that GENCODE release. This is "
                                 f"'not found in GTEx', not 'not expressed'."])
    exact = [r for r in rows if str(r.get("geneSymbol") or "").upper() == symbol.upper()] or rows
    chosen = exact[0]
    warnings: List[str] = []
    if len(rows) == 1 and not [r for r in rows if str(r.get("geneSymbol") or "").upper() == symbol.upper()]:
        warnings.append("GTEx matched %r to %s / %s rather than to that symbol exactly; check it is the gene you "
                        "meant" % (symbol, chosen.get("geneSymbol"), chosen.get("gencodeId")))
    if len(exact) > 1:
        warnings.append(f"GTEx returned {len(exact)} genes for {symbol!r} "
                        f"({', '.join(str(r.get('gencodeId')) for r in exact[:4])}); using {chosen.get('gencodeId')}")
    elif len(rows) > 1:
        warnings.append(f"GTEx matched {len(rows)} genes for {symbol!r} and none by exact symbol; "
                        f"using {chosen.get('geneSymbol')} / {chosen.get('gencodeId')}")
    return Outcome({"gene": chosen.get("geneSymbol") or symbol, "gencode_id": chosen.get("gencodeId"),
                    "dataset": dataset, "gencode_version": version, "found": True,
                    "gene_type": chosen.get("geneType"), "chromosome": chosen.get("chromosome"),
                    "entrez_gene_id": chosen.get("entrezGeneId"), "strand": chosen.get("strand"),
                    "description": chosen.get("description")},
                   sources=[src], warnings=warnings)


def median_expression(gencode_id: str, dataset: str = DEFAULT_DATASET) -> Outcome:
    """Median TPM per tissue for one GENCODE id, as GTEx returns it (all tissues)."""
    if dataset not in DATASETS:
        raise ValueError(f"dataset must be one of {', '.join(sorted(DATASETS))}")
    params = {"gencodeId": gencode_id, "datasetId": dataset, "itemsPerPage": 250}
    resp = get_json(f"{BASE}/expression/medianGeneExpression", source="GTEx medianGeneExpression",
                    params=params, cache_ttl=30 * 86400, timeout=90)
    data = validated_json(resp, "GTEx medianGeneExpression", require="data",
                          refetch=lambda: get_json(f"{BASE}/expression/medianGeneExpression",
                                                   source="GTEx medianGeneExpression", params=params,
                                                   cache_ttl=30 * 86400, timeout=90, refresh=True))
    payload = data.get("data")
    if payload is None:
        payload = []
    if not isinstance(payload, list):
        raise SourceError("GTEx medianGeneExpression", BASE, resp.status,
                          "`data` is a %s, not a list of tissue rows" % type(payload).__name__)
    rows = []
    mismatched = set()
    for r in payload:
        if not isinstance(r, dict):
            continue
        got = r.get("gencodeId")
        if got and str(got) != str(gencode_id):
            mismatched.add(str(got))
            continue
        try:
            tpm = float(r.get("median"))
        except (TypeError, ValueError):
            tpm = None
        if tpm is not None and not math.isfinite(tpm):
            tpm = None  # a NaN would sort to the top and read as a band it is not in
        ontology = r.get("ontologyId")
        rows.append({
            "tissue": r.get("tissueSiteDetailId"),
            "uberon": ontology,  # the contract's field name; see the module docstring
            "ontology_id": ontology,
            "ontology_is_uberon": bool(ontology and str(ontology).startswith("UBERON:")),
            "median_tpm": tpm,
            "unit": r.get("unit") or "TPM",
            **_band(tpm),
        })
    src = record("GTEx medianGeneExpression", f"{gencode_id} {dataset}", resp,
                 note=f"{len(rows)} tissues, median TPM")
    paging = data.get("paging_info") or {}
    warnings: List[str] = []
    if mismatched:
        warnings.append("GTEx returned rows for %s as well as %s; only the asked-for gene's rows were kept"
                        % (", ".join(sorted(mismatched)[:3]), gencode_id))
    missing = sum(1 for r in rows if r["median_tpm"] is None)
    if missing:
        warnings.append("GTEx gave no usable median for %d of the %d tissues returned for %s; those are reported "
                        "as 'not measured', never as zero" % (missing, len(rows), gencode_id))
    if not rows:
        warnings.append(f"GTEx {dataset} returned no expression rows for {gencode_id}: the GENCODE id does not "
                        f"belong to that dataset's release (GENCODE {DATASETS[dataset]}). This is a lookup "
                        f"mismatch, not an absence of expression.")
    elif (paging.get("numberOfPages") or 1) > 1:
        warnings.append(f"GTEx paginated the tissue list for {gencode_id} ({paging.get('totalNumberOfItems')} "
                        f"rows); only the first {len(rows)} were read")
    return Outcome(rows, sources=[src], warnings=warnings)


def top_tissues(symbol: str, n: int = 10, dataset: str = DEFAULT_DATASET) -> Outcome:
    """Median expression per tissue for `symbol`: the top `n`, plus the RNA-test tissues.

    `result["tissues"]` is the top `n` by median TPM; `result["rna_test_tissues"]`
    always carries blood, lymphoblastoid, fibroblasts, skin and skeletal muscle,
    whether or not they are in the top `n`, because those are the ones that decide
    whether an RNA test on this gene is practical.
    """
    try:
        n = int(n)
    except (TypeError, ValueError):
        raise UsageError("--top must be a whole number") from None
    if not 1 <= n <= 60:
        raise UsageError("--top must be between 1 and 60 (GTEx v8 has 54 tissues)")
    gene_out = resolve_gene(symbol, dataset=dataset)
    gene = gene_out.result
    sources = list(gene_out.sources)
    warnings = list(gene_out.warnings)
    base: Dict[str, Any] = {
        "gene": gene.get("gene") or symbol,
        "gencode_id": gene.get("gencode_id"),
        "dataset": dataset,
        "gencode_version": gene.get("gencode_version"),
        "unit": "median TPM across GTEx donors",
        "tissues": [],
        "rna_test_tissues": [],
        "tissues_total": 0,
        "caveats": list(CAVEATS),
    }
    if not gene.get("found"):
        base["note"] = f"GTEx has no gene record for {symbol!r} in {dataset}; nothing was measured, nothing is implied"
        return Outcome(base, sources=sources, warnings=warnings,
                       query={"gene": symbol, "n": n, "dataset": dataset})
    expr_out = attempt("GTEx medianGeneExpression",
                       lambda: median_expression(str(gene["gencode_id"]), dataset=dataset), warnings)
    if expr_out is None:
        base["note"] = ("GTEx returned no expression table for %s in %s; see the warnings. Nothing is implied "
                        "about where this gene is or is not expressed." % (base["gene"], dataset))
        base["rna_test_hint"] = ("the expression table could not be retrieved, so which RNA test is feasible "
                                 "cannot be said from GTEx here")
        base["expression_retrieved"] = False
        return Outcome(base, sources=sources, warnings=warnings,
                       query={"gene": symbol, "n": n, "dataset": dataset})
    sources.extend(expr_out.sources)
    warnings.extend(expr_out.warnings)
    rows = expr_out.result or []
    base["expression_retrieved"] = bool(rows)
    by_name = {r["tissue"]: r for r in rows}
    ranked = sorted(rows, key=lambda r: (r["median_tpm"] is None, -(r["median_tpm"] or 0.0)))
    picks = []
    for name, why in RNA_TEST_TISSUES:
        row = by_name.get(name)
        if row is None:
            picks.append({"tissue": name, "uberon": None, "ontology_id": None, "median_tpm": None,
                          "band": None, "meaning": None, "why_it_matters": why,
                          "note": (f"{name} is not in {dataset}; not measured here" if rows
                                   else "no expression table was returned, so this tissue was not read at all")})
        else:
            picks.append(dict(row, why_it_matters=why))
    base.update(tissues=ranked[:n], tissues_total=len(rows), rna_test_tissues=picks)
    if not rows:
        # No table came back: the GENCODE id did not belong to this dataset's release, or
        # GTEx returned nothing. Saying "no accessible tissue reaches 1 TPM" here would be
        # a negative finding drawn from a lookup failure.
        base["rna_test_hint"] = ("no expression table was returned for %s in %s, so which RNA test is feasible "
                                 "cannot be said from GTEx here. This is a lookup result, not a finding about "
                                 "expression; see the warnings." % (base["gene"], dataset))
        return Outcome(base, sources=sources, warnings=warnings,
                       query={"gene": symbol, "n": n, "dataset": dataset})
    unmeasured = [r["tissue"] for r in picks if r.get("median_tpm") is None]
    best = next((r for r in picks if r.get("median_tpm") is not None and r["median_tpm"] >= 1.0), None)
    if best is not None:
        base["rna_test_hint"] = (f"{best['tissue']} at {best['median_tpm']:.3g} TPM is the most accessible tissue "
                                 f"here with the gene above 1 TPM ({best['why_it_matters']}). Confirm against the "
                                 f"caveats below before planning a test.")
    else:
        base["rna_test_hint"] = ("none of the accessible tissues (blood, lymphoblastoid, fibroblasts, skin, muscle) "
                                 "reaches 1 TPM for this gene in GTEx. An RNA test would need the disease tissue, "
                                 "an induced or differentiated cell model, or a minigene assay instead."
                                 + (" Note that %s had no usable median, so %s not covered by that statement."
                                    % (", ".join(unmeasured), "they are" if len(unmeasured) > 1 else "it is")
                                    if unmeasured else ""))
    if ranked:
        base["highest_tissue"] = {k: ranked[0][k] for k in ("tissue", "uberon", "median_tpm", "band")}
        base["ontology_hint"] = (f"for a tissue-keyed model, the highest-expressing GTEx site is "
                                 f"{ranked[0]['tissue']} ({ranked[0]['uberon']}); pick the tissue the DISEASE acts "
                                 f"in, which is not always the highest-expressing one")
    return Outcome(base, sources=sources, warnings=warnings,
                   query={"gene": symbol, "n": n, "dataset": dataset})
