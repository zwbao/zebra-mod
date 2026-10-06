#!/usr/bin/env python3
"""Capture the reference sequence the `zebra aso` tests pin their arithmetic to.

The offline tests build most of their sequence by hand, so every window position and
dinucleotide can be checked by eye. Two of them instead use real sequence, so the
cryptic-exon geometry is pinned to the genome and not only to the construction:
CFTR intron 22 around c.3718-2477C>T (the 3849+10kb variant, whose 84 nt cryptic exon
is established in RNA), and the same window on the minus strand of a minus-strand gene.

It also captures the raw GTEx and MaveDB API responses the offline source tests parse,
into tests/fixtures/{gtex,mavedb}/ with a `.provenance.json` beside each one.

Run from the repo root:
    python3 tools/aso/fetch_fixtures.py             # everything
    python3 tools/aso/fetch_fixtures.py sequences   # only tests/fixtures/aso/
    python3 tools/aso/fetch_fixtures.py responses   # only the API responses
Every file is written with the URL, the request body, the service's own retrieval time,
the sha256 and the byte count, so a fixture can always be traced back and re-fetched.
Nothing here is hand-edited: if a fixture does not match its sha256, rebuild it.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from zebra.http import request  # noqa: E402

OUT = ROOT / "tests" / "fixtures" / "aso"
# name -> (assembly, chrom, start, end, why)
REGIONS = {
    "cftr_intron22_3849plus10kb": (
        "GRCh38", "7", 117639860, 117639975,
        "CFTR intron 22 around NM_000492.4:c.3718-2477C>T (7-117639961-C-T). SpliceAI reports a gained acceptor "
        "at 117639876 and a gained donor at 117639959: 84 nt apart, the size of the cryptic exon reported for "
        "this variant. The C>T turns the reference GC at 117639960-117639961 into GT, the canonical donor.",
    ),
    "tp53_exon4_boundaries": (
        "GRCh38", "17", 7675960, 7676310,
        "TP53 exon 4 (7675994-7676272) with flanks, on the minus strand: the acceptor is the HIGHER genomic "
        "coordinate and the donor the lower one, which is what the strand handling has to get right.",
    ),
}


# Raw API responses the offline tests parse, captured as they came. Kept next to the
# sequence fixtures so one script rebuilds everything W7's tests pin.
#   name -> (directory, method, url, body or None)
RESPONSES = {
    "reference_gene_cftr_v26": (
        "gtex", "GET",
        "https://gtexportal.org/api/v2/reference/gene?geneId=CFTR&gencodeVersion=v26"
        "&genomeBuild=GRCh38%2Fhg38&itemsPerPage=50", None),
    "median_expression_cftr_v8": (
        "gtex", "GET",
        "https://gtexportal.org/api/v2/expression/medianGeneExpression"
        "?gencodeId=ENSG00000001626.14&datasetId=gtex_v8&itemsPerPage=250", None),
    "reference_gene_missing": (
        "gtex", "GET",
        "https://gtexportal.org/api/v2/reference/gene?geneId=NOTAGENE123&gencodeVersion=v26"
        "&genomeBuild=GRCh38%2Fhg38&itemsPerPage=50", None),
    "score_sets_search_brca1": (
        "mavedb", "POST", "https://api.mavedb.org/api/v1/score-sets/search",
        '{"targets": ["BRCA1"], "published": true}'),
    "calibration_brca1_sge": (
        "mavedb", "GET",
        "https://api.mavedb.org/api/v1/score-calibrations/score-set/urn:mavedb:00000097-0-2/primary", None),
    "scores_brca1_sge_head": (
        "mavedb", "GET",
        "https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000097-0-2/scores?start=0&limit=200", None),
}


def capture_responses() -> None:
    for name, (area, method, url, body) in RESPONSES.items():
        out_dir = ROOT / "tests" / "fixtures" / area
        out_dir.mkdir(parents=True, exist_ok=True)
        accept = "text/csv" if "/scores" in url else "application/json"
        headers = {"Content-Type": "application/json"} if body is not None else None
        resp = request(url, source=f"fixture {name}", method=method, body=body, headers=headers, accept=accept,
                       cache_ttl=0, timeout=180)
        suffix = "csv" if accept == "text/csv" else "json"
        path = out_dir / f"{name}.{suffix}"
        raw = resp.text.encode("utf-8")
        path.write_bytes(raw)
        provenance = out_dir / f"{name}.provenance.json"
        provenance.write_text(json.dumps({
            "name": name, "url": url, "method": method, "request_body": body, "status": resp.status,
            # the service's own retrieval time, not the time this file was written
            "retrieved_at": resp.retrieved_at,
            "written_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            "bytes": len(raw), "characters": len(resp.text),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "built_by": "tools/aso/fetch_fixtures.py",
        }, indent=1, sort_keys=True) + "\n", "utf-8")
        print(f"{path.relative_to(ROOT)}  {len(raw)} bytes  HTTP {resp.status}")


def capture_sequences() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    host = {"GRCh38": "https://rest.ensembl.org", "GRCh37": "https://grch37.rest.ensembl.org"}
    for name, (assembly, chrom, start, end, why) in REGIONS.items():
        region = f"{chrom}:{start}..{end}:1"
        url = f"{host[assembly]}/sequence/region/human/{region}"
        resp = request(url, source="Ensembl sequence", accept="text/plain", cache_ttl=0)
        seq = resp.text.strip().upper()
        payload = {
            "name": name,
            "why": why,
            "assembly": assembly,
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": 1,
            "sequence": seq,
            "length": len(seq),
            "expected_length": end - start + 1,
            "sha256": hashlib.sha256(seq.encode("ascii")).hexdigest(),
            "source": {"db": "Ensembl sequence", "url": url, "retrieved_at": datetime.now(timezone.utc)
                       .replace(microsecond=0).isoformat(), "status": resp.status},
            "built_by": "tools/aso/fetch_fixtures.py",
        }
        path = OUT / f"{name}.json"
        path.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", "utf-8")
        print(f"{path.relative_to(ROOT)}  {len(seq)} bp  sha256 {payload['sha256'][:16]}...")


def main(argv: list) -> int:
    what = (argv[1] if len(argv) > 1 else "all").lower()
    if what not in ("all", "sequences", "responses"):
        print("usage: fetch_fixtures.py [all|sequences|responses]")
        return 2
    if what in ("all", "sequences"):
        capture_sequences()
    if what in ("all", "responses"):
        capture_responses()
    return 0


if __name__ == "__main__":
    os.environ.setdefault("ZEBRA_NO_CACHE", "1")
    raise SystemExit(main(sys.argv))
