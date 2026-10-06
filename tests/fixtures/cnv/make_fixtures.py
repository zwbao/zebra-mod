"""Capture the offline fixtures for tests/test_cnv.py from the live sources (run by hand, network needed).

    python3 tests/fixtures/cnv/make_fixtures.py

Everything written here is a verbatim upstream body or a field subset of one, with
its URL, retrieval time and sha256 recorded in provenance.json:
- ClinGen dosage bulk files (gene and region curation lists, GRCh38 and GRCh37), gzipped as served.
- Ensembl /lookup/symbol?expand=1 for DMD, UBE3A, SMN1, SMN2: the canonical transcript only,
  with its exons and Translation (the coding start/end the frame arithmetic needs).
- Ensembl /overlap/region?feature=gene for the CNV calls the tests use, in the same 4.5 Mb windows
  zebra.cnv walks, keeping the fields zebra.cnv reads (id, external_name, biotype, start, end, strand).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).parent
CLINGEN = "https://ftp.clinicalgenome.org/"
FILES = ["ClinGen_gene_curation_list_GRCh38.tsv", "ClinGen_gene_curation_list_GRCh37.tsv",
         "ClinGen_region_curation_list_GRCh38.tsv", "ClinGen_region_curation_list_GRCh37.tsv"]
LOOKUPS = ["DMD", "UBE3A", "SMN1", "SMN2"]
WINDOW = 4_500_000  # zebra.cnv.OVERLAP_WINDOW
# (build, chrom, start, end): the CNV calls used in tests/test_cnv.py
CALLS = [
    ("GRCh38", "22", 18648855, 21800471),   # 22q11.2 proximal A-D deletion (DiGeorge)
    ("GRCh38", "1", 849466, 6823542),       # 1p36 terminal deletion
    ("GRCh38", "7", 73330452, 74778226),    # Williams-Beuren 7q11.23
    ("GRCh38", "15", 23123715, 28193120),   # Prader-Willi / Angelman 15q11-q13
    ("GRCh38", "17", 16810000, 20450000),   # Smith-Magenis 17p11.2
    ("GRCh38", "16", 29580020, 30180020),   # 16p11.2 BP4-BP5
    ("GRCh38", "13", 55000000, 56000000),   # gene desert (13q21)
    ("GRCh37", "22", 18648855, 21800471),   # the same 22q11.2 call read as GRCh37
]


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def fetch(url: str, accept: str = "application/json") -> bytes:
    req = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "zebra-mod fixture capture"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                body = resp.read()
            time.sleep(0.3)
            return body
        except urllib.error.HTTPError as err:
            if err.code not in (429, 500, 502, 503, 504) or attempt == 3:
                raise
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(url)


def main() -> None:
    prov = {"captured_by": "tests/fixtures/cnv/make_fixtures.py", "files": {}}
    for name in FILES:
        url = CLINGEN + name
        body = fetch(url, "text/tab-separated-values,text/plain,*/*")
        text = body.decode("utf-8")
        rows = sum(1 for line in text.splitlines() if line and not line.startswith("#"))
        header_date = text.splitlines()[1].lstrip("#").strip()
        (HERE / (name + ".gz")).write_bytes(gzip.compress(body, mtime=0))
        prov["files"][name + ".gz"] = {"url": url, "retrieved_at": now(), "sha256": hashlib.sha256(body).hexdigest(),
                                      "bytes": len(body), "data_rows": rows, "file_header_date": header_date}
    lookups = {}
    for symbol in LOOKUPS:
        for build, host in (("GRCh38", "https://rest.ensembl.org"), ("GRCh37", "https://grch37.rest.ensembl.org")):
            if symbol != "DMD" and build == "GRCh37":
                continue
            url = f"{host}/lookup/symbol/homo_sapiens/{symbol}?expand=1"
            body = fetch(url)
            data = json.loads(body)
            canonical = [t for t in data["Transcript"] if t.get("is_canonical")]
            keep = {k: v for k, v in data.items() if k != "Transcript"}
            keep["Transcript"] = [{k: t.get(k) for k in ("id", "display_name", "biotype", "is_canonical", "length",
                                                          "strand", "seq_region_name", "start", "end", "Exon",
                                                          "Translation")} for t in canonical]
            for t in keep["Transcript"]:
                t["Exon"] = [{k: e.get(k) for k in ("id", "start", "end", "strand", "seq_region_name")} for e in t["Exon"]]
                if t.get("Translation"):
                    t["Translation"] = {k: t["Translation"].get(k) for k in ("id", "start", "end", "length")}
            lookups[f"{build}:{symbol}"] = keep
            prov["files"][f"lookups.json:{build}:{symbol}"] = {
                "url": url, "retrieved_at": now(), "sha256": hashlib.sha256(body).hexdigest(),
                "note": f"canonical transcript only ({len(data['Transcript'])} transcripts in the response)"}
    (HERE / "lookups.json").write_text(json.dumps(lookups, indent=1, sort_keys=True), "utf-8")
    overlap = {}
    for build, chrom, start, end in CALLS:
        host = "https://rest.ensembl.org" if build == "GRCh38" else "https://grch37.rest.ensembl.org"
        w = start
        while w <= end:
            we = min(end, w + WINDOW - 1)
            url = f"{host}/overlap/region/human/{chrom}:{w}-{we}?feature=gene"
            body = fetch(url)
            genes = [{k: g.get(k) for k in ("id", "external_name", "biotype", "start", "end", "strand")}
                     for g in json.loads(body)]
            overlap[f"{build}:{chrom}:{w}-{we}"] = genes
            prov["files"][f"overlap.json.gz:{build}:{chrom}:{w}-{we}"] = {
                "url": url, "retrieved_at": now(), "sha256": hashlib.sha256(body).hexdigest(), "genes": len(genes)}
            w = we + 1
    (HERE / "overlap.json.gz").write_bytes(gzip.compress(json.dumps(overlap, sort_keys=True).encode("utf-8"), mtime=0))
    (HERE / "provenance.json").write_text(json.dumps(prov, indent=1, sort_keys=True), "utf-8")
    print(json.dumps({k: {kk: v.get(kk) for kk in ("data_rows", "genes", "file_header_date") if kk in v}
                      for k, v in prov["files"].items()}, indent=1))


if __name__ == "__main__":
    main()
