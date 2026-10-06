"""Build a realistic exome-scale synthetic trio VCF (~25-30k records) from real public variant sites.

    python tests/fixtures/vcf/make_exome_scale.py OUT.vcf [--seed 7]

Sites: real GRCh38 variant sites (position, REF, ALT, gnomAD v2.1.1 exome allele
frequency) pulled from MyVariant.info's query API in four frequency bands. Only
site-level public data is used; no individual's genotypes are fetched.

Genotypes: invented. Mother and father are drawn under Hardy-Weinberg from each
site's gnomAD AF (common and low-frequency bands); each rare-band site is given
to exactly one parent as a heterozygote (an individual carries a few hundred
rare variants, which HWE at AF 1e-4 would almost never produce); "novel"
variants put a different ALT base at a real site (REF stays correct) and are
private to one parent. The child (P, female) inherits one allele per parent
per site; then planted findings: SCN1A c.2134C>T de novo, a mosaic de novo,
CFTR F508del (maternal) + c.350G>A p.Arg117His (paternal) in trans, an ALDH7A1
homozygote with two carrier parents. Depth/GQ/AD are invented.

The output header records the queries, retrieval time, sha256 of the fetched
site list and the counts. The file is written outside the repository (it is
~3 MB); tests do not read it — it is the measurement input for `zebra vcf triage
--prefilter myvariant` and `zebra qc` at exome scale.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))

from zebra.http import get_json, now_iso  # noqa: E402

BASE = "https://myvariant.info/v1/query"
BANDS = (  # (name, lucene range on gnomad_exome.af.af, sites to take)
    ("common", "[0.05 TO 1]", 36000),
    ("uncommon", "[0.01 TO 0.05}", 9000),
    ("low", "[0.001 TO 0.01}", 1500),
    ("rare", "[0.00001 TO 0.001}", 1500),
)
CHROMS = [str(i) for i in range(1, 23)] + ["X"]
PLANTED = [  # chrom, pos, ref, alt, genotypes (M, F, P), what
    ("2", 166042334, "G", "A", ("0/0", "0/0", "0/1"), "SCN1A c.2134C>T p.Arg712* de novo"),
    ("12", 52000000 + 123, "C", "T", ("0/0", "0/0", "0/1:mosaic"), "mosaic de novo (VAF ~0.15), invented site"),
    ("7", 117559590, "ATCT", "A", ("0/1", "0/0", "0/1"), "CFTR F508del, maternal"),
    ("7", 117530975, "G", "A", ("0/0", "0/1", "0/1"), "CFTR c.350G>A p.Arg117His, paternal"),
    ("5", 126554292, "C", "G", ("0/1", "0/1", "1/1"), "ALDH7A1 homozygote, carrier parents"),
]


def fetch_band(rng: str, want: int, log) -> list:
    sites = []
    params = {"q": f"gnomad_exome.af.af:{rng}", "assembly": "hg38", "fields": "vcf,chrom,gnomad_exome.af.af",
              "fetch_all": "true"}
    resp = get_json(BASE, source="MyVariant.info", params=params, cache_ttl=0, timeout=120)
    data = resp.json()
    pages = 1
    while True:
        for h in data.get("hits") or []:
            v = h.get("vcf") or {}
            c = str(h.get("chrom") or "").upper().replace("CHR", "")
            af = ((h.get("gnomad_exome") or {}).get("af") or {}).get("af")
            if c not in CHROMS or not v.get("position") or not v.get("ref") or not v.get("alt") or af is None:
                continue
            if isinstance(af, list):
                af = max(af)
            sites.append((c, int(v["position"]), str(v["ref"]).upper(), str(v["alt"]).upper(), float(af)))
        if len(sites) >= want or not data.get("_scroll_id") or not data.get("hits"):
            break
        resp = get_json(BASE, source="MyVariant.info", params={"scroll_id": data["_scroll_id"], "assembly": "hg38"},
                        cache_ttl=0, timeout=120)
        data = resp.json()
        pages += 1
    log(f"{rng}: {len(sites)} sites in {pages} page(s)")
    return sites[:want]


def gt_field(g: int, rnd: random.Random, mosaic: bool = False) -> str:
    dp = rnd.randint(18, 70)
    gq = rnd.randint(40, 99)
    if mosaic:
        alt = max(3, round(dp * 0.15))
        return f"0/1:{dp - alt},{alt}:{dp}:{gq}"
    if g == 0:
        return f"0/0:{dp},0:{dp}:{gq}"
    if g == 2:
        return f"1/1:0,{dp}:{dp}:{gq}"
    alt = max(1, min(dp - 1, round(dp * rnd.uniform(0.35, 0.62))))
    return f"0/1:{dp - alt},{alt}:{dp}:{gq}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    rnd = random.Random(args.seed)
    log = lambda m: print(m, file=sys.stderr)  # noqa: E731
    retrieved = now_iso()
    bands = {name: fetch_band(rng, want, log) for name, rng, want in BANDS}
    digest = hashlib.sha256(json.dumps(bands, sort_keys=True).encode()).hexdigest()
    rows = {}
    for name, sites in bands.items():
        for c, pos, ref, alt, af in sites:
            key = (c, pos, ref, alt)
            if key in rows:
                continue
            if name in ("common", "uncommon"):
                mom = sum(rnd.random() < af for _ in range(2))
                dad = sum(rnd.random() < af for _ in range(2))
            else:
                mom, dad = (1, 0) if rnd.random() < 0.5 else (0, 1)
            rows[key] = [mom, dad]
    # novel: a different ALT base at real rare-band SNV sites, private to one parent
    novel = 0
    for c, pos, ref, alt, _ in bands["rare"][:600]:
        if len(ref) != 1 or len(alt) != 1:
            continue
        other = [b for b in "ACGT" if b not in (ref, alt)]
        key = (c, pos, ref, rnd.choice(other))
        if key not in rows:
            rows[key] = [1, 0] if rnd.random() < 0.5 else [0, 1]
            novel += 1
        if novel >= 300:
            break

    def child(m: int, f: int) -> int:
        def allele(g: int) -> int:
            return 1 if g == 2 else (0 if g == 0 else int(rnd.random() < 0.5))
        x = allele(m) + allele(f)
        return x

    lines = []
    for (c, pos, ref, alt), (m, f) in rows.items():
        p = child(m, f)
        if c == "X":
            f = 2 if f else 0  # the father's X is hemizygous, written diploid as GATK does
            p = (1 if m == 2 else (int(rnd.random() < 0.5) if m == 1 else 0)) + (1 if f else 0)
        if m == f == p == 0:
            continue
        lines.append((c, pos, ref, alt, [gt_field(m, rnd), gt_field(f, rnd), gt_field(p, rnd)]))
    planted = {(c, pos, ref, alt) for c, pos, ref, alt, _, _ in PLANTED}
    lines = [ln for ln in lines if (ln[0], ln[1], ln[2], ln[3]) not in planted]
    for c, pos, ref, alt, gts, _ in PLANTED:
        fields = []
        for g in gts:
            if g.endswith(":mosaic"):
                fields.append(gt_field(1, rnd, mosaic=True))
            else:
                fields.append(gt_field({"0/0": 0, "0/1": 1, "1/1": 2}[g], rnd))
        lines.append((c, pos, ref, alt, fields))
    order = {c: i for i, c in enumerate(CHROMS)}
    lines.sort(key=lambda r: (order[r[0]], r[1], r[2], r[3]))
    seen = set()
    body = []
    for c, pos, ref, alt, fields in lines:
        if (c, pos, ref, alt) in seen:
            continue
        seen.add((c, pos, ref, alt))
        body.append(f"chr{c}\t{pos}\t.\t{ref}\t{alt}\t100\tPASS\t.\tGT:AD:DP:GQ\t" + "\t".join(fields))
    from zebra.vcf import CHROM_LENGTHS

    head = ["##fileformat=VCFv4.2",
            "##source=zebra-mod make_exome_scale.py: REAL GRCh38 variant sites, INVENTED genotypes (not a person)",
            f"##zebra_sites=MyVariant.info {BASE} q=gnomad_exome.af.af bands "
            + ";".join(f"{n}:{r}:{w}" for n, r, w in BANDS) + f" retrieved {retrieved} sha256 {digest}",
            f"##zebra_counts=sites {sum(len(v) for v in bands.values())}; novel {novel}; planted {len(PLANTED)}; "
            f"records {len(body)}; seed {args.seed}",
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">']
    head += [f"##contig=<ID=chr{c},length={CHROM_LENGTHS['GRCh38'][c]}>" for c in CHROMS + ["Y"]]
    head.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tM\tF\tP")
    with open(args.out, "w") as fh:
        fh.write("\n".join(head + body) + "\n")
    print(json.dumps({"out": args.out, "records": len(body), "novel": novel, "sha256_sites": digest,
                      "bands": {n: len(v) for n, v in bands.items()}}))


if __name__ == "__main__":
    main()
