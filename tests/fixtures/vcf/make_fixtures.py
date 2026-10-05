"""Build the synthetic family VCFs used by tests/test_vcf.py (and optionally re-capture Ensembl responses).

    python tests/fixtures/vcf/make_fixtures.py            # write trio.vcf(.gz), family2.vcf(.gz), bad.vcf, genes.txt
    python tests/fixtures/vcf/make_fixtures.py --capture  # also refresh vep_trio.json, vep_family2.json,
                                                          # lookup_genes.json (network)

Two families, each written once and never rewritten in place: `trio.vcf` is the
original acceptance family, `family2.vcf` the second one added for the triage
defects found in review (common pathogenic recessive allele, de novo plus
inherited pair, diploid male X). `bad.vcf` holds data lines that cannot be read.

Coordinates are real GRCh38 positions; every REF base was checked against Ensembl
/sequence/region (2026-10-05). The variants are real dbSNP/ClinVar alleles
(resolved with Ensembl variant_recoder / overlap); only the genotypes, depths and
the family are invented. Samples: M (mother), F (father), P (proband, female),
S (brother, male). trio.vcf uses chr-prefixed contigs, trio.vcf.gz plain names
(written as BGZF, so it is what bgzip would produce).
"""

from __future__ import annotations

import json
import os
import struct
import sys
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLES = ("M", "F", "P", "S")

# GRCh38 lengths (Ensembl /info/assembly/homo_sapiens, GRCh38.p14)
CONTIGS = [("1", 248956422), ("2", 242193529), ("3", 198295559), ("4", 190214555), ("5", 181538259),
           ("6", 170805979), ("7", 159345973), ("8", 145138636), ("9", 138394717), ("10", 133797422),
           ("11", 135086622), ("12", 133275309), ("13", 114364328), ("14", 107043718), ("15", 101991189),
           ("16", 90338345), ("17", 83257441), ("18", 80373285), ("19", 58617616), ("20", 64444167),
           ("21", 46709983), ("22", 50818468), ("X", 156040895), ("Y", 57227415), ("M", 16569)]

# chrom, pos, id, ref, alts, filter, {sample: (GT, AD, DP, GQ)}, what it tests
RECORDS = [
    ("2", 166036097, "rs2105793395", "A", "C", "PASS",
     {"M": ("0/0", "30,0", 30, 90), "F": ("0/0", "28,0", 28, 84), "P": ("0/1", "3,3", 6, 40), "S": ("0/0", "25,0", 25, 75)},
     "SCN1A ClinVar pathogenic stop: proband DP 6 < 10 -> removed by the quality filter"),
    ("2", 166036116, "rs796052993", "C", "A", "LowQual",
     {"M": ("0/0", "30,0", 30, 90), "F": ("0/0", "26,0", 26, 78), "P": ("0/1", "10,4", 14, 5), "S": ("0/0", "22,0", 22, 66)},
     "SCN1A ClinVar pathogenic stop: FILTER LowQual, GQ 5 -> removed"),
    ("2", 166042334, "rs794726730", "G", "A", "PASS",
     {"M": ("0/0", "31,0", 31, 90), "F": ("0/0", "29,0", 29, 87), "P": ("0/1", "18,17", 35, 99), "S": ("0/0", "25,0", 25, 75)},
     "SCN1A NM_001165963.4:c.2134C>T p.Arg712* de novo in P"),
    ("2", 166053034, "rs3812718", "C", "T,A", "PASS",
     {"M": ("0/1", "15,16,0", 31, 99), "F": ("0/0", "30,0,0", 30, 90), "P": ("0/1", "20,18,0", 38, 99), "S": ("1/1", "0,22,0", 22, 66)},
     "multi-allelic common SCN1A intronic SNP: T common (filtered by AF), A not carried by P"),
    ("2", 166122240, "rs7587026", "C", "A", "PASS",
     {"M": ("0/0", "27,0", 27, 81), "F": ("0/1", "13,14", 27, 99), "P": ("0|1", "14,15", 29, 99), "S": ("0/0", "24,0", 24, 72)},
     "common SCN1A intronic SNP, paternal, phased GT -> filtered by AF"),
    ("2", 178560007, "rs562860372", "C", "T", "PASS",
     {"M": ("0/0", "30,0", 30, 90), "F": ("0/1", "15,15", 30, 99), "P": ("0/1", "16,16", 32, 99), "S": ("0/0", "25,0", 25, 75)},
     "TTN missense VUS, paternal (outside the gene list)"),
    ("2", 178560163, "rs72648224", "T", "A", "PASS",
     {"M": ("0/1", "14,15", 29, 99), "F": ("0/0", "28,0", 28, 84), "P": ("0/1", "17,16", 33, 99), "S": ("0/1", "12,13", 25, 99)},
     "TTN stop, maternal (outside the gene list)"),
    ("3", 25751134, "rs200561967", "G", "A", "PASS",
     {"M": ("0/0", "4,0", 4, 12), "F": ("0/0", "33,0", 33, 99), "P": ("0/1", "20,20", 40, 99), "S": ("0/0", "21,0", 21, 63)},
     "NGLY1 stop; mother DP 4 -> possible_de_novo, not de novo"),
    ("5", 126554292, ".", "C", "G", "PASS",
     {"M": ("0/1", "14,14", 28, 99), "F": ("0/1", "15,13", 28, 99), "P": ("1/1", "0,30", 30, 90), "S": ("0/1", "12,12", 24, 99)},
     "ALDH7A1 NM_001182.5:c.1195G>C p.Gly399Arg homozygous in P, absent from gnomAD"),
    ("7", 117559479, "rs213950", "G", "A", "PASS",
     {"M": ("0/1", "15,15", 30, 99), "F": ("0/1", "14,16", 30, 99), "P": ("1/1", "0,31", 31, 93), "S": ("0/1", "13,13", 26, 99)},
     "CFTR c.1408G>A p.Val470Met common, homozygous in P -> filtered by AF"),
    ("7", 117587806, "rs75527207", "G", "A", "PASS",
     {"M": ("0/1", "16,15", 31, 99), "F": ("0/0", "30,0", 30, 90), "P": ("0/1", "15,17", 32, 99), "S": ("0/0", "27,0", 27, 81)},
     "CFTR c.1652G>A p.Gly551Asp maternal -> comp-het partner"),
    ("7", 117652877, "rs80034486", "C", "G", "PASS",
     {"M": ("0/0", "29,0", 29, 87), "F": ("0/1", "14,14", 28, 99), "P": ("0|1", "16,14", 30, 99), "S": ("0/1", "11,12", 23, 99)},
     "CFTR c.3909C>G p.Asn1303Lys paternal -> comp-het partner"),
    ("X", 31178721, "rs398123832", "G", "A", "PASS",
     {"M": ("0/1", "15,14", 29, 99), "F": ("0", "20,0", 20, 60), "P": ("0/1", "13,14", 27, 99), "S": ("1", "0,24", 24, 72)},
     "DMD c.10171C>T p.Arg3391*: maternal het in P (female), hemizygous in S (male)"),
]

GENES = ["SCN1A", "ALDH7A1", "CFTR", "DMD", "NGLY1", "PCDH19"]

# ---------------------------------------------------------------- family 2
# A second family, so the first one's expectations never move. Samples: M2
# (mother), F2 (father), P2 (proband, male). Every REF base was checked against
# Ensembl /sequence/region on 2026-10-06; the Y SNVs are real dbSNP alleles
# found with Ensembl /overlap/region/human/Y?feature=variation. The caller
# emitted diploid X and Y for every sample, which is the GATK default when the
# ploidy is not set — that is what makes the male X call read as a homozygote.
SAMPLES2 = ("M2", "F2", "P2")
RECORDS2 = [
    ("7", 117559590, "rs113993960", "ATCT", "A", "PASS",
     {"M2": ("0/1", "17,16", 33, 99), "F2": ("0/0", "31,0", 31, 93), "P2": ("0/1", "18,17", 35, 99)},
     "CFTR c.1521_1523del p.Phe508del, maternal: gnomAD grpmax 0.0149 is above the 1 % default, and ClinVar "
     "reports it pathogenic — the AF filter must not delete it, and the significance must attach to an indel"),
    ("7", 117587806, "rs75527207", "G", "A", "PASS",
     {"M2": ("0/0", "30,0", 30, 90), "F2": ("0/0", "29,0", 29, 87), "P2": ("0/1", "16,15", 31, 99)},
     "CFTR c.1652G>A p.Gly551Asp de novo: the comp-het partner of F508del, so it must not be held to the "
     "dominant AF cut-off"),
    ("X", 31178721, "rs398123832", "G", "A", "PASS",
     {"M2": ("0/1", "15,14", 29, 99), "F2": ("0/0", "26,0", 26, 78), "P2": ("1/1", "0,27", 27, 81)},
     "DMD c.10171C>T p.Arg3391* called 1/1 on male X outside the PARs: hemizygous, not a homozygote"),
    ("Y", 2786042, "rs2051121937", "T", "G", "PASS",
     {"M2": ("./.", "0,0", 0, 0), "F2": ("1/1", "0,22", 22, 66), "P2": ("1/1", "0,24", 24, 72)},
     "real Y non-PAR SNV: evidence that the proband is male"),
    ("Y", 2786191, "rs936999469", "C", "T", "PASS",
     {"M2": ("./.", "0,0", 0, 0), "F2": ("1/1", "0,20", 20, 60), "P2": ("1/1", "0,21", 21, 63)},
     "second real Y non-PAR SNV: two Y ALT calls settle the sex inference"),
]

# Data lines that cannot be read, to prove they are counted and reported.
BAD_LINES = [
    "chr1 1000 . A G 500 PASS . GT:DP:GQ 0/1:30:90",          # space-separated
    "chr1\t3e2\t.\tA\tG\t500\tPASS\t.\tGT:DP:GQ\t0/1:30:90",  # POS not an integer
    "chr1\t4000\t.\tA\tG\t500\tPASS\t.",                      # 8 columns, no genotypes (readable)
    "chr1\t5000\t.\tA\tG\t500\tPASS",                         # 7 columns
    "chr1\t0\t.\tA\tG\t500\tPASS\t.\tGT:DP:GQ\t0/1:30:90",    # POS 0 is not a position
    "chr1\t6000\t.\tA\t\t500\tPASS\t.\tGT:DP:GQ\t0/1:30:90",  # empty ALT
    "chr1\t7000\t.\tA\tG\t500\tPASS\t.\tGT:DP:GQ\t0/1:30:90", # readable
]


def vcf_text(chr_prefix: bool, records=None, samples=None) -> str:
    records = RECORDS if records is None else records
    samples = SAMPLES if samples is None else samples
    pre = "chr" if chr_prefix else ""
    lines = [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=LowQual,Description="Low quality">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele">',
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency in this call set (not a population frequency)">',
    ]
    for name, length in CONTIGS:
        cid = (pre + name) if chr_prefix else ("MT" if name == "M" else name)
        lines.append(f"##contig=<ID={cid},length={length}>")
    lines.append("##reference=file:///references/GRCh38_full_analysis_set_plus_decoy_hla.fa")
    lines.append("##source=zebra-mod synthetic fixture: real GRCh38 alleles, invented genotypes")
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))
    for chrom, pos, vid, ref, alts, filt, gts, _ in records:
        n_alt = len(alts.split(","))
        ac = [0] * n_alt
        an = 0
        for s in samples:
            for a in gts[s][0].replace("|", "/").split("/"):
                if a not in (".",):
                    an += 1
                    if a != "0":
                        ac[int(a) - 1] += 1
        info = f"AC={','.join(map(str, ac))};AN={an};AF={','.join(f'{x / an:.3f}' if an else '0.000' for x in ac)}"
        cols = [pre + chrom, str(pos), vid, ref, alts, "500" if filt == "PASS" else "12", filt, info, "GT:AD:DP:GQ"]
        cols += [f"{g}:{ad}:{dp}:{gq}" for g, ad, dp, gq in (gts[s] for s in samples)]
        lines.append("\t".join(cols))
    return "\n".join(lines) + "\n"


def bgzf(data: bytes) -> bytes:
    """BGZF (blocked gzip, as written by bgzip): gzip members of <=64 KB with a 'BC' extra field, plus EOF block."""
    out = b""
    for i in range(0, len(data), 65280):
        block = data[i:i + 65280]
        comp = zlib.compressobj(9, zlib.DEFLATED, -15)
        cdata = comp.compress(block) + comp.flush()
        header = (b"\x1f\x8b\x08\x04" + b"\x00\x00\x00\x00" + b"\x00\xff" + struct.pack("<H", 6) + b"BC"
                  + struct.pack("<H", 2) + struct.pack("<H", len(cdata) + 25))
        out += header + cdata + struct.pack("<II", zlib.crc32(block) & 0xFFFFFFFF, len(block))
    return out + bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def alleles(records=None):
    out = []
    for chrom, pos, _, ref, alts, _, _, _ in (RECORDS if records is None else records):
        for alt in alts.split(","):
            out.append((chrom, pos, ref, alt))
    return out


def _trim_vep(rec):
    keep_tc = ("transcript_id", "gene_id", "gene_symbol", "mane_select", "canonical", "biotype", "impact",
               "consequence_terms", "hgvsc", "hgvsp", "revel", "alphamissense", "cadd_phred", "spliceai")
    keep_cv = ("id", "allele_string", "frequencies", "clin_sig", "clin_sig_allele", "start", "end")
    out = {k: rec.get(k) for k in ("input", "assembly_name", "seq_region_name", "start", "end", "allele_string",
                                   "most_severe_consequence") if k in rec}
    rank = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
    tcs = rec.get("transcript_consequences") or []
    main = [t for t in tcs if t.get("mane_select") or t.get("canonical")]
    floor = max((rank.get(t.get("impact"), 0) for t in main), default=0)
    tcs = main + [t for t in tcs if t not in main and t.get("biotype") == "protein_coding"
                  and rank.get(t.get("impact"), 0) > floor][:3]  # keep 'more severe elsewhere' cases
    out["transcript_consequences"] = [{k: t[k] for k in keep_tc if k in t} for t in tcs]
    out["colocated_variants"] = [{k: c[k] for k in keep_cv if k in c} for c in rec.get("colocated_variants") or []]
    if rec.get("intergenic_consequences"):
        out["intergenic_consequences"] = rec["intergenic_consequences"]
    return out


def capture() -> None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
    from zebra.http import post_json
    from zebra.sources import ensembl

    got = ensembl.vep_batch(alleles(), "GRCh38")
    with open(os.path.join(HERE, "vep_trio.json"), "w") as fh:
        json.dump({"_source": got.sources, "records": [_trim_vep(r) for r in got.result]}, fh, indent=1)
    resp = post_json("https://rest.ensembl.org/lookup/symbol/homo_sapiens", {"symbols": GENES + ["TTN", "NOTAGENE1"]},
                     source="Ensembl lookup", cache_ttl=0, timeout=120)
    keep = ("id", "display_name", "seq_region_name", "start", "end", "strand", "biotype", "assembly_name")
    data = {k: {f: v[f] for f in keep if f in v} for k, v in resp.json().items()}
    with open(os.path.join(HERE, "lookup_genes.json"), "w") as fh:
        json.dump({"_source": {"url": resp.url, "retrieved_at": resp.retrieved_at}, "response": data}, fh, indent=1)
    print(f"captured VEP for {len(got.result)} alleles and lookup for {len(data)} genes")

    got2 = ensembl.vep_batch(alleles(RECORDS2), "GRCh38")
    with open(os.path.join(HERE, "vep_family2.json"), "w") as fh:
        json.dump({"_source": got2.sources, "records": [_trim_vep(r) for r in got2.result]}, fh, indent=1)
    print(f"captured VEP for {len(got2.result)} family-2 alleles")
    capture_cnv()


# Responses for the `zebra cnv` tests. The region is a 300 kb window inside
# 15q11.2-q13.1 that contains UBE3A (ClinGen haploinsufficiency score 3), and
# the transcript is DMD's Ensembl canonical one (the exon numbering NM_004006
# uses), trimmed to the fields zebra reads.
CNV_REGION = ("15", 25200000, 25500000)
CNV_EXON_GENE = "DMD"


def capture_cnv() -> None:
    from zebra.http import get_json
    from zebra.sources import clingen, ensembl

    chrom, start, end = CNV_REGION
    resp = get_json(f"https://rest.ensembl.org/overlap/region/human/{chrom}:{start}-{end}",
                    source="Ensembl overlap", params={"feature": "gene"}, cache_ttl=0, timeout=120)
    keep = ("id", "external_name", "biotype", "start", "end", "strand", "seq_region_name", "assembly_name")
    genes = [{k: g[k] for k in keep if k in g} for g in resp.json()]
    with open(os.path.join(HERE, "overlap_cnv.json"), "w") as fh:
        json.dump({"_source": {"url": resp.url, "retrieved_at": resp.retrieved_at},
                   "region": f"{chrom}:{start}-{end}", "response": genes}, fh, indent=1)

    # the real ClinGen bulk TSV, trimmed to its header and the region's genes:
    # that is what zebra fetches once and parses per gene
    from zebra.http import request as http_request

    wanted = sorted({g.get("external_name") for g in genes if g.get("external_name")
                     and g.get("biotype") == "protein_coding"})
    bulk = http_request(clingen.DOSAGE_TSV["GRCh38"], source="ClinGen dosage",
                        accept="text/tab-separated-values,text/plain,*/*", cache_ttl=0, timeout=120)
    kept = [line for line in bulk.text.splitlines()
            if line.startswith("#Gene Symbol") or line.split("\t")[0].strip() in wanted]
    rows = {gene: clingen.parse_dosage_tsv("\n".join(kept), gene) for gene in wanted}
    with open(os.path.join(HERE, "dosage_cnv.json"), "w") as fh:
        json.dump({"_source": {"url": bulk.url, "retrieved_at": bulk.retrieved_at,
                               "note": "ClinGen dosage bulk TSV, header plus the rows for the region's genes"},
                   "tsv": "\n".join(kept) + "\n", "response": rows}, fh, indent=1)

    data = ensembl.lookup_symbol(CNV_EXON_GENE, "GRCh38", expand=True).result
    tx = [t for t in data["Transcript"] if t.get("is_canonical")][0]
    trimmed = {k: data[k] for k in ("id", "display_name", "seq_region_name", "start", "end", "strand", "biotype",
                                    "assembly_name") if k in data}
    trimmed["Transcript"] = [{
        **{k: tx[k] for k in ("id", "display_name", "biotype", "is_canonical", "length", "strand",
                              "seq_region_name", "start", "end") if k in tx},
        "Exon": [{k: e[k] for k in ("id", "start", "end", "strand", "seq_region_name") if k in e}
                 for e in tx["Exon"]],
    }]
    with open(os.path.join(HERE, f"lookup_{CNV_EXON_GENE.lower()}.json"), "w") as fh:
        json.dump({"_source": "Ensembl /lookup/symbol/homo_sapiens/DMD?expand=1, canonical transcript only",
                   "response": trimmed}, fh, indent=1)
    print(f"captured {len(genes)} genes, {len(rows)} ClinGen dosage rows and {CNV_EXON_GENE}'s canonical transcript")


def main() -> None:
    with open(os.path.join(HERE, "trio.vcf"), "w") as fh:
        fh.write(vcf_text(chr_prefix=True))
    with open(os.path.join(HERE, "trio.vcf.gz"), "wb") as fh:
        fh.write(bgzf(vcf_text(chr_prefix=False).encode("utf-8")))
    with open(os.path.join(HERE, "family2.vcf"), "w") as fh:
        fh.write(vcf_text(chr_prefix=True, records=RECORDS2, samples=SAMPLES2))
    with open(os.path.join(HERE, "family2.vcf.gz"), "wb") as fh:
        fh.write(bgzf(vcf_text(chr_prefix=False, records=RECORDS2, samples=SAMPLES2).encode("utf-8")))
    with open(os.path.join(HERE, "bad.vcf"), "w") as fh:
        head = vcf_text(chr_prefix=True, records=[], samples=("P",)).rstrip("\n")
        fh.write(head + "\n" + "\n".join(BAD_LINES) + "\n")
    with open(os.path.join(HERE, "genes.txt"), "w") as fh:
        fh.write("# genes for the triage acceptance run\n" + "\n".join(GENES) + "\n")
    if "--capture" in sys.argv:
        capture()


if __name__ == "__main__":
    main()
