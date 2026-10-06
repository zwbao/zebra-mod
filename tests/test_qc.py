"""zebra qc: sex, KING-robust kinship, ROH, Mendelian errors / UPD, mosaic de novo — on synthetic families.

The families are simulated here (seeded): founder haplotypes drawn under Hardy-Weinberg at invented allele
frequencies on evenly spaced autosomal positions, children made by transmitting one recombined haplotype from
each parent, then the planted features. No real person's genotypes are involved.

`zebra.qc` is new in v0.2 (CP1-15): on the old code every test here fails at import.
"""

from __future__ import annotations

import json
import random
from typing import Dict, List, Optional, Tuple

import pytest

from zebra import qc as Q
from zebra import vcf as V
from zebra.core import UsageError

LENGTHS = V.CHROM_LENGTHS["GRCh38"]
AUTO = [str(i) for i in range(1, 23)]


class Sim:
    """Haplotypes per individual and chromosome; genotypes written as a joint-called VCF."""

    def __init__(self, seed: int = 7, sites_per_chrom: int = 600, x_sites: int = 120, error: float = 0.001):
        self.r = random.Random(seed)
        self.error = error
        self.sites: Dict[str, List[Tuple[int, float]]] = {}
        for c in AUTO + ["X"]:
            n = x_sites if c == "X" else sites_per_chrom
            lo, hi = (2_800_000, 155_000_000) if c == "X" else (100_000, LENGTHS[c] - 100_000)
            step = (hi - lo) // n
            self.sites[c] = [(lo + i * step + self.r.randrange(step // 2), self.r.uniform(0.08, 0.5))
                             for i in range(n)]
        self.haps: Dict[str, Dict[str, Tuple[List[int], List[int]]]] = {}
        self.sex: Dict[str, str] = {}
        self.extra: List[Tuple[str, int, str, str, Dict[str, str]]] = []  # planted literal rows

    def founder(self, name: str, sex: str) -> None:
        self.sex[name] = sex
        h = {}
        for c, sites in self.sites.items():
            a = [1 if self.r.random() < p else 0 for _, p in sites]
            b = [1 if self.r.random() < p else 0 for _, p in sites]
            if c == "X" and sex == "male":
                b = list(a)  # one X: written as a homozygous diploid call, as GATK does by default
            h[c] = (a, b)
        self.haps[name] = h

    def _gamete(self, parent: str, c: str) -> List[int]:
        a, b = self.haps[parent][c]
        if c == "X" and self.sex[parent] == "male":
            return list(a)
        cut = self.r.randrange(len(a))
        first, second = (a, b) if self.r.random() < 0.5 else (b, a)
        return first[:cut] + second[cut:]

    def child(self, name: str, mother: str, father: str, sex: str) -> None:
        self.sex[name] = sex
        h = {}
        for c in self.sites:
            m = self._gamete(mother, c)
            if c == "X" and sex == "male":
                h[c] = (m, list(m))
            else:
                h[c] = (m, self._gamete(father, c))
        self.haps[name] = h

    def genotypes(self, name: str, c: str) -> List[int]:
        a, b = self.haps[name][c]
        return [x + y for x, y in zip(a, b)]

    def write(self, path, samples: List[str]) -> str:
        rows = []
        for c in AUTO + ["X"]:
            gts = {s: self.genotypes(s, c) for s in samples}
            for i, (pos, _) in enumerate(self.sites[c]):
                calls = []
                for s in samples:
                    g = gts[s][i]
                    if self.r.random() < self.error:  # genotyping error
                        g = (g + self.r.choice((1, 2))) % 3
                    calls.append(g)
                if not any(calls):
                    continue  # joint calling writes only sites where someone carries ALT
                rows.append((c, pos, "A", "G", [_fmt(g) for g in calls]))
        for c, pos, ref, alt, by in self.extra:
            rows.append((c, pos, ref, alt, [by.get(s, "0/0:30,0:30:99") for s in samples]))
        rows.sort(key=lambda r: (AUTO.index(r[0]) if r[0] in AUTO else 99, r[1]))
        head = ["##fileformat=VCFv4.2", "##source=zebra-mod tests: simulated genotypes, invented allele frequencies"]
        head += [f"##contig=<ID=chr{c},length={LENGTHS[c]}>" for c in AUTO + ["X", "Y"]]
        head.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))
        body = [f"chr{c}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT:AD:DP:GQ\t" + "\t".join(f) for c, pos, ref, alt, f in rows]
        path.write_text("\n".join(head + body) + "\n")
        return str(path)


def _fmt(g: int) -> str:
    return {0: "0/0:30,0:30:99", 1: "0/1:15,15:30:99", 2: "1/1:0,30:30:99"}[g]


def _ped(path, rows: List[str]) -> str:
    path.write_text("\n".join(rows) + "\n")
    return str(path)


@pytest.fixture(scope="module")
def family(tmp_path_factory):
    return build_family(tmp_path_factory.mktemp("qcfam"))


def build_family(tmp):
    """Trio M/F/P with a full sibling S; the parents share one ancestral haplotype on chr3 (consanguinity-like
    40 Mb ROH in P); P has maternal isodisomy of chr15 and maternal heterodisomy of chr7; two de novo SNVs in P,
    one at a mosaic allele fraction."""
    sim = Sim(seed=11)
    sim.founder("M", "female")
    sim.founder("F", "male")
    sim.founder("U", "male")  # an unrelated man, for the swap test
    # shared ancestral haplotype on chr3 50-90 Mb in both parents
    s3 = sim.sites["3"]
    region = [i for i, (pos, _) in enumerate(s3) if 50_000_000 <= pos <= 90_000_000]
    for i in region:
        sim.haps["F"]["3"][0][i] = sim.haps["M"]["3"][0][i]
    sim.child("P", "M", "F", "female")
    sim.child("S", "M", "F", "male")
    for i in region:  # P inherits the shared haplotype from both parents there
        sim.haps["P"]["3"][0][i] = sim.haps["M"]["3"][0][i]
        sim.haps["P"]["3"][1][i] = sim.haps["F"]["3"][0][i]
    # UPD: chr15 maternal isodisomy (one maternal haplotype twice), chr7 maternal heterodisomy (both)
    m15 = sim.haps["M"]["15"][0]
    sim.haps["P"]["15"] = (list(m15), list(m15))
    sim.haps["P"]["7"] = (list(sim.haps["M"]["7"][0]), list(sim.haps["M"]["7"][1]))
    sim.extra.append(("12", 40_000_123, "C", "T", {"P": "0/1:40,8:48:99"}))   # mosaic de novo, VAF 0.17
    sim.extra.append(("12", 60_000_321, "G", "A", {"P": "0/1:20,22:42:99"}))  # germline de novo
    vcf = sim.write(tmp / "family.vcf", ["M", "F", "P", "S", "U"])
    ped = _ped(tmp / "family.ped", ["fam1 M 0 0 2 1", "fam1 F 0 0 1 1", "fam1 P F M 2 2", "fam1 S F M 1 1"])
    return {"vcf": vcf, "ped": ped, "tmp": tmp, "sim": sim}


@pytest.fixture(scope="module")
def qc_result(family):
    return Q.run(family["vcf"], ped=family["ped"])


# ------------------------------------------------------------- statistics

def test_cp1_15_king_robust_formula_and_cutoffs():
    # algebra check: phi = 1/2 - sum (x_i - x_j)^2 / (4 N_Aa(i)) with i the sample with fewer hets
    hh, opp, ha, hb = 300, 10, 500, 520
    sq = (ha + hb - 2 * hh) + 4 * opp
    assert Q.king_robust(hh, opp, ha, hb) == pytest.approx(0.5 - sq / (4 * ha))
    assert Q.degree(0.36) == "duplicate/MZ twin" and Q.degree(0.25) == "1st degree"
    assert Q.degree(0.12) == "2nd degree" and Q.degree(0.06) == "3rd degree"
    assert Q.degree(0.01).startswith("unrelated") and Q.king_robust(0, 0, 0, 5) is None


def test_cp1_15_binomial_helpers():
    lo, hi = Q.wilson(8, 48)
    assert 0.08 < lo < 0.10 and 0.29 < hi < 0.31  # 8/48 = 0.167
    assert Q.binom_cdf(5, 10, 0.5) == pytest.approx(0.623046875)
    assert Q.poisson_sf(0, 2.0) == 1.0 and Q.poisson_sf(3, 1.0) == pytest.approx(1 - 2.5 * 2.718281828 ** -1, rel=1e-6)
    assert Q.poisson_sf(30, 2.0) < 1e-20
    m = Q.mosaic_assessment(40, 8)
    assert m["possible_mosaic"] is True and m["vaf"] == 0.167
    assert Q.mosaic_assessment(20, 22)["possible_mosaic"] is False and Q.mosaic_assessment(None, 3) is None


def test_cp1_15_mendelian_table():
    assert Q.ME_TABLE[(1, 0, 0)] == "de_novo_like" and Q.ME_TABLE[(1, 1, 0)] == "ok"
    assert Q.ME_TABLE[(2, 2, 0)] == Q.ME_TABLE[(0, 0, 2)] == Q.ME_TABLE[(2, 1, 0)] == "pat_absent"
    assert Q.ME_TABLE[(0, 2, 0)] == Q.ME_TABLE[(2, 0, 2)] == Q.ME_TABLE[(2, 0, 1)] == "mat_absent"
    assert sum(1 for v in Q.ME_TABLE.values() if v == "ok") == 15


def test_cp1_15_find_roh_on_a_planted_run():
    pos = list(range(1_000_000, 61_000_000, 100_000))  # 600 sites
    hets = [1 if i % 2 else 0 for i in range(600)]
    for i in range(200, 400):  # 20 Mb of homozygous ALT calls with one het error inside
        hets[i] = 0
    hets[300] = 1
    segs = Q.find_roh(pos, hets)
    assert len(segs) == 1
    # the one-het allowance lets a run take in the edge sites of a window (here 2 sites each side)
    assert pos[196] <= segs[0]["start"] <= pos[200] and pos[399] <= segs[0]["end"] <= pos[403]
    assert 19.5 <= segs[0]["mb"] <= 20.8 and segs[0]["hets"] <= 3
    assert Q.find_roh(pos[:10], hets[:10]) == []


# ------------------------------------------------------------- the family

def test_cp1_15_kinship_ped_check_and_sex(qc_result):
    res = qc_result.result
    pairs = {(p["a"], p["b"]): p for p in res["relatedness"]}
    po = pairs[("M", "P")] if ("M", "P") in pairs else pairs[("P", "M")]
    assert 0.18 < po["kinship"] < 0.32 and po["relationship"].endswith("parent-offspring-like")
    sib = pairs.get(("P", "S")) or pairs[("S", "P")]
    assert 0.18 < sib["kinship"] < 0.36 and sib["relationship"].endswith("full-sibling-like")
    mf = pairs[("M", "F")]
    assert mf["kinship"] < 0.0442 and mf["expected"] == "unrelated (parents of one child)" and mf["problem"] is None
    # P carries maternal UPD of chr7 and chr15, so her kinship to F is lower, but still 1st degree
    fp = pairs.get(("F", "P")) or pairs[("P", "F")]
    assert fp["relationship"].startswith("1st degree") and fp["problem"] is None
    assert res["sex"]["M"]["inferred"] == "female" and res["sex"]["F"]["inferred"] == "male"
    assert res["sex"]["S"]["agrees"] is True and res["sex"]["P"]["stated"] == "female"
    assert res["methods"]["relatedness"]["source"].startswith("Manichaikul A")


def test_cp1_15_roh_found_where_planted(qc_result):
    roh = qc_result.result["roh"]["P"]
    chr3 = [s for s in roh["segments"] if s["chrom"] == "3"]
    assert chr3 and any(s["start"] < 55_000_000 and s["end"] > 85_000_000 for s in chr3)
    chr15 = [s for s in roh["segments"] if s["chrom"] == "15"]
    assert chr15 and sum(s["mb"] for s in chr15) > 80  # isodisomy: the whole chromosome
    assert roh["chromosomes"]["chr15"]["covered_span_fraction"] > 0.8
    assert not [s for s in roh["segments"] if s["chrom"] == "7"]  # heterodisomy leaves no run
    # the isodisomic chromosome is taken out of the parental-relatedness reading
    assert roh["excluding_upd"]["chromosomes"] == ["chr15"] and roh["excluding_upd"]["f_roh"] < 0.02
    assert "left out" in roh["reading"]
    assert "AutoMap" in Q.METHODS["roh"]["source"] and "exome" in Q.METHODS["roh"]["assumptions"]


def test_cp1_15_upd_hetero_and_isodisomy_told_apart(qc_result):
    trio = qc_result.result["mendelian"][0]
    assert trio["trio"] == {"child": "P", "mother": "M", "father": "F"}
    upd = {f["chrom"]: f for f in trio["upd"]}
    assert set(upd) == {"7", "15"}, upd
    assert upd["15"]["origin"] == "maternal" and upd["15"]["type"] == "isodisomy"
    assert upd["15"]["imprinting_disorder"] == "Prader-Willi syndrome"
    assert upd["7"]["origin"] == "maternal" and upd["7"]["type"] == "heterodisomy"
    assert upd["7"]["imprinting_disorder"] == "Silver-Russell syndrome"
    assert trio["genome_wide"] is None
    # the chromosome with planted ROH from shared ancestry is not called UPD (no excess one-sided errors)
    assert "3" not in upd
    assert any("possible maternal UPD of chr15" in w and "isodisomy" in w for w in qc_result.warnings)


def test_cp1_15_mosaic_de_novo_flagged(qc_result):
    dn = qc_result.result["de_novo"][0]
    calls = {c["variant"]: c for c in dn["calls"]}
    mosaic, germline = calls["12-40000123-C-T"], calls["12-60000321-G-A"]
    assert mosaic["child"]["possible_mosaic"] is True and "postzygotic mosaic" in mosaic["flags"][0]
    assert germline["child"]["possible_mosaic"] is False and germline["flags"] == []
    assert dn["possible_mosaic"] >= 1


def test_cp1_15_sample_swap_is_named(family, tmp_path):
    """The PED names U as P's father: KING says unrelated, the Mendelian errors say so genome-wide, no UPD call."""
    ped = _ped(tmp_path / "swap.ped", ["fam1 M 0 0 2 1", "fam1 U 0 0 1 1", "fam1 P U M 2 2"])
    got = Q.run(family["vcf"], ped=ped)
    pairs = {frozenset((p["a"], p["b"])): p for p in got.result["relatedness"]}
    up = pairs[frozenset(("U", "P"))]
    assert up["kinship"] < 0.0442 and "non-paternity or a swapped/mislabelled father sample" in up["problem"]
    trio = got.result["mendelian"][0]
    assert trio["genome_wide"] and any("not the biological father" in g for g in trio["genome_wide"])
    assert not [f for f in trio["upd"] if f["origin"] == "maternal"]
    assert any("non-paternity" in w for w in got.warnings)


def test_cp1_15_sex_mismatch_against_the_ped(family, tmp_path):
    ped = _ped(tmp_path / "sex.ped", ["fam1 M 0 0 2 1", "fam1 F 0 0 1 1", "fam1 P F M 1 2"])  # P written male
    got = Q.run(family["vcf"], ped=ped)
    assert got.result["sex"]["P"]["agrees"] is False
    assert any(w.startswith("sex mismatch for P: stated male (PED), genotypes look female") for w in got.warnings)


def test_cp1_15_trio_without_ped_and_cli(family, capsys):
    from zebra.cli import main

    got = Q.run(family["vcf"], proband="P", mother="M", father="F")
    assert got.result["trios"] == [{"child": "P", "mother": "M", "father": "F"}]
    assert {f["chrom"] for f in got.result["mendelian"][0]["upd"]} == {"7", "15"}
    assert main(["--json", "qc", family["vcf"], "--ped", family["ped"]]) == 0
    env = json.loads(capsys.readouterr().out)
    assert env["ok"] and env["command"] == "qc" and env["result"]["trios"][0]["child"] == "P"
    assert env["sources"][0]["db"] == "local VCF"


def test_cp1_15_ped_parsing_errors(tmp_path):
    with pytest.raises(UsageError, match="6 columns"):
        Q.read_ped(_ped(tmp_path / "a.ped", ["fam P F M 2"]))
    with pytest.raises(UsageError, match="appears twice"):
        Q.read_ped(_ped(tmp_path / "b.ped", ["fam P 0 0 2 2", "fam P 0 0 2 2"]))
    with pytest.raises(UsageError, match="father but has sex 2"):
        Q.read_ped(_ped(tmp_path / "c.ped", ["fam M 0 0 2 1", "fam P M 0 2 2"]))
    ped = Q.read_ped(_ped(tmp_path / "d.ped", ["# comment", "f1 M 0 0 2 1", "f1 F 0 0 1 1", "f1 P F M 2 2",
                                               "f1 S F M 1 -9"]))
    assert ped.trios(["M", "F", "P", "S"]) == [("P", "M", "F"), ("S", "M", "F")]
    assert ped.full_sibs("P") == ["S"] and ped.get("S").affected is None and ped.get("P").affected is True
    assert ped.expected("M", "F") == "unrelated (parents of one child)" and ped.expected("P", "M") == "parent-offspring"


def test_cp1_15_qc_refuses_a_mostly_unreadable_file():
    from pathlib import Path

    bad = str(Path(__file__).parent / "fixtures" / "vcf" / "bad.vcf")
    with pytest.raises(UsageError, match="QC refuses to report on a file it cannot read"):
        Q.run(bad)


# ------------------------------------------------ fixes after the adversarial review (W5)

def _subset(src, dst, keep):
    lines = [l for l in open(src) if l.startswith("#") or keep(l)]
    dst.write_text("".join(lines))
    return str(dst)


def test_cp1_15_review_king_hand_example_and_closed_ranges():
    # 10 sites, genotypes a/b: hethet 3 (sites 1-3), opposite hom 1 (site 4: 0 vs 2), het only in a: 2, only in b: 1
    a = [1, 1, 1, 0, 1, 1, 2, 0, 2, 0]
    b = [1, 1, 1, 2, 0, 2, 2, 1, 0, 0]
    hh = sum(1 for x, y in zip(a, b) if x == y == 1)
    opp = sum(1 for x, y in zip(a, b) if {x, y} == {0, 2})
    ha, hb = a.count(1), b.count(1)
    assert (hh, opp, ha, hb) == (3, 2, 5, 4)
    # by hand: i = b (4 hets): (3 - 2*2)/(2*4) + 1/2 - (5+4)/(4*4) = -0.125 + 0.5 - 0.5625 = -0.1875
    assert Q.king_robust(hh, opp, ha, hb) == pytest.approx(-0.1875)
    assert Q.degree(0.177) == "1st degree" and Q.degree(0.0884) == "2nd degree" and Q.degree(0.0442) == "3rd degree"
    assert Q.degree(0.354) == "1st degree" and Q.degree(0.3541) == "duplicate/MZ twin"


def test_cp1_15_review_kinship_is_sampled_across_all_autosomes(family, monkeypatch):
    """A quota per chromosome, not the first N sites of the file (a sorted WGS reached only chr1)."""
    monkeypatch.setattr(Q, "KING_PER_CHROM", 40)
    got = Q.run(family["vcf"], ped=family["ped"])
    assert got.result["scan"]["king_autosomes"] == 22 and got.result["scan"]["king_sites"] == 22 * 40
    pairs = {frozenset((p["a"], p["b"])): p for p in got.result["relatedness"]}
    assert pairs[frozenset(("M", "F"))]["autosomes"] == 22


def test_cp1_15_review_swapped_parent_labels_withhold_the_parent_of_origin(family):
    got = Q.run(family["vcf"], proband="P", mother="F", father="M")
    assert got.result["sex"]["F"]["agrees"] is False and got.result["sex"]["F"]["stated_by"] == "named as a mother"
    upd = got.result["mendelian"][0]["upd"]
    assert {f["chrom"] for f in upd} == {"7", "15"}
    assert all(f["origin"] == "unknown" and f["imprinting_disorder"] is None for f in upd)
    assert not any("Angelman" in w for w in got.warnings)


def test_cp1_15_review_a_parental_deletion_is_segmental_not_heterodisomy(tmp_path):
    sim = Sim(seed=21, sites_per_chrom=2500)
    sim.founder("M", "female")
    sim.founder("F", "male")
    sim.child("P", "M", "F", "female")
    s15 = sim.sites["15"]
    for i, (pos, _) in enumerate(s15):
        if 23_000_000 <= pos <= 28_500_000:  # the paternal copy is missing: P shows only the maternal allele
            sim.haps["P"]["15"][1][i] = sim.haps["P"]["15"][0][i]
    vcf = sim.write(tmp_path / "del.vcf", ["M", "F", "P"])
    got = Q.run(vcf, proband="P", mother="M", father="F")
    f = next(x for x in got.result["mendelian"][0]["upd"] if x["chrom"] == "15")
    assert f["type"].startswith("segmental") and f["child_roh_fraction"] < 0.2
    lo, hi = (int(x) for x in f["region"].split(":")[1].split("-"))
    assert 22_000_000 < lo < 24_000_000 and 27_500_000 < hi < 29_000_000
    assert any("deletion on the father's chromosome" in w for w in got.warnings)


def test_cp1_15_review_child_sample_swap_is_one_finding(family, tmp_path):
    ped = _ped(tmp_path / "child.ped", ["fam1 M 0 0 2 1", "fam1 F 0 0 1 1", "fam1 U F M 1 2"])
    got = Q.run(family["vcf"], ped=ped)
    gw = got.result["mendelian"][0]["genome_wide"]
    assert len(gw) == 1 and "the child's sample is swapped" in gw[0]
    assert got.result["mendelian"][0]["upd"] == []
    assert got.result["de_novo"][0].get("unreliable")


def test_cp1_15_review_a_cohort_vcf_does_not_hide_a_swapped_father(family, tmp_path):
    """Thousands of sites where only other cohort members carry ALT must not dilute the error rate."""
    rows = [l.rstrip("\n") for l in open(family["vcf"])]
    head = [r for r in rows if r.startswith("#")]
    body = [r for r in rows if not r.startswith("#")]
    extra = [f"chr22\t{20_000_000 + i * 7}\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:GQ\t"
             + "\t".join(["0/0:30,0:30:99"] * 5) + "\t0/1:15,15:30:99" for i in range(20000)]
    head[-1] = head[-1] + "\tR"
    body = [b + "\t0/0:30,0:30:99" for b in body]
    p = tmp_path / "cohort.vcf"
    p.write_text("\n".join(head + body + extra) + "\n")
    ped = _ped(tmp_path / "swap.ped", ["fam1 M 0 0 2 1", "fam1 U 0 0 1 1", "fam1 P U M 2 2"])
    got = Q.run(str(p), ped=ped)
    trio = got.result["mendelian"][0]
    assert trio["genome_wide"] and "not the biological father" in trio["genome_wide"][0]
    assert trio["sites"] < 12_000  # the 20,000 cohort-only sites are not informative for this trio


def test_cp1_15_review_nothing_usable_is_refused_not_reassuring(tmp_path):
    p = tmp_path / "lowq.vcf"
    p.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n"
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP\n"
                 + "".join(f"chr1\t{1000 + i}\t.\tA\tG\t5\tLowQual\t.\tGT\t0/1\n" for i in range(50)))
    with pytest.raises(UsageError, match="no record is a usable biallelic autosomal SNV"):
        Q.run(str(p))


def test_cp1_15_review_proband_trio_survives_the_sample_cap(family):
    got = Q.run(family["vcf"], ped=family["ped"], proband="P", max_samples=3)
    assert got.result["samples"] == ["P", "M", "F"] and got.result["trios"][0]["child"] == "P"
    assert any("raise --max-samples" in w for w in got.warnings)


def test_cp1_15_review_par_bounds_are_used_when_the_build_is_unstated(family, tmp_path):
    """A male's PAR1 heterozygous calls must not count as X non-PAR hets when the header names no build."""
    src = [l for l in open(family["vcf"]) if not l.startswith("##contig")]
    hdr = [l for l in src if l.startswith("#")]
    body = [l for l in src if not l.startswith("#")]
    par = [f"chrX\t{100000 + i * 1000}\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP:GQ\t"
           + "\t".join(["0/1:15,15:30:99"] * 5) + "\n" for i in range(10)]
    p = tmp_path / "nocontig.vcf"
    p.write_text("".join(hdr + par + body))
    got = Q.run(str(p), ped=family["ped"])
    assert got.result["assembly"] is None and got.result["sex"]["F"]["inferred"] == "male"


def test_cp1_15_review_gates_and_flags(tmp_path):
    """DP/GQ gates keep a weak parental call out of the de novo set; parental mosaicism needs a real fraction."""
    head = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=248956422>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP\tM\tF"]
    rows = ["chr1\t1000\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99\t0/0:5,0:5:99\t0/0:30,0:30:99",   # DP 5
            "chr1\t2000\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99\t0/0:30,0:30:5\t0/0:30,0:30:99",   # GQ 5
            "chr1\t3000\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99\t0/0:27,3:30:99\t0/0:30,0:30:99",  # 3/30
            "chr1\t4000\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99\t0/0:198,2:200:99\t0/0:30,0:30:99",  # 1%
            "chr1\t5000\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:GQ\t1:0,30:30:99\t0/0:30,0:30:99\t0/0:30,0:30:99"]      # haploid
    p = tmp_path / "g.vcf"
    p.write_text("\n".join(head + rows) + "\n")
    got = Q.run(str(p), proband="P", mother="M", father="F")
    dn = got.result["de_novo"][0]
    calls = {c["variant"]: c for c in dn["calls"]}
    assert set(calls) == {"1-3000-A-G", "1-4000-A-G"}  # 1000 (DP) and 2000 (GQ) gated, 5000 haploid ignored
    assert any("mother has 3 of 30 reads" in f for f in calls["1-3000-A-G"]["flags"])
    assert calls["1-4000-A-G"]["flags"] == []
    assert got.result["roh"]["P"]["assessed"] is False  # 4 informative sites: not assessed, not "no ROH"


def test_cp1_15_review_roh_splits_at_gaps_and_unsorted_input_is_sorted(tmp_path):
    pos = list(range(1_000_000, 4_000_000, 50_000)) + list(range(10_000_000, 13_000_000, 50_000))
    segs = Q.find_roh(pos, [0] * len(pos))
    assert len(segs) == 2  # the 6 Mb gap splits the run
    head = ["##fileformat=VCFv4.2", "##contig=<ID=chr2,length=242193529>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP"]
    import random as _r

    rows = [f"chr2\t{p_}\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:GQ\t1/1:0,30:30:99" for p_ in range(1_000_000, 31_000_000, 20_000)]
    rows += [f"chr2\t{p_}\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99" for p_ in range(40_000_000, 60_000_000, 20_000)]
    _r.Random(1).shuffle(rows)
    p = tmp_path / "shuffled.vcf"
    p.write_text("\n".join(head + rows) + "\n")
    roh = Q.run(str(p), proband="P").result["roh"]["P"]
    assert roh["n_segments"] == 1 and roh["segments"][0]["start"] == 1_000_000


def test_cp1_15_review_po_fs_cut_scales_with_the_unrelated_pairs(qc_result):
    cut = qc_result.result["po_fs_cut"]
    assert cut["basis"].startswith("0.1 x the median IBS0/het of 1 unrelated pair")
    pairs = {frozenset((p["a"], p["b"])): p for p in qc_result.result["relatedness"]}
    assert pairs[frozenset(("M", "F"))]["ibs0_per_het"] == pytest.approx(cut["ibs0_per_het"] * 10, rel=0.01)
    kind, problem = Q._pair_verdict("P", "M", 0.12, 0.0, "parent-offspring", None)
    assert "2nd degree relatives" in problem
    kind, problem = Q._pair_verdict("A", "B", 0.25, 0.015, "full siblings", None, po_cut=0.01)
    assert kind == "1st degree (parent-offspring vs full siblings unclear)" and problem is None


def test_cp1_15_review_few_chromosomes_make_a_verdict_inconclusive(family, tmp_path):
    vcf = _subset(family["vcf"], tmp_path / "chr12.vcf", lambda l: l.startswith(("chr1\t", "chr2\t")))
    ped = _ped(tmp_path / "swap.ped", ["fam1 M 0 0 2 1", "fam1 U 0 0 1 1", "fam1 P U M 2 2"])
    got = Q.run(vcf, ped=ped)
    up = next(p for p in got.result["relatedness"] if {p["a"], p["b"]} == {"U", "P"})
    assert up["autosomes"] == 2 and up["problem"].startswith("inconclusive — kinship from only 2 autosome(s)")


def test_cp1_15_review_sex_without_proband_and_half_sibs(family, tmp_path):
    got = Q.run(family["vcf"], sex="male", ped=family["ped"])
    assert any("--sex states the proband's sex, but no --proband was given" in w for w in got.warnings)
    ped = Q.read_ped(_ped(tmp_path / "h.ped", ["f M 0 0 2 1", "f F 0 0 1 1", "f G 0 0 1 1", "f A F M 1 1",
                                               "f B G M 1 1"]))
    assert ped.expected("A", "B") == "half siblings"
