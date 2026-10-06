"""zebra.aso: the splice-switching antisense feasibility screen.

Most sequence here is hand-built -- filler 'T' with a few bases placed on purpose --
so every window offset, dinucleotide and bystander the screen reports can be checked
by hand. Two tests use captured real sequence (tests/fixtures/aso/, built by
tools/aso/fetch_fixtures.py) so the cryptic-exon geometry is pinned to the genome and
not only to the construction.

CP1-19 is "no ASO design aid": the regression tests for it are the ones that import
`zebra.aso` and run the screen, both of which fail on a tree without this module.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from zebra import aso
from zebra.core import Outcome, UsageError

FIXTURES = Path(__file__).parent / "fixtures" / "aso"


def fixture(name):
    path = FIXTURES / f"{name}.json"
    if not path.exists():  # pragma: no cover - the fixture is committed
        pytest.skip(f"missing fixture {path}; rebuild with tools/aso/fetch_fixtures.py")
    return json.loads(path.read_text("utf-8"))


# ---------------------------------------------------------------- primitives


def test_revcomp_and_gc():
    assert aso.revcomp("ACGTN") == "NACGT"
    assert aso.revcomp("") == ""
    assert aso.gc_fraction("GCGC") == 1.0
    assert aso.gc_fraction("ATAT") == 0.0
    assert aso.gc_fraction("GCAT") == 0.5
    # N is not counted in either the numerator or the denominator
    assert aso.gc_fraction("GCNN") == 1.0
    assert aso.gc_fraction("NNNN") is None
    assert aso.gc_fraction("") is None


def test_hairpin_finds_a_planted_stem():
    # GGGGG ... CCCCC with a 5 nt loop: a 5 bp stem
    planted = aso.hairpin("GGGGG" + "ATATA" + "CCCCC")
    assert planted["stem"] == 5
    assert planted["loop"] >= aso.HAIRPIN_MIN_LOOP
    # no self-complementarity at all
    assert aso.hairpin("AAAAAAAAAAAAAAAA")["stem"] == 0
    # a stem whose arms are closer than the minimum loop is not reported
    assert aso.hairpin("GGCC", min_loop=3)["stem"] == 0


def test_longest_run():
    assert aso.longest_run("AACCCGT") == {"base": "C", "length": 3}
    assert aso.longest_run("ACGT")["length"] == 1
    assert aso.longest_run("")["length"] == 0


def test_parse_lengths():
    assert aso.parse_lengths("18-20") == (18, 19, 20)
    assert aso.parse_lengths("20") == (20,)
    assert aso.parse_lengths("18,25,18") == (18, 25)
    assert aso.parse_lengths(None) == tuple(range(aso.DEFAULT_LENGTHS[0], aso.DEFAULT_LENGTHS[1] + 1))
    for bad in ("25-18", "5-9", "0-50", "abc", "18-200"):
        with pytest.raises(UsageError):
            aso.parse_lengths(bad)


# ---------------------------------------------------------------- Region


def plant(bases, length=200, filler="T"):
    seq = [filler] * length
    for index, base in bases.items():
        seq[int(index)] = base
    return "".join(seq)


def test_region_plus_strand_round_trip():
    region = aso.Region("7", 1001, 1100, 1, plant({}, 100))
    assert len(region) == 100
    assert region.index(1001) == 0
    assert region.genomic(0) == 1001
    assert region.index(1100) == 99
    assert region.premrna == region.forward
    assert region.span(10, 20) == (1011, 1021)
    for pos in (1001, 1050, 1100):
        assert region.genomic(region.index(pos)) == pos


def test_region_minus_strand_is_the_reverse_complement():
    forward = "ACGT" + "T" * 96
    region = aso.Region("7", 1001, 1100, -1, forward)
    assert region.premrna == aso.revcomp(forward)
    # the first pre-mRNA base is the LAST genomic base on a minus-strand gene
    assert region.index(1100) == 0
    assert region.genomic(0) == 1100
    assert region.index(1001) == 99
    # a pre-mRNA slice maps back to a forward-strand span
    assert region.span(0, 9) == (1091, 1100)
    for pos in (1001, 1050, 1100):
        assert region.genomic(region.index(pos)) == pos


def test_region_slice_is_premrna_not_forward():
    forward = "AAAA" + "T" * 16
    region = aso.Region("1", 1, 20, -1, forward)
    # pre-mRNA of a minus-strand gene starts at the far end: 16 A's complement of the T's
    assert region.slice(0, 3) == "AAAA"
    assert region.slice(16, 19) == "TTTT"


# ---------------------------------------------------------------- site windows


def site(kind, index, position=None, delta=0.5):
    return {"kind": kind, "premrna_index": index, "position": position if position is not None else 1000 + index,
            "delta": delta, "spliceai_score": "DS_AG" if kind == "acceptor" else "DS_DG",
            "passed_min_delta": True, "created_by_variant": False}


def test_every_acceptor_window_covers_the_AG_and_the_first_exonic_base():
    region = aso.Region("1", 1, 400, 1, plant({}, 400))
    s = site("acceptor", 200)
    windows = aso._site_windows(region, s, (18, 25), stride=1)
    assert windows
    for w in windows:
        assert w["i0"] <= 200 - 2, w          # the AG is inside
        assert w["i1"] >= 200, w              # the first exonic base is inside
        assert w["i1"] - w["i0"] + 1 == w["length"]


def test_every_donor_window_covers_the_last_exonic_base_and_the_GT():
    region = aso.Region("1", 1, 400, 1, plant({}, 400))
    s = site("donor", 200)
    windows = aso._site_windows(region, s, (18, 25), stride=1)
    assert windows
    for w in windows:
        assert w["i0"] <= 200, w
        assert w["i1"] >= 200 + 2, w


def test_site_windows_honour_the_stride_and_the_region_edges():
    region = aso.Region("1", 1, 400, 1, plant({}, 400))
    dense = aso._site_windows(region, site("donor", 200), (20,), stride=1)
    sparse = aso._site_windows(region, site("donor", 200), (20,), stride=4)
    assert len(sparse) < len(dense)
    # a site at the very start of the region yields no acceptor window (no room for the AG)
    assert aso._site_windows(region, site("acceptor", 1), (20,), stride=1) == []


def test_body_windows_stay_inside_the_cryptic_exon():
    region = aso.Region("1", 1, 400, 1, plant({}, 400))
    pseudoexon = {"acceptor": site("acceptor", 100), "donor": site("donor", 183), "size_bp": 84}
    windows = aso._body_windows(region, pseudoexon, (18, 25), stride=10)
    assert windows
    for w in windows:
        assert 100 <= w["i0"] and w["i1"] <= 183, w
        assert w["region"] == "pseudoexon_body"
        assert w["site_premrna_index"] is None


# ---------------------------------------------------------------- candidates


def test_candidate_antisense_is_the_reverse_complement_of_the_premrna_target():
    forward = plant({40: "A", 41: "C", 42: "G"}, 100)
    region = aso.Region("5", 1, 100, 1, forward)
    window = {"region": "donor_site", "site_kind": "donor", "site_position": 41, "i0": 35, "i1": 54,
              "length": 20, "site_premrna_index": 40}
    c = aso._candidate(region, window, variant_pos=41)
    assert c["target_premrna_5to3"] == forward[35:55]
    assert c["antisense_5to3"] == aso.revcomp(c["target_premrna_5to3"])
    assert c["target_genomic"] == {"chrom": "5", "start": 36, "end": 55, "strand": "+"}
    assert c["contains_patient_variant"] is True
    assert c["variant_offset_in_target"] == 41 - 1 - 35 + 1
    # the site offsets are mirror images across the two strands
    assert c["site_offset_in_target"] == 6
    assert c["site_offset_in_antisense"] == 20 - 5


def test_candidate_on_a_minus_strand_gene_targets_the_premrna_not_the_forward_strand():
    forward = "AAAA" + "T" * 46
    region = aso.Region("X", 1, 50, -1, forward)
    window = {"region": "donor_site", "site_kind": "donor", "site_position": None, "i0": 0, "i1": 19,
              "length": 20, "site_premrna_index": None}
    c = aso._candidate(region, window, variant_pos=1_000_000)
    # pre-mRNA index 0..19 is the forward window's LAST 20 bases, reverse complemented
    assert c["target_premrna_5to3"] == aso.revcomp(forward[-20:])
    assert c["antisense_5to3"] == forward[-20:]
    assert c["target_genomic"]["start"] == 31 and c["target_genomic"]["end"] == 50
    assert c["contains_patient_variant"] is False
    assert "variant_offset_in_target" not in c


def test_candidate_flags_hairpin_and_homopolymer():
    forward = plant({}, 60, filler="A")  # 60 A's -> the antisense is 60 T's: a long run, no hairpin
    region = aso.Region("1", 1, 60, 1, forward)
    window = {"region": "pseudoexon_body", "site_kind": None, "site_position": None, "i0": 0, "i1": 19,
              "length": 20, "site_premrna_index": None}
    c = aso._candidate(region, window, variant_pos=0)
    assert c["homopolymer_flag"] is True and c["longest_homopolymer"] == 20
    assert c["hairpin_flag"] is False
    assert c["gc_fraction"] == 0.0


# ---------------------------------------------------------------- ranking


def _rank_for(**overrides):
    c = {"length": 20, "gc_fraction": 0.5, "hairpin_flag": False, "self_complementary_stem": 2,
         "homopolymer_flag": False, "longest_homopolymer": 2, "homopolymer_base": "A",
         "contains_patient_variant": False, "site_offset_in_target": 10,
         "uniqueness": {"checked": False, "reason": "not checked"}}
    c.update(overrides)
    return aso._score(c)


def test_score_is_exactly_the_sum_of_the_shown_contributions():
    rank = _rank_for()
    assert rank["score"] == pytest.approx(sum(r["contribution"] for r in rank["components"]))
    assert rank["score_max"] == pytest.approx(sum(r["weight"] for r in rank["components"] if r["counted"]))
    assert rank["fraction"] == pytest.approx(rank["score"] / rank["score_max"], abs=5e-4)
    # every component carries its weight and an explanation; nothing is hidden
    assert {r["name"] for r in rank["components"]} == set(aso.RANK_COMPONENTS)
    assert "note" not in rank          # the constant explanation lives once in the result
    for r in rank["components"]:
        assert r["basis"]
        # the explanation is stated once in the result, not duplicated on every candidate
        assert "what_it_means" not in r


def test_unchecked_uniqueness_is_left_out_of_the_score_and_of_the_maximum():
    unchecked = _rank_for()
    unique = _rank_for(uniqueness={"checked": True, "unique": True, "locus_count": 1})
    repeated = _rank_for(uniqueness={"checked": True, "unique": False, "locus_count": 7})
    assert "uniqueness" in unchecked["components_not_counted"]
    assert unchecked["score_max"] == pytest.approx(unique["score_max"] - 2.0)
    assert unique["score"] == pytest.approx(repeated["score"] + 2.0)
    assert unique["score_max"] == repeated["score_max"]


def test_variant_in_target_is_reported_with_weight_zero():
    with_variant = _rank_for(contains_patient_variant=True)
    without = _rank_for(contains_patient_variant=False)
    assert with_variant["score"] == pytest.approx(without["score"])
    row = next(r for r in with_variant["components"] if r["name"] == "variant_in_target")
    assert row["weight"] == 0.0 and row["value"] == 1.0 and row["contribution"] == 0.0


def test_site_centering_rewards_a_centred_site_and_a_body_window_is_not_scored_on_it():
    centred = _rank_for(site_offset_in_target=10)   # length 20 -> centre 10.5
    edge = _rank_for(site_offset_in_target=1)
    body = _rank_for(site_offset_in_target=None)
    def value(rank):
        return next(r for r in rank["components"] if r["name"] == "site_centering")
    assert value(centred)["value"] > value(edge)["value"]
    assert value(edge)["value"] == 0.0
    assert value(body)["value"] is None and value(body)["counted"] is False
    assert body["score_max"] == pytest.approx(centred["score_max"] - 2.0)


def test_gc_band_is_flat_inside_the_band_and_falls_outside_it():
    def value(gc):
        return next(r for r in _rank_for(gc_fraction=gc)["components"] if r["name"] == "gc_in_band")["value"]
    assert value(0.40) == value(0.50) == value(0.60) == 1.0
    assert 0.0 < value(0.30) < 1.0
    assert value(0.20) == 0.0 and value(0.80) == 0.0
    assert value(0.10) == 0.0 and value(0.95) == 0.0


# ---------------------------------------------------------------- event detection


def transcript(deltas, strand="+", refseq="NM_000001.1"):
    return {"gene": "GENE", "transcript": "ENST00000000001", "refseq": refseq, "strand": strand,
            "priority_label": "MANE Select", "biotype": "protein_coding", "deltas": deltas}


def delta(score, value, position, label="gain"):
    return {"score": score, "label": label, "delta": value, "offset": 0, "position": position,
            "ref_prob": 0.0, "alt_prob": value}


def regions_for(forward_reference, variant_index, alt, start=1, strand=1, applied=True):
    """The (patient, reference) pair `fetch_region` would return for a 1-base substitution."""
    patient_seq = forward_reference[:variant_index] + alt + forward_reference[variant_index + 1:]
    end = start + len(forward_reference) - 1
    patient = aso.Region("7", start, end, strand, patient_seq)
    patient.allele_applied = ("yes: %s>%s at 7:%d" % (forward_reference[variant_index], alt, start + variant_index)
                              if applied else "no: not a single-base substitution")
    return patient, aso.Region("7", start, end, strand, forward_reference)


def test_detect_event_calls_a_cryptic_exon_from_a_gained_acceptor_and_donor():
    # acceptor at index 100 (AG at 98-99), donor at index 183 (GT at 184-185): an 84 nt exon
    reference = plant({98: "A", 99: "G", 184: "G", 185: "C"}, 300)
    region, reference_region = regions_for(reference, 185, "T")
    out = aso.detect_event(region, reference_region, transcript([
        delta("DS_AG", 0.35, 101), delta("DS_DG", 0.40, 184),
    ]), "auto", 0.2)
    assert out["type"] == "pseudoexon"
    assert out["pseudoexon"]["size_bp"] == 84
    assert out["pseudoexon"]["in_frame"] is True
    assert out["pseudoexon"]["genomic"] == [101, 184]
    donor = out["pseudoexon"]["donor"]
    assert donor["dinucleotide_patient_allele"] == "GT"
    assert donor["dinucleotide_reference"] == "GC"
    assert donor["created_by_variant"] is True
    acceptor = out["pseudoexon"]["acceptor"]
    assert acceptor["dinucleotide_patient_allele"] == "AG"
    assert acceptor["created_by_variant"] is False
    assert "canonical in the patient's own sequence" in out["sequence_check"]


def test_detect_event_stops_when_no_gain_passes_the_floor():
    reference = plant({}, 300)
    region, reference_region = regions_for(reference, 150, "A")
    out = aso.detect_event(region, reference_region,
                           transcript([delta("DS_AG", 0.03, 151), delta("DS_DG", 0.0, 160)]), "auto", 0.2)
    assert out["type"] is None and out["stop"] is True
    assert "no gained splice site" in out["reason"]
    assert "nothing for a splice-switching antisense oligonucleotide to aim at" in out["reason"]


def test_a_gain_in_the_uninformative_band_is_screened_but_flagged():
    reference = plant({98: "A", 99: "G"}, 300)
    region, reference_region = regions_for(reference, 150, "A")
    out = aso.detect_event(region, reference_region, transcript([delta("DS_AG", 0.162, 101)]), "auto", 0.2)
    assert out["type"] == "cryptic_acceptor"
    assert out["sites"][0]["passed_min_delta"] is False
    assert "BELOW the ClinGen PP3 threshold" in out["confidence"]
    assert out["sites"][0]["delta_band"].startswith("uninformative")


def test_a_pseudoexon_of_an_implausible_size_is_reported_as_two_separate_sites():
    reference = plant({98: "A", 99: "G", 110: "G", 111: "T"}, 300)   # only 10 bp apart
    region, reference_region = regions_for(reference, 150, "A")
    out = aso.detect_event(region, reference_region, transcript([
        delta("DS_AG", 0.4, 101), delta("DS_DG", 0.5, 110)]), "auto", 0.2)
    assert out["type"] == "cryptic_donor"          # the stronger single gain
    assert "pseudoexon_rejected" in out
    assert "outside the" in out["pseudoexon_rejected"]


def test_forcing_pseudoexon_stops_when_only_one_gain_exists():
    reference = plant({98: "A", 99: "G"}, 300)
    region, reference_region = regions_for(reference, 150, "A")
    out = aso.detect_event(region, reference_region, transcript([delta("DS_AG", 0.4, 101)]), "pseudoexon", 0.2)
    assert out["type"] is None and out["stop"] is True
    assert "cryptic exon was asked for" in out["reason"]


def test_forcing_a_cryptic_donor_stops_when_only_an_acceptor_gain_exists():
    reference = plant({98: "A", 99: "G"}, 300)
    region, reference_region = regions_for(reference, 150, "A")
    out = aso.detect_event(region, reference_region, transcript([delta("DS_AG", 0.4, 101)]),
                           "cryptic_donor", 0.2)
    assert out["type"] is None and out["stop"] is True
    assert "no donor gain" in out["reason"]


def test_a_non_canonical_predicted_site_is_called_out_not_hidden():
    reference = plant({}, 300)  # no AG anywhere: the "acceptor" is not a splice site
    region, reference_region = regions_for(reference, 150, "A")
    out = aso.detect_event(region, reference_region, transcript([delta("DS_AG", 0.5, 101)]), "auto", 0.2)
    assert out["type"] == "cryptic_acceptor"
    assert "NOT canonical" in out["sequence_check"]
    assert "unverified" in out["sequence_check"]


def test_created_by_variant_separates_a_new_site_from_an_upgraded_one():
    """GC->GT is how CFTR c.3718-2477C>T makes its donor: an upgrade, not a site from nothing."""
    # the reference has no splice-site dinucleotide here at all
    reference = plant({184: "G", 185: "A"}, 300)
    region, reference_region = regions_for(reference, 185, "T")
    made = aso.detect_event(region, reference_region, transcript([delta("DS_DG", 0.4, 184)]), "auto", 0.2)["sites"][0]
    assert made["created_by_variant"] is True
    assert "not a splice-site dinucleotide" in made["created_by_variant_basis"]

    # the reference has the minor-class GC; the variant makes it the major-class GT
    reference = plant({184: "G", 185: "C"}, 300)
    region, reference_region = regions_for(reference, 185, "T")
    upgraded = aso.detect_event(region, reference_region, transcript([delta("DS_DG", 0.4, 184)]),
                                "auto", 0.2)["sites"][0]
    assert upgraded["created_by_variant"] is True
    assert "minor-class" in upgraded["created_by_variant_basis"]

    # the dinucleotide is unchanged: the variant may cause the site to be USED, not created
    reference = plant({98: "A", 99: "G"}, 300)
    region, reference_region = regions_for(reference, 150, "A")
    unchanged = aso.detect_event(region, reference_region, transcript([delta("DS_AG", 0.4, 101)]),
                                 "auto", 0.2)["sites"][0]
    assert unchanged["created_by_variant"] is False
    assert unchanged["dinucleotide_changed_by_variant"] is False
    assert "it is a site that the variant may cause to be USED" in unchanged["created_by_variant_basis"]


def test_a_minor_class_GC_donor_is_canonical_and_labelled():
    reference = plant({184: "G", 185: "T"}, 300)
    region, reference_region = regions_for(reference, 185, "C")   # GT -> GC on the patient allele
    out = aso.detect_event(region, reference_region, transcript([delta("DS_DG", 0.4, 184)]), "auto", 0.2)
    donor = out["sites"][0]
    assert donor["dinucleotide_patient_allele"] == "GC"
    assert donor["canonical_on_patient_allele"] is True
    assert donor["minor_class_donor"] is True


# ---------------------------------------------------------------- real sequence


def test_cftr_cryptic_exon_geometry_against_captured_sequence():
    """CP1-19: the 84 nt CFTR cryptic exon, from real sequence and real SpliceAI positions."""
    cap = fixture("cftr_intron22_3849plus10kb")
    start = cap["start"]
    variant_pos = 117639961
    region, reference_region = regions_for(cap["sequence"], variant_pos - start, "T", start=start)
    out = aso.detect_event(region, reference_region, transcript([
        delta("DS_AG", 0.104, 117639876), delta("DS_DG", 0.162, 117639959),
    ], refseq="NM_000492.4"), "auto", 0.2)
    assert out["type"] == "pseudoexon"
    assert out["pseudoexon"]["size_bp"] == 84          # the size reported for this variant
    assert out["pseudoexon"]["genomic"] == [117639876, 117639959]
    donor, acceptor = out["pseudoexon"]["donor"], out["pseudoexon"]["acceptor"]
    assert acceptor["dinucleotide_patient_allele"] == "AG"
    assert donor["dinucleotide_reference"] == "GC" and donor["dinucleotide_patient_allele"] == "GT"
    assert donor["created_by_variant"] is True and acceptor["created_by_variant"] is False
    assert "BELOW the ClinGen PP3 threshold" in out["confidence"]
    # the group with the variant-created site comes first
    assert aso._group_order(out)[0] == "donor_site"


def test_minus_strand_exon_boundaries_from_captured_sequence():
    """On a minus-strand gene the acceptor is the higher coordinate and the donor the lower."""
    cap = fixture("tp53_exon4_boundaries")
    start = cap["start"]
    region = aso.Region("17", start, cap["end"], -1, cap["sequence"])
    acceptor_pos, donor_pos = 7676272, 7675994       # TP53 exon 4 on the minus strand
    assert aso._dinucleotide(region, region.index(acceptor_pos), "acceptor") == "AG"
    assert aso._dinucleotide(region, region.index(donor_pos), "donor") == "GT"
    # read the other way round, neither is canonical: the strand handling is load-bearing
    plus = aso.Region("17", start, cap["end"], 1, cap["sequence"])
    assert aso._dinucleotide(plus, plus.index(acceptor_pos), "acceptor") != "AG"


# ---------------------------------------------------------------- the whole screen


def stub_predict(deltas, strand="+", status="ran", reason=None, chrom="7", pos=1000, ref="C", alt="T",
                 pangolin=None):
    tx = transcript(deltas, strand=strand)
    row = {"model": "spliceai", "status": status, "transcript": tx, "distance": 500, "mask": 0,
           "url": "https://example.invalid/spliceai", "reading": "..."}
    if reason:
        row["reason"] = reason
    models = [row]
    if pangolin:
        models.append({"model": "pangolin", "status": "ran", "headline": pangolin,
                       "url": "https://example.invalid/pangolin"})
    variant = {"input": "test", "kind": "vcf", "assembly": "GRCh38", "chrom": chrom, "pos": pos,
               "ref": ref, "alt": alt, "resolved_from": "given as chrom-pos-ref-alt", "ids": [], "hgvs_g": None}
    return Outcome({"variant": variant, "assembly": "GRCh38", "models": models}, sources=[], warnings=[])


def wire(monkeypatch, predicted, forward, start):
    from zebra import s2f
    from zebra.sources import ensembl

    monkeypatch.setattr(s2f, "predict", lambda *a, **k: predicted)

    def sequence(chrom, lo, hi, assembly="GRCh38", strand=1):
        assert lo >= start and hi <= start + len(forward) - 1, (lo, hi, start, len(forward))
        return Outcome(forward[lo - start:hi - start + 1],
                       sources=[{"db": "Ensembl sequence", "record": f"{chrom}:{lo}-{hi}", "url": "x"}])

    monkeypatch.setattr(ensembl, "sequence", sequence)
    monkeypatch.setattr(aso.ensembl, "sequence", sequence)


def test_screen_end_to_end_on_a_planted_cryptic_exon(monkeypatch):
    """CP1-19: `zebra aso` turns a predicted cryptic exon into target windows with provenance."""
    start = 1
    reference = plant({98: "A", 99: "G", 184: "G", 185: "C"}, 400)
    predicted = stub_predict([delta("DS_AG", 0.35, 101), delta("DS_DG", 0.40, 184)], pos=186, ref="C", alt="T",
                             pangolin={"score": "DS_SG", "value": 0.5, "position": 184})
    wire(monkeypatch, predicted, reference, start)
    out = aso.screen("7-186-C-T", lengths=(20,), top=5, with_precedents=False)
    r = out.result
    assert r["event"]["type"] == "pseudoexon"
    assert r["candidate_count"] > 0 and r["candidates_shown"] == 5
    assert {g["region"] for g in r["groups"]} == {"acceptor_site", "donor_site", "pseudoexon_body"}
    assert r["uniqueness_checked"] is False
    for c in r["candidates"]:
        assert c["antisense_5to3"] == aso.revcomp(c["target_premrna_5to3"])
        assert c["uniqueness"]["checked"] is False
        assert "not checked" in c["uniqueness"]["reason"]
        assert "uniqueness" in c["rank"]["components_not_counted"]
    assert "NOT checked" in r["uniqueness_note"]
    assert r["corroboration"]["model"] == "pangolin" and "not independent" in r["corroboration"]["caveat"]
    assert out.text and "cryptic exon" in out.text
    # the screen never claims more than a feasibility screen
    assert "Not a design" in r["what_this_is"]
    assert any("RNA evidence" in p for p in r["preconditions"])
    assert any("no chemistry" in p.lower() or "no dose" in p.lower() for p in r["limits"])


def test_screen_says_so_and_stops_when_no_event_is_predicted(monkeypatch):
    """CP1-19: a variant with no predicted aberrant event gets a refusal, not a target list."""
    reference = plant({}, 400)
    predicted = stub_predict([delta("DS_AG", 0.02, 101), delta("DS_DG", 0.01, 184),
                              {"score": "DS_DL", "delta": 0.3, "position": 150, "label": "donor loss"}],
                             pos=186, ref="T", alt="A")
    wire(monkeypatch, predicted, reference, 1)
    out = aso.screen("7-186-C-T", lengths=(20,), with_precedents=False)
    r = out.result
    assert r["stop"] is True
    assert r["candidates"] == [] and r["candidate_count"] == 0
    assert r["event"]["type"] is None
    assert "no gained splice site" in r["event"]["reason"]
    assert "no antisense target" in out.text
    assert "cannot restore one the variant weakened" in out.text


def test_screen_stops_when_spliceai_did_not_run(monkeypatch):
    reference = plant({}, 400)
    predicted = stub_predict([], status="error", reason="the service refused the request", pos=186,
                             ref="T", alt="A")
    wire(monkeypatch, predicted, reference, 1)
    out = aso.screen("7-186-C-T", with_precedents=False)
    r = out.result
    assert r["stop"] is True and r["event"]["type"] is None
    assert "the service refused the request" in r["event"]["reason"]
    # 'not measured' must never read as 'no aberrant splicing'
    assert "not measured" in r["event"]["reason"]


def test_screen_refuses_a_ref_mismatch(monkeypatch):
    reference = plant({}, 400)  # filler T everywhere, so REF=C at pos 186 is wrong
    predicted = stub_predict([delta("DS_DG", 0.4, 184)], pos=186, ref="C", alt="T")
    wire(monkeypatch, predicted, reference, 1)
    with pytest.raises(UsageError) as err:
        aso.screen("7-186-C-T", with_precedents=False)
    assert "REF mismatch" in str(err.value)


def test_screen_rejects_bad_arguments():
    with pytest.raises(UsageError):
        aso.screen("7-186-C-T", assembly="hg38")
    with pytest.raises(UsageError):
        aso.screen("7-186-C-T", event="nonsense")
    with pytest.raises(UsageError):
        aso.screen("7-186-C-T", lengths=(3,))


def test_screen_warns_and_drops_windows_with_unresolved_bases(monkeypatch):
    reference = plant({98: "A", 99: "G"}, 400, filler="N")
    predicted = stub_predict([delta("DS_AG", 0.4, 101)], pos=186, ref="N", alt="N")
    wire(monkeypatch, predicted, reference, 1)
    out = aso.screen("7-186-N-N", lengths=(20,), with_precedents=False)
    assert out.result["candidates"] == []
    assert any("does not resolve" in w for w in out.warnings)


def test_uniqueness_result_is_carried_into_the_rank(monkeypatch):
    reference = plant({98: "A", 99: "G", 184: "G", 185: "C"}, 400)
    predicted = stub_predict([delta("DS_DG", 0.40, 184)], pos=186, ref="C", alt="T")
    wire(monkeypatch, predicted, reference, 1)
    from zebra.sources import ncbi_blast

    def locate(queries, max_wait=None):
        out = {}
        for i, name in enumerate(queries):
            out[name] = {"checked": True, "unique": i == 0, "locus_count": 1 if i == 0 else 4,
                         "loci": [], "non_primary_alignments": 0, "rid": "TESTRID"}
        return Outcome(out, sources=[{"db": "NCBI BLAST", "record": "TESTRID", "url": "x"}])

    monkeypatch.setattr(ncbi_blast, "locate", locate)
    out = aso.screen("7-186-C-T", lengths=(20,), top=3, uniqueness=True, with_precedents=False)
    r = out.result
    assert r["uniqueness_checked"] is True
    assert all(c["uniqueness"]["checked"] for c in r["candidates"])
    assert [c["uniqueness"]["unique"] for c in r["candidates"]][0] is True   # unique sorts first
    for c in r["candidates"]:
        assert "uniqueness" not in c["rank"]["components_not_counted"]
    assert any(s.get("db") == "NCBI BLAST" for s in out.sources)


def test_a_failed_uniqueness_check_never_reads_as_unique(monkeypatch):
    reference = plant({98: "A", 99: "G", 184: "G", 185: "C"}, 400)
    predicted = stub_predict([delta("DS_DG", 0.40, 184)], pos=186, ref="C", alt="T")
    wire(monkeypatch, predicted, reference, 1)
    from zebra.http import SourceError
    from zebra.sources import ncbi_blast

    def boom(queries, max_wait=None):
        raise SourceError("NCBI BLAST", "https://example.invalid", 503, "service unavailable")

    monkeypatch.setattr(ncbi_blast, "locate", boom)
    out = aso.screen("7-186-C-T", lengths=(20,), top=2, uniqueness=True, with_precedents=False)
    for c in out.result["candidates"]:
        assert c["uniqueness"]["checked"] is False
        assert c["uniqueness"].get("unique") is None
    assert any("NCBI BLAST" in w for w in out.warnings)


def test_exon_skip_uses_the_authentic_sites_of_the_exon(monkeypatch):
    reference = plant({98: "A", 99: "G", 184: "G", 185: "T"}, 400)
    predicted = stub_predict([delta("DS_AG", 0.0, 101)], pos=150, ref="T", alt="A")
    wire(monkeypatch, predicted, reference, 1)
    monkeypatch.setattr(aso, "transcript_exons", lambda tid, assembly="GRCh38": Outcome(
        {"transcript": tid, "strand": 1, "chrom": "7",
         "exons": [{"id": "E1", "start": 101, "end": 184, "number": 1, "length": 84}], "exon_count": 1},
        sources=[{"db": "Ensembl lookup", "record": tid, "url": "x"}]))
    out = aso.screen("7-150-T-A", lengths=(20,), event="exon_skip", top=4, with_precedents=False)
    r = out.result
    assert r["event"]["type"] == "exon_skip"
    assert r["event"]["exon"]["number"] == 1
    assert r["event"]["frame"]["in_frame"] is True
    kinds = {s["kind"] for s in r["event"]["sites"]}
    assert kinds == {"acceptor", "donor"}
    assert all(s["authentic_site"] for s in r["event"]["sites"])
    assert all(s["created_by_variant"] is False for s in r["event"]["sites"])
    assert r["candidate_count"] > 0
    assert "skip exon 1" in out.text


def test_exon_skip_reports_a_frameshifting_exon_as_such(monkeypatch):
    reference = plant({98: "A", 99: "G", 183: "G", 184: "T"}, 400)
    predicted = stub_predict([delta("DS_AG", 0.0, 101)], pos=150, ref="T", alt="A")
    wire(monkeypatch, predicted, reference, 1)
    monkeypatch.setattr(aso, "transcript_exons", lambda tid, assembly="GRCh38": Outcome(
        {"transcript": tid, "strand": 1, "chrom": "7",
         "exons": [{"id": "E1", "start": 101, "end": 183, "number": 1, "length": 83}], "exon_count": 1},
        sources=[]))
    out = aso.screen("7-150-T-A", lengths=(20,), event="exon_skip", with_precedents=False)
    assert out.result["event"]["frame"]["in_frame"] is False
    assert "SHIFTS the frame" in out.result["event"]["frame"]["note"]


def test_exon_skip_stops_when_the_transcript_structure_is_unavailable(monkeypatch):
    reference = plant({}, 400)
    predicted = stub_predict([delta("DS_AG", 0.0, 101)], pos=150, ref="T", alt="A")
    wire(monkeypatch, predicted, reference, 1)
    from zebra.http import SourceError

    def boom(tid, assembly="GRCh38"):
        raise SourceError("Ensembl lookup", "https://example.invalid", 503, "down")

    monkeypatch.setattr(aso, "transcript_exons", boom)
    out = aso.screen("7-150-T-A", event="exon_skip", with_precedents=False)
    assert out.result["stop"] is True
    assert "exon structure could not be retrieved" in out.result["event"]["reason"]
    assert any("Ensembl transcript exons" in w for w in out.warnings)


def test_precedents_cite_nothing_when_europe_pmc_is_unreachable(monkeypatch):
    from zebra.http import SourceError
    from zebra.sources import europepmc

    def boom(*a, **k):
        raise SourceError("Europe PMC", "https://example.invalid", 503, "down")

    monkeypatch.setattr(europepmc, "search", boom)
    monkeypatch.setattr(aso.europepmc, "search", boom)
    out = aso.precedents()
    assert out.result["records"]
    for block in out.result["records"]:
        assert block["hits"] == []
        assert "nothing is cited" in block["note"]
    assert len(out.warnings) == len(aso.PRECEDENT_QUERIES)
    assert "never from memory" in aso.precedents.__doc__


def test_precedents_report_each_retrieved_record_with_its_id(monkeypatch):
    from zebra.sources import europepmc

    def search(query, limit=15, sort=None, result_type="core"):
        return Outcome({"query": query, "hitCount": 1, "hits": [
            {"pmid": "31597037", "title": "Patient-Customized Oligonucleotide Therapy for a Rare Genetic Disease.",
             "journal": "N Engl J Med", "year": "2019", "url": "https://europepmc.org/article/MED/31597037",
             "abstract": "x" * 900, "doi": "10.1056/NEJMoa1813279", "authors": "Kim J et al.",
             "fullTextUrl": None}]}, sources=[{"db": "Europe PMC", "record": query, "url": "u"}])

    monkeypatch.setattr(europepmc, "search", search)
    monkeypatch.setattr(aso.europepmc, "search", search)
    out = aso.precedents()
    keys = [b["key"] for b in out.result["records"]]
    assert keys == [q["key"] for q in aso.PRECEDENT_QUERIES]
    hit = out.result["records"][0]["hits"][0]
    assert hit["pmid"] == "31597037" and hit["year"] == "2019"
    assert len(hit["abstract_excerpt"]) <= aso.ABSTRACT_EXCERPT + 3
    assert "retrieved from Europe PMC" in out.result["rule"]


def test_the_command_is_registered_with_its_flags():
    """CP1-19: `zebra aso` is reachable from the CLI, with the flags the screen needs."""
    import argparse

    from zebra.commands import aso as command

    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command")
    command.register(sub)
    args = parser.parse_args(["aso", "7-186-C-T", "--lengths", "18-25", "--event", "pseudoexon",
                              "--uniqueness", "--top", "5"])
    assert args.variant == "7-186-C-T" and args.lengths == "18-25"
    assert args.event == "pseudoexon" and args.uniqueness is True and args.top == 5
    assert callable(args.func)
    assert parser.parse_args(["aso", "7-186-C-T"]).event == "auto"
    assert parser.parse_args(["aso", "7-186-C-T"]).uniqueness is False


# ---------------------------------------------------------------- live


@pytest.mark.live
def test_live_cftr_pseudoexon_screen():
    out = aso.screen("NM_000492.4:c.3718-2477C>T", lengths=(20, 25), top=4, with_precedents=True)
    r = out.result
    assert r["transcript"]["gene"] == "CFTR"
    assert r["event"]["type"] == "pseudoexon"
    assert r["event"]["pseudoexon"]["size_bp"] == 84
    assert r["event"]["pseudoexon"]["genomic"] == [117639876, 117639959]
    donor = r["event"]["pseudoexon"]["donor"]
    assert donor["created_by_variant"] is True and donor["dinucleotide_patient_allele"] == "GT"
    assert r["candidates"] and r["candidates"][0]["region"] == "donor_site"
    pmids = [h["pmid"] for b in r["precedents"]["records"] for h in b["hits"]]
    assert "31597037" in pmids   # the milasen report, retrieved, not remembered


@pytest.mark.live
def test_live_variant_with_no_predicted_event_stops():
    out = aso.screen("rs1042522", with_precedents=False)
    assert out.result["stop"] is True
    assert out.result["event"]["type"] is None
    assert "no gained splice site" in out.result["event"]["reason"]


@pytest.mark.live
def test_live_transcript_exons_cftr():
    out = aso.transcript_exons("ENST00000003084")
    assert out.result["exon_count"] == 27
    assert out.result["strand"] == 1
    assert out.result["exons"][0]["start"] == 117480025


@pytest.mark.live
def test_live_uniqueness_check():
    from zebra.sources import ncbi_blast

    # a 25-mer from CFTR intron 22 and a poly-A stretch
    out = ncbi_blast.locate({"cftr25": "CAGTATTAAAATGGCGAGTAAGACA", "polya": "A" * 25})
    for name in ("cftr25", "polya"):
        row = out.result[name]
        assert "checked" in row
        if row["checked"]:
            assert isinstance(row["locus_count"], int)
            assert row["criterion"]


# ---------------------------------------------------------------- the uniqueness backend

# The real shape of a `FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular` body, captured from
# RID C9M8N8ES014 on 2026-10-06: the service's HTML preamble, `# blastn`, `# Iteration`,
# `# RID`, `# Database`, a `# Fields:` line and `# N hits found`, then one tab-separated
# line per alignment named by the FASTA title. The first two alignments of cftr_ref25
# are verbatim; `cftr_ref25` is the 25 nt at GRCh38 7:117,639,947-117,639,971, and its
# own locus comes back at exactly those coordinates (NC_060931.1 is T2T-CHM13).
TABULAR = """\
<p><!--
QBlastInfoBegin
\tStatus=READY
QBlastInfoEnd
--><p>
<PRE>
# blastn
# Iteration: 0
# Query: cand01
# RID: C9M8N8ES014
# Database: GPIPE/9606/current/all_top_level
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 2 hits found
cand01\tNC_060931.1\t100.000\t25\t0\t0\t1\t25\t118955310\t118955334\t1.12e-04\t47.3
cand01\tNC_000007.14\t100.000\t25\t0\t0\t1\t25\t117639947\t117639971\t1.12e-04\t47.3
# blastn
# Iteration: 0
# Query: cand02
# RID: C9M8N8ES014
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 4 hits found
cand02\tNC_000007.14\t100.000\t25\t0\t0\t1\t25\t117639947\t117639971\t1.12e-04\t47.3
cand02\tNC_000002.12\t96.000\t25\t1\t0\t1\t25\t500\t476\t3e-04\t42.1
cand02\tNT_187361.1\t100.000\t25\t0\t0\t1\t25\t900\t924\t1.12e-04\t47.3
cand02\tNC_000007.14\t90.000\t20\t2\t0\t1\t20\t900000\t900019\t0.2\t30.0
# blastn
# Iteration: 0
# Query: cand03
# RID: C9M8N8ES014
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 0 hits found
</PRE>
"""


def test_blast_tabular_is_parsed_per_query_using_the_reports_own_field_order():
    from zebra.sources import ncbi_blast as nb

    blocks, order, truncated = nb.parse_tabular(TABULAR)
    assert order == ["cand01", "cand02", "cand03"]
    assert truncated is False
    assert len(blocks["cand02"]) == 4 and blocks["cand03"] == []
    assert len(blocks["cand01"]) == 2
    assert [r["subject"] for r in blocks["cand01"]] == ["NC_060931.1", "NC_000007.14"]
    row = blocks["cand01"][1]
    assert row["align_len"] == 25 and row["mismatches"] == 0 and row["s_start"] == 117639947


def test_blast_tabular_follows_a_reordered_field_line():
    from zebra.sources import ncbi_blast as nb

    swapped = ("# Query: q\n"
               "# Fields: subject acc.ver, query acc.ver, alignment length, % identity, mismatches, gap opens, "
               "s. start, s. end\n"
               "# 1 hits found\n"
               "NC_000001.11\tq\t20\t100.000\t0\t0\t10\t29\n")
    blocks, _, _ = nb.parse_tabular(swapped)
    row = blocks["q"][0]
    assert row["subject"] == "NC_000001.11" and row["align_len"] == 20 and row["identity"] == 100.0
    assert row["s_start"] == 10 and row["s_end"] == 29


def test_blast_loci_count_only_full_length_near_exact_primary_chromosome_hits():
    from zebra.sources import ncbi_blast as nb

    blocks, _, _ = nb.parse_tabular(TABULAR)
    one = nb._loci(blocks["cand01"], 25)
    # the T2T-CHM13 copy of the same locus must not make it look duplicated
    assert one["locus_count"] == 1 and one["unique"] is True
    assert one["loci"][0]["chrom"] == "7" and one["loci"][0]["strand"] == "+"
    assert one["loci"][0]["start"] == 117639947 and one["loci"][0]["end"] == 117639971
    assert one["non_primary_alignments"] == 1
    many = nb._loci(blocks["cand02"], 25)
    # chr7 exact + chr2 one-mismatch = 2 loci; the unplaced scaffold and the 20/25
    # partial alignment are not loci
    assert many["locus_count"] == 2 and many["unique"] is False
    assert {l["chrom"] for l in many["loci"]} == {"7", "2"}
    assert many["non_primary_alignments"] == 1
    assert nb._loci(blocks["cand02"], 25)["loci"][0]["mismatches"] in (0, 1)
    # the minus-strand alignment keeps its coordinates ordered low-to-high
    minus = next(l for l in many["loci"] if l["chrom"] == "2")
    assert minus["start"] < minus["end"] and minus["strand"] == "-"


def test_blast_zero_loci_is_void_not_a_count_of_zero():
    from zebra.sources import ncbi_blast as nb

    warnings = []
    row = dict(nb._loci([], 25), checked=True, rid="R", url="u")
    out = nb._void_if_self_missing(row, warnings, "cand03")
    assert out["checked"] is False
    assert "not even its own locus" in out["reason"]
    assert "unique" not in out
    assert warnings and "void" in warnings[0]


def test_blast_fasta_refuses_input_it_cannot_search_safely():
    from zebra.sources import ncbi_blast as nb

    assert nb.fasta({"a": "acgt"}) == ">a\nACGT\n"
    assert nb.fasta({"a": "AC-GT!"}) == ">a\nACGT\n"          # punctuation stripped
    with pytest.raises(ValueError):
        nb.fasta({})
    with pytest.raises(ValueError):
        nb.fasta({"a": "XXXX"})                                # no DNA bases at all
    with pytest.raises(ValueError):
        nb.fasta({"a b": "ACGT"})                              # a title that would break the FASTA
    with pytest.raises(ValueError):
        nb.fasta({"a": "A" * (nb.MAX_QUERY_LEN + 1)})
    with pytest.raises(ValueError):
        nb.fasta({"q%d" % i: "ACGT" for i in range(nb.MAX_QUERIES + 1)})


def test_blast_parsing_stops_at_the_line_cap():
    from zebra.sources import ncbi_blast as nb

    header = ("# Query: q\n"
              "# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, "
              "q. start, q. end, s. start, s. end, evalue, bit score\n")
    line = "q\tNC_000001.11\t100.000\t20\t0\t0\t1\t20\t1\t20\t1e-5\t40\n"
    blocks, _, truncated = nb.parse_tabular(header + line * (nb.MAX_REPORT_LINES + 10))
    assert truncated is True
    assert len(blocks["q"]) <= nb.MAX_REPORT_LINES


def test_a_donor_gain_upstream_of_an_acceptor_gain_is_not_an_exon():
    reference = plant({98: "A", 99: "G", 50: "G", 51: "T"}, 300)
    region, reference_region = regions_for(reference, 51, "T")
    out = aso.detect_event(region, reference_region, transcript([
        delta("DS_AG", 0.4, 101), delta("DS_DG", 0.5, 50)]), "auto", 0.2)
    assert out["type"] == "cryptic_donor"
    assert "upstream of the gained acceptor" in out["pseudoexon_rejected"]
    assert "cannot bound one exon" in out["pseudoexon_rejected"]


def test_window_enumeration_is_capped_and_says_how_many_it_cut(monkeypatch):
    reference = plant({98: "A", 99: "G", 1099: "G", 1100: "T"}, 1400)
    predicted = stub_predict([delta("DS_AG", 0.4, 101), delta("DS_DG", 0.5, 1099)],
                             pos=600, ref="T", alt="A")
    wire(monkeypatch, predicted, reference, 1)
    monkeypatch.setattr(aso, "MAX_WINDOWS", 25)
    out = aso.screen("7-600-T-A", lengths=(18, 25), top=3, body_stride=1, stride=1, with_precedents=False)
    assert out.result["windows_enumerated"] > 25
    assert out.result["candidate_count"] <= 25
    assert any("thinned evenly across every group and length" in w for w in out.warnings)
    # the cut must not fall entirely on one group or one length
    assert {g["region"] for g in out.result["groups"]} >= {"acceptor_site", "donor_site", "pseudoexon_body"}


def test_blast_rejects_a_body_that_is_not_a_report():
    """The service answers a report it cannot serve with a short message and HTTP 200."""
    from zebra.sources import ncbi_blast as nb

    assert nb._is_tabular(TABULAR) is None
    for body in ("", "   ",
                 "SYSTEM CAN'T PROCESS YOUR REQUEST, PLEASE CONTACT blasthelp.RID: C9M8N8ES014<BR>",
                 "<p><!--\nQBlastInfoBegin\n\tStatus=READY\nQBlastInfoEnd\n--></p>",
                 "<html><body>maintenance</body></html>"):
        reason = nb._is_tabular(body)
        assert reason and "not a BLAST tabular report" in reason, body[:40]


def test_blast_locate_end_to_end_offline(monkeypatch):
    from zebra.sources import ncbi_blast as nb

    monkeypatch.setattr(nb, "FIRST_POLL_AFTER", 0.0)
    monkeypatch.setattr(nb, "POLL_INTERVAL", 0.0)
    monkeypatch.setattr(nb, "submit", lambda q, refresh=False: ("TESTRID", {"db": "NCBI BLAST", "record": "TESTRID",
                                                                            "url": "u"}))
    monkeypatch.setattr(nb, "status", lambda rid: "READY")
    monkeypatch.setattr(nb, "report", lambda rid: (TABULAR, {"db": "NCBI BLAST", "record": rid, "url": "u"}))
    out = nb.locate({"cand01": "CAGTATTAAAATGGCGAGTAAGACA",
                     "cand02": "CAGTATTAAAATGGTGAGTAAGACA",
                     "cand03": "A" * 25})
    assert out.result["cand01"]["checked"] is True and out.result["cand01"]["unique"] is True
    assert out.result["cand02"]["checked"] is True and out.result["cand02"]["locus_count"] == 2
    # the query with no alignment at all is void, not "zero loci"
    assert out.result["cand03"]["checked"] is False
    assert "not even its own locus" in out.result["cand03"]["reason"]
    assert any("void" in w for w in out.warnings)
    assert len(out.sources) == 2


def test_blast_locate_degrades_honestly_when_the_search_does_not_finish(monkeypatch):
    from zebra.sources import ncbi_blast as nb

    monkeypatch.setattr(nb, "FIRST_POLL_AFTER", 0.0)
    monkeypatch.setattr(nb, "POLL_INTERVAL", 0.0)
    monkeypatch.setattr(nb, "submit", lambda q, refresh=False: ("TESTRID", {"db": "NCBI BLAST", "record": "TESTRID"}))
    monkeypatch.setattr(nb, "report", lambda rid: pytest.fail("the report must not be fetched before READY"))
    for state, marker in (("WAITING", "time budget"), ("FAILED", "FAILED"), ("UNKNOWN", "even after re-submitting")):
        monkeypatch.setattr(nb, "status", lambda rid, s=state: s)
        out = nb.locate({"cand01": "ACGT" * 6}, max_wait=0.0 if state == "WAITING" else 5.0)
        row = out.result["cand01"]
        assert row["checked"] is False and row.get("unique") is None
        assert marker in row["reason"], (state, row["reason"])
        assert any("not checked" in w for w in out.warnings)


def test_an_expired_request_id_is_resubmitted_once_not_replayed_from_the_cache(monkeypatch):
    """P1: the submit is cached for hours, so "rerun it" would hand back the same dead id."""
    from zebra.sources import ncbi_blast as nb

    monkeypatch.setattr(nb, "FIRST_POLL_AFTER", 0.0)
    monkeypatch.setattr(nb, "POLL_INTERVAL", 0.0)
    submits = []

    def submit(queries, refresh=False):
        submits.append(bool(refresh))
        return ("RIDFRESH" if refresh else "RIDSTALE"), {"db": "NCBI BLAST", "record": "x"}

    states = {"RIDSTALE": "UNKNOWN", "RIDFRESH": "READY"}
    monkeypatch.setattr(nb, "submit", submit)
    monkeypatch.setattr(nb, "status", lambda rid: states[rid])
    monkeypatch.setattr(nb, "report", lambda rid: (TABULAR, {"db": "NCBI BLAST", "record": rid}))
    out = nb.locate({"cand01": "CAGTATTAAAATGGCGAGTAAGACA"}, max_wait=30.0)
    assert submits == [False, True], submits
    assert out.result["cand01"]["checked"] is True
    assert out.result["cand01"]["rid"] == "RIDFRESH"
    assert any("had expired; a new search was submitted" in w for w in out.warnings)


def test_a_saturated_subject_list_cannot_assert_uniqueness(monkeypatch):
    """P1: HITLIST_SIZE caps SUBJECTS, so further primary loci may never have been reported."""
    from zebra.sources import ncbi_blast as nb

    rows = [{"subject": "NT_%06d.1" % i, "identity": 100.0, "align_len": 25, "mismatches": 0,
             "gaps": 0, "s_start": 1, "s_end": 25} for i in range(nb.HITLIST_SIZE - 1)]
    rows.append({"subject": "NC_000007.14", "identity": 100.0, "align_len": 25, "mismatches": 0,
                 "gaps": 0, "s_start": 117639947, "s_end": 117639971})
    row = nb._loci(rows, 25)
    assert row["locus_count"] == 1
    assert row["unique"] is None                 # not True: the list was saturated
    assert "floor" in row["locus_count_is_floor"]
    assert "not determined" in row["unique_note"]


def test_two_nearby_but_distinct_loci_are_not_merged_and_one_locus_is_not_split():
    """P0: a fixed-width bucket both invents and destroys uniqueness; overlap is the test."""
    from zebra.sources import ncbi_blast as nb

    def rows(*starts):
        return [{"subject": "NC_000001.11", "identity": 100.0, "align_len": 20, "mismatches": 0,
                 "gaps": 0, "s_start": s, "s_end": s + 19} for s in starts]

    # 1-20 and 24-43 do not overlap: two distinct loci, 23 bp apart. A 25-wide bucket
    # would have put both in bucket 0 and called the target unique.
    distinct = nb._loci(rows(1, 24), 20)
    assert distinct["locus_count"] == 2 and distinct["unique"] is False
    # 24-43 and 26-45 overlap: one locus seen as two alignments. A 25-wide bucket would
    # have split them across buckets 0 and 1 and called the target non-unique.
    straddling = nb._loci(rows(24, 26), 20)
    assert straddling["locus_count"] == 1 and straddling["unique"] is True
    assert straddling["loci"][0]["alignments"] == 2
    assert straddling["loci"][0]["start"] == 24 and straddling["loci"][0]["end"] == 45
    # far apart: two loci
    apart = nb._loci(rows(1, 5000), 20)
    assert apart["locus_count"] == 2 and apart["unique"] is False


def test_blast_loci_thresholds_are_each_load_bearing():
    """A short-but-exact alignment and a full-length 2-mismatch alignment are both rejected."""
    from zebra.sources import ncbi_blast as nb

    def row(align_len, mismatches, gaps=0, start=1000):
        return {"subject": "NC_000001.11", "identity": 100.0, "align_len": align_len,
                "mismatches": mismatches, "gaps": gaps, "s_start": start, "s_end": start + align_len - 1}

    assert nb._loci([row(25, 0)], 25)["locus_count"] == 1       # exact, full length
    assert nb._loci([row(24, 0)], 25)["locus_count"] == 1       # one base short is allowed
    assert nb._loci([row(20, 0)], 25)["locus_count"] == 0       # five bases short is not
    assert nb._loci([row(25, 1)], 25)["locus_count"] == 1       # one mismatch is allowed
    assert nb._loci([row(25, 2)], 25)["locus_count"] == 0       # two is not
    assert nb._loci([row(25, 0, gaps=1)], 25)["locus_count"] == 0  # a gap is not


def test_the_published_intervals_cannot_be_shortened_by_the_environment(monkeypatch):
    """The 10 s / 60 s intervals are a condition of using the service, not a preference."""
    from zebra.sources import ncbi_blast as nb

    assert nb._interval("ZEBRA_TEST_BLAST_X", 60.0, 60.0) == 60.0
    monkeypatch.setenv("ZEBRA_TEST_BLAST_X", "0")
    assert nb._interval("ZEBRA_TEST_BLAST_X", 60.0, 60.0) == 60.0
    monkeypatch.setenv("ZEBRA_TEST_BLAST_X", "120")
    assert nb._interval("ZEBRA_TEST_BLAST_X", 60.0, 60.0) == 120.0   # longer is allowed
    monkeypatch.setenv("ZEBRA_TEST_BLAST_X", "nonsense")
    assert nb._interval("ZEBRA_TEST_BLAST_X", 60.0, 60.0) == 60.0
    assert nb.SUBMIT_INTERVAL >= 10.0 and nb.POLL_INTERVAL >= 60.0


def test_blast_locate_respects_the_process_deadline(monkeypatch):
    from zebra.sources import ncbi_blast as nb

    monkeypatch.setenv("ZEBRA_DEADLINE_MS", "1")
    monkeypatch.setattr(nb, "submit", lambda q, refresh=False: ("TESTRID", {"db": "NCBI BLAST", "record": "TESTRID"}))
    monkeypatch.setattr(nb, "status", lambda rid: pytest.fail("no poll may be made with no budget left"))
    out = nb.locate({"cand01": "ACGT" * 6}, max_wait=600.0)
    assert out.result["cand01"]["checked"] is False
    assert "time budget" in out.result["cand01"]["reason"]


@pytest.mark.live
def test_live_uniqueness_check_of_a_real_genomic_target():
    """A 25-mer read off GRCh38 must come back at its own coordinates, or say it could not.

    Verified live on 2026-10-06 (728 s, NCBI's queue): the CFTR intron 22 reference
    25-mer came back unique at NC_000007.14:117,639,947-117,639,971 with 0 mismatches,
    with the T2T-CHM13 copy counted as non-primary. The same window carrying the
    patient's C>T came back with no alignment at all -- there is no 16 nt exact seed on
    either side of the substitution -- which is why the screen searches the reference
    reading of a window and why a zero-alignment answer is void rather than a count.
    """
    from zebra.sources import ncbi_blast

    out = ncbi_blast.locate({"cftr25": "CAGTATTAAAATGGCGAGTAAGACA"}, max_wait=600.0)
    row = out.result["cftr25"]
    if not row["checked"]:
        pytest.skip("NCBI BLAST did not finish in the budget: %s" % row["reason"])
    assert row["locus_count"] == 1 and row["unique"] is True
    assert any(l["chrom"] == "7" and abs(l["start"] - 117639947) <= 2 and l["mismatches"] == 0
               for l in row["loci"]), row["loci"]
    assert row["criterion"]


def test_an_unsettled_uniqueness_answer_is_not_scored_as_repeated():
    rank = _rank_for(uniqueness={"checked": True, "unique": None, "locus_count": 1,
                                 "unique_note": "the subject list was saturated"})
    row = next(r for r in rank["components"] if r["name"] == "uniqueness")
    assert row["counted"] is False and row["value"] is None
    assert "not settled" in row["basis"]
    assert "uniqueness" in rank["components_not_counted"]
    assert rank["score_max"] == pytest.approx(_rank_for()["score_max"])


def test_threshold_boundaries_are_where_they_are_documented():
    """0.2 is the ClinGen SVI cut and 0.1 the floor; both are boundary-inclusive as stated."""
    reference = plant({98: "A", 99: "G"}, 300)
    region, reference_region = regions_for(reference, 150, "A")

    def at(value):
        return aso.detect_event(region, reference_region, transcript([delta("DS_AG", value, 101)]), "auto", 0.2)

    # exactly at the floor: not above it, so no site and no event
    assert at(aso.DELTA_FLOOR)["type"] is None
    assert at(aso.DELTA_FLOOR + 0.001)["type"] == "cryptic_acceptor"
    # exactly at min_delta: passes (>=), as Walker 2023 states the threshold
    assert at(0.2)["sites"][0]["passed_min_delta"] is True
    assert at(0.199)["sites"][0]["passed_min_delta"] is False
    # the strongest gain is reported even when it is below the floor
    assert at(0.09)["strongest_gain"] == pytest.approx(0.09)
    assert "below the floor" in at(0.09)["reason"]


def test_hairpin_and_homopolymer_flags_fire_at_their_documented_values():
    stem = aso.HAIRPIN_STEM
    assert aso.hairpin("G" * stem + "ATATA" + "C" * stem)["stem"] == stem
    region = aso.Region("1", 1, 60, 1, "G" * stem + "ATATA" + "C" * stem + "T" * (60 - 2 * stem - 5))
    window = {"region": "pseudoexon_body", "site_kind": None, "site_position": None,
              "i0": 0, "i1": 2 * stem + 4, "length": 2 * stem + 5, "site_premrna_index": None}
    assert aso._candidate(region, window, 0)["hairpin_flag"] is True
    shorter = aso.Region("1", 1, 60, 1, "G" * (stem - 1) + "ATATA" + "C" * (stem - 1) + "T" * 40)
    window2 = dict(window, i1=2 * (stem - 1) + 4, length=2 * (stem - 1) + 5)
    assert aso._candidate(shorter, window2, 0)["hairpin_flag"] is False
    # HOMOPOLYMER_MAX is the longest run still allowed
    run_ok = aso.Region("1", 1, 40, 1, "A" * aso.HOMOPOLYMER_MAX + "CGCGCGCGCGCGCGCGCG"[:20 - aso.HOMOPOLYMER_MAX])
    assert aso._candidate(run_ok, dict(window, i0=0, i1=19, length=20), 0)["homopolymer_flag"] is False
    run_bad = aso.Region("1", 1, 40, 1, "A" * (aso.HOMOPOLYMER_MAX + 1)
                         + "CGCGCGCGCGCGCGCGCG"[:19 - aso.HOMOPOLYMER_MAX])
    assert aso._candidate(run_bad, dict(window, i0=0, i1=19, length=20), 0)["homopolymer_flag"] is True


def test_confidence_reports_a_pseudoexon_with_one_weak_boundary_as_partly_below():
    """A cryptic exon needs BOTH boundaries, so one passing gain does not make it a finding."""
    reference = plant({98: "A", 99: "G", 184: "G", 185: "C"}, 300)
    region, reference_region = regions_for(reference, 185, "T")
    out = aso.detect_event(region, reference_region, transcript([
        delta("DS_AG", 0.11, 101), delta("DS_DG", 0.90, 184)]), "auto", 0.2)
    assert out["type"] == "pseudoexon"
    assert "PARTLY below the ClinGen PP3 threshold" in out["confidence"]
    assert "needs BOTH boundaries" in out["confidence"]
    # a forced single site must be judged on THAT site, not on the strongest gain anywhere
    forced = aso.detect_event(region, reference_region, transcript([
        delta("DS_AG", 0.15, 101), delta("DS_DG", 0.90, 184)]), "cryptic_acceptor", 0.2)
    assert forced["sites"][0]["kind"] == "acceptor"
    assert forced["sites"][0]["passed_min_delta"] is False
    assert "BELOW the ClinGen PP3 threshold" in forced["confidence"]


def test_a_gain_whose_position_cannot_be_used_is_refused_not_reported_as_absent():
    """zebra.s2f sets position=None when SpliceAI's offset is not an integer."""
    reference = plant({}, 300)
    region, reference_region = regions_for(reference, 150, "A")
    unusable = {"score": "DS_DG", "label": "donor gain", "delta": 0.97, "offset": None, "position": None}
    out = aso.detect_event(region, reference_region, transcript([unusable]), "auto", 0.2)
    assert out["type"] is None and out["stop"] is True
    assert "no usable position" in out["reason"]
    assert "NOT an absence of aberrant splicing" in out["reason"]
    assert out["strongest_gain"] == pytest.approx(0.97)
    assert out["gains_without_a_usable_position"][0]["delta"] == pytest.approx(0.97)


def test_a_malformed_splicing_result_costs_the_command_not_the_process(monkeypatch):
    reference = plant({}, 400)
    for broken in ("not a list", [["nested"]], [{"score": "DS_AG", "delta": "abc", "position": 101}]):
        predicted = stub_predict([], pos=186, ref="T", alt="A")
        predicted.result["models"][0]["transcript"]["deltas"] = broken
        wire(monkeypatch, predicted, reference, 1)
        out = aso.screen("7-186-T-A", with_precedents=False)
        assert out.result["stop"] is True, broken
        assert out.result["event"]["type"] is None


def test_an_unreadable_strand_stops_the_screen_instead_of_assuming_plus(monkeypatch):
    reference = plant({98: "A", 99: "G"}, 400)
    for bad in (None, "", "minus", 0, []):
        predicted = stub_predict([delta("DS_AG", 0.4, 101)], pos=186, ref="T", alt="A")
        predicted.result["models"][0]["transcript"]["strand"] = bad
        wire(monkeypatch, predicted, reference, 1)
        out = aso.screen("7-186-T-A", with_precedents=False)
        assert out.result["stop"] is True, bad
        assert "no readable strand" in out.result["event"]["reason"]
        assert out.result["candidates"] == []
        assert any("no readable transcript strand" in w for w in out.warnings)
    # and the two spellings that ARE readable still work
    for good, expect in (("+", 1), ("1", 1), ("-", -1), ("-1", -1)):
        assert aso._strand_of(good) == expect


def test_a_contig_edge_never_produces_an_empty_target(monkeypatch):
    """P0: a donor predicted beyond the sequence that came back must not yield '' targets."""
    reference = plant({98: "A", 99: "G"}, 260)   # Ensembl serves 260 bases of the 400 asked for

    def short_sequence(chrom, lo, hi, assembly="GRCh38", strand=1):
        return Outcome(reference[:max(0, min(len(reference), hi - lo + 1))],
                       sources=[{"db": "Ensembl sequence", "record": "x", "url": "x"}])

    predicted = stub_predict([delta("DS_AG", 0.4, 101), delta("DS_DG", 0.5, 351)], pos=150, ref="T", alt="A")
    from zebra.sources import ensembl
    monkeypatch.setattr(ensembl, "sequence", short_sequence)
    monkeypatch.setattr(aso.ensembl, "sequence", short_sequence)
    from zebra import s2f
    monkeypatch.setattr(s2f, "predict", lambda *a, **k: predicted)
    out = aso.screen("7-150-T-A", lengths=(20,), top=20, with_precedents=False)
    for c in out.result["candidates"]:
        assert len(c["target_premrna_5to3"]) == c["length"], c
        assert c["target_complete"] is True
        assert c["gc_fraction"] is not None
        assert c["rank"]["fraction"] <= 1.0
    assert any("contig edge" in w for w in out.warnings)


def test_a_window_covering_a_non_canonical_site_is_marked_on_the_candidate(monkeypatch):
    """P1: the sequence-check verdict has to reach --json .candidates, not only the prose."""
    reference = plant({}, 400)   # no AG anywhere: the predicted acceptor is not a splice site
    predicted = stub_predict([delta("DS_AG", 0.5, 101)], pos=186, ref="T", alt="A")
    wire(monkeypatch, predicted, reference, 1)
    out = aso.screen("7-186-T-A", lengths=(20,), top=5, with_precedents=False)
    assert "NOT canonical" in out.result["event"]["sequence_check"]
    assert out.result["candidates"]
    for c in out.result["candidates"]:
        assert c["site_canonical"] is False
        assert "unverified" in c["sequence_check_failed"]


def test_the_allele_applied_field_says_what_actually_happened(monkeypatch):
    reference = plant({98: "A", 99: "G"}, 400)
    # a single-base substitution: applied
    predicted = stub_predict([delta("DS_AG", 0.4, 101)], pos=186, ref="T", alt="A")
    wire(monkeypatch, predicted, reference, 1)
    out = aso.screen("7-186-T-A", lengths=(20,), with_precedents=False)
    assert out.result["sequence_window"]["patient_allele_applied"].startswith("yes:")
    # a deletion: not applied, and the field must not claim otherwise
    predicted = stub_predict([delta("DS_AG", 0.4, 101)], pos=186, ref="TT", alt="T")
    wire(monkeypatch, predicted, reference, 1)
    out = aso.screen("7-186-TT-T", lengths=(20,), with_precedents=False)
    assert out.result["sequence_window"]["patient_allele_applied"].startswith("no:")
    assert any("no allele was applied" in w for w in out.warnings)
    for site in out.result["event"]["sites"]:
        assert site["created_by_variant"] is False
        assert "not determined" in site["created_by_variant_basis"]


def test_a_ref_mismatch_is_caught_for_a_multi_base_ref(monkeypatch):
    reference = plant({}, 400)   # filler T
    predicted = stub_predict([delta("DS_AG", 0.4, 101)], pos=186, ref="GG", alt="A")
    wire(monkeypatch, predicted, reference, 1)
    with pytest.raises(UsageError) as err:
        aso.screen("7-186-GG-A", with_precedents=False)
    assert "REF mismatch" in str(err.value)


def test_top_and_min_delta_are_validated_not_clamped():
    for bad_top in (0, -1, 61, "x"):
        with pytest.raises(UsageError):
            aso.screen("7-186-C-T", top=bad_top)
    for bad_delta in (-5, 0.0, aso.DELTA_FLOOR, 1.5, "x"):
        with pytest.raises(UsageError):
            aso.screen("7-186-C-T", min_delta=bad_delta)


def test_exon_skip_on_a_minus_strand_gene_puts_the_acceptor_at_the_higher_coordinate(monkeypatch):
    """On a minus-strand gene the acceptor is the exon's END and the donor its START."""
    # build the window so the minus-strand pre-mRNA has AG before the acceptor and GT after the donor
    length = 400
    seq = ["T"] * length
    # exon 101..184 on the minus strand: acceptor at 184, donor at 101
    # pre-mRNA index of 184 is end-184; the two bases 5' of it are at 186, 185 -> must be revcomp("AG")
    # pre-mRNA index of 184 is 400-184=216; its two 5' neighbours are indices 214, 215,
    # i.e. genomic 186 and 185 complemented -> forward 186 must be T and 185 must be C
    seq[186 - 1] = "T"
    seq[185 - 1] = "C"
    # pre-mRNA index of 101 is 400-101=299; its two 3' neighbours are indices 300, 301,
    # i.e. genomic 100 and 99 complemented -> forward 100 must be C and 99 must be A
    seq[100 - 1] = "C"
    seq[99 - 1] = "A"
    reference = "".join(seq)
    predicted = stub_predict([delta("DS_AG", 0.0, 150)], strand="-", pos=150, ref="T", alt="A")
    wire(monkeypatch, predicted, reference, 1)
    monkeypatch.setattr(aso, "transcript_exons", lambda tid, assembly="GRCh38": Outcome(
        {"transcript": tid, "strand": -1, "chrom": "7",
         "exons": [{"id": "E1", "start": 101, "end": 184, "number": 1, "length": 84}], "exon_count": 1},
        sources=[]))
    out = aso.screen("7-150-T-A", lengths=(20,), event="exon_skip", top=6, with_precedents=False)
    sites = {s["kind"]: s for s in out.result["event"]["sites"]}
    assert sites["acceptor"]["position"] == 184, sites["acceptor"]
    assert sites["donor"]["position"] == 101, sites["donor"]
    assert sites["acceptor"]["dinucleotide_patient_allele"] == "AG"
    assert sites["donor"]["dinucleotide_patient_allele"] == "GT"
    assert "canonical in the sequence" in out.result["event"]["sequence_check"]
    assert out.result["candidate_count"] > 0
    for c in out.result["candidates"]:
        assert c["target_premrna_5to3"] == aso.revcomp(
            reference[c["target_genomic"]["start"] - 1:c["target_genomic"]["end"]])


def test_exon_skip_names_a_terminal_exon_as_more_than_a_frame_question(monkeypatch):
    reference = plant({98: "A", 99: "G", 184: "G", 185: "T"}, 400)
    predicted = stub_predict([delta("DS_AG", 0.0, 101)], pos=150, ref="T", alt="A")
    wire(monkeypatch, predicted, reference, 1)
    monkeypatch.setattr(aso, "transcript_exons", lambda tid, assembly="GRCh38": Outcome(
        {"transcript": tid, "strand": 1, "chrom": "7",
         "exons": [{"id": "E1", "start": 101, "end": 184, "number": 1, "length": 84}], "exon_count": 3},
        sources=[]))
    out = aso.screen("7-150-T-A", lengths=(20,), event="exon_skip", with_precedents=False)
    assert "FIRST exon" in out.result["event"]["frame"]["note"]
    assert "start codon" in out.result["event"]["frame"]["note"]


def test_body_windows_never_leave_the_sequence_that_came_back():
    """The donor's predicted index can lie past the end of a truncated Ensembl response."""
    region = aso.Region("1", 1, 120, 1, plant({}, 120))   # only 120 bases came back
    pseudoexon = {"acceptor": site("acceptor", 20), "donor": site("donor", 300), "size_bp": 281}
    windows = aso._body_windows(region, pseudoexon, (20,), stride=1)
    assert windows, "windows inside the sequence must still be offered"
    for w in windows:
        assert region.holds(w["i0"]) and region.holds(w["i1"]), w
        assert len(region.slice(w["i0"], w["i1"])) == w["length"]


def test_exon_skip_is_not_refused_because_a_gain_has_no_usable_position(monkeypatch):
    """Exon skipping rests on the transcript structure, not on the predicted gains."""
    reference = plant({98: "A", 99: "G", 184: "G", 185: "T"}, 400)
    predicted = stub_predict([{"score": "DS_DG", "label": "donor gain", "delta": 0.9,
                               "offset": None, "position": None}], pos=150, ref="T", alt="A")
    wire(monkeypatch, predicted, reference, 1)
    monkeypatch.setattr(aso, "transcript_exons", lambda tid, assembly="GRCh38": Outcome(
        {"transcript": tid, "strand": 1, "chrom": "7",
         "exons": [{"id": "E1", "start": 101, "end": 184, "number": 2, "length": 84}], "exon_count": 5},
        sources=[]))
    out = aso.screen("7-150-T-A", lengths=(20,), event="exon_skip", with_precedents=False)
    assert out.result["event"]["type"] == "exon_skip"
    assert out.result["candidate_count"] > 0
    # the same prediction in auto mode refuses, because there it IS the basis
    wire(monkeypatch, predicted, reference, 1)
    auto = aso.screen("7-150-T-A", lengths=(20,), with_precedents=False)
    assert auto.result["stop"] is True and "no usable position" in auto.result["event"]["reason"]


def test_an_unsettled_uniqueness_is_not_rendered_as_a_settled_one(monkeypatch):
    reference = plant({98: "A", 99: "G", 184: "G", 185: "C"}, 400)
    predicted = stub_predict([delta("DS_DG", 0.40, 184)], pos=186, ref="C", alt="T")
    wire(monkeypatch, predicted, reference, 1)
    from zebra.sources import ncbi_blast

    monkeypatch.setattr(ncbi_blast, "locate", lambda q, max_wait=None: Outcome(
        {name: {"checked": True, "unique": None, "locus_count": 1, "loci": [],
                "unique_note": "the subject list was saturated"} for name in q},
        sources=[{"db": "NCBI BLAST", "record": "R", "url": "u"}]))
    out = aso.screen("7-186-C-T", lengths=(20,), top=2, uniqueness=True, with_precedents=False)
    rendered = out.text
    assert " 1?" in rendered, rendered
    for c in out.result["candidates"]:
        assert "uniqueness" in c["rank"]["components_not_counted"]


def test_uniqueness_searches_the_reference_reading_of_the_window(monkeypatch):
    """The patient's allele is not in the reference genome, so searching for it finds nothing."""
    reference = plant({98: "A", 99: "G", 184: "G", 185: "C"}, 400)
    predicted = stub_predict([delta("DS_DG", 0.40, 184)], pos=186, ref="C", alt="T")
    wire(monkeypatch, predicted, reference, 1)
    from zebra.sources import ncbi_blast

    sent = {}

    def locate(queries, max_wait=None):
        sent.update(queries)
        return Outcome({name: {"checked": True, "unique": True, "locus_count": 1, "loci": [],
                               "non_primary_alignments": 0, "rid": "R"} for name in queries},
                       sources=[{"db": "NCBI BLAST", "record": "R", "url": "u"}])

    monkeypatch.setattr(ncbi_blast, "locate", locate)
    out = aso.screen("7-186-C-T", lengths=(20,), top=4, uniqueness=True, with_precedents=False)
    covering = [c for c in out.result["candidates"] if c["contains_patient_variant"]]
    assert covering, "the donor windows here cover the variant"
    for c in covering:
        assert c["differs_from_reference"] is True
        assert c["reference_target_premrna_5to3"] != c["target_premrna_5to3"]
        # what went to the search is the reference reading, not the patient's
        assert c["reference_target_premrna_5to3"] in sent.values()
        assert c["target_premrna_5to3"] not in sent.values()
        assert "reference reading" in c["uniqueness"]["searched_sequence"]
        assert "property of the genome" in c["uniqueness"]["note"]
    # a window that does not cover the variant is identical in both readings
    same = [c for c in out.result["candidates"] if not c["contains_patient_variant"]]
    for c in same:
        assert c["differs_from_reference"] is False
        assert "reference and patient allele are the same" in c["uniqueness"]["searched_sequence"]
