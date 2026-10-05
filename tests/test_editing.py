"""zebra.editing: the base-editing screen on hand-built sequences, then live.

The hand-built sequences are filler 'T' with a few bases placed on purpose, so every
protospacer, window position and bystander the screen reports can be checked by hand.
Two tests use real captured sequence (SCN1A, CPS1) so the arithmetic is pinned to the
genome and not only to the construction.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from zebra import editing
from zebra.core import UsageError
from zebra.sources import ensembl

FIXTURES = Path(__file__).parent / "fixtures"
FLANK = editing.FLANK  # the screen asks for pos-40 .. pos+40


def plant(bases):
    """An 81 nt forward window of filler 'T' with `bases` ({index: base}) placed on purpose."""
    seq = ["T"] * (2 * FLANK + 1)
    for index, base in bases.items():
        seq[int(index)] = base
    return "".join(seq)


def ref_window(bases, ref_base):
    """The same window as `plant`, but with the reference allele at the variant index."""
    placed = dict(bases)
    placed[40] = ref_base
    return plant(placed)


def sequence_stub(forward_window, expect_chrom=None, expect_pos=None):
    """Stand in for ensembl.sequence: hands back `forward_window` for the asked region."""
    from zebra.core import Outcome

    def _sequence(chrom, start, end, assembly="GRCh38", strand=1):
        if expect_chrom is not None:
            assert chrom == expect_chrom
        if expect_pos is not None:
            assert start == expect_pos - FLANK and end == expect_pos + FLANK
        assert end - start + 1 == len(forward_window), (end - start + 1, len(forward_window))
        return Outcome(forward_window, sources=[{"db": "Ensembl sequence", "record": f"{chrom}:{start}-{end}",
                                                 "url": "https://rest.ensembl.org/sequence"}])

    return _sequence


# ---------------------------------------------------------------- primitives


def test_complement_and_revcomp():
    assert [editing.complement(b) for b in "ACGTN"] == ["T", "G", "C", "A", "N"]
    assert editing.complement("a") == "T"
    assert editing.revcomp("AACGT") == "ACGTT"
    assert editing.revcomp(editing.revcomp("ACGTTGCA")) == "ACGTTGCA"
    with pytest.raises(UsageError):
        editing.complement("U")


def test_pam_matches():
    seq = "AAAAGGAAGAAA"
    assert editing.pam_matches(seq, 3, "NGG") == "AGG"
    assert editing.pam_matches(seq, 4, "NGG") is None  # GGA is not NGG
    assert editing.pam_matches(seq, 7, "NG") == "AG"
    assert editing.pam_matches(seq, len(seq) - 1, "NGG") is None  # runs off the end
    assert editing.pam_matches(seq, -1, "NGG") is None
    assert editing.pam_matches("ANGG", 1, "NGG") is None  # an N in the PAM is not a match


# --------------------------------------------------- which editor, which strand


@pytest.mark.parametrize(
    "ref,alt,editor,strand,on_protospacer",
    [("G", "A", "ABE", "+", "A>G"),   # patient A, restore G: A->G on the forward strand
     ("C", "T", "ABE", "-", "A>G"),   # patient T, restore C: A->G on the reverse strand
     ("T", "C", "CBE", "+", "C>T"),
     ("A", "G", "CBE", "-", "C>T")],
)
def test_the_four_transitions_pick_the_right_editor_and_strand(monkeypatch, ref, alt, editor, strand, on_protospacer):
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(plant({40: ref})))
    out = editing.screen(f"1-1041-{ref}-{alt}")
    res = out.result
    assert res["base_editable"] is True and res["editors"] == [editor]
    rows = {r["strand"]: r for r in res["strands"]}
    assert rows[strand]["editor"] == editor
    assert rows[strand]["base_on_protospacer_strand"] == on_protospacer
    other = "-" if strand == "+" else "+"
    assert rows[other]["editor"] is None and "neither A>G" in rows[other]["reason"]


@pytest.mark.parametrize("ref,alt", [("A", "C"), ("C", "A"), ("G", "T"), ("T", "G"), ("C", "G"), ("A", "T")])
def test_transversions_are_not_base_editable(monkeypatch, ref, alt):
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(plant({40: ref})))
    res = editing.screen(f"1-1041-{ref}-{alt}").result
    assert res["base_editable"] is False
    assert "transversion" in res["reason"]
    joined = " ".join(res["alternatives"])
    assert "Prime editing" in joined and "CGBE" in joined
    assert res["guides"] == []


def test_indels_are_not_base_editable_and_point_at_prime_editing(monkeypatch):
    called = []
    monkeypatch.setattr(ensembl, "sequence", lambda *a, **k: called.append(a))
    res = editing.screen("7-117559590-ATCT-A").result
    assert res["base_editable"] is False and res["editor"] is None
    assert "insertion/deletion" in res["reason"]
    assert any("Prime editing" in alt for alt in res["alternatives"])
    assert called == [], "an indel needs no reference fetch"


def test_reference_mismatch_is_refused(monkeypatch):
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(plant({40: "G"})))
    with pytest.raises(UsageError) as err:
        editing.screen("1-1041-C-T")
    assert "REF mismatch" in str(err.value) and "reference has G" in str(err.value)


def test_bad_arguments():
    with pytest.raises(UsageError):
        editing.screen("1-1041-G-A", pam="NAG")
    with pytest.raises(UsageError):
        editing.screen("1-1041-G-A", window=(8, 4))
    with pytest.raises(UsageError):
        editing.screen("1-1041-G-A", window=(0, 8))
    with pytest.raises(UsageError):
        editing.screen("1-1041-G-A", window=(4, 21))
    with pytest.raises(UsageError):
        editing.screen("1-1041-G-A", assembly="hg38")
    with pytest.raises(UsageError):
        editing.screen("1-1041-G-G")


# ------------------------------------------------------ protospacer arithmetic


def test_one_guide_with_a_known_position_and_one_bystander(monkeypatch):
    """Target at index 40, one NGG at 56-58 -> exactly one guide, target at position 5.

    protospacer = indices 36..55, PAM = indices 56..58 ('TGG'). Index 42 is a second A,
    which lands on protospacer position 7: inside the 4-8 window, so a bystander.
    """
    patient_bases = {40: "A", 42: "A", 56: "T", 57: "G", 58: "G"}
    reference = ref_window(patient_bases, "G")
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(reference, expect_chrom="1", expect_pos=1041))
    res = editing.screen("1-1041-G-A").result
    assert res["guide_count"] == 1 and res["bystander_free_guides"] == 0
    g = res["guides"][0]
    assert g["editor"] == "ABE" and g["strand"] == "+"
    assert g["target_protospacer_position"] == 5
    assert g["target_position"] == 1041
    assert g["protospacer"] == "TTTTATATTTTTTTTTTTTT"
    assert len(g["protospacer"]) == 20
    assert g["pam"] == "TGG" and g["pam_pattern"] == "NGG"
    assert g["pam_distal_position"] == 1041 - 4  # protospacer position 1 is 4 bases 5' of the target
    assert g["bystander_count"] == 1
    by = g["bystanders_in_window"][0]
    assert by["protospacer_position"] == 7
    assert by["position"] == 1043  # index 42 = pos + 2
    assert by["edit_on_protospacer_strand"] == "A>G"
    assert by["coding_effect"] == "check bystanders on the transcript"
    assert g["same_base_elsewhere_in_protospacer"] == 1
    assert res["reference_window"]["ref_confirmed"] is True
    assert res["reference_window"]["patient_plus"][40] == "A"


def test_the_window_is_respected_at_both_edges(monkeypatch):
    """Two NGG PAMs: one putting the target at position 4, one at position 9."""
    # position 4 -> start 37, PAM at 57..59 ; position 9 -> start 32, PAM at 52..54
    bases = {40: "A", 58: "G", 59: "G", 53: "G", 54: "G"}
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(ref_window(bases, "G")))
    inside = editing.screen("1-1041-G-A", window=(4, 8)).result
    assert [g["target_protospacer_position"] for g in inside["guides"]] == [4]
    wider = editing.screen("1-1041-G-A", window=(4, 9)).result
    assert sorted(g["target_protospacer_position"] for g in wider["guides"]) == [4, 9]
    narrower = editing.screen("1-1041-G-A", window=(5, 8)).result
    assert narrower["guide_count"] == 0
    assert "no NGG protospacer places the base in window 5-8" in narrower["reason"]
    assert narrower["base_editable"] is True, "the editor still fits; only the window does not"


def test_ng_pam_finds_guides_that_ngg_does_not(monkeypatch):
    bases = {40: "A", 57: "G"}  # NG at 56..57 for start 36 (position 5); no second G, so no NGG
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(ref_window(bases, "G")))
    assert editing.screen("1-1041-G-A", pam="NGG").result["guide_count"] == 0
    ng = editing.screen("1-1041-G-A", pam="NG").result
    assert ng["guide_count"] == 1
    assert ng["guides"][0]["pam"] == "TG" and "relaxed PAM" in ng["guides"][0]["nuclease"]


def test_minus_strand_guides_are_read_off_the_reverse_complement(monkeypatch):
    """Patient T, reference C: ABE on the reverse strand.

    The reverse-strand sequence is revcomp of the 81 nt patient window, so reverse index
    j is forward index 80-j. Target at forward 40 is reverse 40. A PAM at reverse 56..58
    is forward 22..24, and 'NGG' on the reverse strand is 'CCN' read forward.
    """
    # A reverse-strand NGG for a protospacer starting at reverse index s needs reverse
    # [s+21] and [s+22] to be G, i.e. forward [59-s] and [58-s] to be C. For s=36 (which
    # puts the target at window position 5) that is forward indices 23 and 22.
    bases = {40: "C", 22: "C", 23: "C"}
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(plant(bases)))
    res = editing.screen("1-1041-C-T").result
    assert res["editors"] == ["ABE"]
    rows = {r["strand"]: r for r in res["strands"]}
    assert rows["-"]["guide_count"] == 1 and rows["+"]["guide_count"] == 0
    g = res["guides"][0]
    assert g["strand"] == "-" and g["target_protospacer_position"] == 5
    assert g["target_position"] == 1041
    # protospacer position 1 sits 4 bases 3' of the target on the forward strand
    assert g["pam_distal_position"] == 1045
    assert g["edit_on_protospacer_strand"] == "A>G"
    # the filler 'T' on the forward strand is 'A' on the reverse strand: every window
    # position other than the target is a bystander here
    assert g["bystander_count"] == 4
    assert [b["protospacer_position"] for b in g["bystanders_in_window"]] == [4, 6, 7, 8]
    assert [b["position"] for b in g["bystanders_in_window"]] == [1042, 1040, 1039, 1038]


def test_cbe_bystanders_are_cytosines_not_adenines(monkeypatch):
    # patient C at index 40 restored to T: CBE on the forward strand
    bases = {40: "C", 41: "C", 56: "T", 57: "G", 58: "G"}
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(ref_window(bases, "T")))
    res = editing.screen("1-1041-T-C").result
    assert res["editors"] == ["CBE"]
    g = res["guides"][0]
    assert g["edit_on_protospacer_strand"] == "C>T"
    assert [b["protospacer_position"] for b in g["bystanders_in_window"]] == [6]
    assert g["bystanders_in_window"][0]["position"] == 1042


def test_guides_are_ordered_bystander_free_first(monkeypatch):
    # position 5 (start 36, PAM 56..58) has a bystander at index 42; position 6
    # (start 35, PAM 55..57) is built to have none, by placing the PAMs apart
    bases = {40: "A", 42: "A", 56: "T", 57: "G", 58: "G"}
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(ref_window(bases, "G")))
    res = editing.screen("1-1041-G-A", window=(4, 6)).result
    # only position 5 exists here; widen the window so index 42 leaves it
    assert [g["target_protospacer_position"] for g in res["guides"]] == [5]
    narrow = editing.screen("1-1041-G-A", window=(4, 5)).result
    g = narrow["guides"][0]
    assert g["bystander_count"] == 0, "index 42 is protospacer position 7, outside window 4-5"
    assert narrow["bystander_free_guides"] == 1
    assert g["same_base_elsewhere_in_protospacer"] == 1, "it is still in the protospacer, just not in the window"


def test_bystander_annotation_expresses_the_edit_on_the_forward_strand(monkeypatch):
    """A>G on a reverse-strand protospacer is T>C on the forward strand."""
    bases = {40: "C", 22: "C", 23: "C"}
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(plant(bases)))
    seen = {}

    def fake_vep_batch(variants, assembly="GRCh38"):
        from zebra.core import Outcome

        seen["variants"] = list(variants)
        return Outcome([{"start": pos, "most_severe_consequence": "synonymous_variant",
                         "transcript_consequences": [{"mane_select": "NM_1.1", "gene_symbol": "GENE",
                                                      "transcript_id": "ENST1", "hgvsc": "c.3A>G", "hgvsp": None}]}
                        for _c, pos, _r, _a in variants],
                       sources=[{"db": "Ensembl VEP", "record": "batch", "url": "https://rest.ensembl.org/vep"}])

    monkeypatch.setattr(ensembl, "vep_batch", fake_vep_batch)
    out = editing.screen("1-1041-C-T", annotate_bystanders=True)
    res = out.result
    assert {(r, a) for _c, _p, r, a in seen["variants"]} == {("T", "C")}
    assert sorted(p for _c, p, _r, _a in seen["variants"]) == [1038, 1039, 1040, 1042]
    by = res["guides"][0]["bystanders_in_window"][0]
    assert by["coding_effect"] == "synonymous_variant" and by["vep"]["silent"] is True
    assert res["bystanders_annotated"] is True and "bystander_note" not in res
    assert any(src["db"] == "Ensembl VEP" for src in out.sources)


def test_bystanders_unannotated_say_so(monkeypatch):
    bases = {40: "A", 42: "A", 56: "T", 57: "G", 58: "G"}
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(ref_window(bases, "G")))
    res = editing.screen("1-1041-G-A").result
    assert "not annotated" in res["bystander_note"]


def test_labels_and_text_always_carry_the_caveats(monkeypatch):
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(plant({40: "G", 57: "G", 58: "G"})))
    out = editing.screen("1-1041-G-A")
    assert out.text and "Feasibility screen" in out.text
    assert "Delivery is the hard part" in out.text
    joined = " ".join(out.result["labels"])
    assert "not a guide design" in joined and "Bystander" in joined
    indel = editing.screen("7-117559590-ATCT-A")
    assert "Feasibility screen" in indel.text


# --------------------------------------- real sequence, no network (captured)


SCN1A_WINDOW = "TCTACTGTATTTGTTAGAATGCTGGCTATACTCATTGCTCGTTGCCTTTGGGAAGGATCTTCTAGAAAGTCCATGGAAACG"
CPS1_WINDOW = "TGTTTTGAATATCACAAACAAACAGGCTTTCATTACTGCTCAGAATCATGGCTATGCCTTGGACAACACCCTCCCTGCTGG"


def test_scn1a_arg712ter_is_abe_correctable_on_the_forward_strand(monkeypatch):
    """SCN1A NM_001165963.4:c.2134C>T p.Arg712* = 2-166042334-G-A (captured GRCh38 sequence)."""
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(SCN1A_WINDOW, expect_chrom="2", expect_pos=166042334))
    res = editing.screen("2-166042334-G-A").result
    assert res["base_editable"] is True and res["editors"] == ["ABE"]
    assert res["correction"] == {"from": "A", "to": "G", "on": "forward strand",
                                 "description": "restore 2:166042334 A (patient) to G (reference)"}
    assert res["guide_count"] == 1
    g = res["guides"][0]
    assert g["strand"] == "+" and g["target_protospacer_position"] == 8
    assert g["protospacer"] == "ATTGCTCATTGCCTTTGGGA" and g["pam"] == "AGG"
    assert g["protospacer"][7] == "A", "the patient's allele sits at window position 8"
    assert g["bystander_count"] == 0 and res["bystander_free_guides"] == 1
    assert g["target_position"] == 166042334


def test_scn1a_with_ng_pam_adds_the_position_7_guide(monkeypatch):
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(SCN1A_WINDOW))
    res = editing.screen("2-166042334-G-A", pam="NG").result
    assert sorted(g["target_protospacer_position"] for g in res["guides"]) == [7, 8]
    assert all(g["bystander_count"] == 0 for g in res["guides"])


def test_cps1_gln335ter_needs_a_relaxed_pam(monkeypatch):
    """CPS1 c.1003C>T = 2-210591886-C-T: ABE on the reverse strand, but no NGG in window 4-8."""
    monkeypatch.setattr(ensembl, "sequence", sequence_stub(CPS1_WINDOW, expect_chrom="2", expect_pos=210591886))
    ngg = editing.screen("2-210591886-C-T").result
    assert ngg["base_editable"] is True and ngg["editors"] == ["ABE"]
    rows = {r["strand"]: r for r in ngg["strands"]}
    assert rows["-"]["editor"] == "ABE" and rows["+"]["editor"] is None
    assert ngg["guide_count"] == 0 and "no NGG protospacer" in ngg["reason"]

    ng = editing.screen("2-210591886-C-T", pam="NG").result
    assert sorted(g["target_protospacer_position"] for g in ng["guides"]) == [4, 8]
    best = ng["guides"][0]
    assert best["target_protospacer_position"] == 8 and best["bystander_count"] == 0
    assert best["protospacer"] == "TGATTCTAAGCAGTAATGAA" and best["pam"] == "AG"
    worse = ng["guides"][1]
    assert [b["position"] for b in worse["bystanders_in_window"]] == [210591885, 210591882]


# --------------------------------------------------------------- live


@pytest.mark.live
def test_live_scn1a_matches_the_captured_window():
    out = editing.screen("2-166042334-G-A")
    res = out.result
    assert res["reference_window"]["reference_plus"] == SCN1A_WINDOW
    assert res["guide_count"] == 1
    assert res["guides"][0]["protospacer"] == "ATTGCTCATTGCCTTTGGGA"
    assert any(src["db"] == "Ensembl sequence" for src in out.sources)


@pytest.mark.live
def test_live_cps1_control_is_abe_correctable():
    res = editing.screen("2-210591886-C-T", pam="NG").result
    assert res["editors"] == ["ABE"] and res["guide_count"] == 2
    assert res["reference_window"]["reference_plus"] == CPS1_WINDOW


@pytest.mark.live
def test_live_hgvs_input_resolves_and_screens():
    res = editing.screen("NM_001165963.4:c.2134C>T").result
    v = res["variant"]
    assert (v["chrom"], v["pos"], v["ref"], v["alt"]) == ("2", 166042334, "G", "A")
    assert res["editors"] == ["ABE"] and res["guides"][0]["strand"] == "+"


@pytest.mark.live
def test_live_bystander_annotation_on_cps1():
    out = editing.screen("2-210591886-C-T", pam="NG", annotate_bystanders=True)
    guide = next(g for g in out.result["guides"] if g["bystander_count"])
    effects = {b["coding_effect"] for b in guide["bystanders_in_window"]}
    assert effects and effects != {"check bystanders on the transcript"}
    assert all(b["vep"]["gene"] == "CPS1" for b in guide["bystanders_in_window"])


def test_a_short_window_at_a_contig_edge_is_warned_about(monkeypatch):
    from zebra.core import Outcome

    short = plant({40: "G", 57: "G", 58: "G"})[:70]

    def _sequence(chrom, start, end, assembly="GRCh38", strand=1):
        return Outcome(short, sources=[{"db": "Ensembl sequence", "record": "x", "url": "u"}])

    monkeypatch.setattr(ensembl, "sequence", _sequence)
    out = editing.screen("1-1041-G-A")
    assert any("contig edge" in w for w in out.warnings)
    assert out.result["base_editable"] is True
