"""zebra repeat: D4Z4 units and DMPK classes by rule (0.3.5)."""

import pytest

from zebra import repeats
from zebra.core import UsageError


@pytest.mark.parametrize("kb,digest,units", [
    (17.0, "EcoRI", 3),      # the case that made this command: 4 by memory, 3 by the formula
    (10.0, "EcoRI", 1),
    (38.0, "EcoRI", 9),      # the conventional upper bound of FSHD1 sizes sits at 9-10 units
    (14.0, "EcoRI/BlnI", 3), # the same allele after BlnI: ~3 kb shorter, same units
    (15.0, "EcoRI/HindIII", 3),
    (16.3, "EcoRI/HindIII", 4),  # (16.3 - 4.7) / 3.3 = 3.52 (NCT04369209); the old 4.8 flank gave 3
])
def test_d4z4_units_by_digest(kb, digest, units):
    out = repeats.d4z4(kb, digest=digest)
    assert out.result["units"] == units
    assert out.sources and all(s.get("url") for s in out.sources)


def test_d4z4_band_and_warnings():
    out = repeats.d4z4(17.0, haplotype="4qA", methylation_pct=16.2)
    assert out.result["band"].startswith("1-3")
    assert "retinal" in out.result["meaning"]
    assert any("10q" in w for w in out.warnings)  # EcoRI alone cannot place the array on 4q
    assert out.result["methylation_pct"] == 16.2
    assert "4qB" in " ".join(repeats.d4z4(17.0, digest="EcoRI/BlnI", haplotype="4qB").warnings)
    assert repeats.d4z4(60.0).result["band"].startswith("11")


def test_d4z4_refuses_nonsense():
    with pytest.raises(UsageError):
        repeats.d4z4(17.0, digest="HindIII only")
    with pytest.raises(UsageError):
        repeats.d4z4(0.5)
    with pytest.raises(UsageError):
        repeats.d4z4(17.0, methylation_pct=150)


@pytest.mark.parametrize("ctg,band", [(20, "normal"), (40, "premutation"), (80, "mild"), (120, "mild / classic overlap"),
                                      (600, "classic"), (1500, "congenital")])
def test_dmpk_classes(ctg, band):
    assert repeats.dmpk(ctg).result["band"].startswith(band)


def test_d4z4_hindiii_flank_from_the_protocol():
    out = repeats.d4z4(16.3, digest="EcoRI/HindIII")
    assert out.result["flank_kb"] == 4.7 and out.result["units_exact"] == 3.52
    assert out.result["band"].startswith("4-6")
    assert any("4.7" in s.get("record", "") for s in out.sources)


@pytest.mark.parametrize("digest", ["HindIII", "hindiii", "BlnI", "blni", " B "])
def test_d4z4_refuses_single_enzyme_digests(digest):
    with pytest.raises(UsageError, match="double digests"):
        repeats.d4z4(17.0, digest=digest)


def test_d4z4_10q_warning_for_hindiii_too():
    def has_10q(**kw):
        return any("10q" in w for w in repeats.d4z4(20.0, **kw).warnings)
    assert has_10q(digest="EcoRI")
    assert has_10q(digest="EcoRI/HindIII")
    assert not has_10q(digest="EcoRI/BlnI")
    assert not has_10q(digest="EcoRI/HindIII", locus="4q")
    assert not has_10q(digest="EcoRI", locus="4q35")
    assert any("not FSHD1" in w for w in repeats.d4z4(20.0, locus="10q").warnings)
    with pytest.raises(UsageError):
        repeats.d4z4(20.0, locus="chr7")


@pytest.mark.parametrize("hap", ["4qB", "B", "4q B", "4qb", " 4QB163 "])
def test_d4z4_haplotype_b_normalised(hap):
    out = repeats.d4z4(17.0, digest="EcoRI/BlnI", haplotype=hap)
    assert out.result["haplotype_read"] == "4qB"
    assert any("4qB" in w for w in out.warnings)


@pytest.mark.parametrize("hap", ["4qA", "A", "4q A", "4qa"])
def test_d4z4_haplotype_a_normalised(hap):
    out = repeats.d4z4(17.0, digest="EcoRI/BlnI", haplotype=hap)
    assert out.result["haplotype_read"] == "4qA"
    assert not any("4qB" in w for w in out.warnings)


def test_d4z4_methylation_recorded_not_interpreted():
    out = repeats.d4z4(60.0, digest="EcoRI/BlnI", haplotype="4qA", methylation_pct=12.0)
    assert out.result["methylation_interpreted"] is False
    note = out.result["methylation_note"]
    assert "not interpreted" in note and "laboratory" in note
    assert all(g in note for g in ("FSHD2", "SMCHD1", "DNMT3B", "LRIF1"))
    assert "not interpreted" in out.text
    assert any("FSHD2" in w for w in out.warnings)  # non-contracted 4qA allele with a methylation value


def test_d4z4_rounds_half_up():
    # (6.8 + 3.3 * 2.5) = 15.05 kb is exactly 2.5 units: banker's rounding gave 2
    assert repeats.d4z4(6.8 + 3.3 * 2.5).result["units"] == 3
    assert repeats._half_up(4.5) == 5 and repeats._half_up(3.49) == 3


def test_d4z4_boundary_mosaic_and_sex_notes():
    out = repeats.d4z4(38.0)
    joined = " ".join(out.warnings)
    assert "10-38 kb" in joined and "one unit lower" in joined
    assert "mosaic" in joined and "women" in joined
    assert "mosaic" in " ".join(repeats.d4z4(20.0, digest="EcoRI/BlnI").warnings)


def test_dmpk_congenital_wording_is_genereviews():
    out = repeats.dmpk(1500)
    assert "childhood" not in out.result["meaning"]
    assert "congenital" in out.result["meaning"] and "NBK1165" in out.result["meaning"]
