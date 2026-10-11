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


@pytest.mark.parametrize("ctg,band", [(20, "normal"), (40, "premutation"), (80, "mild"), (600, "classic"), (1500, "large")])
def test_dmpk_classes(ctg, band):
    assert repeats.dmpk(ctg).result["band"].startswith(band)
