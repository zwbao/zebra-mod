import pytest

from zebra import acmg


@pytest.mark.parametrize(
    "codes,points,cls",
    [
        (["PVS1", "PS2", "PM2_Supporting"], 13, "Pathogenic"),
        (["PVS1", "PM2_Supporting"], 9, "Likely pathogenic"),
        (["PS3", "PM2_Supporting", "PP3"], 6, "Likely pathogenic"),
        (["PM2_Supporting", "PP3"], 2, "Uncertain significance"),
        (["BS1", "BP4"], -5, "Likely benign"),
        (["BS1", "BS2"], -8, "Benign"),
        (["PS3", "BS1"], 0, "Uncertain significance"),
    ],
)
def test_points(codes, points, cls):
    r = acmg.classify(codes)
    assert r["points"] == points
    assert r["classification"] == cls


def test_ba1_standalone():
    r = acmg.classify(["BA1"])
    assert r["classification"] == "Benign" and r["points"] is None


def test_richards_disagreement_is_reported():
    r = acmg.classify(["PVS1", "PM2_Supporting"])
    assert r["classification_richards_2015"] == "Uncertain significance"
    assert r["agree"] is False and "note" in r


def test_richards_rules():
    assert acmg.classify(["PS1", "PS3"])["classification_richards_2015"] == "Pathogenic"
    assert acmg.classify(["PM1", "PM2", "PM5"])["classification_richards_2015"] == "Likely pathogenic"
    # a lone BS1 meets no benign rule, so the 2015 rules see no contradiction; points weigh it
    r = acmg.classify(["PVS1", "PS3", "BS1"])
    assert r["classification_richards_2015"] == "Pathogenic" and r["classification"] == "Likely pathogenic"
    assert acmg.classify(["PS1", "PM1", "BS1", "BP4"])["classification_richards_2015"] == "Uncertain significance"


def test_pm1_pp3_cap():
    r = acmg.classify(["PM1", "PP3_Strong"])
    assert r["points"] == 4
    assert any("cap" in w for w in r["warnings"])


def test_warnings():
    r = acmg.classify(["PVS1", "PP3", "PP5", "PM2"])
    text = " ".join(r["warnings"])
    assert "double counting" in text and "PP5" in text and "PM2_Supporting" in text


def test_duplicate_counts_once():
    r = acmg.classify(["PM2_Supporting", "PM2"])
    assert r["points"] == 1


@pytest.mark.parametrize("bad", ["PX1", "PVS2", "PM2_Huge", "BA1_Strong", "PP3_StandAlone"])
def test_bad_codes(bad):
    with pytest.raises(ValueError):
        acmg.classify([bad])


def test_suggest_revel_and_frequency():
    s = acmg.suggest_from_data({"revel": 0.95, "consequence": ["missense_variant"], "gnomad_ac": 0, "spliceai_max": 0.0})
    codes = {x["code"] for x in s}
    assert "PP3_Strong" in codes and "PM2_Supporting" in codes and "BP4" not in codes


def test_suggest_benign_side():
    s = acmg.suggest_from_data({"revel": 0.01, "consequence": "missense_variant", "grpmax_af": 0.08, "grpmax_an": 30000})
    codes = {x["code"] for x in s}
    assert "BP4_Strong" in codes and "BA1" in codes


def test_suggest_bs1_from_max_af():
    s = acmg.suggest_from_data({"faf95": 0.001, "max_credible_af": 0.00004, "gnomad_ac": 30})
    assert [x["code"] for x in s] == ["BS1"]


def test_spliceai_not_with_truncating():
    s = acmg.suggest_from_data({"spliceai_max": 0.6, "consequence": ["stop_gained"]})
    assert all(x["code"] != "PP3" for x in s)
