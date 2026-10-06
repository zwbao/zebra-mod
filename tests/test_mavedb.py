"""MaveDB score sets and one variant's score (CP1-20: no calibrated functional data layer).

Offline tests parse captured MaveDB responses (tests/fixtures/mavedb/, built by
tools/aso/fetch_fixtures.py), including the BRCA1 SGE score set's own primary
calibration -- the thing that makes a score readable at all.

Two rules the tests pin, because breaking either returns a wrong answer under an ACMG
label:
  * a score is never given a functional class this module invented;
  * a variant is never matched across numbering systems -- not across transcripts or
    proteins (different accession), and not into a score set whose target is a
    construct (its position 10 is not the gene's position 10).
"""

from __future__ import annotations

import json
import os

import pytest

from zebra.http import Response, SourceError
from zebra.sources import mavedb

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "mavedb")
SGE = "urn:mavedb:00000097-0-2"            # accession-numbered (NM_007294.3)
CONSTRUCT = "urn:mavedb:00001222-a-2"      # a 191 aa BRCT construct
NO_CALIB_404 = '{"detail":"No primary score calibrations found for score set urn:mavedb:x"}'


def text(name):
    path = os.path.join(FIX, name)
    if not os.path.exists(path):  # pragma: no cover - the fixtures are committed
        pytest.skip("missing fixture %s; rebuild with tools/aso/fetch_fixtures.py" % path)
    with open(path) as fh:
        return fh.read()


def wire(monkeypatch, *, search=None, calibration=None, calibration_404=NO_CALIB_404,
         scores=None, scores_pages=None, scores_error=None):
    """Route the three MaveDB calls to captured bodies."""
    search_body = text(search) if search else json.dumps({"scoreSets": [], "numScoreSets": 0})
    calls = {"search": 0, "calibration": 0, "scores": 0}

    def fake_post(url, payload, source, **kw):
        calls["search"] += 1
        assert "score-sets/search" in url
        return Response(url, 200, search_body, "2026-10-06T00:00:00+00:00", True)

    def fake_get(url, source, **kw):
        calls["calibration"] += 1
        assert "score-calibrations" in url, url
        if calibration is None:
            return Response(url, 404, calibration_404, "2026-10-06T00:00:00+00:00", True)
        body = calibration if calibration.startswith("{") else text(calibration)
        return Response(url, 200, body, "2026-10-06T00:00:00+00:00", True)

    def fake_request(url, source, **kw):
        calls["scores"] += 1
        if scores_error is not None:
            raise scores_error
        params = kw.get("params") or {}
        if scores_pages is not None:
            body = scores_pages(int(params.get("start") or 0), int(params.get("limit") or 0))
        else:
            body = text(scores) if scores else "accession,hgvs_nt,hgvs_pro,score\n"
        bad = (kw.get("validate") or (lambda t: None))(body)
        if bad:
            raise SourceError("MaveDB scores", url, 200, bad)
        return Response("%s?start=%s" % (url, params.get("start")), 200, body,
                        "2026-10-06T00:00:00+00:00", True)

    monkeypatch.setattr(mavedb, "post_json", fake_post)
    monkeypatch.setattr(mavedb, "get_json", fake_get)
    monkeypatch.setattr(mavedb, "request", fake_request)
    return calls


def band(label, lo, hi, criterion=None, *, lower_inclusive=False, upper_inclusive=False,
         classification="normal", strength="STRONG", odds=None):
    return {"label": label, "functional_classification": classification, "range": [lo, hi],
            "range_as_given": [lo, hi], "is_score_band": lo is not None or hi is not None, "class": None,
            "inclusive_lower_bound": lower_inclusive, "inclusive_upper_bound": upper_inclusive,
            "acmg_criterion": criterion, "acmg_evidence_strength": strength if criterion else None,
            "odds_path": odds, "positive_likelihood_ratio": None, "variant_count": 1, "description": None}


def calib(*bands, **kw):
    out = {"title": kw.get("title", "test calibration"), "research_use_only": kw.get("research_use_only", False),
           "baseline_score": 0.0, "baseline_score_description": None, "notes": None,
           "classifications": list(bands)}
    out["score_bands"] = sum(1 for c in out["classifications"] if c["is_score_band"])
    return out


# ---------------------------------------------------------------- HGVS normalisation


def test_protein_hgvs_spellings_all_normalise_to_mavedbs_own():
    for spelling in ("p.R1699W", "p.Arg1699Trp", "NP_009225.1:p.Arg1699Trp", "p.(Arg1699Trp)",
                     "p.(R1699W)", "p.r1699w", "p.ARG1699TRP", "R1699W", "p. R1699W"):
        assert mavedb.normalise_hgvs(spelling) == "p.Arg1699Trp", spelling
    # the legacy one-letter stop and the modern one agree
    for spelling in ("p.R1699X", "p.R1699*", "p.Arg1699Ter", "p.R1699Ter"):
        assert mavedb.normalise_hgvs(spelling) == "p.Arg1699Ter", spelling


def test_coding_hgvs_keeps_its_prefix_and_ignores_case_and_accession():
    for spelling in ("NM_007294.3:c.5565A>T", "c.5565A>T", "c.5565a>t", "C.5565A>T"):
        assert mavedb.normalise_hgvs(spelling) == "c.5565A>T", spelling
    assert mavedb.normalise_hgvs("g.100A>T") == "g.100A>T"
    assert mavedb.normalise_hgvs(None) is None
    assert mavedb.normalise_hgvs("") is None


def test_an_hgvs_that_identifies_no_single_variant_is_refused():
    """p.= means "the protein is unchanged": it matches any number of rows and names none."""
    for spelling in ("p.=", "p.(=)", "NP_1:p.=", "p."):
        assert mavedb.normalise_change(mavedb.split_hgvs(spelling)[1]) is None, spelling


def test_split_hgvs_separates_the_accession():
    assert mavedb.split_hgvs("NM_007294.3:c.5565A>T") == ("NM_007294.3", "c.5565A>T")
    assert mavedb.split_hgvs("c.5565A>T") == (None, "c.5565A>T")
    assert mavedb.split_hgvs("") == (None, None)


def test_matching_refuses_a_different_transcript_or_gene():
    """Coding numbering is accession-specific: the same c. change on another transcript is another variant."""
    cell = "NM_007294.3:c.5565A>T"
    assert mavedb.match_hgvs(cell, "NM_007294.3:c.5565A>T") == "accession"
    assert mavedb.match_hgvs(cell, "NM_007294.1:c.5565A>T") == "accession"     # version-insensitive
    assert mavedb.match_hgvs(cell, "NM_007300.4:c.5565A>T") is None            # other BRCA1 transcript
    assert mavedb.match_hgvs(cell, "ENST00000357654.9:c.5565A>T") is None
    assert mavedb.match_hgvs(cell, "NM_000059.4:c.5565A>T") is None            # BRCA2
    # with no accession on one side the match stands but is labelled, not asserted
    assert mavedb.match_hgvs(cell, "c.5565A>T") == "change_only"
    assert mavedb.match_hgvs("p.Glu10Lys", "p.E10K") == "change_only"
    assert mavedb.match_hgvs("p.Glu10Lys", "p.E11K") is None
    assert mavedb.match_hgvs("NA", "p.E10K") is None


# ---------------------------------------------------------------- score sets


def test_search_records_whether_a_target_is_accession_or_construct_numbered(monkeypatch):
    """P0: a 191 aa BRCT construct numbers from its own position 1, not BRCA1's."""
    wire(monkeypatch, search="score_sets_search_brca1.json")
    out = mavedb.search_score_sets("BRCA1")
    by_urn = {d["urn"]: d for d in out.result}
    assert by_urn[SGE]["numbering"] == "accession"
    assert by_urn[SGE]["target_accession"] == "NM_007294.3"
    assert by_urn[SGE]["searchable_by_hgvs"] is True
    assert by_urn[CONSTRUCT]["numbering"] == "construct"
    assert by_urn[CONSTRUCT]["searchable_by_hgvs"] is False
    assert by_urn[CONSTRUCT]["target_sequence_length"] == 191
    assert "numbered from position 1 of that construct" in by_urn[CONSTRUCT]["numbering_note"]
    for d in out.result:
        assert d["target"] == "BRCA1"
        assert d["url"].startswith("https://mavedb.org/score-sets/")
    counts = [d["num_variants"] or 0 for d in out.result]
    assert counts == sorted(counts, reverse=True)


def test_search_drops_a_score_set_whose_other_target_matched(monkeypatch):
    body = json.dumps({"scoreSets": [
        {"urn": "urn:mavedb:1-a-1", "title": "many targets", "numVariants": 10,
         "targetGenes": [{"name": "BRAF"}, {"name": "KRAS"}]},
        {"urn": "urn:mavedb:2-a-1", "title": "the right one", "numVariants": 20,
         "targetGenes": [{"name": "BRCA1", "targetAccession": {"accession": "NM_007294.3"}}]},
    ], "numScoreSets": 2})
    monkeypatch.setattr(mavedb, "post_json",
                        lambda url, payload, source, **kw: Response(url, 200, body, "t", True))
    out = mavedb.search_score_sets("BRCA1")
    assert [d["urn"] for d in out.result] == ["urn:mavedb:2-a-1"]


def test_search_warns_when_the_count_it_declared_does_not_match(monkeypatch):
    body = json.dumps({"scoreSets": [
        {"urn": "urn:mavedb:2-a-1", "title": "t", "numVariants": 20,
         "targetGenes": [{"name": "BRCA1", "targetAccession": {"accession": "NM_007294.3"}}]}], "numScoreSets": 9})
    monkeypatch.setattr(mavedb, "post_json",
                        lambda url, payload, source, **kw: Response(url, 200, body, "t", True))
    out = mavedb.search_score_sets("BRCA1")
    assert any("said it had 9" in w for w in out.warnings)


def test_search_rejects_an_empty_gene():
    with pytest.raises(ValueError):
        mavedb.search_score_sets("  ")


# ---------------------------------------------------------------- calibration


def test_calibration_is_read_as_the_dataset_stated_it(monkeypatch):
    wire(monkeypatch, calibration="calibration_brca1_sge.json")
    out = mavedb.calibration(SGE)
    c = out.result
    assert c["title"]
    labels = {b["label"] for b in c["classifications"]}
    assert {"Functional", "Non-functional"} <= labels
    functional = next(b for b in c["classifications"] if b["label"] == "Functional")
    nonfunctional = next(b for b in c["classifications"] if b["label"] == "Non-functional")
    assert functional["acmg_criterion"] == "BS3"
    assert nonfunctional["acmg_criterion"] == "PS3"
    assert functional["acmg_evidence_strength"] and nonfunctional["acmg_evidence_strength"]
    assert nonfunctional["odds_path"] and nonfunctional["odds_path"] > 1
    assert c["score_bands"] >= 2
    assert c["research_use_only"] is False


def test_a_score_set_without_a_calibration_gives_none_not_a_guess(monkeypatch):
    wire(monkeypatch, calibration=None)
    out = mavedb.calibration("urn:mavedb:00000081-a-2")
    assert out.result is None
    assert out.sources and out.sources[0]["db"] == "MaveDB score calibration"


def test_a_404_about_the_score_set_is_a_failure_not_an_absence_of_calibration(monkeypatch):
    """MaveDB answers a missing calibration, a missing score set and a bad route all with 404."""
    wire(monkeypatch, calibration=None, calibration_404='{"detail":"score set with URN urn:x not found"}')
    with pytest.raises(SourceError):
        mavedb.calibration("urn:x")
    wire(monkeypatch, calibration=None, calibration_404='{"detail":"Not Found"}')
    with pytest.raises(SourceError):
        mavedb.calibration("urn:x")


def test_calibration_keeps_only_finite_bounds(monkeypatch):
    body = json.dumps({"title": "t", "functionalClassifications": [
        {"label": "numeric", "range": [-1.0, 1.0], "inclusiveLowerBound": True},
        {"label": "categorical", "range": None, "class": "LoF",
         "acmgClassification": {"criterion": "PS3", "evidenceStrength": "STRONG"}},
        {"label": "short range", "range": [], "acmgClassification": {"criterion": "BS3"}},
    ]})
    wire(monkeypatch, calibration=body)
    c = mavedb.calibration("urn:x").result
    by_label = {b["label"]: b for b in c["classifications"]}
    assert by_label["numeric"]["is_score_band"] is True
    assert by_label["categorical"]["is_score_band"] is False and by_label["categorical"]["class"] == "LoF"
    assert by_label["short range"]["is_score_band"] is False
    assert c["score_bands"] == 1


# ---------------------------------------------------------------- score -> class


def test_interpret_respects_the_bounds_the_calibration_declared():
    c = calib(
        band("Functional", -0.748, None, "BS3", lower_inclusive=True, classification="normal", odds=0.02),
        band("Intermediate", -1.328, -0.748, None, classification="not_specified"),
        band("Non-functional", None, -1.328, "PS3", upper_inclusive=True, classification="abnormal", odds=52.4),
    )
    assert mavedb.interpret(0.0, c)["label"] == "Functional"
    assert mavedb.interpret(-0.748, c)["label"] == "Functional"       # inclusive lower bound
    assert mavedb.interpret(-1.0, c)["label"] == "Intermediate"
    assert mavedb.interpret(-1.328, c)["label"] == "Non-functional"   # inclusive upper bound
    assert mavedb.interpret(-5.0, c)["label"] == "Non-functional"
    assert mavedb.interpret(-1.0, c)["acmg_criterion"] is None        # no criterion for the middle band
    assert mavedb.interpret(0.0, c)["source"] == "MaveDB primary score calibration"


def test_interpret_invents_nothing_it_cannot_read():
    c = calib(band("X", 0.0, 1.0, "BS3", lower_inclusive=True, upper_inclusive=True))
    assert mavedb.interpret(-3.0, None) is None                        # no calibration at all
    # a score outside every band is said to be outside them, not pushed into the nearest
    out = mavedb.interpret(5.0, c)
    assert out["label"] is None and "no band" in out["note"]


def test_a_band_with_no_finite_bound_cannot_contain_a_score():
    """A categorical calibration must not label every score with its first class."""
    c = calib(band("LoF", None, None, "PS3", classification="abnormal"),
              band("WT-like", None, None, "BS3", classification="normal"))
    c["classifications"][0]["class"] = "LoF"
    for score in (0.5, -99.0, 0.0):
        out = mavedb.interpret(score, c)
        assert out["label"] is None
        assert "no numeric score band" in out["note"]
        assert out.get("acmg_criterion") is None


def test_a_non_finite_score_is_not_placed_in_any_band():
    c = calib(band("Functional", -0.748, None, "BS3", lower_inclusive=True),
              band("Non-functional", None, -1.328, "PS3", upper_inclusive=True))
    for score in (float("nan"), float("inf"), float("-inf"), None):
        out = mavedb.interpret(score, c)
        assert out["label"] is None, score
        assert out.get("acmg_criterion") is None
        assert "not a finite number" in out["note"]


def test_bands_that_overlap_and_disagree_give_no_class():
    c = calib(band("Normal", -1.0, None, "BS3", lower_inclusive=True, classification="normal"),
              band("Abnormal", None, 0.0, "PS3", upper_inclusive=True, classification="abnormal"))
    out = mavedb.interpret(-0.5, c)
    assert out["label"] is None and out.get("acmg_criterion") is None
    assert "do not agree" in out["note"]
    # reversing the order gives the same refusal, not the other verdict
    reversed_calib = calib(*reversed(c["classifications"]))
    assert mavedb.interpret(-0.5, reversed_calib)["label"] is None


def test_a_research_use_only_calibration_withholds_its_criterion():
    c = calib(band("Functional", -1.0, None, "BS3", lower_inclusive=True), research_use_only=True)
    out = mavedb.interpret(0.0, c)
    assert out["label"] == "Functional"
    assert out["acmg_criterion"] is None
    assert out["acmg_criterion_withheld"] == "BS3"
    assert "research-use-only" in out["note"]


# ---------------------------------------------------------------- the variant's score


def test_lookup_finds_a_coding_variant_and_labels_it_with_the_datasets_calibration(monkeypatch):
    """CP1-20: a calibrated MAVE score reaches the user with the dataset's own class."""
    wire(monkeypatch, search="score_sets_search_brca1.json", calibration="calibration_brca1_sge.json",
         scores="scores_brca1_sge_head.csv")
    out = mavedb.lookup("BRCA1", hgvs_c="NM_007294.3:c.5565A>T")
    r = out.result
    assert r["datasets"] and r["scores_searched"] is True
    assert r["datasets_searched"] == [SGE]              # the construct-numbered sets are not searched
    assert len(r["scores"]) == 1
    hit = r["scores"][0]
    assert hit["urn"] == SGE and hit["hgvs"] == "NM_007294.3:c.5565A>T"
    assert hit["matched_column"] == "hgvs_nt" and hit["accession_checked"] is True
    assert hit["numbering"] == "accession"
    assert hit["score"] == pytest.approx(-0.0153, abs=1e-3)
    interp = hit["interpretation"]
    assert interp["label"] == "Functional" and interp["acmg_criterion"] == "BS3"
    assert interp["source"] == "MaveDB primary score calibration"
    assert "Brnich" in r["how_to_use"]
    assert "A raw score is not a functional class" in r["how_to_use"]


def test_a_construct_numbered_score_set_is_never_searched_by_default(monkeypatch):
    """P0: `p.E10K` in a 191 aa BRCT construct is BRCA1 residue 1682, not residue 10."""
    wire(monkeypatch, search="score_sets_search_brca1.json", calibration="calibration_brca1_sge.json",
         scores="scores_brca1_sge_head.csv")
    out = mavedb.lookup("BRCA1", hgvs_p="p.E10K")
    r = out.result
    assert r["scores"] == []
    assert CONSTRUCT in [d["urn"] for d in r["datasets"]]
    assert CONSTRUCT in [item["urn"] for item in r["datasets_not_searched"]]
    assert any("construct-numbered" in item["reason"] for item in r["datasets_not_searched"])
    assert any("construct-numbered and were not searched" in w for w in out.warnings)
    assert "position 1 of THAT construct" in r["numbering_rule"]


def test_a_construct_hit_is_reported_without_a_class_when_asked_for(monkeypatch):
    rows = ("accession,hgvs_nt,hgvs_splice,hgvs_pro,score\n"
            "urn:mavedb:00001222-a-2#82,NA,NA,p.Glu10Lys,-0.031\n")
    wire(monkeypatch, search="score_sets_search_brca1.json", calibration="calibration_brca1_sge.json",
         scores_pages=lambda start, limit: rows if start == 0 else "accession,hgvs_nt,hgvs_pro,score\n")
    out = mavedb.lookup("BRCA1", hgvs_p="p.E10K", include_construct_numbered=True)
    hits = [s for s in out.result["scores"] if s["numbering"] == "construct"]
    assert hits, "the construct set should be searched when the caller opts in"
    for hit in hits:
        assert hit["interpretation"] is None
        assert "construct's residue, not the gene's" in hit["interpretation_note"]


def test_a_match_without_an_accession_on_both_sides_is_labelled(monkeypatch):
    wire(monkeypatch, search="score_sets_search_brca1.json", calibration="calibration_brca1_sge.json",
         scores="scores_brca1_sge_head.csv")
    out = mavedb.lookup("BRCA1", hgvs_c="c.5565A>T")
    hit = out.result["scores"][0]
    assert hit["accession_checked"] is False
    assert "matched on the change alone" in hit["match_note"]
    assert "NM_007294.3" in hit["match_note"]
    assert any("without an accession on both sides" in w for w in out.warnings)


def test_a_variant_on_another_transcript_is_not_matched(monkeypatch):
    wire(monkeypatch, search="score_sets_search_brca1.json", calibration="calibration_brca1_sge.json",
         scores="scores_brca1_sge_head.csv")
    for other in ("NM_007300.4:c.5565A>T", "NM_000059.4:c.5565A>T", "ENST00000357654.9:c.5565A>T"):
        out = mavedb.lookup("BRCA1", hgvs_c=other)
        assert out.result["scores"] == [], other


def test_a_score_from_an_uncalibrated_set_carries_no_class(monkeypatch):
    wire(monkeypatch, search="score_sets_search_brca1.json", calibration=None,
         scores="scores_brca1_sge_head.csv")
    out = mavedb.lookup("BRCA1", hgvs_c="c.5565A>T")
    hit = out.result["scores"][0]
    assert hit["interpretation"] is None
    assert "publishes no primary calibration" in hit["interpretation_note"]
    assert hit["score"] is not None     # the number is still reported


def test_a_missing_score_is_not_reported_as_a_missing_calibration(monkeypatch):
    rows = ("accession,hgvs_nt,hgvs_pro,score\n"
            "urn:x#1,NM_007294.3:c.5565A>T,NA,NA\n")
    wire(monkeypatch, search="score_sets_search_brca1.json", calibration="calibration_brca1_sge.json",
         scores_pages=lambda start, limit: rows if start == 0 else "accession,hgvs_nt,hgvs_pro,score\n")
    out = mavedb.lookup("BRCA1", hgvs_c="NM_007294.3:c.5565A>T")
    hit = out.result["scores"][0]
    assert hit["score"] is None
    assert hit["interpretation"] is None
    assert "no usable score" in hit["interpretation_note"]
    assert "publishes no primary calibration" not in hit["interpretation_note"]


def test_a_non_finite_score_cell_is_read_as_no_score(monkeypatch):
    for cell in ("nan", "NaN", "inf", "-inf", "1e400", "", "NA"):
        rows = "accession,hgvs_nt,hgvs_pro,score\nurn:x#1,NM_007294.3:c.5565A>T,NA,%s\n" % cell
        wire(monkeypatch, search="score_sets_search_brca1.json", calibration="calibration_brca1_sge.json",
             scores_pages=lambda start, limit, b=rows: b if start == 0 else "accession,hgvs_nt,hgvs_pro,score\n")
        out = mavedb.lookup("BRCA1", hgvs_c="NM_007294.3:c.5565A>T")
        hit = out.result["scores"][0]
        assert hit["score"] is None, cell
        assert hit["interpretation"] is None, cell


def test_a_json_error_body_is_not_read_as_an_empty_score_table(monkeypatch):
    """MaveDB answers some failures with a JSON body and HTTP 200; `score` appears in it."""
    bad = '{"detail":"score set not found or scores unavailable"}'
    wire(monkeypatch, search="score_sets_search_brca1.json", calibration="calibration_brca1_sge.json",
         scores_pages=lambda start, limit: bad)
    out = mavedb.lookup("BRCA1", hgvs_c="NM_007294.3:c.5565A>T")
    assert out.result["scores"] == []
    assert any("not the scores CSV" in w for w in out.warnings)
    assert "no score set was searched in full" in out.result["note"]


def test_a_csv_without_a_score_column_is_rejected(monkeypatch):
    bad = "accession,hgvs_nt,hgvs_pro,score_rep1\nurn:x#1,NM_007294.3:c.5565A>T,NA,0.5\n"
    wire(monkeypatch, search="score_sets_search_brca1.json", scores_pages=lambda start, limit: bad)
    out = mavedb.lookup("BRCA1", hgvs_c="NM_007294.3:c.5565A>T")
    assert out.result["scores"] == []
    assert any("no `score` column" in w for w in out.warnings)


def test_a_variant_the_assay_did_not_cover_is_not_assayed_not_no_effect(monkeypatch):
    wire(monkeypatch, search="score_sets_search_brca1.json", calibration="calibration_brca1_sge.json",
         scores="scores_brca1_sge_head.csv")
    out = mavedb.lookup("BRCA1", hgvs_p="NP_009225.1:p.Arg1699Trp")
    assert out.result["scores"] == []
    assert "not assayed there" in out.result["note"]
    assert "no effect" in out.result["note"]


def test_a_truncated_scan_is_excluded_from_the_not_assayed_claim(monkeypatch):
    header = "accession,hgvs_nt,hgvs_pro,score\n"

    def pages(start, limit):
        return header + "\n".join("urn:x#%d,NM_1.1:c.%dA>T,NA,0.5" % (start + i, start + i)
                                  for i in range(limit)) + "\n"

    wire(monkeypatch, search="score_sets_search_brca1.json", scores_pages=pages)
    monkeypatch.setattr(mavedb, "MAX_ROWS", 25)
    monkeypatch.setattr(mavedb, "PAGE_ROWS", 10)
    out = mavedb.lookup("BRCA1", hgvs_c="NM_007294.3:c.999999A>T")
    assert out.result["scores"] == []
    assert out.result["datasets_incomplete"] == [SGE]
    assert any("stopped after reading" in w for w in out.warnings)
    assert "searched only partly or failed" in out.result["note"]
    assert "not assayed there" not in out.result["note"]


def test_pagination_actually_advances(monkeypatch):
    header = "accession,hgvs_nt,hgvs_pro,score\n"
    seen = []

    def pages(start, limit):
        seen.append(start)
        if start >= 20:
            return header
        return header + "\n".join("urn:x#%d,NM_1.1:c.%dA>T,NA,0.5" % (start + i, start + i)
                                  for i in range(limit)) + "\n"

    wire(monkeypatch, search="score_sets_search_brca1.json", scores_pages=pages)
    monkeypatch.setattr(mavedb, "PAGE_ROWS", 10)
    mavedb.variant_rows(SGE, ["NM_1.1:c.15A>T"], max_rows=100)
    assert seen[:3] == [0, 10, 20], seen


def test_rows_marked_NA_are_not_matched(monkeypatch):
    header = "accession,hgvs_nt,hgvs_splice,hgvs_pro,score\n"
    body = header + "urn:x#1,NA,NA,NA,0.5\n"
    wire(monkeypatch, scores_pages=lambda start, limit: body if start == 0 else header)
    out = mavedb.variant_rows(SGE, ["NA", "na"], max_rows=100)
    assert out.result == []


def test_one_failing_dataset_does_not_lose_the_scores_from_the_others(monkeypatch):
    """P1: a 503 on one score set must not discard a hit already found in another."""
    rows = "accession,hgvs_nt,hgvs_pro,score\nurn:x#1,NM_007294.3:c.5565A>T,NA,-0.0153\n"
    header = "accession,hgvs_nt,hgvs_pro,score\n"
    state = {"calls": 0}

    def pages(start, limit):
        state["calls"] += 1
        if state["calls"] > 1 and start == 0:
            raise SourceError("MaveDB scores", "https://example.invalid", 503, "service unavailable")
        return rows if start == 0 else header

    body = json.dumps({"scoreSets": [
        {"urn": "urn:a", "title": "one", "numVariants": 100,
         "targetGenes": [{"name": "BRCA1", "targetAccession": {"accession": "NM_007294.3"}}]},
        {"urn": "urn:b", "title": "two", "numVariants": 50,
         "targetGenes": [{"name": "BRCA1", "targetAccession": {"accession": "NM_007294.3"}}]},
    ], "numScoreSets": 2})
    wire(monkeypatch, calibration=None, scores_pages=pages)
    monkeypatch.setattr(mavedb, "post_json",
                        lambda url, payload, source, **kw: Response(url, 200, body, "t", True))
    out = mavedb.lookup("BRCA1", hgvs_c="NM_007294.3:c.5565A>T")
    assert len(out.result["scores"]) == 1
    assert out.result["scores"][0]["urn"] == "urn:a"
    assert "urn:b" in out.result["datasets_incomplete"]
    assert any("MaveDB scores urn:b unavailable" in w for w in out.warnings)


def test_no_score_set_for_the_gene_is_an_absence_of_an_assay(monkeypatch):
    wire(monkeypatch)
    out = mavedb.lookup("NOTAGENE123", hgvs_p="p.R1W")
    r = out.result
    assert r["datasets"] == [] and r["scores"] == []
    assert "not evidence about any variant" in r["note"]


def test_lookup_without_a_variant_lists_only_the_score_sets(monkeypatch):
    wire(monkeypatch, search="score_sets_search_brca1.json")
    out = mavedb.lookup("BRCA1")
    assert out.result["scores_searched"] is False
    assert out.result["datasets"]
    assert "no variant given" in out.result["note"]


def test_lookup_refuses_an_hgvs_that_names_no_variant(monkeypatch):
    wire(monkeypatch, search="score_sets_search_brca1.json")
    out = mavedb.lookup("BRCA1", hgvs_p="p.(=)")
    assert out.result["scores"] == [] and out.result["scores_searched"] is False
    assert "identifies no single variant" in out.result["note"]


def test_lookup_says_which_score_sets_it_did_not_search(monkeypatch):
    body = json.dumps({"scoreSets": [
        {"urn": "urn:%d" % i, "title": "t%d" % i, "numVariants": 100 - i,
         "targetGenes": [{"name": "BRCA1", "targetAccession": {"accession": "NM_007294.3"}}]}
        for i in range(4)], "numScoreSets": 4})
    wire(monkeypatch, calibration=None, scores_pages=lambda start, limit: "accession,hgvs_nt,hgvs_pro,score\n")
    monkeypatch.setattr(mavedb, "post_json",
                        lambda url, payload, source, **kw: Response(url, 200, body, "t", True))
    out = mavedb.lookup("BRCA1", hgvs_c="c.5565A>T", max_datasets=2)
    assert out.result["datasets_searched"] == ["urn:0", "urn:1"]
    assert any("the variant was looked for in the 2 largest" in w for w in out.warnings)


def test_render_names_the_criterion_the_framework_and_the_numbering_rule(monkeypatch):
    wire(monkeypatch, search="score_sets_search_brca1.json", calibration="calibration_brca1_sge.json",
         scores="scores_brca1_sge_head.csv")
    out = mavedb.lookup("BRCA1", hgvs_c="NM_007294.3:c.5565A>T")
    rendered = mavedb.render(out.result)
    assert SGE in rendered
    assert "BS3" in rendered
    assert "Brnich" in rendered
    assert "construct-numbered: not searched by HGVS" in rendered


def test_render_survives_a_dataset_with_no_urn():
    rendered = mavedb.render({"gene": "X", "dataset_count": 1, "scores": [],
                              "datasets": [{"urn": None, "title": None, "num_variants": None}],
                              "numbering_rule": "n", "how_to_use": "h"})
    assert "MaveDB X" in rendered


# ---------------------------------------------------------------- live


@pytest.mark.live
def test_live_mavedb_brca1_score_sets_and_a_calibrated_score():
    out = mavedb.lookup("BRCA1", hgvs_c="NM_007294.3:c.5074G>A")
    r = out.result
    assert r["datasets"], "MaveDB should have published BRCA1 score sets"
    assert any(d["numbering"] == "accession" for d in r["datasets"])
    assert any(d["numbering"] == "construct" for d in r["datasets"])
    hits = [s for s in r["scores"] if (s.get("interpretation") or {}).get("acmg_criterion")]
    assert hits, "the BRCA1 SGE set publishes a calibration; the score should carry its class"
    assert hits[0]["interpretation"]["acmg_criterion"] in ("PS3", "BS3")
    assert hits[0]["accession_checked"] is True


@pytest.mark.live
def test_live_mavedb_refuses_the_construct_numbering_false_positive():
    """p.E10K must not come back as the BRCT construct's residue 10 under an ACMG label."""
    out = mavedb.lookup("BRCA1", hgvs_p="p.E10K")
    assert out.result["scores"] == []
    assert out.result["datasets_not_searched"]


@pytest.mark.live
def test_live_mavedb_refuses_a_cross_gene_accession():
    out = mavedb.lookup("BRCA1", hgvs_c="NM_000059.4:c.5565A>T")   # a BRCA2 transcript
    assert out.result["scores"] == []


@pytest.mark.live
def test_live_mavedb_gene_with_no_score_set():
    out = mavedb.lookup("NOTAGENE123", hgvs_p="p.R1W")
    assert out.result["datasets"] == []
    assert "not evidence about any variant" in out.result["note"]


def test_in_range_refuses_a_band_with_no_bound_at_all():
    """A band with neither bound would otherwise contain every score, including NaN's neighbours."""
    boundless = band("any", None, None, "PS3")
    assert boundless["is_score_band"] is False
    for score in (0.0, -99.0, 1e9):
        assert mavedb._in_range(score, boundless) is False, score
    # a one-sided band still works
    assert mavedb._in_range(0.5, band("upper", None, 1.0)) is True
    assert mavedb._in_range(2.0, band("upper", None, 1.0)) is False
    assert mavedb._in_range(2.0, band("lower", 1.0, None)) is True
    assert mavedb._in_range(0.5, band("lower", 1.0, None)) is False
