"""zebra.phenopacket: reading GA4GH Phenopackets and writing a zebra case as one."""

import copy
import json

import pytest

from zebra import phenopacket as P

# Modelled on phenopacket-store 0.1.27 (PMID_36996813_Individual_11.json, trimmed)
STORE_CASE = {
    "id": "PMID_36996813_Individual_11",
    "subject": {"id": "Individual 11", "timeAtLastEncounter": {"gestationalAge": {"weeks": 32, "days": 0}},
                "sex": "UNKNOWN_SEX"},
    "phenotypicFeatures": [
        {"type": {"id": "HP:0000319", "label": "Smooth philtrum"}, "excluded": True},
        {"type": {"id": "HP:0002170", "label": "Intracranial hemorrhage"},
         "onset": {"ontologyClass": {"id": "HP:0030674", "label": "Antenatal onset"}}},
        {"type": {"id": "HP:0000238", "label": "Hydrocephalus"}, "onset": {"age": {"iso8601duration": "P2M"}}},
        {"type": {"id": "MAXO:0000001", "label": "not a phenotype"}},
    ],
    "interpretations": [{
        "id": "x1", "progressStatus": "SOLVED",
        "diagnosis": {
            "disease": {"id": "OMIM:620371", "label": "Neurodevelopmental disorder with intracranial hemorrhage"},
            "genomicInterpretations": [{
                "subjectOrBiosampleId": "Individual 11", "interpretationStatus": "CAUSATIVE",
                "variantInterpretation": {
                    "acmgPathogenicityClassification": "PATHOGENIC",
                    "variationDescriptor": {
                        "id": "v", "geneContext": {"valueId": "HGNC:17474", "symbol": "ESAM"},
                        "expressions": [{"syntax": "hgvs.c", "value": "NM_138961.3:c.35T>A"},
                                        {"syntax": "hgvs.p", "value": "NP_620411.2:p.(Leu12Ter)"}],
                        "vcfRecord": {"genomeAssembly": "hg38", "chrom": "chr11", "pos": 124762120, "ref": "A",
                                      "alt": "T"},
                        "allelicState": {"id": "GENO:0000136", "label": "homozygous"}}}}]}}],
    "diseases": [{"term": {"id": "OMIM:620371", "label": "NDD"}, "onset": {"gestationalAge": {"weeks": 21}}}],
    "metaData": {"phenopacketSchemaVersion": "2.0.2",
                 "externalReferences": [{"id": "PMID:36996813", "description": "ESAM"}]},
}


def test_read_a_phenopacket_store_case():
    r = P.read(STORE_CASE)
    assert [t["id"] for t in r["present"]] == ["HP:0002170", "HP:0000238"]
    assert [t["id"] for t in r["excluded"]] == ["HP:0000319"]
    assert r["present"][0]["onset"] == "HP:0030674 Antenatal onset"
    assert r["present"][1]["onset"] == "P2M"
    assert r["diseases"] == [{"id": "OMIM:620371", "label": "Neurodevelopmental disorder with intracranial hemorrhage",
                              "excluded": False, "from": "interpretation", "status": "SOLVED", "onset": None}]
    assert r["genes"] == [{"symbol": "ESAM", "id": "HGNC:17474", "status": "CAUSATIVE"}]
    v = r["variants"][0]
    assert v["hgvs_c"] == "NM_138961.3:c.35T>A" and v["vcf"] == "chr11-124762120-A-T" and v["zygosity"] == "homozygous"
    assert v["acmg"] == "PATHOGENIC" and v["kind"] == "small"
    assert r["references"] == ["PMID:36996813"]
    assert r["subject"]["age"] == "gestational 32w0d"
    assert any("MAXO:0000001" in w and "not an HPO term" in w for w in r["warnings"])


def test_read_from_a_file_path(tmp_path):
    path = tmp_path / "p.json"
    path.write_text(json.dumps(STORE_CASE), "utf-8")
    assert P.read(str(path))["present"][0]["id"] == "HP:0002170"
    assert P.read(path)["id"] == "PMID_36996813_Individual_11"


def test_read_v1_negated_and_top_level_variants():
    v1 = {"id": "v1", "subject": {"id": "p"},
          "phenotypicFeatures": [{"type": {"id": "HP:0001250", "label": "Seizure"}},
                                 {"type": {"id": "HP:0001252", "label": "Hypotonia"}, "negated": True}],
          "genes": [{"id": "HGNC:10585", "symbol": "SCN1A"}],
          "variants": [{"hgvsAllele": {"hgvs": "NM_001165963.4:c.2134C>T"},
                        "zygosity": {"id": "GENO:0000135", "label": "heterozygous"}}],
          "diseases": [{"term": {"id": "OMIM:607208", "label": "Dravet"}}]}
    r = P.read(v1)
    assert [t["id"] for t in r["present"]] == ["HP:0001250"] and [t["id"] for t in r["excluded"]] == ["HP:0001252"]
    assert r["genes"][0]["symbol"] == "SCN1A" and r["variants"][0]["hgvs_c"].endswith("c.2134C>T")
    assert r["diseases"][0]["from"] == "diseases"


def test_read_conflicting_duplicate_feature_is_reported_not_silently_merged():
    pp = {"phenotypicFeatures": [{"type": {"id": "HP:0001250"}}, {"type": {"id": "HP:0001250"}, "excluded": True}]}
    r = P.read(pp)
    assert [t["id"] for t in r["present"]] == ["HP:0001250"] and r["excluded"] == []
    assert any("both present and excluded" in w for w in r["warnings"])


def test_read_excluded_disease_is_kept_with_its_flag():
    pp = {"diseases": [{"term": {"id": "OMIM:1"}, "excluded": True}], "phenotypicFeatures": []}
    assert P.read(pp)["diseases"][0]["excluded"] is True


@pytest.mark.parametrize("bad, msg", [
    ([1, 2], "file path or a parsed JSON object"),
    ({"foo": 1}, "not a phenopacket"),
])
def test_read_refuses_what_is_not_a_phenopacket(bad, msg):
    with pytest.raises(P.PhenopacketError, match=msg):
        P.read(bad)


def test_read_refuses_bad_files(tmp_path):
    with pytest.raises(P.PhenopacketError, match="cannot read"):
        P.read(tmp_path / "missing.json")
    bad = tmp_path / "bad.json"
    bad.write_text("{not json", "utf-8")
    with pytest.raises(P.PhenopacketError, match="not JSON"):
        P.read(bad)
    arr = tmp_path / "arr.json"
    arr.write_text("[]", "utf-8")
    with pytest.raises(P.PhenopacketError, match="JSON object"):
        P.read(arr)
    big = tmp_path / "big.json"
    big.write_bytes(b" " * (P.MAX_BYTES + 1))
    with pytest.raises(P.PhenopacketError, match="at most"):
        P.read(big)


def test_read_tolerates_wrong_shapes_inside():
    pp = {"phenotypicFeatures": [None, "x", {"type": "HP:0001250"}, {"type": {"id": 5}}],
          "interpretations": [None, {"diagnosis": "x"}, {"diagnosis": {"genomicInterpretations": [None, {}]}}],
          "metaData": {"externalReferences": [None, {"id": 3}]}}
    r = P.read(pp)
    assert r["present"] == [] and r["diseases"] == [] and r["references"] == []


def test_read_family_message_reads_the_proband():
    fam = {"id": "fam", "proband": copy.deepcopy(STORE_CASE)}
    r = P.read(fam)
    assert r["present"][0]["id"] == "HP:0002170" and any("Family" in w for w in r["warnings"])


# ---------------------------------------------------------------- from_case

CASE = {
    "schema": "zebra.case/1", "id": "lily-zhang", "title": "Lily Zhang 张丽丽", "role": "family", "language": "zh",
    "proband": {"sex": "女", "age": "3岁", "ancestry": "Han", "consanguinity": False},
    "phenotypes": [
        {"id": "HP:0002373", "label": "Febrile seizure", "status": "present", "onset": "6 months",
         "source": "records/张丽丽-neuro.pdf", "note": "first at 6 months"},
        {"id": "HP:0001252", "label": "Hypotonia", "status": "excluded", "onset": None, "note": "EEG normal"},
        {"id": "HP:0001263", "label": "Global developmental delay", "status": "present", "onset": "大约一岁时"},
        {"id": "HP:0003593", "label": "Infantile onset", "status": "present", "onset": "婴儿期"},
    ],
    "variants": [
        {"id": "v1", "kind": "small", "gene": "SCN1A", "hgvs_c": "NM_001165963.4:c.2134C>T", "hgvs_p": "p.Arg712*",
         "vcf": "2-166002594-G-A", "assembly": "GRCh38", "zygosity": "heterozygous",
         "classification_lab": "Pathogenic", "source": "records/report.pdf"},
        {"id": "v2", "kind": "cnv", "region": "22:18648855-21800471", "cnv_type": "loss", "zygosity": "het"},
        {"id": "v3", "kind": "repeat_expansion", "gene": "FMR1", "motif": "CGG", "repeat_count": 230},
    ],
    "hypotheses": [
        {"id": "h1", "disease": "Dravet syndrome", "status": "confirmed", "ids": {"ORPHA": "33069", "OMIM": "607208"}},
        {"id": "h2", "disease": "GEFS+", "status": "excluded", "ids": {"ORPHA": "36387"}},
        {"id": "h3", "disease": "Something", "status": "leading", "ids": {}},
        {"id": "h4", "disease": "Considered thing", "status": "considered", "ids": {"OMIM": "100100"}},
    ],
    "family": {"members": [{"name": "Zhang Wei"}], "notes": "mother 13800138000"},
    "questions": ["ask about 13800138000"],
    "privacy": {"identifiers": ["张丽丽", "13800138000"]},
}


@pytest.fixture(autouse=True)
def _no_local_hpo(tmp_path, monkeypatch):
    """Labels come from the HPO release when one is installed; these tests pin the no-release path."""
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(tmp_path / "no-hpo"))


def test_from_case_writes_a_v2_phenopacket():
    warnings = []
    pp = P.from_case(CASE, created="2026-10-06T00:00:00Z", hpo_version="2026-09-01", warnings=warnings)
    assert pp["metaData"]["phenopacketSchemaVersion"] == "2.0"
    assert pp["metaData"]["createdBy"].startswith("zebra-mod ")
    assert {r["namespacePrefix"] for r in pp["metaData"]["resources"]} >= {"HP", "GENO", "SO", "OMIM", "ORPHA"}
    hp = next(r for r in pp["metaData"]["resources"] if r["namespacePrefix"] == "HP")
    assert hp["version"] == "2026-09-01"
    assert pp["subject"] == {"id": "proband", "sex": "FEMALE", "timeAtLastEncounter": {"age": {"iso8601duration": "P3Y"}}}
    feats = {f["type"]["id"]: f for f in pp["phenotypicFeatures"]}
    assert feats["HP:0001252"]["excluded"] is True
    assert feats["HP:0002373"]["onset"] == {"age": {"iso8601duration": "P6M"}}
    # an onset that does not parse is never turned into a made-up age, and free text is not exported
    assert "onset" not in feats["HP:0001263"] and "description" not in feats["HP:0001263"]
    assert any("HP:0001263: onset is free text" in w for w in warnings)
    assert feats["HP:0003593"]["onset"] == {"ontologyClass": {"id": "HP:0003593", "label": "Infantile onset"}}
    interp = pp["interpretations"]
    assert interp[0]["progressStatus"] == "SOLVED"
    assert interp[0]["diagnosis"]["disease"]["id"] == "OMIM:607208"
    assert interp[0]["id"] == "interpretation-1"
    gis = interp[0]["diagnosis"]["genomicInterpretations"]
    assert len(gis) == 3 and all(g["interpretationStatus"] == "CANDIDATE" for g in gis)
    small = gis[0]["variantInterpretation"]
    assert small["acmgPathogenicityClassification"] == "PATHOGENIC"
    vd = small["variationDescriptor"]
    assert vd["vcfRecord"] == {"genomeAssembly": "GRCh38", "chrom": "2", "pos": 166002594, "ref": "G", "alt": "A"}
    assert vd["allelicState"]["id"] == "GENO:0000135" and vd["geneContext"]["symbol"] == "SCN1A"
    cnv = gis[1]["variantInterpretation"]["variationDescriptor"]
    assert cnv["structuralType"]["id"] == "SO:0001743" and "22:18648855-21800471" in cnv["description"]
    rep = gis[2]["variantInterpretation"]["variationDescriptor"]
    assert rep["structuralType"]["id"] == "SO:0002162"
    assert {"term": {"id": "ORPHA:36387", "label": "GEFS+"}, "excluded": True} in pp["diseases"]
    assert any("leading hypothesis has no OMIM/ORPHA/MONDO id" in w for w in warnings)


def test_from_case_carries_no_identifying_text():
    text = json.dumps(P.from_case(CASE, created="2026-10-06T00:00:00Z", hpo_version="x"), ensure_ascii=False)
    for secret in ("lily", "Lily", "张丽丽", "Zhang Wei", "13800138000", "records/", "EEG normal", "first at 6 months",
                   "Han"):
        assert secret not in text, secret
    a, b = P.from_case(CASE, created="t", hpo_version="x")["id"], P.from_case(CASE, created="t", hpo_version="x")["id"]
    assert a.startswith("zebra-") and a != b, "random, so it cannot be reversed to the case (folder) name"


def test_from_case_round_trips_through_read():
    pp = P.from_case(CASE, created="2026-10-06T00:00:00Z", hpo_version="2026-09-01")
    back = P.read(json.loads(json.dumps(pp)))
    assert [t["id"] for t in back["present"]] == ["HP:0002373", "HP:0001263", "HP:0003593"]
    assert [t["id"] for t in back["excluded"]] == ["HP:0001252"]
    assert back["present"][0]["onset"] == "P6M"
    assert back["present"][1]["onset"] is None and "description" not in back["present"][1]
    dis = {d["id"]: d for d in back["diseases"]}
    assert dis["OMIM:607208"]["status"] == "SOLVED" and dis["OMIM:607208"]["excluded"] is False
    assert dis["ORPHA:36387"]["excluded"] is True
    assert [g["symbol"] for g in back["genes"]] == ["SCN1A", "FMR1"]
    kinds = [v["kind"] for v in back["variants"]]
    assert kinds == ["small", "cnv", "repeat_expansion"]
    assert back["variants"][0]["hgvs_c"] == "NM_001165963.4:c.2134C>T"
    assert back["variants"][0]["zygosity"] == "heterozygous"
    assert back["subject"]["sex"] == "FEMALE" and back["subject"]["age"] == "P3Y"


def test_from_case_variants_without_a_diagnosis_still_export():
    case = {"id": "c", "phenotypes": [{"id": "HP:0001250", "label": "Seizure"}],
            "variants": [{"id": "v1", "gene": "SCN1A", "hgvs_c": "c.1A>G"}], "hypotheses": []}
    warnings = []
    pp = P.from_case(case, created="t", hpo_version="x", warnings=warnings)
    assert pp["interpretations"][0]["progressStatus"] == "IN_PROGRESS"
    # Diagnosis.disease is required in v2: the root class "disease" says no more than is known ...
    assert pp["interpretations"][0]["diagnosis"]["disease"] == {"id": "MONDO:0000001", "label": "disease"}
    assert any("without a diagnosis" in w for w in warnings)
    back = P.read(pp)
    assert back["variants"][0]["hgvs_c"] == "c.1A>G"
    assert back["diseases"] == [], "... and read() does not report it as a diagnosis"


def test_from_case_minimal_and_bad_input():
    pp = P.from_case({"id": "x"}, created="t", hpo_version="x")
    assert pp["phenotypicFeatures"] == [] and pp["subject"]["sex"] == "UNKNOWN_SEX"
    with pytest.raises(P.PhenopacketError):
        P.from_case([1])
    warnings = []
    P.from_case({"id": "x", "phenotypes": [{"id": "Seizure"}], "proband": {"age": "about three"}},
                created="t", hpo_version="x", warnings=warnings)
    assert warnings
    assert any("not an HPO id" in w for w in warnings) and any("not an unambiguous age" in w for w in warnings)


@pytest.mark.parametrize("text, iso", [
    ("P1Y6M", "P1Y6M"), ("6 months", "P6M"), ("3岁", "P3Y"), ("18 mo", "P18M"), ("1.5 years", "P18M"),
    ("6个月", "P6M"), ("12 weeks old", "P12W"), ("2.5 months", None), ("about 6 months", None), ("", None),
    ("P", None), (3, None), (True, None), (-1, None),
    ("6月", None), ("2月", None), ("6月龄", "P6M"),  # R-A-P1-11: 6月 is usually June, not six months
])
def test_parse_age_never_estimates(text, iso):
    assert P.parse_age(text) == iso



# ---------------------------------------------------------------- review findings (R-A)

def test_r_a_p0_1_a_protected_identifier_anywhere_refuses_the_export():
    base = {"id": "c", "privacy": {"identifiers": ["张三", "S12345", "2023-05-12"]},
            "phenotypes": [{"id": "HP:0001250", "label": "Seizure", "onset": "2023-05-12 张三 first seizure at 协和医院"}],
            "hypotheses": [{"disease": "Dravet", "status": "confirmed", "ids": {"OMIM": "607208"}}],
            "variants": [{"id": "v1", "kind": "cnv", "region": "1:100-200", "description": "S12345 华大",
                          "method": "协和 CMA"}]}
    warnings = []
    pp = P.from_case(base, created="t", hpo_version="x", warnings=warnings)  # free text is not exported at all
    text = json.dumps(pp, ensure_ascii=False)
    for leaked in ("张三", "S12345", "2023-05-12", "华大", "协和"):
        assert leaked not in text, leaked
    # a label or hypothesis name the person typed is checked, and a hit refuses the export
    for field in ({"phenotypes": [{"id": "HP:0001250", "label": "张三的抽搐"}]},
                  {"hypotheses": [{"disease": "S12345 Dravet", "status": "confirmed", "ids": {"OMIM": "607208"}}]}):
        bad = dict(base, variants=[], **field)
        with pytest.raises(P.PhenopacketError, match="protected identifier"):
            P.from_case(bad, created="t", hpo_version="x")
    with pytest.raises(P.PhenopacketError, match="resident ID"):
        P.from_case({"id": "c", "hypotheses": [{"disease": "x 330106201903140021", "status": "confirmed",
                                                 "ids": {"OMIM": "607208"}}]}, created="t", hpo_version="x")


def test_r_a_p1_5_hand_edited_fields_are_normalised_or_refused():
    warnings = []
    case = {"id": "c", "phenotypes": [{"id": "HP:0001250", "label": "Seizure", "status": "absent"},
                                      {"id": "HP:0001263", "label": "GDD", "status": "excluded"}],
            "variants": [{"id": "v1", "kind": "SNV", "gene": "SCN1A", "hgvs_c": "c.1A>G"}],
            "hypotheses": [{"disease": "x", "status": "confirmed", "ids": {"OMIM": None}},
                           {"disease": "y", "status": "confirmed", "ids": {"OMIM": "60720"}}]}
    pp = P.from_case(case, created="t", hpo_version="x", warnings=warnings)
    assert [f["type"]["id"] for f in pp["phenotypicFeatures"]] == ["HP:0001263"]
    assert pp["phenotypicFeatures"][0]["excluded"] is True
    assert any("'absent' is neither present nor excluded" in w for w in warnings)
    vd = pp["interpretations"][0]["diagnosis"]["genomicInterpretations"][0]["variantInterpretation"]["variationDescriptor"]
    assert "structuralType" not in vd and vd["expressions"][0]["value"] == "c.1A>G"  # SNV is a small variant
    assert "OMIM:None" not in json.dumps(pp) and pp["interpretations"][0]["diagnosis"]["disease"]["id"] == "MONDO:0000001"
    assert any("60720 is not a well-formed OMIM id" in w for w in warnings)


@pytest.mark.parametrize("lab, want", [("P", "PATHOGENIC"), ("LP", "LIKELY_PATHOGENIC"), ("VUS", "UNCERTAIN_SIGNIFICANCE"),
                                       ("LB", "LIKELY_BENIGN"), ("B", "BENIGN"), ("致病性", "PATHOGENIC"),
                                       ("可能致病性", "LIKELY_PATHOGENIC"), ("whatever", "NOT_PROVIDED")])
def test_r_a_p1_6_acmg_abbreviations_are_read(lab, want):
    assert P._acmg({"classification_lab": lab}) == want


def test_r_a_p1_9_unreadable_inputs_raise_phenopacket_error(tmp_path):
    with pytest.raises(P.PhenopacketError):
        P.read(tmp_path)  # a directory
    deep = tmp_path / "deep.json"
    deep.write_text('{"phenotypicFeatures": ' + "[" * 50000 + "]" * 50000 + "}", "utf-8")
    with pytest.raises(P.PhenopacketError):
        P.read(deep)
    nested = {"start": {}}
    cur = nested
    for _ in range(2000):
        cur["ageRange"] = {"start": {}}
        cur = cur["ageRange"]["start"]
    r = P.read({"phenotypicFeatures": [{"type": {"id": "HP:0001250"}, "onset": {"ageRange": nested}}]})
    assert r["present"][0]["id"] == "HP:0001250"


def test_read_is_strict_about_flags_and_shapes():
    r = P.read({"phenotypicFeatures": [{"type": {"id": "HP:0001250"}, "excluded": "false"},
                                       {"type": {"id": "HP:0001252"}, "excluded": float("nan")},
                                       {"type": {"id": "HP:0001263"}, "excluded": "true"}],
                "interpretations": [{"diagnosis": {"genomicInterpretations": [{"variantInterpretation": {
                    "variationDescriptor": {"id": "v", "vcfRecord": {"chrom": "1", "pos": 5}}}}]}}],
                "subject": {"id": 5, "sex": ["x"], "timeAtLastEncounter": {"gestationalAge": {"weeks": {"a": 1}}}}})
    assert [t["id"] for t in r["present"]] == ["HP:0001250", "HP:0001252"]
    assert [t["id"] for t in r["excluded"]] == ["HP:0001263"]
    assert "vcf" not in r["variants"][0]
    assert r["subject"] == {"id": None, "sex": None, "age": None}
    with pytest.raises(P.PhenopacketError, match="Cohort"):
        P.read({"id": "cohort", "members": [STORE_CASE]})


def test_from_case_does_not_guess_assembly_or_loss_or_zygosity():
    warnings = []
    case = {"id": "c", "variants": [
        {"id": "v1", "kind": "small", "gene": "SCN1A", "vcf": "2-166002594-G-A", "zygosity": "mosaic"},
        {"id": "v2", "kind": "copy_number", "gene": "DMD", "copy_number": 1},
        {"id": "v3", "kind": "small", "gene": "X", "hgvs_c": "c.1A>G", "zygosity": "unknown"}]}
    pp = P.from_case(case, created="t", hpo_version="x", warnings=warnings)
    gis = pp["interpretations"][0]["diagnosis"]["genomicInterpretations"]
    vd = [g["variantInterpretation"]["variationDescriptor"] for g in gis]
    assert "vcfRecord" not in vd[0] and any("no genome assembly" in w for w in warnings)
    assert "allelicState" not in vd[0] and any("mosaic" in w for w in warnings)
    assert vd[1]["structuralType"]["id"] == "SO:0001019"  # one copy is normal for X in a male: not "loss"
    assert vd[2]["allelicState"]["id"] == "GENO:0000137"



def test_r_a_p0_1_labels_come_from_the_hpo_release_when_installed(tmp_path, monkeypatch):
    """With a release installed, a typed label is replaced by the HPO label (so it is exported, not
    refused) and a disease the HPOA knows gets its curated name."""
    import os
    import shutil

    from zebra import hpo_local

    fix = os.path.join(os.path.dirname(__file__), "fixtures", "hpo")
    d = tmp_path / "hpo"
    d.mkdir()
    for f in hpo_local.FILES:
        shutil.copy(os.path.join(fix, f), d / f)
    with open(d / "phenotype.hpoa", "a", encoding="utf-8") as fh:
        fh.write("OMIM:607208\tDravet syndrome (curated)\t\tHP:0001250\tPMID:1\tPCS\t\t\t\t\tP\tHPO:x\n")
    monkeypatch.setenv("ZEBRA_HPO_DIR", str(d))
    case = {"id": "c", "privacy": {"identifiers": ["张三"]},
            "phenotypes": [{"id": "HP:0001250", "label": "张三的抽搐"}],
            "hypotheses": [{"disease": "张三's Dravet", "status": "confirmed", "ids": {"OMIM": "607208"}}]}
    pp = P.from_case(case, created="t")
    assert pp["phenotypicFeatures"][0]["type"]["label"] == "Seizure"
    assert pp["interpretations"][0]["diagnosis"]["disease"] == {"id": "OMIM:607208", "label": "Dravet syndrome (curated)"}
    hp = next(r for r in pp["metaData"]["resources"] if r["namespacePrefix"] == "HP")
    assert hp["version"] == "fixture-1"


def test_r_a_p0_1_numbers_are_checked_as_well_as_text():
    case = {"id": "c", "privacy": {"identifiers": ["88812345"]},
            "variants": [{"id": "v1", "kind": "repeat_expansion", "gene": "FMR1", "repeat_count": 88812345}]}
    with pytest.raises(P.PhenopacketError, match="protected identifier"):
        P.from_case(case, created="t", hpo_version="x")
