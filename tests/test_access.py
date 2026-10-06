"""`zebra access` and the bundled China access tables (CP0-4): approvals, 医保, China trials, 协作网 hospitals."""

from __future__ import annotations

import json

import pytest

from zebra import cli
from zebra.commands import access, china
from zebra.core import Outcome
from zebra.sources import china_access as ca
from zebra.sources import ctgov, opentargets as ot, regulators

APPROVALS = {
    "schema": 1,
    "title": "test copy",
    "provenance": {"coverage": "TEST: two rows from one document; not a complete register",
                   "sources": {"vat2019": {"title": "第一批罕见病药品清单", "document_no": "财税〔2019〕24号",
                                           "url": "https://www.gov.cn/x", "retrieved_at": "2026-10-06T00:00:00Z",
                                           "sha256": "a" * 64}}},
    "drugs": [
        {"drug_zh": "诺西那生钠注射液", "inn": "nusinersen", "inn_source": "vat2019", "brand": None,
         "indication_zh": "用于治疗5q脊髓性肌萎缩症", "indication_en": None,
         "list_diseases": [{"list": 1, "no": 110, "how": "indication names 脊髓性肌萎缩症"}],
         "approval": {"date": "2019-02-21", "year": 2019, "type": "进口", "source_id": "vat2019"},
         "evidence": [{"source_id": "vat2019", "what": "listed in 第一批罕见病药品清单"}]},
        {"drug_zh": None, "inn": "Risdiplam", "inn_source": "vat2019", "brand": None,
         "indication_zh": "脊髓性肌萎缩症", "indication_en": None, "list_diseases": [], "approval": None,
         "evidence": [{"source_id": "vat2019", "what": "临床急需境外新药名单 — 待申报"}]},
        {"drug_zh": "依达拉奉注射液", "inn": None, "inn_source": None, "brand": None,
         "indication_zh": "肌萎缩侧索硬化", "indication_en": None, "list_diseases": [], "approval": None,
         "evidence": [{"source_id": "vat2019", "what": "listed"}]},
    ],
}


@pytest.fixture
def bundled(tmp_path, monkeypatch):
    p = tmp_path / "approvals.json"
    p.write_text(json.dumps(APPROVALS, ensure_ascii=False), "utf-8")
    monkeypatch.setattr(ca, "APPROVALS_FILE", str(p))
    monkeypatch.setattr(ca, "_cache", {})
    return p


@pytest.fixture
def offline(monkeypatch):
    """No network: Open Targets, the agency checks and ClinicalTrials.gov answer from here."""
    calls = {"ctgov": [], "check": []}

    def search(text, entities=("disease", "target"), limit=10):
        hits = []
        if text.lower().startswith("spinal muscular atrophy"):
            hits = [{"id": "MONDO_0001516", "name": "spinal muscular atrophy", "entity": "disease"}]
        if text.lower() == "nusinersen":
            hits = [{"id": "CHEMBL3833342", "name": "NUSINERSEN", "entity": "drug"}]
        return Outcome({"query": text, "total": len(hits), "hits": hits}, sources=[{"db": "Open Targets search"}])

    def disease(efo_id, drug_limit=40, n_targets=10):
        rows = [{"drug": "NUSINERSEN", "stage": "APPROVAL"}, {"drug": "RISDIPLAM", "stage": "APPROVAL"},
                {"drug": "SALANERSEN", "stage": "PHASE_3"}]
        return Outcome({"id": efo_id, "name": "spinal muscular atrophy", "synonyms": ["SMA", "spinal muscular atrophies"],
                        "drugs": {"rows": rows}}, sources=[{"db": "Open Targets"}])

    def check(drug, terms, with_fda_date=True):
        calls["check"].append((drug, list(terms)))
        return Outcome({"drug": drug, "FDA": {"status": "approved_for_disease", "reading": "label", "approved_on": "2016-12-23"},
                        "EMA": {"status": "approved_for_disease", "reading": "EU", "date": "2017-05-30"},
                        "approved_in": ["FDA (United States)", "EMA (European Union)"], "refused_or_withdrawn_in": [],
                        "headline": "approved"}, sources=[{"db": "openFDA drug label", "record": drug, "url": "u"}])

    def trials(condition, term=None, country=None, status="RECRUITING", limit=20, **kw):
        calls["ctgov"].append((condition, country, status, limit))
        return Outcome({"condition": condition, "country": "China", "status": status, "total": 1, "returned": 1,
                        "studies": [{"nct_id": "NCT05187260", "title": "ASO for SMA", "phases": ["PHASE1"],
                                     "status": "RECRUITING", "sites_in_country": 1,
                                     "sites": [{"facility": "Peking Union Medical College Hospital", "city": "Beijing"}],
                                     "url": "https://clinicaltrials.gov/study/NCT05187260"}],
                        "filtered": [{"nct_id": "NCT07190300"}]},
                       sources=[{"db": "ClinicalTrials.gov"}], warnings=["1 trial(s) set aside"])

    def by_ingredient(base):
        return Outcome({"ingredient": base, "applications": [
            {"application": "NDA209531", "brands": ["SPINRAZA"], "approved_on": "2016-12-23", "url": "u"}]},
            sources=[{"db": "openFDA Drugs@FDA", "record": base}])

    def ema(drug, terms=()):
        return Outcome({"records": [{"medicine": "Spinraza", "status": "Authorised", "authorised_on": "2017-05-30",
                                     "url": "https://www.ema.europa.eu/en/medicines/human/EPAR/spinraza"}]},
                       sources=[{"db": "EMA medicines (EPAR register export)"}])

    monkeypatch.setattr(ot, "search", search)
    monkeypatch.setattr(ot, "disease", disease)
    monkeypatch.setattr(regulators, "check", check)
    monkeypatch.setattr(regulators, "fda_by_ingredient", by_ingredient)
    monkeypatch.setattr(regulators, "ema", ema)
    monkeypatch.setattr(ctgov, "search", trials)
    return calls


def _access(argv, capsys):
    code = cli.main(["--json", "access", *argv])
    return code, json.loads(capsys.readouterr().out)


def test_cp0_4_a_chinese_disease_name_answers_all_four_questions(bundled, offline, capsys):
    code, env = _access(["脊髓性肌萎缩症", "--province", "浙江"], capsys)
    assert code == 0, env
    r = env["result"]
    assert r["kind"] == "disease" and r["national_rare_disease_list"]["status"] == "on_list"
    # approvals: FDA/EMA from the agency check, China from the bundled documents
    fe = {d["drug"]: d for d in r["approvals"]["fda_ema"]["drugs"]}
    assert fe["NUSINERSEN"]["FDA"]["status"] == "approved_for_disease"
    assert fe["NUSINERSEN"]["china_bundled"].startswith("approved in China per the bundled documents")
    # named in the documents but not approved there: never "in the list" as if approved
    assert fe["RISDIPLAM"]["china_bundled"].startswith("named in the bundled documents")
    assert "status unknown here" in fe["RISDIPLAM"]["china_bundled"]
    cn = r["approvals"]["china_nmpa"]
    assert cn["status"] == "approved_in_china" and cn["rows"][0]["drug_zh"] == "诺西那生钠注射液"
    assert any(r["approved_in_china"] is False for r in cn["rows"])
    assert cn["rows"][0]["approved_in_china"] is True and cn["rows"][0]["approval"]["date"] == "2019-02-21"
    assert "TEST" in cn["coverage"]
    # 医保 from the real bundled NRDL: the restriction names the disease
    nr = r["nrdl_医保"]
    assert nr["status"] == "listed" and "医保发〔2025〕33号" in nr["edition"]
    names = {x["name_zh"]: x for x in nr["rows"]}
    assert names["诺西那生钠注射液"]["restriction"] == "限5q脊髓性肌萎缩症。"
    assert "利司扑兰口服溶液用散" in names
    # trials in China, with the ChiCTR note
    tr = r["trials_china"]
    assert tr["status"] == "ok" and tr["trials"][0]["nct_id"] == "NCT05187260" and "HTTP 405" in tr["chictr"]
    assert offline["ctgov"][0][1] == "China"
    # 协作网 hospitals of the province
    ho = r["hospitals_协作网"]
    assert ho["status"] == "ok" and ho["count"] == 15 and ho["province"] == "浙江"
    # every block cites a source
    dbs = {s["db"] for s in env["sources"]}
    assert {"国家医保药品目录 (NHSA)", "China drug approvals (NMPA/CDE/MOF documents)", "全国罕见病诊疗协作网 (NHC)",
            "ClinicalTrials.gov"} <= dbs
    text = access._render(r)
    assert "not approved" not in text.lower().replace("non-approval", "")


def test_cp0_4_missing_china_tables_are_not_bundled_never_not_approved(offline, monkeypatch, capsys):
    monkeypatch.setattr(ca, "APPROVALS_FILE", "/nonexistent/a.json")
    monkeypatch.setattr(ca, "NRDL_FILE", "/nonexistent/n.json")
    monkeypatch.setattr(ca, "_cache", {})
    code, env = _access(["脊髓性肌萎缩症"], capsys)
    r = env["result"]
    assert r["approvals"]["china_nmpa"]["status"] == "not_bundled"
    assert r["nrdl_医保"]["status"] == "not_bundled"
    assert any("not evidence of anything about approval" in w for w in env["warnings"])
    assert any("医保 status not checked" in w for w in env["warnings"])


def test_cp0_4_a_chinese_drug_name_is_a_drug_with_its_医保_restriction(bundled, offline, capsys):
    code, env = _access(["诺西那生钠"], capsys)
    r = env["result"]
    assert r["kind"] == "drug" and r["inn"] == "nusinersen"
    assert r["approvals"]["china_nmpa"]["status"] == "approved_in_china"
    rows = r["nrdl_医保"]["rows"]
    assert rows and rows[0]["name_zh"] == "诺西那生钠注射液" and rows[0]["section"].startswith("协议期内谈判药品")
    assert r["approvals"]["FDA"]["status"] == "approved" and r["approvals"]["EMA"]["status"] == "Authorised"


def test_cp0_4_an_english_drug_with_no_chinese_name_says_医保_was_not_checked(offline, monkeypatch, capsys):
    monkeypatch.setattr(ca, "APPROVALS_FILE", "/nonexistent/a.json")
    monkeypatch.setattr(ca, "_cache", {})
    monkeypatch.setattr(regulators, "fda_by_ingredient", lambda b: Outcome({"ingredient": b, "applications": []}))
    monkeypatch.setattr(regulators, "ema", lambda d, terms=(): Outcome({"records": []}))
    code, env = _access(["nusinersen"], capsys)
    r = env["result"]
    assert r["kind"] == "drug"
    assert r["nrdl_医保"]["status"] == "not_checked" and "Chinese generic name" in r["nrdl_医保"]["note"]
    assert r["approvals"]["china_nmpa"]["status"] == "not_bundled"
    assert r["approvals"]["FDA"]["status"] == "no_record" and "not proof it is unapproved" in r["approvals"]["FDA"]["note"]


def test_cp0_4_a_drug_not_in_the_bundled_list_is_worded_as_such(bundled, offline, monkeypatch, capsys):
    monkeypatch.setattr(regulators, "fda_by_ingredient", lambda b: Outcome({"ingredient": b, "applications": []}))
    monkeypatch.setattr(regulators, "ema", lambda d, terms=(): Outcome({"records": []}))
    code, env = _access(["--as", "drug", "givinostat"], capsys)
    cn = env["result"]["approvals"]["china_nmpa"]
    assert cn["status"] == "not_in_bundled_list"
    assert "not in the bundled list" in cn["note"] and "not evidence of anything about approval" in cn["note"]


def test_cp0_4_a_qualified_disease_is_not_called_on_the_list(bundled, offline, capsys):
    code, env = _access(["地中海贫血", "--trials", "0"], capsys)
    r = env["result"]
    assert r["national_rare_disease_list"]["status"] == "qualified"
    assert r["trials_china"]["status"] in ("skipped", "not_checked")


def test_cp0_4_an_unresolved_name_is_said_so(bundled, offline, capsys):
    code, env = _access(["某某病"], capsys)
    assert env["result"]["kind"] == "unresolved" and "not answered" in env["result"]["note"]


def test_cp0_4_bundled_nrdl_is_the_current_edition_with_provenance():
    d = ca.load_nrdl()
    p = d["provenance"]
    assert p["edition"].startswith("2025") and p["document_no"] == "医保发〔2025〕33号" and p["effective"] == "2026-01-01"
    s = p["sources"]["nrdl2025"]
    assert s["url"].startswith("https://") and len(s["sha256"]) == 64 and s["retrieved_at"]
    assert p["counts"]["west_part_numbers"] == 1446 and p["counts"]["negotiated_and_bidding_west"] == 411
    assert p["builder"].startswith("tools/china_access/build_nrdl.py")
    assert "医保基金不予支付" in d["commercial_innovative_list_label"]


@pytest.mark.parametrize("a,b,same", [("诺西那生钠", "诺西那生钠注射液", True), ("诺西那生", "诺西那生钠注射液", True),
                                      ("利司扑兰", "利司扑兰口服溶液用散", True), ("阿加糖酶α", "阿加糖酶β注射用浓溶液", False),
                                      ("沙丙蝶呤", "盐酸沙丙蝶呤片", True), ("维拉苷酶β", "注射用维拉苷酶β", True),
                                      ("钠", "诺西那生钠注射液", False)])
def test_cp0_4_chinese_drug_names_match_on_the_moiety(a, b, same):
    assert ca.same_drug_zh(a, b) is same


def test_cp0_4_nrdl_by_disease_finds_the_restriction_and_labels_the_commercial_list():
    rows = ca.nrdl_for_disease(["戈谢病"])
    lists = {r["list"] for r in rows}
    assert "CIDL" in lists  # 维拉苷酶β is on the commercial list only
    nrdl = ca.load_nrdl()
    c = ca.compact_nrdl(next(r for r in rows if r["list"] == "CIDL"), nrdl)
    assert "医保 does NOT pay" in c["list"]


@pytest.mark.live
def test_live_access_dravet(capsys):
    code = cli.main(["--json", "access", "Dravet syndrome"])
    env = json.loads(capsys.readouterr().out)
    r = env["result"]
    assert code == 0 and r["kind"] == "disease"
    fe = {d["drug"]: d for d in r["approvals"]["fda_ema"]["drugs"]}
    assert fe["CANNABIDIOL"]["FDA"]["status"] == "approved_for_disease"
    assert r["national_rare_disease_list"]["status"] == "on_list"
    assert r["hospitals_协作网"]["count"] == 33


def test_cp0_4_a_near_miss_is_never_answered_for(bundled, offline, monkeypatch, capsys):
    """`access x` used to answer for Open Targets' first disease hit (Waldenstrom macroglobulinemia)."""
    def search(text, entities=("disease", "target"), limit=10):
        return Outcome({"hits": [{"id": "MONDO_0100280", "name": "Waldenstrom macroglobulinemia", "entity": "disease"}]})

    monkeypatch.setattr(ot, "search", search)
    code, env = _access(["x"], capsys)
    r = env["result"]
    assert r["kind"] == "unresolved" and r["search_hits"][0]["id"] == "MONDO_0100280"
    assert "approvals" not in r and "nrdl_医保" not in r
    code, env = _access(["--as", "disease", "x"], capsys)
    assert env["result"]["kind"] == "unresolved"


def test_cp0_4_a_disease_id_is_answered_with_its_open_targets_name(bundled, offline, capsys):
    code, env = _access(["MONDO_0001516", "--trials", "0"], capsys)
    r = env["result"]
    assert r["kind"] == "disease" and r["name"] == "spinal muscular atrophy"
    assert r["national_rare_disease_list"]["status"] == "on_list"
    assert r["approvals"]["china_nmpa"]["status"] == "approved_in_china"


def test_cp0_4_negative_trials_is_a_usage_error(capsys):
    code, env = _access(["x", "--trials", "-1"], capsys)
    assert code == 2 and env["error"]["type"] == "UsageError"


# ---------------------------------------------------------------- adversarial review (W4 round), access

@pytest.mark.parametrize("a,b", [("胰岛素", "胰岛素瘤"), ("维拉苷酶", "注射用维拉苷酶β"), ("葡萄糖", "葡萄糖酸钙注射液"),
                                 ("西妥昔单抗", "西妥昔单抗β"), ("甲状腺素", "甲状腺片")])
def test_rev_p0_4_a_different_product_is_not_the_same_drug(a, b):
    assert ca.same_drug_zh(a, b) is False


@pytest.mark.parametrize("a,b", [("司美替尼", "硫酸氢司美替尼胶囊"), ("尼达尼布", "乙磺酸尼达尼布软胶囊"),
                                 ("诺西那生", "诺西那生钠注射液")])
def test_rev_p1_salt_and_form_words_are_stripped(a, b):
    assert ca.same_drug_zh(a, b) is True


def test_rev_p0_4_a_disease_shaped_chinese_name_is_not_a_drug(bundled, offline, capsys):
    code, env = _access(["胰岛素瘤", "--trials", "0"], capsys)
    assert env["result"]["kind"] != "drug"


def test_rev_p0_5_the_users_inn_is_used_for_fda_and_ema(offline, monkeypatch, tmp_path, capsys):
    data = json.loads(json.dumps(APPROVALS))
    data["drugs"].append({"drug_zh": None, "inn": "Galafold（Migalastat hydrochloride）", "list_diseases": [],
                          "approval": None, "evidence": [{"source_id": "vat2019", "what": "x"}]})
    p = tmp_path / "a.json"
    p.write_text(json.dumps(data, ensure_ascii=False), "utf-8")
    monkeypatch.setattr(ca, "APPROVALS_FILE", str(p))
    monkeypatch.setattr(ca, "_cache", {})
    asked = []
    monkeypatch.setattr(regulators, "fda_by_ingredient", lambda b: asked.append(b) or Outcome({"applications": []}))
    monkeypatch.setattr(regulators, "ema", lambda d, terms=(): asked.append(d) or Outcome({"records": []}))
    code, env = _access(["migalastat"], capsys)
    assert asked == ["migalastat", "migalastat"]
    assert env["result"]["approvals"]["china_nmpa"]["status"] == "named_not_approved"
    assert ca.clean_inn("Galafold（Migalastat hydrochloride）") == "Migalastat hydrochloride"


def test_rev_p0_6_a_secondary_condition_is_not_the_disease():
    assert not ca.mentions("限中危或高危原发性骨髓纤维化(PMF)、真性红细胞增多症继发性骨髓纤维化(PPV-MF)", "真性红细胞增多症")
    assert ca.mentions("单药适用于既往接受羟基脲治疗效果不佳的真性红细胞增多症成人患者。", "真性红细胞增多症")
    rows = ca.nrdl_for_disease(["真性红细胞增多症"])
    assert not [r for r in rows if "PPV-MF" in (r.get("restriction") or "") and "真性红细胞增多症患者" not in
                (r.get("restriction") or "")]


def test_rev_p1_nmpa_status_says_approved_named_or_absent(bundled, offline, monkeypatch, capsys):
    monkeypatch.setattr(regulators, "fda_by_ingredient", lambda b: Outcome({"applications": []}))
    monkeypatch.setattr(regulators, "ema", lambda d, terms=(): Outcome({"records": []}))
    code, env = _access(["--as", "drug", "risdiplam"], capsys)
    cn = env["result"]["approvals"]["china_nmpa"]
    assert cn["status"] == "named_not_approved" and "status unknown here" in cn["status_meaning"]
