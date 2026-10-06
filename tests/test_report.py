"""zebra report export: Markdown → HTML / Word / PDF for a family handout (CP1-13, 代操作)."""

from __future__ import annotations

import json
import os
import zipfile
import xml.dom.minidom

import pytest

from zebra import case as case_mod
from zebra import report_export as rx
from zebra.cli import main

LETTER = """# 给家属的信

这封信说明**看了什么**、*发现了什么* [E12]。

- 表型 [E3, E4]
  - 首次发作在 6 个月
- 变异 `NM_001165963.4:c.2134C>T`
1. 检索 ClinVar
2. 对照表型

| 项目 | 结果 |
|---|---|
| 变异 | SCN1A c.2134C>T \\| 杂合 |

> 钠通道阻滞剂可能加重发作 [PMID 31010818]。

见 [GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1318/) 与 [本地](file:///etc/passwd)。
"""


def _case_with_letter(tmp_path, text=LETTER, identifiers=("王小雨", "2019-03-02", "MZ0012345")):
    d = tmp_path / "case"
    case_mod.init(str(d), title="t")
    case_mod.add_identifiers(str(d), list(identifiers))
    (d / "reports").mkdir(exist_ok=True)
    md = d / "reports" / "family-letter.md"
    md.write_text(text, "utf-8")
    return d, md


def test_cp1_13_parse_keeps_bullets_and_numbers_apart():
    blocks = rx.parse(LETTER)
    lists = [b for b in blocks if b[0] == "list"]
    assert len(lists) == 2
    assert [o for _, o, _, _ in lists[0][1]] == [False, False, False]
    assert [o for _, o, _, _ in lists[1][1]] == [True, True]
    table = next(b for b in blocks if b[0] == "table")
    assert table[2][0] == ["变异", "SCN1A c.2134C>T | 杂合"]  # escaped pipe stays in the cell


def test_cp1_13_html_and_docx_are_well_formed_and_chinese(tmp_path):
    _d, md = _case_with_letter(tmp_path)
    got = rx.export(str(md), ["html", "docx"], identifiers=[], engine=None)
    files = {f["format"]: f["path"] for f in got["files"]}
    page = open(files["html"], encoding="utf-8").read()
    assert '<html lang="zh-CN">' in page and "<ol>" in page and "<table>" in page
    assert 'class="ev">[E12]' in page
    assert 'href="file:' not in page  # only http(s) links become links
    with zipfile.ZipFile(files["docx"]) as z:
        for name in z.namelist():
            xml.dom.minidom.parseString(z.read(name))  # every part is well-formed XML
        doc = z.read("word/document.xml").decode("utf-8")
        rels = z.read("word/_rels/document.xml.rels").decode("utf-8")
    assert "给家属的信" in doc and "1." in doc and "2." in doc and "3.\t" not in doc
    assert "https://www.ncbi.nlm.nih.gov/books/NBK1318/" in rels and "file:" not in rels
    # rPr children in schema order: rStyle before b, rFonts before b, color before shd
    assert "<w:b/><w:bCs/><w:rStyle" not in doc


def test_cp1_13_export_refuses_a_report_with_a_protected_identifier(tmp_path):
    for leak in ("王 小雨 的随访", "生日 2019年3月2日", "病历号 MZ-0012345", "身份证 330106201903021234"):
        d, md = _case_with_letter(tmp_path / leak[:2], LETTER + "\n" + leak + "\n")
        with pytest.raises(PermissionError):
            rx.export(str(md), ["docx"], identifiers=case_mod.load(str(d))["privacy"]["identifiers"])
        assert not (md.parent / "family-letter.docx").exists()
    # an internal copy may keep it, with a warning
    got = rx.export(str(md), ["html"], identifiers=["王小雨"], allow_identifiers=True)
    assert any("internal" in w for w in got["warnings"])


def test_cp1_13_research_words_are_not_identifiers(tmp_path):
    # a short pinyin name must not match inside words; numbers in HGVS are not record numbers
    hit, _ = rx.identifier_hits("SCN1A c.2134C>T; PMID 31010818; Li-Fraumeni", ["Li", "2134"])
    assert hit is None


def test_cp1_13_cli_finds_the_case_and_reports_a_missing_pdf_engine(tmp_path, capsys, monkeypatch):
    _d, md = _case_with_letter(tmp_path)
    monkeypatch.setattr(rx, "find_pdf_engine", lambda: None)
    assert main(["--json", "report", "export", str(md), "--to", "word,pdf"]) == 0
    out = json.loads(capsys.readouterr().out)
    assert [f["format"] for f in out["result"]["files"]] == ["docx"]
    assert out["result"]["skipped"][0]["format"] == "pdf"
    assert any("PDF not made" in w for w in out["warnings"])
    md.write_text(LETTER + "\n王小雨\n", "utf-8")
    assert main(["--json", "report", "export", str(md)]) == 2
    assert "protected identifier" in json.loads(capsys.readouterr().out)["error"]["message"]


@pytest.mark.skipif(os.environ.get("ZEBRA_TEST_PDF") != "1" or rx.find_pdf_engine() is None,
                    reason="set ZEBRA_TEST_PDF=1 on a machine with Chrome/Edge/LibreOffice")
def test_cp1_13_pdf_is_printed_by_a_local_browser(tmp_path):
    _d, md = _case_with_letter(tmp_path)
    got = rx.export(str(md), ["pdf"], identifiers=[])
    pdf = next(f for f in got["files"] if f["format"] == "pdf")
    assert open(pdf["path"], "rb").read(5) == b"%PDF-"
    assert not (md.parent / "family-letter.html").exists()  # the HTML was only a print source


def test_cp1_13_a_record_number_never_matches_inside_a_scientific_id(tmp_path):
    text = "| HP:0012345 | Seizure |\nPMID 30012345, rs1060500100, NM_0012345.1, chr2:166012345\n"
    assert rx.identifier_hits(text, ["MZ0012345", "30012345"])[0] is None
    assert rx.identifier_hits(text + "病历号 MZ-001-2345\n", ["MZ0012345"])[0] is not None
