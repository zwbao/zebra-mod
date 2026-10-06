"""Turn a Markdown report into HTML, Word (.docx) and PDF that a family can open.

Standard library only. The Markdown is the subset the report skills write:
headings, paragraphs, bullet and numbered lists (nested by indentation),
pipe tables, block quotes, code blocks, rules, and inline bold, italic, code,
links and evidence tags ([E12], [PMID 12345678]).

HTML and Word are written here. PDF is printed from the HTML by a browser
already on the machine (Chrome, Edge, Chromium, Brave) in headless mode with a
throw-away profile, or converted from the Word file by LibreOffice; with
neither, PDF is skipped and the reason says so.

Before anything is written the report is checked for the case's protected
identifiers and for Chinese resident ID numbers: a handout leaves the machine
the moment it is printed or sent, and it never needs them.
"""

from __future__ import annotations

import html
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unicodedata
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

MAX_BYTES = 5 * 1024 * 1024
PDF_TIMEOUT = 90.0

# ------------------------------------------------------------------ markdown

Block = Tuple[Any, ...]

_FENCE = re.compile(r"^\s*(```|~~~)")
_HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")
_RULE = re.compile(r"^\s*([-*_])(\s*\1){2,}\s*$")
_ITEM = re.compile(r"^(\s*)([-*+]|\d{1,3}[.)])\s+(.*)$")
_TABLE_SEP = re.compile(r"^\s*\|?\s*:?-{2,}:?\s*(\|\s*:?-{2,}:?\s*)*\|?\s*$")
_CJK = re.compile(r"[\u3000-\u303f\u3400-\u9fff\uf900-\ufaff\uff00-\uffef]")


def _join(a: str, b: str) -> str:
    """Join two wrapped lines: no space between two CJK characters, one otherwise."""
    if not a:
        return b
    if _CJK.match(a[-1:]) and _CJK.match(b[:1]):
        return a + b
    return a + " " + b


def _cells(line: str) -> List[str]:
    s = line.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|") and not s.endswith("\\|"):
        s = s[:-1]
    out, cur, i = [], "", 0
    while i < len(s):
        if s[i] == "\\" and i + 1 < len(s) and s[i + 1] == "|":
            cur += "|"
            i += 2
            continue
        if s[i] == "|":
            out.append(cur.strip())
            cur = ""
        else:
            cur += s[i]
        i += 1
    out.append(cur.strip())
    return out


def parse(text: str) -> List[Block]:
    """Blocks: ("h", level, text) ("p", text) ("list", [(depth, ordered, marker, text)])
    ("table", header, rows) ("quote", text) ("code", text) ("hr",)."""
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    blocks: List[Block] = []
    para = ""
    i = 0

    def flush() -> None:
        nonlocal para
        if para.strip():
            blocks.append(("p", para.strip()))
        para = ""

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if stripped.startswith("<!--"):
            flush()
            while i < len(lines) and "-->" not in lines[i]:
                i += 1
            i += 1
            continue
        if _FENCE.match(line):
            flush()
            fence = _FENCE.match(line).group(1)  # type: ignore[union-attr]
            body: List[str] = []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith(fence):
                body.append(lines[i])
                i += 1
            blocks.append(("code", "\n".join(body)))
            i += 1
            continue
        if not stripped:
            flush()
            i += 1
            continue
        m = _HEADING.match(stripped)
        if m:
            flush()
            blocks.append(("h", len(m.group(1)), m.group(2)))
            i += 1
            continue
        if _RULE.match(line):
            flush()
            blocks.append(("hr",))
            i += 1
            continue
        if "|" in stripped and i + 1 < len(lines) and _TABLE_SEP.match(lines[i + 1]) and "-" in lines[i + 1]:
            flush()
            header = _cells(line)
            rows: List[List[str]] = []
            i += 2
            while i < len(lines) and "|" in lines[i] and lines[i].strip():
                row = _cells(lines[i])
                rows.append((row + [""] * len(header))[: max(len(header), 1)])
                i += 1
            blocks.append(("table", header, rows))
            continue
        if stripped.startswith(">"):
            flush()
            quote = ""
            while i < len(lines) and lines[i].strip().startswith(">"):
                quote = _join(quote, lines[i].strip()[1:].strip()) if lines[i].strip()[1:].strip() else quote + "\n"
                i += 1
            blocks.append(("quote", quote.strip()))
            continue
        m = _ITEM.match(line)
        if m:
            flush()
            items: List[Tuple[int, bool, str, str]] = []
            base = len(m.group(1).expandtabs(4))
            while i < len(lines):
                m = _ITEM.match(lines[i])
                if m:
                    depth = max(0, (len(m.group(1).expandtabs(4)) - base) // 2)
                    marker = m.group(2)
                    ordered = marker[0].isdigit()
                    # a top-level switch between bullets and numbers starts a new list
                    if depth == 0 and items and next((o for d, o, _, _ in items if d == 0), ordered) != ordered:
                        break
                    items.append((min(depth, 4), ordered, marker, m.group(3).strip()))
                    i += 1
                    continue
                nxt = lines[i]
                if nxt.strip() and (nxt.startswith("  ") or nxt.startswith("\t")) and items:
                    d, o, mk, t = items[-1]
                    items[-1] = (d, o, mk, _join(t, nxt.strip()))
                    i += 1
                    continue
                break
            blocks.append(("list", items))
            continue
        para = _join(para, stripped)
        i += 1
    flush()
    return blocks


# inline runs: (text, flags) with flags from b, i, code, ev, and ("href", url)
Run = Tuple[str, Dict[str, Any]]
_INLINE = re.compile(
    r"(?P<code>`[^`]+`)"
    r"|(?P<bold>\*\*(?P<bt>.+?)\*\*|__(?P<bt2>.+?)__)"
    r"|(?P<ital>(?<![\w*])\*(?!\s)(?P<it>[^*]+?)(?<!\s)\*(?!\w)|(?<![\w_])_(?!\s)(?P<it2>[^_]+?)(?<!\s)_(?![\w]))"
    r"|(?P<link>\[(?P<lt>[^\]]+)\]\((?P<lu>[^)\s]+)\))"
    r"|(?P<auto><(?P<au>https?://[^>\s]+)>)"
    r"|(?P<ev>\[(?:E\d+(?:\s*[,，、;；]\s*E\d+)*|PMID[:：]?\s*\d+(?:\s*[,，、;；]\s*(?:PMID[:：]?\s*)?\d+)*)\])"
)


def inline(text: str, flags: Optional[Dict[str, Any]] = None) -> List[Run]:
    flags = dict(flags or {})
    runs: List[Run] = []
    pos = 0
    for m in _INLINE.finditer(text):
        if m.start() > pos:
            runs.append((text[pos:m.start()], dict(flags)))
        if m.group("code"):
            runs.append((m.group("code")[1:-1], {**flags, "code": True}))
        elif m.group("bold"):
            runs.extend(inline(m.group("bt") or m.group("bt2") or "", {**flags, "b": True}))
        elif m.group("ital"):
            runs.extend(inline(m.group("it") or m.group("it2") or "", {**flags, "i": True}))
        elif m.group("link"):
            url = m.group("lu")
            safe = url if re.match(r"^https?://", url, re.I) else None
            runs.extend(inline(m.group("lt"), {**flags, **({"href": safe} if safe else {})}))
            if not safe:
                runs.append((f" ({url})", dict(flags)))
        elif m.group("auto"):
            runs.append((m.group("au"), {**flags, "href": m.group("au")}))
        elif m.group("ev"):
            runs.append((m.group("ev"), {**flags, "ev": True}))
        pos = m.end()
    if pos < len(text):
        runs.append((text[pos:], dict(flags)))
    return [r for r in runs if r[0]]


def title_of(blocks: Sequence[Block], fallback: str) -> str:
    for b in blocks:
        if b[0] == "h" and b[1] == 1:
            return "".join(t for t, _ in inline(b[2]))
    return fallback


# ------------------------------------------------------------------ privacy

_RESIDENT_ID = re.compile(r"(?<!\d)\d{6}(?:18|19|20)\d{2}(?:0[1-9]|1[0-2])(?:0[1-9]|[12]\d|3[01])\d{3}[\dXx](?!\d)")
_MOBILE = re.compile(r"(?<!\d)1[3-9]\d{9}(?!\d)")
_SCIENTIFIC_IDS = re.compile("|".join([
    r"\b(?:hp|omim|mim|orpha|orphanet|mondo|doid|ncit|efo|uberon|go|hgnc|cl|chebi|mp)\s*[:_]\s*\d+",
    r"\brs\d+", r"\b(?:nm|nr|np|nc|ng|nt|xm|xp|enst|ensg|ensp|ccds|lrg)_?\d+(?:\.\d+)?",
    r"\bpmid\s*:?\s*\d+", r"\bpmc\d+", r"\bnct\d{8}\b", r"\bchictr-?\w+",
    r"\b(?:chr)?(?:\d{1,2}|x|y|mt?)\s*[:-]\s*\d+(?:\s*[-_]\s*\d+)?", r"\b[cgpnmr]\.\S+", r"\b10\.\d{4,9}/\S+",
]), re.I)
_DATE = re.compile(r"^\s*(\d{4})\s*[-/.年]\s*(\d{1,2})\s*[-/.月]\s*(\d{1,2})\s*日?\s*$")


def _fold(s: str) -> str:
    return unicodedata.normalize("NFKC", s).casefold()


def _date_forms(value: str) -> List[str]:
    m = _DATE.match(value)
    if not m:
        return []
    y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
    return [f"{y}-{mo:02d}-{d:02d}", f"{y}-{mo}-{d}", f"{y}/{mo:02d}/{d:02d}", f"{y}/{mo}/{d}", f"{y}.{mo}.{d}",
            f"{y}.{mo:02d}.{d:02d}", f"{y}年{mo}月{d}日", f"{y}年{mo:02d}月{d:02d}日", f"{y}{mo:02d}{d:02d}",
            f"{d:02d}/{mo:02d}/{y}", f"{mo:02d}/{d:02d}/{y}"]


def identifier_hits(text: str, identifiers: Sequence[str]) -> Tuple[Optional[str], List[str]]:
    """(why the text must not be exported, or None; warnings). Checked line by line."""
    warnings: List[str] = []
    lines = [_fold(line) for line in text.split("\n")]
    tight_lines = [re.sub(r"\s+", "", line) for line in lines]
    # a record number hides in digits, but never in a scientific identifier's (HP:0012345, rs…, PMID …)
    digit_lines = [re.sub(r"\D", "", _SCIENTIFIC_IDS.sub(" ", line)) for line in lines]
    for n, raw in enumerate(identifiers, 1):
        value = _fold(str(raw or "").strip())
        if len(value) < 2:
            continue
        label = f"protected identifier #{n} of the case"
        forms = [re.sub(r"\s+", "", f) for f in _date_forms(value)]
        if forms and any(f in t for t in tight_lines for f in forms):
            return f"{label} (a date)", warnings
        digits = re.sub(r"\D", "", value)
        if len(digits) >= 6 and len(digits) >= len(re.sub(r"[\s\-_/]", "", value)) - 3:
            if any(digits in d for d in digit_lines):
                return f"{label} (a number)", warnings
            continue
        if _CJK.search(value):
            if any(re.sub(r"\s+", "", value) in t for t in tight_lines):
                return label, warnings
            continue
        tokens = [t for t in re.split(r"[\s,;]+", value) if len(t) >= 2]
        if not tokens or (len(tokens) == 1 and len(tokens[0]) < 3):
            continue
        pats = [re.compile(r"(?<![\w])" + re.escape(t) + r"(?![\w])") for t in tokens]
        if any(all(p.search(line) for p in pats) for line in lines):
            return label, warnings
    for line in lines:
        if _RESIDENT_ID.search(line):
            return "a Chinese resident ID number", warnings
    if any(_MOBILE.search(line) for line in lines):
        warnings.append("the report contains what looks like a mobile phone number: keep it only if it is an "
                        "organisation's public number, never the family's")
    return None, warnings


# ------------------------------------------------------------------ HTML

_CSS = """
:root { color-scheme: light; }
@page { size: A4; margin: 18mm 16mm; }
body { margin: 0 auto; max-width: 46rem; padding: 2rem 1.25rem; color: #1d2329; background: #fff;
  font: 15px/1.75 "PingFang SC", "Hiragino Sans GB", "Noto Sans CJK SC", "Source Han Sans SC", "Microsoft YaHei",
  -apple-system, "Segoe UI", Roboto, sans-serif; }
h1 { font-size: 1.6rem; line-height: 1.35; margin: 0 0 .4rem; }
h2 { font-size: 1.22rem; margin: 1.8rem 0 .5rem; padding-bottom: .25rem; border-bottom: 1px solid #d9dee3; }
h3 { font-size: 1.05rem; margin: 1.3rem 0 .4rem; }
h4, h5, h6 { font-size: 1rem; margin: 1rem 0 .3rem; }
p { margin: .55rem 0; }
ul, ol { margin: .4rem 0 .6rem; padding-left: 1.5rem; }
li { margin: .2rem 0; }
.table { overflow-x: auto; margin: .7rem 0; }
table { border-collapse: collapse; width: 100%; font-size: .92rem; }
th, td { border: 1px solid #c9d0d6; padding: .35rem .55rem; text-align: left; vertical-align: top; }
th { background: #f1f4f6; font-weight: 600; }
blockquote { margin: .8rem 0; padding: .5rem .9rem; border-left: 4px solid #8aa4b8; background: #f5f8fa; }
code { font: .88em/1.4 ui-monospace, "SF Mono", Menlo, Consolas, monospace; background: #f1f3f5; padding: .05em .3em;
  border-radius: 3px; }
pre { background: #f1f3f5; padding: .7rem .9rem; overflow-x: auto; border-radius: 4px; }
pre code { background: none; padding: 0; }
hr { border: 0; border-top: 1px solid #d9dee3; margin: 1.4rem 0; }
a { color: #1f5f8b; }
.ev { color: #6b7780; font-size: .8em; white-space: nowrap; }
.meta { color: #6b7780; font-size: .85rem; margin-bottom: 1.2rem; }
@media print { body { max-width: none; padding: 0; } a { color: inherit; text-decoration: none; }
  a.ext::after { content: " (" attr(href) ")"; font-size: .85em; color: #4a5560; word-break: break-all; }
  h2, h3 { break-after: avoid; } tr, blockquote { break-inside: avoid; } }
"""


def _html_runs(runs: Sequence[Run]) -> str:
    out = []
    for text, f in runs:
        s = html.escape(text)
        if f.get("code"):
            s = f"<code>{s}</code>"
        if f.get("ev"):
            s = f'<span class="ev">{s}</span>'
        if f.get("i"):
            s = f"<em>{s}</em>"
        if f.get("b"):
            s = f"<strong>{s}</strong>"
        if f.get("href"):
            # on paper a link is only useful with its address: print it after the text
            ext = ' class="ext"' if text.strip() != f["href"] else ""
            s = f'<a href="{html.escape(f["href"], quote=True)}"{ext}>{s}</a>'
        out.append(s)
    return "".join(out)


def _html_list(items: Sequence[Tuple[int, bool, str, str]]) -> str:
    out: List[str] = []
    stack: List[str] = []
    for depth, ordered, _marker, text in items:
        tag = "ol" if ordered else "ul"
        while len(stack) > depth + 1:
            out.append(f"</li></{stack.pop()}>")
        if len(stack) == depth + 1 and stack[-1] != tag:
            out.append(f"</li></{stack.pop()}>")  # same level, the other kind of list
        elif len(stack) == depth + 1:
            out.append("</li>")
        while len(stack) < depth + 1:
            stack.append(tag)
            out.append(f"<{tag}>")
        out.append(f"<li>{_html_runs(inline(text))}")
    while stack:
        out.append(f"</li></{stack.pop()}>")
    return "".join(out)


def to_html(blocks: Sequence[Block], title: str, lang: str, generated: str) -> str:
    body: List[str] = []
    for b in blocks:
        kind = b[0]
        if kind == "h":
            body.append(f"<h{b[1]}>{_html_runs(inline(b[2]))}</h{b[1]}>")
        elif kind == "p":
            body.append(f"<p>{_html_runs(inline(b[1]))}</p>")
        elif kind == "list":
            body.append(_html_list(b[1]))
        elif kind == "table":
            head = "".join(f"<th>{_html_runs(inline(c))}</th>" for c in b[1])
            rows = "".join("<tr>" + "".join(f"<td>{_html_runs(inline(c))}</td>" for c in r) + "</tr>" for r in b[2])
            body.append(f'<div class="table"><table><thead><tr>{head}</tr></thead><tbody>{rows}</tbody></table></div>')
        elif kind == "quote":
            body.append("<blockquote>" + "".join(f"<p>{_html_runs(inline(p))}</p>" for p in b[1].split("\n") if p.strip())
                        + "</blockquote>")
        elif kind == "code":
            body.append(f"<pre><code>{html.escape(b[1])}</code></pre>")
        elif kind == "hr":
            body.append("<hr>")
    meta = ("生成于 " if lang.startswith("zh") else "Generated ") + generated + " · zebra-mod"
    return (f'<!doctype html>\n<html lang="{lang}">\n<head>\n<meta charset="utf-8">\n'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">\n'
            f"<title>{html.escape(title)}</title>\n<style>{_CSS}</style>\n</head>\n<body>\n"
            + "\n".join(body) + f'\n<p class="meta">{html.escape(meta)}</p>\n</body>\n</html>\n')


# ------------------------------------------------------------------ Word

_W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
_R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_XML_BAD = re.compile("[\x00-\x08\x0b\x0c\x0e-\x1f\ufffe\uffff]")


def _x(s: str) -> str:
    return html.escape(_XML_BAD.sub("", s), quote=True)


class _Docx:
    def __init__(self) -> None:
        self.links: List[str] = []

    def run(self, text: str, f: Dict[str, Any]) -> str:
        # rPr children in the order the schema requires (Word rejects others as unreadable)
        props = ""
        if f.get("href"):
            props += '<w:rStyle w:val="Hyperlink"/>'
        if f.get("code"):
            props += '<w:rFonts w:ascii="Consolas" w:hAnsi="Consolas" w:cs="Consolas"/>'
        if f.get("b"):
            props += "<w:b/><w:bCs/>"
        if f.get("i"):
            props += "<w:i/><w:iCs/>"
        if f.get("ev"):
            props += '<w:color w:val="6B7780"/><w:sz w:val="16"/><w:szCs w:val="16"/>'
        if f.get("code"):
            props += '<w:shd w:val="clear" w:color="auto" w:fill="F1F3F5"/>'
        body = "<w:br/>".join(f'<w:t xml:space="preserve">{_x(part)}</w:t>' for part in text.split("\n"))
        r = f"<w:r>{f'<w:rPr>{props}</w:rPr>' if props else ''}{body}</w:r>"
        if f.get("href"):
            self.links.append(f["href"])
            r = f'<w:hyperlink r:id="rIdL{len(self.links)}" w:history="1">{r}</w:hyperlink>'
            if text.strip() != f["href"]:
                r += self.run(f" ({f['href']})", {"ev": True})  # a printed letter needs the address itself
        return r

    def para(self, runs: Sequence[Run], style: Optional[str] = None, ppr: str = "", bold: bool = False) -> str:
        pp = (f'<w:pStyle w:val="{style}"/>' if style else "") + ppr
        body = "".join(self.run(t, {**f, "b": True} if bold else f) for t, f in runs)
        return f"<w:p>{f'<w:pPr>{pp}</w:pPr>' if pp else ''}{body}</w:p>"

    def document(self, blocks: Sequence[Block], footer_note: str) -> str:
        out: List[str] = []
        for b in blocks:
            kind = b[0]
            if kind == "h":
                out.append(self.para(inline(b[2]), f"Heading{min(b[1], 4)}" if b[1] > 1 else "Title"))
            elif kind == "p":
                out.append(self.para(inline(b[1])))
            elif kind == "list":
                counters: Dict[int, Tuple[bool, int]] = {}
                for depth, ordered, _marker, text in b[1]:
                    for d in [d for d in counters if d > depth]:
                        del counters[d]
                    was, n = counters.get(depth, (ordered, 0))
                    counters[depth] = (ordered, n + 1 if was == ordered else 1)
                    bullet = f"{counters[depth][1]}." if ordered else ("•" if depth % 2 == 0 else "◦")
                    left = 420 + 360 * depth
                    out.append(self.para([(bullet + "\t", {})] + inline(text), "ListParagraph",
                                         f'<w:tabs><w:tab w:val="left" w:pos="{left}"/></w:tabs>'
                                         f'<w:ind w:left="{left}" w:hanging="300"/>'))
            elif kind == "table":
                ncols = max(len(b[1]), 1)
                width = 9638 // ncols
                grid = "".join(f'<w:gridCol w:w="{width}"/>' for _ in range(ncols))

                def row(cells: Sequence[str], header: bool) -> str:
                    tr = '<w:trPr><w:tblHeader/></w:trPr>' if header else ""
                    tcs = "".join(
                        f'<w:tc><w:tcPr><w:tcW w:w="{width}" w:type="dxa"/>'
                        + ('<w:shd w:val="clear" w:color="auto" w:fill="F1F4F6"/>' if header else "")
                        + f"</w:tcPr>{self.para(inline(c), 'TableText', bold=header)}</w:tc>"
                        for c in cells)
                    return f"<w:tr>{tr}{tcs}</w:tr>"

                # borders on the table itself, not only in its style: WPS, Pages and previewers
                # that ignore table styles still draw the grid
                borders = "".join(f'<w:{e} w:val="single" w:sz="4" w:space="0" w:color="C9D0D6"/>'
                                  for e in ("top", "left", "bottom", "right", "insideH", "insideV"))
                out.append('<w:tbl><w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="5000" w:type="pct"/>'
                           f'<w:tblBorders>{borders}</w:tblBorders><w:tblLayout w:type="autofit"/></w:tblPr>'
                           f"<w:tblGrid>{grid}</w:tblGrid>{row(b[1], True)}"
                           + "".join(row(r, False) for r in b[2]) + "</w:tbl>")
                out.append(self.para([], ppr='<w:spacing w:after="0"/>'))
            elif kind == "quote":
                for p in [p for p in b[1].split("\n") if p.strip()]:
                    out.append(self.para(inline(p), "Quote"))
            elif kind == "code":
                out.append(self.para([(b[1], {"code": True})], "Code"))
            elif kind == "hr":
                out.append(self.para([], ppr='<w:pBdr><w:bottom w:val="single" w:sz="6" w:space="1" w:color="D9DEE3"/></w:pBdr>'))
        out.append(self.para([(footer_note, {"ev": True})]))
        sect = ('<w:sectPr><w:pgSz w:w="11906" w:h="16838"/>'
                '<w:pgMar w:top="1134" w:right="1134" w:bottom="1134" w:left="1134" w:header="567" w:footer="567" w:gutter="0"/>'
                "</w:sectPr>")
        return (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<w:document xmlns:w="{_W}" xmlns:r="{_R}">'
                f"<w:body>{''.join(out)}{sect}</w:body></w:document>")


def _styles(east_asia: str) -> str:
    def heading(sid: str, name: str, size: int, before: int, outline: Optional[int]) -> str:
        lvl = f'<w:outlineLvl w:val="{outline}"/>' if outline is not None else ""
        return (f'<w:style w:type="paragraph" w:styleId="{sid}"><w:name w:val="{name}"/><w:basedOn w:val="Normal"/>'
                f'<w:next w:val="Normal"/><w:qFormat/><w:pPr><w:keepNext/><w:spacing w:before="{before}" w:after="120"/>{lvl}</w:pPr>'
                f'<w:rPr><w:b/><w:bCs/><w:sz w:val="{size}"/><w:szCs w:val="{size}"/></w:rPr></w:style>')

    return (
        f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<w:styles xmlns:w="{_W}">'
        '<w:docDefaults><w:rPrDefault><w:rPr>'
        f'<w:rFonts w:ascii="Calibri" w:hAnsi="Calibri" w:eastAsia="{east_asia}" w:cs="Calibri"/>'
        '<w:sz w:val="22"/><w:szCs w:val="22"/><w:lang w:val="en-US" w:eastAsia="zh-CN"/></w:rPr></w:rPrDefault>'
        '<w:pPrDefault><w:pPr><w:spacing w:after="120" w:line="300" w:lineRule="auto"/></w:pPr></w:pPrDefault></w:docDefaults>'
        '<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:qFormat/></w:style>'
        + heading("Title", "Title", 36, 0, 0)
        + heading("Heading2", "heading 2", 28, 360, 1)
        + heading("Heading3", "heading 3", 24, 240, 2)
        + heading("Heading4", "heading 4", 22, 200, 3)
        + '<w:style w:type="paragraph" w:styleId="ListParagraph"><w:name w:val="List Paragraph"/><w:basedOn w:val="Normal"/>'
          '<w:pPr><w:spacing w:after="60"/></w:pPr></w:style>'
        + '<w:style w:type="paragraph" w:styleId="TableText"><w:name w:val="Table Text"/><w:basedOn w:val="Normal"/>'
          '<w:pPr><w:spacing w:after="0" w:line="260" w:lineRule="auto"/></w:pPr><w:rPr><w:sz w:val="20"/><w:szCs w:val="20"/></w:rPr></w:style>'
        + '<w:style w:type="paragraph" w:styleId="Quote"><w:name w:val="Quote"/><w:basedOn w:val="Normal"/>'
          '<w:pPr><w:pBdr><w:left w:val="single" w:sz="18" w:space="8" w:color="8AA4B8"/></w:pBdr>'
          '<w:shd w:val="clear" w:color="auto" w:fill="F5F8FA"/><w:ind w:left="240"/></w:pPr></w:style>'
        + '<w:style w:type="paragraph" w:styleId="Code"><w:name w:val="Code"/><w:basedOn w:val="Normal"/>'
          '<w:pPr><w:shd w:val="clear" w:color="auto" w:fill="F1F3F5"/><w:spacing w:after="120" w:line="240" w:lineRule="auto"/></w:pPr>'
          '<w:rPr><w:rFonts w:ascii="Consolas" w:hAnsi="Consolas" w:cs="Consolas"/><w:sz w:val="18"/></w:rPr></w:style>'
        + '<w:style w:type="character" w:styleId="Hyperlink"><w:name w:val="Hyperlink"/>'
          '<w:rPr><w:color w:val="1F5F8B"/><w:u w:val="single"/></w:rPr></w:style>'
        + '<w:style w:type="table" w:styleId="TableGrid"><w:name w:val="Table Grid"/><w:tblPr><w:tblBorders>'
        + "".join(f'<w:{s} w:val="single" w:sz="4" w:space="0" w:color="C9D0D6"/>'
                  for s in ("top", "left", "bottom", "right", "insideH", "insideV"))
        + '</w:tblBorders><w:tblCellMar><w:left w:w="100" w:type="dxa"/><w:right w:w="100" w:type="dxa"/></w:tblCellMar>'
          "</w:tblPr></w:style></w:styles>"
    )


def to_docx(blocks: Sequence[Block], title: str, path: Path, generated: str, lang: str) -> None:
    d = _Docx()
    note = ("生成于 " if lang.startswith("zh") else "Generated ") + generated + " · zebra-mod"
    document = d.document(blocks, note)
    links = "".join(f'<Relationship Id="rIdL{i}" Type="{_R}/hyperlink" Target="{_x(u)}" TargetMode="External"/>'
                    for i, u in enumerate(d.links, 1))
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    files = {
        "[Content_Types].xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
            '<Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>'
            '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
            "</Types>"),
        "_rels/.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'<Relationship Id="rId1" Type="{_R}/officeDocument" Target="word/document.xml"/>'
            '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>'
            "</Relationships>"),
        "word/_rels/document.xml.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'<Relationship Id="rIdS" Type="{_R}/styles" Target="styles.xml"/>{links}</Relationships>'),
        "word/document.xml": document,
        "word/styles.xml": _styles("Microsoft YaHei"),
        "docProps/core.xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
            'xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '
            'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
            f"<dc:title>{_x(title)}</dc:title><dc:creator>zebra-mod</dc:creator>"
            f'<dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>'
            f'<dcterms:modified xsi:type="dcterms:W3CDTF">{now}</dcterms:modified></cp:coreProperties>'),
    }
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for name, text in files.items():
            z.writestr(name, text.encode("utf-8"))
    os.replace(tmp, path)


# ------------------------------------------------------------------ PDF

_BROWSERS_MAC = [
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
    "/Applications/Chromium.app/Contents/MacOS/Chromium",
    "/Applications/Brave Browser.app/Contents/MacOS/Brave Browser",
]
_BROWSERS_WIN = [
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
]
_BROWSER_NAMES = ["google-chrome", "google-chrome-stable", "chromium", "chromium-browser", "microsoft-edge",
                  "microsoft-edge-stable", "msedge", "brave-browser", "chrome"]
_OFFICE = ["/Applications/LibreOffice.app/Contents/MacOS/soffice", r"C:\Program Files\LibreOffice\program\soffice.exe"]


def find_pdf_engine() -> Optional[Tuple[str, str]]:
    """("browser" | "office", executable) for the first PDF engine on this machine, or None."""
    env = os.environ.get("ZEBRA_BROWSER")
    if env and os.path.isfile(env) and os.access(env, os.X_OK):
        return ("office" if "office" in os.path.basename(env).lower() else "browser", env)
    for p in (_BROWSERS_MAC if sys.platform == "darwin" else _BROWSERS_WIN if os.name == "nt" else []):
        if os.path.isfile(p):
            return "browser", p
    for name in _BROWSER_NAMES:
        found = shutil.which(name)
        if found:
            return "browser", found
    for p in _OFFICE:
        if os.path.isfile(p):
            return "office", p
    for name in ("soffice", "libreoffice"):
        found = shutil.which(name)
        if found:
            return "office", found
    return None


def _pdf_complete(path: Path) -> bool:
    try:
        data = path.read_bytes()
    except OSError:
        return False
    return len(data) > 100 and data[:5] == b"%PDF-" and b"%%EOF" in data[-1024:]


def _kill_group(proc: "subprocess.Popen[bytes]") -> None:
    if proc.poll() is not None:
        return
    if os.name == "posix":
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except OSError:
            pass
    else:
        proc.kill()
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        pass


def _run_isolated(argv: List[str], timeout: float, done: Optional[Path] = None) -> Tuple[int, str]:
    """Run argv in its own process group; the whole group is killed when it overruns, or as soon as
    `done` holds a complete PDF (a headless browser often keeps running after it has printed)."""
    kw: Dict[str, Any] = {"stdout": subprocess.DEVNULL, "stderr": subprocess.PIPE, "stdin": subprocess.DEVNULL}
    if os.name == "posix":
        kw["start_new_session"] = True
    proc = subprocess.Popen(argv, **kw)
    deadline = time.monotonic() + timeout
    last_size = -1
    try:
        while proc.poll() is None:
            if time.monotonic() > deadline:
                _kill_group(proc)
                return -1, f"timed out after {timeout:.0f} s"
            if done is not None and done.exists():
                size = done.stat().st_size
                if size == last_size and _pdf_complete(done):
                    _kill_group(proc)
                    return 0, ""
                last_size = size
            time.sleep(0.3)
        err = proc.stderr.read().decode("utf-8", "replace")[-400:] if proc.stderr else ""
        return proc.returncode, err
    finally:
        _kill_group(proc)
        if proc.stderr:
            proc.stderr.close()


def to_pdf(html_path: Path, docx_path: Optional[Path], pdf_path: Path,
           engine: Optional[Tuple[str, str]] = None) -> Tuple[bool, str]:
    """(made it, how or why not)."""
    engine = engine or find_pdf_engine()
    if engine is None:
        return False, ("no PDF engine on this machine (Chrome, Edge, Chromium, Brave or LibreOffice); "
                       "open the HTML or Word file and print to PDF, or set ZEBRA_BROWSER")
    kind, exe = engine
    with tempfile.TemporaryDirectory(prefix="zebra-pdf-") as tmp:
        if kind == "browser":
            out = Path(tmp) / "out.pdf"
            argv = [exe, "--headless=new", "--disable-gpu", "--no-first-run", "--no-default-browser-check",
                    "--disable-extensions", "--disable-background-networking", f"--user-data-dir={tmp}/profile",
                    "--no-pdf-header-footer", "--print-to-pdf-no-header", f"--print-to-pdf={out}",
                    html_path.resolve().as_uri()]
        else:
            source = docx_path or html_path
            argv = [exe, "--headless", "--norestore", f"-env:UserInstallation={Path(tmp, 'profile').as_uri()}",
                    "--convert-to", "pdf", "--outdir", tmp, str(source.resolve())]
            out = Path(tmp) / (source.stem + ".pdf")
        code, err = _run_isolated(argv, PDF_TIMEOUT, out if kind == "browser" else None)
        if not _pdf_complete(out):
            return False, f"{os.path.basename(exe)} did not produce a PDF (exit {code}): {err.strip()[-200:]}"
        shutil.copyfile(out, pdf_path)
    return True, f"printed by {os.path.basename(exe)}"


# ------------------------------------------------------------------ export

def export(md_path: str, formats: Sequence[str], out_dir: Optional[str] = None, title: Optional[str] = None,
           identifiers: Sequence[str] = (), allow_identifiers: bool = False,
           engine: Optional[Tuple[str, str]] = None) -> Dict[str, Any]:
    src = Path(md_path).expanduser()
    if not src.is_file():
        raise FileNotFoundError(f"no such report: {src}")
    if src.stat().st_size > MAX_BYTES:
        raise ValueError(f"{src.name} is over {MAX_BYTES // (1024 * 1024)} MB; a report is text")
    text = src.read_text("utf-8", errors="replace")
    hit, warnings = identifier_hits(text, identifiers)
    if hit and not allow_identifiers:
        raise PermissionError(
            f"{src.name} contains {hit}. A handout leaves this machine once it is printed or sent and never needs "
            "it: remove it from the report and export again (an internal copy that must keep it: --allow-identifiers)")
    if hit:
        warnings.append(f"exported with {hit} in it (--allow-identifiers): keep this copy internal")
    blocks = parse(text)
    name = title or title_of(blocks, src.stem)
    lang = "zh-CN" if len(_CJK.findall(text)) > len(text) // 20 else "en"
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    dest = Path(out_dir).expanduser() if out_dir else src.parent
    dest.mkdir(parents=True, exist_ok=True)
    stem = src.stem
    made: List[Dict[str, Any]] = []
    skipped: List[Dict[str, str]] = []
    html_path = dest / f"{stem}.html"
    docx_path = dest / f"{stem}.docx"
    want = [f for f in ("html", "docx", "pdf") if f in formats]
    need_html = "html" in want or "pdf" in want
    tmp_html: Optional[Path] = None
    if need_html:
        page = to_html(blocks, name, lang, generated)
        if "html" in want:
            html_path.write_text(page, "utf-8")
            made.append({"format": "html", "path": str(html_path), "bytes": html_path.stat().st_size})
        else:
            fd, raw = tempfile.mkstemp(prefix="zebra-", suffix=".html")
            os.close(fd)
            tmp_html = Path(raw)
            tmp_html.write_text(page, "utf-8")
    if "docx" in want:
        to_docx(blocks, name, docx_path, generated, lang)
        made.append({"format": "docx", "path": str(docx_path), "bytes": docx_path.stat().st_size})
    if "pdf" in want:
        pdf_path = dest / f"{stem}.pdf"
        try:
            ok, how = to_pdf(tmp_html or html_path, docx_path if "docx" in want else None, pdf_path, engine)
        finally:
            if tmp_html is not None:
                tmp_html.unlink()
        if ok:
            made.append({"format": "pdf", "path": str(pdf_path), "bytes": pdf_path.stat().st_size, "how": how})
        else:
            skipped.append({"format": "pdf", "reason": how})
            warnings.append(f"PDF not made: {how}")
    return {"title": name, "language": lang, "files": made, "skipped": skipped, "warnings": warnings,
            "blocks": len(blocks)}
