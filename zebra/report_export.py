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
                rows.append(_cells(lines[i]))
                i += 1
            width = max([len(header)] + [len(r) for r in rows])
            header = header + [""] * (width - len(header))
            rows = [r + [""] * (width - len(r)) for r in rows]
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
    r"|(?P<link>\[(?P<lt>[^\]]+)\]\((?P<lu>(?:[^()\s]|\([^()\s]*\))+)\))"
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

_RESIDENT_ID = re.compile(r"(?<!\d)[1-9]\d{5}[ -]?(?:18|19|20)\d{2}[ -]?(?:0[1-9]|1[0-2])[ -]?(?:0[1-9]|[12]\d|3[01])"
                          r"[ -]?\d{3}[\dXx](?!\d)")
_MOBILE = re.compile(r"(?<!\d)1[3-9]\d{9}(?!\d)")
# identifiers of science, not of people: removed before looking for numbers and dates in them.
# A coordinate needs "chr" or a colon, so a record number like 12-0042317 is not taken for one.
_SCIENTIFIC_IDS = re.compile("|".join([
    r"\b(?:hp|omim|mim|orpha|orphanet|mondo|doid|ncit|efo|uberon|go|hgnc|cl|chebi|mp)\s*[:_]\s*\d+",
    r"\brs\d+", r"\b(?:nm|nr|np|nc|ng|nt|xm|xp|enst|ensg|ensp|ccds|lrg)_?\d+(?:\.\d+)?",
    r"\bpmid\s*:?\s*\d+", r"\bpmc\d+", r"\bnct\d{8}\b", r"\bchictr-?\w+",
    r"\bchr(?:\d{1,2}|x|y|mt?)\s*[:-]\s*[\d,]+(?:\s*[-_]\s*[\d,]+)?",
    r"(?<![\w-])(?:\d{1,2}|x|y|mt?)\s*:\s*[\d,]+(?:\s*[-_]\s*[\d,]+)?",
    r"\b[cgpnmr]\.\S+", r"\b10\.\d{4,9}/\S+",
]), re.I)
_MONTHS = ["january", "february", "march", "april", "may", "june", "july", "august", "september", "october",
           "november", "december"]


def _fold(s: str) -> str:
    """NFKC, lower case, tone marks and other diacritics removed (Zhāng → zhang)."""
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    return unicodedata.normalize("NFKC", s).casefold()


def _date_parts(value: str) -> List[Tuple[int, int, int]]:
    t = value.strip()
    m = re.match(r"^(\d{4})\s*[-/.年]\s*(\d{1,2})\s*[-/.月]\s*(\d{1,2})\s*日?$", t)
    if m:
        return [(int(m.group(1)), int(m.group(2)), int(m.group(3)))]
    m = re.match(r"^(\d{1,2})[-/.](\d{1,2})[-/.](\d{4})$", t)
    if m:
        a, b, y = int(m.group(1)), int(m.group(2)), int(m.group(3))
        return [(y, b, a), (y, a, b)]  # day-month-year and month-day-year
    return []


def _date_matchers(value: str) -> List["re.Pattern[str]"]:
    """A registered date in the ways people write one, with digit boundaries (2019-3-25 is not 2019-03-02)."""
    out = []
    for y, m, d in _date_parts(value):
        if not (1 <= m <= 12 and 1 <= d <= 31):
            continue
        mo, da, yy = f"0?{m}", f"0?{d}", str(y)[2:]
        sep = r"\s*[-/.]\s*"
        loose = r"[\s\-/.,]*"
        name = f"(?:{_MONTHS[m - 1]}|{_MONTHS[m - 1][:3]}\\.?)"
        ord_ = r"(?:st|nd|rd|th)?"
        out += [re.compile(p) for p in (
            rf"(?<!\d){y}\s*[-/.年]\s*{mo}\s*[-/.月]\s*{da}(?!\d)",
            rf"(?<!\d){y}{m:02d}{d:02d}(?!\d)",
            rf"(?<!\d){da}{sep}{mo}{sep}{y}(?!\d)",
            rf"(?<!\d){mo}{sep}{da}{sep}{y}(?!\d)",
            rf"(?<!\d){yy}[-/.]{mo}[-/.]{da}(?!\d)",
            rf"(?<![a-z]){name}{loose}{da}{ord_}{loose}{y}(?!\d)",
            rf"(?<!\d){da}{ord_}{loose}(?:of\s*)?{name}{loose}{y}(?!\d)",
        )]
    return out


def _paragraphs(text: str) -> List[str]:
    """What a reader sees: lines joined within a paragraph (a soft line break is a space in HTML/Word)."""
    out, cur = [], []
    for line in text.split("\n"):
        if line.strip():
            cur.append(line.strip())
        elif cur:
            out.append(" ".join(cur))
            cur = []
    if cur:
        out.append(" ".join(cur))
    return out


def identifier_hits(text: str, identifiers: Sequence[str]) -> Tuple[Optional[str], List[str]]:
    """(why the text must not be exported, or None; warnings).

    Checked paragraph by paragraph as a reader sees them (soft line breaks joined), and for a
    Chinese name across the whole text with all spacing removed. Latin names are matched as
    words whose boundaries are ASCII letters and digits, so 患儿Zhang Xiaoming still matches.
    """
    warnings: List[str] = []
    paras = [_fold(p) for p in _paragraphs(text)]
    scrubbed = [_SCIENTIFIC_IDS.sub(" ", p) for p in paras]  # NM_…, rs…, chr…: never a person's number
    names = [p.replace("_", " ") for p in paras]  # Zhang_Xiaoming is a name; NM_000492 is not touched above
    tight_all = re.sub(r"[\s_]+", "", "".join(paras))
    digit_paras = [re.sub(r"\D", "", p) for p in scrubbed]
    for n, raw in enumerate(identifiers, 1):
        value = _fold(str(raw or "").strip())
        if len(value) < 2:
            continue
        label = f"protected identifier #{n} of the case"
        dates = _date_matchers(value)
        if dates:
            if any(d.search(p) for d in dates for p in scrubbed):
                return f"{label} (a date)", warnings
            continue
        digits = re.sub(r"\D", "", value)
        if len(digits) >= 6 and len(digits) >= len(re.sub(r"[\s\-_/]", "", value)) - 3:
            if any(digits in d for d in digit_paras):
                return f"{label} (a number)", warnings
            continue
        if _CJK.search(value):
            if re.sub(r"\s+", "", value) in tight_all:
                return label, warnings
            continue
        tokens = [t for t in re.split(r"[\s,;]+", value) if len(t) >= 2]
        if not tokens or (len(tokens) == 1 and len(tokens[0]) < 3):
            continue
        pats = [re.compile(r"(?<![a-z0-9])" + re.escape(t) + r"(?![a-z0-9])") for t in tokens]
        if any(all(p_.search(par) for p_ in pats) for par in names):
            return label, warnings
    for par in scrubbed:
        if _RESIDENT_ID.search(par):
            return "a Chinese resident ID number", warnings
    if any(_MOBILE.search(par) for par in scrubbed):
        warnings.append("the report contains what looks like a mobile phone number: keep it only if it is an "
                        "organisation's public number, never the family's")
    return None, warnings


# ------------------------------------------------------------------ HTML

_CSS = """
:root { color-scheme: light;
  --ink: #1b2430; --muted: #5d6b78; --faint: #8b97a3; --rule: #dfe5ea; --rule-strong: #c5cfd8;
  --accent: #0f5c6e; --accent-soft: #e8f2f4; --paper: #ffffff; --wash: #f6f8f9; --stripe: #f8fafb;
  --ok: #1d6b45; --ok-bg: #e7f4ec; --src: #2a5d8f; --src-bg: #e8f0f8; --warn: #8a5a00; --warn-bg: #fdf3dc;
  --sans: "PingFang SC", "Hiragino Sans GB", "Noto Sans CJK SC", "Source Han Sans SC", "Microsoft YaHei",
          -apple-system, "Segoe UI", Roboto, "Helvetica Neue", sans-serif; }
* { box-sizing: border-box; }
html { background: var(--wash); }
body { margin: 0; color: var(--ink); background: var(--wash); font: 15px/1.72 var(--sans);
  -webkit-font-smoothing: antialiased; text-rendering: optimizeLegibility; }
.page { max-width: 60rem; margin: 2rem auto; background: var(--paper); padding: 2.6rem 3rem 2rem;
  border: 1px solid var(--rule); border-radius: 10px; box-shadow: 0 1px 2px rgba(16,24,40,.04); }
.doc-head { border-bottom: 3px solid var(--accent); padding-bottom: 1rem; margin-bottom: 1.4rem; }
.kicker { font-size: .74rem; letter-spacing: .12em; text-transform: uppercase; color: var(--accent);
  font-weight: 600; margin-bottom: .35rem; }
h1 { font-size: 1.75rem; line-height: 1.3; margin: 0; font-weight: 700; letter-spacing: -.01em; }
.subtitle { color: var(--muted); font-size: .92rem; margin: .5rem 0 0; }
h2 { font-size: 1.22rem; line-height: 1.4; margin: 2.2rem 0 .7rem; padding: 0 0 .35rem; color: var(--accent);
  border-bottom: 1px solid var(--rule-strong); font-weight: 700; }
h3 { font-size: 1.04rem; margin: 1.5rem 0 .45rem; font-weight: 700; }
h4, h5, h6 { font-size: .97rem; margin: 1.1rem 0 .3rem; font-weight: 600; color: var(--muted); }
p { margin: .5rem 0; }
ul, ol { margin: .35rem 0 .7rem; padding-left: 1.4rem; }
li { margin: .18rem 0; }
li::marker { color: var(--accent); }
strong { font-weight: 650; }
a { color: var(--src); text-decoration: none; border-bottom: 1px solid rgba(42,93,143,.25); overflow-wrap: anywhere; }
.table { margin: .8rem 0 1.1rem; overflow-x: auto; }
table { border: 1px solid var(--rule-strong); border-collapse: collapse; width: 100%; font-size: .86rem; line-height: 1.55; table-layout: fixed; }
th, td { padding: .42rem .6rem; text-align: left; vertical-align: top; overflow-wrap: anywhere; word-break: normal; }
thead th { background: var(--accent-soft); color: #134451; font-weight: 650; border-bottom: 1px solid var(--rule-strong);
  font-size: .8rem; }
tbody td { border-top: 1px solid var(--rule); }
tbody tr:first-child td { border-top: 0; }
tbody tr:nth-child(even) td { background: var(--stripe); }
td + td, th + th { border-left: 1px solid var(--rule); }
table.wide { font-size: .8rem; }
blockquote { margin: 1rem 0; padding: .7rem 1rem; border-left: 4px solid var(--accent); background: var(--accent-soft);
  border-radius: 0 6px 6px 0; font-size: .92rem; }
blockquote p { margin: .25rem 0; }
code { font: .86em/1.4 ui-monospace, "SF Mono", Menlo, Consolas, monospace; background: #eef1f4; padding: .05em .3em;
  border-radius: 3px; overflow-wrap: anywhere; }
pre { background: #eef1f4; padding: .7rem .9rem; overflow-x: auto; border-radius: 6px; white-space: pre-wrap; }
pre code { background: none; padding: 0; }
hr { border: 0; border-top: 1px solid var(--rule); margin: 1.6rem 0; }
.ev { color: var(--faint); font-size: .78em; white-space: nowrap; }
.tag { display: inline-block; font-size: .72em; line-height: 1.5; padding: 0 .42em; margin: 0 .12em; border-radius: 4px;
  font-weight: 600; white-space: nowrap; vertical-align: .08em; }
.tag-ok { color: var(--ok); background: var(--ok-bg); }
.tag-src { color: var(--src); background: var(--src-bg); }
.tag-warn { color: var(--warn); background: var(--warn-bg); }
.tag-note { color: var(--faint); font-size: .8em; overflow-wrap: anywhere; }
.toc { background: var(--wash); border: 1px solid var(--rule); border-radius: 8px; padding: .8rem 1.1rem; margin: 0 0 1.6rem;
  font-size: .88rem; }
.toc-title { font-weight: 650; color: var(--muted); font-size: .76rem; letter-spacing: .08em; margin-bottom: .3rem; }
.toc ol { list-style: none; margin: 0; padding: 0; columns: 2; column-gap: 2rem; }
.toc li { margin: .1rem 0; break-inside: avoid; }
.toc a { border: 0; color: var(--ink); }
.toc a:hover { color: var(--accent); }
.doc-foot { margin-top: 1.6rem; color: var(--faint); font-size: .8rem; }
.urgent { margin: 1.6rem 0; padding: .2rem 1.1rem .6rem; border: 1px solid #f0c4c0; border-left: 5px solid #b42318;
  background: #fef4f3; border-radius: 0 8px 8px 0; break-inside: avoid; }
.urgent h2 { color: #b42318; border-bottom-color: #f0c4c0; margin-top: .9rem; }
.urgent li::marker { color: #b42318; }
@media (max-width: 720px) {
  .page { margin: 0; border: 0; border-radius: 0; padding: 1.4rem 16px; }
  table { table-layout: auto; min-width: 36rem; } .toc ol { columns: 1; } h1 { font-size: 1.45rem; } }
@page { size: A4; margin: 16mm 14mm 18mm; }
@media print {
  html, body { background: #fff; }
  body { font-size: 10.2pt; line-height: 1.62; }
  .page { max-width: none; margin: 0; padding: 0; border: 0; box-shadow: none; }
  .toc { display: none; }
  h1 { font-size: 19pt; } h2 { font-size: 13pt; margin-top: 1.4rem; } h3 { font-size: 11pt; }
  table { font-size: 8.6pt; } table.wide { font-size: 8pt; } thead th { font-size: 8pt; }
  th, td { padding: .3rem .45rem; }
  thead { display: table-header-group; }
  tr, blockquote, .tag { break-inside: avoid; }
  h2, h3, h4 { break-after: avoid; }
  p, li { orphans: 2; widows: 2; }
  a { color: inherit; border: 0; }
  .table { overflow: visible; border-radius: 0; }
  tbody tr:nth-child(even) td { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
  thead th, blockquote, .tag, .doc-head, .urgent { -webkit-print-color-adjust: exact; print-color-adjust: exact; } }
"""

# provenance labels the report skills write in square brackets: [已检索：source] [凭记忆] [原始报告] ...
_TAGS = [
    (re.compile(r"^(已检索|retrieved|checked this session)", re.I), "tag-ok"),
    (re.compile(r"^(原始报告|报告原文|original report|from the report)", re.I), "tag-src"),
    (re.compile(r"^(凭记忆|未核实|需核实|待核实|from memory|unverified|to verify|not checked)", re.I), "tag-warn"),
]
_TAG_RE = re.compile(r"\[(已检索|原始报告|报告原文|凭记忆|未核实|需核实|待核实|retrieved|original report|from the report|"
                     r"from memory|unverified|to verify|not checked this session)([：:，,；;][^\[\]]*)?\]", re.I)
_BARE_URL = re.compile(r"(?<![(<\"'=])\bhttps?://[^\s\]\)）\"'<>，。；]+", re.I)


def _short_url(url: str) -> str:
    s = re.sub(r"^https?://(www\.)?", "", url, flags=re.I).rstrip("/")
    return s if len(s) <= 48 else s[:30] + "…" + s[-14:]


def _text_html(text: str) -> str:
    """Escape plain text, turning bare URLs into short links and provenance labels into tags."""
    out: List[str] = []
    pos = 0
    for m in _TAG_RE.finditer(text):
        out.append(_urls_html(text[pos:m.start()]))
        label = m.group(1)
        cls = next((c for rx, c in _TAGS if rx.match(label)), "tag-warn")
        note = (m.group(2) or "")[1:].strip()
        out.append(f'<span class="tag {cls}">{html.escape(label)}</span>'
                   + (f'<span class="tag-note">{_urls_html(note)}</span>' if note else ""))
        pos = m.end()
    out.append(_urls_html(text[pos:]))
    return "".join(out)


def _urls_html(text: str) -> str:
    out: List[str] = []
    pos = 0
    for m in _BARE_URL.finditer(text):
        url = m.group(0).rstrip(".,;:")
        out.append(html.escape(text[pos:m.start()]))
        out.append(f'<a href="{html.escape(url, quote=True)}">{html.escape(_short_url(url))}</a>')
        pos = m.start() + len(url)
    out.append(html.escape(text[pos:]))
    return "".join(out)


def _html_runs(runs: Sequence[Run]) -> str:
    out = []
    for text, f in runs:
        s = html.escape(text) if (f.get("code") or f.get("href")) else _text_html(text)
        if f.get("code"):
            s = f"<code>{s}</code>"
        if f.get("ev"):
            s = f'<span class="ev">{s}</span>'
        if f.get("i"):
            s = f"<em>{s}</em>"
        if f.get("b"):
            s = f"<strong>{s}</strong>"
        if f.get("href"):
            shown = _short_url(text) if text.strip() == f["href"] else None
            s = f'<a href="{html.escape(f["href"], quote=True)}">{html.escape(shown) if shown else s}</a>'
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


def _plain(cell: str) -> str:
    return "".join(t for t, _ in inline(cell))


def _display_len(s: str) -> float:
    """Width in Latin-character units: a CJK character is about two, a URL breaks anywhere."""
    s = _BARE_URL.sub("x" * 24, s)
    return sum(2.0 if _CJK.match(ch) else 1.0 for ch in s)


def col_widths(header: Sequence[str], rows: Sequence[Sequence[str]]) -> List[float]:
    """Column widths in percent from the content: wide enough for the header and the longest
    unbreakable word, the rest shared by how much text each column holds."""
    n = len(header)
    if n == 0:
        return []
    floors, weights = [], []
    for c in range(n):
        cells = [_plain(r[c]) for r in rows if c < len(r)]
        head = _plain(header[c])
        word = max([len(w) for x in cells + [head] for w in re.findall(r"[A-Za-z0-9.+\-/%≤≥<>]+", x)] + [0])
        cjk_head = _display_len(head)
        floor = max(min(word, 16), min(cjk_head, 12), 3.5) + 2.0
        lens = sorted(_display_len(x) for x in cells) or [0.0]
        mean = sum(lens) / len(lens)
        p90 = lens[min(len(lens) - 1, int(len(lens) * .9))]
        if lens[-1] <= 16:  # a column of short labels: wide enough that none of them wraps
            floor = max(floor, lens[-1] + 2.0)
        floors.append(floor)
        weights.append(max(floor, .55 * mean + .45 * p90))
    # each column first gets its floor (as a share of a ~95-character line), then the rest of the
    # width is shared by how much more text it holds, damped so a long column cannot starve the others
    line = 95.0 if n < 5 else 105.0
    base = [100.0 * f / line for f in floors]
    if sum(base) >= 100.0:
        pct = [100.0 * x / sum(base) for x in base]
    else:
        extra = [max(0.0, w - f) ** .75 for w, f in zip(weights, floors)]
        room = 100.0 - sum(base)
        total = sum(extra) or 1.0
        pct = [x + room * e / total for x, e in zip(base, extra)]
    return [round(p, 1) for p in pct]


def _table_html(header: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    widths = col_widths(header, rows)
    cols = "".join(f'<col style="width:{w}%">' for w in widths)
    head = "".join(f"<th>{_html_runs(inline(c))}</th>" for c in header)
    body = "".join("<tr>" + "".join(f"<td>{_html_runs(inline(c))}</td>" for c in r) + "</tr>" for r in rows)
    cls = ' class="wide"' if len(header) >= 5 else ""
    return (f'<div class="table"><table{cls}><colgroup>{cols}</colgroup><thead><tr>{head}</tr></thead>'
            f"<tbody>{body}</tbody></table></div>")


# a section that tells the reader when to seek care now is drawn as a warning box
_URGENT = re.compile(r"马上就医|立即就医|立刻就医|急诊|紧急|危险信号|红旗|emergency|urgent|red flags|seek care now|call 120|call 911",
                     re.I)


def _slug(i: int) -> str:
    return f"s{i}"


def _css_string(s: str) -> str:
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ") + '"'


def to_html(blocks: Sequence[Block], title: str, lang: str, generated: str) -> str:
    zh = lang.startswith("zh")
    blocks = list(blocks)
    # the H1 and a short line right under it (date, audience) make the header
    head_title = title
    subtitle = None
    if blocks and blocks[0][0] == "h" and blocks[0][1] == 1:
        head_title = blocks[0][2]
        blocks = blocks[1:]
        if blocks and blocks[0][0] == "p" and len(blocks[0][1]) <= 240:
            subtitle = blocks[0][1]
            blocks = blocks[1:]
    sections = [(i, b[2]) for i, b in enumerate(blocks) if b[0] == "h" and b[1] == 2]
    body: List[str] = []
    in_urgent = False
    for i, b in enumerate(blocks):
        kind = b[0]
        if kind == "h" and b[1] <= 2 and in_urgent:
            body.append("</section>")
            in_urgent = False
        if kind == "h":
            anchor = f' id="{_slug(i)}"' if b[1] == 2 else ""
            if b[1] == 2 and _URGENT.search(_plain(b[2])):
                body.append('<section class="urgent">')
                in_urgent = True
            body.append(f"<h{b[1]}{anchor}>{_html_runs(inline(b[2]))}</h{b[1]}>")
        elif kind == "p":
            body.append(f"<p>{_html_runs(inline(b[1]))}</p>")
        elif kind == "list":
            body.append(_html_list(b[1]))
        elif kind == "table":
            body.append(_table_html(b[1], b[2]))
        elif kind == "quote":
            body.append("<blockquote>" + "".join(f"<p>{_html_runs(inline(p))}</p>" for p in b[1].split("\n") if p.strip())
                        + "</blockquote>")
        elif kind == "code":
            body.append(f"<pre><code>{html.escape(b[1])}</code></pre>")
        elif kind == "hr":
            # a rule just before a section heading only doubles the heading's own line
            nxt = blocks[i + 1] if i + 1 < len(blocks) else None
            if not (nxt and nxt[0] == "h" and nxt[1] <= 2):
                if in_urgent:
                    body.append("</section>")
                    in_urgent = False
                body.append("<hr>")
    if in_urgent:
        body.append("</section>")
    toc = ""
    if len(sections) >= 4:
        items = "".join(f'<li><a href="#{_slug(i)}">{_html_runs(inline(t))}</a></li>' for i, t in sections)
        toc = f'<nav class="toc"><div class="toc-title">{"目录" if zh else "CONTENTS"}</div><ol>{items}</ol></nav>'
    kicker = "zebra-mod"
    plain_title = "".join(t for t, _ in inline(head_title))
    page_word = ("第 ", " 页 / 共 ", " 页") if zh else ("Page ", " of ", "")
    margin_css = (
        "@page { @top-left { content: " + _css_string(plain_title[:60]) + "; font: 7.5pt sans-serif; color: #8b97a3; }"
        " @bottom-left { content: " + _css_string(("生成于 " if zh else "Generated ") + generated) + "; font: 7.5pt sans-serif; color: #8b97a3; }"
        f' @bottom-right {{ content: "{page_word[0]}" counter(page) "{page_word[1]}" counter(pages) "{page_word[2]}";'
        " font: 7.5pt sans-serif; color: #8b97a3; } }"
        " @page :first { @top-left { content: none; } }"
    )
    head = (f'<header class="doc-head"><div class="kicker">{html.escape(kicker)}</div>'
            f"<h1>{_html_runs(inline(head_title))}</h1>"
            + (f'<p class="subtitle">{_html_runs(inline(subtitle))}</p>' if subtitle else "") + "</header>")
    meta = ("生成于 " if zh else "Generated ") + generated + " · zebra-mod"
    return (f'<!doctype html>\n<html lang="{lang}">\n<head>\n<meta charset="utf-8">\n'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">\n'
            f"<title>{html.escape(plain_title)}</title>\n<style>{_CSS}\n{margin_css}</style>\n</head>\n<body>\n"
            f'<main class="page">\n{head}\n{toc}\n' + "\n".join(body)
            + f'\n<footer class="doc-foot meta">{html.escape(meta)}</footer>\n</main>\n</body>\n</html>\n')


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
                pcts = col_widths(b[1], b[2]) or [100.0]
                widths = [max(400, int(9638 * p / 100)) for p in pcts]
                grid = "".join(f'<w:gridCol w:w="{w}"/>' for w in widths)

                def row(cells: Sequence[str], header: bool) -> str:
                    tr = '<w:trPr><w:tblHeader/></w:trPr>' if header else ""
                    tcs = "".join(
                        f'<w:tc><w:tcPr><w:tcW w:w="{w}" w:type="dxa"/>'
                        + ('<w:shd w:val="clear" w:color="auto" w:fill="E8F2F4"/>' if header else "")
                        + f"</w:tcPr>{self.para(inline(c), 'TableText', bold=header)}</w:tc>"
                        for c, w in zip(cells, widths))
                    return f"<w:tr>{tr}{tcs}</w:tr>"

                # borders on the table itself, not only in its style: WPS, Pages and previewers
                # that ignore table styles still draw the grid
                borders = "".join(f'<w:{e} w:val="single" w:sz="4" w:space="0" w:color="C9D0D6"/>'
                                  for e in ("top", "left", "bottom", "right", "insideH", "insideV"))
                out.append('<w:tbl><w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="5000" w:type="pct"/>'
                           f'<w:tblBorders>{borders}</w:tblBorders><w:tblLayout w:type="fixed"/></w:tblPr>'
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
    # the title and the file name are published too: in <title>, the Word metadata and the file names
    hit, warnings = identifier_hits("\n\n".join([text, title or "", src.stem]), identifiers)
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
