#!/usr/bin/env python3
"""Build zebra/data/china_rare_network_hospitals.json: the member hospitals of
全国罕见病诊疗协作网 (National Rare Disease Diagnosis and Treatment Collaboration
Network), parsed from the National Health Commission's own documents.

Rerun it from ~/zebra-mod:

    python3 -I tools/china_access/build_hospitals.py --dir <download-dir> --fetch
    python3 -I tools/china_access/build_hospitals.py --dir <download-dir>      # offline, files already there

    options: --out <json>   (default: zebra/data/china_rare_network_hospitals.json, relative to the repo)

--fetch downloads every source that a script can download (at least 1.2 s between requests to the
same host) into --dir and logs url / final url / HTTP status / UTC time in <dir>/fetch_log.json:

    govcn_2019_zhengceku.html        gov.cn 政策库 copy of 国卫办医函〔2019〕157号 (the 2019 table is inline)
    nhc_2024_64_att1_hospitals.docx  附件1 of 国卫办医政函〔2024〕64号, served by nhc.gov.cn without a challenge
    govcn_2024_yaowen.html           gov.cn news page stating the 419 total

nhc.gov.cn *pages* (.shtml) answer scripts with HTTP 412 (JavaScript challenge), so these were saved from
a real browser (ego-browser, 2026-10-06) and are read from --dir if present (their retrieval time is the
file mtime). Missing ones are reported in provenance.not_retrieved, and the check that needs them is
recorded as not performed:

    nhc_2019_157_notice_rendered.html       NHC original of the 2019 notice: second copy of the table
    nhc_2024_64_notice_rendered.html        NHC notice page carrying the 2024 docx (document number, link)
    nhc_2024_tia_02403_reply_rendered.html  NHC reply to CPPCC proposal 02403 (2024): "目前已有419家医院"
    nhc_2026_305_notice_rendered.html       国卫办医政函〔2026〕305号, checked to contain no hospital list
    nhc_yzygj_policy_list_2021-2026.json    titles/urls of NHC 医政司 政策文件 (search evidence, read only)

What the build asserts (it stops on any failure): the 2019 table numbers 1..324 with 1/32/291 roles as the
notice text states; the two tables on the gov.cn page are identical; the NHC original equals gov.cn row
by row; the 2024 table numbers 1..419 and totals 419 as NHC states; no duplicate name within a list; both
lists use the same 32 province strings; every role string is known.

Rows: the 2024 list (current membership) in its published order, then the 2019 hospitals whose name is
not in the 2024 list, in 2019 order. `name` is always exactly as published in the row's `source_id`.
A 2019 and a 2024 name are treated as the same hospital only on string evidence inside the two lists
(same province; identical, or one name / bracketed alias of at least 4 characters contained in the
other, or one name an in-order subsequence of the other of at least 6 characters), and only when the
pairing is one-to-one. No pairing comes from outside knowledge.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import re
import sys
import time
import urllib.request
import zipfile
from collections import Counter, OrderedDict
from datetime import datetime, timezone
from html.parser import HTMLParser
from xml.etree import ElementTree as ET

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
DEFAULT_OUT = os.path.join(REPO, "zebra", "data", "china_rare_network_hospitals.json")
UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/128.0 Safari/537.36")

URL_2019_GOVCN = "https://www.gov.cn/zhengce/zhengceku/2019-10/08/content_5436962.htm"
URL_2019_NHC = "https://www.nhc.gov.cn/yzygj/c100068/201902/39c01916f98a45c58418a839d82fda0c.shtml"
URL_2024_NHC = "https://www.nhc.gov.cn/yzygj/c100068/202403/e0328f505bcc47619e464148b2304dc0.shtml"
URL_2024_DOCX = ("https://www.nhc.gov.cn/yzygj/c100068/202403/e0328f505bcc47619e464148b2304dc0/"
                 "files/1732871019315_14550.docx")
URL_2024_TIA = "https://www.nhc.gov.cn/wjw/tia/202408/bb5a1989852b4c64a7f876f86c0cb93a.shtml"
URL_2024_GOVNEWS = "https://www.gov.cn/yaowen/liebiao/202410/content_6981620.htm"
URL_2026_NHC = "https://www.nhc.gov.cn/yzygj/c100068/202609/1bf21303ba184697963d65dd3b5bc1ad.shtml"
URL_POLICY_LIST = ("https://www.nhc.gov.cn/search/2e198e5e7e4647e48f135466e8b0f18f?_isAgg=true&_isJson=true"
                   "&_pageSize=300&_template=index&page=1")

F_2019_GOVCN = "govcn_2019_zhengceku.html"
F_2019_NHC = "nhc_2019_157_notice_rendered.html"
F_2024_DOCX = "nhc_2024_64_att1_hospitals.docx"
F_2024_NHC = "nhc_2024_64_notice_rendered.html"
F_2024_TIA = "nhc_2024_tia_02403_reply_rendered.html"
F_2024_GOVNEWS = "govcn_2024_yaowen.html"
F_2026_NHC = "nhc_2026_305_notice_rendered.html"
F_POLICY_LIST = "nhc_yzygj_policy_list_2021-2026.json"

FETCHABLE = [(URL_2019_GOVCN, F_2019_GOVCN), (URL_2024_DOCX, F_2024_DOCX), (URL_2024_GOVNEWS, F_2024_GOVNEWS)]

ROLE = {"国家级牵头医院": "national_lead", "省级牵头医院": "provincial_lead",
        "成员医院": "member", "成员单位": "member"}
HEADER = ["编号", "省份", "类型", "协作网医院"]
# published province string -> short name (strip 省/市; autonomous regions lose 自治区 and the
# ethnic qualifier); 新疆生产建设兵团 is not a province-level government and is kept whole.
AUTONOMOUS = {"内蒙古自治区": "内蒙古", "广西壮族自治区": "广西", "西藏自治区": "西藏",
              "宁夏回族自治区": "宁夏", "新疆维吾尔自治区": "新疆"}
PAREN = re.compile(r"[（(]([^）)]*)[）)]")


# ---------------------------------------------------------------- download / file bookkeeping

def utc_now():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def fetch_all(d):
    os.makedirs(d, exist_ok=True)
    log_path = os.path.join(d, "fetch_log.json")
    log = json.load(open(log_path, encoding="utf-8")) if os.path.exists(log_path) else {}
    last = {}
    for url, fname in FETCHABLE:
        host = url.split("/")[2]
        wait = 1.2 - (time.time() - last.get(host, 0))
        if wait > 0:
            time.sleep(wait)
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=60) as r:
            body = r.read()
            status, final = r.status, r.geturl()
        last[host] = time.time()
        if status != 200 or not body:
            sys.exit(f"fetch failed: {url} -> HTTP {status}, {len(body)} bytes")
        with open(os.path.join(d, fname), "wb") as fh:
            fh.write(body)
        log[fname] = {"url": url, "final_url": final, "status": status, "retrieved_at": utc_now(),
                      "bytes": len(body)}
        print(f"fetched {fname}: HTTP {status}, {len(body)} bytes")
    with open(log_path, "w", encoding="utf-8") as fh:
        json.dump(log, fh, ensure_ascii=False, indent=1)


def file_meta(d, fname, log):
    p = os.path.join(d, fname)
    b = open(p, "rb").read()
    ent = log.get(fname)
    if ent:
        when, how = ent["retrieved_at"], "downloaded by this script (fetch_log.json)"
    else:
        when = datetime.fromtimestamp(os.path.getmtime(p), timezone.utc).replace(microsecond=0)
        when, how = when.isoformat().replace("+00:00", "Z"), "saved from a browser; time is the file mtime"
    return b, {"retrieved_at": when, "sha256": hashlib.sha256(b).hexdigest(), "bytes": len(b), "how": how}


def visible_text(raw: bytes) -> str:
    t = raw.decode("utf-8", errors="replace")
    t = re.sub(r"(?is)<(script|style)\b.*?</\1>", "", t)
    t = html.unescape(re.sub(r"<[^>]+>", "", t))
    return re.sub(r"[\s　\xa0]+", "", t)


# ---------------------------------------------------------------- table parsing

def clean(s: str) -> str:
    # cell text as published, minus layout whitespace (the source cells contain no internal spaces)
    return re.sub(r"[\s　\xa0]+", "", html.unescape(s))


class _Tables(HTMLParser):
    """Collects every top-level <table> as rows of {text, rowspan, colspan} cells."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tables, self.depth, self.cur, self.row, self.cell = [], 0, None, None, None

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "table":
            self.depth += 1
            if self.depth == 1:
                self.cur = []
        elif tag == "tr" and self.cur is not None and self.depth == 1:
            self.row = []
        elif tag in ("td", "th") and self.row is not None:
            self.cell = {"text": "", "rowspan": int(a.get("rowspan") or 1), "colspan": int(a.get("colspan") or 1)}

    def handle_endtag(self, tag):
        if tag in ("td", "th") and self.cell is not None:
            self.row.append(self.cell)
            self.cell = None
        elif tag == "tr" and self.row is not None:
            self.cur.append(self.row)
            self.row = None
        elif tag == "table":
            self.depth -= 1
            if self.depth == 0 and self.cur is not None:
                self.tables.append(self.cur)
                self.cur = None

    def handle_data(self, data):
        if self.cell is not None:
            self.cell["text"] += data


def _expand(rows):
    """Resolve rowspan/colspan into a rectangular grid of strings."""
    grid, pending = [], {}
    for r in rows:
        line, col, cells = [], 0, list(r)
        while cells or col in pending:
            if col in pending:
                txt, left = pending[col]
                line.append(txt)
                if left > 1:
                    pending[col] = (txt, left - 1)
                else:
                    del pending[col]
                col += 1
                continue
            c = cells.pop(0)
            txt = clean(c["text"])
            for _ in range(c["colspan"]):
                line.append(txt)
                if c["rowspan"] > 1:
                    pending[col] = (txt, c["rowspan"] - 1)
                col += 1
        grid.append(line)
    return grid


def html_member_tables(raw: bytes):
    p = _Tables()
    p.feed(raw.decode("utf-8", errors="replace"))
    out = []
    for t in p.tables:
        g = _expand(t)
        if g and g[0] == HEADER:
            out.append(g[1:])
    return out


def docx_member_table(raw: bytes, path_for_errors: str):
    import io
    W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    root = ET.fromstring(zipfile.ZipFile(io.BytesIO(raw)).read("word/document.xml"))
    tables = list(root.iter(W + "tbl"))
    if len(tables) != 1:
        sys.exit(f"{path_for_errors}: expected 1 table, found {len(tables)}")
    paras = [clean("".join(x.text or "" for x in p.iter(W + "t"))) for p in root.find(W + "body").findall(W + "p")]
    grid, last = [], {}
    for tr in tables[0].findall(W + "tr"):
        line = []
        for i, tc in enumerate(tr.findall(W + "tc")):
            pr = tc.find(W + "tcPr")
            vm = pr.find(W + "vMerge") if pr is not None else None
            if pr is not None and pr.find(W + "gridSpan") is not None:
                sys.exit(f"{path_for_errors}: unexpected gridSpan")
            txt = clean("".join(x.text or "" for x in tc.iter(W + "t")))
            if vm is not None and vm.get(W + "val", "continue") == "continue":
                txt = last[i]          # vertically merged cell: repeats the cell above
            last[i] = txt
            line.append(txt)
        grid.append(line)
    if grid[0] != HEADER:
        sys.exit(f"{path_for_errors}: unexpected header {grid[0]}")
    return [p for p in paras if p], grid[1:]


def to_records(grid, label):
    recs = []
    for i, row in enumerate(grid, 1):
        if len(row) != 4:
            sys.exit(f"{label}: row {i} has {len(row)} cells: {row}")
        no, prov, role_zh, name = row
        if not no.isdigit() or int(no) != i:
            sys.exit(f"{label}: numbering breaks at row {i}: {row}")
        if role_zh not in ROLE:
            sys.exit(f"{label}: unknown role string {role_zh!r} at 编号 {no}")
        if not prov or not name:
            sys.exit(f"{label}: empty cell at 编号 {no}: {row}")
        recs.append({"no": i, "province_published": prov, "role_zh": role_zh, "name": name})
    return recs


def short_province(p: str) -> str:
    if p in AUTONOMOUS:
        return AUTONOMOUS[p]
    if p == "新疆生产建设兵团":
        return p
    if p.endswith("省") or p.endswith("市"):
        return p[:-1]
    sys.exit(f"unmapped province string {p!r}")


# ---------------------------------------------------------------- 2019 <-> 2024 pairing

def _aliases(name):
    return {PAREN.sub("", name)} | set(PAREN.findall(name))


def _is_subsequence(short, long_):
    it = iter(long_)
    return all(ch in it for ch in short)


def pair_basis(a: str, b: str):
    """String evidence that 2019 name `a` and 2024 name `b` denote one hospital, or None."""
    for x in _aliases(a):
        if len(x) >= 4 and x in b:
            return "2019_name_in_2024_name" if x == PAREN.sub("", a) else "2019_alias_in_2024_name"
    for y in _aliases(b):
        if len(y) >= 4 and y in a:
            return "2024_name_in_2019_name" if y == PAREN.sub("", b) else "2024_alias_in_2019_name"
    s, l = (a, b) if len(a) <= len(b) else (b, a)
    if len(s) >= 6 and _is_subsequence(s, l):
        return "2019_name_subsequence_of_2024" if s == a else "2024_name_subsequence_of_2019"
    return None


def match_lists(l19, l24):
    by_exact = {}
    idx24 = {(r["province_published"], r["name"]): r for r in l24}
    for r in l19:
        k = (r["province_published"], r["name"])
        if k in idx24:
            by_exact[r["no"]] = idx24[k]["no"]
    used19 = set(by_exact)
    used24 = set(by_exact.values())
    cand = {}
    for a in l19:
        if a["no"] in used19:
            continue
        for b in l24:
            if b["no"] in used24 or b["province_published"] != a["province_published"]:
                continue
            basis = pair_basis(a["name"], b["name"])
            if basis:
                cand.setdefault(a["no"], []).append((b["no"], basis))
    back = Counter(b for v in cand.values() for b, _ in v)
    renamed, ambiguous = {}, []
    for a_no, v in cand.items():
        if len(v) == 1 and back[v[0][0]] == 1:
            renamed[a_no] = v[0]
        else:
            ambiguous.append((a_no, v))
    return by_exact, renamed, ambiguous


# ---------------------------------------------------------------- build

def counts_of(rows):
    roles = Counter(r["role"] for r in rows)
    byp = OrderedDict()
    for r in rows:
        byp[r["province"]] = byp.get(r["province"], 0) + 1
    return OrderedDict([("total", len(rows)), ("national_lead", roles["national_lead"]),
                        ("provincial_lead", roles["provincial_lead"]), ("member", roles["member"]),
                        ("by_province", byp)])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="download directory (one dir only for these files)")
    ap.add_argument("--fetch", action="store_true", help="download the script-fetchable sources first")
    ap.add_argument("--out", default=DEFAULT_OUT)
    args = ap.parse_args(argv)
    d = os.path.abspath(args.dir)
    if args.fetch:
        fetch_all(d)
    log_path = os.path.join(d, "fetch_log.json")
    log = json.load(open(log_path, encoding="utf-8")) if os.path.exists(log_path) else {}
    have = lambda f: os.path.exists(os.path.join(d, f))
    for f in (F_2019_GOVCN, F_2024_DOCX):
        if not have(f):
            sys.exit(f"missing {f} in {d}; run with --fetch")
    not_retrieved = []

    # ---- 2019: gov.cn copy (parsed) + NHC original (cross-check)
    raw19, m19 = file_meta(d, F_2019_GOVCN, log)
    t19 = visible_text(raw19)
    for must in ("国卫办医函〔2019〕157号", "324家医院", "1家国家级牵头医院、32家省级牵头医院和291家协作网成员医院",
                 "附件：第一批罕见病诊疗协作网医院名单"):
        assert must in t19, f"gov.cn 2019 page lacks {must!r}"
    tabs = html_member_tables(raw19)
    assert len(tabs) == 2 and tabs[0] == tabs[1], "gov.cn page: expected the same member table twice"
    l19 = to_records(tabs[0], "2019 gov.cn")
    roles19 = Counter(ROLE[r["role_zh"]] for r in l19)
    assert len(l19) == 324 and roles19 == {"national_lead": 1, "provincial_lead": 32, "member": 291}, roles19
    if have(F_2019_NHC):
        rawn, mn = file_meta(d, F_2019_NHC, log)
        tn = html_member_tables(rawn)
        assert len(tn) == 1, f"NHC 2019 page: {len(tn)} member tables"
        diffs = [(a, b) for a, b in zip(tabs[0], tn[0]) if a != b]
        assert len(tn[0]) == len(tabs[0]) and not diffs, f"2019 copies differ: {diffs[:5]}"
        xc19 = (f"gov.cn page carries the attachment table twice: both copies identical (324 rows). "
                f"Compared row by row (编号, 省份, 类型, 名称) with the NHC original page {URL_2019_NHC} "
                f"(saved from a browser, sha256 {mn['sha256']}): 324/324 rows identical, so per-province and "
                f"per-role counts agree. Parsed roles 1/32/291 equal the notice text "
                f"'1家国家级牵头医院、32家省级牵头医院和291家协作网成员医院'.")
    else:
        mn = None
        xc19 = ("gov.cn duplicate table identical; NHC original not in --dir, so the second-copy comparison "
                "was not performed in this run. Parsed roles 1/32/291 equal the notice text.")
        not_retrieved.append({"what": "NHC original of 国卫办医函〔2019〕157号 (second copy of the table)",
                              "url": URL_2019_NHC, "reason": f"{F_2019_NHC} not in --dir (nhc.gov.cn pages answer scripts with HTTP 412)"})

    # ---- 2024: NHC docx (parsed)
    raw24, m24 = file_meta(d, F_2024_DOCX, log)
    paras24, g24 = docx_member_table(raw24, F_2024_DOCX)
    assert paras24[:2] == ["附件1", "全国罕见病诊疗协作网成员医院名单"], paras24[:3]
    l24 = to_records(g24, "2024 docx")
    roles24 = Counter(ROLE[r["role_zh"]] for r in l24)
    assert len(l24) == 419, len(l24)
    stated_419 = []
    if have(F_2024_TIA):
        rawt, mt = file_meta(d, F_2024_TIA, log)
        assert "目前已有419家医院纳入协作网范围" in visible_text(rawt)
        stated_419.append("nhc_reply_2024_02403")
    else:
        mt = None
        not_retrieved.append({"what": "NHC reply to CPPCC proposal 02403 (states 419)", "url": URL_2024_TIA,
                              "reason": f"{F_2024_TIA} not in --dir (HTTP 412 to scripts)"})
    if have(F_2024_GOVNEWS):
        rawg, mg = file_meta(d, F_2024_GOVNEWS, log)
        assert "全国罕见病诊疗协作网医院总数达到419家" in visible_text(rawg)
        stated_419.append("govcn_news_2024_419")
    else:
        mg = None
    assert stated_419, "no retrieved source states the 419 total"
    if have(F_2024_NHC):
        rawp, mp = file_meta(d, F_2024_NHC, log)
        tp = visible_text(rawp)
        assert "国卫办医政函〔2024〕64号" in tp and "全国罕见病诊疗协作网成员医院名单" in tp and "2024年2月28日" in tp
        assert URL_2024_DOCX.split("/files/")[1] in rawp.decode("utf-8", "replace")
    else:
        mp = None
        not_retrieved.append({"what": "NHC notice page of 国卫办医政函〔2024〕64号", "url": URL_2024_NHC,
                              "reason": f"{F_2024_NHC} not in --dir (HTTP 412 to scripts); docx used without it"})
    if have(F_2026_NHC):
        raw26, m26 = file_meta(d, F_2026_NHC, log)
        t26 = visible_text(raw26)
        assert "国卫办医政函〔2026〕305号" in t26
        assert "<table" not in raw26.decode("utf-8", "replace").lower() and "附件" not in t26, \
            "2026 notice now carries a table/attachment: check it for a membership list"
    else:
        m26 = None
        not_retrieved.append({"what": "国卫办医政函〔2026〕305号 (checked for a membership list)", "url": URL_2026_NHC,
                              "reason": f"{F_2026_NHC} not in --dir (HTTP 412 to scripts); the coverage statement "
                                        f"about it was not re-verified in this run"})
    if not have(F_POLICY_LIST):
        not_retrieved.append({"what": "NHC 医政司 政策文件 list (search for later membership updates)", "url": URL_POLICY_LIST,
                              "reason": f"{F_POLICY_LIST} not in --dir (browser-only); not re-searched in this run"})
    policy_hits = []
    if have(F_POLICY_LIST):
        pl = json.load(open(os.path.join(d, F_POLICY_LIST), encoding="utf-8"))
        policy_hits = [f"{x['date'][:10]} {x['title']}" for x in pl if re.search("罕见|协作网", x["title"])]
        pl_span = (min(x["date"] for x in pl)[:10], max(x["date"] for x in pl)[:10], len(pl))

    # ---- provinces
    p19 = list(OrderedDict.fromkeys(r["province_published"] for r in l19))
    p24 = list(OrderedDict.fromkeys(r["province_published"] for r in l24))
    assert set(p19) == set(p24) and len(p24) == 32, (set(p19) ^ set(p24))
    pmap = OrderedDict((p, short_province(p)) for p in p24)
    assert len(set(pmap.values())) == len(pmap), "province short names collide"
    for lst, lab in ((l19, "2019"), (l24, "2024")):
        dup = [k for k, v in Counter((r["province_published"], r["name"]) for r in lst).items() if v > 1]
        assert not dup, f"duplicate names in {lab}: {dup}"
        dupn = [k for k, v in Counter(r["name"] for r in lst).items() if v > 1]
        assert not dupn, f"same name under two provinces in {lab}: {dupn}"

    # ---- pairing
    exact, renamed, ambiguous = match_lists(l19, l24)
    assert not ambiguous, f"ambiguous rename candidates (not one-to-one): {ambiguous}"
    rev = {v: k for k, v in exact.items()}
    rev.update({v[0]: k for k, v in renamed.items()})
    basis_of = {v[0]: v[1] for v in renamed.values()}
    by19 = {r["no"]: r for r in l19}

    S19, S24 = "nhc2019_157", "nhc2024_64"
    rows, role_changes = [], []
    for r in l24:
        row = OrderedDict([("name", r["name"]), ("province", pmap[r["province_published"]]),
                           ("role", ROLE[r["role_zh"]]), ("role_zh", r["role_zh"]), ("source_id", S24)])
        o = by19.get(rev.get(r["no"]))
        if o is None:
            row["status"] = "added_2024"
            row["no"] = r["no"]
        else:
            row["status"] = "listed_2019" if o["name"] == r["name"] else "renamed_2024"
            row["no"] = r["no"]
            row["no_2019"] = o["no"]
            if o["name"] != r["name"]:
                row["name_2019"] = o["name"]
                row["matched_by"] = basis_of[r["no"]]
            if ROLE[o["role_zh"]] != ROLE[r["role_zh"]]:
                row["role_zh_2019"] = o["role_zh"]
                role_changes.append(f"{row['province']}: {r['name']} {o['role_zh']}→{r['role_zh']}")
        rows.append(row)
    kept19 = set(rev.values())
    for o in l19:
        if o["no"] in kept19:
            continue
        rows.append(OrderedDict([("name", o["name"]), ("province", pmap[o["province_published"]]),
                                 ("role", ROLE[o["role_zh"]]), ("role_zh", o["role_zh"]),
                                 ("source_id", S19), ("status", "removed_2024"), ("no", o["no"])]))

    st = Counter(r["status"] for r in rows)
    assert st["listed_2019"] + st["renamed_2024"] + st["added_2024"] == 419
    assert st["listed_2019"] + st["renamed_2024"] + st["removed_2024"] == 324
    current = [r for r in rows if r["source_id"] == S24]
    c24 = counts_of(current)
    c19 = counts_of([{"role": ROLE[o["role_zh"]], "province": pmap[o["province_published"]]} for o in l19])
    renamed_list = [f"{r['province']}: {r['name_2019']} → {r['name']} ({r['matched_by']})"
                    for r in rows if r["status"] == "renamed_2024"]
    unusual_roles = [f"{pmap[r['province_published']]} 编号{r['no']} {r['name']}: {r['role_zh']}"
                     for r in l24 if r["role_zh"] == "成员单位"]

    built = utc_now()
    sources = OrderedDict()
    sources[S19] = OrderedDict([
        ("title", "国家卫生健康委办公厅关于建立全国罕见病诊疗协作网的通知"),
        ("document_no", "国卫办医函〔2019〕157号"), ("issuer", "国家卫生健康委办公厅"),
        ("published", "2019-02-12"), ("url", URL_2019_GOVCN), ("origin", URL_2019_NHC),
        ("attachment_url", None), ("retrieved_at", m19["retrieved_at"]), ("sha256", m19["sha256"]),
        ("bytes", m19["bytes"]), ("rows", len(l19)), ("rows_in_file", st["removed_2024"]),
        ("stated_total", 324),
        ("stated_by_role", OrderedDict([("national_lead", 1), ("provincial_lead", 32), ("member", 291)])),
        ("cross_check", xc19),
        ("note", "附件 第一批罕见病诊疗协作网医院名单 is inline on the page as an HTML table (编号/省份/类型/协作网医院). "
                 "Dated 2019-02-12; NHC posted it 2019-02-15, gov.cn 政策库 2019-10-08. All 324 rows are parsed; "
                 "only the hospitals whose name is absent from the 2024 list are rows of this file with this "
                 "source_id. The others are carried by their 2024 row (no_2019 gives their 2019 编号, name_2019 "
                 "the 2019 spelling when it differs)."),
    ])
    if mn:
        sources["nhc2019_157_nhc_copy"] = OrderedDict([
            ("title", "国家卫生健康委办公厅关于建立全国罕见病诊疗协作网的通知"), ("document_no", "国卫办医函〔2019〕157号"),
            ("issuer", "国家卫生健康委办公厅"), ("published", "2019-02-12"), ("url", URL_2019_NHC), ("origin", URL_2019_NHC),
            ("attachment_url", None), ("retrieved_at", mn["retrieved_at"]), ("sha256", mn["sha256"]),
            ("bytes", mn["bytes"]), ("rows", len(tabs[0])), ("rows_in_file", 0), ("stated_total", 324),
            ("cross_check", f"identical to {S19} row by row"),
            ("note", "Original issuing page. nhc.gov.cn answers scripts with HTTP 412, so it was opened in a real "
                     "browser (ego-browser) and the rendered DOM saved; used only as the second copy for the cross-check."),
        ])
    xc24 = (f"Parsed 419 rows, 编号 1..419 continuous. The notice and the docx state no total; the 419 total is "
            f"stated by NHC itself ({'nhc_reply_2024_02403' if mt else 'n/a'}: '目前已有419家医院纳入协作网范围') and by "
            f"the gov.cn news page ({'govcn_news_2024_419' if mg else 'n/a'}): parsed 419 = stated 419. No official "
            f"per-role total exists for 2024; parsed roles are {roles24['national_lead']} national lead / "
            f"{roles24['provincial_lead']} provincial lead / {roles24['member']} member. No second official copy "
            f"of the 2024 list was found (gov.cn 政策库 has no entry for this notice), so no row-level cross-check. "
            f"The same docx bytes were obtained via the browser and via urllib (same sha256). "
            f"Against 2019: {st['listed_2019']} same name, {st['renamed_2024']} renamed on string evidence, "
            f"{st['added_2024']} new names, {st['removed_2024']} 2019 names absent.")
    sources[S24] = OrderedDict([
        ("title", "国家卫生健康委办公厅关于调整全国罕见病诊疗协作网成员医院和办公室人员的通知"),
        ("document_no", "国卫办医政函〔2024〕64号"), ("issuer", "国家卫生健康委办公厅"),
        ("published", "2024-02-28"), ("url", URL_2024_NHC), ("origin", URL_2024_NHC),
        ("attachment_url", URL_2024_DOCX), ("retrieved_at", m24["retrieved_at"]), ("sha256", m24["sha256"]),
        ("bytes", m24["bytes"]), ("rows", len(l24)), ("rows_in_file", len(current)),
        ("stated_total", None),   # neither the notice nor the docx states a number; 419 is stated by the sources below
        ("cross_check", xc24),
        ("note", "Dated 2024-02-28, posted by NHC 医政司 2024-03-18. 附件1 全国罕见病诊疗协作网成员医院名单 is the full "
                 "adjusted list (it replaces the 2019 list; it does not mark what changed). sha256/bytes are of the "
                 "docx, which nhc.gov.cn serves to scripts without the 412 challenge; the notice page itself was read "
                 "in a browser" + (f" (rendered DOM sha256 {mp['sha256']})" if mp else "") + ". 附件2 (协作网办公室人员名单) "
                 "is not used. Two 西藏 rows are typed 成员单位 instead of 成员医院; role=member, role_zh as published."),
    ])
    if mt:
        sources["nhc_reply_2024_02403"] = OrderedDict([
            ("title", "关于政协第十四届全国委员会第二次会议第02403号（医疗卫生类173号）提案答复的函"), ("document_no", None),
            ("issuer", "国家卫生健康委员会"), ("published", "2024-08"), ("url", URL_2024_TIA), ("origin", URL_2024_TIA),
            ("attachment_url", None), ("retrieved_at", mt["retrieved_at"]), ("sha256", mt["sha256"]), ("bytes", mt["bytes"]),
            ("rows", 0), ("rows_in_file", 0), ("stated_total", 419),
            ("cross_check", f"states '2024年…对全国罕见病诊疗协作网成员医院进行调整，目前已有419家医院纳入协作网范围' = {S24} parsed 419"),
            ("note", "Used only for the stated total. Read in a browser (HTTP 412 to scripts); month from the URL path."),
        ])
    if mg:
        sources["govcn_news_2024_419"] = OrderedDict([
            ("title", "诊疗水平再提升！全国罕见病诊疗协作网医院达419家"), ("document_no", None),
            ("issuer", "新华社 (中国政府网 要闻)"), ("published", "2024-10-19"), ("url", URL_2024_GOVNEWS),
            ("origin", URL_2024_GOVNEWS), ("attachment_url", None), ("retrieved_at", mg["retrieved_at"]),
            ("sha256", mg["sha256"]), ("bytes", mg["bytes"]), ("rows", 0), ("rows_in_file", 0), ("stated_total", 419),
            ("cross_check", f"states '全国罕见病诊疗协作网医院总数达到419家' = {S24} parsed 419"),
            ("note", "News item on gov.cn; used only for the stated total."),
        ])
    if m26:
        sources["nhc2026_305"] = OrderedDict([
            ("title", "关于进一步加强罕见病防治工作的通知"), ("document_no", "国卫办医政函〔2026〕305号"),
            ("issuer", "国家卫生健康委办公厅、国家中医药局综合司"), ("published", "2026-09-14"), ("url", URL_2026_NHC),
            ("origin", URL_2026_NHC), ("attachment_url", None), ("retrieved_at", m26["retrieved_at"]),
            ("sha256", m26["sha256"]), ("bytes", m26["bytes"]), ("rows", 0), ("rows_in_file", 0), ("stated_total", None),
            ("cross_check", "page has no table and no attachment: it changes no membership"),
            ("note", "Posted 2026-09-20. Refers to the network as 全国罕见病诊疗协作网络 (协作网络) and directs that "
                     "有条件的中医医院 be brought into it, but lists no hospitals. Read in a browser (HTTP 412 to scripts)."),
        ])

    not_retrieved += [
        {"what": "nhc.gov.cn notice pages by script (2019/157, 2024/64, 2026/305, 2024 proposal reply)",
         "url": "https://www.nhc.gov.cn/", "reason": "HTTP 412 JavaScript challenge to curl/urllib; read in a real browser "
                                                      "instead and the rendered DOM saved (see sources)"},
        {"what": "a second official copy of the 2024 list (国卫办医政函〔2024〕64号 附件1)",
         "url": "https://sousuo.www.gov.cn/search-gov/data?t=zhengcelibrary_bm",
         "reason": "not found: gov.cn 政策库 returns only the 2019/157 and 2020/2 notices for 罕见病诊疗协作网; no provincial "
                   "卫健委 repost with the full attachment turned up in web search"},
        {"what": "provincial implementation lists (各省罕见病诊疗协作网 实施方案/名单)", "url": None,
         "reason": "out of scope for the national list; not retrieved"},
    ]

    search_note = ("gov.cn 政策库 API (sousuo.www.gov.cn/search-gov/data, t=zhengcelibrary_bm and zhengcelibrary, q=罕见病诊疗协作网 "
                   "and q=调整全国罕见病诊疗协作网成员医院)")
    if policy_hits:
        search_note += (f"; NHC 医政司 政策文件 list ({pl_span[2]} most recent documents, {pl_span[0]}..{pl_span[1]}), "
                        f"titles with 罕见/协作网: " + "; ".join(policy_hits))
    search_note += ("; web search (exa, WebSearch) for 罕见病诊疗协作网 成员医院 调整/增补/新增 2020-2026 and for reposts of "
                    "国卫办医政函〔2024〕64号")
    coverage = (
        f"Hospitals of the national 全国罕见病诊疗协作网 as published by the National Health Commission: the first list "
        f"(国卫办医函〔2019〕157号, 324 hospitals) and the adjusted full list that replaced it (国卫办医政函〔2024〕64号 附件1, "
        f"419 hospitals: {roles24['national_lead']} national lead, {roles24['provincial_lead']} provincial lead, "
        f"{roles24['member']} member). Current membership = the {len(current)} rows with source_id {S24}; "
        f"rows with source_id {S19} are 2019 names absent from the 2024 list. As of {built[:10]} the 2024 list is the only "
        f"later official membership revision found; the 2023 第二批罕见病目录 notice carries a disease list only, and "
        f"国卫办医政函〔2026〕305号 (2026-09-14) asks for 有条件的中医医院 to be added but names none. Searched: {search_note}. "
        f"NOT covered: provincial or city-level rare-disease networks, 罕见病门诊 designations required by the 2026 notice, "
        f"国家罕见病医学中心, which diseases or departments each hospital covers, addresses or contacts, Hong Kong/Macao/"
        f"Taiwan (none are listed), and any membership change NHC made after the 2024 list without publishing a list. "
        f"A 2019 name counts as renamed only on string evidence inside the two lists; a hospital that changed its name "
        f"beyond that (for example a medical school upgraded to a university) appears as removed_2024 plus added_2024, "
        f"so removed_2024 means 'this 2019 name is not in the 2024 list', not proof the hospital left the network.")

    out = OrderedDict()
    out["schema"] = 1
    out["title"] = "全国罕见病诊疗协作网 member hospitals"
    prov = OrderedDict()
    prov["built_at"] = built
    prov["builder"] = "tools/china_access/build_hospitals.py"
    prov["sources"] = sources
    c = OrderedDict([("total", c24["total"]), ("national_lead", c24["national_lead"]),
                     ("provincial_lead", c24["provincial_lead"]), ("member", c24["member"]),
                     ("by_province", c24["by_province"]),
                     ("rows_in_file", len(rows)),
                     ("by_status", OrderedDict((k, st[k]) for k in ("listed_2019", "renamed_2024", "added_2024", "removed_2024"))),
                     ("list_2019", c19)])
    prov["counts"] = c
    prov["status_definitions"] = OrderedDict([
        ("listed_2019", "in the 2019 list and in the 2024 list under the same name (no_2019 = 2019 编号)"),
        ("renamed_2024", "in the 2024 list under a name that differs from its 2019 entry; name_2019 is the 2019 spelling "
                         "and matched_by the string evidence used (no outside knowledge)"),
        ("added_2024", "in the 2024 list; its name is not in the 2019 list of the same province (may include renamings "
                       "the string rules cannot see)"),
        ("removed_2024", "in the 2019 list; this name is not in the 2024 list (left the network, or renamed beyond what "
                         "the string rules can see)"),
    ])
    prov["fields"] = ("no = 编号 in the source table of source_id; role_zh_2019 appears only where the 2019 role differs. "
                      "province is the short form; provinces maps it to the spelling published in both lists.")
    prov["provinces"] = OrderedDict((v, k) for k, v in pmap.items())
    prov["changes_2019_to_2024"] = OrderedDict([
        ("role_changes", role_changes), ("renamed_on_string_evidence", renamed_list),
        ("role_zh_variants", unusual_roles)])
    prov["coverage"] = coverage
    prov["not_retrieved"] = not_retrieved
    out["provenance"] = prov
    out["hospitals"] = rows

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    tmp = args.out + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
        fh.write("\n")
    os.replace(tmp, args.out)
    print(f"wrote {args.out}: {len(rows)} rows; current {len(current)} "
          f"({c24['national_lead']}/{c24['provincial_lead']}/{c24['member']}); status {dict(st)}")
    print("renamed:", *renamed_list, sep="\n  ")
    print("role changes:", *role_changes, sep="\n  ")
    return 0


if __name__ == "__main__":
    sys.exit(main())
