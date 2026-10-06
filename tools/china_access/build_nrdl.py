#!/usr/bin/env python3
"""Build zebra/data/china_nrdl.json: the current 国家医保药品目录 (NRDL), Western-medicine part.

Edition: 《国家基本医疗保险、生育保险和工伤保险药品目录（2025年）》, 医保发〔2025〕33号
(国家医保局 人力资源社会保障部, 成文 2025-12-05, published 2025-12-07, in force from 2026-01-01; it
repealed the 2024 edition). The same notice issued the first 《商业健康保险创新药品目录（2025年）》, which is
bundled as a separate list: it is NOT 医保 reimbursement (the notice: "商保创新药目录内药品医保基金不予支付").

What is parsed (every row comes from the official PDF attachment, nothing is typed in by hand):
  * 西药部分 (all rows, including the "★(n)" repeat rows for extra dosage forms);
  * 协议期内谈判药品部分 (一)西药 and the 竞价药品部分 printed after it (the 凡例 counts 竞价 inside the 411
    Western 谈判 drugs);
  * the 凡例 clauses (verbatim) and the two 凡例 lists that say what a "◇" class row covers;
  * 商业健康保险创新药品目录 (all 19 rows).
  NOT parsed: 中成药部分, 协议期内谈判药品部分 (二)中成药, 中药饮片部分.

Checks (the build aborts if one fails): numbering is continuous; parsed counts equal the totals the 凡例
states (西药部分 1446, 西药甲类 393, 谈判 472 = 西药 411 + 中成药 61) and the totals in NHSA's statistics sheet
(西药 1857 = 1446 + 411; 协议期内谈判品种 472); the gov.cn and nhsa.gov.cn copies of each attachment are
byte-identical; a second, independent table extractor (PyMuPDF) yields the same cell text for every row.

Rerun (from ~/zebra-mod; needs python3.11 with pdfplumber and PyMuPDF; build tooling only, not used by zebra):
  # fresh download into a NEW empty directory, then parse:
  python3 -I tools/china_access/build_nrdl.py --download --dl-dir /tmp/nrdl-dl
  # re-parse files already downloaded by a previous --download run (uses its manifest.json):
  python3 -I tools/china_access/build_nrdl.py --dl-dir /tmp/nrdl-dl
  Output: zebra/data/china_nrdl.json (override with --out).
"""
import argparse
import datetime as dt
import hashlib
import html
import json
import os
import re
import sys
import time
import urllib.request

import fitz  # PyMuPDF, used only for the independent cross-check extraction
import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DEFAULT_OUT = os.path.join(REPO, 'zebra', 'data', 'china_nrdl.json')
UA = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126 Safari/537.36'

NOTICE_GOVCN = 'https://www.gov.cn/zhengce/zhengceku/202512/content_7050533.htm'
NOTICE_NHSA = 'https://www.nhsa.gov.cn/art/2025/12/7/art_104_18970.html'
STATS_PAGE = 'https://www.nhsa.gov.cn/art/2026/8/5/art_257_21672.html'
DOWNLOADS = [  # (key, url, referer, filename)
    ('notice_govcn', NOTICE_GOVCN, None, 'notice_govcn.html'),
    ('nrdl_pdf_govcn', 'https://www.gov.cn/zhengce/zhengceku/202512/P020251207694651387946.pdf', NOTICE_GOVCN,
     'nrdl2025_govcn.pdf'),
    ('cidl_pdf_govcn', 'https://www.gov.cn/zhengce/zhengceku/202512/P020251207694651942582.pdf', NOTICE_GOVCN,
     'cidl2025_govcn.pdf'),
    ('notice_nhsa', NOTICE_NHSA, None, 'notice_nhsa.html'),
    ('nrdl_pdf_nhsa', 'https://www.nhsa.gov.cn/module/download/downfile.jsp?classid=0&filename='
     'a32f9f2f3fc046afaf08471f87456ce3.pdf', NOTICE_NHSA, 'nrdl2025_nhsa.pdf'),
    ('cidl_pdf_nhsa', 'https://www.nhsa.gov.cn/module/download/downfile.jsp?classid=0&filename='
     'da298734646c4a4ca1b5ca612ee23fb1.pdf', NOTICE_NHSA, 'cidl2025_nhsa.pdf'),
    ('stats_page', STATS_PAGE, None, 'stats_page.html'),
    ('stats_pdf', 'https://www.nhsa.gov.cn/module/download/downfile.jsp?classid=0&filename='
     '49e01eeefa02436aa748c9f174fdbb7c.pdf', STATS_PAGE, 'stats_2000_2025.pdf'),
]


def now_utc():
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace('+00:00', 'Z')


def sha256_file(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def download(dl_dir):
    if os.path.isdir(dl_dir) and os.listdir(dl_dir):
        raise SystemExit(f'--download needs a new empty directory: {dl_dir} is not empty')
    os.makedirs(dl_dir, exist_ok=True)
    manifest, last_hit = {}, {}
    for key, url, referer, fn in DOWNLOADS:
        host = url.split('/')[2]
        wait = 1.5 - (time.time() - last_hit.get(host, 0))
        if wait > 0:
            time.sleep(wait)
        req = urllib.request.Request(url, headers={'User-Agent': UA, **({'Referer': referer} if referer else {})})
        with urllib.request.urlopen(req, timeout=180) as r:
            data, status, ctype, final = r.read(), r.status, r.headers.get('Content-Type'), r.geturl()
        last_hit[host] = time.time()
        if status != 200 or not data:
            raise SystemExit(f'download failed: {url} HTTP {status}')
        if fn.endswith('.pdf') and not data.startswith(b'%PDF'):
            raise SystemExit(f'not a PDF: {url} ({ctype})')
        path = os.path.join(dl_dir, fn)
        open(path, 'wb').write(data)
        manifest[key] = {'url': url, 'final_url': final, 'file': fn, 'http_status': status, 'content_type': ctype,
                         'retrieved_at': now_utc(), 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
        print(f'  {key}: HTTP {status} {len(data)} bytes', file=sys.stderr)
    json.dump(manifest, open(os.path.join(dl_dir, 'manifest.json'), 'w'), ensure_ascii=False, indent=1)
    return manifest


def load_manifest(dl_dir):
    manifest = json.load(open(os.path.join(dl_dir, 'manifest.json'), encoding='utf-8'))
    for key, m in manifest.items():  # the files must still be the bytes that were downloaded
        assert sha256_file(os.path.join(dl_dir, m['file'])) == m['sha256'], f'{key}: file changed since download'
    return manifest


def html_text(path):
    t = open(path, encoding='utf-8', errors='replace').read()
    t = re.sub(r'<script.*?</script>|<style.*?</style>|<!--.*?-->', '', t, flags=re.S)
    s = html.unescape(re.sub(r'<[^>]+>', '\n', t))
    s = re.sub(r'[ \t　\xa0]+', ' ', s)
    return re.sub(r'\n\s*\n+', '\n', s), t


def parse_notice(path):
    """Metadata of the notice; asserts that it is the document and attachments we think it is."""
    s, raw = html_text(path)
    flat = re.sub(r'\s+', '', s)
    meta = {'document_no': re.search(r'医保发〔\d{4}〕\d+号', flat).group(0)}
    assert meta['document_no'] == '医保发〔2025〕33号', meta
    m = re.search(r'国家医保局人力资源社会保障部关于印发《国家基本医疗保险、生育保险和工伤保险药品目录》以及'
                  r'《商业健康保险创新药品目录》（2025年）的通知', flat)
    assert m, 'notice title not found'
    meta['title'] = ('国家医保局 人力资源社会保障部关于印发《国家基本医疗保险、生育保险和工伤保险药品目录》以及'
                     '《商业健康保险创新药品目录》（2025年）的通知')
    assert '新版药品目录自2026年1月1日起正式执行' in flat
    assert '《国家基本医疗保险、工伤保险和生育保险药品目录（2024年）》（医保发〔2024〕33号）同时废止' in flat
    assert '商保创新药目录内药品医保基金不予支付' in flat
    meta['attachments'] = [(html.unescape(u), re.sub(r'\s+', '', html.unescape(n)))
                           for u, n in re.findall(r'href="([^"]+)"[^>]*>([^<]*)</a>', raw)
                           if re.search(r'\.pdf|downfile', u) and '目录' in n]
    m = re.search(r'成文日期：(\d{4})年(\d{2})月(\d{2})日', flat)
    if m:
        meta['signed'] = '-'.join(m.groups())
    m = re.search(r'发布日期：(\d{4}-\d{2}-\d{2})', flat)
    if m:
        meta['published'] = m.group(1)
    return meta


# ---------------------------------------------------------------- cell text reconstruction
ASCII_END = re.compile(r'[A-Za-z0-9.,;:%)\]]$')
ASCII_START = re.compile(r'^[A-Za-z(\[]')


def cell_lines(page, bbox):
    sub = page.within_bbox(bbox)
    return sub.extract_text_lines(strip=True, return_chars=True)


OPENERS = '(（[【〔“‘<《'
CLOSERS = ')）]】〕”’>》，。、；：？！,.;:%'
CJK = re.compile(r'[\u3400-\u9fff\uf900-\ufaff\u3000-\u303f\uff00-\uffef]')


def first_unit_width(line):
    """Width of the first unbreakable unit of a line: leading opening brackets + one CJK character or a
    whole run of non-CJK characters (a Latin word wraps as a whole) + closing punctuation that may not start
    a line. A previous line ended early if this unit would not have fitted in the room left on it."""
    chs = [c for c in line['chars'] if c['text'].strip()]
    if not chs:
        return 0.0
    i = 0
    while i < len(chs) - 1 and chs[i]['text'] in OPENERS:
        i += 1
    if CJK.match(chs[i]['text']) and chs[i]['text'] not in OPENERS + CLOSERS:
        i += 1
    else:
        while i < len(chs) and not CJK.match(chs[i]['text']) and (
                i == 0 or chs[i]['x0'] - chs[i - 1]['x1'] < 1.0):
            i += 1
    if 0 < i < len(chs) and chs[i - 1]['text'] in OPENERS:  # '11(' cannot end a line: take what follows
        i += 1
    while i < len(chs) and chs[i]['text'] in CLOSERS:
        i += 1
    i = max(i, 1)
    return chs[min(i, len(chs)) - 1]['x1'] - chs[0]['x0']


def join_cell(page, bbox, breaks_log, where, split_items):
    """Text of one table cell. Lines are soft wraps of one text and are joined with nothing (Chinese) or one
    space (between two Latin words). Only when split_items is set (the drug-name column), a line that ended
    although the next line's first unit would have fitted (several names stacked under one number) starts a
    new item; the items are returned as a list."""
    lines = cell_lines(page, bbox)
    if not lines:
        return [], ''
    right = bbox[2] - 2.0  # the cell's text area ends ~2 pt inside the ruling line
    items, cur = [], lines[0]['text']
    raw = ''.join(l['text'] for l in lines)
    for a, b in zip(lines, lines[1:]):
        room = right - a['x1']
        unit = first_unit_width(b)
        latin_pair = bool(ASCII_END.search(a['text']) and ASCII_START.search(b['text']))
        if latin_pair:
            unit += 3.0  # the space between two Latin words was dropped at the wrap
        hard = room >= unit + 2.0 and b['text'][0] not in OPENERS + CLOSERS
        if hard:
            breaks_log.append(('hard' if split_items else 'hard-joined', where, a['text'], b['text'],
                               round(room, 1), round(unit, 1)))
        if hard and split_items:
            items.append(cur)
            cur = b['text']
            continue
        sep = ' ' if latin_pair else ''
        if sep:
            breaks_log.append(('space-join', where, a['text'][-14:], b['text'][:14]))
        cur += sep + b['text']
    items.append(cur)
    return items, raw


# ---------------------------------------------------------------- PDF table parsing
HEADINGS = {'西药部分': 'west', '中成药部分': 'tcm', '协议期内谈判药品部分': 'neg', '(一)西药': 'neg_west',
            '(二)中成药': 'neg_tcm', '竞价药品部分': 'bid', '中药饮片部分': 'herb',
            '商业健康保险创新药品目录': 'cidl'}
NO_RE = re.compile(r'^(\d+)$|^★[(（](\d+)[)）]$')
CODE_RE = re.compile(r'^[XZ][A-Z]\d{2}[A-Z]{0,2}$|^[XZ][A-Z]$')


def norm_heading(t):
    return re.sub(r'\s+', '', t).replace('（', '(').replace('）', ')')


def page_headings(page):
    out = []
    for ln in page.extract_text_lines(strip=True):
        key = HEADINGS.get(norm_heading(ln['text']))
        if key:
            out.append((ln['top'], key))
    return out


def parse_tables(pdf_path, layout, wanted):
    """Walk every page, follow section headings, and return parsed rows of the wanted sections.
    layout: {section: [column names in order]} for the sections whose tables are parsed."""
    rows, cats, breaks, skipped, footnotes, sec_pages = [], {}, [], [], [], {}
    pdf = pdfplumber.open(pdf_path)
    section, cat = None, None
    for pno, page in enumerate(pdf.pages, start=1):
        heads = page_headings(page)
        tables = sorted(page.find_tables(), key=lambda t: t.bbox[1])
        for tb in tables:
            for top, key in heads:
                # '协议期内谈判药品部分' is always followed by its own '(一)西药' / '(二)中成药' line
                if top < tb.bbox[1] and key != 'neg':
                    section = key
            sec_pages.setdefault(section, []).append(pno)
            if section not in wanted:
                continue
            cols = layout[section]
            for ri, row in enumerate(tb.rows):
                cells = row.cells
                if len(cells) != len(cols):
                    # a different (small) table, e.g. the tiered-price table under 备注2: keep it verbatim
                    skipped.append({'pdf_page': pno, 'section': section, 'ncols': len(cells),
                                    'rows': [[re.sub(r'\s+', ' ', c or '').strip() for c in r]
                                             for r in tb.extract()]})
                    break
                txt, raw = {}, {}
                for name, bbox in zip(cols, cells):
                    if bbox is None:
                        txt[name], raw[name] = None, ''
                        continue
                    items, r = join_cell(page, bbox, breaks, (pno, ri, name), name == 'name')
                    txt[name], raw[name] = items, r
                flat = {k: (''.join(v) if v else '') for k, v in txt.items()}
                if flat.get('no') == '编号' or flat.get('code', '').startswith('药品分'):
                    continue  # header row
                # vertical extent of the row: from its code cell (category row) or number cell (drug row);
                # the empty category columns on the left can be merged over several rows
                key_cell = cells[0] if flat['code'] else cells[cols.index('no')]
                assert key_cell, (pno, ri, flat)
                top, bottom = key_cell[1], key_cell[3]
                band = re.sub(r'\s+', '', page.within_bbox((tb.bbox[0], top, tb.bbox[2], bottom)).extract_text())
                code = flat['code']
                if code:
                    assert CODE_RE.match(code), (pno, code)
                    # merged category-name cells are not always recovered as cells: read the whole band
                    assert band.startswith(code), (pno, ri, band)
                    name = band[len(code):]
                    assert name, (pno, ri, 'category row without a name', flat)
                    cats.setdefault(code, {}).setdefault(name, []).append(section)
                    cat = code
                    assert not flat.get('no'), (pno, 'number on a category row', flat)
                    continue
                # every printed character of the row must have landed in some cell
                got = sum(len(re.sub(r'\s+', '', r)) for r in raw.values())
                assert got == len(band), (pno, ri, 'text outside the cells', band, raw)
                no = flat.get('no', '')
                m = NO_RE.match(no)
                assert m, (pno, ri, 'row without code or valid number', flat)
                rows.append({'section': section, 'pdf_page': pno, 'row': ri, 'category_code': cat,
                             'primary_no': int(m.group(1) or m.group(2)), 'items': txt, 'raw': raw})
        # footnotes printed under the 谈判/竞价 tables ("备注1：…")
        if section in ('neg_west', 'neg_tcm', 'bid'):
            for ln in page.extract_text_lines(strip=True):
                if re.match(r'^备注\d?：', ln['text']) and not any(
                        t.bbox[1] <= ln['top'] <= t.bbox[3] for t in tables):
                    footnotes.append({'pdf_page': pno, 'section': section, 'text': ln['text']})
    return rows, cats, breaks, skipped, footnotes, sec_pages


def fitz_rows(pdf_path, pages, ncols):
    """Independent extraction with PyMuPDF: per page, the drug rows' cell texts (whitespace removed)."""
    doc = fitz.open(pdf_path)
    out = {}
    for pno in sorted(set(pages)):
        page_rows = []
        for tb in doc[pno - 1].find_tables().tables:
            for r in tb.extract():
                if len(r) != ncols:
                    break
                cells = [re.sub(r'\s+', '', c or '') for c in r]
                if NO_RE.match(cells[6]):
                    page_rows.append(cells)
        out[pno] = page_rows
    return out


# ---------------------------------------------------------------- 凡例 (general rules)
def parse_fanli(pdf_path):
    """凡例 (printed pages 1-7): stated totals, the 17 clauses verbatim, and the three tables in it."""
    pdf = pdfplumber.open(pdf_path)
    pages, tables = [], []
    for p in pdf.pages[1:8]:  # PDF pages 2-8 = printed pages 1-7
        tbs = p.find_tables()
        tables += [(p.page_number, tb) for tb in tbs]
        inside = lambda o, tbs=tbs: any(  # noqa: E731
            tb.bbox[0] <= (o['x0'] + o['x1']) / 2 <= tb.bbox[2] and tb.bbox[1] <= (o['top'] + o['bottom']) / 2
            <= tb.bbox[3] for tb in tbs)
        prose = p.filter(lambda o: o.get('object_type') != 'char' or not inside(o))
        lines = [l['text'] for l in prose.extract_text_lines(strip=True)
                 if not re.fullmatch(r'第\s*\d+\s*页', l['text'])]
        pages.append('\n'.join(lines))
    flat = re.sub(r'\s+', '', '\n'.join(pages))  # line breaks here are only page-width wraps
    stated = {
        'west_part': int(re.search(r'其中西药部分(\d+)个', flat).group(1)),
        'tcm_part': int(re.search(r'中成药部分(\d+)个（含民族药', flat).group(1)),
        'negotiated_total': int(re.search(r'协议期内谈判药品部分(\d+)个（含西药', flat).group(1)),
        'negotiated_west': int(re.search(r'（含西药(\d+)个、中成药', flat).group(1)),
        'negotiated_tcm': int(re.search(r'个、中成药(\d+)个），共计', flat).group(1)),
        'grand_total': int(re.search(r'共计(\d+)个', flat).group(1)),
        'west_class_jia': int(re.search(r'西药甲类药品(\d+)个', flat).group(1)),
    }
    clauses = {}
    marks = list(re.finditer(r'[（(]([一二三四五六七八九十]+)[）)]', flat))
    for i, m in enumerate(marks):
        body = flat[m.end():marks[i + 1].start() if i + 1 < len(marks) else len(flat)]
        body = re.sub(r'[一二三四五]、(目录构成|编排与分类|名称与剂型|限定支付范围|其他)$', '', body)
        clauses['（%s）' % m.group(1)] = body
    assert list(clauses)[0] == '（一）' and list(clauses)[-1] == '（十七）' and len(clauses) == 17, list(clauses)
    # table under (八): merged dosage-form groups; table under (十)4: the 47 names covered by 西药 1303
    forms, otc_1303 = {}, {}
    for pno, tb in tables:
        for r in tb.extract():
            r = [re.sub(r'\s+', '', c or '') for c in r]
            if len(r) == 2 and r[0] != '合并归类的剂型':
                forms[r[0]] = r[1]
            elif len(r) == 4 and r[0] != '序号':
                for n, name in ((r[0], r[1]), (r[2], r[3])):
                    if n:
                        otc_1303[int(n)] = name
    assert sorted(otc_1303) == list(range(1, len(otc_1303) + 1)), 'list under 1303 not continuous'
    assert '口服常释剂型' in forms and '注射剂' in forms, forms
    m = re.search(r'西药部分第207号“缓解消化道不适症状的复方OTC制剂”包括：([^。]+)。', flat)
    list_207 = m.group(1).split('、')
    ranges = {}
    for num, a, b in re.findall(r'[（(](十[四五])[）)]参保人员使用西药部分第(\d+)-(\d+)号', flat):
        ranges['（%s）' % num] = (int(a), int(b))
    assert len(ranges) == 2, ranges
    return stated, clauses, forms, list_207, [otc_1303[k] for k in sorted(otc_1303)], ranges


def parse_stats(pdf_path):
    text = pdfplumber.open(pdf_path).pages[0].extract_text()
    lines = text.split('\n')
    hdr = next(l for l in lines if l.startswith('版本'))
    assert hdr.split()[-1] == '2025', hdr

    def last(prefix):
        ln = next(l for l in lines if l.startswith(prefix))
        return int(ln.split()[-1])
    return {'total': last('药品数量'), 'west': last('西药 '), 'tcm': last('中成药 '),
            'negotiated': last('其中：协议期内谈判品种'), 'newly_negotiated_2025': last('新增谈判准入药品'),
            'antineoplastic_west': last('抗肿瘤药物（西药）'), 'rare_disease': last('罕见病治疗药物')}


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--dl-dir', required=True, help='directory holding (or receiving) the downloaded files')
    ap.add_argument('--download', action='store_true', help='download the sources into --dl-dir first')
    ap.add_argument('--out', default=DEFAULT_OUT)
    a = ap.parse_args()
    t0 = time.time()
    manifest = download(a.dl_dir) if a.download else load_manifest(a.dl_dir)
    P = lambda k: os.path.join(a.dl_dir, manifest[k]['file'])  # noqa: E731

    # 1. the notice (both official copies) and its attachments
    n_gov, n_nhsa = parse_notice(P('notice_govcn')), parse_notice(P('notice_nhsa'))
    for n in (n_gov, n_nhsa):
        names = [x[1] for x in n['attachments']]
        assert '国家基本医疗保险、生育保险和工伤保险药品目录（2025年）' in names, names
        assert '商业健康保险创新药品目录（2025年）' in names, names
        assert not [u for u, _ in n['attachments'] if re.search(r'\.(xlsx?|et|csv)\b', u, re.I)]
    for k in ('nrdl', 'cidl'):
        assert manifest[f'{k}_pdf_govcn']['sha256'] == manifest[f'{k}_pdf_nhsa']['sha256'], \
            f'{k}: gov.cn and nhsa.gov.cn copies differ'
    nrdl_pdf, cidl_pdf = P('nrdl_pdf_govcn'), P('cidl_pdf_govcn')

    # 2. 凡例 and the stated totals
    stated, clauses, fanli_forms, list_207, list_1303, fanli_ranges = parse_fanli(nrdl_pdf)
    stats = parse_stats(P('stats_pdf'))

    # 3. the tables
    WEST = ['code', 'c1', 'c2', 'c3', 'c4', 'class', 'no', 'name', 'form', 'remark']
    NEG = ['code', 'c1', 'c2', 'c3', 'c4', 'class', 'no', 'name', 'price', 'remark', 'period']
    rows, cats, breaks, skipped, footnotes, sec_pages = parse_tables(
        nrdl_pdf, {'west': WEST, 'neg_west': NEG, 'neg_tcm': NEG, 'bid': NEG}, {'west', 'neg_west', 'neg_tcm', 'bid'})
    by_sec = {s: [r for r in rows if r['section'] == s] for s in ('west', 'neg_west', 'neg_tcm', 'bid')}
    print('section pages:', {k: (min(v), max(v)) for k, v in sec_pages.items() if k}, file=sys.stderr)

    # hard breaks: only allowed in the name column (several names stacked under one number)
    hard = [b for b in breaks if b[0] == 'hard']
    bad = [b for b in hard if b[1][2] != 'name']
    assert not bad, ('unexpected hard line breaks', bad[:10])

    # 4. counts against the stated totals
    w = by_sec['west']
    prim = [r for r in w if not r['items']['no'][0].startswith('★')]
    assert [r['primary_no'] for r in prim] == list(range(1, len(prim) + 1)), '西药部分 numbering not continuous'
    assert len(prim) == stated['west_part'] == 1446, (len(prim), stated)
    assert all(r['primary_no'] <= len(prim) for r in w)
    # 甲/乙 is printed per row (dosage form); the 凡例 total counts numbers with at least one 甲 row
    jia = len({r['primary_no'] for r in w if r['items']['class'][0] == '甲'})
    assert jia == stated['west_class_jia'], ('甲类', jia, stated['west_class_jia'])
    for s in ('neg_west', 'neg_tcm', 'bid'):
        assert [r['primary_no'] for r in by_sec[s]] == list(range(1, len(by_sec[s]) + 1)), s
    nw, nt, nb = len(by_sec['neg_west']), len(by_sec['neg_tcm']), len(by_sec['bid'])
    assert nw + nb == stated['negotiated_west'] == 411, (nw, nb, stated)
    assert nt == stated['negotiated_tcm'] and nw + nb + nt == stated['negotiated_total'] == 472
    assert stats['west'] == len(prim) + nw + nb == 1857, stats
    assert stats['negotiated'] == nw + nb + nt and stats['total'] == stated['grand_total']

    # 5. independent re-extraction (PyMuPDF) must give the same cell text, row by row
    mism = []
    for s, cols, ncol in (('west', WEST, 10), ('neg_west', NEG, 11), ('bid', NEG, 11)):
        fz = fitz_rows(nrdl_pdf, [r['pdf_page'] for r in by_sec[s]], ncol)
        for pno in sorted(set(r['pdf_page'] for r in by_sec[s])):
            mine = [r for r in by_sec[s] if r['pdf_page'] == pno]
            theirs = fz[pno]
            if len(mine) != len(theirs):
                mism.append((s, pno, 'row count', len(mine), len(theirs)))
                continue
            for r, f in zip(mine, theirs):
                for ci, c in enumerate(cols):
                    if c.startswith('c') or c == 'code':
                        continue
                    if re.sub(r'\s+', '', r['raw'][c]) != f[ci]:
                        mism.append((s, pno, r['items']['no'], c, r['raw'][c], f[ci]))
    assert not mism, ('pdfplumber vs PyMuPDF differ', mism[:10])

    # 6. assemble rows
    def one(items):
        return None if not items else '、'.join(items)

    def single(items, what, r):
        if not items:
            return None
        assert len(items) == 1, (what, r['pdf_page'], items)
        return items[0]

    SEC_LABEL = {'west': '西药部分', 'neg_west': '协议期内谈判药品部分（一）西药', 'bid': '竞价药品部分'}
    PERIOD_RE = re.compile(r'^\d{4}年\d{1,2}月\d{1,2}日至\d{4}年\d{1,2}月\d{1,2}日$')
    drugs, stacked = [], []
    for s in ('west', 'neg_west', 'bid'):
        for r in by_sec[s]:
            it = r['items']
            if len(it['name']) > 1:
                stacked.append({'section': SEC_LABEL[s], 'no': it['no'][0], 'names': it['name']})
            d = {'no': single(it['no'], 'no', r), 'name_zh': one(it['name']),
                 'dosage_form': single(it.get('form'), 'form', r) if s == 'west' else None,
                 'section': SEC_LABEL[s], 'class': single(it['class'], 'class', r),
                 'restriction': single(it['remark'], 'remark', r),
                 'agreement_period': single(it.get('period'), 'period', r) if s != 'west' else None}
            if s != 'west':
                d['payment_standard'] = single(it['price'], 'price', r)
                assert d['agreement_period'] and PERIOD_RE.match(d['agreement_period']), (r['pdf_page'], d)
            assert d['class'] in ('甲', '乙') and d['name_zh'], d
            assert not (d['restriction'] and '日至' in d['restriction']), d
            d['primary_no'] = r['primary_no']
            d['category_code'] = r['category_code']
            d['pdf_page'] = r['pdf_page']
            if s == 'west' and d['primary_no'] == 207:
                d['includes'] = list_207
            if s == 'west' and d['primary_no'] == 1303:
                d['includes'] = list_1303
            for clause, (lo, hi) in fanli_ranges.items():
                if s == 'west' and lo <= d['primary_no'] <= hi:
                    d['see_fanli'] = clause
            m = re.search(r'见(备注\d)', d.get('payment_standard') or '')
            if m:  # e.g. 治疗用碘[131I]化钠胶囊: the price is a tiered scheme printed under the table
                fn = [f['text'] for f in footnotes if f['text'].startswith(m.group(1))]
                tbl = [t['rows'] for t in skipped if t['pdf_page'] == r['pdf_page']]
                assert len(fn) == 1 and len(tbl) == 1, (m.group(1), fn, tbl)
                d['payment_standard_note'] = {'text': fn[0], 'table': tbl[0]}
            d['source_id'] = 'nrdl2025'
            drugs.append(d)
    neg_class = {d['class'] for d in drugs if d['section'] != '西药部分'}
    assert neg_class == {'乙'}, neg_class
    assert [t['pdf_page'] for t in skipped] == [d['pdf_page'] for d in drugs if 'payment_standard_note' in d], skipped

    # extra stated figure (reported, not asserted: its definition is the publisher's)
    xl01 = sorted({(d['section'] if d['section'] == '西药部分' else 'neg', d['primary_no']) for d in drugs
                   if (d['category_code'] or '').startswith('XL01')})

    # 7. 商业健康保险创新药品目录
    CID = ['code', 'c1', 'c2', 'c3', 'c4', 'no', 'name', 'brand', 'indication', 'mah', 'licensee', 'period']
    crow, ccats, cbreaks, cskipped, _, _ = parse_tables(cidl_pdf, {'cidl': CID}, {'cidl'})
    assert not cskipped and not [b for b in cbreaks if b[0] == 'hard'], ([b for b in cbreaks if b[0] == 'hard'],)
    assert [r['primary_no'] for r in crow] == list(range(1, len(crow) + 1))
    cidl = []
    for r in crow:
        it = r['items']
        e = {'no': it['no'][0], 'name_zh': it['name'][0], 'brand_name_zh': one(it['brand']),
             'indication': one(it['indication']), 'marketing_authorization_holder': one(it['mah']),
             'authorized_company': one(it['licensee']), 'validity': one(it['period']),
             'category_code': r['category_code'], 'pdf_page': r['pdf_page'], 'source_id': 'cidl2025'}
        assert PERIOD_RE.match(e['validity']), e
        cidl.append(e)
    for k, v in ccats.items():
        for name, secs in v.items():
            cats.setdefault(k, {}).setdefault(name, []).extend(secs)
    SEC_ALL = {**SEC_LABEL, 'neg_tcm': '协议期内谈判药品部分（二）中成药', 'cidl': '商业健康保险创新药品目录'}
    categories, cat_variants = {}, {}
    for k in sorted(cats):
        names = sorted(cats[k], key=lambda n: (0 if 'west' in cats[k][n] else 1, n))
        categories[k] = names[0]
        if len(names) > 1:  # the same code printed with slightly different wording in different tables
            cat_variants[k] = {n: sorted({SEC_ALL[x] for x in cats[k][n]}) for n in names}

    # 8. provenance + write
    M = manifest
    counts = {
        'west_part_rows': len(w), 'west_part_numbers': len(prim), 'west_part_class_jia_numbers': jia,
        'negotiated_west': nw, 'bidding_west': nb, 'negotiated_and_bidding_west': nw + nb,
        'west_numbers_total': len(prim) + nw + nb, 'drugs_rows_total': len(drugs),
        'rows_with_restriction': sum(1 for d in drugs if d['restriction']),
        'negotiated_tcm_not_bundled': nt, 'commercial_innovative_list': len(cidl),
        'xl01_antineoplastic_numbers_parsed': len(xl01),
        'xl01_antineoplastic_stated_by_nhsa_stats': stats['antineoplastic_west'],
        'xl01_note': ('reported, not asserted: unique numbers filed under XL01 in 西药部分 + 谈判 + 竞价; NHSA\'s '
                      'figure follows its own definition (注释5) and differs by %d' % (len(xl01) - stats['antineoplastic_west'])),
    }
    src_nrdl = {
        'title': '国家基本医疗保险、生育保险和工伤保险药品目录（2025年）',
        'document_no': n_gov['document_no'], 'issuer': '国家医保局 人力资源社会保障部',
        'published': n_nhsa.get('published', '2025-12-07'), 'signed': n_gov.get('signed'),
        'url': M['notice_govcn']['url'], 'origin': M['notice_nhsa']['url'],
        'attachment_url': M['nrdl_pdf_govcn']['url'], 'attachment_url_origin': M['nrdl_pdf_nhsa']['url'],
        'retrieved_at': M['nrdl_pdf_govcn']['retrieved_at'], 'sha256': M['nrdl_pdf_govcn']['sha256'],
        'bytes': M['nrdl_pdf_govcn']['bytes'], 'rows': len(drugs),
        'stated_total': stated['west_part'] + stated['negotiated_west'],
        'stated_totals_fanli': stated,
        'cross_check': (
            f"The gov.cn and nhsa.gov.cn copies of the PDF are byte-identical (sha256 {M['nrdl_pdf_govcn']['sha256']}); "
            "no Excel/other copy is attached to either notice. Parsed with pdfplumber and re-extracted independently "
            "with PyMuPDF: same rows and identical cell text (whitespace ignored) for all "
            f"{len(drugs)} rows. Counts equal the 凡例: 西药部分 {len(prim)} numbers (stated {stated['west_part']}), "
            f"甲类 {jia} (stated {stated['west_class_jia']}), 谈判 西药 {nw} + 竞价 {nb} = {nw + nb} "
            f"(stated 含西药 {stated['negotiated_west']}), 谈判 中成药 {nt} (stated {stated['negotiated_tcm']}); "
            f"and NHSA's statistics sheet: 西药 {len(prim) + nw + nb} (stated {stats['west']}), "
            f"协议期内谈判品种 {nw + nb + nt} (stated {stats['negotiated']})."),
        'note': ('Rows are as printed. "★(n)" rows repeat number n for another dosage form; class can differ from '
                 'the primary row, and the 凡例 甲类 total (393) counts numbers with at least one 甲 row. dosage_form '
                 'is null where the 目录 does not list it separately: every 谈判/竞价 row, and 西药部分 rows that '
                 'carry the full approved name (凡例（七）: 2025 additions and drugs transferred from the 谈判 part; '
                 'so a null dosage_form does not by itself mean "new in 2025"). Line breaks inside a cell are not '
                 'kept (they are wraps; joined with nothing, or one space between two Latin words). In the '
                 '谈判/竞价 tables "*" in 医保支付标准 means the company asked for the price to be kept confidential '
                 '(备注1). Several names printed one per line under a single number are joined with "、". '
                 '"no" is as printed (including full-width "★（406）") and restarts in each section, so (section, no) '
                 'is the key; primary_no is its number as an integer. pdf_page is the 1-based page of the PDF file; '
                 'the printed "第N页" is pdf_page - 1 (page 1 is the unnumbered contents page).'),
        'footnotes': [f['text'] for f in footnotes if f['section'] in ('neg_west', 'bid')],
        'stacked_name_rows': stacked,
    }
    src_cidl = {
        'title': '商业健康保险创新药品目录（2025年）', 'document_no': n_gov['document_no'],
        'issuer': '国家医保局 人力资源社会保障部', 'published': src_nrdl['published'],
        'url': M['notice_govcn']['url'], 'origin': M['notice_nhsa']['url'],
        'attachment_url': M['cidl_pdf_govcn']['url'], 'attachment_url_origin': M['cidl_pdf_nhsa']['url'],
        'retrieved_at': M['cidl_pdf_govcn']['retrieved_at'], 'sha256': M['cidl_pdf_govcn']['sha256'],
        'bytes': M['cidl_pdf_govcn']['bytes'], 'rows': len(cidl), 'stated_total': None,
        'cross_check': ('gov.cn and nhsa.gov.cn copies byte-identical; numbering 1..%d continuous; no stated '
                        'total in the notice or the list itself' % len(cidl)),
        'note': ('NOT 医保 reimbursement. The notice (第三部分（八）): "商保创新药目录内药品医保基金不予支付", '
                 'recommended for commercial health insurance / 医疗互助 reference.'),
    }
    src_stats = {
        'title': '2000年-2025年国家基本医保药品目录相关数据', 'document_no': None, 'issuer': '国家医疗保障局',
        'published': '2026-08-05', 'url': M['stats_page']['url'], 'attachment_url': M['stats_pdf']['url'],
        'retrieved_at': M['stats_pdf']['retrieved_at'], 'sha256': M['stats_pdf']['sha256'],
        'bytes': M['stats_pdf']['bytes'], 'rows': 0, 'stated_total': stats['west'],
        'stated_2025': stats,
        'cross_check': 'used only for its stated 2025 totals (西药 1857 = 1446 + 411; 谈判 472), which match',
        'note': 'no rows taken from this source',
    }
    out = {
        'schema': 1,
        'title': '国家医保药品目录 (NRDL) 2025年版, Western medicines',
        'provenance': {
            'built_at': now_utc(), 'builder': 'tools/china_access/build_nrdl.py',
            'edition': '2025年版', 'document_no': n_gov['document_no'], 'effective': '2026-01-01',
            'replaces': '国家基本医疗保险、工伤保险和生育保险药品目录（2024年）（医保发〔2024〕33号）',
            'sources': {'nrdl2025': src_nrdl, 'cidl2025': src_cidl, 'nhsa_stats_2000_2025': src_stats},
            'counts': counts,
            'coverage': (
                '国家基本医疗保险、生育保险和工伤保险药品目录（2025年）, in force since 2026-01-01 and still the current '
                'edition on 2026-10-06 (the latest 2026-adjustment item on nhsa.gov.cn: the 形式审查公示 解读 of 2026-06-29). Included: the '
                f'whole 西药部分 ({len(w)} rows = {len(prim)} numbers plus "★" repeat rows for other dosage forms), '
                f'the 协议期内谈判药品部分 (一)西药 ({nw}) and the 竞价药品部分 ({nb}), each with number, name, dosage '
                'form, 甲/乙 class, the 备注 restriction verbatim, and for 谈判/竞价 the 医保支付标准 and 协议有效期 '
                '(竞价: 支付标准有效期); the 凡例 clauses verbatim (fanli) and the two 凡例 lists of what rows 207 '
                'and 1303 cover. Separately, the 商业健康保险创新药品目录（2025年） (19 drugs), which is NOT paid by '
                '医保. NOT included: 中成药部分 (1335), 协议期内谈判药品部分 (二)中成药 (61), 中药饮片部分 (892 and the '
                'excluded list), English/INN names (the documents give none), specifications/manufacturers (the '
                '目录 does not list them), provincial 双通道/单独支付 lists and local reimbursement ratios. '
                'Restrictions apply to 基本医疗保险/生育保险 only (凡例（十一）: 工伤保险 is not bound by them); '
                'rows 292-309 and 1367-1381 carry an extra rule stated only in the 凡例 (see_fanli).'),
            'not_retrieved': [
                {'what': 'Excel (or any non-PDF) copy of the 2025 目录', 'url': M['notice_nhsa']['url'],
                 'reason': 'not published: both the NHSA and gov.cn notices attach only the two PDFs'},
                {'what': '2026年版 目录', 'url': 'https://www.nhsa.gov.cn/art/2026/6/29/art_105_21130.html',
                 'reason': 'not issued yet on 2026-10-06 (2026 adjustment in progress)'},
                {'what': 'previous (2024) edition, for an optional derived "not in the 2024 谈判 list" flag',
                 'url': 'https://www.gov.cn/zhengce/zhengceku/202411/P020241128820415409368.pdf',
                 'reason': ('optional, not bundled: the PDF downloads (HTTP 200, 1459033 bytes, 2026-10-06) but its '
                            'tables are laid out differently (e.g. "★（1060）" split across two cells on p.64, text '
                            'crossing row rules on p.187), so this parser\'s integrity checks reject it')},
                {'what': 'per-drug list of the 2025 additions (新增)',
                 'url': M['notice_nhsa']['url'],
                 'reason': ('neither the notice nor the 目录 marks which drugs were added in 2025, and no official '
                            'per-drug list was found (gov.cn policy-library search); NHSA\'s statistics sheet gives '
                            'only the number (新增谈判准入药品 2025: %d)' % stats['newly_negotiated_2025'])},
            ],
        },
        'categories': categories,
        'category_name_variants': cat_variants,
        'fanli': clauses,
        'fanli_dosage_form_groups': fanli_forms,
        'drugs': drugs,
        'commercial_innovative_list_label': ('商业健康保险创新药品目录（2025年）：不属于基本医疗保险支付范围，'
                                             '医保基金不予支付（医保发〔2025〕33号 三、（八）：“商保创新药目录内药品医保基金不予支付”）'),
        'commercial_innovative_list': cidl,
    }
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, 'w', encoding='utf-8') as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
        f.write('\n')
    print(json.dumps({'out': a.out, 'bytes': os.path.getsize(a.out), 'counts': counts,
                      'space_joins': sum(1 for b in breaks + cbreaks if b[0] == 'space-join'),
                      'stacked_name_rows': stacked, 'skipped_tables': skipped,
                      'wall_s': round(time.time() - t0, 1)}, ensure_ascii=False, indent=1), file=sys.stderr)
    return breaks + cbreaks


if __name__ == '__main__':
    main()
