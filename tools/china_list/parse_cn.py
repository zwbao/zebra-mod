import re, html, json, sys

def cells_from_html(path):
    t = open(path, encoding='utf-8', errors='replace').read()
    rows = []
    for tr in re.findall(r'<tr[^>]*>(.*?)</tr>', t, flags=re.S | re.I):
        cells = re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', tr, flags=re.S | re.I)
        cells = [re.sub(r'\s+', ' ', html.unescape(re.sub(r'<[^>]+>', '', c))).strip() for c in cells]
        if cells:
            rows.append(cells)
    return rows

def list1(path):
    out = []
    for r in cells_from_html(path):
        r = [c for c in r if c != '']
        if len(r) >= 3 and re.fullmatch(r'\d{1,3}', r[0]):
            out.append({'no': int(r[0]), 'zh': r[1], 'en': ' '.join(r[2:])})
    return out

for p in sys.argv[1:]:
    L = list1(p)
    print(p, len(L), L[:2], L[-2:])
    nums = [x['no'] for x in L]
    print(' numbers ok:', nums == list(range(1, len(L) + 1)))
