import re, json, hashlib
exec(open('parse_cn.py').read().split("for p in sys.argv")[0])
CJK = re.compile(r'[㐀-鿿（）]')

def list2(path):
    t = open(path, encoding='utf-8').read().replace('\x07', '\n')
    lines = [l.strip() for l in t.split('\n') if l.strip()]
    i = lines.index('疾病名称（英文）') + 1
    out, cur = [], None
    for l in lines[i:]:
        m = re.fullmatch(r'(\d{1,3})', l)
        m2 = re.match(r'^(\d{1,3})(?=[^\d\-])(.*)$', l)
        nxt = (cur['no'] + 1) if cur else 1
        if m and int(m.group(1)) == nxt:
            cur = {'no': nxt, 'parts': []}; out.append(cur); continue
        if m2 and int(m2.group(1)) == nxt and not l[len(m2.group(1)):].startswith(('型', '-')):
            cur = {'no': nxt, 'parts': [m2.group(2)]}; out.append(cur); continue
        cur['parts'].append(l)
    res = []
    for e in out:
        zh, en = [], []
        for p in e['parts']:
            # a part may hold both: Chinese name immediately followed by English
            m = re.match(r'^(.*[㐀-鿿）】])\s*([A-Za-z][^㐀-鿿]*)$', p)
            if m and CJK.search(m.group(1)) and not CJK.search(m.group(2)):
                zh.append(m.group(1)); en.append(m.group(2))
            elif CJK.search(p):
                zh.append(p)
            else:
                en.append(p)
        res.append({'no': e['no'], 'zh': ''.join(zh).strip(), 'en': ' '.join(en).strip()})
    return res

L1 = list1('cn1_zhengceku.html')[:121]
L1g = list1('cn1_gongbao.html')
L2 = list2('cn2_govcn.txt'); L2m = list2('cn2_miit.txt')
assert [x['no'] for x in L1] == list(range(1, 122)), 'list1 numbering'
assert [x['no'] for x in L2] == list(range(1, 87)), ('list2 numbering', [x['no'] for x in L2])
for x in L2:
    assert x['zh'] and x['en'], x
diff2 = [(a, b) for a, b in zip(L2, L2m) if a != b]
print('list2 gov.cn vs MIIT diffs:', diff2)
zh_diff1 = [(a['no'], a['zh'], b['zh']) for a, b in zip(L1, L1g) if re.sub(r'[\s（）()]', '', a['zh']) != re.sub(r'[\s（）()]', '', b['zh'])]
print('list1 zh diffs (ignoring spaces/brackets):', zh_diff1)
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
data = {
  'schema': 'zebra.china_rare/1',
  'title': '国家罕见病目录 (China national rare disease lists)',
  'provenance': {
    'retrieved_at': '2026-10-05T15:01:53Z',
    'counts': {'1': len(L1), '2': len(L2), 'total': len(L1) + len(L2)},
    'lists': {
      '1': {'name_zh': '第一批罕见病目录', 'document': '关于公布第一批罕见病目录的通知 国卫医发〔2018〕10号', 'date': '2018-05-11',
            'issuers': '国家卫生健康委员会 科学技术部 工业和信息化部 国家药品监督管理局 国家中医药管理局',
            'source_url': 'https://www.gov.cn/zhengce/zhengceku/2018-12/31/content_5435167.htm',
            'source_sha256': sha('cn1_zhengceku.html'),
            'cross_checked_against': 'https://www.gov.cn/gongbao/content/2018/content_5338244.htm',
            'cross_check': 'same 121 numbered entries and Chinese names (whitespace/bracket style and one 症 suffix differ); the State Council Gazette copy has English typos (e.g. "21-Hydroxyulase", "Ldrenoleuko Dystrophy") that the policy-library copy used here does not',
            'origin': 'NHC (www.nhc.gov.cn/yzygj/c100068/201806/bd1611850ff14bc8888c149567fe0a55.shtml); nhc.gov.cn serves a JavaScript challenge to scripts, so the gov.cn reproduction (来源：卫生健康委网站) was used'},
      '2': {'name_zh': '第二批罕见病目录', 'document': '关于公布第二批罕见病目录的通知 国卫医政发〔2023〕26号', 'date': '2023-09-18',
            'issuers': '国家卫生健康委 科技部 工业和信息化部 国家药监局 国家中医药局 中央军委后勤保障部',
            'source_url': 'https://www.gov.cn/zhengce/zhengceku/202309/content_6905273.htm',
            'attachment_url': 'https://www.gov.cn/zhengce/zhengceku/202309/P020260429440050482485.doc',
            'source_sha256': sha('cn2_govcn.doc'),
            'cross_checked_against': 'https://www.miit.gov.cn/jgsj/xfpgys/wjfb/art/2023/art_f25169aef87b47afa03862549f2125ed.html (attachment d1b717ee57c14f7bb674822c0383a459.doc)',
            'cross_check': 'identical 86 entries except no. 64, where the MIIT copy gives the English name "Primary ciliary dyskinesia" against the Chinese 原发性生长激素缺乏症; the gov.cn attachment (used here) reads "Primary growth hormone deficiency"'},
    },
    'method': 'list 1 parsed from the HTML table of the gov.cn policy-library page; list 2 from the gov.cn .doc attachment converted to text (macOS textutil); names stored exactly as published',
  },
  'diseases': [{'list': 1, 'no': x['no'], 'name_zh': x['zh'], 'name_en': x['en']} for x in L1] +
              [{'list': 2, 'no': x['no'], 'name_zh': x['zh'], 'name_en': x['en']} for x in L2],
}
json.dump(data, open('/Users/baozhiwei/zebra-mod/zebra/data/china_rare_diseases.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
print('written', data['provenance']['counts'])
for x in L2[:10] + L2[15:20] + L2[83:]: print(x)
