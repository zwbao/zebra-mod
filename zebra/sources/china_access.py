"""Bundled China-access tables: NMPA rare-disease approvals and the national reimbursement list (医保).

Both files are built from official documents by scripts in tools/china_access/
and carry their provenance (URL, retrieval time, sha256, counts, coverage):

  zebra/data/china_rare_drug_approvals.json   NMPA / CDE / MOF documents naming
      rare-disease drugs approved or listed in China (build_drug_approvals.py).
      There is no complete official register of NMPA rare-disease approvals, so
      the file states its own coverage: a drug missing from it may still be
      approved in China.
  zebra/data/china_nrdl.json   国家基本医疗保险、生育保险和工伤保险药品目录 —
      the Western-medicine part in full, the negotiated (协议期内谈判) and bidding
      (竞价) tables with their payment restrictions verbatim, and separately the
      商业健康保险创新药品目录, which the 医保 fund does NOT pay for (build_nrdl.py).

Reading only, standard library only; nothing here touches the network. A
missing or unreadable file raises `NotBundled`, which callers report as "not
bundled", never as "not approved" or "not covered".
"""

from __future__ import annotations

import json
import os
import re
import unicodedata
from typing import Any, Dict, Iterable, List, Optional, Sequence

from zebra.sources import record as source_record

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
APPROVALS_FILE = os.path.join(DATA_DIR, "china_rare_drug_approvals.json")
NRDL_FILE = os.path.join(DATA_DIR, "china_nrdl.json")
_cache: Dict[str, Any] = {}
_CJK = re.compile(r"[㐀-鿿]")

# dosage-form words that end (or, for 注射用, begin) an NRDL / NMPA product name
_FORM_PREFIX = ("注射用",)
# salt words written before the moiety in Chinese names (盐酸沙丙蝶呤 = sapropterin hydrochloride)
_SALT_PREFIX = ("盐酸", "硫酸氢", "硫酸", "磷酸氢", "磷酸", "枸橼酸", "马来酸", "甲磺酸", "二甲磺酸", "乙磺酸", "醋酸",
                "酒石酸", "富马酸", "琥珀酸", "苯磺酸", "氢溴酸", "乳酸", "葡萄糖酸", "甲苯磺酸", "对甲苯磺酸", "门冬氨酸",
                "二盐酸", "双盐酸", "依地酸", "羟乙磺酸")
# words after a disease name that make the text about a different, secondary condition
# (真性红细胞增多症继发的骨髓纤维化 is myelofibrosis, not polycythaemia vera)
_SECONDARY = re.compile(r"^(继发|后|相关|所致|引起|导致|伴发|转化)")
# a Chinese query ending like a disease name is never read as a drug name (胰岛素瘤 is not 胰岛素)
DISEASE_TAIL = re.compile(r"(瘤|症|病|缺乏|综合征|综合症|血症|炎|癌|畸形|障碍|不全|发育不良|萎缩)$")
_FORM_SUFFIX = tuple(sorted((
    "口服溶液用散", "注射用浓溶液", "口服混悬液", "肠溶干混悬剂", "干混悬剂", "口服溶液", "注射液", "混悬注射液",
    "混悬液", "肠溶胶囊", "缓释胶囊", "软胶囊", "缓释片", "肠溶片", "分散片", "咀嚼片", "口崩片", "胶囊", "颗粒", "片",
    "散", "注射剂", "滴眼液", "眼用凝胶", "乳膏", "凝胶", "吸入溶液", "吸入粉雾剂", "吸入剂", "气雾剂", "鼻喷雾剂",
    "喷雾剂", "植入剂", "贴剂", "栓", "丸", "糖浆", "口服液", "溶液剂", "溶液", "粉雾剂", "乳剂", "外用溶液"),
    key=len, reverse=True))
# the only characters a longer name may add to a shorter one and still name the same moiety
_SALT_SUFFIX = ("二钠", "三钠", "葡胺", "钠", "钾", "钙", "镁", "锌", "锂", "铵")


class NotBundled(Exception):
    pass


def _load(key: str, path: str, rows_key: str) -> Dict[str, Any]:
    if key not in _cache:
        if not os.path.exists(path):
            raise NotBundled(f"not bundled: {os.path.basename(path)} is missing")
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
        except (OSError, ValueError) as err:
            raise NotBundled(f"{os.path.basename(path)} unreadable: {err}") from None
        if not isinstance(data, dict) or not data.get(rows_key):
            raise NotBundled(f"{os.path.basename(path)} holds no {rows_key}")
        _cache[key] = data
    return _cache[key]


def load_approvals() -> Dict[str, Any]:
    return _load("approvals", APPROVALS_FILE, "drugs")


def load_nrdl() -> Dict[str, Any]:
    return _load("nrdl", NRDL_FILE, "drugs")


def zh_norm(text: str) -> str:
    t = unicodedata.normalize("NFKC", text or "").lower()
    return re.sub(r"[\s\-_'\",.;:()\[\]{}/、，。（）·]+", "", t)


def zh_generic(name: str) -> str:
    """The generic part of a Chinese product name: 诺西那生钠注射液 -> 诺西那生钠, 注射用维拉苷酶β -> 维拉苷酶β."""
    t = zh_norm(name)
    for p in _FORM_PREFIX:
        if t.startswith(p) and len(t) > len(p) + 1:
            t = t[len(p):]
    for p in sorted(_SALT_PREFIX, key=len, reverse=True):
        if t.startswith(p) and len(t) > len(p) + 1:
            t = t[len(p):]
            break
    changed = True
    while changed:
        changed = False
        for s in _FORM_SUFFIX:
            if t.endswith(s) and len(t) > len(s) + 1:
                t = t[:-len(s)]
                changed = True
                break
    return t


def _latin_words(text: str) -> List[str]:
    return [w for w in re.split(r"[^a-z0-9]+", (text or "").lower()) if w]


def same_drug_zh(a: str, b: str) -> bool:
    """Two Chinese drug names name the same moiety: equal generics, or one generic inside the other (salt forms)."""
    ga, gb = zh_generic(a), zh_generic(b)
    if not ga or not gb or not _CJK.search(ga) or not _CJK.search(gb):
        return False
    if ga == gb:
        return True
    short, long_ = sorted((ga, gb), key=len)
    # 诺西那生 ~ 诺西那生钠 (a salt). Anything else added is another product: 胰岛素 is not 胰岛素瘤,
    # 维拉苷酶 is not 维拉苷酶β, 葡萄糖 is not 葡萄糖酸钙 (adversarial review P0-4)
    return len(short) >= 3 and long_.startswith(short) and long_[len(short):] in _SALT_SUFFIX


def clean_inn(text: Optional[str]) -> Optional[str]:
    """The INN in a bundled English name: 'Galafold（Migalastat hydrochloride）' -> 'Migalastat hydrochloride'."""
    t = unicodedata.normalize("NFKC", text or "").strip()
    if not t:
        return None
    m = re.search(r"\(([^()]+)\)", t)
    if m and re.search(r"[a-z]", m.group(1)):
        t = m.group(1)
    t = re.sub(r"\b(for injection|injection|tablets?|capsules?|oral solution|for oral suspension|solution|"
               r"concentrate for solution for infusion|powder)\b", " ", t, flags=re.I)
    t = " ".join(t.split())
    return t or None


def mentions(text: str, name: str) -> bool:
    """`text` names `name` as the condition itself, not as the cause of a secondary one."""
    t, n = zh_norm(text), zh_norm(name)
    if not n:
        return False
    start = t.find(n)
    while start >= 0:
        if not _SECONDARY.match(t[start + len(n):]):
            return True
        start = t.find(n, start + 1)
    return False


def same_drug_latin(query: str, inn: Optional[str]) -> bool:
    """`query` names the moiety of `inn` (word-bounded): 'nusinersen' ~ 'nusinersen sodium'."""
    from zebra.sources.regulators import moiety

    q = moiety(query)
    i = moiety(clean_inn(inn) or "")
    if not q or not i:
        return False
    return q == i or re.search(r"(?<![a-z0-9])" + re.escape(q) + r"(?![a-z0-9])", i) is not None or \
        re.search(r"(?<![a-z0-9])" + re.escape(i) + r"(?![a-z0-9])", q) is not None


def approvals_for_drug(names: Iterable[str]) -> List[Dict[str, Any]]:
    """Rows of the bundled NMPA/CDE/MOF table for a drug given by Chinese and/or English (INN) names."""
    data = load_approvals()
    names = [n for n in names if n]
    out = []
    for row in data["drugs"]:
        hit = None
        zh_names = [row.get("drug_zh") or ""] + list(row.get("names_zh") or [])
        for n in names:
            if _CJK.search(n):
                z = next((x for x in zh_names if x and same_drug_zh(n, x)), None)
                if z:
                    hit = f"Chinese name {z}"
            elif row.get("inn") and same_drug_latin(n, row["inn"]):
                hit = f"English name {row['inn']}"
            if hit:
                break
        if hit:
            out.append(dict(row, matched_by=hit))
    return out


def approvals_for_disease(list_keys: Sequence[tuple], names_zh: Iterable[str] = (),
                          names_en: Iterable[str] = ()) -> List[Dict[str, Any]]:
    """Rows linked to a national-list entry, or whose published indication names the disease."""
    data = load_approvals()
    keys = {tuple(k) for k in list_keys}
    zh = [zh_norm(n) for n in names_zh if n and len(zh_norm(n)) >= 3]
    en = [n.lower() for n in names_en if n and len(n) >= 5]
    out = []
    for row in data["drugs"]:
        why = None
        texts = [row.get("indication_zh") or ""] + [e.get("indication_zh") or "" for e in row.get("evidence") or []]
        texts = [t for t in texts if t]
        for d in row.get("list_diseases") or []:
            if (d.get("list"), d.get("no")) in keys:
                name = d.get("name_zh") or ""
                # the builder linked on the literal name; a text that names it only as the cause of
                # another condition (真性红细胞增多症继发的骨髓纤维化) does not make it this drug's indication
                if name and texts and any(zh_norm(name) in zh_norm(t) for t in texts) \
                        and not any(mentions(t, name) for t in texts):
                    continue
                why = f"linked to national list entry {d.get('list')}#{d.get('no')} ({d.get('how') or 'by the builder'})"
                break
        if not why:
            ind_zh = row.get("indication_zh") or ""
            ind_en = (row.get("indication_en") or "").lower()
            z = next((n for n in zh if mentions(ind_zh, n)), None)
            e = next((n for n in en if re.search(r"(?<![a-z0-9])" + re.escape(n) + r"(?![a-z0-9])", ind_en)), None)
            if z or e:
                why = f"its published indication names '{z or e}'"
        if why:
            out.append(dict(row, matched_by=why))
    return out


def nrdl_for_drug(names_zh: Iterable[str]) -> List[Dict[str, Any]]:
    """NRDL rows (and 商保创新药目录 rows, labelled) for Chinese drug names."""
    data = load_nrdl()
    names = [n for n in names_zh if n and _CJK.search(n)]
    out = []
    for row in data["drugs"]:
        if any(same_drug_zh(n, row.get("name_zh") or "") for n in names):
            out.append(dict(row, list="NRDL"))
    for row in data.get("commercial_innovative_list") or []:
        if any(same_drug_zh(n, row.get("name_zh") or "") for n in names):
            out.append(dict(row, list="CIDL"))
    return out


def nrdl_for_disease(names_zh: Iterable[str]) -> List[Dict[str, Any]]:
    """NRDL rows whose payment restriction names the disease, and CIDL rows whose indication does."""
    data = load_nrdl()
    names = [zh_norm(n) for n in names_zh if n and _CJK.search(n) and len(zh_norm(n)) >= 3]
    out = []
    if not names:
        return out
    for row in data["drugs"]:
        hit = next((n for n in names if mentions(row.get("restriction") or "", n)), None)
        if hit:
            out.append(dict(row, list="NRDL", matched_by=f"its 医保 restriction names '{hit}'"))
    for row in data.get("commercial_innovative_list") or []:
        hit = next((n for n in names if mentions(row.get("indication") or "", n)), None)
        if hit:
            out.append(dict(row, list="CIDL", matched_by=f"its indication names '{hit}'"))
    return out


def compact_nrdl(row: Dict[str, Any], nrdl: Dict[str, Any]) -> Dict[str, Any]:
    """The fields a family needs from one NRDL / CIDL row, with what the list means."""
    p = nrdl.get("provenance") or {}
    if row.get("list") == "CIDL":
        return {"list": "商业健康保险创新药品目录 (commercial insurance list — 医保 does NOT pay)",
                "name_zh": row.get("name_zh"), "brand_zh": row.get("brand_name_zh"),
                "indication": row.get("indication"), "validity": row.get("validity"), "no": row.get("no"),
                "matched_by": row.get("matched_by"), "source_id": row.get("source_id")}
    return {"list": f"国家医保药品目录 {p.get('edition', '')} ({p.get('document_no', '')}, in force {p.get('effective', '')})",
            "section": row.get("section"), "no": row.get("no"), "name_zh": row.get("name_zh"),
            "dosage_form": row.get("dosage_form"), "class": row.get("class"),
            "restriction": row.get("restriction") or "no restriction printed (备注 empty)",
            "agreement_period": row.get("agreement_period"), "matched_by": row.get("matched_by"),
            "source_id": row.get("source_id")}


def sources_for(data: Dict[str, Any], source_ids: Iterable[str], db: str) -> List[Dict[str, Any]]:
    """Provenance rows for the documents behind the rows used (all of them when none matched)."""
    p = data.get("provenance") or {}
    srcs = p.get("sources") or {}
    ids = [s for s in dict.fromkeys(source_ids) if s in srcs] or list(srcs)
    recs = []
    for sid in ids:
        s = srcs[sid]
        recs.append(source_record(db, s.get("document_no") or s.get("title") or sid, url=s.get("url"),
                                  note=f"{s.get('title', '')}; {s.get('issuer', '')} {s.get('published', '')}; "
                                       f"bundled copy retrieved {s.get('retrieved_at')}, sha256 "
                                       f"{str(s.get('sha256') or '')[:12]}…"
                                       + (f"; origin {s['origin']}" if s.get("origin") and s.get("origin") != s.get("url") else "")))
    return recs


def coverage(data: Dict[str, Any]) -> Optional[str]:
    return (data.get("provenance") or {}).get("coverage")
