#!/usr/bin/env python3
"""Second pass over zebra/data/china_disease_aliases.json: hand-entered folk names.

Belongs in tools/china_list/ next to build_china_aliases.py, which does the
first pass (the retrieved layer: official list names, their bracketed parts, and
Orphanet's English and Chinese preferred terms and synonyms for each linked
ORPHAcode). This pass adds only names that come from no dataset, so every one is
marked `kind: "folk_name"` and carries the official `name_zh` it maps to.

Idempotent: run it again and nothing changes. The build fails if a target entry
is not on the list, so a name can never be attached to a disease that is not
there.

    python3 add_folk_names.py ../../zebra/data/china_disease_aliases.json

Why each name is here:
  瑞特综合征 / 雷特综合征  the usual Chinese transliterations of Rett; the first
      pass left 瑞特综合征 matching 白塞病 (Behçet) at 0.42 bigram overlap, which
      is the defect this alias layer exists to fix
  蝴蝶宝贝 / 蝴蝶宝宝     "butterfly children", epidermolysis bullosa
  黏宝宝 / 粘宝宝         mucopolysaccharidosis (both spellings of 黏/粘 are used)
  月亮孩子               "moon children", albinism
  不食人间烟火的孩子      phenylketonuria
  玻璃人                 haemophilia
"""

from __future__ import annotations

import json
import sys

SOURCE = ("folk name or Chinese transliteration, hand-entered; not from a retrieved dataset "
          "(second builder pass, 2026-10-06)")

# (list, number) -> the folk names that map to that entry
SUPPLEMENT = {
    (2, 72): ["瑞特综合征", "雷特综合征"],      # Rett综合征 / Rett syndrome
    (1, 39): ["蝴蝶宝贝", "蝴蝶宝宝"],          # 遗传性大疱性表皮松解症
    (1, 73): ["黏宝宝", "粘宝宝"],              # 黏多糖贮积症
    (1, 2): ["月亮孩子"],                       # 白化病
    (1, 90): ["不食人间烟火的孩子"],             # 苯丙酮尿症
    (1, 36): ["玻璃人"],                        # 血友病
}

NOTE = ("second builder pass (tools/china_list/add_folk_names.py) 2026-10-06 added {n} hand-entered "
        "folk names and Chinese transliterations (瑞特/雷特综合征 for Rett综合征 — the first pass left "
        "瑞特综合征 matching 白塞病, the defect this layer exists to fix; 蝴蝶宝贝/蝴蝶宝宝, 黏宝宝/粘宝宝, "
        "月亮孩子, 不食人间烟火的孩子, 玻璃人). Each is marked folk_name and names the official list entry "
        "it maps to; none comes from a retrieved dataset.")


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    path = argv[1]
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    by_key = {(e["list"], e["no"]): e for e in data["entries"]}
    added = 0
    for key, texts in SUPPLEMENT.items():
        if key not in by_key:
            raise SystemExit(f"list entry {key} is not on the national lists: refusing to add {texts}")
        entry = by_key[key]
        have = {(a["text"], a["kind"]) for a in entry["aliases"]}
        for text in texts:
            if (text, "folk_name") in have:
                continue
            entry["aliases"].append({"text": text, "kind": "folk_name",
                                     "maps_to": entry["name_zh"], "source": SOURCE})
            added += 1
    if not added:
        print(f"{path}: already up to date")
        return 0
    prov = data["provenance"]
    prov["counts"]["aliases"] += added
    prov["counts"]["folk_names"] += added
    prov.setdefault("notes", []).append(NOTE.format(n=added))
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=1)
        fh.write("\n")
    print(f"{path}: added {added} folk names; counts {prov['counts']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
