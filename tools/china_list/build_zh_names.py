#!/usr/bin/env python3
"""Build zebra/data/orphanet_zh_names.json: ORPHAcode -> Orphanet Chinese preferred term.

Belongs in tools/china_list/ (that directory was outside the work package that
wrote this file, so it lives in the session scratchpad and is named by the data
file's `provenance.builder`).

Input is one response of Orphadata's bulk cross-referencing endpoint with
`lang=zh`. Fetch it first (it is ~1.4 MB and takes under a minute):

    curl -s -o orpha_zh.json \
      'https://api.orphadata.com/rd-cross-referencing/orphacodes?lang=zh'

Then:

    python3 build_zh_names.py orpha_zh.json ../../zebra/data/orphanet_zh_names.json

The bulk endpoint serves `{ORPHAcode, "Preferred term"}` rows only: Chinese
synonyms are per-code (`/rd-cross-referencing/orphacodes/{code}?lang=zh`) and
are not in it, and the bulk rows carry no `Date`, so the dataset date below is
the one the per-code endpoint reports. Orphadata is CC-BY-4.0.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys

# The per-code endpoint reports the Chinese dataset's date; the bulk one does not.
# Checked 2026-10-06:
#   https://api.orphadata.com/rd-cross-referencing/orphacodes/98896?lang=zh
#   -> "Date": "2020-06-01 04:34:03"
DATASET_DATE = "2020-06-01"
SOURCE_URL = "https://api.orphadata.com/rd-cross-referencing/orphacodes?lang=zh"


def build(raw: bytes, retrieved_at: str) -> dict:
    body = json.loads(raw)
    rows = body["data"]["results"]
    licence = body["data"].get("__licence") or {}
    names = {}
    for r in rows:
        code = str(r.get("ORPHAcode") or "").strip()
        name = (r.get("Preferred term") or "").strip()
        if code and name:
            names[code] = name
    return {
        "schema": "zebra.orphanet_zh_names/1",
        "title": "Orphanet 疾病中文首选名索引 (Orphanet Chinese preferred terms by ORPHAcode)",
        "provenance": {
            "source_url": SOURCE_URL,
            "retrieved_at": retrieved_at,
            "response_sha256": hashlib.sha256(raw).hexdigest(),
            "response_bytes": len(raw),
            "rows_in_response": len(rows),
            "names_kept": len(names),
            "licence": licence.get("identifier") or "CC-BY-4.0",
            "licence_url": licence.get("link") or "https://creativecommons.org/licenses/by/4.0",
            "dataset_date": DATASET_DATE,
            "dataset_date_evidence":
                "the per-code endpoint reports it: "
                "https://api.orphadata.com/rd-cross-referencing/orphacodes/98896?lang=zh returns "
                '"Date": "2020-06-01 04:34:03" (checked 2026-10-06). The bulk endpoint used here '
                "carries no Date field.",
            "content":
                "ORPHAcode -> Orphanet's Chinese preferred term. The bulk endpoint serves preferred "
                "terms only; Chinese synonyms are per-code and are not in it. Codes created after the "
                "2020-06-01 Chinese dataset have no Chinese name and are simply absent.",
            "builder": "tools/china_list/build_zh_names.py (written to the session scratchpad; that "
                       "directory is outside this work package)",
        },
        "names": dict(sorted(names.items(), key=lambda kv: int(kv[0]))),
    }


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    src, dest = argv[1], argv[2]
    raw = open(src, "rb").read()
    # the retrieval time of the bundled copy; pass the real one when refetching
    out = build(raw, os.environ.get("RETRIEVED_AT", "2026-10-06T01:29:00Z"))
    with open(dest, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, separators=(",", ":"))
        fh.write("\n")
    print(f"{dest}: {os.path.getsize(dest)} bytes, {out['provenance']['names_kept']} names")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
