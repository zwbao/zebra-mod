# Rebuilding zebra/data/china_rare_diseases.json

`parse_cn.py` and `build_cn.py` rebuild the bundled list from the official sources named in the JSON's
`provenance` block: list 1 from the gov.cn policy-library HTML table, list 2 from the gov.cn `.doc` attachment
converted to text with macOS `textutil -convert txt`. Download the pages into this folder first; the scripts
record each source's sha256. Names are stored exactly as published; differences between official copies are
noted in `provenance.lists.*.cross_check`.

## The alias layer

`build_china_aliases.py` rebuilds `zebra/data/china_disease_aliases.json` from the bundled lists and
Orphadata (CC BY 4.0), keeping every response it fetches in a cache directory so a rebuild is reproducible:

    python3 -I tools/china_list/build_china_aliases.py --cache <dir> [--compare <old json>]

Pass 1 links each list entry to one ORPHAcode by exact rules only (rules a–d in the file's `provenance.rules`)
and takes the official names, their bracket/slash parts and the linked code's Orphanet terms as aliases; the
2026-10-06 rebuild reproduced the earlier (lost) builder's 149 links and every alias exactly
(`provenance.rebuild_check`). Pass 2 records what a match on each alias means: the unqualified head of the three
qualified entries (地中海贫血（重型）, 帕金森病（青年型、早发型）, 糖原累积病（I型、Ⅱ型）) is `scope: wider`; Orphanet
subtypes inside a qualifier carry their own ORPHAcode (Pompe disease, ORPHA:365, inside 糖原累积病 Ⅱ型); group
entries get their Orphanet descendants as members; short all-capital aliases are `acronym: true`. Pass 3 adds the
hand-entered folk names (marked `folk_name`, with the official name they map to). `build_zh_names.py` rebuilds
`zebra/data/orphanet_zh_names.json` (the 2020-06-01 Chinese dataset, so ORPHAcodes created later have no
Chinese name). `add_folk_names.py` was the second pass of the earlier builder; its table is now part of
`build_china_aliases.py`.

The China access tables (NMPA rare-disease approvals, the national reimbursement list, the 协作网 hospitals)
are built by the scripts in `tools/china_access/`; each script's docstring says how to rerun it.
