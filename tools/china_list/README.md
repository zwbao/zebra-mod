# Rebuilding zebra/data/china_rare_diseases.json

`parse_cn.py` and `build_cn.py` rebuild the bundled list from the official sources named in the JSON's
`provenance` block: list 1 from the gov.cn policy-library HTML table, list 2 from the gov.cn `.doc` attachment
converted to text with macOS `textutil -convert txt`. Download the pages into this folder first; the scripts
record each source's sha256. Names are stored exactly as published; differences between official copies are
noted in `provenance.lists.*.cross_check`.

## The alias layer

`build_zh_names.py` rebuilds `zebra/data/orphanet_zh_names.json` from Orphadata's Chinese dataset
(CC BY 4.0; the bundled copy is the 2020-06-01 release, so ORPHAcodes created later have no Chinese
name). `add_folk_names.py` adds the folk names families use (瓷娃娃, 渐冻症, 小胖威利, 快乐木偶) to
`zebra/data/china_disease_aliases.json`; it is idempotent, and every alias carries its `kind` and
`source`, with folk names marked `folk_name` and pointing at the official entry they map to.
Known limits are recorded in the data files themselves: 58 of the 207 list entries have no single
ORPHAcode, and the list→ORPHA link is by name rather than through Orphanet's classification.
