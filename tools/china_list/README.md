# Rebuilding zebra/data/china_rare_diseases.json

`parse_cn.py` and `build_cn.py` rebuild the bundled list from the official sources named in the JSON's
`provenance` block: list 1 from the gov.cn policy-library HTML table, list 2 from the gov.cn `.doc` attachment
converted to text with macOS `textutil -convert txt`. Download the pages into this folder first; the scripts
record each source's sha256. Names are stored exactly as published; differences between official copies are
noted in `provenance.lists.*.cross_check`.
