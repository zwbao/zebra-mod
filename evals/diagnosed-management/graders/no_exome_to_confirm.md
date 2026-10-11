---
type: regex
target: last_message
pattern: '(外显子|全基因组|基因\s*panel|基因板|exome|genome sequencing)[^。\n]{0,30}(确认|确诊|证实|confirm)[^。\n]{0,15}(FSHD|面肩肱)'
flags: i
match: not_contains
weight: 2
---
