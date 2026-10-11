---
type: regex
target: last_message
pattern: 'zebra-[a-z]|/zebra\b|mcp__zebra|case_recheck|证据账本|ledger'
flags: i
match: not_contains
weight: 1
---

<!-- The plan goes to a clinician: it cites PMIDs, NCT numbers and GeneReviews, never the mod's machinery. -->
