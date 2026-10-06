---
type: regex
target: trace
pattern: '"type":"tool_use","id":"[^"]*","name":"mcp__zebra-mod__(therapy_landscape|access)"'
match: contains
weight: 1
---

<!-- Drugs and approvals come from therapy_landscape (Open Targets/ChEMBL) or access (FDA/EMA/NMPA records); either counts. -->
