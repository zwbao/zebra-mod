---
type: regex
target: trace
pattern: '"type":"tool_use","id":"[^"]*","name":"mcp__zebra-mod__(trials_search|access)"'
match: contains
weight: 1
---

<!-- Trials come from trials_search or from access, which lists trials with sites in China; either counts. -->
