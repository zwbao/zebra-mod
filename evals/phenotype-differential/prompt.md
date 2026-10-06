---
description: "Clinician gives HPO terms with an excluded one: ranking from phenotype_rank, scores not invented, sources kept apart"
tags: [journey-clinician, diagnosis]
runs: 3
max_turns: 30
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, "mcp__zebra-mod__*"]
---

Phenotype-driven differential, please. Present: HP:0002373, HP:0007359, HP:0001263, HP:0002133. Excluded: HP:0001252. Give me the top candidates and which genes to test.
