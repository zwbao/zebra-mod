---
description: "Chinese lay phrases become verified HPO ids: 头围小, 听力下降, 肌张力低下, 走路晚 (regressions D-P1-4, CP1-7)"
tags: [zh, hpo, regression]
runs: 3
max_turns: 15
timeout_seconds: 300
allowed_tools: [Read, Glob, Grep, Skill, "mcp__zebra-mod__*"]
---

帮我把孩子的这些情况转成标准的HPO术语，给出HPO编号：头围小，听力下降，肌张力低下，孩子走路晚。
