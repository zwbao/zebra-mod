---
description: "Chinese parent shares a clinic note and asks what it could be: Chinese reply, HPO ids from hpo_search, ranking from phenotype_rank, no identifiers sent out"
tags: [journey-family, zh, diagnosis]
runs: 3
max_turns: 30
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, "mcp__zebra-mod__*"]
---

这是我女儿的门诊病历：

门诊病历（摘录）
姓名：王小雨    性别：女    年龄：2岁3个月    门诊号：MZ20240518007
联系电话：13912345678
主诉：反复发热时抽搐1年余。
现病史：患儿6月龄起出现发热时全身抽搐，持续10-20分钟，此后多次发作，部分发作为一侧肢体抽动。
1岁后出现不发热时的抽搐。走路晚，1岁8个月才会独走，现在说话少，只会叫“爸爸”“妈妈”。
查体：肌张力正常。
辅助检查：头颅MRI未见异常；脑电图（1岁）：背景正常。
家族史：无类似疾病。

她这两年老是发烧就抽，现在不发烧也会抽，走路也晚。医生说可能是癫痫，但一直查不出原因。这可能是什么病？我们下一步该问医生什么？
