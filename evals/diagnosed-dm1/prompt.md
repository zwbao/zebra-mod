---
description: "Clinician asks for a management plan for a woman with genetically confirmed DM1: the diagnosis is taken as established, the conduction disease and the planned surgery are handled as urgent safety issues, pregnancy plans get genetic counselling, a concrete clinician plan"
tags: [journey-clinician, zh, diagnosed, care]
runs: 4
max_turns: 40
timeout_seconds: 900
allowed_tools: [Read, Glob, Grep, Skill, WebSearch, WebFetch, "mcp__zebra-mod__*"]
---

我是神经内科医生。下面是一位患者的资料摘要（已去标识化），请给出具体的管理、预后、营养等尽可能详细的方案。

患者：女，32 岁，已婚，计划明年怀孕。
基因检测（2025）：DMPK 基因 3'UTR CTG 重复，一个等位约 600 次，另一个 13 次。报告结论：符合强直性肌营养不良 1 型（DM1）。
症状：握拳后松手慢 5 年；近 2 年下楼梯费力、常绊倒；白天嗜睡明显，开会时会睡着；30 岁做过双眼白内障手术。
查体（外院简要记录）：颞肌、胸锁乳突肌萎缩，叩击性肌强直（+），足背屈 4 级，握力减弱。
2026-08：CK 310 U/L，GGT 95 U/L（<45），ALT 48 U/L，空腹血糖 6.4 mmol/L，HbA1c 6.3%，TSH 正常，睾酮未查。
心电图：窦律，PR 236 ms，QRS 126 ms（右束支传导阻滞型）。未做 Holter 和心超。
肺功能：FVC 78% 预计值。
另：外科建议下个月做腹腔镜胆囊切除（胆囊结石反复发作）。
家族史：母亲 50 多岁有白内障，"手没力气"，未检查；一个弟弟 28 岁，无症状。
