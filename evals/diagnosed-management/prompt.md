---
description: "Clinician brings an FSHD1 patient already diagnosed by D4Z4 sizing + methylation and asks for a management plan: the diagnosis is taken as established, no exome 'to confirm FSHD', missing exam not counted against it, a concrete plan in clinician register"
tags: [journey-clinician, zh, diagnosed, care]
runs: 4
max_turns: 40
timeout_seconds: 900
allowed_tools: [Read, Glob, Grep, Skill, WebSearch, "mcp__zebra-mod__*"]
---

我是神经内科医生，这是我的一个患者的资料摘要（已去标识化），请根据这个病例的情况，给出具体的管理、预后、营养等尽可能详细的方案。

患者：男，25 岁。
基因检测（2023，送检时临床诊断写的是"FSHD 待排除"）：4q35 EcoRI/p13E-11 片段一条缩短至 17.0 kb（正常 38–300 kb），4qA 型；另一条 >38 kb，4qA 型。4q35-D4Z4/PAS 甲基化 16.2%（正常 >69.5%）。实验室结论：结合临床表型，符合 FSHD1 的诊断。
2023-12：CK 693 U/L，LDH 319，肌酐 40 µmol/L，尿酸 583 µmol/L。肺功能：FVC 62% 预计值，TLC 73%，RV/TLC 190% 预计值，DLCO/VA 正常，结论中度限制性通气障碍。膈肌超声正常。心脏超声 LVEF 67%，各腔室正常。心电图：窦性心律不齐，不完全性右束支传导阻滞。肌肉 MRI：臀肌、腘绳肌脂肪浸润，左侧重；右小腿胫骨前肌、腓肠肌内侧头脂肪浸润，左侧同名肌 STIR 高信号。吸烟。
2026-06：CK 407 U/L，肌酐 34 µmol/L，胱抑素 C 0.92 mg/L，肌红蛋白 85 ng/mL，hs-TnT 0.083 ng/mL（≤0.014），CK-MB mass 15.3 ng/mL，NT-proBNP 19 pg/mL，IL-6 8.4 pg/mL，同型半胱氨酸 18.2 µmol/L，AFP 8.6 ng/mL（<7），甲状腺球蛋白抗体 311 U/mL（<115），TSH 正常。双侧大腿、小腿 MRI：萎缩、脂肪浸润，STIR 高信号。颈椎 MRI：颈 2–3 椎间盘水平脊髓可疑片状高信号。心电图窦律，HR 85。
手环心率 2 个多月：夜间均值 66，最高 143。
资料里没有查体记录和家族史。
