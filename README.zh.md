# 🦓 zebra-mod

**zebra-mod 是一个 Claude Code mod，把 Claude Code 改造成罕见病研究与分析工作台**，服务对象包括正在寻求诊断的患者和家属、临床医生与科研人员。

> 医学教育有一句格言："听到马蹄声，先想到马，而不是斑马。"罕见病患者自称"斑马"，因为他们的答案往往不在常见路径上。zebra-mod 的目标正是在这些路径之外寻找答案，并为每一步给出可核查的证据。

[English → README.md](README.md)

---

## 功能

| 任务 | 做法 |
|---|---|
| **鉴别诊断**：根据表型推断可能的疾病 | 病历 → HPO 表型档案（术语经核验，区分存在/排除，记录起病时间与出处页码）→ 三种相互独立的排序方法（基于 HPO 注释的离线 Resnik 相似度、Monarch 语义相似度、PubCaseFinder），各自保留排名，以方法间的一致性作为信号 → 鉴别特征表 → 每个候选疾病应采用何种检测（外显子组、基因组、染色体微阵列、重复扩增检测、甲基化、线粒体 DNA、代谢筛查） |
| **变异解读**：ACMG/AMP 框架 | VEP（MANE 转录本、AlphaMissense、REVEL、CADD、SpliceAI）+ gnomAD（grpmax、faf95）+ ClinVar（审核星级）+ ClinGen 基因-疾病有效性与剂量敏感性 + LitVar → 按 ClinGen 校准阈值给出数据驱动的证据条目 → 由模型论证的判断类条目 → **由代码计算积分**（Tavtigian 2020），并与 2015 年组合规则的结论并列呈现 |
| **序列到功能（S2F）**：剪接、非编码与调控变异 | SpliceAI/Pangolin（Broad 在线查询），AlphaGenome、Evo 2 与 GPN-MSA 通过 [s2f-penguin](https://github.com/zwbao/s2f-penguin) 的 `s2f` 命令行调用；逐轴解读并标注结论上限，以逻辑一致性综合，不做数值平均；由机制推导验证方案（RNA 检测） |
| **外显子组/基因组重分析**，在本机完成 | VCF 分诊：质量过滤、遗传模式（新发、纯合、复合杂合、X 连锁半合子）、按表型基因限定范围，仅对剩余候选变异调用 VEP 注释；VCF 文件不离开本机 |
| **基因型到治疗（G2T）** | 先判定致病机制（功能丧失、功能获得、显性负效应、剪接异常、重复扩增）→ 已上市与在研药物（Open Targets / ChEMBL）、正在招募的临床试验（ClinicalTrials.gov，含中国）、文献证据、个体化治疗可行性初筛（反义寡核苷酸、基因替代、**碱基编辑可行性**）→ 按 A–E 分级，并逐条核对机制方向 |
| **经典统计方法** | 共分离似然比与 PP1 强度、最大可信等位基因频率（BS1）、携带者频率与遗传患病率、再发风险（含 X 连锁贝叶斯计算）、Fisher 精确检验与负荷检验、新发变异富集、Kaplan–Meier 自然史分析与对数秩检验、N-of-1 试验设计与分析 |
| **文献** | Europe PMC、PubTator3、LitVar2；PMID 仅使用工具返回的结果，结论引用原文语句 |
| **面向家庭** | 通俗书面语解释（中文用户默认中文）、如实说明意义不明变异（VUS）、复诊问题清单、再发风险与家系筛查、患者组织、国家罕见病目录 |
| **报告** | 临床医生摘要与家属信，每条结论附出处，定稿前由**独立的对抗式审查**核对 |

## 原则（由机制保证，而非仅写在文档中）

1. **证据优先，不凭记忆。** 每个结论引用本次会话中检索到的证据账本编号、PMID 或数据库记录。每个病例维护一份只追加的**证据账本**（`E1、E2……`），记录用到的全部来源。
2. **代码计分，模型解读。** 排名、ACMG 积分与统计量均由 `zebra` 引擎计算；模型不编写、不修改任何数值。
3. **信息缺失时先澄清，绝不编造。** 基因组版本、转录本、合子状态、遗传方式、性别等关键信息缺失时，先询问，或以条件句给出结论。
4. **研究级分析，而非临床报告。** 不以确定口吻告知诊断，不提供用药剂量；后续步骤以"供医疗团队考虑的问题"形式呈现。
5. **隐私：如实说明边界。** 病例文件只写在本机；zebra 向公共数据库发出的查询只携带生物学信息（HPO 编号、基因、变异、疾病名），不含姓名、出生日期或病历号。mod 内置**隐私闸门**：外发调用若包含病例中登记的受保护标识，或身份证号、手机号、电子邮箱格式的内容，一律拒绝；病例目录内的文件或原始基因组文件（VCF/BAM/CRAM/FASTQ）外传前须经确认；若读不到病例的标识清单，闸门会关闭并拒绝所有外发调用。**闸门做不到的事**：你让 Claude 阅读的每一份病历（PDF、照片、报告）都会作为对话内容发送给模型服务商，这与在 Claude Code 中打开任何文件相同；闸门检查的是工具调用，不是对话本身。建议尽早登记受保护标识（`zebra case identifiers --add`），如有顾虑请在提供文件前自行脱敏。标识匹配已覆盖常见编码与写法变体，但属尽力而为，不构成保证。

## mod 为 Claude Code 增加的能力

- **15 个工具**，模型可直接调用（`mcp__zebra-mod__*`）：`case_status`、`case_update`、`hpo_search`、`phenotype_rank`、`gene_card`、`variant_card`、`disease_card`、`acmg`、`s2f_predict`、`therapy_landscape`、`trials_search`、`literature_search`、`rare_stats`、`edit_check`、`china_rare`。公共数据库查询已预先授权，只读研究操作无需逐次确认。
- **研究守则**写入系统提示词（即上述原则），并附当前病例信息。
- **病例看板**面板（`/zebra board`）与状态栏：表型、变异及其研究级分类、诊断假设、治疗线索、待解决问题、证据条数，随病例更新实时刷新。
- **`/zebra`** 命令：`new <目录> [标题]`、`case <目录>`、`board`、`ledger`、`doctor`、`close`。
- **12 个技能**（由 `/zebra-mod:zebra-start` 统一分流）：`zebra-safety`（急症红旗与特定罕见病的用药、麻醉、操作禁忌）、`zebra-intake`、`zebra-diagnose`、`zebra-variant`、`zebra-reanalysis`、`zebra-s2f`、`zebra-therapy`、`zebra-stats`、`zebra-literature`、`zebra-family`、`zebra-report`。
- **6 个子代理**：`phenotype-curator`、`variant-curator`、`s2f-analyst`、`therapy-scout`、`literature-scout`、`evidence-auditor`。
- **`zebra` 命令行工具**（仅依赖 Python 标准库，Python ≥ 3.9）：所有工具背后的唯一实现，可在任意终端、笔记本或其他 agent 中使用。

## 安装

前置条件：Claude Code ≥ 2.1.289（支持函数钩子 mod）；`PATH` 中有 Python ≥ 3.9（名为 `python3`，或在插件选项 `python` 中指定）。无需 pip 安装任何依赖。

```bash
# 从 GitHub 安装（私有仓库：需具备访问权限，并配置好 SSH 或 token）
claude plugin marketplace add zwbao/zebra-mod
claude plugin install zebra-mod@zebra-mod

# 或从本地目录加载，仅对当前会话生效
git clone https://github.com/zwbao/zebra-mod ~/zebra-mod
claude --plugin-dir ~/zebra-mod
```

可选，一次性下载约 80 MB 的 HPO 数据，以启用离线表型排序与即时术语检索：

```bash
~/zebra-mod/bin/zebra hpo fetch
```

可选的重型 S2F 模型：安装 [s2f-penguin](https://github.com/zwbao/s2f-penguin) 的 `s2f` 命令行（`uv tool install "git+https://github.com/zwbao/s2f-penguin"`），并设置 `ALPHAGENOME_API_KEY`（AlphaGenome：仅限非商业用途，不得用于临床决策）和/或 `NVCF_RUN_KEY`（NVIDIA 托管的 Evo 2）。可用 `/zebra doctor` 检查全部环境。

插件选项（`/config` → zebra-mod，或 settings 中的 `pluginConfigs`）：`python`（解释器）、`doctrine`（`always` | `case` | `off`）、`privacyGate`（默认开启）。

## 快速上手

```text
/zebra new ~/cases/lily "Lily：6 月龄起反复抽搐"
# 将检查报告、化验单、基因检测报告放入 ~/cases/lily/records/
/zebra-mod:zebra-intake
/zebra-mod:zebra-diagnose
/zebra-mod:zebra-variant   NM_001165963.4:c.2134C>T
/zebra-mod:zebra-therapy
/zebra-mod:zebra-report    中文家属信
```

也可以直接用自己的话提问，例如："我女儿 6 个月开始发热抽搐，基因报告显示 SCN1A 有一个变异，这意味着什么？"分流技能会接手后续流程。

## 数据来源与使用条款

zebra-mod 实时查询公共资源并保留出处，各资源适用其自身许可（见 [docs/DATA-SOURCES.md](docs/DATA-SOURCES.md)）。需特别注意：AlphaGenome 的输出仅限非商业用途且不得用于临床决策；SpliceAI 模型权重采用 CC BY-NC 许可；OMIM 内容不做转载分发。

## 非医疗器械

zebra-mod 提供研究级分析，帮助使用者提出更好的问题。它不做诊断、不开处方，也不能替代临床医生或遗传咨询师。

## 渊源

本项目吸收了此前两个项目的经验：g2t-harness（私有仓库——同心式 harness：模型不写排名、先澄清不编造、每个技能写明证伪条件）与 [s2f-penguin](https://github.com/zwbao/s2f-penguin)（带运行凭证的序列到功能模型、结论上限、以逻辑而非算术综合多轴证据）。重型序列模型由 zebra-mod 调用 s2f-penguin 的 `s2f` 命令行完成。

MIT 许可 © 2026 鲍志炜
