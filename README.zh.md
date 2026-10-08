# 🦓 zebra-mod

**zebra-mod 是一个 Claude Code mod，把 Claude Code 改造成罕见病研究与分析工作台**，服务对象包括正在寻求诊断的患者和家属、临床医生与科研人员。

> 医学教育有一句格言："听到马蹄声，先想到马，而不是斑马。"罕见病患者自称"斑马"，因为他们的答案往往不在常见路径上。zebra-mod 的目标正是在这些路径之外寻找答案，并为每一步给出可核查的证据。

[English → README.md](README.md)

---

## 功能

| 任务 | 做法 |
|---|---|
| **鉴别诊断**：根据表型推断可能的疾病 | 病历 → HPO 表型档案（术语经核验，区分存在/排除，记录起病时间与出处页码）→ 三种相互独立的排序方法（基于 HPO 注释的离线 Resnik 相似度、Monarch 语义相似度、PubCaseFinder），各自保留排名，以方法间的一致性作为信号 → 鉴别特征表 → 每个候选疾病应采用何种检测（外显子组、基因组、染色体微阵列、重复扩增检测、甲基化、线粒体 DNA、代谢筛查） |
| **变异解读**：ACMG/AMP 框架 | VEP（MANE 转录本、AlphaMissense、REVEL、CADD）+ ±4,999 nt 范围的 SpliceAI + gnomAD（grpmax、faf95；线粒体变异的异质性计数）+ ClinVar（四星制审核星级）+ ClinGen 基因-疾病有效性与剂量敏感性 + MaveDB 功能实验分数（按数据集自身校准）+ LitVar（排除同一位点的其他等位变异）+ 中国人群频率 → 按 ClinGen 校准阈值给出数据驱动的证据条目，以及 PVS1、PS1/PM5 的**判断依据** → 由模型论证的判断类条目 → **由代码计算积分**（Tavtigian 2020），并与 2015 年组合规则的结论并列呈现；GRCh37 输入自动换算至 GRCh38 |
| **序列到功能（S2F）**：剪接、非编码与调控变异 | SpliceAI/Pangolin（Broad 在线查询），AlphaGenome、Evo 2 与 GPN-MSA 通过 [s2f-penguin](https://github.com/zwbao/s2f-penguin) 的 `s2f` 命令行调用；逐轴解读并标注结论上限，以逻辑一致性综合，不做数值平均；由机制推导验证方案（RNA 检测，依据 GTEx 选择该基因有表达的组织）；剪接调控型反义寡核苷酸初筛（`aso_screen`：在患者序列上核对异常剪接事件、列出候选靶向窗口与已发表的 N-of-1 先例） |
| **外显子组/基因组重分析**，在本机完成 | `zebra qc`：性别核对、KING 亲缘系数（样本对调、非亲生）、纯合片段、孟德尔错误与单亲二体、嵌合新发变异。VCF 分诊：质量过滤、遗传模式（新发、纯合、复合杂合、X 连锁半合子）、支持 PED 家系文件、按表型基因限定范围、对靠前的剪接/非编码候选运行 SpliceAI；全外显子组可经 MyVariant 频率预过滤。VCF 文件留在本机；变异位点会发往注释服务——默认只发过滤后的候选，开启预过滤则发送全部位点（会先征得同意） |
| **基因型到治疗（G2T）** | 先判定致病机制（功能丧失、功能获得、显性负效应、剪接异常、重复扩增）→ 已上市与在研药物（Open Targets / ChEMBL）、**可及性**：FDA 与 EMA 以监管机构自身记录为准，国内批准情况依据国家药监局/药审中心官方文件，2025 版国家医保目录含限定支付范围原文，以及在中国设有中心的临床试验；文献证据、个体化治疗可行性初筛（反义寡核苷酸、基因替代、**碱基编辑可行性**）→ 按 A–E 分级，并逐条核对机制方向 |
| **经典统计方法** | 共分离似然比与 PP1 强度、最大可信等位基因频率（BS1）、携带者频率与遗传患病率、再发风险（含 X 连锁贝叶斯计算）、Fisher 精确检验与负荷检验、新发变异富集、Kaplan–Meier 自然史分析与对数秩检验、N-of-1 试验设计与分析 |
| **文献** | Europe PMC、PubTator3、LitVar2；PMID 仅使用工具返回的结果，结论引用原文语句 |
| **面向家庭** | 通俗书面语解释（中文用户默认中文）、如实说明意义不明变异（VUS）、复诊问题清单、再发风险与家系筛查、患者组织、国家罕见病目录及各省罕见病诊疗协作网医院 |
| **长期随访** | `case_recheck` 重新查询病例的关键问题——各变异的 ClinVar 分类与星级、相关基因的 ClinGen 有效性、各待定诊断的在招试验与新发表文献——报告自上次复查以来的变化，并记入病例时间线 |
| **报告** | 临床医生摘要、家属信与就诊准备单，每条结论附出处，定稿前由**独立的对抗式审查**核对，再导出为中文排版的 **Word 与 PDF**——医生或志愿者可替家庭完成分析，把文件交给家属 |

## 原则（由机制保证，而非仅写在文档中）

1. **证据优先，不凭记忆。** 每个结论引用本次会话中检索到的证据账本编号、PMID 或数据库记录。每个病例维护一份只追加的**证据账本**（`E1、E2……`），记录用到的全部来源。
2. **代码计分，模型解读。** 排名、ACMG 积分与统计量均由 `zebra` 引擎计算；模型不编写、不修改任何数值。
3. **信息缺失时先澄清，绝不编造。** 基因组版本、转录本、合子状态、遗传方式、性别等关键信息缺失时，先询问，或以条件句给出结论。
4. **研究级分析，而非临床报告。** 不以确定口吻告知诊断，不提供用药剂量；后续步骤以"供医疗团队考虑的问题"形式呈现。
5. **隐私：如实说明边界。** 病例文件只写在本机；zebra 向公共数据库发出的查询只携带生物学信息（HPO 编号、基因、变异、疾病名），不含姓名、出生日期或病历号。mod 内置**隐私闸门**：外发调用若包含病例中登记的受保护标识，或身份证号、手机号、电子邮箱格式的内容，一律拒绝；病例目录内的文件或原始基因组文件（VCF/BAM/CRAM/FASTQ）外传前须经确认；若读不到病例的标识清单，闸门会关闭并拒绝所有外发调用。**闸门做不到的事**：你让 Claude 阅读的每一份病历（PDF、照片、报告）都会作为对话内容发送给模型服务商，这与在 Claude Code 中打开任何文件相同；闸门检查的是工具调用，不是对话本身。建议尽早登记受保护标识（病历中出现患者姓名时，模型会通过 `case_update` 的 `identifiers` 登记），如有顾虑请在提供文件前自行脱敏。导出 Word/PDF 时，报告中若仍含已登记的标识或身份证号，导出会被拒绝。标识匹配已覆盖常见编码与写法变体，但属尽力而为，不构成保证。

## 效果如何

以测量结果为准，不做宣称——方法、数据校验和与置信区间见 [docs/BENCHMARK.md](docs/BENCHMARK.md)。

- **表型排序**：在 GA4GH phenopacket-store 0.1.27（10,374 个已发表病例、780 种疾病）上按疾病划分开发集与留出集。对留出集中"病例本身的论文不是 HPO 注释来源"的病例（最接近新患者的公平测试），正确疾病进入本地排序前 10 的比例为 **28%**（0.1.0 为 15%），致病基因为 40%。全部留出病例为 70%，只能视为上限：多数已发表病例的疾病注释正是从该论文整理而来。在 100 例联网抽样中，本地 + Monarch + PubCaseFinder 合并后前 10 命中 75%。排序结果是待检验的假设清单，不是答案。
- 一项如实披露的代价：PRD 中的示例查询（热性惊厥、局灶性与强直阵挛发作、发育落后；排除肌张力低下）在本地排序中，Dravet 综合征现排第 16（0.1.0 为第 6），原因是排除表型改为标注而不计分——这一选择在留出集上测得更好。
- `evals/` 收录 10 个端到端评测用例，供 `claude plugin eval` 使用（中文家长病历、急症优先、VUS、CNV、线粒体变异、治疗与试验、隐私、超出范围的请求）。

## mod 为 Claude Code 增加的能力

- **21 个工具**，模型可直接调用（`mcp__zebra-mod__*`）：`case_status`、`case_update`、`case_recheck`、`hpo_search`、`phenotype_rank`、`gene_card`、`variant_card`、`disease_card`、`acmg`、`cnv_interpret`、`s2f_predict`、`therapy_landscape`、`trials_search`、`literature_search`、`rare_stats`、`edit_check`、`china_rare`、`access`、`expression`、`aso_screen`、`report_export`。每次调用都经过你的权限规则与隐私闸门；只读查询默认无需确认（除非你的规则另有要求），写入病例遵循当前权限模式（可选择"本次会话内允许"）。
- **研究守则**写入系统提示词（即上述原则），并附当前病例信息。
- **看得见它在工作**（全部为增量挂载，与罕见病无关的工作照常显示）：它的工具调用有自己的行（`🦓 变异卡  NM_001165963.4:c.2134C>T`，随后 `⎿ 证据 15 条（E39–E53）· 来源 Ensembl VEP、gnomAD、ClinVar、LitVar2 · 3 条来自缓存`），取代原始 JSON；登记身份信息的 case_update 行只显示条数，不显示内容；等待动画写明正在查询哪些数据库；查询进行时，输入框上方有一匹奔跑的小斑马，隐私闸门拦下调用后，这里改为盾牌提示；底部状态区有 `🦓` 标签（打开病例时显示病例标题）；用到 zebra-mod 的一轮结束时，留下一行小结：调用次数、数据源、新增证据、被拦截的调用。中英文随病例语言或你的提问切换。选项 `interface`：`full`（默认）、`quiet`（无动画、不改等待动画文字）或 `off`。
- **病例看板**面板（`/zebra board`）与状态栏：表型、变异及其研究级分类、诊断假设、治疗线索、待解决问题、证据条数，随病例更新实时刷新。
- **`/zebra`** 命令：`new <目录> [标题]`、`case <目录>`、`board`、`ledger`、`doctor`、`close`。
- **12 个技能**（由 `/zebra-mod:zebra-start` 统一分流）：`zebra-safety`（急症红旗与特定罕见病的用药、麻醉、操作禁忌）、`zebra-intake`、`zebra-diagnose`、`zebra-variant`、`zebra-reanalysis`、`zebra-s2f`、`zebra-therapy`、`zebra-stats`、`zebra-literature`、`zebra-family`、`zebra-report`。
- **6 个子代理**：`phenotype-curator`、`variant-curator`、`s2f-analyst`、`therapy-scout`、`literature-scout`、`evidence-auditor`。
- **`zebra` 命令行工具**（仅依赖 Python 标准库，Python ≥ 3.9）：所有工具背后的唯一实现，可在任意终端、笔记本或其他 agent 中使用。

## 安装

前置条件：Claude Code ≥ 2.1.289（支持函数钩子 mod）；Python ≥ 3.9。无需 sudo，无需 pip 安装任何依赖。

**一句话安装。** 在 Claude Code 里说：

```text
帮我安装 https://github.com/zwbao/zebra-mod
```

Claude Code 会读取仓库中的 [INSTALL.md](INSTALL.md) 并逐步完成：检查版本与 Python、注册并安装插件、下载 HPO 数据（约 80 MB）、运行 `zebra doctor` 自检。装好后在新会话中输入 `/zebra-mod:zebra-start` 即可开始（已打开的会话可输入 `/reload-plugins`）。

**一条命令安装**（效果相同）：

```bash
curl -fsSL https://raw.githubusercontent.com/zwbao/zebra-mod/main/install.sh | sh
# 仓库为私有时，先克隆再运行：
git clone https://github.com/zwbao/zebra-mod ~/zebra-mod && ~/zebra-mod/install.sh
```

**手动安装**：

```bash
claude plugin marketplace add zwbao/zebra-mod
claude plugin install zebra-mod@zebra-mod
python3 <安装目录>/bin/zebra hpo fetch     # 可选：离线表型排序与即时术语检索
# 仅对当前会话试用：claude --plugin-dir ~/zebra-mod
```

可选的重型 S2F 模型：安装 [s2f-penguin](https://github.com/zwbao/s2f-penguin) 的 `s2f` 命令行（`uv tool install "git+https://github.com/zwbao/s2f-penguin"`），并设置 `ALPHAGENOME_API_KEY`（AlphaGenome：仅限非商业用途，不得用于临床决策）和/或 `NVCF_RUN_KEY`（NVIDIA 托管的 Evo 2）。可用 `/zebra doctor` 检查全部环境。

插件选项（`/config` → zebra-mod，或 settings 中的 `pluginConfigs`）：`python`（解释器）、`doctrine`（默认 `auto`：打开病例时注入完整规则，否则只放一段仅对罕见病问题生效的简短说明，你在 Claude Code 里的其它工作不受影响；也可选 `always` | `case` | `off`）、`privacyGate`（默认开启；未打开病例时，其它工具调用中出现邮箱或手机号会先询问，而不直接拒绝）、`interface`（默认 `full`：工具专属行、等待动画文字、奔跑的斑马与每轮小结；`quiet` 不显示动画和等待动画文字；`off` 完全使用 Claude Code 自己的显示）。

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

本项目的序列到功能（S2F）能力建立在 [s2f-penguin](https://github.com/zwbao/s2f-penguin) 之上：带运行凭证的序列到功能模型、结论上限、以逻辑而非算术综合多轴证据。重型序列模型由 zebra-mod 调用 s2f-penguin 的 `s2f` 命令行完成。

MIT 许可 © 2026 zwbao
