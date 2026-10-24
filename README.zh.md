<div align="center">

# 🦓 zebra-mod

**装进 Claude Code 的罕见病研究工作站**

从症状到鉴别诊断，从变异到 ACMG 分级，从序列到功能，从基因型到治疗——每一条结论都附有可以打开核对的来源。

[![License: MIT](https://img.shields.io/badge/license-MIT-0f5c6e.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.3.9-0f5c6e.svg)](CHANGELOG.md)
[![Claude Code](https://img.shields.io/badge/Claude%20Code-%E2%89%A5%202.1.289-1b2430.svg)](https://claude.com/claude-code)
[![Python](https://img.shields.io/badge/python-%E2%89%A5%203.9%20·%20stdlib%20only-1b2430.svg)](zebra/)

[English](README.md) · **中文**

</div>

> *"听到马蹄声，先想到马，而不是斑马。"* 医学教育让人优先考虑常见病；罕见病患者自称"斑马"，因为答案往往不在常规路径上。zebra-mod 专门往那里找，并且每一步都给出证据。

zebra-mod 是一个 [Claude Code](https://claude.com/claude-code) mod，面向**在求诊路上的患者与家庭、临床医生和研究者**。它为 Claude 提供 22 个工具，实时查询公共数据库（HPO、Monarch、Orphanet、ClinVar、gnomAD、Ensembl VEP、ClinGen、PanelApp、Open Targets、ClinicalTrials.gov、Europe PMC 等），提供带证据账本的病例工作区，以及诊断、变异解读、治疗检索和报告撰写的规范流程。你像平常一样用中文或英文与 Claude 对话，它会自行调用这些工具。

## 主要能力

- **从病历到鉴别诊断。** 门诊记录、化验单和照片整理为经核验的 HPO 表型；三套独立的排序（离线、Monarch、PubCaseFinder）并列呈现，与 Claude 自己的临床推理放在一起，而不是取而代之。
- **变异解读，ACMG 分值由代码计算。** VEP、gnomAD、ClinVar、ClinGen、SpliceAI、AlphaMissense、REVEL 与 MaveDB 提供证据；证据条目由你判断，分值由 zebra 计算（Tavtigian 2020 分值与 2015 组合规则并列）。
- **不止单核苷酸变异。** 拷贝数变异与芯片结果对照 ClinGen 剂量敏感性，外显子缺失给出读码框计算，SMN1/SMN2 拷贝数，重复扩增（D4Z4、DMPK、FMR1、HTT 等）按规则换算。
- **序列到功能模型在 mod 内运行。** SpliceAI、Pangolin、GPN-MSA、Evo 2 与 AlphaGenome 各自按所属维度报告，并标注结论上限；以一致性综合，从不取平均。
- **从基因型到治疗，并说明可及性。** 按机制检索已批准与在研药物，FDA/EMA 状态取自监管机构记录，国内获批与医保目录情况，有中国研究中心的在招试验，以及个体化治疗路径（反义寡核苷酸、碱基编辑）。
- **家属拿得到手的报告。** 临床医生摘要与通俗版家属资料，经独立子代理审查后，导出为中文排版的 HTML、Word 与 PDF。
- **隐私内建。** 病例文件只保存在本机；隐私闸门阻止已登记的姓名、出生日期和病历号出现在任何对外调用中。

## 快速开始

**1. 安装。** 在 Claude Code 中输入：

```text
安装 https://github.com/zwbao/zebra-mod
```

Claude 会按 [INSTALL.md](INSTALL.md) 完成其余步骤：检查 Claude Code（≥ 2.1.289）与 Python（≥ 3.9），安装插件，下载 HPO 数据（约 80 MB），并运行环境自检。也可以在终端执行：

```bash
curl -fsSL https://raw.githubusercontent.com/zwbao/zebra-mod/main/install.sh | sh
```

**2. 试用演示病例。** 新开一个会话，输入 `/zebra demo`。它会打开一个合成病例（一份门诊记录和一份基因报告），第一个问题已填好；按回车，即可看到它读病历、做鉴别诊断、解读变异的全过程。

**3. 问你自己的情况。** 无需记任何命令：

```text
我女儿6个月开始发热抽搐，基因报告说 SCN1A NM_001165963.4:c.2134C>T，
这是什么意思？下次见神经科医生该问什么？
```

`/zebra` 查看 mod 能做什么；`/zebra new <文件夹>` 为一位患者新建病例，把资料集中管理。

## 可以这样问

| 你是 | 示例 |
|---|---|
| **患者或家属** | "跑了好几家医院都没查出来，这是孩子的检查单，可能是什么病？" · "VUS 对我儿子意味着什么？" · "国内有没有 FSHD 的临床试验？" · "帮我列一页下次复诊要问的问题。" |
| **临床医生** | "按这些 HPO 术语给出鉴别诊断排序，以及每个候选该用什么检测确认。" · "结合这个家系的共分离，对该变异做 ACMG 解读。" · "DM1 合并传导阻滞患者的管理与随访方案。" · "血管型 EDS 患者能做导管血管造影吗？" |
| **研究者** | "结合表型重新分析这份三人家系 VCF。" · "对这个深部内含子变异跑 SpliceAI、Pangolin 和 AlphaGenome。" · "计算这些减数分裂的共分离似然比。" · "DMD 第 45–50 外显子缺失，跳读哪个外显子能恢复读码框？" |

每项能力的详细说明（按工具逐一列出）见 [docs/FEATURES.md](docs/FEATURES.md)。

## 工作原理

```
你
└─► Claude Code
      ├─ 研究准则（系统提示）· 13 个技能 · 6 个子代理
      └─► zebra-mod 钩子
            ├─ 隐私闸门
            └─► 22 个工具（mcp__zebra-mod__*）
                  └─► zebra 引擎（Python，仅标准库）
                        ├─► 公共数据库与模型：HPO · Monarch · ClinVar · gnomAD · VEP · ClinGen …
                        └─► 本机病例文件夹：case.json · evidence/ledger.jsonl（E1、E2…）· reports/
```

- **代码计算，模型解读。** 排序、ACMG 分值、重复单位数与统计量均由 `zebra` 引擎算出；模型从不编写或改动数字。当排序器看不到的证据才是决定因素时，模型会公开地推翻排序并说明理由。
- **每个来源都留档。** 每次工具结果连同数据库、记录号和检索时间写入病例的证据账本，报告引用这些来源。
- **欢迎 Claude 自身的知识，但要标注。** 本次检索所得与凭记忆的内容分开标示，读者能分辨哪些已核实。
- **mod 只做加法，不做减法。** 它不设任何先决关卡：先回答问题，病例整理等工作随后或同时进行。

同一个引擎也是命令行工具，可在任何终端、笔记本或其他智能体中使用：

```bash
zebra phenotype rank --present HP:0002373 HP:0001263 --sources local,monarch,pubcasefinder
zebra acmg suggest NM_001165963.4:c.2134C>T --inheritance AD
zebra acmg classify PVS1 PS2 PM2_Supporting
zebra s2f predict 12-25245350-C-T --models spliceai,pangolin,gpn_msa
zebra vcf triage trio.vcf.gz --proband P --mother M --father F --sex female --hpo-genes
zebra stats maxaf --prevalence 0.0000625 --allelic 0.05 --penetrance 0.9
```

## 效果如何

均有实测，方法与置信区间见 [docs/BENCHMARK.md](docs/BENCHMARK.md)：

- **表型排序**（GA4GH phenopacket-store，10,374 例已发表病例，按疾病划分开发集与留出集）。在"本身论文未被用于构建 HPO 注释"的留出病例上（这才是对新患者的公平检验），正确疾病进入离线排序前 10 的比例为 **28%**，致病基因为 **40%**。排序只是待验证的假设清单，所以 Claude 自己的鉴别诊断与它并列。
- **装与不装 mod 的对比**（6 个留出的端到端病例，由未接触 mod 说明的独立代理编写，每组 8 次，按指南逐项评分）。合并后两组持平（差值 −0.02，95% 置信区间 −0.16 至 +0.12）：长病历中夹带急症、VUS 共分离两例中装 mod 更好，另有两例更差。"装了不会更差"尚未得到证实，基准文档如实记录。
- [`evals/`](evals/) 中有 **19 个端到端评测**，用 `claude plugin eval` 运行，涵盖从中文家长上传的病历，到应按肿瘤学而非罕见病思路解读的体细胞肿瘤报告。

## 隐私与安全

- **只留本机：** 病例文件、病历和证据账本只写在你的电脑上。
- **zebra 从不发送：** 姓名、出生日期和病历号。数据库查询只携带生物学信息（HPO 编号、基因、变异）；外发调用若含已登记的标识，或身份证号、手机号、电子邮箱格式的内容，隐私闸门会拒绝；病例文件或原始基因组文件外传前会先询问。
- **闸门做不到的事：** 你让 Claude 阅读的每份病历都会作为对话内容发送给模型服务商，与在 Claude Code 中打开任何文件相同。如有顾虑，请先自行脱敏。
- **急症优先：** 需要当天就医的症状，会在任何分析之前先说明。
- **研究级：** zebra-mod 产出供医生审阅的分析草稿。它不是医疗器械，不自行做出诊断或处方，不能替代临床医生或遗传咨询师。

## 配置

在 `/config` → zebra-mod 中设置（或 settings 中的 `pluginConfigs`）：

| 选项 | 默认 | 作用 |
|---|---|---|
| `doctrine` | `auto` | 打开病例时启用完整研究准则，否则只在罕见病问题上生效的简版；也可设为 `always`、`case`、`off` |
| `privacyGate` | 开 | 上文所述的外发调用检查 |
| `interface` | `full` | 工具行、提示文字与奔跑的斑马；`quiet` 不显示动画；`off` 只用 Claude Code 自身的界面 |
| `python` | `python3` | 引擎使用的 Python 解释器 |

可选的序列模型（`/zebra doctor` 会逐项检查）：

| 模型 | 需要 |
|---|---|
| SpliceAI、Pangolin、GPN-MSA | 无需设置 |
| Evo 2 | `NVCF_RUN_KEY`（NVIDIA 托管服务） |
| AlphaGenome | `ALPHAGENOME_API_KEY` 与 `pip install alphagenome`（或设置 `ZEBRA_ALPHAGENOME_PYTHON`）；仅限非商业用途，不得用于临床决策 |

## 文档

| | |
|---|---|
| [docs/FEATURES.md](docs/FEATURES.md) | 全部能力、工具、技能与子代理 |
| [docs/BENCHMARK.md](docs/BENCHMARK.md) | 排序基准与装/不装 mod 的对比评测 |
| [docs/DATA-SOURCES.md](docs/DATA-SOURCES.md) | 各数据库与模型的使用条款和许可 |
| [docs/PRD.md](docs/PRD.md) · [docs/ROADMAP.md](docs/ROADMAP.md) | 设计与后续计划 |
| [INSTALL.md](INSTALL.md) · [CHANGELOG.md](CHANGELOG.md) | 安装步骤与版本说明 |

## 开发

```
hooks/        mod 本体：工具、准则、隐私闸门、病例面板、/zebra      （TypeScript）
skills/       13 个流程，由 zebra-start 分派
agents/       6 个子代理（整理、检索、证据审查）
zebra/        引擎：数据源、ACMG、统计、VCF、S2F、报告导出          （Python，仅标准库）
evals/        `claude plugin eval` 端到端评测
tests/        pytest（默认离线；设置 ZEBRA_NETWORK_TESTS=1 运行联网测试）
```

```bash
python3 -m pytest -q            # 引擎
claude plugin test .            # mod 钩子
claude plugin validate .        # 清单、技能、子代理
claude --plugin-dir .           # 用当前代码运行一次会话
```

欢迎提交 issue 和 pull request。修改准则或技能时，请附上装/不装 mod 的对比评测（`claude plugin eval --ablation with-without`）；`evals/heldout-*` 中的留出评分标准不得为迎合结果而修改。

## 数据来源与许可

zebra-mod 实时查询公共资源并记录每个答案的出处，各资源保留其自身条款（见 [docs/DATA-SOURCES.md](docs/DATA-SOURCES.md)）。特别注意：AlphaGenome 输出仅限非商业用途且不得用于临床决策；SpliceAI 权重为 CC BY-NC；不转载 OMIM 内容。

## 引用

```bibtex
@software{zebra_mod,
  author  = {zwbao},
  title   = {zebra-mod: a rare-disease research workstation for Claude Code},
  year    = {2026},
  version = {0.3.9},
  url     = {https://github.com/zwbao/zebra-mod}
}
```

## 致谢

本项目建立在 HPO、Monarch、Orphanet、ClinGen、ClinVar、gnomAD、Ensembl、PanelApp、Open Targets、Europe PMC、PubTator、LitVar 与 PubCaseFinder 团队的工作之上，也感谢 SpliceAI、Pangolin、GPN-MSA、Evo 2 与 AlphaGenome 的作者。输入框上方奔跑的斑马，取自埃德沃德·迈布里奇 1878 年的连续摄影《运动中的马》（公有领域）。

MIT 许可 © 2026 zwbao。
