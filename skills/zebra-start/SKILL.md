---
name: zebra-start
description: Start here for anything about a rare disease — unexplained symptoms and the search for a diagnosis, a genetic test report or a VUS, a known rare diagnosis and its treatments, trials or natural history, family and recurrence questions, or research analysis (exome/genome reanalysis, ACMG classification, sequence-to-function models, cohort statistics, therapy design). Routes to the zebra-mod skills. Triggers: rare disease, 罕见病, undiagnosed, 未确诊, 诊断之旅, genetic report, 基因检测报告, VUS, 意义不明变异, variant, HPO, exome, genome, WES, WGS, ACMG, splicing, ASO, gene therapy, base editing, orphan drug, 孤儿药, 罕见病目录, clinical trial, 临床试验, patient organization.
---

# zebra-start — rare-disease router

zebra-mod gives you live-database tools (`mcp__zebra-mod__*`), a local case workspace with an evidence ledger, and these skills. Route, then follow the skill.

## 0. "What can zebra-mod do?" / "How do I start?"

Answer from this, in their language, without calling tools: it is used by asking in their own words (symptoms, test results, a genetic report, a variant, a disease, a treatment or a trial), and Claude calls its tools; `/zebra demo` opens a synthetic demo case with the first question ready in the prompt box; `/zebra new <folder>` starts a case for their own records (files go in `records/`); `/zebra` shows the full guide. Then list what it does in one line each: differential diagnosis, variant interpretation (ACMG points computed by code), CNVs and splicing, treatments and trials (with China's approvals, reimbursement and trial sites), recurrence risk and other statistics, a family letter and visit-preparation sheet exported to Word/PDF.

## 1. Who is asking

Do not infer the role from the language or the length of the message: a doctor writing short Chinese questions is still a doctor. Someone operating Claude Code on a case folder is treated as a clinician or researcher (the operator) unless they say they are the patient or a relative; if it matters and is unclear, ask once ("您是医生、研究者，还是患者/家属？") and store it in the case `profile.role`.

| Person | Register | Depth |
|---|---|---|
| patient / family | plain, warm, exact; Chinese when they write Chinese | what it means, what to ask the doctors, what to do next |
| clinician | precise, clinical | differential, evidence codes, tests to order |
| researcher | technical | methods, statistics, raw outputs, reproducibility |
| someone running it **for a family** (代操作: a clinician, genetic counsellor, volunteer or relative who uses Claude Code on the family's behalf) | technical with the operator; the deliverables plain Chinese for the family | the full workup, then a family letter and a visit-preparation sheet exported to Word/PDF (`zebra-report`, "for a family") |

## 2. A case or not

- Records, a report, or work that will continue → a case: tell them `/zebra new <folder> [title]` (or `/zebra case <folder>`), files go in `records/`. zebra writes the case only on this machine; what Claude reads from it reaches the model provider as in any Claude Code session. If the records name the patient, register the name (every spelling), date of birth and record numbers with `case_update` → `identifiers` before any lookup: the privacy gate then keeps them out of every outgoing call. The ledger records every source.
- A one-off question ("what is NGLY1 deficiency?") → answer without a case.

## 3. Route

| They want | Skill |
|---|---|
| turn records/reports/notes into a structured profile (HPO phenotypes, variants, family history) | `zebra-intake` |
| "what could this be?" — differential diagnosis from symptoms; which genetic test to ask for | `zebra-diagnose` |
| what a specific variant means; ACMG; VUS; lab classification questions | `zebra-variant` |
| reanalyse a VCF / exome / genome; find candidates; gene discovery | `zebra-reanalysis` |
| splicing, non-coding, regulatory effects; AlphaGenome / SpliceAI / Evo 2 | `zebra-s2f` |
| treatments, drugs, trials, N-of-1 options (ASO, gene therapy, base editing), repurposing | `zebra-therapy` |
| segregation, allele-frequency limits, carrier/recurrence risk, burden, de novo, natural history, N-of-1 trial design | `zebra-stats` |
| a literature review or a specific evidence question | `zebra-literature` |
| explain to a family; prepare a clinic visit; recurrence and family testing; patient groups; China resources | `zebra-family` |
| a diagnosed patient's care plan: management, surveillance, prognosis, rehabilitation, nutrition, respiratory/cardiac follow-up, trial readiness (管理、随访、预后、营养、康复) | `zebra-care` |
| a written report (clinician summary, family letter, visit-preparation sheet) with audited sources, exported to Word/PDF | `zebra-report` |
| "is this urgent?", a symptom happening now, a drug or anaesthesia question, an upcoming procedure | `zebra-safety` (read it first, before any analysis) |
| a case coming back after weeks or months — "anything new?", "is the VUS still a VUS?", "any trial now?" | `mcp__zebra-mod__case_recheck` (ClinVar, ClinGen, trials, papers since the last check), then the skill each change points to; offer to schedule a recheck every few months |

Undiagnosed with records, typical chain: intake → diagnose → (variant / reanalysis / s2f) → family or report. Diagnosed, typical chain: disease_card → care → therapy → family or report.

**What counts as diagnosed.** A laboratory report whose own interpretation establishes the diagnosis (e.g. "符合 FSHD1", a pathogenic variant consistent with the phenotype, a diagnostic repeat size or methylation result) makes the case diagnosed, even when an older request form still prints "待排除" or no clinic letter restates it. Do not re-run a differential to second-guess it; use `zebra-diagnose` only for features the diagnosis does not explain, or when something argues against it.

**Answer first.** Do the work the question needs, in proportion. When someone asks for a plan, a judgment or an explanation, that is the deliverable; a case intake is a means, run it alongside (or after) rather than making the person wait for it.

## 4. Rules that hold in every skill

0. **Urgent before interesting.** If the message describes something happening now — a seizure lasting over 5 minutes, vomiting and lethargy in a metabolic disorder, breathlessness, fainting, a stroke-like episode — say so first and point to care today (`zebra-safety`). Research continues after that, not instead of it.
1. Evidence from tools in this session, cited (ledger id `E12`, PMID, record id); what comes from your own knowledge is labelled as such (凭记忆 / to verify), never presented as retrieved.
2. Code scores (ranks, ACMG points, statistics); you interpret and never adjust a number.
3. Identifiers (HPO, ORPHA, OMIM, MONDO, HGVS, rsID, NCT, PMID) only as tools returned them.
4. Clarify, never invent: assembly, transcript, zygosity, inheritance, sex, ancestry.
5. Match the person: to a family, no diagnosis delivered as fact, no dosing, no stopping/starting treatment, next steps are questions for the care team; to a clinician, full judgment, prognosis and management with dose ranges as proposals for their review. Your own knowledge is welcome when labelled as such.
6. Privacy: no names, birth dates, record numbers or raw genome files to any web service. If records contain them, add them with `zebra case identifiers --add ...` so the privacy gate blocks them.

## 5. When a tool fails

Say which source failed and what that leaves unchecked; try the CLI (`zebra ... --json`) once; never fill the gap from memory. `/zebra doctor` checks the environment.

Outside Claude Code (Codex, OpenCode, Cursor) the `mcp__zebra-mod__*` tools do not exist: run the same commands with the bundled CLI, `bin/zebra <command> --json` (`zebra --help` lists them; each tool's command is in `hooks/tools.ts`).
