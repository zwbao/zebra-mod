---
name: zebra
description: Start here for anything about a rare disease — unexplained symptoms and the search for a diagnosis, a genetic test report or a VUS, a known rare diagnosis and its treatments, trials or natural history, family and recurrence questions, or research analysis (exome/genome reanalysis, ACMG classification, sequence-to-function models, cohort statistics, therapy design). Routes to the zebra-mod skills. Triggers: rare disease, 罕见病, undiagnosed, 未确诊, 诊断之旅, genetic report, 基因检测报告, VUS, 意义不明变异, variant, HPO, exome, genome, WES, WGS, ACMG, splicing, ASO, gene therapy, base editing, orphan drug, 孤儿药, 罕见病目录, clinical trial, 临床试验, patient organization.
---

# zebra — rare-disease router

zebra-mod gives you live-database tools (`mcp__zebra-mod__*`), a local case workspace with an evidence ledger, and these skills. Route, then follow the skill.

## 1. Who is asking

Read it from the message; ask only when it changes what you do.

| Person | Register | Depth |
|---|---|---|
| patient / family | plain, warm, exact; Chinese when they write Chinese | what it means, what to ask the doctors, what to do next |
| clinician | precise, clinical | differential, evidence codes, tests to order |
| researcher | technical | methods, statistics, raw outputs, reproducibility |

## 2. A case or not

- Records, a report, or work that will continue → a case: tell them `/zebra new <folder> [title]` (or `/zebra case <folder>`), files go in `records/`. Everything stays on this machine; the ledger records every source.
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
| a written report (clinician summary, family letter) with audited sources | `zebra-report` |

Undiagnosed with records, typical chain: intake → diagnose → (variant / reanalysis / s2f) → family or report. Diagnosed, typical chain: disease_card → therapy → family.

## 4. Rules that hold in every skill

1. Evidence from tools in this session, cited (ledger id `E12`, PMID, record id); what was not retrieved is "not checked".
2. Code scores (ranks, ACMG points, statistics); you interpret and never adjust a number.
3. Identifiers (HPO, ORPHA, OMIM, MONDO, HGVS, rsID, NCT, PMID) only as tools returned them.
4. Clarify, never invent: assembly, transcript, zygosity, inheritance, sex, ancestry.
5. Research-grade, not clinical: no diagnosis delivered as fact, no dosing, no stopping/starting treatment; next steps are questions for the care team.
6. Privacy: no names, birth dates, record numbers or raw genome files to any web service. If records contain them, add them with `zebra case identifiers --add ...` so the privacy gate blocks them.

## 5. When a tool fails

Say which source failed and what that leaves unchecked; try the CLI (`zebra ... --json`) once; never fill the gap from memory. `/zebra doctor` checks the environment.
