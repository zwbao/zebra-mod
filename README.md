# 🦓 zebra-mod

**A Claude Code mod that turns Claude Code into a rare-disease research workstation** — for patients and families on a diagnostic odyssey, for clinicians, and for researchers.

> *"When you hear hoofbeats, think horses, not zebras."* Medicine is taught to expect the common. People with rare diseases call themselves zebras because their answers lie off the beaten path. zebra-mod is built to look there — and to show its evidence for every step.

[中文说明 → README.zh.md](README.zh.md)

---

## What it does

| Job | How zebra-mod does it |
|---|---|
| **"What could this be?"** — phenotype-driven differential diagnosis | Records → HPO profile (verified terms, present/excluded, onset, source page) → three independent rankings (offline Resnik over HPO annotations, Monarch semantic similarity, PubCaseFinder) kept separate, agreement as signal → discriminating features → which test finds each candidate (exome, genome, CMA, repeat assays, methylation, mtDNA, metabolic) |
| **"What does this variant mean?"** — ACMG/AMP interpretation | VEP (MANE, AlphaMissense, REVEL, CADD, SpliceAI) + gnomAD (grpmax, faf95) + ClinVar (review stars) + ClinGen validity and dosage + LitVar → data-driven codes with ClinGen-calibrated thresholds → your judged codes → **points computed by code** (Tavtigian 2020) and the 2015 combining rules side by side |
| **Sequence-to-function (S2F)** — splicing, non-coding, regulatory | SpliceAI/Pangolin (Broad lookup), AlphaGenome, Evo 2 and GPN-MSA through the [s2f-penguin](https://github.com/zwbao/s2f-penguin) `s2f` CLI; axis by axis with claim ceilings, combined by agreement, never averaged; mechanism → how to confirm with RNA |
| **Exome/genome reanalysis** on your own machine | VCF triage: quality, inheritance models (de novo, homozygous, compound het, X-hemizygous), phenotype-gene restriction, VEP annotation of survivors only — the VCF never leaves the computer |
| **Genotype → therapy (G2T)** | Mechanism first (LoF / GoF / DN / splicing / repeat) → approved and investigational drugs (Open Targets / ChEMBL), recruiting trials (ClinicalTrials.gov, incl. China), literature, N-of-1 screens (antisense, gene replacement, **base editing feasibility**) → leads tiered A–E with a mechanism-fit check |
| **Classical statistics** | Cosegregation LR → PP1, maximum credible allele frequency (BS1), carrier frequency and genetic prevalence, recurrence risk (incl. Bayesian X-linked), Fisher/burden, de novo enrichment, Kaplan–Meier natural history + log-rank, N-of-1 trial design and analysis |
| **Literature** | Europe PMC, PubTator3, LitVar2 — PMIDs only as returned, findings quoted from the paper |
| **For families** | Plain-language explanations (Chinese by default for Chinese speakers), VUS explained honestly, visit preparation, recurrence and cascade testing, patient organisations, China's national rare disease lists |
| **Reports** | Clinician summary and family letter, every claim cited, then an **independent adversarial audit** before it is final |

## Principles (enforced, not just written)

1. **Evidence, not recall.** Every answer cites a ledger id, PMID or database record retrieved in the session. Each case keeps an append-only **evidence ledger** (`E1, E2, …`) of every source used.
2. **Code scores, the model interprets.** Rankings, ACMG points and statistics are computed by the `zebra` engine; the model never writes or adjusts a number.
3. **Clarify, never invent.** Assembly, transcript, zygosity, inheritance, sex — missing and material means ask, or conclude conditionally.
4. **Research-grade, not clinical.** No diagnosis delivered as fact, no dosing; next steps are questions for the care team.
5. **Privacy by construction.** Case files stay local. A **privacy gate** in the mod refuses any outgoing call carrying a case's protected identifiers (names, birth dates, record numbers) or ID-number/phone/email patterns, and asks before raw genome files (VCF/BAM/CRAM/FASTQ) leave the machine — even in bypass-permissions mode.

## What the mod adds to Claude Code

- **15 tools** the model calls directly (`mcp__zebra-mod__*`): `case_status`, `case_update`, `hpo_search`, `phenotype_rank`, `gene_card`, `variant_card`, `disease_card`, `acmg`, `s2f_predict`, `therapy_landscape`, `trials_search`, `literature_search`, `rare_stats`, `edit_check`, `china_rare` — public-database lookups are pre-approved; no permission prompts for read-only research.
- **A research doctrine** in the system prompt (the rules above), with the active case.
- **A case board** pane (`/zebra board`) and a status line: phenotypes, variants and their research class, hypotheses, therapy leads, open questions, evidence count — live as the case changes.
- **`/zebra`** command: `new <dir> [title]`, `case <dir>`, `board`, `ledger`, `doctor`, `close`.
- **11 skills** (`/zebra-mod:zebra` routes): `zebra-intake`, `zebra-diagnose`, `zebra-variant`, `zebra-reanalysis`, `zebra-s2f`, `zebra-therapy`, `zebra-stats`, `zebra-literature`, `zebra-family`, `zebra-report`.
- **6 subagents**: `phenotype-curator`, `variant-curator`, `s2f-analyst`, `therapy-scout`, `literature-scout`, `evidence-auditor`.
- **The `zebra` CLI** (Python standard library only, Python ≥ 3.9): the single implementation behind every tool, usable from any shell, notebook or other agent.

## Install

Requirements: Claude Code ≥ 2.1.289 (function-hook mods), Python ≥ 3.9 on `PATH` as `python3` (or set the plugin's `python` option). No pip installs.

```bash
# from GitHub (needs access to the repository)
claude plugin marketplace add zwbao/zebra-mod
claude plugin install zebra-mod@zebra-mod

# or from a local checkout, for one session
git clone https://github.com/zwbao/zebra-mod ~/zebra-mod
claude --plugin-dir ~/zebra-mod
```

Optional, once (≈ 80 MB, enables offline phenotype ranking and instant HPO search):

```bash
~/zebra-mod/bin/zebra hpo fetch
```

Optional heavy S2F models: install the `s2f` CLI from [s2f-penguin](https://github.com/zwbao/s2f-penguin) (`uv tool install "git+https://github.com/zwbao/s2f-penguin"`), and set `ALPHAGENOME_API_KEY` (AlphaGenome: non-commercial, not for clinical decisions) and/or `NVCF_RUN_KEY` (Evo 2 on NVIDIA). Check everything with `/zebra doctor`.

Options (`/config` → zebra-mod, or `pluginConfigs` in settings): `python` (interpreter), `doctrine` (`always` | `case` | `off`), `privacyGate` (on by default).

## Quick start

```text
/zebra new ~/cases/lily "Lily — seizures since 6 months"
# put reports, lab sheets, the genetic report in ~/cases/lily/records/
/zebra-mod:zebra-intake
/zebra-mod:zebra-diagnose
/zebra-mod:zebra-variant   NM_001165963.4:c.2134C>T
/zebra-mod:zebra-therapy
/zebra-mod:zebra-report    family letter in Chinese
```

Or just ask in your own words — "我女儿6个月开始发热抽搐，基因报告说SCN1A有个变异，这是什么意思？" — the router skill takes it from there.

Researchers can drive the engine directly:

```bash
zebra hpo rank HP:0002373 HP:0007359 HP:0002133 HP:0001263 --exclude HP:0001252
zebra phenotype rank --present HP:0002373 HP:0001263 --sources local,monarch,pubcasefinder
zebra variant NM_000492.4:c.1652G>A --json
zebra acmg suggest NM_001165963.4:c.2134C>T --inheritance AD
zebra acmg classify PVS1 PS2 PM2_Supporting
zebra s2f predict 7-117559590-ATCT-A --models spliceai,pangolin
zebra vcf triage trio.vcf.gz --proband P --mother M --father F --sex female --hpo-genes --case ~/cases/lily
zebra stats maxaf --prevalence 0.0000625 --allelic 0.05 --penetrance 0.9 --faf95 0.00012
zebra edit 1-12345678-T-C --assembly GRCh38
```

## Layout

```
.claude-plugin/plugin.json      manifest (+ marketplace.json)
hooks/register.tsx              the mod: tools, doctrine, case board, /zebra, privacy gate
hooks/{tools,doctrine,privacy}.ts
types/index.d.ts                the mod's state contract
skills/<name>/SKILL.md          11 skills (vercel-labs/skills layout)
agents/*.md                     6 subagents
zebra/                          Python engine (stdlib only): sources/, commands/, acmg, stats, hpo_local, vcf, s2f, editing, case
bin/zebra                       CLI launcher
tests/                          pytest (offline + live-marked) and mod tests
docs/                           PRD, architecture, data sources and terms
```

## Data sources and terms

zebra-mod queries public resources live and keeps their provenance; each keeps its own licence (see [docs/DATA-SOURCES.md](docs/DATA-SOURCES.md)). Notably: AlphaGenome outputs are non-commercial and not for clinical decision-making; SpliceAI weights are CC BY-NC; OMIM content is not redistributed.

## Not a medical device

zebra-mod produces research-grade analyses to help people ask better questions. It does not diagnose, prescribe, or replace a clinician or genetic counsellor.

## Lineage

Built on lessons from [g2t-harness](https://github.com/zwbao/g2t-harness) (concentric harness: the model never writes a rank; clarify, never invent; falsification per skill) and [s2f-penguin](https://github.com/zwbao/s2f-penguin) (sequence-to-function models with receipts, claim ceilings, triangulation by logic not arithmetic).

MIT licence © 2026 Zhiwei Bao.
