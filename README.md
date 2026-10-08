# 🦓 zebra-mod

**A Claude Code mod that turns Claude Code into a rare-disease research workstation** — for patients and families on a diagnostic odyssey, for clinicians, and for researchers.

> *"When you hear hoofbeats, think horses, not zebras."* Medicine is taught to expect the common. People with rare diseases call themselves zebras because their answers lie off the beaten path. zebra-mod is built to look there — and to show its evidence for every step.

[中文说明 → README.zh.md](README.zh.md)

---

## What it does

| Job | How zebra-mod does it |
|---|---|
| **"What could this be?"** — phenotype-driven differential diagnosis | Records → HPO profile (verified terms, present/excluded, onset, source page) → three independent rankings (offline Resnik over HPO annotations, Monarch semantic similarity, PubCaseFinder) kept separate, agreement as signal → discriminating features → which test finds each candidate (exome, genome, CMA, repeat assays, methylation, mtDNA, metabolic) |
| **"What does this variant mean?"** — ACMG/AMP interpretation | VEP (MANE, AlphaMissense, REVEL, CADD) + SpliceAI at ±4,999 nt + gnomAD (grpmax, faf95; mtDNA heteroplasmy counts) + ClinVar (n of 4 stars) + ClinGen validity and dosage + MaveDB scores under their own calibration + LitVar (other alleles at the position excluded) + Chinese-cohort frequencies → data-driven codes with ClinGen-calibrated thresholds, PVS1 and PS1/PM5 **inputs** → your judged codes → **points computed by code** (Tavtigian 2020) and the 2015 combining rules side by side; GRCh37 input mapped to GRCh38 |
| **Sequence-to-function (S2F)** — splicing, non-coding, regulatory | SpliceAI/Pangolin (Broad lookup), AlphaGenome, Evo 2 and GPN-MSA through the [s2f-penguin](https://github.com/zwbao/s2f-penguin) `s2f` CLI; axis by axis with claim ceilings, combined by agreement, never averaged; mechanism → how to confirm with RNA, in a tissue GTEx shows the gene is expressed in; a splice-switching antisense screen (`aso_screen`: the aberrant event checked in the patient sequence, candidate target windows, published N-of-1 precedents) |
| **The other report forms** — microarray/CNV, exon-level deletions, SMN1 copy number, repeat expansions | `cnv_interpret`: genes covered, ClinGen dosage sensitivity, ACMG/ClinGen CNV evidence inputs (never a classification), and for an out-of-frame exon deletion the frame arithmetic that says which exon skipping would restore it |
| **Exome/genome reanalysis** on your own machine | `zebra qc`: sex check, KING kinship (swaps, non-paternity), runs of homozygosity, Mendelian errors and uniparental disomy, mosaic de novo calls. VCF triage: quality, inheritance models (de novo, homozygous, compound het, X-hemizygous), PED input, phenotype-gene restriction, SpliceAI on the top splice/non-coding candidates; a whole exome through a MyVariant frequency prefilter. The VCF file stays on the computer; variant positions go to the annotation services — by default only the filtered candidates, with the prefilter every variant (the mod asks first) |
| **Genotype → therapy (G2T)** | Mechanism first (LoF / GoF / DN / splicing / repeat) → approved and investigational drugs (Open Targets / ChEMBL), **access**: FDA and EMA status from the agencies' own records, approval in China from official NMPA/CDE documents, China's 2025 reimbursement list with its restriction text, trials with sites in China; literature, N-of-1 screens (antisense, gene replacement, **base editing feasibility**) → leads tiered A–E with a mechanism-fit check |
| **Classical statistics** | Cosegregation LR → PP1, maximum credible allele frequency (BS1), carrier frequency and genetic prevalence, recurrence risk (incl. Bayesian X-linked), Fisher/burden, de novo enrichment, Kaplan–Meier natural history + log-rank, N-of-1 trial design and analysis |
| **Literature** | Europe PMC, PubTator3, LitVar2 — PMIDs only as returned, findings quoted from the paper |
| **For families** | Plain-language explanations (Chinese by default for Chinese speakers), VUS explained honestly, visit preparation, recurrence and cascade testing, patient organisations, China's national rare disease lists and the collaboration-network hospitals by province |
| **Over months and years** | `case_recheck` asks the case's questions again — ClinVar class and stars of each variant, ClinGen validity of its genes, recruiting trials and new papers for each open hypothesis — and reports what changed since the last check, into the case timeline |
| **Reports** | Clinician summary, family letter and visit-preparation sheet, every claim cited, an **independent adversarial audit** before it is final, then exported to **Word and PDF** with Chinese typography — so a clinician or volunteer can run zebra-mod for a family and hand them the files |

## Principles (enforced, not just written)

1. **Evidence, not recall.** Every answer cites a ledger id, PMID or database record retrieved in the session. Each case keeps an append-only **evidence ledger** (`E1, E2, …`) of every source used.
2. **Code scores, the model interprets.** Rankings, ACMG points and statistics are computed by the `zebra` engine; the model never writes or adjusts a number.
3. **Clarify, never invent.** Assembly, transcript, zygosity, inheritance, sex — missing and material means ask, or conclude conditionally.
4. **Research-grade, not clinical.** No diagnosis delivered as fact, no dosing; next steps are questions for the care team.
5. **Privacy, stated exactly.** Case files are written only on your machine, and zebra's own database queries carry biology (HPO ids, genes, variants, disease names), never a name, a date of birth or a record number: a **privacy gate** in the mod refuses an outgoing call that carries the case's registered identifiers or an ID-number/phone/email pattern, asks before a file from the case folder or a raw genome file (VCF/BAM/CRAM/FASTQ) leaves the machine, and closes (refusing outgoing calls) if it cannot read the case's identifier list. **What the gate cannot do:** every record you ask Claude to read — every PDF, photo and report — is sent to the model provider as part of the conversation, like any other file you open in Claude Code; the gate inspects tool calls, not the conversation. Register identifiers early (the model does it with `case_update` → `identifiers` when records name the patient), and redact before sharing if that matters to you. Exports to Word/PDF refuse a report that still holds a registered identifier or an ID number. Matching is best effort over encodings and spellings, not a guarantee.

## How good is it

Measured, not claimed — [docs/BENCHMARK.md](docs/BENCHMARK.md) has the method, the data digests and the confidence intervals.

- **Phenotype ranking** on GA4GH phenopacket-store 0.1.27 (10,374 published cases, 780 diseases), split by disease into development and held-out halves. On held-out cases whose own paper is **not** one of the HPO annotation sources — the fair test for a new patient — the correct disease is in the local top 10 for **28%** (0.1.0: 15%) and the causal gene for 40%. Across all held-out cases it is 70%, an upper bound: for most published cases the disease's annotations were curated from that very paper. On a 100-case web sample, local + Monarch + PubCaseFinder together put the correct disease in the top 10 for 75%. A ranking is a list of hypotheses to test.
- One honest cost: the PRD's own example (febrile, focal and tonic-clonic seizures, developmental delay; hypotonia excluded) now ranks Dravet syndrome 16th locally (6th in 0.1.0), because excluded terms are flagged rather than scored — the choice that measured better on the held-out half.
- `evals/` holds 10 end-to-end cases for `claude plugin eval` (a Chinese parent's records, urgent symptoms first, a VUS, a CNV, mtDNA, therapy and trials, privacy, an out-of-scope request).

## What the mod adds to Claude Code

- **21 tools** the model calls directly (`mcp__zebra-mod__*`): `case_status`, `case_update`, `case_recheck`, `hpo_search`, `phenotype_rank`, `gene_card`, `variant_card`, `disease_card`, `acmg`, `cnv_interpret`, `s2f_predict`, `therapy_landscape`, `trials_search`, `literature_search`, `rare_stats`, `edit_check`, `china_rare`, `access`, `expression`, `aso_screen`, `report_export`. Every call goes through your permission rules and the privacy gate; read-only lookups skip the prompt unless a rule of yours says otherwise, and case writes follow your permission mode (with an "allow for this session" choice).
- **A research doctrine** in the system prompt (the rules above), with the active case.
- **You can see it at work** (all additive; unrelated work draws as before): its tool calls get their own rows (`🦓 Variant card  NM_001165963.4:c.2134C>T`, then `⎿ 15 evidence rows (E39–E53) · from Ensembl VEP, gnomAD, ClinVar, LitVar2 · 3 cached`) instead of a JSON envelope, and a case_update row never shows the identifiers it registers; the spinner names the databases being queried; a small galloping zebra above the prompt while a query runs, a shield there after the privacy gate stops a call; a `🦓` label in the footer (with the open case's title); and one line at the end of a turn that used it — calls, databases, evidence rows, calls stopped. Chinese or English, from the case's language or your prompts. Option `interface`: `full` (default), `quiet` (no animation, no spinner text) or `off`.
- **A case board** pane (`/zebra board`) and a status line: phenotypes, variants and their research class, hypotheses, therapy leads, open questions, evidence count — live as the case changes.
- **`/zebra`** command: `new <dir> [title]`, `case <dir>`, `board`, `ledger`, `doctor`, `close`.
- **12 skills** (`/zebra-mod:zebra-start` routes): `zebra-safety` (urgent red flags and disease-specific drug, anaesthesia and procedure hazards), `zebra-intake`, `zebra-diagnose`, `zebra-variant`, `zebra-reanalysis`, `zebra-s2f`, `zebra-therapy`, `zebra-stats`, `zebra-literature`, `zebra-family`, `zebra-report`.
- **6 subagents**: `phenotype-curator`, `variant-curator`, `s2f-analyst`, `therapy-scout`, `literature-scout`, `evidence-auditor`.
- **The `zebra` CLI** (Python standard library only, Python ≥ 3.9): the single implementation behind every tool, usable from any shell, notebook or other agent.

## Install

Requirements: Claude Code ≥ 2.1.289 (function-hook mods) and Python ≥ 3.9. No sudo, no pip installs.

**One sentence.** In Claude Code, say:

```text
Install https://github.com/zwbao/zebra-mod
```

Claude Code reads [INSTALL.md](INSTALL.md) and does the rest: checks versions and Python, registers and installs the plugin, downloads the HPO release (~80 MB), and runs `zebra doctor`. Then start a new session (or type `/reload-plugins`) and run `/zebra-mod:zebra-start`.

**One command** (the same steps):

```bash
curl -fsSL https://raw.githubusercontent.com/zwbao/zebra-mod/main/install.sh | sh
# while the repository is private, clone it first:
git clone https://github.com/zwbao/zebra-mod ~/zebra-mod && ~/zebra-mod/install.sh
```

**By hand:**

```bash
claude plugin marketplace add zwbao/zebra-mod
claude plugin install zebra-mod@zebra-mod
python3 <installPath>/bin/zebra hpo fetch     # optional: offline phenotype ranking, instant HPO search
# one session only, no install: claude --plugin-dir ~/zebra-mod
```

Optional heavy S2F models: install the `s2f` CLI from [s2f-penguin](https://github.com/zwbao/s2f-penguin) (`uv tool install "git+https://github.com/zwbao/s2f-penguin"`), and set `ALPHAGENOME_API_KEY` (AlphaGenome: non-commercial, not for clinical decisions) and/or `NVCF_RUN_KEY` (Evo 2 on NVIDIA). Check everything with `/zebra doctor`.

Options (`/config` → zebra-mod, or `pluginConfigs` in settings): `python` (interpreter), `doctrine` (`auto` by default: the full rules while a case is open, otherwise a short section that only applies to rare-disease questions, so the rest of your work in Claude Code is untouched; or `always` | `case` | `off`), `privacyGate` (on by default; with no case open, an email or phone number in another tool's call is asked about rather than refused), `interface` (`full` by default: its own tool rows, the spinner text, the galloping zebra and the end-of-turn line; `quiet` without the animation and spinner text; `off` for Claude Code's own drawing only).

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
skills/<name>/SKILL.md          12 skills (vercel-labs/skills layout)
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

zebra-mod's sequence-to-function (S2F) layer builds on [s2f-penguin](https://github.com/zwbao/s2f-penguin): sequence-to-function models with run receipts, claim ceilings, and triangulation by logic rather than arithmetic. zebra-mod calls its `s2f` CLI for the heavy models.

MIT licence © 2026 zwbao.
