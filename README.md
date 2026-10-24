<div align="center">

# 🦓 zebra-mod

**A rare-disease research workstation inside Claude Code**

From symptoms to a differential, from a variant to an ACMG class, from a sequence to its function, from a genotype to a therapy — every claim tied to a source you can open.

[![License: MIT](https://img.shields.io/badge/license-MIT-0f5c6e.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.3.9-0f5c6e.svg)](CHANGELOG.md)
[![Claude Code](https://img.shields.io/badge/Claude%20Code-%E2%89%A5%202.1.289-1b2430.svg)](https://claude.com/claude-code)
[![Python](https://img.shields.io/badge/python-%E2%89%A5%203.9%20·%20stdlib%20only-1b2430.svg)](zebra/)

**English** · [中文](README.zh.md)

</div>

> *"When you hear hoofbeats, think horses, not zebras."* Medicine is taught to expect the common. People with rare diseases call themselves zebras because their answers lie off the beaten path. zebra-mod is built to look there, and to show its evidence at every step.

zebra-mod is a [Claude Code](https://claude.com/claude-code) mod for **families on a diagnostic odyssey, clinicians, and researchers**. It gives Claude 22 tools over live public databases (HPO, Monarch, Orphanet, ClinVar, gnomAD, Ensembl VEP, ClinGen, PanelApp, Open Targets, ClinicalTrials.gov, Europe PMC and more), a case workspace with an evidence ledger, and a set of procedures for diagnosis, variant interpretation, therapy search and reporting. You talk to Claude as usual, in English or Chinese; it reaches for the tools itself.

## Highlights

- **Differential diagnosis from records.** Clinic notes, lab sheets and photos become a verified HPO profile; three independent rankers (offline, Monarch, PubCaseFinder) are reported side by side, next to Claude's own clinical reasoning, never instead of it.
- **Variant interpretation with code-computed ACMG points.** VEP, gnomAD, ClinVar, ClinGen, SpliceAI, AlphaMissense, REVEL and MaveDB feed the evidence; you justify the codes, zebra does the arithmetic (Tavtigian 2020 points beside the 2015 rules).
- **Beyond SNVs.** CNVs and microarray results against ClinGen dosage maps, exon deletions with reading-frame arithmetic, SMN1/SMN2 copy number, repeat expansions (D4Z4, DMPK, FMR1, HTT…) sized by rule.
- **Sequence-to-function, run in the mod.** SpliceAI, Pangolin, GPN-MSA, Evo 2 and AlphaGenome, each reported on its own axis with the highest claim it supports, combined by agreement and never averaged.
- **Genotype to therapy, with access.** Approved and investigational drugs by mechanism, FDA/EMA status from the agencies' records, approval and reimbursement in China, recruiting trials with sites in China, N-of-1 routes (antisense, base editing).
- **Reports a family can hold.** Clinician summaries and plain-language family sheets, audited by an independent subagent, exported to HTML, Word and PDF with Chinese typography.
- **Private by construction.** Case files stay on your machine; a privacy gate keeps registered names, birth dates and record numbers out of every outgoing call.

## Quick start

**1. Install.** In Claude Code, say:

```text
Install https://github.com/zwbao/zebra-mod
```

Claude reads [INSTALL.md](INSTALL.md) and does the rest: checks Claude Code (≥ 2.1.289) and Python (≥ 3.9), installs the plugin, fetches the HPO release (~80 MB) and runs a health check. Or from a shell:

```bash
curl -fsSL https://raw.githubusercontent.com/zwbao/zebra-mod/main/install.sh | sh
```

**2. Try the demo.** Start a new session and type `/zebra demo`. It opens a synthetic case (a clinic note and a genetic report) with the first question ready; press Enter and watch it read the records, build the differential and interpret the variant.

**3. Ask about your own case.** No commands to learn:

```text
My daughter started having fever-triggered seizures at 6 months. The genetic report says
SCN1A NM_001165963.4:c.2134C>T. What does it mean, and what should we ask the neurologist?
```

`/zebra` shows what the mod can do; `/zebra new <folder>` starts a case for a patient whose records you want to keep together.

## What you can ask

| You are | For example |
|---|---|
| **A patient or family** | "跑了好几家医院都没查出来，这是孩子的检查单，可能是什么病？" · "What does a VUS mean for my son?" · "Are there trials for FSHD in China?" · "Write a one-page list of questions for our next visit." |
| **A clinician** | "Ranked differential for these HPO terms, and which test finds each candidate." · "ACMG reading of this variant with the segregation in this pedigree." · "Management and surveillance plan for a patient with DM1 and conduction disease." · "Is a catheter angiogram safe in vascular EDS?" |
| **A researcher** | "Re-analyse this trio VCF with the phenotype." · "SpliceAI, Pangolin and AlphaGenome for this deep-intronic variant." · "Cosegregation likelihood ratio for these meioses." · "Which exon skip restores the frame for a DMD exon 45–50 deletion?" |

A full list of what each capability does, tool by tool, is in [docs/FEATURES.md](docs/FEATURES.md).

## How it works

```
you
└─► Claude Code
      ├─ research doctrine (system prompt) · 13 skills · 6 subagents
      └─► zebra-mod hooks
            ├─ privacy gate
            └─► 22 tools (mcp__zebra-mod__*)
                  └─► zebra engine (Python, stdlib only)
                        ├─► public databases and models: HPO · Monarch · ClinVar · gnomAD · VEP · ClinGen …
                        └─► case folder on your machine: case.json · evidence/ledger.jsonl (E1, E2, …) · reports/
```

- **The engine computes, the model interprets.** Rankings, ACMG points, repeat sizes and statistics come from the `zebra` engine; the model never writes or adjusts a number, and openly overrides a ranking when evidence the ranker cannot see decides the question.
- **Every source is kept.** Each tool result is appended to the case's evidence ledger with its database, record and retrieval time; reports cite those sources.
- **Claude's own knowledge is welcome, labelled.** What was retrieved this session and what comes from recall are marked apart, so a reader can tell which is which.
- **The mod adds, it never subtracts.** Nothing in it is a gate the person must pass before getting an answer: the answer comes first, the case work around it.

The same engine is a command-line tool you can use from any shell, notebook or agent:

```bash
zebra phenotype rank --present HP:0002373 HP:0001263 --sources local,monarch,pubcasefinder
zebra acmg suggest NM_001165963.4:c.2134C>T --inheritance AD
zebra acmg classify PVS1 PS2 PM2_Supporting
zebra s2f predict 12-25245350-C-T --models spliceai,pangolin,gpn_msa
zebra vcf triage trio.vcf.gz --proband P --mother M --father F --sex female --hpo-genes
zebra stats maxaf --prevalence 0.0000625 --allelic 0.05 --penetrance 0.9
```

## How well it works

Measured, with the method and confidence intervals in [docs/BENCHMARK.md](docs/BENCHMARK.md):

- **Phenotype ranking** (GA4GH phenopacket-store, 10,374 published cases, split by disease). On held-out cases whose own paper was not used to build the HPO annotations — the fair test for a new patient — the correct disease is in the offline top 10 for **28%** and the causal gene for **40%**. A ranking is a list of hypotheses, which is why Claude's own differential stands beside it.
- **With vs without the mod** (six held-out end-to-end cases written without access to the mod's instructions, 8 runs per arm, graded item by item against guideline ground truth). Pooled, the two arms are level (difference −0.02, 95% CI −0.16 to +0.12): the mod was better on an emergency buried in long records and on a VUS with segregation, worse on two cases. "Never worse" is not yet established, and the benchmark says so.
- **19 end-to-end evals** in [`evals/`](evals/) run with `claude plugin eval`, from a Chinese parent's records to a somatic tumour report that must be read on oncology's terms, not as a rare disease.

## Privacy and safety

- **Stays local:** case files, records and the evidence ledger are written only on your machine.
- **Never sent by zebra:** names, birth dates and record numbers. Database queries carry biology (HPO ids, genes, variants); a privacy gate refuses an outgoing call that contains a registered identifier or an ID-number, phone or email pattern, and asks before a case file or raw genome file leaves the machine.
- **What the gate cannot do:** any record you ask Claude to read is sent to the model provider as part of the conversation, as with any file in Claude Code. Redact first if that matters to you.
- **Urgent first:** a symptom that needs care today is said first, before any analysis.
- **Research-grade:** zebra-mod drafts analyses for a clinician's review. It is not a medical device, does not diagnose or prescribe on its own authority, and does not replace a clinician or genetic counsellor.

## Configuration

Options live under `/config` → zebra-mod (or `pluginConfigs` in settings):

| Option | Default | What it does |
|---|---|---|
| `doctrine` | `auto` | Full research rules while a case is open, a short rare-disease-only section otherwise; or `always`, `case`, `off` |
| `privacyGate` | on | The outgoing-call check described above |
| `interface` | `full` | Tool rows, spinner text and the galloping zebra; `quiet` without animation; `off` for Claude Code's own drawing |
| `python` | `python3` | Interpreter for the engine |

Optional sequence models (`/zebra doctor` checks them all):

| Model | Needs |
|---|---|
| SpliceAI, Pangolin, GPN-MSA | nothing |
| Evo 2 | `NVCF_RUN_KEY` (NVIDIA hosted) |
| AlphaGenome | `ALPHAGENOME_API_KEY` and `pip install alphagenome` (or `ZEBRA_ALPHAGENOME_PYTHON`); non-commercial, not for clinical decisions |

## Documentation

| | |
|---|---|
| [docs/FEATURES.md](docs/FEATURES.md) | Every capability, tool, skill and subagent |
| [docs/BENCHMARK.md](docs/BENCHMARK.md) | Ranking benchmark and with/without-mod evaluations |
| [docs/DATA-SOURCES.md](docs/DATA-SOURCES.md) | Each database and model, its terms and licence |
| [docs/PRD.md](docs/PRD.md) · [docs/ROADMAP.md](docs/ROADMAP.md) | Design and what comes next |
| [INSTALL.md](INSTALL.md) · [CHANGELOG.md](CHANGELOG.md) | Installation steps and release notes |

## Development

```
hooks/        the mod: tools, doctrine, privacy gate, case board, /zebra   (TypeScript)
skills/       13 procedures, routed by zebra-start
agents/       6 subagents (curators, scouts, evidence auditor)
zebra/        the engine: sources, ACMG, statistics, VCF, S2F, reports    (Python, stdlib only)
evals/        end-to-end cases for `claude plugin eval`
tests/        pytest (offline; live network with ZEBRA_NETWORK_TESTS=1)
```

```bash
python3 -m pytest -q            # engine
claude plugin test .            # mod hooks
claude plugin validate .        # manifest, skills, agents
claude --plugin-dir .           # run this checkout for one session
```

Issues and pull requests are welcome. A change to the doctrine or a skill should come with a with/without-mod eval run (`claude plugin eval --ablation with-without`); held-out rubrics in `evals/heldout-*` are never edited to fit a result.

## Data sources and licences

zebra-mod queries public resources live and records where each answer came from; every resource keeps its own terms ([docs/DATA-SOURCES.md](docs/DATA-SOURCES.md)). Notably: AlphaGenome outputs are non-commercial and not for clinical decision-making, SpliceAI weights are CC BY-NC, and OMIM content is not redistributed.

## Citation

```bibtex
@software{zebra_mod,
  author  = {zwbao},
  title   = {zebra-mod: a rare-disease research workstation for Claude Code},
  year    = {2026},
  version = {0.3.9},
  url     = {https://github.com/zwbao/zebra-mod}
}
```

## Acknowledgements

Built on the work of the HPO, Monarch, Orphanet, ClinGen, ClinVar, gnomAD, Ensembl, PanelApp, Open Targets, Europe PMC, PubTator, LitVar and PubCaseFinder teams, and of the authors of SpliceAI, Pangolin, GPN-MSA, Evo 2 and AlphaGenome. The galloping zebra above the prompt is traced from Eadweard Muybridge's *The Horse in Motion* (1878, public domain).

MIT licence © 2026 zwbao.
