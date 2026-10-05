# zebra-mod — product requirements and design

Version 0.1 · 2026-10-05

## 1. Problem

Hundreds of millions of people live with one of more than 6,000 rare diseases (263–446 million by the Orphanet-based estimate of Nguengang Wakap et al., Eur J Hum Genet 2020). Diagnosis often takes years; exome sequencing leaves most tested patients without a molecular diagnosis; variants of uncertain significance are common answers; and most rare diseases have no approved treatment. The knowledge exists but is scattered across ontologies (HPO, Orphanet, Mondo, OMIM), variant resources (ClinVar, gnomAD, ClinGen), prediction models (SpliceAI, AlphaMissense, AlphaGenome, Evo 2), literature, trial registries and drug databases, each with its own identifiers and caveats.

General-purpose AI assistants answer rare-disease questions from memory — fluent, confident, and unverifiable. For this domain that is the wrong tool: a wrong HPO id, an outdated ClinVar classification or a mis-tiered therapy lead does real harm.

## 2. Users and jobs

| User | Jobs |
|---|---|
| Patient / family, undiagnosed | organise records; understand what the phenotype could mean; know which test to ask for; prepare visits |
| Patient / family, with a report | understand a variant or VUS; disease natural history; treatments, trials, N-of-1 routes; recurrence; support |
| Clinician | differential diagnosis; variant interpretation with evidence; reanalysis; next tests; concise summaries |
| Researcher | exome/genome reanalysis; ACMG with calibrated predictors; S2F mechanism analysis; cohort statistics; therapy hypotheses; novel gene candidates |

## 3. Product principles

1. **Evidence, not recall** — every claim traces to a source retrieved in the session; each case keeps an append-only evidence ledger.
2. **Code scores, the model interprets** — rankings, ACMG points, statistics come from deterministic code (g2t-harness: "the LLM never writes a rank").
3. **Clarify, never invent** — missing assembly/transcript/zygosity/inheritance → ask or condition.
4. **Judgement in prose, verification in code** — the skills hold the judgement criteria (ACMG codes, mechanism fit, test choice); code holds arithmetic, structure and provenance.
5. **Axes combine by logic, not arithmetic** — S2F models and phenotype rankers are reported separately; agreement is the signal (s2f-penguin).
6. **Research-grade, honest** — no diagnosis as fact, no dosing; every lead tiered.
7. **Privacy by construction** — local cases, a privacy gate on outgoing calls, genome files never uploaded silently.
8. **Zero-install core** — Python standard library only; heavy models optional.

## 4. Architecture

```
Claude Code
 ├─ zebra-mod (function-hook mod: hooks/register.tsx)
 │   ├─ prompt.compose  → research doctrine + active case
 │   ├─ $.tool.register → 15 tools  ──┐
 │   ├─ tool.check      → auto-allow own read-only tools; privacy gate (deny identifiers/PII; ask before genome uploads)
 │   ├─ /zebra command, case board pane, status line ($.state, $.store)
 │   └─ PATH += bin/, ZEBRA_CASE env for Bash
 ├─ skills/ (11)  agents/ (6)            │
 └─────────────────────────────────────── ▼
   zebra CLI (Python stdlib) — the single implementation
    ├─ sources/: hpo, monarch, pubcasefinder, orphanet, ols, genereviews, ensembl (VEP+plugins),
    │            gnomad, clinvar, clingen, panelapp, uniprot, gene, variant, europepmc, pubtator,
    │            litvar, ctgov, opentargets
    ├─ hpo_local (offline Resnik BMA ranking from HPO release files)
    ├─ acmg (Tavtigian points + Richards rules + ClinGen-calibrated suggestions)
    ├─ stats (segregation, maxaf, carrier, recurrence, Fisher/burden, de novo, KM/log-rank, N-of-1)
    ├─ vcf (inspect, trio triage), s2f (SpliceAI/Pangolin lookup; s2f-penguin bridge), editing (ABE/CBE screen)
    └─ case (case.json, evidence/ledger.jsonl, records/, reports/)
```

Every command returns one envelope: `{ok, command, query, result, sources[], warnings[], ledger[]}`. A failing source becomes a named warning, never a silent gap.

## 5. Components and acceptance

| Component | Acceptance |
|---|---|
| Case workspace + ledger | init/apply/summary; HPO ids verified on write; sequential ledger ids under a file lock; identifiers never echoed |
| Phenotype ranking | Dravet-like profile → Dravet in the top 5 of ≥ 2 independent methods |
| Variant card | CFTR F508del, SCN1A R712*, a missense with REVEL/AlphaMissense, a GRCh37 input — consistent normalisation, frequencies, ClinVar stars |
| ACMG | points and 2015 rules agree on textbook combinations; SVI warnings (PM2_Supporting, PVS1+PP3, PM1+PP3 cap, PP5/BP6) |
| S2F | SpliceAI/Pangolin deltas for a splice variant; heavy models ran or `not_run` with a reason |
| VCF triage | synthetic trio: the SCN1A de novo ranked first; VCF never sent |
| Therapy / trials / literature | Dravet, SMA, CFTR return known drugs and trials with sources |
| Mod | loads in Claude Code; tools callable; doctrine present; privacy gate denies identifiers even with permissions bypassed |
| End to end | ≥ 2 cases, ≥ 2 disease classes, run by a fresh Claude Code session from a user-style prompt; audited by an independent reviewer |

## 6. Out of scope (v0.1)

Clinical reporting; CNV calling from BAM; long-read/repeat calling; local deep-learning model execution (delegated to s2f-penguin); writing to external matchmaking services (summaries are prepared for humans to submit).

## 7. Roadmap

Kept current in [ROADMAP.md](ROADMAP.md), ordered P0/P1/P2.
