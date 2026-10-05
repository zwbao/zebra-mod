# Changelog

## 0.1.0 — 2026-10-05

First release.

- **Mod** (`hooks/register.tsx`): 15 model-callable tools backed by the `zebra` CLI; research doctrine in the system prompt (`always` / `case` / `off`); case board pane and status line; `/zebra` command (`new`, `case`, `board`, `ledger`, `doctor`, `close`); auto-approval of the mod's own read-only tools; privacy gate (denies outgoing calls carrying a case's protected identifiers or ID-number/phone/email patterns, including Bash calls to the zebra CLI; asks before raw genome files leave the machine).
- **Engine** (`zebra/`, Python standard library, ≥ 3.9): HTTP layer with caching, retries and provenance; case workspace with an append-only evidence ledger; ACMG/AMP points (Tavtigian 2020) beside the 2015 combining rules, ClinGen SVI warnings and calibrated suggestions (REVEL per Pejaver 2022, SpliceAI per Walker 2023, BA1/BS1/PM2_Supporting); classical statistics; offline HPO ranking (Resnik best-match average) with English and official Chinese label search; Ensembl VEP client with AlphaMissense/REVEL/CADD/SpliceAI; sources for HPO, Monarch, PubCaseFinder, Orphanet, OLS, GeneReviews, gnomAD, ClinVar, ClinGen, PanelApp, UniProt, Europe PMC, PubTator, LitVar, ClinicalTrials.gov, Open Targets; S2F bridge (SpliceAI/Pangolin lookup; AlphaGenome, Evo 2, GPN-MSA via s2f-penguin); base-editing feasibility screen; VCF inspection and trio triage; `zebra doctor`.
- **Skills** (11) and **subagents** (6).
