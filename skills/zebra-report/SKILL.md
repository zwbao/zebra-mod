---
name: zebra-report
description: Write the case up — a clinician summary (phenotype profile, variants with ACMG codes and points, differential with method scores, sequence-to-function evidence, tiered therapy leads, open questions, methods and sources) and/or a family letter in plain language — every claim tied to a ledger id or PMID, then audited by an independent evidence-auditor subagent before it is final. Triggers: 写报告, write a report, summary for the doctor, 给医生的总结, family letter, 给家属的信, case summary, 病例总结.
---

# zebra-report — write, audit, fix, deliver

## Steps

1. **Gather**: `mcp__zebra-mod__case_status`; `zebra case ledger --case <case>` (Bash) for the evidence rows; re-run a tool if a needed fact is missing — never fill from memory.
2. **Draft** to `<case>/reports/<kind>-<YYYY-MM-DD>.md`:

   **Clinician summary** (language per case; terminology exact)
   1. Patient (no identifiers): sex, age band, ancestry if relevant, consanguinity, key history.
   2. Phenotype profile: HPO table (id, label, onset, source).
   3. Genetic findings: variant table (gene, transcript HGVS, protein, zygosity, inheritance, lab class, zebra research class with points and codes).
   4. Differential: candidates with per-method rank/score, fit and misfit features, decisive tests.
   5. Sequence-to-function evidence: axis matrix with models, scores, ceilings.
   6. Therapy and trials: tiered leads (A–E) with mechanism fit.
   7. Open questions and recommended next steps (as considerations for the care team).
   8. Methods and sources: tools, database versions/dates, ledger ids; "research-grade analysis, not a clinical report".

   **Family letter**: what was looked at, what was found, what it means and does not mean, what to ask, where to find support. Plain language; no ids except where they help (disease name, gene).

   Every factual sentence ends with its evidence (`[E12]`, `[PMID 12345678]`).
3. **Audit**: dispatch `zebra-mod:evidence-auditor` with the draft path and the case path. It works independently and returns P0/P1/P2 findings.
4. **Fix** every P0 and P1 (re-query, correct, or remove the claim); note P2s you leave and why. Re-audit if P0s were found.
5. **Deliver**: final path(s), what changed after audit, what remains uncertain.

## Never

A number without a source; an identifier from the records; a classification presented as clinical; a lead without its tier.
