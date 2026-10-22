---
name: zebra-report
description: Write the case up and export it to Word/PDF — a clinician summary (phenotype profile, variants with ACMG codes and points, differential with method scores, sequence-to-function evidence, tiered therapy leads, open questions, methods and sources) and/or a family letter in plain language — every claim tied to a ledger id or PMID, then audited by an independent evidence-auditor subagent before it is final; also the family letter and visit-preparation sheet a clinician or volunteer hands to a family (代操作). Triggers: 导出, export, PDF, Word, 打印, 代操作, 写报告, write a report, summary for the doctor, 给医生的总结, family letter, 给家属的信, case summary, 病例总结.
---

# zebra-report — write, audit, fix, deliver

## Steps

1. **Gather**: `mcp__zebra-mod__case_status`; `zebra case ledger --case <case>` (Bash) for the evidence rows; re-run a tool if a needed fact is missing; anything still from memory is labelled as such, never passed off as retrieved.
2. **Draft** to `<case>/reports/<kind>-<YYYY-MM-DD>.md`:

   **Clinician summary** (language per case; terminology exact)
   1. Patient (no identifiers): sex, age band, ancestry if relevant, consanguinity, key history.
   2. Phenotype profile: HPO table (id, label, onset, source).
   3. Genetic findings: variant table (gene, transcript HGVS, protein, zygosity, inheritance, lab class, zebra research class with points and codes).
   4. Differential: candidates with per-method rank/score, fit and misfit features, decisive tests.
   5. Sequence-to-function evidence: axis matrix with models, scores, ceilings.
   6. Therapy and trials: tiered leads (A–E) with mechanism fit.
   7. Open questions and recommended next steps (for a clinician reader: the management plan and dose ranges as proposals for their review; for a family-facing document: considerations for the care team).
   8. Methods and sources: tools, database versions/dates, ledger ids; which statements are from clinical knowledge rather than retrieved; "research-grade analysis, not a clinical report" — except for findings an accredited laboratory has already reported, which are stated as the laboratory reported them.

   **Family letter**: what was looked at, what was found, what it means and does not mean, what to ask, where to find support. Plain language; no ids except where they help (disease name, gene).

   Every retrieved fact ends with its evidence (`[E12]`, `[PMID 12345678]`); clinical reasoning, guideline content and management proposals from your own knowledge are kept and labelled (凭记忆 / clinical judgment, to verify) — a report is never thinner than the answer given in conversation.
3. **Audit**: dispatch `zebra-mod:evidence-auditor` with the draft path and the case path. It works independently and returns P0/P1/P2 findings.
4. **Fix** every P0 and P1: correct a wrong claim; for an unsourced one, re-query, or keep it labelled as clinical knowledge (to verify) — remove a claim only when it is wrong or cannot be supported even as knowledge. Dose ranges and management proposals in a clinician document stay (marked as proposals). Note P2s you leave and why. Re-audit if P0s were found.
5. **Clean** before export: no zebra skill or tool names, `/zebra` commands or bare ledger ids (E12) in the text — replace each with the source it stands for (PMID, NCT, database record).
6. **Export** what goes to people outside the terminal: `mcp__zebra-mod__report_export` with the report path (Word and PDF by default, `formats: ["html"]` for a page to open in a browser). It refuses a report that still holds a protected identifier — remove it, never pass around the check. A missing PDF engine is reported; Word is still made.
7. **Deliver**: final path(s) of the .md and the exported files, what changed after audit, what remains uncertain.

## For a family (代操作)

When the person running zebra-mod does it on a family's behalf, the family receives files, never a terminal:

1. Work the case through the usual skills (intake → diagnose / variant / reanalysis → therapy → safety), with the operator.
2. Write, in plain Chinese unless the family reads another language:
   - `reports/family-letter-<date>.md` — 给家属的信: what was looked at, what was found, what it means and does not mean, in short paragraphs; one table at most; evidence ids kept but small (the export greys them).
   - `reports/visit-prep-<date>.md` — 就诊准备单 (`zebra-family` → visit preparation): the 3–5 questions to ask, tests to ask about and why, what to bring, and the hazards from `zebra-safety` for this disease.
3. Audit both (`evidence-auditor`), fix, then export both with `report_export` (Word and PDF).
4. Tell the operator which files to hand over, and what the family should hear in person rather than read: a possible diagnosis is discussed with the care team, never delivered by paper.

## Never

A number without a source; an identifier from the records; a classification presented as clinical; a lead without its tier; an exported handout that skipped the audit.
