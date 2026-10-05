---
name: zebra-intake
description: Turn a patient's records (PDF reports, photos of lab sheets, clinic letters, genetic test reports, the family's own notes) into a structured zebra case — HPO phenotypes present and excluded with onset and source page, variants exactly as reported, family history, tests already done, and protected identifiers. Use before diagnosis or variant work, or when new records arrive. Triggers: 整理病历, 整理资料, organize records, intake, phenotype profile, HPO coding, 表型.
---

# zebra-intake — records → case

Output: an updated case (via `mcp__zebra-mod__case_update`) and a short summary of what is recorded and what is missing. Never interpret beyond what the records say.

## Steps

1. **Case.** `mcp__zebra-mod__case_status`. No case → ask the person to run `/zebra new <folder>` and put files in `<folder>/records/`.
2. **Profile.** Set the case `profile` (via `case_update`): `role` and `language` from how the person writes (a mother writing Chinese → family, zh; a fellow writing English → clinician or researcher, en), proband sex and age band. No names or birth dates.
3. **Inventory.** List `records/`. More than 6 files → fan out: one `zebra-mod:phenotype-curator` subagent per ~5 files, all in one message, each returning its JSON; otherwise read them yourself (Read handles PDF and images).
4. **Protect identifiers first.** Patient and relatives' names, birth dates, ID/record/insurance numbers, phone, address → `zebra case identifiers --add "<value>" ...` via Bash with `--case <case dir>`. Never repeat them in chat or in any other tool call.
5. **Phenotypes.** For each clinical finding:
   - phrase it in English clinical terms, `mcp__zebra-mod__hpo_search`, choose the most specific term the record supports (not more specific);
   - present vs excluded: excluded only when the record states absence or a normal result for that feature (normal brain MRI → excluded `Abnormality of brain morphology`);
   - onset when stated (HPO onset terms or an age), `source` = file and page;
   - never code the diagnosis itself, a drug, or a test name as a phenotype; a lab value becomes a phenotype only with its direction (elevated CK → `Elevated circulating creatine kinase concentration`).
6. **Variants.** Copy exactly as printed: gene, transcript HGVS (`NM_...:c.`), protein change, genomic coordinates **with the assembly the report names** (Chinese reports often use hg19/GRCh37), zygosity, inheritance if parents were tested, the lab's classification. Do not normalise here; `zebra-variant` does that.
7. **Family and tests.** Consanguinity, affected relatives (who, what, age), parental testing; tests done and results (CMA, panel, exome/genome — singleton or trio, year, lab), biochemical/metabolic tests, imaging, EEG/EMG, biopsy. Put these in the case's notes via questions or variants' `source` fields; record "test not done" gaps as questions.
8. **Write once.** One `case_update` call with all phenotypes, variants and questions.
9. **Report back.** Board summary, then the gaps that would change the analysis (e.g. no parental samples, exome from 2019 never reanalysed, no metabolic screen, onset ages unknown) as questions for the care team.

## Quality bar

- Every phenotype has a source; every HPO id came from `hpo_search`.
- 5–15 specific phenotypes beat 40 vague ones; include the striking, unusual ones (they discriminate).
- When a record is ambiguous, record the finding with a note, or ask — do not guess.

What would make this intake wrong: a phenotype that no record states; an excluded term inferred from silence; a variant transcribed with a different transcript or assembly than printed.
