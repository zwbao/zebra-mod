---
name: zebra-intake
description: Turn a patient's records (PDF reports, photos of lab sheets, clinic letters, genetic test reports, the family's own notes) into a structured zebra case — HPO phenotypes present and excluded with onset and source page, variants exactly as reported, family history, tests already done, and protected identifiers. Use before diagnosis or variant work, or when new records arrive. Triggers: 整理病历, 整理资料, organize records, intake, phenotype profile, HPO coding, 表型.
---

# zebra-intake — records → case

Output: an updated case (via `mcp__zebra-mod__case_update`) and a short summary of what is recorded and what is missing. Record what the records say; your interpretation belongs in the answer, labelled as yours.

**Proportion.** Intake serves the person's question; it never replaces it. If they asked something (a diagnosis, a plan, what a result means), answer it in the same turn — from the records you read — and let the structured case catch up. Keep the tool traffic small: register only the identifiers that point to a person, search HPO for the 5–15 findings that matter (send several `hpo_search` calls in one message), and write the case once.

## Steps

1. **Case.** `mcp__zebra-mod__case_status`. No case → ask the person to run `/zebra new <folder>` and put files in `<folder>/records/`.
2. **Profile.** Set the case `profile` (via `case_update`): `role` from what the person says about themselves — not from their language or writing style. Someone operating Claude Code on a case folder is a clinician/researcher unless they say they are the patient or a relative ("我妈妈的报告", "我的孩子"); when unclear and it matters, ask once. `language` from how they write; proband sex and age band. No names or birth dates.
3. **Inventory.** List `records/`. More than 6 files → fan out: one `zebra-mod:phenotype-curator` subagent per ~5 files, all in one message, each returning its JSON; otherwise read them yourself (Read handles PDF and images).
4. **Protect identifiers first** — before any lookup. What points to a person: patients' and relatives' names (汉字 and pinyin, each spelling), birth dates, national ID numbers, phone, address, insurance numbers, hospital record/admission/outpatient numbers → `case_update` with `identifiers: [...]` (the result shows only a count). Exam numbers, barcodes and report serials do not identify a person on their own and need not be registered. From the shell the same is `zebra case identifiers --add "<value>" ... --case <case dir>`. Never send them in any other tool call; in replies to a family do not repeat them.
5. **Phenotypes.** For each clinical finding:
   - phrase it in English clinical terms, `mcp__zebra-mod__hpo_search`, choose the most specific term the record supports (not more specific);
   - present vs excluded: excluded only when the record states absence or a normal result for that feature (normal brain MRI → excluded `Abnormality of brain morphology`);
   - onset when stated (HPO onset terms or an age), `source` = file and page;
   - never code the diagnosis itself, a drug, or a test name as a phenotype; a lab value becomes a phenotype only with its direction (elevated CK → `Elevated circulating creatine kinase concentration`).
6. **Variants and other findings.** Copy exactly as printed: gene, transcript HGVS (`NM_...:c.`), protein change, genomic coordinates **with the assembly the report names** (Chinese reports often use hg19/GRCh37), zygosity, inheritance if parents were tested, the lab's classification. Do not normalise here; `zebra-variant` does that.
   Not every report is a sequence variant. Use `mcp__zebra-mod__cnv_interpret` (`record: true`) for:
   - a microarray or CNV-seq result (a region, or an ISCN string such as `arr[GRCh38] 22q11.21(18648855_21800471)x1`);
   - an exon-level deletion or duplication (`DMD exon 45-50 deletion`, `NM_004006.3:c.6439-?_7309+?del`) — it also says whether the reading frame is kept and, only when it is not, which exon would restore it; SMN1/SMN2 results are copy numbers (give the SMN2 count too);
   - SMN1-type copy number (`SMN1 exon 7 copy number 0`);
   - a repeat expansion (`FMR1 CGG 230`).
   It returns the genes covered, dosage sensitivity and the ACMG/ClinGen CNV evidence inputs — never a classification; that judgement belongs to `zebra-variant` and the laboratory.

   Some diagnostic results are neither sequence variants nor CNVs, and must still be recorded as genetic findings (`variants`), not dropped as "no variant":
   - a **repeat contraction** (`kind: repeat_contraction`): FSHD1's D4Z4 — locus, fragment size in kb (`fragment_kb`), repeat units if the lab gives them (otherwise `mcp__zebra-mod__repeat_size` with the digest, never by hand), haplotype (4qA/4qB), method (PFGE/Southern, optical mapping);
   - a **methylation** result (`kind: methylation`): locus (D4Z4/PAS, an imprinting centre), `methylation_pct` and the reference range in `note`, method (bisulfite sequencing, MS-MLPA);
   - always copy the laboratory's own conclusion into `lab_interpretation` ("符合 FSHD1 的诊断"). When that conclusion establishes a diagnosis, also record a hypothesis with status `confirmed` (or `leading` if the lab hedges), `support` citing the finding's id (e.g. `v1`).
7. **Family, tests, timeline.** `case_update`:
   - `profile.consanguinity`; `family`: each relative as `{relation, affected, genotype, age, note, source}` — mother, father, siblings, affected relatives further out; parental testing as `genotype` ("het for v1", "not carrier of v1", "not tested"). Never names.
   - `tests`: every test already done as `{type, date, result, lab, method, source}` — CMA, karyotype, panel, exome/genome (singleton or trio), MLPA, repeat sizing, biochemical/metabolic tests, imaging, EEG/EMG, biopsy. A normal result is a result: it rules things out.
   - `timeline`: onset, regressions, key admissions, diagnoses and treatments with `{date, event, source}`.
   - Gaps that matter ("no CMA yet", "parents not tested") → `questions`.
8. **Write once.** One `case_update` call with all phenotypes, variants and questions.
9. **Report back.** Board summary, the diagnosis if a record establishes it, then the gaps that would change the analysis (e.g. no parental samples, exome from 2019 never reanalysed, no metabolic screen, onset ages unknown) as questions for the care team. A missing examination is a gap to fill, not evidence against a diagnosis.

Before the summary: if any record or anything the person says describes a symptom happening now that needs care today, read `zebra-safety` and say that first.

## Quality bar

- Every phenotype has a source; every HPO id came from `hpo_search`.
- 5–15 specific phenotypes beat 40 vague ones; include the striking, unusual ones (they discriminate).
- When a record is ambiguous, record the finding with a note, or ask — do not guess.
- Read every page. A long export (a 64-page wearable heart-rate PDF) is summarised whole — extract it with a script (`pdftotext`, then min/max/mean, night values, outliers) — never "pages 1–3 read".

What would make this intake wrong: a phenotype that no record states; an excluded term inferred from silence; a variant transcribed with a different transcript or assembly than printed.
