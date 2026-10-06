# Changelog

## 0.2.0 — 2026-10-06

Making the claim literally true — a rare-disease workstation for families on a diagnostic odyssey, clinicians and researchers. Three independent reviews of 0.1.0 set the list; every change below carries a regression test that fails on 0.1.0.

**Privacy and permissions**
- The mod's own tools now go through the permission chain and the privacy gate (in 0.1.0 they bypassed both: a registered name could reach Europe PMC). Read-only lookups skip the prompt unless a rule says otherwise; case writes follow the permission mode, with "allow for this session".
- The gate reads commands as the shell will (quotes, `bash -c`, backticks, process substitution, `/dev/tcp`, runners such as `uv run` and `Rscript`), matches dates of birth in every common order and spelling, never matches a record number inside a scientific id (HP:, rs, NM_, PMID, coordinates), checks names in free-form keys, asks before case files or genome data leave by scp/rsync/aws/gsutil/gh gist/pipes or a synced folder (Dropbox, iCloud, OneDrive, Nutstore, Baidu), asks before the MyVariant whole-exome prefilter, and refuses outbound calls if it crashes.
- Identifiers can be registered before a case exists (held for the session, never written, added to the next case); with a case open they are refused in every outgoing call. Without a case, an email or phone number in another tool's call is asked about, not refused.
- Installed, the mod no longer changes unrelated work: the default doctrine mode `auto` adds the full rules only while a case is open, otherwise a short section that applies to rare-disease questions only.

**Diagnosis**
- Benchmark on GA4GH phenopacket-store (10,374 cases, held-out split, leakage strata): correct disease in the local top 10 for 28% of held-out cases whose paper is not an annotation source (0.1.0: 15%), 70% across all held-out cases (an upper bound). docs/BENCHMARK.md.
- Excluded terms are flagged (`excluded_hits`), not scored; ties reported; curated NOT annotations pooled correctly; ~50 Chinese lay phrases, inside short sentences, negation-aware, also without the local HPO files.
- Phenopacket v2 export (no free text, no identifiers) and v1/v2 import; `evals/` with 10 end-to-end cases for `claude plugin eval`.

**Variants, CNVs, S2F**
- ACMG: PM2 re-derived from ClinGen SVI PM2 v1.0, Whiffin formula chosen by inheritance (BS1 for F508del no longer offered), PM1+PP3 cap kept at 4 points, BS4 with its strength; SpliceAI at ±4,999 nt in `acmg suggest`; PVS1 and PS1/PM5 inputs; GRCh37 input mapped to GRCh38; mtDNA through gnomAD's mitochondrial data, heteroplasmy and MITOMAP; LitVar counts exclude other alleles; "SCN1A c.2134C>T" read on the MANE transcript with a warning.
- CNVs: ClinGen curated regions with HI/TS and the ACMG section-2 row (22q11.2 → ISCA-37446, HI 3, 2A); every gene checked; signed intronic offsets; SMN1/SMN2 clinical exon numbering and copy-number reading; FMR1/HTT/DMPK/FXN/C9orf72 size bands with quoted sources.
- `aso_screen` (splice-switching antisense feasibility), `expression` (GTEx tissues), MaveDB scores under their own calibration.

**Reanalysis and QC**
- Whole-exome triage through a MyVariant frequency prefilter (cold 556 s, warm 3 s, resumable); SpliceAI on the top splice/non-coding candidates; PED input; CNV-only VCFs handed to `zebra cnv`.
- `zebra qc`: sex check, KING kinship, runs of homozygosity, Mendelian errors and uniparental disomy, mosaic de novo calls.

**China**
- `access`: FDA/EMA status from the agencies' records, NMPA/CDE approval from bundled official documents, the 2025 national reimbursement list with restriction text, trials with sites in China, collaboration-network hospitals by province; Chinese-cohort allele frequencies (NyuWa, WBBC, 1000G East Asian, Taiwan Biobank); stricter Chinese name matching (糖尿病 is no longer maple syrup urine disease).

**Families and follow-up**
- `report_export`: Word, PDF (local browser or LibreOffice) and HTML with Chinese typography, refusing a report that still holds an identifier; the for-a-family (代操作) flow for family letters and visit-preparation sheets.
- `case_recheck`: what changed since the last check — ClinVar, ClinGen validity, recruiting trials, new papers — into the case timeline.
- Case model: family members, tests already done, timeline, identifiers.

**Install and engine**
- One sentence ("install https://github.com/zwbao/zebra-mod") through INSTALL.md, or `install.sh`; `zebra doctor` marks optional keys as optional.
- HTTP: an error body sent with 200 is evicted, accepted not-found answers live a day, GraphQL error 500s are not retried, per-host pacing for the new sources, PubCaseFinder's hourly and daily limits counted across processes; results are trimmed before provenance.

## 0.1.0 — 2026-10-05

First release.

- **Mod** (`hooks/register.tsx`): 15 model-callable tools backed by the `zebra` CLI; research doctrine in the system prompt (`always` / `case` / `off`); case board pane and status line; `/zebra` command (`new`, `case`, `board`, `ledger`, `doctor`, `close`); auto-approval of the mod's own read-only tools; privacy gate (denies outgoing calls carrying a case's protected identifiers or ID-number/phone/email patterns, including Bash calls to the zebra CLI; asks before raw genome files leave the machine).
- **Engine** (`zebra/`, Python standard library, ≥ 3.9): HTTP layer with caching, retries and provenance; case workspace with an append-only evidence ledger; ACMG/AMP points (Tavtigian 2020) beside the 2015 combining rules, ClinGen SVI warnings and calibrated suggestions (REVEL per Pejaver 2022, SpliceAI per Walker 2023, BA1/BS1/PM2_Supporting); classical statistics; offline HPO ranking (Resnik best-match average) with English and official Chinese label search; Ensembl VEP client with AlphaMissense/REVEL/CADD/SpliceAI; sources for HPO, Monarch, PubCaseFinder, Orphanet, OLS, GeneReviews, gnomAD, ClinVar, ClinGen, PanelApp, UniProt, Europe PMC, PubTator, LitVar, ClinicalTrials.gov, Open Targets; S2F bridge (SpliceAI/Pangolin lookup; AlphaGenome, Evo 2, GPN-MSA via s2f-penguin); base-editing feasibility screen; VCF inspection and trio triage; `zebra doctor`.
- **Skills** (12, including `zebra-safety`: urgent red flags and disease-specific drug, anaesthesia and procedure hazards) and **subagents** (6).
