---
name: zebra-variant
description: Interpret a genetic variant for a rare disease with the ACMG/AMP framework and ClinGen refinements — anchor it (transcript, assembly, consequence), weigh population frequency, computational and splicing predictions, gene–disease validity and mechanism, de novo/segregation/phase/functional evidence, then classify with points (zebra computes) and state what would change the class. Research-grade. Triggers: 变异解读, what does this variant mean, VUS, 意义不明, ACMG, pathogenic, likely pathogenic, 致病性, reclassification, 重新分类, ClinVar.
---

# zebra-variant — evidence → ACMG codes → points

Several variants (a compound heterozygote, a reanalysis shortlist) → one `zebra-mod:variant-curator` subagent per variant, all in one message; then reconcile phase and inheritance across them yourself.

You justify codes; `zebra` does the arithmetic. Every code cites a ledger id or names the missing evidence. The result is research-grade until an accredited laboratory confirms it.

## Steps

1. **Anchor.** `mcp__zebra-mod__variant_card` with transcript HGVS (`NM_...:c.`) or coordinates **plus assembly**. Check:
   - the normalised HGVS matches the report (same transcript, same change); a 3′-shift or a different transcript is fine only if you say so; a mismatch → stop and report it;
   - the MANE Select consequence; NMD relevance for truncating variants (last exon / last 50 nt of the penultimate exon escape NMD).
2. **Gene and mechanism.** `mcp__zebra-mod__gene_card`:
   - ClinGen validity for the disease in question — Disputed/Refuted/No known: no P/LP for that disease; Limited: say so;
   - inheritance; mechanism (loss of function: haploinsufficiency score 3, low LOEUF; gain of function / dominant negative: missense clustering, few LoF in patients) — PVS1 needs LoF to be the established mechanism.
3. **Codes from numbers.** `mcp__zebra-mod__acmg` `mode: suggest` (with `inheritance`; for BS1 also `prevalence`, `allelic`, `genetic`, `penetrance` — AD/AR pick the Whiffin formula themselves, XLR/XLD/unknown need `inheritance_mode`; a VCEP's gene-specific PM2 ceiling goes in `pm2_max_af`). It returns BA1 / BS1 / PM2_Supporting and PP3 / BP4 / BP7 with their calibrated thresholds (REVEL: Pejaver 2022; SpliceAI at ±4,999 nt from the SpliceAI lookup: Walker 2023), `codes_after_cap` (PM1 + PP3 capped at 4 points), and two blocks of **inputs, never codes**: `pvs1_inputs` (NMD prediction on the real exon structure, fraction of the protein lost, canonical-site exon and whether skipping it keeps the frame, ClinGen HI, LOEUF) and `ps1_pm5_inputs` (ClinVar P/LP records at the same codon, with stars). PM2 comes only from absence at a covered site, or a recessive frequency under the stated convention; BS1 and PM2 are never offered together. GRCh37 input is mapped to GRCh38 and both builds are shown. Check coverage before trusting absence.
4. **Codes from judgement** — apply only with evidence in hand:
   - **PVS1** null variant + LoF mechanism; ClinGen decision tree (Abou Tayoun 2018) walked with `pvs1_inputs`: NMD-escaping or last-exon truncations, in-frame exon skipping and alternative start codons lower it to Strong/Moderate/Supporting. Do not add PP3 for the same LoF effect.
   - **PS1 / PM5** same amino-acid change / another change at the residue classified P/LP — from `ps1_pm5_inputs` (ClinVar review status is reported as "n of 4 stars"; ≥ 2 stars is the usual floor).
   - **PS2 / PM6** de novo: confirmed parentage (PS2) or assumed (PM6); phenotype consistent with the gene.
   - **PS3 / BS3** well-established functional assay (Brnich 2019); default Supporting unless the assay is calibrated. `variant_card` returns MaveDB scores where a dataset covers the variant, labelled only by that dataset's own published calibration (with its OddsPath); a construct-numbered dataset or one without a calibration gives no class.
   - **PS4** enrichment in cases or several unrelated affected probands (literature, `zebra-literature`).
   - **PM1** hotspot or critical domain without benign variation (cap PM1 + PP3 at Strong — zebra applies it).
   - **PM3 / BP2** recessive: in trans with a P/LP variant (needs phase: parents or reads); in cis argues benign.
   - **PM4 / BP3** in-frame length change outside / inside a repeat.
   - **PP1 / BS4** segregation → `rare_stats` `segregation` (counted meioses: `ad_meioses`, `xlr_male_meioses`, `ar_affected_sibs`); an affected relative without the variant → `nonsegregations` (BS4); unaffected relatives count only with `full_penetrance: true`.
   - **PP4** phenotype highly specific for the gene (use sparingly).
   - **BS2** observed in healthy adults incompatible with penetrance.
   - **BP7** synonymous or deep intronic with SpliceAI ≤ 0.1.
   Splice region, deep intronic, UTR or promoter variants → `zebra-s2f` before deciding PP3/BP4/PVS1.
   **mtDNA** (`m.3243A>G 35%`, or `heteroplasmy` in percent): gnomAD's mitochondrial counts (homoplasmic/heteroplasmic, maximum heteroplasmy) and MITOMAP's disease association come back; the nuclear PM2/BS1/BA1 rules are not applied — use the mitochondrial specifications (McCormick 2020) and say the heteroplasmy in the tested tissue matters.
   **Literature counts** (LitVar) exclude records spelled as a different substitution at the same position and say how many were excluded: a PMID count is for this allele.
5. **Classify.** `mcp__zebra-mod__acmg` `mode: classify` with your codes. Report both readings (points; 2015 rules) when they differ.
6. **Record.** `case_update` → `acmg: [{variant_id, codes, note}]`.
7. **Present** a table: code | strength | evidence (ledger id) | why. Then:
   - the class, points, and the single piece of evidence most likely to move it (parental testing, RNA study, segregation, a functional assay);
   - if it differs from the laboratory's class: "worth asking the laboratory to review", listing the evidence — never "the lab is wrong"; labs reclassify on request.

## Not a sequence variant

A CNV, an exon-level deletion or duplication, a copy-number result or a repeat expansion is scored by a different framework (ACMG/ClinGen CNV 2019, or the gene's own repeat thresholds), not by the codes above. Run `mcp__zebra-mod__cnv_interpret` for the genes covered, dosage sensitivity (every gene checked against ClinGen), the ClinGen curated regions it overlaps with their HI/TS scores and the section-2 row they point to (e.g. complete overlap of an HI 3 region → 2A), section 3 by gene count, and recessive genes whose loss means carrier status; then do the scoring yourself (sections 4–5 need the family and the literature) and say which framework you used. Repeat sizes come back placed in the published bands with their sources — a referenced reading, not a classification; methylation, AGG interruptions or the other allele may be needed first. For an out-of-frame exon deletion, the frame arithmetic it returns is also the exon-skipping question — carry it to `zebra-therapy`.

## Families

A VUS is not a diagnosis and not a negative result: it is "not enough evidence yet". Say what evidence could change it and how (parents' samples, other affected relatives, new publications). Do not let a family act on a VUS.

## What would make this wrong

A transcript or assembly mismatch; a code applied twice for one fact; PVS1 in a gene without established LoF mechanism; absence in gnomAD at a site with poor coverage; PM2 for a recessive allele that is common in the patient's ancestry group.
