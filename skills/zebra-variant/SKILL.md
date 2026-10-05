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
3. **Codes from numbers.** `mcp__zebra-mod__acmg` `mode: suggest` (with `inheritance`). It returns BA1 / BS1 / PM2_Supporting and PP3 / BP4 with their calibrated thresholds (REVEL: Pejaver 2022; SpliceAI: Walker 2023). For BS1 with a known prevalence, run `rare_stats` `maxaf`. Check coverage before trusting absence.
4. **Codes from judgement** — apply only with evidence in hand:
   - **PVS1** null variant + LoF mechanism; ClinGen decision tree (Abou Tayoun 2018): NMD-escaping or last-exon truncations, in-frame exon skipping and alternative start codons lower it to Strong/Moderate/Supporting. Do not add PP3 for the same LoF effect.
   - **PS1 / PM5** same amino-acid change / another change at the residue classified P/LP (ClinVar ≥ 2 stars, from `variant_card` or literature).
   - **PS2 / PM6** de novo: confirmed parentage (PS2) or assumed (PM6); phenotype consistent with the gene.
   - **PS3 / BS3** well-established functional assay (Brnich 2019); default Supporting unless the assay is calibrated.
   - **PS4** enrichment in cases or several unrelated affected probands (literature, `zebra-literature`).
   - **PM1** hotspot or critical domain without benign variation (cap PM1 + PP3 at Strong — zebra applies it).
   - **PM3 / BP2** recessive: in trans with a P/LP variant (needs phase: parents or reads); in cis argues benign.
   - **PM4 / BP3** in-frame length change outside / inside a repeat.
   - **PP1 / BS4** segregation → `rare_stats` `segregation` (counted meioses); non-segregation in an affected relative → BS4.
   - **PP4** phenotype highly specific for the gene (use sparingly).
   - **BS2** observed in healthy adults incompatible with penetrance.
   - **BP7** synonymous or deep intronic with SpliceAI ≤ 0.1.
   Splice region, deep intronic, UTR or promoter variants → `zebra-s2f` before deciding PP3/BP4/PVS1.
5. **Classify.** `mcp__zebra-mod__acmg` `mode: classify` with your codes. Report both readings (points; 2015 rules) when they differ.
6. **Record.** `case_update` → `acmg: [{variant_id, codes, note}]`.
7. **Present** a table: code | strength | evidence (ledger id) | why. Then:
   - the class, points, and the single piece of evidence most likely to move it (parental testing, RNA study, segregation, a functional assay);
   - if it differs from the laboratory's class: "worth asking the laboratory to review", listing the evidence — never "the lab is wrong"; labs reclassify on request.

## Not a sequence variant

A CNV, an exon-level deletion or duplication, a copy-number result or a repeat expansion is scored by a different framework (ACMG/ClinGen CNV 2019, or the gene's own repeat thresholds), not by the codes above. Run `mcp__zebra-mod__cnv_interpret` for the genes covered, dosage sensitivity (ClinGen haploinsufficiency/triplosensitivity) and the section 1–5 evidence inputs, then do the scoring yourself and say which framework you used. For an out-of-frame exon deletion, the frame arithmetic it returns is also the exon-skipping question — carry it to `zebra-therapy`.

## Families

A VUS is not a diagnosis and not a negative result: it is "not enough evidence yet". Say what evidence could change it and how (parents' samples, other affected relatives, new publications). Do not let a family act on a VUS.

## What would make this wrong

A transcript or assembly mismatch; a code applied twice for one fact; PVS1 in a gene without established LoF mechanism; absence in gnomAD at a site with poor coverage; PM2 for a recessive allele that is common in the patient's ancestry group.
