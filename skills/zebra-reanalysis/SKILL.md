---
name: zebra-reanalysis
description: Reanalyse exome or genome data (VCF, singleton or trio) for a rare-disease patient on this machine — check build and samples, filter by inheritance model (de novo, homozygous, compound heterozygous, hemizygous X, dominant), population frequency and consequence, rank candidates by phenotype fit, then hand the best to variant interpretation; for researchers, prepare novel-gene candidates for matchmaking. Triggers: 重分析, reanalysis, VCF, exome, genome, WES, WGS, 全外显子, 全基因组, trio, 家系, candidate genes, 候选基因, GeneMatcher, novel gene.
---

# zebra-reanalysis — VCF → candidates

The VCF never leaves the machine: only the candidate variants are queried (VEP, gnomAD). A VCF identifies the patient and relatives — do not upload it, attach it, or paste it.

## Steps

1. **Inspect**: `zebra vcf inspect <vcf>` (Bash). Samples, variant count, build guessed from contig lengths/header, whether genotypes are present. Build unknown → ask; GRCh37 and GRCh38 are both handled, never mixed.
2. **Pedigree**: which sample is the proband, mother, father; sexes; affected status of others. Missing parents → singleton analysis (de novo cannot be called; say so).
3. **Triage**: 
   ```
   zebra --case <case> vcf triage <vcf> --proband P --mother M --father F --sex female \
       --max-af 0.01 --hpo-genes --assembly GRCh38 --out <case>/reports/triage.tsv
   ```
   It filters by quality and inheritance model, annotates the survivors with VEP (consequence, MANE, REVEL/SpliceAI where available) and gnomAD frequency, scores phenotype fit from the case's HPO profile, and writes a ranked TSV with every score component (a transparent heuristic, not an ACMG class). `--max-af` (default 0.01) applies to recessive classes; dominant classes (de novo, inherited heterozygous) use `--max-af-dominant` (default min(max-af, 0.0001)). Without a case, pass phenotypes directly with `--hpo HP:... HP:...`.
4. **Read the top candidates** (typically ≤ 20): gene–disease fit with the phenotype (`gene_card`, `disease_card`), inheritance fit, then `zebra-variant` for the best 1–5, `zebra-s2f` for splice/non-coding ones.
5. **What SNV/indel VCFs miss**: CNVs (unless a CNV VCF is given), repeat expansions, mtDNA (unless called), mosaic variants at low allele fraction, regions with poor coverage. List the ones relevant to the differential as next steps.
6. **Novel gene candidates (research)**: rare damaging variants in genes without a disease association — constraint (`gene_card`: LOEUF, missense Z), expression in the affected tissue, model-organism phenotypes, paralogs with disease. Write a matchmaking summary (gene, variant, inheritance, key HPO terms, no identifiers) for the clinician to submit to GeneMatcher / Matchmaker Exchange.
7. **Record** candidates as variants and hypotheses in the case.

## What would make this wrong

Mixed builds; parents mislabelled (check that inherited variants are shared); a de novo call from low-depth parents; filtering away a recessive allele common in the patient's ancestry; treating "no candidate" as "no genetic cause".
