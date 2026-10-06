---
name: zebra-reanalysis
description: Reanalyse exome or genome data (VCF, singleton or trio) for a rare-disease patient on this machine — check build and samples, filter by inheritance model (de novo, homozygous, compound heterozygous, hemizygous X, dominant), population frequency and consequence, rank candidates by phenotype fit, then hand the best to variant interpretation; for researchers, prepare novel-gene candidates for matchmaking. Triggers: 重分析, reanalysis, VCF, exome, genome, WES, WGS, 全外显子, 全基因组, trio, 家系, candidate genes, 候选基因, GeneMatcher, novel gene.
---

# zebra-reanalysis — VCF → candidates

The VCF file never leaves the machine, but variant positions do: by default the candidates that survive local filtering go to VEP and gnomAD, and up to `--s2f-top` (default 5) splice/non-coding candidates to the SpliceAI lookup. `--prefilter myvariant` (needed for a whole exome or genome) sends **every** quality-passing variant to myvariant.info — taken together that is the person's genome, so the mod asks first; say so before using it. `result.sent_off_machine` reports what reached each host. A VCF identifies the patient and relatives — do not upload it, attach it, or paste it.

## Steps

1. **Inspect**: `zebra vcf inspect <vcf>` (Bash). Samples, variant count, build guessed from contig lengths/header, whether genotypes are present. Build unknown → ask; GRCh37 and GRCh38 are both handled, never mixed.
2. **Pedigree and QC**: which sample is the proband, mother, father; sexes; affected status of others (a PED file: `--ped`). Then `zebra qc <vcf> --ped <ped>` (or `--proband/--mother/--father`): sex from genotypes against the stated sex, KING kinship (sample swaps, non-paternity, related parents), runs of homozygosity (consanguinity, regions for recessive candidates), Mendelian errors and uniparental disomy by chromosome with the imprinting disorder it would mean, mosaic de novo calls. A swap or parentage problem stops the analysis until it is resolved. Missing parents → singleton analysis (de novo cannot be called; say so).
3. **Triage**: 
   ```
   zebra --case <case> vcf triage <vcf> --proband P --mother M --father F --sex female \
       --max-af 0.01 --hpo-genes --assembly GRCh38 --out <case>/reports/triage.tsv
   # a whole exome or genome: add --prefilter myvariant (asks first: it sends every variant position)
   ```
   A VCF holding only CNV/SV calls is refused with a ready `zebra cnv "chrN:a-b loss"` command for each call → `mcp__zebra-mod__cnv_interpret`.
   It filters by quality and inheritance model, annotates the survivors with VEP (consequence, MANE, REVEL/SpliceAI where available) and gnomAD frequency, scores phenotype fit from the case's HPO profile, and writes a ranked TSV with every score component (a transparent heuristic, not an ACMG class). `--max-af` (default 0.01) applies to recessive classes; dominant classes (de novo, inherited heterozygous) use `--max-af-dominant` (default min(max-af, 0.0001)). Without a case, pass phenotypes directly with `--hpo HP:... HP:...`.
4. **Read the top candidates** (typically ≤ 20): gene–disease fit with the phenotype (`gene_card`, `disease_card`), inheritance fit, then `zebra-variant` for the best 1–5, `zebra-s2f` for splice/non-coding ones.
5. **What SNV/indel VCFs miss**: CNVs (from a CMA or a CNV caller — interpret them with `cnv_interpret`), repeat expansions, mtDNA (unless called), mosaic variants at low allele fraction, regions with poor coverage. List the ones relevant to the differential as next steps.
6. **Novel gene candidates (research)**: rare damaging variants in genes without a disease association — constraint (`gene_card`: LOEUF, missense Z), expression in the affected tissue, model-organism phenotypes, paralogs with disease. Write a matchmaking summary (gene, variant, inheritance, key HPO terms, no identifiers) for the clinician to submit to GeneMatcher / Matchmaker Exchange.
7. **Record** candidates as variants and hypotheses in the case.

## What would make this wrong

Mixed builds; parents mislabelled (check that inherited variants are shared); a de novo call from low-depth parents; filtering away a recessive allele common in the patient's ancestry; treating "no candidate" as "no genetic cause".
