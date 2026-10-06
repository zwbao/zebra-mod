---
name: zebra-stats
description: Classical statistics for rare-disease genetics, computed offline with the model stated beside every number — cosegregation likelihood ratio and PP1 strength, maximum credible allele frequency (BS1), carrier frequency and genetic prevalence, recurrence risk (including Bayesian X-linked carrier risk), Fisher/burden tests for cohorts, de novo enrichment, Kaplan–Meier natural history with log-rank, and N-of-1 trial design and analysis. Triggers: segregation, 共分离, LOD, allele frequency, 等位基因频率, carrier frequency, 携带者频率, prevalence, 患病率, recurrence risk, 再发风险, burden test, de novo enrichment, natural history, 自然史, survival, N-of-1, 单病例试验.
---

# zebra-stats — the right test, its assumptions, its number

Run `mcp__zebra-mod__rare_stats` (`method`, `params` = the flags below without dashes) or `zebra stats <method> ...` in Bash. Always print the model line the tool returns next to the number.

| Question | Method | Key params | Assumption to state |
|---|---|---|---|
| Does the variant track with disease in the family? | `segregation` | `ad_meioses`, `xlr_male_meioses`, or `ar_affected_sibs` / `ar_unaffected_sibs`; against it: `nonsegregations` (affected relatives without the variant → BS4), `unaffected_carriers`; `full_penetrance: true` before any unaffected relative counts | no phenocopies; count only informative meioses; points per ClinGen 2024 (Biesecker) → PP1 strength; BS4 from non-segregation |
| Too common to cause this disease? | `maxaf` | `prevalence`, `allelic`, `genetic`, `penetrance`, `inheritance`, `faf95`, `an` | Whiffin 2017; prevalence and heterogeneity estimates drive it — cite where they came from |
| How many carriers / how common genetically? | `carrier` | `prevalence` or `allele_freqs` | Hardy–Weinberg, random mating, full penetrance; consanguinity and founder effects break it |
| Risk for the next child? | `recurrence` | `mode` (AR, AD-inherited, AD-de-novo, XLR-carrier-mother, XLR-bayes), `penetrance`, `mosaic`, `prior`, `unaffected_sons` | parental genotypes confirmed; de novo risk is germline mosaicism, gene-dependent |
| Enriched in my cases? | `burden` / `fisher` | carrier counts and sizes | comparable sequencing and calling in cases and controls; say so when controls are gnomAD |
| More de novo hits than chance? | `denovo` | `observed`, `trios`, `mu` | per-gene mutation rate for the variant class (cite the table); exome-wide threshold ~2.6e-6 per class |
| Natural history, event-free survival | `km` | `csv`, `time_col`, `event_col`, `group_col`, `event_coding` (`0/1` default, `1/2` as R's survival package) | non-informative censoring; small n → wide intervals; report n at risk |
| Does treatment help this one patient? | `nof1` | design: `effect`, `sd_diff`; analysis: `treatment`, `control` | exchangeable periods, washout, no carryover, stable disease |

Rules:
- Inputs come from retrieved sources (prevalence from `disease_card`, frequencies from `variant_card`) or from the person; say which.
- Small numbers: report exact tests and intervals, not just p-values; with n < 10 talk about what the data cannot exclude.
- Genetic risk numbers for families go through `zebra-family` wording and a recommendation to see a genetic counsellor.
