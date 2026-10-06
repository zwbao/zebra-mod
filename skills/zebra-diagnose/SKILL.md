---
name: zebra-diagnose
description: Phenotype-driven differential diagnosis for a rare disease — rank diseases and genes from HPO phenotypes with three independent methods, build a discriminating-features table for the leading candidates, and say which genetic or biochemical test would detect each one (including what a negative exome misses). Use when symptoms are unexplained or a diagnosis is in doubt. Triggers: 可能是什么病, differential diagnosis, 鉴别诊断, undiagnosed, 未确诊, what test should we ask for, 该做什么检查, exome negative, 全外显子阴性.
---

# zebra-diagnose — phenotypes → differential → tests

You interpret rankings; you never re-rank by intuition. A disease the methods miss enters only with retrieved evidence.

## Steps

1. **Profile.** `mcp__zebra-mod__case_status`. Fewer than 3 specific present phenotypes → run `zebra-intake` or ask for the striking features first. Note sex and age of onset; they gate candidates.
2. **Rank, three ways.** `mcp__zebra-mod__phenotype_rank` with `from_case: true`, `sources: ["local","monarch","pubcasefinder"]`, `top: 20`. Each source ranks on its own:
   - candidates high in two or three sources are the strongest signal;
   - read `matches`: which query terms the disease explains exactly, which only through a broad ancestor (weak), which not at all;
   - a term id that came back without a label: verify it with `hpo_search` before naming or dropping it; genes and diseases you mention come from these outputs and the disease cards, never from memory;
   - `excluded_hits`: diseases that usually have a feature the patient is recorded as *not* having. The local score does not count them (on the benchmark that ranked better on average), so you must: a hallmark feature excluded (seen in most patients with the disease, at this age) argues strongly against that candidate and goes in your discriminating table; an occasional feature excluded argues little;
   - `ties_at_top`: diseases with the same score are not ordered by the ranking — never present one of them as "first";
   - how good the ranking is (docs/BENCHMARK.md, phenopacket-store): for new patients the correct disease was in the local top 10 about 28% of the time and the three sources together do better — the list is a set of hypotheses to test, not an answer.
3. **Characterise the leading 5–10.** `mcp__zebra-mod__disease_card` for each: inheritance, onset, prevalence, genes, hallmark features. Inheritance and sex must fit the family (an X-linked recessive disease in a girl needs an explanation).
4. **Discriminating table.** Rows = candidates; columns = features that separate them (present / absent / unknown in this patient, from the case and the disease cards). Unknown cells that separate the top candidates become questions for the care team (`case_update` → `questions`).
5. **Which test finds each candidate.** Match the test to the usual mechanism:

   | Mechanism | Detected by | Missed by |
   |---|---|---|
   | SNVs/small indels across many genes | exome or genome, trio best | single-gene tests, CMA |
   | CNVs, aneuploidy | CMA, genome | most exome pipelines (small CNVs) |
   | repeat expansions (FMR1, DMPK, HTT, FXN, ATXN*, C9orf72, RFC1) | targeted repeat assay, long-read | exome, CMA |
   | imprinting / methylation (Prader-Willi, Angelman, Beckwith-Wiedemann, Silver-Russell) | methylation-specific MLPA | sequencing alone |
   | mitochondrial DNA | mtDNA sequencing (tissue matters: muscle, urine) | blood exome |
   | deep intronic, structural, regulatory | genome ± RNA-seq of an expressing tissue | exome |
   | somatic mosaic (e.g. PIK3CA overgrowth, many neurocutaneous) | deep sequencing of affected tissue | blood tests |
   | inborn errors of metabolism | plasma amino acids, acylcarnitines, urine organic acids, enzyme assays — often faster than genomics | — |

   Exome negative: reanalysis if data ≥ 1 year old (`zebra-reanalysis`), trio if singleton, genome, RNA-seq, long-read, and re-phenotyping (new features appear with age).
6. **Record.** `case_update` → `hypotheses` for each candidate kept: disease name, ids as the tools returned them, status (`leading` / `considered` / `excluded`), `support` and `against` ledger ids, a one-line note.
7. **Say it.**
   - Family: "possibilities the doctors may want to consider", why each fits, what test would answer it, which questions to bring. Never "your child has X".
   - Clinician: the table, scores per source, matched/unmatched features, tests, evidence ids.

A worsening on a specific drug (carbamazepine in Dravet, valproate in a urea-cycle disorder) is both a diagnostic clue and a hazard: see `zebra-safety`.

## What would falsify a leading hypothesis

A hallmark feature confirmed absent; inheritance incompatible with the family; the decisive test negative with adequate coverage for that mechanism. State it for each leading candidate.
