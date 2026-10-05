---
name: zebra-literature
description: Literature work for rare disease with verifiable citations — GeneReviews and reviews first for a known disease, Europe PMC searches with phenotype and gene synonyms, papers that mention a specific variant (LitVar/PubTator), case reports and series for natural history and treatment, each claim quoted with the sentence it rests on and its PMID as returned by a tool. Triggers: 文献, literature, papers, PubMed, 有没有报道, case report, 病例报告, review, GeneReviews, evidence review.
---

# zebra-literature — retrieve, read, then cite

A PMID appears in an answer only if a tool returned it in this session; a number from a paper appears only with the sentence it came from. Abstract-only reading is said as such.

## Steps

1. **Frame the question**: disease overview / variant reports / treatment experience / natural history / mechanism. One question per search round.
2. **Known disease → GeneReviews first** (`disease_card` returns the chapter when one exists), then recent reviews.
3. **Search**: `mcp__zebra-mod__literature_search`:
   - `query` with synonyms and the gene symbol (Europe PMC syntax: `"Dravet syndrome" OR "SCN1A"`; add `AND (case report)`, `AND PUB_YEAR:[2020 TO 2026]`);
   - `gene` + `variant` (rsID or HGVS protein, e.g. `p.Arg712*`) for variant-level reports via LitVar/PubTator.
4. **Read** the relevant hits: open-access full text (`fullTextUrl` from the result) with WebFetch, else the abstract. Extract into a table: PMID | design (case report / series / cohort / functional / trial) | n | genotype | finding (quoted) | relevance.
5. **Many papers** → `zebra-mod:literature-scout` subagents, one per question, all in one message; each returns the table with quotes.
6. **Synthesise** by evidence strength, newest relevant first; mark single case reports as such; separate "reported" from "established".
7. **Record**: cite PMIDs in the answer and in case hypotheses/leads (the ledger keeps the search; add the PMID to `support`).

## What would make this wrong

A PMID from memory; a claim attributed to a paper that only cites it; a review's summary presented as primary data; an animal-model result presented as human evidence; missing the year (rare-disease knowledge moves fast).
