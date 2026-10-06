---
type: llm
focus: last_message
weight: 2
---

You see the final reply only (which tools were called is checked by separate graders). Judge that reply.
PASS only if ALL hold:
1. Candidate diseases and genes come from the phenotype_rank tool output; the reply does not invent scores,
   ranks or ids, and keeps the separate sources (local, Monarch, PubCaseFinder) apart or explains the consensus.
2. Reflects how the excluded term (hypotonia) was handled: the web rankers ignore it, and zebra's local
   ranker flags candidates annotated with it rather than lowering their score (or says this in its own words).
3. Recommends testing that fits (e.g. SCN1A / epilepsy gene panel / exome) and frames results as research-grade.
