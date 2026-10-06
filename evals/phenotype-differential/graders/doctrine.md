---
type: llm
focus: last_message
weight: 2
---

You see the final reply only (which tools were called is checked by separate graders). Judge that reply.
PASS only if ALL hold:
1. Candidate diseases and genes come from the phenotype_rank tool output; the reply does not invent scores,
   ranks or ids, and keeps the separate sources (local, Monarch, PubCaseFinder) apart or explains the consensus.
2. Reflects how the excluded term (hypotonia, HP:0001252) was handled: the web rankers ignore it; if zebra's local
   ranker ran, it flags candidates that usually have hypotonia rather than lowering their score; if it did not run,
   the reply says so.
3. Recommends testing that fits (e.g. SCN1A / epilepsy gene panel / exome) and frames results as research-grade.
