---
type: regex
target: last_message
pattern: 'HP:\d{7}'
match: not_contains
weight: 1
---

<!-- A somatic tumour report is not a rare-disease phenotype problem: no HPO terms in the reply.
     AMP/ASCO/CAP somatic guideline, Li et al. J Mol Diagn 2017;19:4-23. -->
