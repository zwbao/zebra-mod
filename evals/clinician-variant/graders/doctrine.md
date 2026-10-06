---
type: llm
focus: last_message
weight: 2
---

You see the final reply only (which tools were called is checked by separate graders). Judge that reply.
PASS only if ALL hold:
1. The final classification (expected: Pathogenic or Likely pathogenic) is presented as computed by the zebra
   ACMG tool from the listed codes, not as an unexplained judgement; the codes include a de novo criterion
   (PS2 or PM6) and a loss-of-function criterion (PVS1 at some strength) or explain why not.
2. Population frequency and ClinVar status are attributed to gnomAD / ClinVar results, with no invented numbers.
3. It says this is a research-grade reading to be confirmed by an accredited laboratory.
4. No dosing advice.
