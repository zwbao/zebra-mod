---
type: llm
focus: last_message
weight: 2
---

You see the final reply only (which tools were called is checked by separate graders). Judge that reply.
PASS only if ALL hold:
1. The reply does not invent scores, ranks or ids for the rankers, and keeps the separate sources (local, Monarch,
   PubCaseFinder) apart or explains the consensus. Candidates added from clinical reasoning are allowed (and good when
   they fit) provided they are labelled as such and not given ranker scores. Dravet syndrome / SCN1A is named.
2. Reflects how the excluded term (hypotonia, HP:0001252) was handled: the web rankers ignore it; if zebra's local
   ranker ran, it flags candidates that usually have hypotonia rather than lowering their score; if it did not run,
   the reply says so.
3. Recommends testing that fits (e.g. SCN1A / epilepsy gene panel / exome).
