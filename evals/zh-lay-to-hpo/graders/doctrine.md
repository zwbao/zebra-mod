---
type: llm
focus: last_message
weight: 2
---

You see the final reply only (that hpo_search was called, and which ids came back, is checked by separate graders).
PASS only if ALL hold:
1. Each phrase is given an HPO id (HP:nnnnnnn) with its label.
2. The final reply is written in Chinese.
3. Each phrase maps to a term with the same meaning: 头围小 (small head circumference) to decreased head
   circumference or microcephaly, never an increased/large term; 听力下降 (hearing loss) to a hearing term, never a
   visual term; 肌张力低下 (low muscle tone) to hypotonia; 孩子走路晚 (walking late) to delayed ability to walk
   (motor delay may be added).
Mentioning an alternative term as an option is fine.
