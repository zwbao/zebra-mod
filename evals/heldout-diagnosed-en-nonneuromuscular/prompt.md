---
description: "Held-out: an English-speaking physician asks for a management and surveillance plan for a 29-year-old man with molecularly confirmed vascular EDS (COL3A1 glycine substitution), a proposed catheter angiogram and contact sport, plus an incidental 2.2 cm TI-RADS 4 thyroid nodule. Ground truth: diagnosis taken as established; celiprolol and BP control; non-invasive imaging only; avoid catheter angiography and contact sport; emergency card and planning; cascade testing; the nodule handled on its own terms (FNA)."
tags: [heldout, en, clinician, diagnosed, care]
runs: 8
max_turns: 30
timeout_seconds: 900
allowed_tools: [Read, Glob, Grep, Skill, WebSearch, WebFetch, "mcp__zebra-mod__*"]
---

I'm a general internist taking over the care of a young man with vascular Ehlers-Danlos syndrome and would like a concrete management and surveillance plan. De-identified summary:

- Man, 29. Molecular diagnosis 2024: heterozygous COL3A1 glycine substitution in the triple-helical domain, classified pathogenic by the lab.
- History: spontaneous sigmoid colon perforation at 24 (Hartmann's procedure, reversal at 25 complicated by wound dehiscence). Easy bruising since childhood, thin translucent skin on the chest, two spontaneous pneumothoraces in his late teens.
- Family: mother died at 43 of a "ruptured aorta"; one sister, 26, untested; he has no children but he and his partner are thinking about it.
- Current: no medications. BP 134/86 on repeated clinic readings. Plays amateur rugby at weekends and lifts weights.
- CT angiography neck-to-pelvis last month: 11 mm fusiform dilatation of the left common iliac artery, no dissection; otherwise normal. The vascular surgeon has suggested a diagnostic catheter angiogram "to characterise it better".
- Incidental on the same CT: a thyroid nodule. Dedicated ultrasound: 2.2 cm, solid, hypoechoic, smooth margins, wider-than-tall, no echogenic foci (reported ACR TI-RADS TR4). TSH normal. No neck irradiation, no family history of thyroid cancer.

Please give me the management and surveillance plan, and tell me what to do about the thyroid nodule.
