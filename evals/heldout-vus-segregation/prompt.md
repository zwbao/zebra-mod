---
description: "Held-out: a cardiologist asks for an ACMG/AMP reading of a novel MYH7 rod-domain missense VUS (NM_000257.4:c.4322C>A p.Ala1441Asp) that segregates with HCM in three affected relatives, with one unaffected adult sister who does not carry it, and whether it can be used for predictive testing. Ground truth: PP1 at supporting (or at most moderate) from 2-3 informative meioses, unaffected non-carrier not counted for an age-dependent-penetrance disease, PM2 at supporting, no PM1 outside the head domain; remains VUS; no predictive testing on a VUS; relatives stay under clinical screening."
tags: [heldout, en, clinician, variant]
runs: 8
max_turns: 30
timeout_seconds: 900
allowed_tools: [Read, Glob, Grep, Skill, WebSearch, WebFetch, "mcp__zebra-mod__*"]
---

Cardiologist (inherited cardiac conditions clinic). I'd like a careful ACMG/AMP reading of a variant and advice on what we can and cannot do with it in the family.

Proband: man, 38, hypertrophic cardiomyopathy, asymmetric septal hypertrophy 21 mm, no LVOT obstruction at rest, no hypertension. Diagnosed at 31.

Genetic test (clinical HCM panel, 2025): MYH7 NM_000257.4:c.4322C>A p.(Ala1441Asp), heterozygous, reported as a variant of uncertain significance. GRCh38 chr14:23417534 G>T. No other P/LP/VUS in the sarcomere genes on the panel. The lab report says the variant is not observed in population databases and not previously reported in the literature.

Family, all by targeted testing since:
- Father, 64: HCM (septum 18 mm), carries the variant.
- Brother, 33: HCM (septum 16 mm, diagnosed on family screening), carries the variant.
- Paternal aunt (father's sister), 67: HCM diagnosed in her 50s, carries the variant.
- Sister, 44: normal ECG, echo and cardiac MRI last year; does NOT carry the variant.
- Mother, 65: unaffected, not tested. Paternal grandparents deceased, cause unknown.

Questions:
1. Which ACMG/AMP criteria apply and at what strength, and what is your overall classification? Please show how you counted the segregation evidence.
2. Can we offer predictive (cascade) genetic testing with this variant to the proband's children (8 and 11) and the brother's child (4)?
3. Can the sister be discharged from cardiac screening because she doesn't carry it?
4. What would move this variant to likely pathogenic?
