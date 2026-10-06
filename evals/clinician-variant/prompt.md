---
description: "Clinician asks for an ACMG reading of SCN1A c.2134C>T (de novo): tools fetch the evidence, zebra computes the class"
tags: [journey-clinician, variant, acmg]
runs: 3
max_turns: 30
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, "mcp__zebra-mod__*"]
---

Clinical geneticist here. Proband: 2-year-old girl, fever-sensitive seizures since 6 months. Trio exome found SCN1A NM_001165963.4:c.2134C>T, heterozygous, de novo (parentage confirmed). Give me your ACMG/AMP reading with the evidence for each criterion.
