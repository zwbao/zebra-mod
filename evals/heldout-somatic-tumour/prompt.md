---
description: "Held-out: an oncologist pastes a tissue NGS report for stage IV lung adenocarcinoma (EGFR exon 19 deletion VAF 27%, TP53 R273H VAF 29%, BRCA2 frameshift VAF 49% at 40% tumour content; PD-L1 TPS 60%; small brain metastases) and asks what is actionable. Ground truth: somatic framing (AMP/ASCO/CAP or ESCAT tiers), first-line osimertinib-based therapy, no first-line immunotherapy despite PD-L1, germline follow-up for BRCA2 only, no germline ACMG of somatic variants, no rare-disease phenotype ranking."
tags: [heldout, en, clinician, oncology, somatic]
runs: 8
max_turns: 30
timeout_seconds: 900
allowed_tools: [Read, Glob, Grep, Skill, WebSearch, WebFetch, "mcp__zebra-mod__*"]
---

Medical oncologist. Could you go through this tumour NGS report and tell me what is actionable and how you'd treat first line? De-identified.

Patient: woman, 61, never-smoker. Lung adenocarcinoma, stage IVA (cT2a N2 M1a: malignant pleural effusion) plus two asymptomatic brain metastases (6 mm and 8 mm) on MRI, so overall stage IVB. ECOG 1. No prior treatment.

Tissue NGS, 520-gene DNA panel + RNA fusion panel, tumour-only; estimated tumour content 40%.

| Gene | Variant | VAF |
|---|---|---|
| EGFR | NM_005228.5:c.2235_2249del p.(Glu746_Ala750del), exon 19 | 27% |
| TP53 | NM_000546.6:c.818G>A p.(Arg273His) | 29% |
| BRCA2 | NM_000059.4:c.5946del p.(Ser1982Argfs*22) | 49% |

- No KRAS, BRAF, ERBB2 or MET exon 14 alterations. No ALK, ROS1, RET or NTRK fusions on RNA.
- TMB 3 mut/Mb. MSI stable.
- PD-L1 (22C3) TPS 60%.

Family history: her sister had breast cancer at 46; father had prostate cancer in his 70s.

What is actionable, what would you give first line, and is there anything else I should do?
