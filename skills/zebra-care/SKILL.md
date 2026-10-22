---
name: zebra-care
description: Care plan for a patient with an established rare diagnosis — what the records already show, severity and prognosis, guideline surveillance (respiratory, cardiac, eyes, hearing, bone), rehabilitation and exercise prescription, nutrition and metabolic issues, incidental findings that need their own plan, medicines and anaesthesia to avoid, genetics and family, psychosocial support, trial readiness and a follow-up calendar. Written for the clinician who will review and edit it, with a plain-language page for the patient. Triggers: 管理, 管理方案, 随访, 随访计划, 预后, 预测, 营养, 饮食, 康复, 运动处方, 护理, 长期管理, care plan, management, surveillance, prognosis, nutrition, rehabilitation, follow-up.
---

# zebra-care — diagnosis → care plan

The deliverable is a plan someone can act on. Use your full clinical knowledge; use the tools to check and cite what can be checked; label the rest (已检索 / 凭记忆). For a clinician this is a draft for their review: give judgments, thresholds and dose ranges as proposals, never hedge them into nothing. For a patient or family, follow `zebra-family` for register (no prescription doses; the typical course as population ranges, not a forecast), but the content still has to be useful: an adult patient asking about their own nutrition, exercise or supplements gets concrete general ranges with "confirm with your doctor".

Read `zebra-safety` first: anything urgent in the records comes before the plan.

## 1. Ground it in the records

`mcp__zebra-mod__case_status`, then read the records the plan depends on (or the case's `tests` and `variants`). Before writing, list for yourself:
- the diagnosis and what establishes it (the genetic finding and the laboratory's interpretation); genotype features that bear on severity (repeat size or units, methylation level, variant class);
- every abnormal result with its date and trend between visits (CK, creatinine, troponin, lung function, imaging activity), and results that are normal and reassuring;
- what is **not** in the records and would change the plan (examination, onset, family history, symptoms) — these become a fill-in list, not reasons to withhold the plan;
- which values you read from a photo or OCR and which you checked against the original image — mark the former "待核对".

Genotype numbers by rule: repeat units, repeat classes and risks come from `mcp__zebra-mod__repeat_size` / `rare_stats` when zebra implements the rule, never from mental arithmetic; any other derived number shows its formula and inputs. If the tool fails, give your estimate labelled unverified with the formula you used.

## 2. Disease knowledge, checked

`mcp__zebra-mod__disease_card` for the disease (inheritance, natural history, hallmark complications); `mcp__zebra-mod__literature_search` for the current guideline or consensus statement and anything decision-relevant you are unsure of (surveillance intervals, screening triggers, exercise evidence). Guideline content you know but did not retrieve is fine to use — label it.

## 3. The plan — sections

Order by what matters for this patient; drop sections that do not apply; put a one-screen **decision summary** (the 3–7 actions, each with an urgency tag: urgent / 1–3 months / routine) first.

1. **Severity and prognosis.** Factors pointing to a more and a less severe course, from genotype, age at involvement, organs involved and trend. Give prognosis as ranges and as the objective measures to track (lung function, timed walks, strength scores, quantitative muscle MRI, biomarkers), not a single percentage.
2. **Respiratory.** What to measure (the measures the disease's guideline names, adapted to what the patient can perform), thresholds for non-invasive ventilation and cough assistance, vaccines, smoking.
3. **Cardiac.** Baseline and interval ECG/echo/Holter for the disease; any abnormal marker gets its own differential (a raised troponin T in a myopathy may be skeletal — check cardiac troponin I).
4. **Exercise and rehabilitation.** Gate it on any open cardiac question (start low until it is answered). A concrete prescription (type, frequency, duration, intensity, anchored to the person's own measurements where the records have them), what to avoid, orthoses, surgery criteria and their contraindications, falls, pain, fatigue.
5. **Nutrition and metabolism.** Energy and protein targets in grams, vitamin D/calcium, weight and body composition, and a plan for each abnormal metabolic result (uric acid, homocysteine, glucose, lipids). Every target names its guideline and the tier of that guideline this patient is in, chosen from what the records document — not from what the patient might have. When the deciding facts are unknown (blood pressure, BMI, renal function, prior attacks), give the default tier for an undocumented patient and say which finding would move them to the stricter one. Quoting the strictest tier for a patient the records do not place there is an error, not caution. Supplements with their evidence level, including "not worth paying for".
6. **Findings the diagnosis does not explain.** Each gets its own row: what it could be, the next step, the urgency. Never fold an unexplained finding (a spinal cord signal, a tumour marker) into the main diagnosis. Think of the common and the treatable before the exotic, and of the dangerous before the benign — and of the patient's age, sex and exposures: a marker or imaging finding outside the disease's usual pattern gets the work-up it would get in anyone else, at the pace its worst plausible cause demands; an assay that cross-reacts with the diseased tissue gets the organ-specific alternative before it is attributed to the disease.
7. **Screening the genotype calls for** (eyes, hearing, others the guideline names) and why.
8. **Medicines and procedures.** Anaesthesia hazards, drugs to avoid or use with care, an emergency card.
9. **Genetics and family.** Inheritance and recurrence (`zebra-stats` for numbers), testing of relatives and what it would change, reproductive options.
10. **Psychosocial, work, insurance and rare-disease registration** (China: 罕见病目录, 医保, network hospitals — from the tools, `zebra-family` resources).
11. **Treatments and trials.** Hand off to `zebra-therapy` for approved drugs, trials (with China sites) and repurposing; summarise its tiers here and add trial-readiness steps (baseline measures, records to collect, things that would exclude the patient later).
12. **Follow-up calendar.** Now / 1–3 months / every 6 months / yearly / every 2–3 years / see a doctor at once if.
13. **Fill-in list.** What the patient or clinician should supply (examination, history, originals of photographed reports), at the end.

## 4. Deliver

A recommendation class or level (IIa, Class I) is stated only with its source, and an ECG or lab threshold only as the guideline writes it; when unsure, write "consider" rather than an invented class.

The deliverable cites what the clinician can check — PMIDs, NCT numbers, GeneReviews, guideline names — and never zebra's machinery (skill or tool names, `/zebra` commands, ledger ids like E12): those mean nothing outside this terminal. That includes the closing lines: an offer to save the case or export the plan is made in words ("I can export this as Word/PDF"), never as a command.

- Clinician: the plan as written, technical, editable. Offer `zebra-report` to export it (Word/PDF/HTML); for a page the clinician will edit, HTML with the uncertain values highlighted works well.
- Patient page: one page of plain-language actions (what to do, what to watch for, when to call), in their language.
- Record new questions in the case (`case_update` → `questions`) and treatment leads (`leads`).

## What would make this plan wrong

A value copied from OCR without a mark; a recommendation that contradicts a guideline without saying so; an unexplained finding attributed to the diagnosis; a section of generic advice that ignores this patient's numbers; withholding a judgment the clinician asked for.
