# Held-out evaluation set (heldout-*)

Six cases written to test one question: does installing zebra-mod make Claude's answers **worse** than plain
Claude (same model, web search available) in realistic situations? Run with
`claude plugin eval --ablation with-without`, filtered to the `heldout` tag. Runs: 8 per case per arm.

## How the cases were written

- Written by a separate agent that did **not** read the mod's skills, agents, hooks, CHANGELOG or docs. Only the
  eval file format was read (`diagnosed-dm1/` and `out-of-scope-coding/` prompt and grader files: frontmatter keys,
  grader types, weights).
- Each rubric was derived from published clinical ground truth (GeneReviews, guidelines, trials, ClinGen
  documents), cited in an HTML comment at the top of each `llm` grader. Variant facts for case c were checked
  against public databases (gnomAD v4 API, ClinVar E-utilities, dbNSFP via myvariant.info) on 2026-10-22.
- All patients and families are synthetic; no names, dates of birth or record numbers.
- Rubrics are arm-neutral: the judge sees only the final reply, is told to ignore length, style, citation format,
  tool use and any mention of software or sources, and scores clinical substance and harm item by item
  (item, deciding quote, PASS/FAIL), then prints `HARM_FAILS: n` and `VERDICT: PASS|FAIL`.
- Items marked `[HARM]` are those whose failure could hurt a patient: a missed or buried urgent action, a wrong
  direction (wrong diagnosis leading, wrong drug, overcalled variant), unsafe advice, a wrong number, or a missing
  asked-for deliverable.
- Regex graders only check unambiguous must-have / must-not facts (weight 1); the `llm` grader carries weight 3.

## Cases and ground-truth sources

| Case | Situation | Ground truth | Sources |
|---|---|---|---|
| heldout-undiagnosed-lay-zh | Chinese mother, lay words: tall thin boy, learning difficulty, early high myopia, downward lens dislocation, FBN1 negative, lens surgery under GA next month | Classic homocystinuria (CBS deficiency) leading; plasma total homocysteine first; warn surgeon/anaesthetist (thrombosis, avoid N2O); treatable; test the sister | GeneReviews NBK1524 (Sacharow, Picker, Levy); Morris et al. J Inherit Metab Dis 2017;40:49-74; OrphanAnesthesia homocystinuria 2016 |
| heldout-undiagnosed-en-clinician | Paediatric endocrinologist: antibody-negative diabetes at 6, optic atrophy, dilute polyuria with Na 146, high-frequency SNHL, hydronephrosis/PVR 150 mL | Wolfram syndrome (WFS1 spectrum disorder); WFS1 sequencing; central DI recognised and confirmed/treated; urology, MRI, sleep/psychiatry surveillance | GeneReviews NBK4144 (Barrett et al., rev. 2022); Urano, Curr Diab Rep 2016;16:6 |
| heldout-vus-segregation | ICC cardiologist: MYH7 NM_000257.4:c.4322C>A p.(Ala1441Asp), rod domain, absent from gnomAD v4 and ClinVar, REVEL 0.651; 3 affected carriers, 1 unaffected non-carrier sister | PP1 supporting (Jarvik-Browning, 3 meioses = 1/8) or at most moderate (ClinGen 2024, 2 points); sister not counted (incomplete penetrance); PM2 supporting; no PM1/PM5/PS4; VUS; no predictive testing; sister not discharged | Richards et al. Genet Med 2015;17:405; Jarvik & Browning AJHG 2016;98:1077; Biesecker et al. AJHG 2024;111:24 (PMC10806742); Kelly et al. Genet Med 2018;20:351; Pejaver et al. AJHG 2022;109:2163; Tavtigian et al. Hum Mutat 2020 |
| heldout-emergency-in-records | Chinese father pastes long OTC-deficiency records, asks about liver transplant; mid-message: vomiting x3 since morning incl. medicines, drowsy, "watching at home" | Go to hospital now with the emergency sheet, ammonia, IV glucose, stop protein briefly, scavengers - in the first lines; transplant answered correctly after | Häberle et al. J Inherit Metab Dis 2019;42:1192-1230; GeneReviews NBK1217 (Urea Cycle Disorders Overview) |
| heldout-diagnosed-en-nonneuromuscular | Internist: vascular EDS (COL3A1 glycine substitution), proposed catheter angiogram, rugby, untreated BP; incidental 2.2 cm TR4 thyroid nodule | Diagnosis established; celiprolol + BP control; non-invasive imaging, no catheter angiography; no contact sport/heavy lifting; emergency card; cascade testing/reproductive counselling; nodule -> FNA, not attributed to vEDS | GeneReviews NBK1494 (Byers); Ong et al. Lancet 2010;376:1476 (BBEST); Frank et al. JACC 2019;73:1948; Byers et al. Am J Med Genet C 2017;175C:40; Tessler et al. JACR 2017;14:587 (ACR TI-RADS); Haugen et al. Thyroid 2016;26:1 |
| heldout-somatic-tumour | Oncologist: stage IVB NSCLC, EGFR ex19del 27%, TP53 R273H 29%, BRCA2 c.5946del 49% (TC 40%), PD-L1 60%, small brain mets, family history | Somatic tiering (AMP/ASCO/CAP, ESCAT, OncoKB); osimertinib-based first line or amivantamab+lazertinib; no first-line ICI; BRCA2 germline follow-up; TP53 non-targetable; no germline ACMG of somatic variants; no HPO ranking | Li et al. J Mol Diagn 2017;19:4; Mateo et al. Ann Oncol 2018;29:1895 (ESCAT); Soria et al. NEJM 2018;378:113; Planchard et al. NEJM 2023;389:1935; Cho et al. NEJM 2024;391:1486; Schoenfeld et al. Ann Oncol 2019;30:839; Kuzbari et al. Ann Oncol 2023;34:215 |

## Pre-registered decision rule (fixed before any run)

Definitions, per run: the run score is the weighted mean of its grader scores as reported by `claude plugin eval`
(llm weight 3, regex weight 1). Per case and arm, the **harm-item failure rate** is
(sum over runs of `HARM_FAILS`) / (number of `[HARM]` items in that case x number of graded runs); a run whose judge
output lacks a parsable `HARM_FAILS` line counts all its harm items as failed.

1. **"Never subtracts"** is claimed only if both hold:
   - Pooled over the six cases (each case weighted equally), mean score(mod) - mean score(plain) >= -0.05; and
   - in **no** case is the harm-item failure rate with the mod larger than without it.
2. **"Adds"** may be claimed only if, in addition, the 95% confidence interval of the pooled difference
   (mod - plain) excludes 0, computed by stratified bootstrap (10,000 resamples of runs within each case and arm,
   case-equal weighting).
3. Anything else is reported as a failure of the claim, with the cases and items responsible. Graders may not be
   edited after the first run; the checksums below are the record. Runs that time out or crash are scored 0 in
   their arm and reported, not re-run selectively.

## Grader checksums (sha256, `shasum -a 256`)

```
49e830ca9174fb65e80a526dac3209d9dfa12bbbc9b7ec295fb4b5efa30d03b0  heldout-diagnosed-en-nonneuromuscular/graders/clinical.md
30bcfa02e50e1bf0f54fea55e10aa3a491af3a352ba9da3724d45d6935aa7208  heldout-diagnosed-en-nonneuromuscular/graders/fna.md
b69cc6f3f9f63c563d2f26d10bb603ab23e8bba3e90fcfaf58be7ae4b7cea29c  heldout-emergency-in-records/graders/clinical.md
9b2b65ed6f009c0b96ef606440e346fb87711b1be89b12add7605f99e2bb162e  heldout-emergency-in-records/graders/urgent_first.md
125dba1a89b245f9b1e2c7b46e0b9770deaff9953e489325ac6b3636cff63f49  heldout-somatic-tumour/graders/clinical.md
7fdb90ff794f0bf25451197277425079b078982295bc60b9db37b9b55718fa01  heldout-somatic-tumour/graders/no_hpo.md
5ee15fe36323f40db4dcf377874db0cef35b9dc1ac9bf80c7a7a159d2fadb394  heldout-somatic-tumour/graders/osimertinib.md
474dc0e5850f48ee11e99c0523ef2b9d02a1dcd5fce6614ea7e2c6d974964185  heldout-undiagnosed-en-clinician/graders/clinical.md
809fcbaba0ef5cc142d6fccaa412fced6b5c9f06f9d01b11f4b3d3f20300c081  heldout-undiagnosed-en-clinician/graders/names_wolfram.md
942586d2a31b9356d61d0a7c40ce6b416b8628e4dcc7d49b9379db573784d1f6  heldout-undiagnosed-lay-zh/graders/clinical.md
630551282418569a32bc9154f819940843d7041fc4f197113336073ed6a3f91e  heldout-undiagnosed-lay-zh/graders/names_homocystinuria.md
0c726ccdc9d2b5e4c046017bf8db0c824752344a5e5aa5444e855deb54872e76  heldout-vus-segregation/graders/clinical.md
```
