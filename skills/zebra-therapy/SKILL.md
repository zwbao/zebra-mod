---
name: zebra-therapy
description: Genotype-to-therapy for a rare disease — start from the disease mechanism, then map approved and investigational drugs, recruiting trials (worldwide and in China), repurposing leads with evidence, and individualized (N-of-1) routes: splice-switching or knockdown antisense oligonucleotides, gene replacement, base or prime editing — each lead tiered by evidence and checked for mechanism fit. Leads to discuss with the care team, never a prescription. Triggers: 治疗, treatment, 有什么药, drug, 孤儿药, orphan drug, clinical trial, 临床试验, repurposing, 老药新用, ASO, 反义寡核苷酸, gene therapy, 基因治疗, base editing, 碱基编辑, N-of-1, 个体化治疗.
---

# zebra-therapy — mechanism → leads → tiers

A lead enters only with retrieved evidence and passes the mechanism check. No dosing, no advice to start or stop anything.

Before any lead: read `zebra-safety` for this disease's drug, anaesthesia and procedure hazards, and raise them as questions for the care team.

## 1. Mechanism first

From the case and `gene_card` / `disease_card` (and `zebra-variant` / `zebra-s2f` when the variant matters):

| Mechanism | Rational directions | Wrong direction |
|---|---|---|
| loss of function — haploinsufficiency | raise the remaining allele (ASO against a poison exon or uORF, CRISPRa), gene replacement | inhibitors of the gene product |
| loss of function — biallelic | replacement (enzyme, AAV gene therapy), substrate reduction, chaperones, read-through for nonsense | — |
| gain of function / toxic | allele-specific or total knockdown (ASO, siRNA), inhibitors | replacement or upregulation |
| dominant negative | allele-specific silencing (± replacement) | replacement alone |
| splicing defect | splice-switching ASO (pseudoexon, cryptic site), exon skipping to restore frame | — |
| repeat expansion | knockdown, repeat-targeting | — |

## 2. Gather leads (all in parallel when possible)

- `mcp__zebra-mod__therapy_landscape` (disease and gene): known drugs with phase and mechanism (Open Targets / ChEMBL), tractability.
- `mcp__zebra-mod__trials_search`: recruiting first; with `country: "China"` for Chinese families as well as worldwide. It returns the eligibility criteria (`full_eligibility: true` for the whole text) and contacts: read the genotype and age criteria against this patient instead of guessing, and repeat a status the tool flagged as stale as "last updated <date>, worth checking". It does not cover ChiCTR (which blocks scripts), so for China add that the Chinese registry should also be searched; trials whose conditions do not name the disease are moved to `filtered` (`keep_unrelated: true` keeps them).
- `mcp__zebra-mod__literature_search`: case reports and series of treatment in this disease/gene (off-label use, compassionate use), reviews.
- Access, worldwide and in China: `mcp__zebra-mod__access` for each drug or for the disease — FDA/EMA status from the agencies' own records (a withdrawal or refusal shows as such), approval in China from the bundled NMPA/CDE documents (`not_in_bundled_list` is not "not approved"), the 2025 national reimbursement list entry with its restriction text verbatim, trials with sites in China, and the collaboration-network hospitals of a province. ChiCTR cannot be queried by a script: give the person its search address instead of saying there are no Chinese trials.
- Large searches → one `zebra-mod:therapy-scout` subagent per angle (approved/investigational, trials, literature/repurposing, N-of-1), all in one message.

## 3. N-of-1 feasibility screens (research)

- **Antisense**: variant class (splice / pseudoexon / GoF allele / haploinsufficiency with a targetable poison exon), tissue reachable (CNS by intrathecal dosing, eye, liver; muscle is harder), gene expressed in it (`mcp__zebra-mod__expression`), onset vs disease stage. For a splice or pseudoexon variant run `mcp__zebra-mod__aso_screen`; for an out-of-frame exon deletion, `aso_screen` with `event: "exon_skip"`. It retrieves the N=1 Collaborative / n-Lorem precedents from Europe PMC in the run — cite those, and if the block is empty say the search failed, never paraphrase from memory.
- **Exon skipping**: for an out-of-frame exon deletion, `mcp__zebra-mod__cnv_interpret` gives the frame arithmetic (coding bases, intronic breakpoints read with their sign) and which added exon restores it; that is the exon the antisense drug must skip. An in-frame deletion returns "already in frame" and no skip target; SMN1/SMN2 is never an exon-skipping target (it is a copy-number finding — nusinersen modulates SMN2 splicing instead). Check an approved or trial drug for exactly that exon, never for the gene in general.
- **Base editing**: `mcp__zebra-mod__edit_check` (SNV revertible by ABE/CBE? protospacer with the base in the window? bystanders?). No NGG protospacer → rerun with `pam: "NG"` (relaxed-PAM Cas9 variants); `annotate_bystanders: true` checks whether bystander edits change the protein. Delivery to the relevant tissue is the hard part — say so.
- **Gene replacement**: coding sequence vs AAV capacity (~4.7 kb including regulatory elements); dosage sensitivity (overexpression toxicity, e.g. MECP2).

## 4. Tier and check every lead

| Tier | Meaning |
|---|---|
| A | approved for this disease (where; source) |
| B | in clinical trials for this disease or genotype |
| C | clinical case reports / off-label in this disease |
| D | preclinical evidence in models of this disease/gene |
| E | mechanistic hypothesis only |

Mechanism fit: does the lead act on the disease gene or its pathway, in the right direction, in a tissue that matters? A lead that fails is dropped or labelled "mechanism mismatch".

## 5. Record and present

`case_update` → `leads` (kind: approved / trial / repurposing / n-of-1 / supportive) with evidence ids. Families: what exists, what is being studied, what is only an idea, and what to ask the doctors — hope stated honestly. Clinicians/researchers: the tiered table with sources.

What would make a lead wrong: a mechanism mismatch; a trial already closed or excluding this genotype; a "drug for this disease" that the source shows was for another indication.
