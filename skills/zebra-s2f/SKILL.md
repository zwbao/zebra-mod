---
name: zebra-s2f
description: Sequence-to-function analysis of a variant for rare disease — splicing (SpliceAI, Pangolin), regulation and expression by tissue (AlphaGenome), zero-shot sequence likelihood (Evo 2), conservation (GPN-MSA), missense function (AlphaMissense, REVEL) — read axis by axis with claim ceilings and combined by agreement, never averaged; then what the mechanism implies (NMD, exon skipping, pseudoexon, therapy angle) and how to confirm it in the lab. Triggers: splice, 剪接, deep intronic, 深内含子, non-coding, 非编码, regulatory, promoter, UTR, AlphaGenome, SpliceAI, Evo2, sequence-to-function, S2F, 序列到功能.
---

# zebra-s2f — predictions as hypotheses, axis by axis

Three rules carried over from s2f-agent:
1. **You own the judgement**: which axes matter for this variant and this disease.
2. **Code owns the numbers**: a number enters only as a tool returned it, with its model and version.
3. **Axes combine by logic, not arithmetic**: no averaging, no voting, no "consensus score". A model that was not run is `not_run`, never "no effect".

## Axes

| Axis | Question | Models (zebra) | Ceiling |
|---|---|---|---|
| splicing | does it create/destroy a splice site; which exon/intron? | SpliceAI, Pangolin (`s2f_predict`); AlphaGenome splice tracks | molecular |
| regulation | does it change expression / accessibility in the relevant tissue? | AlphaGenome (needs key + s2f CLI) | cellular |
| constraint | is the position conserved across vertebrates? | GPN-MSA (s2f CLI), phyloP/GERP when VEP returns them | molecular |
| sequence likelihood | is ALT less likely than REF under a genome language model? | Evo 2 (needs key) | molecular |
| missense function | does the residue change damage the protein? | AlphaMissense, REVEL, CADD (`variant_card`) | molecular |
| evidence | what is observed in people? | ClinVar, gnomAD, literature (`variant_card`, `zebra-literature`) | the only route to clinical claims |

No prediction reaches "causes the disease". The highest a model speaks to is a cellular change.

## Steps

1. **Anchor**: `mcp__zebra-mod__variant_card` → consequence on the MANE transcript, exon/intron number, the intronic offset in the c. notation (c.3718-2477 lies 2,477 nt upstream of the exon that begins at c.3718), gene, assembly.
2. **Choose axes** by variant class: splice region / deep intronic / synonymous → splicing first; promoter / UTR / enhancer → regulation (the disease tissue's ontology: take it from `mcp__zebra-mod__expression`, which lists each GTEx tissue with its UBERON/EFO id and where the gene is expressed — e.g. brain `UBERON:0000955`, liver `UBERON:0002107`, skeletal muscle `UBERON:0001134`, heart `UBERON:0000948`; ask if unclear); missense → missense function + splicing (exonic variants can break splicing too).
3. **Run**: `mcp__zebra-mod__s2f_predict` (`models` as chosen). Several slow axes (Evo 2, AlphaGenome in several tissues) or several variants → one `zebra-mod:s2f-analyst` subagent per variant or axis, all in one message, each returning its agreement rows. Heavy models report `not_run` with the reason when keys or the `s2f` CLI are missing — say so; offer setup (below) instead of guessing.
4. **Read each output** before combining:
   - SpliceAI Δ: ≥ 0.2 supports a splice effect (PP3, Walker 2023), ≥ 0.5 high confidence, ≤ 0.1 supports no effect (BP4/BP7). Read *which* score (acceptor/donor gain/loss) and the position → predict the transcript: exon skipping, intron retention, cryptic exon (pseudoexon), shifted site; is the change a multiple of 3 (in-frame) or does it create a premature stop → NMD?
   - AlphaGenome: a window-wide mean near 0 is not "no effect"; look at the local change in the relevant tissue's tracks.
   - Evo 2 `delta_loglik` (nats, ALT − REF): negative = ALT less likely; a hypothesis about constraint, not pathogenicity.
5. **Agreement matrix** in your answer: axis × {supports, opposes, silent, not_run}. Two independent axes agreeing is stronger than one model with a high score.
6. **Mechanism → implications**: loss of function via NMD; in-frame skip (may be partially functional); pseudoexon inclusion, cryptic acceptor or donor → `mcp__zebra-mod__aso_screen` (the aberrant event with its boundaries checked in the patient sequence, candidate target windows, published N-of-1 precedents); the milasen route continues in `zebra-therapy`. A missense variant: `variant_card` carries MaveDB scores where a calibrated dataset covers it.
7. **How to confirm**: RNA from a tissue that expresses the gene — `mcp__zebra-mod__expression` says whether blood, lymphoblastoid cells, fibroblasts, skin or muscle reach 1 TPM, so the test is asked for in a tissue that can show it; an NMD-targeted transcript may need NMD inhibition (cycloheximide/puromycin) to be seen. RT-PCR or RNA-seq; minigene assay; allele-specific expression. Predictions do not replace this; RNA evidence outranks every model.

## Setting up the heavy models

What each needs (a missing piece gives `not_run` with the reason, never a guess):
- AlphaGenome: `ontology` is required — the tissue the disease acts in; without it zebra does not run it (the upstream script would silently use transverse colon). ~5 s.
- GPN-MSA: GRCh38 SNVs only; needs `tabix` (htslib). ~1 min.
- Evo 2: hosted, 3–10 minutes per variant; ask before running it, or run it in the background via Bash (`zebra s2f predict ... --models evo2 --timeout 900`).
- GRCh37 inputs: SpliceAI/Pangolin only.

- `s2f` CLI from s2f-penguin: `uv tool install "git+https://github.com/zwbao/s2f-penguin"`, set `S2F_DEPLOY_ROOT` (where its model environments live; and `S2F_S2F_CORE_PYTHON` for a standalone install), then `s2f doctor`. zebra finds it on `PATH` or via `S2F_BIN`.
- AlphaGenome: `ALPHAGENOME_API_KEY` (Google DeepMind; non-commercial use only, not for clinical decision-making). Evo 2: `NVCF_RUN_KEY` (NVIDIA). Keys go in the environment, never in a file in the case.
- SpliceAI/Pangolin come from the Broad SpliceAI-lookup service (research use; the SpliceAI weights are CC BY-NC).

## What would make this wrong

A residue number without its transcript; reading `not_run` or an empty track as "no effect"; comparing raw scores across models; using the same axis twice as two pieces of evidence; a tissue the gene is not expressed in.
