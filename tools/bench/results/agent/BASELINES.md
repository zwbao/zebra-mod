# Published baselines: rare-disease diagnosis AI tools

Compiled 2026-10-25. Every number is quoted from a primary source saved under `sources/`. Nothing here was computed or estimated. "—" means the source does not report the value. "not retrieved" means I could not get the value from a primary source. All values are percentages.

Abbreviations: R@k = Recall@k (also called top-k or hit@k). "LLM judge" = an LLM decides whether a predicted disease name matches the gold diagnosis. "ID match" = matching ontology IDs (OMIM, ORPHA or MONDO).

Source keys:
- **[DR-SI]**: DeepRare, Nature 2026 (doi 10.1038/s41586-025-10097-9), Supplementary Information PDF (`sources/deeprare/MOESM1_ESM.pdf`), §1.2 **"Table 1 | Rare Disease Diagnosis Result (HPO Input only)"**, p.4. The Nature main text calls this "Supplementary Table 1". In arXiv 2506.20430v3 the same table is "Extended Data Table 2" (§13.2), and the arXiv text calls it "Supplementary Table 2". The numbers are identical in the two versions.
- **[DR-SI-T2]**: same SI PDF, Table 2 (ablation on the RareBench subsets with GPT-4o as host), p.6.
- **[RB]**: RareBench, KDD 2024 (arXiv 2402.06341v2, `sources/rarebench/rarebench_arxiv.pdf`).
- **[HPOQ]**: HPOQuest, arXiv 2609.18431v1 (described as MICCAI 2026), Table 1, p.6.
- **[BioNLP]**: Elmofty & Leser, "When Does Retrieval Beat Direct LLM Diagnosis in Rare Disease?", BioNLP 2026 (aclanthology 2026.bionlp-1.41), Table 2, p.512.
- **[RDG]**: RDguru, IEEE JBHI (doi 10.1109/jbhi.2024.3464555), §III.B.3 / Fig. 6. Taken from the open-access HTML full text (`sources/rdguru/jbhi_exa_capture_excerpts.txt`). The IEEE PDF was blocked to curl.
- **[REESE]**: Reese et al., EJHG 2026; 34:498–504 (doi 10.1038/s41431-026-02054-5), Results "Evaluation on 5213 cases", p.501 (`sources/reese/article.pdf`).
- **[META]**: Nguyen et al., medRxiv 2026 (doi 10.64898/2026.03.26.26349194), Table 1 (p.6), Fig. 2 forest plot, Table 2 (`sources/metaanalysis/full.pdf`).

## How each source judged correctness

| Source | Judge | Ranked list | Input |
|---|---|---|---|
| DeepRare [DR-SI] | **LLM judge.** Main text (Methods): "For each diagnostic task, we generated five most probable diagnostic predictions. The position of the correct diagnosis within these predictions was determined using GPT-4o under Prompt 2 (more detailed reliability validation can be found in Supplementary Methods)". The SI does not contain a judge-reliability section (I searched for "reliab", "agreement" and "kappa"). | top-5 (no R@10) | HPO terms, given to LLMs as term names |
| RareBench [RB] | Paper: "deeming a diagnosis correct if the actual disease appears within the top-k predictions" (§4.3.1). The paper does not name the matching method. Code (`utils/evaluation.py`, `main.py`) uses an **LLM judge** with `--judge_model` default `"chatgpt"` (choices `gpt4`, `chatgpt`). | top-10 | HPO terms (public subsets) |
| HPOQuest [HPOQ] | "GPT-4o-mini with DeepRare's judging prompt … accepting predictions verified by either the LLM judge or a MONDO synonym" (§3) | top-5 | full HPO profile, or 3 seed HPO terms plus active questioning |
| BioNLP 2026 [BioNLP] | "All our metrics use strict OMIM ID matching" (Table 2 caption) | top-5 reported | HPO profile / narrative prompt |
| RDguru [RDG] | ID match: GPT-4 output mapped to ORPHA IDs with "regular expression matching and ORPHA retriever". The ground truth is ORPHA (`summary.csv` has an ORPHA column). | top-10 | HPO terms ("All present phenotypes in each test case") |
| Reese [REESE] | Ontology grounding. Each item was converted to MONDO with OAK (exact label or synonym), and the rest with CurateGPT above a threshold ("94.0% of items … were assigned Mondo terms"). An item counts as correct if its ID "matched the identifier for gold standard diagnosis … exactly, or was mapped as equivalent", or if the gold diagnosis is a close MONDO descendant. Scoring used the PhEval framework. | top-10 | prompts generated from phenopackets by phenopacket2prompt (HPO-derived text) |

**DeepRare Prompt 2, verbatim** ([DR-SI] §1.8, p.8, "Prompt 2. Prompt for Diagnosis Result Evaluation"):

```
You are a specialist in the field of rare diseases.

I will now give you five predicted diseases. Please identify the rank of the following gold-standard
diagnosis.

Please output the predicted rank; otherwise, output "No". Only output "No" or "1-5" numbers. If the
predicted disease has multiple conditions, only output the top rank.
Output only "No" or one number, no additional output.

Predicted diseases: {predict_diagnosis}
Standard diagnosis: {golden_diagnosis}
```

The released code (`github.com/MAGIC-AI4Med/DeepRare` `eval.py`, saved at `sources/deeprare/github/eval.py`) uses slightly different wording. Its system prompt is "You are a specialist in the field of rare diseases." and its user prompt is: `'I will now give you five predicted diseases if the predicted diagnosis is in the standard diagnosis. Please output the predicted rank, otherwise output "No", only output "No" or "1-5" numbers, if the predicted disease has multiple conditions, only output the top rank. Output only "No" or one number, no additional output.'` followed by the predicted and standard diagnoses. The judge model defaults to `--openai_model gpt-4o`. RareBench's judge prompt is the same text with "ten" and "1-10".

For context, DeepRare's baseline-LLM prompt is Prompt 1 ([DR-SI] p.7). It says "Enumerate the top 5 most likely diagnoses … (at least 5). The top 10 diagnoses are:", so the prompt itself is internally inconsistent.

---

## 1. RareBench-MME (n = 40)

| System | R@1 | R@3 | R@5 | R@10 | n | Judge | Source + location |
|---|---|---|---|---|---|---|---|
| DeepRare (gemini-2), best config by R@1 | 77.50 | 85.00 | 90.00 | — | 40 | GPT-4o, Prompt 2 | [DR-SI] Table 1 |
| DeepRare (GPT-4o) | 72.50 | 77.50 | 82.50 | — | 40 | GPT-4o, Prompt 2 | [DR-SI] Table 1 |
| DeepRare (claude-3.7) | 72.50 | 85.00 | 87.50 | — | 40 | " | [DR-SI] Table 1 |
| DeepRare (DS-V3) | 70.00 | 77.50 | 77.50 | — | 40 | " | [DR-SI] Table 1 |
| DeepRare (DS-R1) | 67.57 | 81.08 | 83.78 | — | 40 | " | [DR-SI] Table 1 |
| DeepRare GPT-4o host, case search off ("LLM + Tool-calling + Reflection") | 42.50 | 57.50 | 57.70 | — | 40 | " | [DR-SI-T2] |
| RAG (DS-V3) | 42.50 | 45.00 | 47.50 | — | 40 | " | [DR-SI] Table 1 |
| PubCaseFinder | 40.00 | 57.50 | 57.50 | — | 40 | " | [DR-SI] Table 1 |
| DeepSeek R1 | 35.00 | 60.00 | 65.00 | — | 40 | " | [DR-SI] Table 1 |
| Claude 3.7 think | 27.50 | 50.00 | 57.50 | — | 40 | " | [DR-SI] Table 1 |
| PhenoBrain | 25.00 | 45.00 | 80.00 | — | 40 | " | [DR-SI] Table 1 |
| Claude-3.7-Sonnet | 7.50 | 27.50 | 32.50 | — | 40 | " | [DR-SI] Table 1 |
| o3mini | 2.50 | 7.50 | 7.50 | — | 40 | " | [DR-SI] Table 1 |
| GPT-4o | 2.50 | 5.00 | 10.00 | — | 40 | " | [DR-SI] Table 1 |
| Retrieve+rerank, Gemma-4 31B (best config chosen on the test set) | 70.0 | — | 90.0 | — | 40 | strict OMIM ID | [BioNLP] Table 2 |
| PhenoDP retriever | 65.0 | — | 85.0 | — | 40 | strict OMIM ID | [BioNLP] Table 2 |
| Llama-3.3 70B, unrestricted | 5.0 | — | 10.0 | — | 40 | strict OMIM ID | [BioNLP] Table 2 |
| HPOQuest (Qwen-32B), 3 seed HPO terms | 35.0 | 42.5 | 50.0 | — | 40 | GPT-4o-mini with DeepRare prompt, or MONDO synonym | [HPOQ] Table 1 |
| DeepRare (Qwen-32B), full profile, case retrieval disabled | 17.5 | 22.5 | 35.0 | — | 40 | " | [HPOQ] Table 1 |

## 2. RareBench-HMS (n = 88)

| System | R@1 | R@3 | R@5 | R@10 | n | Judge | Source + location |
|---|---|---|---|---|---|---|---|
| DeepRare (DS-V3), best config by R@1 | 56.97 | 65.12 | 70.93 | — | 88 | GPT-4o, Prompt 2 | [DR-SI] Table 1 |
| DeepRare (claude-3.7) | 56.82 | 63.64 | 64.77 | — | 88 | " | [DR-SI] Table 1 |
| DeepRare (GPT-4o) | 54.54 | 65.90 | 69.31 | — | 88 | " | [DR-SI] Table 1 |
| DeepRare (DS-R1) | 51.76 | 60.00 | 67.06 | — | 88 | " | [DR-SI] Table 1 |
| DeepRare (gemini-2) | 40.91 | 54.25 | 57.95 | — | 88 | " | [DR-SI] Table 1 |
| GPT-4o (best non-DeepRare baseline on HMS) | 47.73 | 60.23 | 65.91 | — | 88 | " | [DR-SI] Table 1 |
| Claude 3.7 think | 44.32 | 57.95 | 62.50 | — | 88 | " | [DR-SI] Table 1 |
| Claude-3.7-Sonnet | 43.18 | 60.23 | 63.64 | — | 88 | " | [DR-SI] Table 1 |
| DeepSeek R1 | 38.64 | 53.41 | 60.23 | — | 88 | " | [DR-SI] Table 1 |
| o3mini | 30.68 | 44.32 | 54.55 | — | 88 | " | [DR-SI] Table 1 |
| PhenoBrain | 22.35 | 32.94 | 42.35 | — | 88 | " | [DR-SI] Table 1 |
| PubCaseFinder | 7.32 | 13.41 | 15.85 | — | 88 | " | [DR-SI] Table 1 |
| Med42 70B, unrestricted | 20.6 | — | 27.9 | — | **68** | strict OMIM ID | [BioNLP] Table 2 (OMIM-evaluable subset) |
| Retrieve+rerank, Llama-3.3 70B | 16.2 | — | 36.8 | — | 68 | " | [BioNLP] Table 2 |
| PhenoDP retriever | 7.4 | — | 25.0 | — | 68 | " | [BioNLP] Table 2 |
| HPOQuest (Qwen-32B), 3 seeds | 15.9 | 23.9 | 25.0 | — | 88 | GPT-4o-mini judge or MONDO synonym | [HPOQ] Table 1 |

## 3. RareBench-LIRICAL (n = 370)

| System | R@1 | R@3 | R@5 | R@10 | n | Judge | Source + location |
|---|---|---|---|---|---|---|---|
| DeepRare (DS-V3), best config by R@1 | 55.95 | 65.16 | 67.57 | — | 370 | GPT-4o, Prompt 2 | [DR-SI] Table 1 |
| DeepRare (DS-R1) | 52.43 | 63.51 | 67.84 | — | 370 | " | [DR-SI] Table 1 |
| DeepRare (GPT-4o) | 51.62 | 62.70 | 65.40 | — | 370 | " | [DR-SI] Table 1 |
| DeepRare (gemini-2) | 48.91 | 63.24 | 68.91 | — | 370 | " | [DR-SI] Table 1 |
| DeepRare (claude-3.7) | 40.00 | 51.08 | 54.86 | — | 370 | " | [DR-SI] Table 1 |
| Claude 3.7 think (best non-DeepRare baseline) | 38.11 | 50.27 | 52.43 | — | 370 | " | [DR-SI] Table 1 |
| PubCaseFinder | 37.77 | 47.01 | 48.91 | — | 370 | " | [DR-SI] Table 1 |
| DeepSeek R1 | 36.76 | 45.14 | 47.57 | — | 370 | " | [DR-SI] Table 1 |
| Claude-3.7-Sonnet | 33.24 | 46.49 | 50.00 | — | 370 | " | [DR-SI] Table 1 |
| GPT-4o | 31.08 | 42.16 | 48.92 | — | 370 | " | [DR-SI] Table 1 |
| PhenoBrain | 27.10 | 44.72 | 57.72 | — | 370 | " | [DR-SI] Table 1 |
| o3mini | 25.68 | 35.41 | 40.00 | — | 370 | " | [DR-SI] Table 1 |
| Retrieve+rerank, Gemma-4 31B / Llama-3.3 70B | 54.1 | — | 75.9 | — | 370 | strict OMIM ID | [BioNLP] Table 2 |
| PhenoDP retriever | 52.7 | — | 70.8 | — | 370 | " | [BioNLP] Table 2 |
| HPOQuest (Qwen-32B), 3 seeds | 28.6 | 40.3 | 45.1 | — | 370 | GPT-4o-mini judge or MONDO synonym | [HPOQ] Table 1 |

## 4. RareBench-RAMEDIS (n = 624)

| System | R@1 | R@3 | R@5 | R@10 | n | Judge | Source + location |
|---|---|---|---|---|---|---|---|
| DeepRare (GPT-4o), best config by R@1 | 73.23 | 81.57 | 85.26 | — | 624 | GPT-4o, Prompt 2 | [DR-SI] Table 1 |
| DeepRare (DS-V3) | 72.60 | 82.85 | 84.62 | — | 624 | " | [DR-SI] Table 1 |
| DeepRare (claude-3.7) | 70.03 | 81.57 | 84.46 | — | 624 | " | [DR-SI] Table 1 |
| DeepRare (DS-R1) | 62.17 | 76.28 | 79.32 | — | 624 | " | [DR-SI] Table 1 |
| DeepRare (gemini-2) | 49.67 | 69.55 | 75.64 | — | 624 | " | [DR-SI] Table 1 |
| Claude-3.7-Sonnet (best non-DeepRare baseline) | 42.47 | 68.59 | 73.88 | — | 624 | " | [DR-SI] Table 1 |
| Claude 3.7 think | 41.83 | 63.78 | 72.92 | — | 624 | " | [DR-SI] Table 1 |
| o3mini | 39.42 | 64.42 | 67.47 | — | 624 | " | [DR-SI] Table 1 |
| GPT-4o | 35.47 | 52.24 | 61.06 | — | 624 | " | [DR-SI] Table 1 |
| DeepSeek R1 | 33.17 | 47.60 | 53.21 | — | 624 | " | [DR-SI] Table 1 |
| PubCaseFinder | 25.84 | 32.58 | 39.97 | — | 624 | " | [DR-SI] Table 1 |
| PhenoBrain | 20.67 | 38.30 | 55.77 | — | 624 | " | [DR-SI] Table 1 |
| Med42 70B, unrestricted | 46.8 | — | 64.1 | — | 624 | strict OMIM ID | [BioNLP] Table 2 |
| Llama-3.3 70B, unrestricted | 46.0 | — | 66.7 | — | 624 | " | [BioNLP] Table 2 |
| PhenoDP retriever | 10.3 | — | 28.7 | — | 624 | " | [BioNLP] Table 2 |
| HPOQuest (Qwen-32B), 3 seeds | 27.6 | 45.4 | 49.0 | — | 624 | GPT-4o-mini judge or MONDO synonym | [HPOQ] Table 1 |
| DxGPT (Claude 3 Opus) | 55.0 | not retrieved | not retrieved | not retrieved | **200** | not retrieved | [META] Table 1 only (do Olmo 2024); not checked against the primary paper |

**RareBench paper, per subset:** RareBench itself reports Task 4 only **pooled over 2,185 cases**. That pool is the public subsets plus private PUMCH cases. The paper has no per-subset numbers for MME, HMS, LIRICAL or RAMEDIS. Best pooled result in [RB] Table 4, Task 4: **GPT-4 CoT R@1 33.2, R@3 46.4, R@10 59.5, MR 4.0**. GPT-4 0-shot: 32.3 / 45.4 / 58.9, MR 5.0. GPT-3.5-Turbo: 21.1 / 34.2 / 48.2. Gemini Pro: 14.6 / 22.7 / 33.0. Other [RB] numbers come from the private PUMCH-75 cases: Table 5 (dynamic 3-shot GLM3-Turbo hit@1 72.0) and Table 6 (GPT-4 with phenotypes hit@1 52.0). Those cases are not public, apart from the 75 PUMCH_ADM cases on HuggingFace.

**RareBench sizes, confirmed:** the HF `data.zip` jsonl files have line counts MME 40, LIRICAL 370, HMS 88, RAMEDIS 624, which matches the HF README table. The README gives disease counts of 17, 39, 252 and 74 for MME, HMS, LIRICAL and RAMEDIS. Each record has `Phenotype` (HPO IDs) and `RareDisease` (OMIM / Orphanet / CCRD codes).

## 5. MyGene2 (n = 146 in DeepRare)

| System | R@1 | R@3 | R@5 | R@10 | n | Judge | Source + location |
|---|---|---|---|---|---|---|---|
| DeepRare (GPT-4o), best config by R@1 | 75.86 | 80.00 | 80.69 | — | 146 | GPT-4o, Prompt 2 | [DR-SI] Table 1 |
| DeepRare (DS-R1) | 75.00 | 80.55 | 81.94 | — | 146 | " | [DR-SI] Table 1 |
| DeepRare (claude-3.7) | 73.79 | 80.68 | 80.69 | — | 146 | " | [DR-SI] Table 1 |
| DeepRare (DS-V3) | 72.60 | 76.03 | 76.03 | — | 146 | " | [DR-SI] Table 1 |
| DeepRare (gemini-2) | 69.18 | 77.40 | 79.45 | — | 146 | " | [DR-SI] Table 1 |
| PhenoBrain (best non-DeepRare baseline) | 45.21 | 58.90 | 61.64 | — | 146 | " | [DR-SI] Table 1 |
| Gemini 2 think | 39.72 | 53.42 | 53.42 | — | 146 | " | [DR-SI] Table 1 |
| DeepSeek R1 | 35.62 | 46.58 | 47.95 | — | 146 | " | [DR-SI] Table 1 |
| o3mini | 34.48 | 46.21 | 49.66 | — | 146 | " | [DR-SI] Table 1 |
| Claude 3.7 think | 34.25 | 37.67 | 39.72 | — | 146 | " | [DR-SI] Table 1 |
| PubCaseFinder | 20.55 | 26.71 | 28.08 | — | 146 | " | [DR-SI] Table 1 |
| Claude-3.7-Sonnet | 19.44 | 36.11 | 39.58 | — | 146 | " | [DR-SI] Table 1 |
| GPT-4o | 19.31 | 35.17 | 37.93 | — | 146 | " | [DR-SI] Table 1 |
| Retrieve+rerank, Gemma-4 31B | 48.9 | — | 61.8 | — | **131** | strict OMIM ID | [BioNLP] Table 2 (PhenoDP-compatible subset) |
| PhenoDP retriever | 38.9 | — | 61.8 | — | 131 | " | [BioNLP] Table 2 |
| Llama-3.3 70B, unrestricted | 3.8 | — | 15.3 | — | 131 | " | [BioNLP] Table 2 |

## 6. DDD (n = 2,283 in DeepRare; disease-level G2P profiles, not patients)

| System | R@1 | R@3 | R@5 | R@10 | n | Judge | Source + location |
|---|---|---|---|---|---|---|---|
| DeepRare (gemini-2), best config by R@1 | 47.59 | 59.82 | 62.70 | — | 2283 | GPT-4o, Prompt 2 | [DR-SI] Table 1 |
| DeepRare (claude-3.7) | 45.49 | 58.85 | 63.76 | — | 2283 | " | [DR-SI] Table 1 |
| DeepRare (DS-R1) | 44.27 | 57.09 | 61.27 | — | 2283 | " | [DR-SI] Table 1 |
| DeepRare (DS-V3) | 42.97 | 53.48 | 56.72 | — | 2283 | " | [DR-SI] Table 1 |
| DeepRare (GPT-4o) | 39.85 | 51.36 | 55.88 | — | 2283 | " | [DR-SI] Table 1 |
| RAG (DS-V3) (best non-DeepRare baseline) | 39.42 | 48.84 | 51.77 | — | 2283 | " | [DR-SI] Table 1 |
| PhenoBrain | 38.92 | 54.68 | 63.11 | — | 2283 | " | [DR-SI] Table 1 |
| DeepSeek R1 | 38.24 | 49.43 | 53.01 | — | 2283 | " | [DR-SI] Table 1 |
| o3mini | 34.09 | 42.96 | 45.50 | — | 2283 | " | [DR-SI] Table 1 |
| Claude 3.7 think | 34.01 | 45.27 | 48.83 | — | 2283 | " | [DR-SI] Table 1 |
| GPT-4o | 33.30 | 42.28 | 45.18 | — | 2283 | " | [DR-SI] Table 1 |
| Claude-3.7-Sonnet | 32.28 | 42.76 | 45.93 | — | 2283 | " | [DR-SI] Table 1 |
| PubCaseFinder | 26.92 | 32.32 | 34.05 | — | 2283 | " | [DR-SI] Table 1 |
| PhenoDP retriever / reranker | 59.0 | — | 70.2 | — | **2248** | strict OMIM ID | [BioNLP] Table 2 |

## 7. MIMIC-IV-Rare (n = 1,875; credentialed access)

| System | R@1 | R@3 | R@5 | R@10 | n | Judge | Source + location |
|---|---|---|---|---|---|---|---|
| DeepRare (DS-V3), best config by R@1 | 29.19 | 36.60 | 39.06 | — | 1875 | GPT-4o, Prompt 2 | [DR-SI] Table 1 |
| DeepRare (claude-3.7) | 24.32 | 30.30 | 32.12 | — | 1875 | " | [DR-SI] Table 1 |
| DeepRare (DS-R1) | 23.27 | 31.16 | 34.79 | — | 1875 | " | [DR-SI] Table 1 |
| DeepRare (GPT-4o) | 21.44 | 28.89 | 31.71 | — | 1875 | " | [DR-SI] Table 1 |
| DS-R1-search (best non-DeepRare baseline) | 19.95 | 28.75 | 32.69 | — | 1875 | " | [DR-SI] Table 1 |
| Gemini 2 think | 15.47 | 20.53 | 25.39 | — | 1875 | " | [DR-SI] Table 1 |
| DeepSeek R1 | 13.81 | 19.92 | 22.93 | — | 1875 | " | [DR-SI] Table 1 |
| GPT-4o | 13.76 | 22.35 | 25.60 | — | 1875 | " | [DR-SI] Table 1 |
| Claude 3.7 think | 13.70 | 21.28 | 24.32 | — | 1875 | " | [DR-SI] Table 1 |
| Claude-3.7-Sonnet | 12.11 | 18.24 | 20.80 | — | 1875 | " | [DR-SI] Table 1 |
| o3mini | 10.88 | 19.52 | 23.90 | — | 1875 | " | [DR-SI] Table 1 |
| PhenoBrain | 8.11 | 10.24 | 11.41 | — | 1875 | " | [DR-SI] Table 1 |
| PubCaseFinder | 1.06 | 1.39 | 1.60 | — | 1875 | " | [DR-SI] Table 1 |
| MIMIC-RD Llama-3.3 70B | 19.3 | not retrieved | not retrieved | not retrieved | **145** | not retrieved | [META] Table 1 (Al Din 2026); different MIMIC subset; not checked against the primary paper |

## 8. Xinhua Hospital (n = 975; in-house, controlled access)

Only local models were run, "because of privacy considerations".

| System | R@1 | R@3 | R@5 | R@10 | n | Judge | Source |
|---|---|---|---|---|---|---|---|
| DeepRare (DS-V3) | 58.27 | 71.04 | 74.35 | — | 975 | GPT-4o, Prompt 2 | [DR-SI] Table 1 |
| DeepRare (DS-R1) | 54.30 | 66.90 | 70.20 | — | 975 | " | [DR-SI] Table 1 |
| Baichuan-14B | 43.35 | 54.40 | 58.89 | — | 975 | " | [DR-SI] Table 1 |
| DeepSeek V3 | 42.43 | 54.91 | 59.41 | — | 975 | " | [DR-SI] Table 1 |
| DeepSeek R1 | 39.16 | 49.85 | 54.98 | — | 975 | " | [DR-SI] Table 1 |
| MMedS-Llama 3 | 29.55 | 40.18 | 45.30 | — | 975 | " | [DR-SI] Table 1 |

DeepRare's 7-dataset average (excluding Xinhua), from [DR-SI] Table 1 "Avg." column: DeepRare (DS-V3) **57.18 / 65.25 / 67.49**. Best baseline averages: Claude 3.7 think 33.39 / 46.60 / 51.17, DeepSeek R1 33.03 / 46.01 / 49.99, RAG (DS-V3) 31.31 / 40.62 / 43.56. These averages are unweighted means over datasets of very different sizes.

## 9. Phenopacket-store

| System | R@1 | R@3 | R@5 | R@10 | n | Judge | Source + location |
|---|---|---|---|---|---|---|---|
| Exomiser (phenotype-only) | 35.5 | 46.3 | — | 58.5 | 5213 (v0.1.14) | MONDO grounding via OAK/CurateGPT; exact or equivalent ID | [REESE] p.501 Results |
| o1-preview (best LLM) | 23.6 | 31.2 | — | 36.8 | 5213 | " | [REESE] p.501 |
| Meditron3-70B (best medical LLM) | 19 | 28 | — | 31 | 5213 | " | [REESE] p.501 ("19%, 28%, and 31%") |
| Medfound-176B | 0.1 | 0.23 | — | 0.25 | 5213 | " | [REESE] p.501 |
| GPT-4o, o1-mini, Gemini 2.0 Flash, Meditron-70B | not retrieved (shown only as bars in Fig. 3; no numbers in text) | | | | 5213 | " | [REESE] Fig. 3. Per-case data in the Zenodo `correct_results.tsv` (not re-scored here) |
| PhenoDP retriever | 51.2 | — | 67.2 | — | **6901** ("Phenopackets (Wen et al., 2025)") | strict OMIM ID | [BioNLP] Table 2 |
| Retrieve+rerank, Gemma-4 31B | 54.3 | — | 67.6 | — | 6901 | " | [BioNLP] Table 2 |
| Llama-3.3 70B, unrestricted | 14.3 | — | 25.6 | — | 6901 | " | [BioNLP] Table 2 |
| GPT-4o (English) | 19.9 | not retrieved | not retrieved | not retrieved | 4917 | not retrieved | [META] Table 1 (Chimirri 2025); not checked against the primary paper |
| RareSeek-R1 + GraphRAG (phenotype only) | 41.4 (0.414 ± 0.014) | not retrieved | not retrieved | not retrieved | 5213 (version not stated) | Orphanet code match (+ parent terms) | RareSeek arXiv 2511.14638v1, text near Fig. 4b; matches [META] Table 2 |
| Exomiser (in RareSeek) | 27.0 (0.270 ± 0.011) | — | — | — | 5213 | " | same. Compare Reese's Exomiser 35.5: a different judge/pipeline on the same n |

## 10. RDguru 238 published cases (literature-derived; public)

| System | R@1 | R@3 | R@5 | R@10 | n | Judge | Source + location |
|---|---|---|---|---|---|---|---|
| RDguru MixDiagDQN | 41.60 (99/238) | — | 63.87 (152/238) | 69.33 (165/238) | 238 | ORPHA ID (GPT-4 names mapped by regex + ORPHA retriever); candidate set 4,257 ORPHA RDs | [RDG] §III.B.3, Fig. 6 |
| PheLR | 41.60 | — | 58.40 | 68.49 | 238 | " | [RDG] §III.B.3 |
| Phenotype matching (Orphadata frequency) | 25.63 | — | 44.96 | 52.52 | 238 | " | [RDG] §III.B.3 |
| GPT-4 | 26.05 | — | 42.02 | 45.80 | 238 | " | [RDG] §III.B.3 |

Caveat: **the model was selected on the test set.** Algorithm 1, steps 22–23: "Repeat steps 9–15 to test all testing cases using Q(s,a|θ) with ε-greedy set to 0. Save the best model". §II.D.3: "followed by the execution of the testing process, and the parameters of the best-performing model are saved". Training used 10,000 simulated cases whose phenotype counts were matched to the 238 test cases. Tables I and II of [RDG] did not render in the HTML capture.

## 11. Meta-analysis index (medRxiv 2026, 19 system–dataset R@1 entries)

Reproduced from [META] Table 1 (p.6) and the Fig. 2 forest plot. [META] marks some values as "Approximate value (extracted from figure or calculated from reported averages)" (†).

| System (study) | Benchmark | N | R@1 (META) | Verified in primary source? |
|---|---|---|---|---|
| MEDDxAgent GPT-4o few-shot (Rose 2025) | RareBench | 100 | 82% (forest 0.82) | not verified here |
| MEDDxAgent GPT-4o zero-shot (Rose 2025) | RareBench | 100 | 46.0% | not verified here |
| RareSeek-R1 (Yang 2025b) | EHR-Internal | 4,306 | 68.4% | not checked here (private data). RareSeek on public RareBench-1,197: 0.392 (§12) |
| DeepRare (Zhao 2026), DS-V3 | 7 public benchmarks avg | 5,426 | 57.2%† | **verified**: [DR-SI] Avg 57.18 |
| GPT-4o standalone (Zhao 2026) | 7-bench avg | 5,426 | 25.6%† | **mismatch**: [DR-SI] GPT-4o Avg R@1 = **26.16**. 25.60 is GPT-4o's **R@5** on MIMIC-IV-Rare. The Nature main text does say "average Recall@1 … improves substantially from 25.60% to 54.67%" (5 public datasets), so [META] took its value from that sentence |
| RareAgents (Chen 2024a/2025) | RareBench | 1,197 | 55.9% (baseline 41.7%) | **verified**: RareAgents Table 3, 0.5589; GPT-4o 0.4169 |
| Deep-DxSearch (Zheng 2025) | 3 rare datasets avg | 1,426 | 46.1%† | not verified here |
| MINT (Wang 2026) | GMDB test | 386 | 43.0%† | not verified here |
| MAC (Chen 2025) | MAC Medline | 302 | 34.1% | not verified here |
| RareScale augmented / baseline (Schumacher 2025) | Simulated conversations | 3,403 | 33.1% / 28.6% | not verified here |
| GPT-4 (RareBench, Chen 2024b) | RareBench | 2,185 | 32.3% | **verified**: [RB] Table 4, Task 4, GPT-4 0-shot k=1 = 32.3 (pooled, includes private PUMCH) |
| Llama-3.3 (Al Din 2026) | MIMIC-IV | 145 | 19.3% | not verified here |
| GPT-4o (Chimirri 2025) | Phenopacket Store | 4,917 | 19.9% | not verified here |
| GPT-4-turbo (Ilić 2025) | Skeletal dysplasia | 45 | 55.6% | not verified here |
| Claude 3.5 Sonnet (Ao 2025) | Chinese RD List | 152 | 64.5% | not verified here |
| o1-preview (Reese 2026) | Phenopacket Store | 5,213 | 23.6% | **verified**: [REESE] p.501 |
| DxGPT Claude 3 Opus (do Olmo 2024) | RAMEDIS | 200 | 55.0% | not verified here |

[META] Table 2 also lists per-benchmark augmented bests, all attributed to DeepRare: RAMEDIS 72.6% (matches [DR-SI] DS-V3 72.60), LIRICAL 56.0% (matches DS-V3 55.95), RareBench agg. 65.8%.

## 12. RareBench pooled public sets (RareAgents, RareSeek, PhenoBrain)

None of these three papers reports MME, HMS, LIRICAL or RAMEDIS separately. All report pooled sets, so they **cannot be compared directly** with the per-subset tables above.

| System | R@1 | R@3 | R@5 | R@10 | n (pool) | Judge | Source + location |
|---|---|---|---|---|---|---|---|
| RareAgents (Llama-3.1-70B) | 0.5589 | 0.6867 | — | 0.7811 (MR 1.0) | 1,197 "RareBench-Public" (5 subsets incl. PUMCH_ADM 75) | GPT-4o, top-10 rank prompt ("I will now give you ten predicted diseases … '1-10'") | RareAgents arXiv 2412.12475v3, Table 3 |
| RareAgents (Llama-3.1-8B) | 0.4511 | 0.5647 | — | 0.7377 | 1,197 | " | same, Table 3 |
| GPT-4o zero-shot CoT (in RareAgents) | 0.4169 | 0.5815 | — | 0.7068 | 1,197 | " | same, Table 3 |
| MDAgents Llama-3.1-70B (in RareAgents) | 0.4042 | 0.5640 | — | 0.6586 | 1,197 | " | same, Table 3 |
| PhenoBrain (in RareAgents) | 0.2857 | 0.4670 | — | 0.6341 | 1,197 | " | same, Table 4 |
| RareSeek-R1 | 0.392 ± 0.025 | not retrieved (Fig. 2d only) | not retrieved | 0.577 ± 0.027 | 1,197 (MME + LIRICAL + HMS + RAMEDIS + PUMCH_ADM) | Orphanet code match via MONDO API; parent terms also accepted; up to 20 predictions | RareSeek arXiv 2511.14638v1, Results text (Fig. 2d) |
| GPT-5 (in RareSeek) | 0.353 ± 0.049 | not retrieved | not retrieved | not retrieved | 1,197 | " | same |
| Exomiser (in RareSeek) | 0.232 ± 0.016 | not retrieved | not retrieved | 0.575 ± 0.022 | 1,197 | " | same |
| PhenoBrain ensemble (own paper) | 0.304 (0.273–0.334) | 0.483 | — | 0.640 | 873 "Public Test Set" (60% of RAMEDIS = 375, plus MME, LIRICAL, HMS) | exact disease in top-k | Mao et al., npj Digit Med 2025 (doi 10.1038/s41746-025-01452-1), Table 2 |
| Best of 12 baselines (PhenoBrain paper) | 0.257 | 0.424 | — | 0.593 | 873 | " | same, Table 2 |

**RareSeek on phenopacket-store** (n = 5,213, 378 diseases; version not stated). Source: arXiv 2511.14638v1, text near Fig. 4b.
- Phenotype only, RareSeek-R1 + GraphRAG: **R@1 0.414 ± 0.014**.
- Phenotype only, RareSeek-R1 without GraphRAG: 0.210 ± 0.007.
- Phenotype only, Exomiser: 0.270 ± 0.011.
- Variant only, RareSeek-R1 + GraphRAG: 0.615 ± 0.019.
- Phenotype + variant, RareSeek-R1 + GraphRAG: R@1 0.734 ± 0.010, R@10 0.910 ± 0.007.

Note: the paper also says "Because raw VCFs were unavailable…".

Availability:
- PhenoBrain test cases: https://zenodo.org/records/10774650. Code: https://github.com/xiaohaomao/timgroup_disease_diagnosis.
- RareSeek: no code repository; the EHR data are not public.
- RareAgents: the PDF refers to "the code repository" but gives no URL.

---

## Data availability

| Dataset | Where | Notes |
|---|---|---|
| RareBench MME/HMS/LIRICAL/RAMEDIS (+ PUMCH_ADM 75) | https://huggingface.co/datasets/chenxz/RareBench (`data.zip` → `data/*.jsonl`; mapping files in `mapping/`) · code https://github.com/chenxz1111/RareBench | Apache-2.0. HPO IDs plus OMIM/Orphanet/CCRD truth. Downloaded to `sources/rarebench/` |
| **MyGene2, 146 cases** | **Harvard Dataverse doi:10.7910/DVN/TZTPFL**, file **`mygene2_5.7.22.txt`** (datafile id 6689035): https://dataverse.harvard.edu/api/access/datafile/6689035 | Cited by DeepRare's Data availability statement. The dataset is the SHEPHERD release (Alsentzer, Li et al.), CC0, not restricted. My download (`sources/mygene2_dataverse/`) has **146 JSON lines**. Fields: `positive_phenotypes` (HPO IDs), `true_genes` (Ensembl), `true_diseases` (internal/MONDO idx), `omim`, `disease_name`, `orpha_id`, `orpha_category`. So it is **HPO-coded with OMIM and ORPHA truth**. DeepRare (arXiv v3, Datasets) says it used "preprocessed data (146 patients spanning 55 MONDO diseases) from [42 = Alsentzer et al.]". The file holds 44 distinct `omim` values (my count). DeepRare's Extended Data Table 1 says 58 diseases. DeepRare's own processed `dataset/mygene_test.csv` is **not** in its GitHub repo or HF repo, so I could not confirm that the 146 records are identical |
| DeepRare code / knowledge base | https://github.com/MAGIC-AI4Med/DeepRare · https://huggingface.co/datasets/Angelakeke/DeepRare | The HF repo holds only the KB/mapping/embedding files. **No test sets** (code expects `dataset/{mygene,ddd,mimic,xinhua}_test.csv`, which are not released) |
| DDD | DeepRare cites https://www.deciphergenomics.org/ddd/ddgenes. arXiv v3 says it was built from the G2P download https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/2025_05_28/ | Disease-level profiles, not patients. The exact 2,283 split is not released ([BioNLP]: "the original split is not publicly available") |
| MIMIC-IV-Rare | https://physionet.org/content/mimiciv/3.1/ (credentialed) | The 1,875-case derived set is not released |
| Xinhua / Hunan | NGDC PRJCA052720; GVM001237/GVM001238; OMIX013512; GSA-Human HRA015470 | Controlled access via the Xinhua DAC |
| RDguru 238 cases | https://github.com/1678492347/rdguru-test (PMIDs, proband IDs, extracted HPO, text reports; `summary.csv` has 238 rows with an ORPHA column) · MixDiagDQN code https://github.com/1678492347/mix-diag-dqn | Saved under `sources/rdguru/github_rdguru-test/` |
| Phenopacket-store / Reese | Zenodo https://zenodo.org/records/15324355 (`phenopackets.tar.gz`, `prompts.tar.gz`, `correct_results.tsv`, `all_models_responses.zip`), CC-BY · https://github.com/monarch-initiative/phenopacket2prompt · https://github.com/monarch-initiative/pheval.llm | Saved under `sources/reese/zenodo_15324355/` (responses zip not downloaded) |
| BioNLP 2026 | https://github.com/mofty8/BioNLP26 ("Code, prompts, and per-case outputs") | |
| HPOQuest | no code URL found in arXiv v1 text | |

---

## Discrepancies / not retrieved

1. **MME R@1, 78% vs 70%: both figures are in the Nature version.** Main text: "in the RareBench-MME evaluation, DeepRare achieves exceptional scores of 78%, and 85% for Recall@1 and Recall@3". That matches the **gemini-2** row (77.50/85.00). The Fig. 2a caption says "RareBench-MME (70.0% top-1 accuracy)", which matches the **DS-V3** row (70.00). arXiv v3 has both sentences too (lines for "78%, and 85%" and "70.0% top-1"). So the discrepancy is between configurations inside one paper, not between the arXiv and Nature versions. The SI table values are identical in the arXiv and Nature versions.
2. **MyGene2 headline:** main text "74% and 81%". No SI row gives exactly 74/81. The closest is claude-3.7 at 73.79/80.68. The best R@1 row is GPT-4o at 75.86/80.00.
3. **Table numbering:** the full table is "Supplementary Table 1" (Nature main text and SI) but "Supplementary Table 2" or "Extended Data Table 2" in arXiv v3. In Nature, SI Tables 2 and 3 are the ablations.
4. **Disease count:** Nature says 2,919 diseases. arXiv v3 main text says "3,134 rare diseases". RareBench public cases: RareBench paper says "1,122 cases spanning over 362 rare diseases" (40+88+370+624 = 1,122). DeepRare arXiv says "public subset comprising 1,114 patient cases".
5. **Multimodal headline:** Nature reports "69.1% … Exomiser's 55.9% on 168 cases". The DeepRare GitHub README says "70.60% at Recall@1 compared to Exomiser's 53.20% in 109 cases".
6. **MDAgents (DS-V3) MME R@5** is printed as "0.025" in [DR-SI] Table 1. This is probably a typo, but it is reproduced as printed in the source file and omitted from the tables above.
7. **Judge-reliability validation:** the Nature text points to "Supplementary Methods" for it. I found no such validation in the SI PDF.
8. **BioNLP's copy of DeepRare numbers:** [BioNLP] Table 2 lists DeepRare GPT-4o LIRICAL R@1 as **56.0**, but [DR-SI] gives **51.62** (56.0 matches DS-V3 55.95). It also lists DeepRare GPT-4o DDD as **43.0**, where [DR-SI] gives **39.85** (43.0 matches DS-V3 42.97). Its other DeepRare cells match the SI after rounding.
9. **HPOQuest's "DeepRare (GPT-4o)" row** (MME 42.5/57.5/57.7; HMS 50.0/62.5/64.8; LIRICAL 39.5/49.2/52.4; RAMEDIS 36.9/46.4/55.4) equals [DR-SI-T2] "LLM + Tool-calling + Reflection", i.e. DeepRare with case search off. It is not DeepRare's headline result. HPOQuest's judge is GPT-4o-mini, not GPT-4o.
10. **Subset sizes differ between papers:** BioNLP uses HMS n=68 (of 88), MyGene2 n=131 (of 146), DDD n=2,248 (of 2,283), Phenopackets n=6,901. It also selects reranker configs on each test set ("upper-bound"). [META] lists RareAgents on RareBench n=1,197, DxGPT on RAMEDIS n=200, and MIMIC n=145.
11. **R@10 is not available** for DeepRare or any of its baselines, because the evaluation is top-5 only.
12. **RareBench per-subset numbers are not published** in RareBench itself (pooled 2,185 only). RareBench's judge model is not stated in the paper. The code defaults to `chatgpt`.
13. **Reese:** GPT-4o, o1-mini, Gemini 2.0 Flash and Meditron-70B values are only in Fig. 3 bars. Not retrieved as text.
14. **RDguru:** Table I (case summary) was not retrieved, because the IEEE PDF was blocked and Table I did not render in the HTML. No top-3 is reported. Model selection used the test set.
15. **RareSeek, RareAgents and PhenoBrain** report only pooled RareBench sets (1,197 or 873 cases), with no per-subset numbers. RareSeek's top-3/5 values appear only in figures and were not retrieved. It is unclear whether RareSeek's headline numbers use exact or parent-term Orphanet matching.
16. **Exomiser on the same 5,213 phenopackets:** 35.5 R@1 in Reese, 27.0 in RareSeek. Pipelines and judges differ.
17. **HPOQuest venue:** I could not confirm MICCAI 2026. The arXiv v1 has no venue comment, and the papers.miccai.org 2026 listing has no match.
18. **Not verified against primary papers:** 11 of the 19 [META] entries (§11): MEDDxAgent, Deep-DxSearch, MINT, MAC, RareScale, Llama-3.3 MIMIC, Chimirri, Ilić, Ao, DxGPT, and RareSeek-EHR (private data). Their PDFs are partly saved under `sources/others/` but were not read for this table.
