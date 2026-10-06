# Phenotype-ranking benchmark

How often does `zebra phenotype rank` put the right diagnosis, and the right gene, near the top?
Measured on the GA4GH phenopacket-store, with the numbers, their uncertainty, what they do not show,
and how to reproduce them. Every number below is in `tools/bench/results/benchmark.json` (written by
`tools/bench/summary.py`); per-case positions for every case are in
`tools/bench/results/local_cases.jsonl.gz`.

## Headline

The offline `local` ranker that `zebra phenotype rank` uses by default in 0.2, on the 5,195 held-out
cases (363 diseases never used for any choice), compared with zebra 0.1.0:

| correct disease in the top … | top-1 | top-3 | top-10 | top-25 |
|---|---|---|---|---|
| zebra 0.2, held-out | **49.4%** (39.1–62.2) | 59.7% (48.3–73.0) | **69.6%** (57.6–82.2) | 76.5% (66.3–87.0) |
| zebra 0.1.0, held-out | 24.9% (19.5–32.3) | 31.9% (26.0–40.0) | 39.5% (33.5–47.5) | 46.6% (41.5–53.8) |
| zebra 0.2, held-out, **paper not an annotation source** (1,785 cases) | 6.8% (2.6–16.9) | 14.7% (8.1–27.3) | **28.3%** (16.0–43.9) | 41.6% (29.1–56.2) |
| zebra 0.1.0, same cases | 3.6% (1.8–7.2) | 7.7% (4.6–12.0) | 14.6% (9.0–20.8) | 24.9% (15.7–31.8) |

Causal gene in the top 10, held-out: 72.0% (61.7–82.1) in 0.2, 44.2% in 0.1.0; on the cases whose paper
was not an annotation source: 39.6% and 26.5%.

95% intervals from a bootstrap that resamples whole diseases (cases of one disease, often one family in one
paper, are not independent). The truth is an OMIM id and only OMIM diseases are ranked for these numbers;
the product list also holds Orphanet entries (held-out top-10 on that list, joined to the truth through
MONDO exact matches: 71.4% now, 41.3% in 0.1.0).

**Two ways to read it, and the second is the one to quote for a new patient.**
- HPO's disease annotations are partly built from these same publications: for 74% of all cases (7,719
  of 10,374; 66% of the held-out cases) the truth disease's HPOA rows cite the case's own PMID. On those
  cases the ranker is partly scored against its own source (held-out top-10 91.2%). The held-out figure
  above mixes the two groups and is an **upper bound**.
- On held-out cases whose paper is *not* among the annotation sources, top-10 is **28.3%** (per case) or
  45.3% averaged per disease. Three single-paper cohorts make up 53% of these 1,785 cases (STXBP1 462,
  ANKRD11 333, SATB2 157), so the per-case and per-disease figures differ; `benchmark.json` carries both
  (`macro_by_disease`).

## What changed, and why: excluded phenotypes

The largest single effect is how excluded ("not present") phenotypes are used. zebra 0.1.0 lowered the
score of every candidate annotated with a phenotype the patient does not have, by that term's information
content times its annotation frequency, summed over the excluded terms. On this dataset that penalty
costs accuracy at every number of excluded terms, in cases whose paper fed the annotations and in cases
whose paper did not:

Held-out top-10 (correct disease), Resnik ranking with the 0/n fix, excluded terms penalised (0.1.0
weight) vs not scored:

| excluded terms | paper was a source: n | penalised | not scored | paper not a source: n | penalised | not scored |
|---|---|---|---|---|---|---|
| 0 | 275 | 93.8% | 93.8% | 601 | 14.6% | 14.6% |
| 1–3 | 612 | 80.7% | 94.6% | 404 | 20.8% | 32.9% |
| 4–10 | 1,141 | 56.0% | 91.8% | 350 | 17.7% | 46.6% |
| >10 | 1,382 | 30.7% | 88.6% | 430 | 4.9% | 28.1% |

Why, as far as the data shows: the excluded terms recorded in these case reports are, overwhelmingly,
features of the true disease that this patient lacks. 88.3% of the excluded terms in the "paper was a
source" cases are annotated to the truth disease (a random OMIM disease: 3.5%); the HPOA frequencies
there (3/10, 0/2 …) were computed from the very table the negatives were copied from. That part is
circular by construction and makes the effect look larger than it will be prospectively; most of the
30-point held-out gain (25 of 30 points) comes from those cases. In the other cases the share is still
38.3% (random disease: 3.9%), consistent with how negatives are recorded in practice: authors and
clinicians note the absence of features they were looking for, which are the features of the disease they
suspect. Either way, a penalty that treats each such absence as evidence against a disease pushes the
right answer down. The size of the effect for a new patient is not known from this benchmark; its
direction held in every stratum, in both halves, leaked or not.

A frequency-aware likelihood-ratio penalty, a cap, weights of 0.1–0.5, and a penalty restricted to
features present in at least 80% of patients were all measured; every one was worse on the development
half than not scoring excluded terms (tuning log below).

So in 0.2 the local ranker **records and checks excluded terms but does not score them**: an excluded
term that contradicts a present one (the same term, or an ancestor such as "seizure" excluded while
"febrile seizure" is present) is dropped with a warning; each candidate annotated with an excluded feature
is flagged in `excluded_hits`; the score is not lowered. `--excluded-weight 1` restores the 0.1.0 penalty.
The library default in `zebra.hpo_local.rank` is unchanged (weight 1), so other callers (the `zebra vcf`
phenotype fit) behave as before; whether they should follow is a separate decision.

This has a visible cost on single cases. The PRD's acceptance query (febrile seizure, focal-onset seizure,
generalized tonic-clonic seizure, global developmental delay; hypotonia excluded) put Dravet syndrome 6th
locally in 0.1.0, because the hypotonia penalty pushed several developmental and epileptic encephalopathies
down; in 0.2 Dravet is 16th locally (6th with `--excluded-weight 1`). The default was chosen on the
development half and confirmed on the held-out half, not on one query.

## Dataset

- GA4GH phenopacket-store release **0.1.27** (published 2026-06-09), asset `all_phenopackets.zip`,
  19,431,098 bytes, sha256 `d0c70005bb09b87035087516252b6b539fbc0415be7c664b73c9cde7eedf205a`
  (equal to the digest GitHub publishes for the asset), retrieved 2026-10-06T09:46:01Z.
  https://github.com/monarch-initiative/phenopacket-store/releases/tag/0.1.27
- 10,377 phenopackets; 10,374 used (3 have no present phenotype). 780 diseases; every diagnosis is a
  single OMIM id; 10,305 cases name a causal gene (CAUSATIVE or CONTRIBUTORY). Present terms per case:
  1–3 in 2,039 cases, 4–6 in 2,466, 7–10 in 2,655, more than 10 in 3,214. Many cases list more than ten
  excluded terms (from the papers' supplementary tables).
- HPO release **v2026-09-01** (`hp.json` sha256 `a7b3a012…f833501`, `phenotype.hpoa` `e89aa39c…ece1f72`,
  whose header says version 2026-09-02, `genes_to_phenotype.txt` `507a17bf…6ab248c`; all three equal to
  the digests GitHub publishes for that release's assets). 29 cases have a truth disease with no HPO
  annotation at all; no present term of any case failed to resolve in this release. 179 cases name a
  causal gene that `genes_to_phenotype.txt` does not list for the truth disease, so no gene list built from
  it can credit them.
- MONDO SSSOM mappings at tag v2026-09-01 (`src/ontology/mappings/mondo.sssom.tsv`, 13,113,756 bytes,
  sha256 `317ac4bf93be1deff4406b4d7d18bb730ead3c4f5c7c814f75d24be7f6ef302f`), used only to join Orphanet
  and MONDO hits to the OMIM truth. 7 truth ids (58 cases) have no MONDO exact match.

## Method

**Cases.** Each phenopacket is read with `zebra.phenopacket.read`: present and excluded HPO terms as
recorded, the diagnosis, the causal genes, the case report's PMID.

**Split.** Diseases are shuffled with a fixed seed (20261006) and each disease goes, with all its cases, to
whichever half has fewer cases so far: development 5,179 cases / 417 diseases, held-out 5,195 cases / 363
diseases. No disease is in both halves (24 papers that report two diseases do span both, as do 34 genes
with two diseases; only global parameters were tuned). Every parameter was chosen on the development half;
the held-out half was scored three times in total: the 0.1.0 code, the 0/n fix alone, the chosen
configuration.

**Leakage flag.** A case is "leaked" when the truth disease's rows in `phenotype.hpoa` cite the case's own
PMID (HPOA's `reference` column). It misses indirect leakage (78 cases whose PMID is cited under another
OMIM id), so the "paper not a source" rows are, if anything, still slightly optimistic.

**Scoring a list.** The position of the first correct item, 1-based; beyond 50 is a miss.
- Disease, primary: only OMIM diseases are ranked, so correct means the truth's OMIM id.
- Disease, product list (OMIM + Orphanet, what `zebra phenotype rank` shows): correct when the hit is the
  truth or is joined to it by MONDO exactMatch, through zebra's own consensus join
  (`zebra.commands.phenotype.resolve_keys`). A broader or narrower class is a miss, as is any id with no
  exact link.
- Gene: the first causal gene symbol in the ranker's gene list.
- Web rankers: Monarch returns MONDO classes, scored through the same joins and the OMIM xrefs Monarch
  returns with each hit (86 of 9,164 of those xrefs are not SSSOM exact matches; none affected a result
  here); PubCaseFinder's OMIM list by id equality.
- Ties: positions are as the ranker orders them (zebra breaks ties by exactly matched terms, then
  annotation specificity; PubCaseFinder serves tied ranks 1, 1, 3). `benchmark.json` also gives
  `ties_fair`, the expected rate if the truth were equally likely anywhere in its tie block; the two differ
  by at most 1.2 points for the headline rows. In 0.2 the truth shares its score with another disease in
  33.0% of the held-out cases where it is found.

**Intervals.** Two 95% intervals per rate (both in `benchmark.json`): Wilson (cases as independent) and a
percentile bootstrap resampling diseases (1,000 replicates, seed 20261006). The tables quote the bootstrap.

## Local ranker by subset (0.1.0 → 0.2)

Correct disease (OMIM-only ranking) and causal gene:

| subset | n | top-1 | top-10 | gene top-10 |
|---|---|---|---|---|
| all cases | 10,374 | 28.1 → 58.0 | 43.5 → 76.1 | 46.2 → 77.3 |
| held-out | 5,195 | 24.9 → 49.4 | 39.5 → 69.6 | 44.2 → 72.0 |
| held-out, paper not a source | 1,785 | 3.6 → 6.8 | 14.6 → 28.3 | 26.5 → 39.6 |
| held-out, paper a source | 3,410 | – | 52.6 → 91.2 | – |
| held-out, 1–3 present terms | 1,085 | 8.7 → 24.3 | 18.2 → 66.6 | 23.9 → 65.6 |
| held-out, 4–6 present terms | 1,232 | 17.9 → 47.8 | 36.4 → 68.6 | 39.9 → 69.1 |
| held-out, 7–10 present terms | 1,289 | 28.3 → 54.0 | 45.2 → 70.1 | 50.6 → 74.4 |
| held-out, >10 present terms | 1,589 | 38.5 → 63.9 | 51.9 → 71.9 | 56.9 → 76.7 |
| held-out, autosomal recessive | 1,443 | 43.0 → 70.2 | 58.4 → 82.5 | 56.5 → 80.7 |
| held-out, autosomal dominant | 3,432 | 17.6 → 41.7 | 32.4 → 66.4 | 39.3 → 68.7 |
| held-out, X-linked | 140 | 17.1 → 32.1 | 33.6 → 46.4 | 42.9 → 69.3 |

Held-out top-10 in 0.2 by the disease's dominant organ system (the top-level HPO branch holding most of
its annotations): nervous system 67.2% (1,911 cases, 135 diseases), head or neck 48.7% (1,208, 48),
musculoskeletal 88.1% (586, 61), neoplasm 99.5% (423, 4 diseases), cardiovascular 79.6% (221, 16), eye
69.0% (200, 17), immune 85.8% (113, 18). With few diseases per subset the intervals are wide (nervous
system 45.9–89.2).

## Tuning log (development half only)

Correct disease (OMIM-only), 5,179 cases. `resnik` is the 0.1.0 scoring with the 0/n fix.

| candidate | top-1 | top-10 | gene top-10 |
|---|---|---|---|
| 0.1.0 code | 31.4 | 47.4 | 48.2 |
| resnik (0/n case counts no longer curated NOT, D-P1-3) | 32.1 | 48.2 | 48.9 |
| **resnik, excluded terms not scored** | **66.7** | **82.7** | **82.5** |
| resnik, excluded weight 0.1 / 0.25 / 0.5 | 60.5 / 54.1 / 44.7 | 78.1 / 72.5 / 63.0 | 76.4 / 71.2 / 62.3 |
| resnik, penalty only if frequency ≥ 80%, weight 1 / 0.25 | 55.1 / 64.2 | 71.4 / 80.4 | 71.6 / 78.8 |
| resnik, excluded not scored, frequency tie-break | 64.8 | 82.1 | 82.8 |
| lr (likelihood ratio, defaults) | 44.5 | 61.0 | 63.2 |
| lr, unknown frequency 0.75 / 0.5 | 46.2 / 46.9 | 63.3 / 64.0 | 64.8 / 65.3 |
| lr, excluded cap 0.5 / 0.99 | 46.1 / 41.5 | 64.8 / 56.9 | 66.7 / 59.9 |
| lr, unexplained term −1 / −2 | 44.6 / 44.7 | 61.4 / 61.4 | 63.4 / 63.6 |
| lr, frequency floor 0.1 | 44.7 | 61.3 | 63.5 |
| lr, penalty only if frequency ≥ 80% / 95% | 54.0 / 54.9 | 70.5 / 71.6 | 72.2 / 72.3 |
| lr, ≥ 80% and cap 0.5 | 56.8 | 75.2 | 76.3 |
| lr, excluded not scored | 53.0 | 74.9 | 76.0 |
| lr, excluded not scored, unknown frequency 0.5 | 54.7 | 76.7 | 77.5 |

Selection rule (written in `tools/bench/tune.py` before the held-out half was scored): the best
development top-10 overall; it becomes the default only if its held-out top-10 beats the shipped code's.
Held-out top-10: 0.1.0 39.5%, 0/n fix alone 39.9%, chosen configuration 69.6%. Kept. The 0/n fix
(D-P1-3) is a correctness fix kept regardless; on its own it moved 0.3 points.

**The likelihood-ratio method** (`--local-method lr`, kept as an option): in the spirit of LIRICAL
(Robinson et al., AJHG 2020), not LIRICAL itself — no genotype, onset or sex model, no pretest
probabilities. For each present term q, ln LR = max over the disease's annotations t of
[IC(MICA(q, t)) + ln f(t)], where f is the HPOA frequency (pooled case counts, Orphanet classes,
percentages); each excluded term adds ln(1 − f) with f capped; a curated NOT on a present term costs
ln 0.1. Frequencies did not help here: diseases annotated from cohorts carry fractional frequencies that
cost ln f per matched term, while classic annotations with no stated frequency count as 1.

## Web rankers and consensus (held-out sample)

Monarch (semsim) on 240 held-out cases and PubCaseFinder on 100 of them (the first 25 of each present-term
stratum), one case per disease, 60 per stratum, seed 20261007 (`tools/bench/results/sample.json`). Run
2026-10-06 against the live services. PubCaseFinder publishes limits of 10 requests a minute, 100 an hour
and 1,000 a day: the runner made 200 requests over about four hours; the request log shows at most 3 in
any minute and 59 in any hour, from all zebra use on the machine. The consensus is rebuilt exactly as the
CLI builds it (top 15 of each source, at least two sources, exact joins), with the 0.2 local ranker.

Correct disease (OMIM truth), percent (Wilson 95% interval):

| list | n | top-1 | top-3 | top-10 | top-25 | gene top-10 |
|---|---|---|---|---|---|---|
| local 0.2 | 240 | 63.3 (57–69) | 73.3 | 81.2 (76–86) | 89.6 | 82.9 |
| local 0.1.0 | 240 | 38.8 (33–45) | 47.9 | 54.6 (48–61) | 60.0 | 57.9 |
| Monarch | 240 | 40.0 (34–46) | 49.6 | 54.2 (48–60) | 58.3 | 47.5 |
| consensus local + Monarch (0.2) | 240 | 47.5 | 52.5 | 53.8 | 53.8 | 47.1 |
| PubCaseFinder | 100 | 42.0 (33–52) | 54.0 | 64.0 (54–73) | 68.0 | 62.0 |
| local 0.2, same 100 | 100 | 64.0 (54–73) | 72.0 | 81.0 (72–87) | 90.0 | 84.0 |
| Monarch, same 100 | 100 | 45.0 | 53.0 | 54.0 | 61.0 | 45.0 |
| consensus of all three (0.2) | 100 | 51.0 (41–61) | 71.0 | 75.0 (66–82) | 76.0 | 74.0 |
| consensus of all three (0.1.0 local) | 100 | 40.0 | 57.0 | 62.0 | 62.0 | 62.0 |

On the 40 sampled cases whose paper was not an annotation source: local 0.2 top-10 45.0%, Monarch 30.0%,
consensus local + Monarch 27.5%; on the 14 such PubCaseFinder cases: local 35.7%, PubCaseFinder 35.7%,
Monarch 14.3%. One case per disease gives 83% leaked cases in this sample (66% in the held-out half), so
the comparison between local and the web rankers is tilted towards local, whose annotations include these
cohorts; Monarch's and PubCaseFinder's knowledge may lag or lead the 2026-09 HPO release.

Differences from the CLI, kept on purpose: PubCaseFinder's Orphanet list was not requested (two requests
per case instead of three), so the consensus has one PubCaseFinder list where the CLI has two; joins come
from MONDO's SSSOM file rather than live Monarch/OLS calls (obsolete MONDO classes are not followed).

## What this does and does not measure

Measured: whether the right OMIM disease and causal gene come back near the top when the phenotype terms
are exactly as curators recorded them from published cases, for 780 Mendelian diseases in phenopacket-store
0.1.27.

Not measured:
- Phenotypes as families or clinicians enter them in zebra (lay phrases mapped through `hpo_search`, fewer
  and less specific terms, Chinese entry). The 1–3-term stratum is the closest proxy.
- Diseases outside phenopacket-store's selection (recent, gene-defined disorders), non-OMIM diagnoses, CNV
  syndromes, mitochondrial disorders (10 held-out cases).
- Anything genotype-driven (`vcf triage`, ACMG). This is phenotype-only ranking.
- The circularity above: for most cases the annotations were derived partly from the case's own
  publication. The "paper not a source" rows are the better estimate for a new patient, and they are low.
- What Claude says with these rankings; `evals/` covers that separately.

## Reproduce

From the repository root, Python 3.9+ and network access (about 35 MB, plus the HPO release if it is not
installed). The local runs take about 4 minutes on 12 cores (20 for the 0.1.0 code); the web sample takes
hours because of PubCaseFinder's limits.

```
python3 -I tools/bench/fetch_store.py --tag 0.1.27 --mondo-tag v2026-09-01 --hpo-tag v2026-09-01
export ZEBRA_HPO_DIR=~/.cache/zebra-mod/bench/hpo/v2026-09-01   # the pinned, digest-checked HPO release
python3 -I tools/bench/cases.py                   # case list, split, leakage flags
mkdir -p /tmp/hpo010 && ln -sf $ZEBRA_HPO_DIR/* /tmp/hpo010/    # the 0.1.0 code keeps its own index there
git show 54c0b4e:zebra/hpo_local.py > /tmp/hpo_local_010.py
python3 -I tools/bench/run_local.py --label base010 --module /tmp/hpo_local_010.py --hpo-dir /tmp/hpo010
python3 -I tools/bench/run_local.py --label v02 --options '{"method": "resnik", "params": {"excluded_weight": 0.0}}'
python3 -I tools/bench/tune.py dev                # every candidate, development half only
python3 -I tools/bench/tune.py heldout resnik resnik_excl0
python3 -I tools/bench/run_remote.py sample --per-bucket 60
python3 -I tools/bench/run_remote.py monarch --per-bucket 60
python3 -I tools/bench/run_remote.py pubcasefinder --per-bucket 25   # ≤ 90/h and 900/day, counted before sending
python3 -I tools/bench/remote_eval.py --label after --json after.json   # the CLI's own local configuration
python3 -I tools/bench/remote_eval.py --module /tmp/hpo_local_010.py --hpo-dir /tmp/hpo010 --label before --json before.json
python3 -I tools/bench/summary.py --before base010 --after v02 --remote-before before.json --remote-after after.json
python3 -I tools/bench/metrics.py base010 v02     # the tables above, as markdown
```

`tests/test_bench_live.py` (`pytest -m live`) re-scores a fixed slice of held-out cases with the shipped
configuration and fails if any recorded position moves, so a change to the ranker has to come with a new
benchmark run. The web results depend on the services' state on the day; answers are cached for seven days.
