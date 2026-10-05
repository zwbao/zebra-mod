# zebra-mod backlog

From three independent adversarial reviews of 0.1.0 (engineering, scientific correctness, first-principles product and safety), 2026-10-05/06. Full reports with reproductions are kept outside the repo; every item below names the file and the defect so it can be reproduced from the code.

Status: `[ ]` open · `[x]` fixed in this tree · `[>]` deliberately deferred with a reason.

## P0 — wrong result, unsafe statement, data leak, or a blocked core user

### Scientific correctness
- [x] **S1** `acmg.py` BA1 fires on the gnomAD point estimate with no filtering AF and no exception list → HFE p.Cys282Tyr (ClinVar Pathogenic) was called stand-alone Benign. Require faf95 (or AN ≥ 2000 and the grpmax estimate), and never offer BA1 for a variant with a ClinVar P/LP assertion without saying so.
- [x] **S2** `acmg.py` BP4 offered for missense from SpliceAI alone; Walker 2023 allows that only once protein impact is excluded.
- [x] **S3** `vcf.py` gnomAD allele **counts** parsed as frequencies (`"af" in "afr"`), so every annotated variant was dropped and triage answered "no candidate".
- [x] **S4** `sources/variant.py` + `commands/acmg.py`: the BS1 branch was unreachable — `max_credible_af` was never passed, so a variant too common to be fully penetrant never got BS1.
- [x] **S5** `acmg.py` `BP4_VeryStrong` alone classified a variant Benign from one computational score (Pejaver 2022: predictors alone cannot classify).
- [x] **S6** `s2f.py` emitted "≤ 0.1 supports no effect (BP4/BP7)" for every variant class, handing BP7 to missense variants.
- [x] **S7** `sources/gnomad.py` grpmax dropped the Middle Eastern group outside its proper scope, under-reporting founder alleles (MEFV p.Met694Val).

### Engineering
- [x] **E1** Any session in any folder adopted the last case used anywhere (global `$.store`) and wrote evidence into that patient's ledger.
- [x] **E2** Updating a hypothesis without `status` demoted `leading`/`confirmed` to `considered`.
- [x] **E3** `case_update` stored wrong-typed values as character lists (`support:"E12"` → `E,1,2`), and a non-string `disease` broke every later write.
- [x] **E4** The privacy gate was bypassed by URL-encoding, `+`-joining, re-spacing, reordering, `\u` escapes or full-width characters.
- [x] **E5** One "case identifiers" anywhere in a Bash command exempted the whole command.
- [x] **E6** File contents were never checked (`curl -d @records/…`), and `Artifact`, remote `Agent` and `SendMessage` were not gated at all.
- [x] **E7** `sources/variant.py`: a variant given in the wrong genome build got a fabricated allele, then read as "absent from gnomAD at a covered site" (feeding PM2) and "no ClinVar record".
- [x] **E8** `NCBI_API_KEY` was placed in URLs that flow into sources, ledger rows and reports.
- [x] **E9** Without local HPO files, `hpo_search "seizure"` never returned HP:0001250.
- [x] **E10** `stats km` hung forever on a NaN time or non-0/1 event codes, and miscounted tied events.
- [x] **E11** `vcf triage` filtered out common pathogenic recessive alleles (an F508del homozygote got "no candidate").

### Product and safety
- [x] **P1a** `/zebra` was refused at registration because the router skill reserved the name, so case management was unreachable in a real session (observed). The router skill is now `zebra-start`.
- [x] **P1b** README claimed records "stay on this machine" while every record the model reads is sent to the model provider. Said plainly now.
- [x] **P1c** No emergency red flags and no contraindication warnings anywhere (sodium-channel blockers in Dravet, anaesthesia risk in DMD, HLA-B\*15:02 and carbamazepine in Han Chinese patients).
- [x] **P1d** The commonest Chinese report forms could not be entered: CNV/CMA results, DMD exon deletions, SMN1 copy number, repeat expansions.
- [x] **P1e** Variants in mitochondrial and non-coding genes were annotated to the wrong gene (m.3243A>G → MT-ND1 instead of MT-TL1): transcript choice ignored overlap.
- [x] **P1f** `zebra china` and `zebra disease` answered "not on the list" for SMA, DMD, Prader-Willi, Dravet and others under their common Chinese names, all of which are on the 2018 list.
- [x] **P1g** "Approved" was not trustworthy: a withdrawn EMA authorisation showed as APPROVAL, and approvals were not split by jurisdiction.
- [x] **P1h** The repository did not exist, so every install command in the README failed for anyone else.

## P1 — materially weakens the claim

- [x] **F1** `tool.check` auto-allowed all 15 tools without consulting the user's rules, plan mode or deny rules, `case_update` writes included.
- [x] **F2** Privacy-gate false positives: emails, `ssh user@host`, `git@gitlab.com`, two-letter pinyin identifiers matching inside JSON keys, years, long rsIDs.
- [x] **F3** The gate failed open when `case.json` could not be read, and had a 0–4 s window after `identifiers --add`.
- [x] **F4** The gate denied local `case_update` writes with a false "would send … off this machine".
- [x] **F5** Results over 60 k characters lost their warnings, sources and ledger ids, and the JSON was left invalid.
- [x] **F6** The CLI's retry budget exceeded the tool timeout, so one hung source lost every other source's answer.
- [x] **F7** Non-JSON stdout surfaced as "ProcessError: SyntaxError"; positional values were not separated with `--`, so `variant: "--help"` printed help; tracebacks were cut at the wrong end.
- [x] **F8** `case_update` could crash part-way and leave the case half-written; unknown keys were ignored while reporting success.
- [x] **F9** `/zebra new "my case"` created a folder named `"my`; paths with spaces were impossible; a bare `~` resolved to `cwd/~`.
- [x] **F10** HTTP 200 bodies that were errors, HTML or empty were cached for 7–30 days.
- [x] **F11** `core.attempt` did not catch AttributeError, so one upstream shape change killed a whole card.
- [x] **F12** A deterministic HTML 404 was retried with backoff and the same call made twice (62 s for one disease query).
- [x] **F13** `therapy` built a landscape for the first fuzzy match ("SMA" → proximal SMA) without a warning.
- [x] **F14** `phenotype rank` with no usable terms still queried Monarch with an empty set and reported success.
- [x] **F15** `gene` built a card for a junk symbol, and multi-word input produced a false "previous symbol" claim.
- [x] **F16** Triage never attached ClinVar significance to indels.
- [x] **F17** NaN/Infinity were emitted as invalid JSON, and `faf95 nan` produced a wrong verdict.
- [x] **F18** `stats` accepted out-of-range inputs (negative frequencies, penetrance 1.5) and showed raw tracebacks.
- [x] **F19** `stats denovo` p-values lost precision or underflowed to 1.0.
- [x] **F20** `s2f_predict` ran models sequentially at 540 s each under a 600 s cap, losing every result and orphaning the child process.
- [x] **F21** A truncated `.vcf.gz` produced a traceback with no envelope.
- [x] **F22** Malformed VCF data lines were dropped silently.
- [x] **F23** `acmg suggest` could offer BP4 alongside PP3, BP4 with no splice prediction, or two PP3s.
- [x] **F24** `acmg suggest` dropped BA1/BS1 silently when the allele number or the maximum credible AF was missing.
- [x] **F25** PM2 came from AC = 0 with no coverage check; no PM2 for dominant disorders unless AC was exactly 0.
- [x] **F26** Any code could be modified to any strength (`PP1_VeryStrong` → 8 points), contradicting zebra's own segregation ceiling.
- [x] **F27** The PM1 + PP3 cap was not applied to the 2015 combining-rule reading.
- [x] **F28** Segregation → PP1 under-counted against ClinGen's current guidance, and the cited threshold table was not the one implemented; no BS4 at all, though two skills instruct it.
- [x] **F29** The HPO excluded-term penalty was documented as a sum and implemented as a maximum; HPOA `NOT` annotations were discarded.
- [x] **F30** One editing window was used for both ABE and CBE, with no CBE by-product or guide-independent deamination caveat.
- [x] **F31** Triage: a de novo candidate's compound-heterozygous partner was held to the dominant AF cut-off; missing parental genotypes read as "not carried"; a diploid male X was misread without `--sex`; phase was parsed but unused.
- [x] **F32** S2F: the 0.1–0.2 uninformative band was unnamed; the PP3 supporting cap unstated; the ±500 window differs from the calibration's ±4,999; `--mask` help contradicted the thresholds; Pangolin's headline was an undisclosed maximum over tissues; SpliceAI and Pangolin were treated as independent.
- [x] **F33** "Covered site" was judged on mean depth while the fraction over 20× was fetched and ignored.
- [x] **F34** The doctrine was unenforced beyond HPO ids and ACMG arithmetic: a non-existent OMIM id, an evidence id with an empty ledger and a "confirmed" status were all accepted without warning.
- [x] **F35** Offline phenotype ranking produced large ties broken by id order, and triage used its top 200 genes as a hard filter.
- [x] **F36** HPO search was literal: common Chinese lay phrases found nothing or the wrong term; a fresh install returned a raw error for Chinese input.
- [x] **F37** Transcript and build handling: legacy transcript versions errored; a bare `GENE:c.` switched transcript silently; a wrong build produced another gene's card behind one warning.
- [x] **F38** Trials: no eligibility criteria or contacts, no ChiCTR, and stale "recruiting" statuses unflagged.
- [x] **F39** The family skill asked for patient organisations and expert centres that no source returned; no GeneReviews text, no plain-language source.
- [x] **F40** Literature results were not relevance-checked (an unrelated oncology paper for an SCN1A variant).
- [>] **F41** Researcher infrastructure: response snapshots, freeze/replay, batch mode, PED/cohort input, Phenopacket exchange. Deferred to 0.3 — each is a feature, not a defect; `--json` plus the ledger covers provenance meanwhile.
- [>] **F42** Deeper ACMG automation (PVS1 decision tree, PM1 hotspots, PS1/PM5, gene-specific VCEP rules). Deferred to 0.3: these are judgement codes the skills ask the model to justify, and a half-built decision tree is worse than none.

## P2 — polish and next

Kept in the reports: cache temp-file races, ledger ids from line counts, `fsync` on case writes, board text escaping, `uploadsGenome` and `isOutboundShell` gaps (`nc`, `tar | ssh`, `gh gist`, `mail`, DNS), deferring tool schemas until a case is active, registering tools before the environment checks, two interpreters for one CLI, structural validation of a hand-edited `case.json`, Chinese UI strings, PDF/Word report export, Chinese population frequencies, ChiCTR trials, and a Phen2Gene-class second gene ranker.
