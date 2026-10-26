# Changelog

## Unreleased — 2026-10-26 (benchmark only; the plugin is unchanged)

**End to end against published systems.** `tools/bench/agent/`: build HPO-only cases from RareBench, the RDguru
test set, MyGene2 and the phenopacket sample (`prepare.py`); run each through headless Claude Code with this
checkout's plugin, no case open, zebra tools only (`run.py`, `supervise.sh` resumes across usage-limit windows,
`--provider minimax` for a MiniMax arm); score by id, by MONDO descendant, by an LLM judge and by DeepRare's own
judge prompt (`score.py`); compare the tools' lists with the final answer (`tool_vs_final.py`). Opus 5.5 on 1,348
cases: MME 82.5, HMS 53.4, LIRICAL 58.6, RAMEDIS 44.4, MyGene2 67.1 R@1 (DeepRare 77.5 / 57.0 / 56.0 / 73.2 / 75.9);
RDguru set 57.5 vs 41.6. Measured on 0.3.6; 0.3.7's differential doctrine (rankings as a weak input) is not re-measured. Published baselines with sources in `tools/bench/results/agent/BASELINES.md`; reading and
caveats in `docs/BENCHMARK.md`.

## 0.3.9 — 2026-10-24

**Reports that look like documents.** Exported PDFs had tables running off the page (five-column trial tables lost their last column), raw URLs that would not wrap, large blank gaps and no page numbers.
- Tables are sized from their content: every column first gets room for its header, its longest word and, for a column of short labels, the whole label; the rest of the width goes to the columns that hold the text. Fixed layout, words and URLs wrap inside their cell, so nothing overflows; the same widths go into the Word file.
- New look for HTML and PDF: a header with the title and the line under it, numbered teal section headings, striped tables with repeated header rows, provenance labels ([已检索], [原始报告], [凭记忆] and their English forms) drawn as small coloured tags with their source in grey, bare URLs shortened to their address, a section about when to seek care now drawn as a red warning box, page header and "第 n 页 / 共 m 页" footer in print.
- HTML is exported by default with Word and PDF: it opens in any browser or phone, with a contents list for reports of four or more sections. The FSHD clinician summary went from 14 pages with clipped tables to 10.

## 0.3.8 — 2026-10-24

**Sequence-to-function models run inside zebra.** The heavy models no longer need an external command-line tool; `zebra/s2f_backends.py` runs them and every run writes its own receipt (`<case>/evidence/s2f/<run_id>/receipt.json`: inputs, model, version, scores, preflight checks).
- GPN-MSA: the authors' published hg38 table read with a pure-Python tabix reader (HTTP byte ranges into the BGZF file through its index, the 2.8 MB index cached once). No key and no `tabix`; checked live against the published KRAS G12D, BRAF V600E and TP53 R175H values. Now in the default model set for hg38 SNVs. A position missing from the table is `absent`, a transport failure an error, neither a benign call.
- Evo 2: NVIDIA's hosted `forward` endpoint called directly (NVCF queueing handled); the .npy logits decoded in the standard library (float32/float16/bfloat16); the score is the summed causal log-likelihood, ALT minus REF, over the window. evo2-7b first, evo2-40b if it fails, and the model that answered is the one reported.
- AlphaGenome: the official client in a child interpreter (`pip install alphagenome`; `ZEBRA_ALPHAGENOME_PYTHON` to name another one), REF checked against Ensembl before the call; besides the window-wide mean it reports the largest change within ±500 bp, its track and offset.
- `zebra doctor` checks the AlphaGenome client instead of the external CLI and `tabix`. Docs, the `zebra-s2f` skill and the tool description updated.

## 0.3.7 — 2026-10-22

**Never subtracts, beyond the diagnosed case.** A review of 0.3.6 found the 0.3.4 fix reached the diagnosed-patient path only: the undiagnosed path, failed lookups, report export and several skills still held rules that leave a person with less than Claude gives without the mod, and the evidence for "never subtracts" was two co-tuned development cases with p≈0.2–0.5.
- Differential: a ranking is one weak input, not the boundary of the differential. Claude writes its own differential before reading the rankings, adds diagnoses the rankers miss (gestalt, a treatable disorder, a hallmark finding) labelled as its own and checked with `disease_card`, and never drops a candidate only because no ranker listed it. `zebra-diagnose` no longer claims the sources together beat the local ranker (BENCHMARK shows the opposite).
- A failed or timed-out lookup is said in one line and the answer continues from knowledge marked unverified (doctrine rule 4; `zebra-start`, `zebra-therapy`, `zebra-family`, `zebra-literature`, `therapy-scout`).
- Identifiers: a variant copied verbatim from the person's report, or a well-known variant name, may be used directly; a paper or trial known but not fetched is named in words. PMIDs, NCT numbers and ontology ids still only as tools returned them.
- Derived numbers: a zebra tool when zebra implements the rule; otherwise the formula and inputs are shown.
- Roles: a clinician is the care team — hazards, acute management and dose ranges directly (`zebra-safety` sections 1–3 are family wording; textbook contraindications are stated at once, labelled "to verify" when unchecked). Patients and families: no prescription doses, but general nutrition, exercise and supplement ranges and the typical course of an established diagnosis as population ranges. With no clue to the role, a clinical answer plus a plain summary.
- Report export: labelled clinical knowledge and clinician dose ranges survive the audit (`zebra-report`, `evidence-auditor`); a claim is removed only when wrong.
- Answer first: no case, intake, fan-out (now only above ~20 files) or ledger write before the answer. Somatic tumour reports, polygenic risk, pharmacogenomics and NIPT are out of scope. `zebra-safety` loads for urgent symptoms and treatment/procedure questions, not every history. A laboratory-classified pathogenic CNV, exon deletion or SMN1 result is a diagnosis.
- Case-specific facts the FSHD/DM1 rubrics check (17 kb → 3 units, urate tiers, AFP work-up, the mask and wearable lines) removed from `zebra-care` in favour of general rules.
- Repeats: EcoRI/HindIII flank 4.7 kb as the cited protocol gives (was 4.8: 16.3 kb read 3 units instead of 4); HindIII- or BlnI-only sizes refused; 4q/10q warning for EcoRI/HindIII too, `--locus` (`array_locus` in the tool); methylation marked as recorded, not interpreted, with the FSHD2 pattern; haplotype read in any spelling; half-up rounding, band-boundary and 38–40 kb notes, mosaicism and sex not modelled; DMPK >1000 "congenital" only.
- Case: a missing role loads as `unspecified`, `role_confirmed` (cases from 0.3.3 and earlier start unconfirmed; the doctrine asks); one batch may cite tests and relatives it records; `support_reasons` ledger ids validated, `e3` read as E3; removing a finding removes its citations with a note; `case init` language from the locale; `/zebra new` passes the session's language.
- Privacy gate: a two-character Chinese name inside a longer word (张力 in 肌张力) asks instead of denying; standalone, next to a name word or with another identifier it is still denied. An unreadable case file still closes zebra's own tools and anything touching the case, and asks for other calls.
- Tools: `phenotype_rank` described as a weak input; `trials_search` reminds to read `filtered` and other statuses before "no trials"; `case_update` description shorter; `literature_search` no longer deferred, `acmg` and `cnv_interpret` deferred.
- Evals: all cases allow WebSearch/WebFetch so the without-mod arm is not crippled; the clinician no-dosing grader and the "candidates only from phenotype_rank" rubric item removed. New held-out set (`evals/heldout-*`, `evals/HELDOUT.md`): six cases written without access to the mod's instructions (undiagnosed lay Chinese, undiagnosed English clinician, VUS with segregation, an emergency inside long records, a non-neuromuscular diagnosed case, a somatic tumour report), 8 runs, grader hashes and a pre-registered non-inferiority rule. FSHD/DM1 kept as regression cases; BENCHMARK says what their ablations do and do not show.
- Tests: live-network tests skip unless `ZEBRA_NETWORK_TESTS=1` (the full suite no longer hangs): 1467 passed, 145 skipped; hooks 111 passed.

## 0.3.6 — 2026-10-11

**Targets from the documented tier; rubrics that judge substance.**
- `zebra-care`: a guideline target is taken from the tier the records document, not the strictest one (with the default tier when deciding facts are unknown, and what would move the patient); a recommendation class (IIa, Class I) only with its source, otherwise "consider"; closing offers in words, never `/zebra` commands.
- Evals: both diagnosed rubrics rewritten from an itemised re-grade of 0.3.5 — the judge writes a verdict and a quote per item before its verdict; diagnosis wording, D4Z4 units, urate tiers, NIV thresholds, DM1 conduction criteria (QRS ≥120 ms; PR 236 is first-degree block below 240), anaesthesia and GGT spelled out; a new FSHD item and regex grader for the AFP germ-cell work-up; 4 runs per arm.
- Ablation (Opus 5.5, judge Sonnet 5.5, 4 runs per arm, $11.70): FSHD 0.83 with the mod vs 0.61 without (pass 2/4 vs 0/4; the AFP work-up 4/4 vs 2/4); DM1 0.70 vs 0.55 (pass 2/4 vs 1/4). With the mod every run called `repeat_size` and loaded `zebra-care`; no run leaked internal names. The mod's runs cost about 2.7× and take about 1.7× as long.

## 0.3.5 — 2026-10-11

**Numbers by rule, not by recall.** In the 0.3.4 ablation both arms made the same recall errors: an EcoRI fragment of 17 kb read as 4 D4Z4 units (it is 3: the 3 kb flank belongs to the EcoRI/BlnI digest), a uric-acid target too strict for asymptomatic hyperuricaemia, and the mod's replies leaked its own tool names into a clinician's plan.
- New tool `repeat_size` (`zebra repeat d4z4 | dmpk`): D4Z4 units from a 4q35 fragment size by digest (EcoRI (size − 6.8)/3.3, per the NCT04369209 protocol; EcoRI/BlnI and EcoRI/HindIII offsets), with the sizing-error range, the FSHD1 band and its meaning; DMPK CTG classes for DM1. Offline, with sources.
- Doctrine: a number that follows from a rule comes from a tool; a guideline target names its guideline and population; deliverables cite PMIDs, NCT numbers and GeneReviews, never zebra's machinery.
- `zebra-care`: genotype numbers from `repeat_size`; unexplained findings — common, treatable and dangerous causes first (B12 / nitrous oxide for a cord signal with raised homocysteine, germ-cell work-up for AFP in a young man, cardiac troponin I before calling troponin T skeletal); exercise gated on an open cardiac question; no internal names in the plan. `zebra-report` cleans them before export.
- Evals: the FSHD rubric no longer rewards a source-labelling style (that item drove most of 0.3.4's gap); it checks factual accuracy instead, and a regex grader catches internal names. New case `evals/diagnosed-dm1` (DM1, ~600 CTG, PR 236/QRS 126, surgery next month, pregnancy planned) so the gains are not FSHD-specific.
- Ablation (Opus 5.5, judge Sonnet 5.5, 3 runs per arm, $8.45): DM1 0.93 with the mod vs 0.60 without; FSHD 0.625 in both arms (every reply failed the judge). An itemised re-grade found: in FSHD the mod's replies all computed 3 units (the plain ones said 3–4, 4, and 3 with a hedge) and all worked up the AFP for a germ-cell tumour within a month (the plain ones repeated it in 3–6 months), but two of three still gave a uric-acid target of <360 for an asymptomatic patient; in DM1 all six replies were clinically sound, so that gap is mostly judge noise. The rubrics were themselves ambiguous in places (fixed in 0.3.6).
- Tests: `tests/test_repeats.py` (units by digest, bands, warnings, refusals, DMPK classes); a tool-argv test for `repeat_size`.

## 0.3.4 — 2026-10-11

**The mod adds, it never subtracts.** A clinician gave a real case (an adult patient; D4Z4 EcoRI fragment 17 kb, 4qA, D4Z4/PAS methylation 16.2%, the laboratory's conclusion "符合 FSHD1") to Claude Code with and without the mod. Without it, Claude read the report as diagnostic and wrote a management plan. With it, the case was recorded as "no variants", the role was set to `family` because the person wrote short Chinese questions, phenotype ranking put myosin-storage, TTN and GNE myopathies above FSHD, the undocumented facial and shoulder examination was counted against FSHD, an exome was recommended "to confirm FSHD" (an exome cannot see a D4Z4 contraction), and the requested management/prognosis/nutrition plan was never written — rules 1, 2 and 5 forbade it. Intake took 15 minutes (5 subagents, 109 identifiers, 30 HPO searches). The intake's extraction was better than the plain session's in places (an incomplete right bundle branch block the plain session missed), so the structured reading stays; what changes is everything that narrowed the answer.
- Doctrine: a new principle 0 (the mod adds, never subtracts); rule 1 prefers retrieved evidence and asks for recall to be labelled, instead of forbidding it; rule 2 lets evidence the ranker cannot see override a ranking, openly; rule 4 adds that a finding the records do not mention is unknown, never evidence against; rule 5 is matched to the person — clinicians get diagnostic judgment, prognosis, a management plan and dose ranges as proposals for review, families keep the protective register; rule 8 registers only identifiers that point to a person; rule 10 says to answer the question first. The short no-case section says the same in brief.
- Role: never inferred from the language; someone operating Claude Code on a case is a clinician/researcher unless they say otherwise; new cases start `unspecified` (was `family`), and the doctrine says how to treat that.
- Case schema: variant kinds `repeat_contraction` (FSHD1 D4Z4: `fragment_kb`, `haplotype`, `locus`) and `methylation` (`methylation_pct`), and `lab_interpretation` for any finding. A hypothesis's `support`/`against` take ledger ids, case item ids (`v1`, `x2`, `f1` — so a laboratory result can support a diagnosis without a web lookup) and reasons in words (kept in `support_reasons`/`against_reasons` instead of failing the write); status `confirmed` when a lab result establishes the diagnosis.
- `zebra-diagnose`: step 0 — an established diagnosis is recorded and the work moves to care; the ranking runs only for what it does not explain. The test table gains FSHD1 (D4Z4 sizing by PFGE/Southern or optical mapping; missed by exome, genome short-read, CMA, panels) and FSHD2.
- New skill `zebra-care`: a diagnosed patient's plan — decision summary, severity and prognosis as measures to track, respiratory, cardiac, exercise prescription, nutrition, findings the diagnosis does not explain, screening, medicines and anaesthesia, genetics, psychosocial, trials and trial readiness, a follow-up calendar, a fill-in list.
- `zebra-intake` and `phenotype-curator`: proportion (answer first, batch HPO searches, write once); non-sequence genetic results recorded, never "no variant"; long PDFs read whole by script, never "pages 1–3 of 64". `zebra-safety`, `zebra-therapy`, `zebra-start` and `evidence-auditor` follow the role rule; the auditor also flags the opposite failure (an established diagnosis treated as unconfirmed, absence of documentation counted against, the answer replaced by process).
- `zebra case init` warns when the local HPO files are missing (phenotype ranking would skip its local ranker).
- Eval: `evals/diagnosed-management` (de-identified): FSHD1 treated as established, no exome to confirm it, missing examination not counted against, a concrete clinician plan, troponin T handled as skeletal vs cardiac (cTnI), the cord signal handled on its own.
- Tests: `tests/test_case.py` — non-sequence findings, hypotheses citing case items and reasons, a missing item refused, value ranges; the evidence-id test now expects a reason in words to be recorded.

## 0.3.3 — 2026-10-08

- The galloping zebra no longer jumps up and down. 0.3.2 aligned the eleven frames on each panel's ground line, but Muybridge's photographs come from eleven cameras, each framed and lined a little differently: the back jumped by up to 6 dots (a cell and a half) from one frame to the next, in no rhythm a horse has. The frames are now aligned on the horse itself, the line of the back from the croup to the withers: the body stays level (0 dots of movement), the legs and the head move, and in the gathered frames the hooves leave the ground as they do in a gallop.

## 0.3.2 — 2026-10-08

**A zebra that looks like one.** The 0.3.0 sprite was 28 × 12 half-block pixels, hand-placed: one-pixel legs that read as dots, a dog's head, barcode stripes, a grey outline that cluttered dark backgrounds, and only the legs moving across four frames.
- Now eleven positions of one gallop stride traced from Eadweard Muybridge's *The Horse in Motion* (1878, public domain): each silhouette cut from its panel, the rider removed along the line of the back, the frames aligned on the ground and the body's centre so the body rises and falls as a galloping horse's does.
- Drawn in braille (2 × 4 dots to a cell: 44 × 24 dots in 22 × 6 cells, about twice the detail of the half blocks in the same band), with stripes one dot wide that lean back on the rump and forward on the neck, a mane, an eye, solid legs; one colour, the terminal's own text colour, so it reads on light and dark themes alike; the ground is a dim dotted line that runs backwards.
- One stride in about 0.8 s (70 ms a frame). The still pose of the welcome card and the shield is the same zebra.
- `tools/zebra_sprite/build.py` rebuilds the frames from the photograph (numpy, scipy, pillow; a build step, not a dependency of the mod). The sprite test checks every frame decodes to braille cells in the terminal's colour.

## 0.3.1 — 2026-10-08

**Onboarding: after installing, a person knows what to do.** A person who installed 0.3.0 saw "installed" and nothing after it: the installer's last message was mostly the doctor report, the welcome toast waited for a first prompt, and nothing said that zebra-mod is used by asking in plain words, or offered something to try it on.
- `/zebra demo` (and `zebra case demo [dir] --lang zh|en`): a bundled synthetic demo case — a clinic note and a trio exome report (SCN1A c.2134C>T, de novo), Chinese (`~/zebra-cases/demo-xiaoyu`) or English (`demo-lily`) — created or reopened, its identifiers registered so the privacy gate is armed from the first prompt, made the active case with its board open, and its first question put in the prompt box: press Enter and the whole flow runs. Records are copied, never linked, and put back if deleted; a folder holding another case is refused.
- The installer creates the demo case (INSTALL.md step 7, `install.sh`) and ends with a guide written for the person: just ask, with example questions; try `/zebra demo`; what it can do; `/zebra new` for their own records; the privacy line. The doctor report shrinks to one or two lines.
- The first interactive session after installing or updating shows a toast at once (not after the first prompt any more) and a card in the band above the prompt — standing zebra, "just describe symptoms or test results", `/zebra demo`, `/zebra` — until the first prompt or until a case is open.
- The open case moves from a pinned status line to the footer label (`🦓 <title> · HPO n · E n`): Claude Code draws pinned lines as notices with a warning sign, which read as an alarm on a family's case. With `interface: off` it is pinned as before; a case that cannot be read is still pinned (that one is a warning).
- `/zebra` is the guide: how to use it (just ask, three examples), the demo, what it does, the case commands, privacy; Chinese or English from the case or `LANG` (`/zebra zh`, `/zebra en`).
- The router skill answers "what can zebra-mod do / how do I start" the same way.
- Other: the background work `session.start` starts (CLI check, adopting a remembered case, the board's poll) no longer surfaces a failure as an unhandled rejection; it was already best effort.
- Tests: 9 for the mod (`hooks/onboarding.test.tsx`: the demo verb, the guide in both languages, the toast and card at session start and their end) and 4 for the CLI (`tests/test_case_demo.py`: created, copied not linked, identifiers registered and never printed, reused and repaired, another case refused, the home default); checked in a live interactive session: the card at start, `/zebra demo` opening the case and board with the question in the prompt box.

## 0.3.0 — 2026-10-08

**You can see zebra-mod at work.** Until now the only signs were a status line and a pane opened on request; a person could not tell whether a reply used the mod. All of the below is additive and scoped to zebra-mod's own activity, so unrelated work in Claude Code draws exactly as before (`hooks/ui.tsx`, `hooks/ui-sprite.ts`, registered first in `register.tsx` so they wrap its hooks without changing them).
- Each zebra tool call is its own row: `🦓`, what it is (变异卡 / Variant card …), what it is about (the databases being queried are named by the spinner and the band; the row adds them too where the engine draws the row while it runs). Its result is one line — evidence rows and their ledger range, the databases actually queried, how many answers came from the cache, notes — instead of the JSON envelope. A `case_update` row counts the identifiers it registers and never shows them (the engine's own row printed the input, names included).
- The spinner says which databases are being queried (`🦓 正在查询 Europe PMC、PubTator3、LitVar2…`).
- A 28 × 12 pixel zebra gallops in the band above the prompt while a query runs (a `Raster` of half blocks, four frames, repainted with `$.ui.blit` inside the call's own dispatch); after the privacy gate stops a call, the band shows a standing zebra and says why (an identifier in the call, or a case file it could not read) and that the call was not sent, until the next turn. Text only on surfaces without `Raster`; nothing while a survey holds the band.
- A footer label, `🦓 zebra` or `🦓 <case title>` (`🛡` after a refusal), so it is always clear the mod is installed.
- One line at the end of a turn that used the mod: calls, databases, new evidence rows, calls not run, calls the gate stopped (a notice: the model never reads it).
- One toast at the first turn after installing: what the mod does and `/zebra`.
- Chinese or English from the case's language, then from the prompts (any Chinese → Chinese), then from `LANG`.
- New option `interface`: `full` (default), `quiet` (no animation, no spinner text; a toast instead of the band after a refusal) or `off`.
- 19 new tests (`hooks/ui.test.tsx`): rows on the terminal and desktop surfaces, pass-through for other tools, identifiers never drawn, the result line, the spinner, the footer, the band with and without `Raster`, `interface: off`, and a whole turn (start → zebra call → complete) leaving its line; checked in a live interactive session as well.

## 0.2.1 — 2026-10-06

Fixes from an independent adversarial review of 0.2.0 (3 P0, 10 P1 and the cheap P2s), each with a regression test that fails on 0.2.0, and from the first full runs of the eval suite.

**Privacy gate**
- A shell command that reaches the network is scanned whole: query strings after `&`, heredocs, `echo … | curl`, variables set before the call and a script's source no longer slip past (P0). Quotes are respected when splitting commands.
- Every Artifact action is scanned (database writes, comment replies, string replacements), not only a page title (P0).
- Identifiers registered before a case exists are refused in outgoing calls, as with a case open; a command naming another case with `--case` is checked against that case's identifiers.
- Record numbers written with a hyphen (12-0042317), birth dates as 02-Mar-2019 / Mar-02-2019 / 19-03-02, names in HTML entities, split by tags or with tone marks, query-string paths (FHIR `/Patient?name=…`), case folders archived or opened by a script, `gh release upload`, argparse abbreviations of `--prefilter`, and `zebra case $(…)` are all caught; genomic coordinates are no longer read as ID numbers or birth dates; pages are scanned paragraph by paragraph; the email check is linear (a 200 kB sequence took 15 s).
- Report export checks what a reader sees (soft line breaks joined, a name next to Chinese characters, underscores), the title and the file name; links with parentheses and wider table rows render whole.

**Recheck**
- A source that failed keeps its last known answer, so a later change is still reported (P0); a ClinVar or gnomAD outage is "not checked", never "record gone"; paper windows follow each question's last successful check and never overlap; a recruiting list longer than one page is not diffed (its count is).

**Behaviour seen in the eval runs**
- `access` keeps each FDA label's indication text, so "amenable to exon 51 skipping" is quoted, not recalled; the doctrine forbids writing a patient's identifiers back in replies; the skills add heteroplasmy tissue guidance for m.3243A>G, trial eligibility "decided by the trial team", verification of an unlabeled HPO id, and rescue medication in first aid.
- Eval graders judge the final reply (the trace given to a judge is cut at 100,000 characters), leak checks ignore the local identifier registration, the clinic note is in the prompt, and the tools that now answer a question (`access`) count.

**Other**
- Case defaults are copied deeply (one case's relatives never appear in another loaded in the same process); `install.sh` piped through curl no longer treats the current directory as the checkout, and an explicit `ZEBRA_SOURCE` replaces a registered marketplace; INSTALL.md passes every option.

## 0.2.0 — 2026-10-06

Making the claim literally true — a rare-disease workstation for families on a diagnostic odyssey, clinicians and researchers. Three independent reviews of 0.1.0 set the list; every change below carries a regression test that fails on 0.1.0.

**Privacy and permissions**
- The mod's own tools now go through the permission chain and the privacy gate (in 0.1.0 they bypassed both: a registered name could reach Europe PMC). Read-only lookups skip the prompt unless a rule says otherwise; case writes follow the permission mode, with "allow for this session".
- The gate reads commands as the shell will (quotes, `bash -c`, backticks, process substitution, `/dev/tcp`, runners such as `uv run` and `Rscript`), matches dates of birth in every common order and spelling, never matches a record number inside a scientific id (HP:, rs, NM_, PMID, coordinates), checks names in free-form keys, asks before case files or genome data leave by scp/rsync/aws/gsutil/gh gist/pipes or a synced folder (Dropbox, iCloud, OneDrive, Nutstore, Baidu), asks before the MyVariant whole-exome prefilter, and refuses outbound calls if it crashes.
- Identifiers can be registered before a case exists (held for the session, never written, added to the next case); with a case open they are refused in every outgoing call. Without a case, an email or phone number in another tool's call is asked about, not refused.
- Installed, the mod no longer changes unrelated work: the default doctrine mode `auto` adds the full rules only while a case is open, otherwise a short section that applies to rare-disease questions only.

**Diagnosis**
- Benchmark on GA4GH phenopacket-store (10,374 cases, held-out split, leakage strata): correct disease in the local top 10 for 28% of held-out cases whose paper is not an annotation source (0.1.0: 15%), 70% across all held-out cases (an upper bound). docs/BENCHMARK.md.
- Excluded terms are flagged (`excluded_hits`), not scored; ties reported; curated NOT annotations pooled correctly; ~50 Chinese lay phrases, inside short sentences, negation-aware, also without the local HPO files.
- Phenopacket v2 export (no free text, no identifiers) and v1/v2 import; `evals/` with 10 end-to-end cases for `claude plugin eval`.

**Variants, CNVs, S2F**
- ACMG: PM2 re-derived from ClinGen SVI PM2 v1.0, Whiffin formula chosen by inheritance (BS1 for F508del no longer offered), PM1+PP3 cap kept at 4 points, BS4 with its strength; SpliceAI at ±4,999 nt in `acmg suggest`; PVS1 and PS1/PM5 inputs; GRCh37 input mapped to GRCh38; mtDNA through gnomAD's mitochondrial data, heteroplasmy and MITOMAP; LitVar counts exclude other alleles; "SCN1A c.2134C>T" read on the MANE transcript with a warning.
- CNVs: ClinGen curated regions with HI/TS and the ACMG section-2 row (22q11.2 → ISCA-37446, HI 3, 2A); every gene checked; signed intronic offsets; SMN1/SMN2 clinical exon numbering and copy-number reading; FMR1/HTT/DMPK/FXN/C9orf72 size bands with quoted sources.
- `aso_screen` (splice-switching antisense feasibility), `expression` (GTEx tissues), MaveDB scores under their own calibration.

**Reanalysis and QC**
- Whole-exome triage through a MyVariant frequency prefilter (cold 556 s, warm 3 s, resumable); SpliceAI on the top splice/non-coding candidates; PED input; CNV-only VCFs handed to `zebra cnv`.
- `zebra qc`: sex check, KING kinship, runs of homozygosity, Mendelian errors and uniparental disomy, mosaic de novo calls.

**China**
- `access`: FDA/EMA status from the agencies' records, NMPA/CDE approval from bundled official documents, the 2025 national reimbursement list with restriction text, trials with sites in China, collaboration-network hospitals by province; Chinese-cohort allele frequencies (NyuWa, WBBC, 1000G East Asian, Taiwan Biobank); stricter Chinese name matching (糖尿病 is no longer maple syrup urine disease).

**Families and follow-up**
- `report_export`: Word, PDF (local browser or LibreOffice) and HTML with Chinese typography, refusing a report that still holds an identifier; the for-a-family (代操作) flow for family letters and visit-preparation sheets.
- `case_recheck`: what changed since the last check — ClinVar, ClinGen validity, recruiting trials, new papers — into the case timeline.
- Case model: family members, tests already done, timeline, identifiers.

**Install and engine**
- One sentence ("install https://github.com/zwbao/zebra-mod") through INSTALL.md, or `install.sh`; `zebra doctor` marks optional keys as optional.
- HTTP: an error body sent with 200 is evicted, accepted not-found answers live a day, GraphQL error 500s are not retried, per-host pacing for the new sources, PubCaseFinder's hourly and daily limits counted across processes; results are trimmed before provenance.

## 0.1.0 — 2026-10-05

First release.

- **Mod** (`hooks/register.tsx`): 15 model-callable tools backed by the `zebra` CLI; research doctrine in the system prompt (`always` / `case` / `off`); case board pane and status line; `/zebra` command (`new`, `case`, `board`, `ledger`, `doctor`, `close`); auto-approval of the mod's own read-only tools; privacy gate (denies outgoing calls carrying a case's protected identifiers or ID-number/phone/email patterns, including Bash calls to the zebra CLI; asks before raw genome files leave the machine).
- **Engine** (`zebra/`, Python standard library, ≥ 3.9): HTTP layer with caching, retries and provenance; case workspace with an append-only evidence ledger; ACMG/AMP points (Tavtigian 2020) beside the 2015 combining rules, ClinGen SVI warnings and calibrated suggestions (REVEL per Pejaver 2022, SpliceAI per Walker 2023, BA1/BS1/PM2_Supporting); classical statistics; offline HPO ranking (Resnik best-match average) with English and official Chinese label search; Ensembl VEP client with AlphaMissense/REVEL/CADD/SpliceAI; sources for HPO, Monarch, PubCaseFinder, Orphanet, OLS, GeneReviews, gnomAD, ClinVar, ClinGen, PanelApp, UniProt, Europe PMC, PubTator, LitVar, ClinicalTrials.gov, Open Targets; S2F layer (SpliceAI/Pangolin lookup; AlphaGenome, Evo 2, GPN-MSA); base-editing feasibility screen; VCF inspection and trio triage; `zebra doctor`.
- **Skills** (12, including `zebra-safety`: urgent red flags and disease-specific drug, anaesthesia and procedure hazards) and **subagents** (6).
