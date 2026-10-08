# Changelog

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
- **Engine** (`zebra/`, Python standard library, ≥ 3.9): HTTP layer with caching, retries and provenance; case workspace with an append-only evidence ledger; ACMG/AMP points (Tavtigian 2020) beside the 2015 combining rules, ClinGen SVI warnings and calibrated suggestions (REVEL per Pejaver 2022, SpliceAI per Walker 2023, BA1/BS1/PM2_Supporting); classical statistics; offline HPO ranking (Resnik best-match average) with English and official Chinese label search; Ensembl VEP client with AlphaMissense/REVEL/CADD/SpliceAI; sources for HPO, Monarch, PubCaseFinder, Orphanet, OLS, GeneReviews, gnomAD, ClinVar, ClinGen, PanelApp, UniProt, Europe PMC, PubTator, LitVar, ClinicalTrials.gov, Open Targets; S2F bridge (SpliceAI/Pangolin lookup; AlphaGenome, Evo 2, GPN-MSA via s2f-penguin); base-editing feasibility screen; VCF inspection and trio triage; `zebra doctor`.
- **Skills** (12, including `zebra-safety`: urgent red flags and disease-specific drug, anaesthesia and procedure hazards) and **subagents** (6).
