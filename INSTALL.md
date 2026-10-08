# Installing zebra-mod

This file is written so that Claude Code itself can follow it. A person can give Claude Code
the repository's address and say "install this" ("帮我安装 https://github.com/zwbao/zebra-mod");
Claude Code reads this file and runs the steps below. A person can also run them by hand.

## For Claude Code: install it for the person

Do the steps in order. Run each command with Bash, read its output, and stop at the first failure
with the reason in plain words (in the person's language). Do not use sudo and do not pip-install
anything: zebra-mod needs neither.

1. **Claude Code version.** `claude --version` must print 2.1.289 or later. If older: tell the
   person to run `claude update`, then start again.

2. **Python 3.9 or later.** Find one: `python3 --version` (or `python --version`). If there is none:
   - with `uv` available: `uv python install 3.12`, then use the path from `uv python find 3.12`;
   - macOS without uv: ask the person to run `xcode-select --install` (it installs python3), or
     `brew install python`;
   - Linux: ask them to install python3 with their package manager;
   - Windows: https://www.python.org/downloads/ (tick "Add python.exe to PATH").

3. **Register and install.** If the repository has been cloned to a folder, use that folder as the
   source; otherwise use `zwbao/zebra-mod`.
   ```bash
   claude plugin marketplace add zwbao/zebra-mod        # or: claude plugin marketplace add /path/to/clone
   claude plugin install zebra-mod@zebra-mod --scope user \
     --config python=python3 --config doctrine=auto --config privacyGate=true --config interface=full
   ```
   If the Python from step 2 is not `python3` on PATH, add `--config python=/full/path/to/python`
   to the install command. If the marketplace is already registered, run
   `claude plugin marketplace update zebra-mod` and `claude plugin update zebra-mod@zebra-mod`
   instead. A private repository needs the person's GitHub access: if `marketplace add` fails with an
   authentication error, say so and ask them to clone the repository first.

4. **Find where it was installed.** `claude plugin list --json` → the entry with
   `"id": "zebra-mod@zebra-mod"` → its `installPath`. Call it ROOT below.

5. **Offline phenotype data (recommended, ~80 MB, once).** `python3 ROOT/bin/zebra hpo fetch`.
   If it fails (network), zebra-mod still works online; tell the person they can run it later.

6. **Check.** `python3 ROOT/bin/zebra doctor`. Keep what you tell the person about it to one or two
   lines: everything passed, or which sources failed and what that leaves unavailable (a failure
   usually means a firewall or proxy). Optional keys (AlphaGenome, Evo 2, NCBI, OMIM, `s2f`) are
   not failures: mention them only if the person does research.

7. **The demo case.** `python3 ROOT/bin/zebra case demo --lang zh` for a person who writes Chinese,
   `--lang en` otherwise. It creates `~/zebra-cases/demo-xiaoyu` (or `demo-lily`): two synthetic
   records (a clinic note and a genetic report) with their identifiers already registered, so the
   person can try the whole flow before using real records. Running it again reuses the folder.

8. **Tell the person how to use it** — this is the part they read, so write it for them, in their
   language, in this order, and keep the installation details short:
   1. It is installed and takes effect in a **new** Claude Code session (exit and run `claude`
      again), or in this one after `/reload-plugins`.
   2. **How to use it: just ask.** No commands to learn. They describe symptoms, paste test results
      or a genetic report, or ask about a variant, a disease, a treatment or a trial, in their own
      words; Claude calls zebra-mod's tools by itself. Give two or three example questions (for
      instance: 「孩子 6 个月开始发热抽搐，基因报告说 SCN1A c.2134C>T，这是什么意思？」,
      「Dravet 综合征现在有哪些获批的药和在招募的临床试验？我们在中国」).
   3. **Try it now:** in the new session type `/zebra demo`. It opens the demo case from step 7 and
      puts the first question in the prompt box; they press Enter and watch it work.
   4. **What it can do**, one line each: differential diagnosis from symptoms; variant
      interpretation with ACMG points computed by code; CNVs, exon deletions and splicing; approved
      treatments and recruiting trials (with China's approvals, reimbursement and trial sites);
      recurrence risk and other statistics; a family letter and a visit-preparation sheet exported to
      Word/PDF.
   5. **Their own records:** `/zebra new ~/cases/<name>`, put the files in its `records/` folder, and
      ask. `/zebra` shows this guide again at any time.
   6. Once, plainly: what Claude reads from records reaches the model provider as in any Claude Code
      session; zebra-mod keeps registered names and record numbers out of every outgoing database
      query. It is research-grade and does not replace a clinician.

## By hand

```bash
curl -fsSL https://raw.githubusercontent.com/zwbao/zebra-mod/main/install.sh | sh
```

or, from a clone (needed while the repository is private):

```bash
git clone https://github.com/zwbao/zebra-mod ~/zebra-mod
~/zebra-mod/install.sh
```

`install.sh` does steps 1–7 above and prints the guide of step 8. Options: `ZEBRA_PYTHON=/path/to/python3`,
`ZEBRA_SKIP_DATA=1` (skip the HPO download), `ZEBRA_SOURCE=<owner/repo, git URL or folder>`.

## Try it without installing

```bash
claude --plugin-dir ~/zebra-mod
```

loads the plugin for that one session only.

## Uninstall

```bash
claude plugin uninstall zebra-mod@zebra-mod
claude plugin marketplace remove zebra-mod
rm -rf ~/.cache/zebra-mod          # downloaded HPO data and the HTTP cache
```

Case folders you created are yours and are not touched.
