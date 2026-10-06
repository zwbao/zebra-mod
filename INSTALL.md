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
     --config python=python3 --config doctrine=auto --config privacyGate=true
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

6. **Check.** `python3 ROOT/bin/zebra doctor`. Report which checks passed. Network sources that
   fail here usually mean a firewall or proxy; heavy sequence-to-function models (AlphaGenome,
   Evo 2) are optional and need their own keys — mention them only if the person does research.

7. **Tell the person** (in their language) that zebra-mod is installed and takes effect in a **new**
   Claude Code session (exit and run `claude` again), or in this one after `/reload-plugins`; then
   type `/zebra-mod:zebra-start`. To work on a
   patient's records: `/zebra new ~/cases/<name>`, put the files in its `records/` folder, and ask.
   Mention once that what Claude reads from records reaches the model provider as in any Claude
   Code session, and that zebra-mod keeps registered names and record numbers out of every outgoing
   database query.

## By hand

```bash
curl -fsSL https://raw.githubusercontent.com/zwbao/zebra-mod/main/install.sh | sh
```

or, from a clone (needed while the repository is private):

```bash
git clone https://github.com/zwbao/zebra-mod ~/zebra-mod
~/zebra-mod/install.sh
```

`install.sh` does steps 1–6 above. Options: `ZEBRA_PYTHON=/path/to/python3`,
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
