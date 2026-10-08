#!/bin/sh
# zebra-mod installer.
#
#   from GitHub:     curl -fsSL https://raw.githubusercontent.com/zwbao/zebra-mod/main/install.sh | sh
#   from a checkout: ./install.sh
#
# It checks Claude Code and Python, registers the zebra-mod marketplace and installs the
# plugin for your user, downloads the HPO release once (~80 MB, for offline phenotype
# ranking), and runs `zebra doctor`. Nothing needs sudo; nothing is installed with pip.
#
# Environment: ZEBRA_SOURCE   marketplace source (owner/repo, git URL or a local folder)
#              ZEBRA_PYTHON   the Python 3.9+ to use
#              ZEBRA_SKIP_DATA=1  skip the HPO download
set -eu

MIN_CC="2.1.289"

say() { printf '%s\n' "$*"; }
step() { printf '\n==> %s\n' "$*"; }
die() { printf '\nzebra-mod install stopped: %s\n' "$*" >&2; exit 1; }

version_ge() { # is $1 >= $2 (dotted numbers)
  awk -v a="$1" -v b="$2" 'BEGIN { n = split(a, x, "."); split(b, y, ".");
    for (i = 1; i <= 3; i++) { if ((x[i] + 0) > (y[i] + 0)) exit 0; if ((x[i] + 0) < (y[i] + 0)) exit 1 } exit 0 }'
}

python_ok() { # a working Python 3.9+
  "$1" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)' >/dev/null 2>&1
}

# ---------------------------------------------------------------- Claude Code
step "Claude Code"
command -v claude >/dev/null 2>&1 || die "the 'claude' command was not found. Install Claude Code first: https://docs.claude.com/en/docs/claude-code"
CC_VERSION=$(claude --version 2>/dev/null | awk '{ print $1 }')
[ -n "$CC_VERSION" ] || die "'claude --version' printed nothing"
version_ge "$CC_VERSION" "$MIN_CC" || die "Claude Code $CC_VERSION is too old: zebra-mod needs $MIN_CC or later. Run: claude update"
say "Claude Code $CC_VERSION"

# ---------------------------------------------------------------- Python
step "Python 3.9+"
PY=""
for c in ${ZEBRA_PYTHON:-} python3 python3.13 python3.12 python3.11 python3.10 python3.9 python; do
  if command -v "$c" >/dev/null 2>&1 && python_ok "$c"; then PY="$c"; break; fi
done
if [ -z "$PY" ] && command -v uv >/dev/null 2>&1; then
  say "No Python 3.9+ on PATH; installing one with uv (no sudo needed)…"
  uv python install 3.12 >/dev/null 2>&1 || true
  c=$(uv python find 3.12 2>/dev/null || true)
  if [ -n "$c" ] && python_ok "$c"; then PY="$c"; fi
fi
if [ -z "$PY" ]; then
  case "$(uname -s 2>/dev/null)" in
    Darwin) hint="run 'xcode-select --install' (it includes python3) or 'brew install python', then run this again" ;;
    Linux) hint="install it with your package manager (e.g. 'sudo apt install python3'), then run this again" ;;
    *) hint="install it from https://www.python.org/downloads/ and run this again" ;;
  esac
  die "no Python 3.9 or later was found: $hint"
fi
PY_PATH=$(command -v "$PY" 2>/dev/null || printf '%s' "$PY")
say "$("$PY" -c 'import sys; print("Python %d.%d.%d" % sys.version_info[:3])') at $PY_PATH"

# ---------------------------------------------------------------- plugin
step "Register and install the plugin"
# $0 is this script only when it was run from a file; piped (curl … | sh) it is the shell, and the
# current directory says nothing about where zebra-mod is
HERE=""
if [ -f "$0" ] && [ "$(basename -- "$0")" = "install.sh" ]; then
  HERE=$(CDPATH= cd -- "$(dirname -- "$0")" 2>/dev/null && pwd || true)
fi
if [ -n "${ZEBRA_SOURCE:-}" ]; then
  SOURCE="$ZEBRA_SOURCE"
elif [ -n "$HERE" ] && [ -f "$HERE/.claude-plugin/marketplace.json" ]; then
  SOURCE="$HERE" # run from a checkout: install that copy
else
  SOURCE="zwbao/zebra-mod"
fi
say "source: $SOURCE"

if [ -n "${ZEBRA_SOURCE:-}" ] && claude plugin marketplace list 2>/dev/null | grep -q "zebra-mod"; then
  claude plugin marketplace remove zebra-mod >/dev/null 2>&1 || true  # an explicit source replaces the old one
fi
if claude plugin marketplace list 2>/dev/null | grep -q "zebra-mod"; then
  claude plugin marketplace update zebra-mod || die "could not update the zebra-mod marketplace"
else
  claude plugin marketplace add "$SOURCE" || die "could not add the marketplace from '$SOURCE'. If the repository is private you need access to it: clone it with git and run ./install.sh inside the clone"
fi

INSTALLED=$(claude plugin list --json 2>/dev/null | "$PY" -c 'import json, sys
try:
    rows = json.load(sys.stdin)
except ValueError:
    rows = []
print(next((r.get("installPath", "") for r in rows if r.get("id") == "zebra-mod@zebra-mod"), ""))' || true)
if [ -n "$INSTALLED" ]; then
  claude plugin update zebra-mod@zebra-mod || die "could not update zebra-mod"
else
  # every option given, so the install does not end on "options not yet set": these are the defaults,
  # except python when the interpreter found is not plain python3
  PY_OPT="python3"
  [ "$PY" = "python3" ] || PY_OPT="$PY_PATH"
  claude plugin install zebra-mod@zebra-mod --scope user --config "python=$PY_OPT" --config "doctrine=auto" \
    --config "privacyGate=true" --config "interface=full" || die "could not install zebra-mod"
fi

ROOT=$(claude plugin list --json 2>/dev/null | "$PY" -c 'import json, sys
rows = json.load(sys.stdin)
print(next((r.get("installPath", "") for r in rows if r.get("id") == "zebra-mod@zebra-mod"), ""))')
[ -n "$ROOT" ] && [ -f "$ROOT/bin/zebra" ] || die "zebra-mod is not listed by 'claude plugin list' after installing"
say "installed at $ROOT"

# ---------------------------------------------------------------- data and check
if [ "${ZEBRA_SKIP_DATA:-0}" != "1" ]; then
  step "HPO release for offline phenotype ranking (~80 MB, once)"
  "$PY" "$ROOT/bin/zebra" hpo fetch || say "(the download did not finish: zebra-mod still works online; run '$PY $ROOT/bin/zebra hpo fetch' later)"
fi

step "zebra doctor"
"$PY" "$ROOT/bin/zebra" doctor || say "(some checks did not pass: the lines above say which and why)"

step "demo case (synthetic records to try it on)"
DEMO_LANG="en"
case "${LC_ALL:-${LANG:-}}" in zh*) DEMO_LANG="zh" ;; esac
"$PY" "$ROOT/bin/zebra" case demo --lang "$DEMO_LANG" || say "(the demo case was not created: /zebra demo creates it later)"

cat <<EOF

zebra-mod is installed. It takes effect in a new Claude Code session (run: claude), or in an open
one after /reload-plugins.

How to use it: just ask. Describe symptoms, paste test results or a genetic report, or ask about a
variant, a disease, a treatment or a trial, in your own words; Claude calls zebra-mod by itself.
  e.g. "My daughter has had seizures with fever since 6 months; her report says SCN1A c.2134C>T.
        What does it mean?"
Try it now:   /zebra demo                  opens a synthetic demo case, the first question ready
Your records: /zebra new ~/cases/<name>    then put the files in its records/ folder and ask
Everything:   /zebra

zebra-mod 已安装，在新的 Claude Code 会话（运行 claude）或已打开会话中输入 /reload-plugins 后生效。

怎么用：直接提问。用自己的话描述症状、贴检查结果或基因报告，或者问某个变异、疾病、药物、临床试验，
Claude 会自己调用 zebra-mod。
  例如：「孩子 6 个月开始发热抽搐，基因报告说 SCN1A c.2134C>T，这是什么意思？」
先试一下：/zebra demo                 打开一个合成示例病例，第一个问题已放进输入框
你自己的资料：/zebra new ~/cases/<名字>  把病历放进它的 records/ 文件夹再提问
全部功能：/zebra
EOF
