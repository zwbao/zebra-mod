#!/usr/bin/env python3
"""Run each case through Claude Code with zebra-mod loaded, headless, and keep the transcript.

    python3 -I tools/bench/agent/run.py --cases cases.jsonl --out RUN_DIR [--only MME] [--limit 10]
        [--jobs 6] [--model claude-opus-5-5]

Each case runs in its own empty directory under RUN_DIR/cases/<id>/ with ZEBRA_CASE unset, so no
case is open and nothing is written to a real case. The plugin is loaded from this checkout
(--plugin-dir) and the user-scope install is disabled, so the version under test is the one in
plugin.json. Tools: zebra-mod's MCP tools and Skill only (no web search, no shell). The model is
given the HPO ids and labels and nothing else. Resumable: a case with a result.json is skipped.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

PROMPT = """I am a clinical geneticist. A patient with a suspected rare disease has these phenotypes (HPO):

Present:
{present}
{excluded}
Please give a ranked differential diagnosis of the 10 most likely rare diseases.

When you are done, end your reply with a fenced ```json block, exactly in this shape (10 entries, most likely first; "id" is an OMIM or ORPHA identifier as a tool returned it, or null if you have none):
{{"differential": [{{"rank": 1, "name": "...", "id": "OMIM:123456"}}]}}"""


def prompt_for(case: dict) -> str:
    fmt = lambda ts: "\n".join(f"- {t['id']} {t['label'] or ''}".rstrip() for t in ts)
    exc = f"\nExcluded (looked for and absent):\n{fmt(case['excluded'])}\n" if case["excluded"] else ""
    return PROMPT.format(present=fmt(case["present"]), excluded=exc)


def parse_final(text: str):
    import re
    blocks = re.findall(r"```json\s*(\{.*?\})\s*```", text or "", re.S)
    for b in reversed(blocks):
        try:
            d = json.loads(b)
            if isinstance(d.get("differential"), list):
                return d["differential"]
        except Exception:
            continue
    return None


PROVIDERS = {
    # Claude Code pointed at MiniMax's Anthropic-compatible endpoint; key from $MINIMAX_API_KEY.
    # Every model slot (main, small/fast, subagents) goes to the same MiniMax model.
    "minimax": {"ANTHROPIC_BASE_URL": "https://api.minimaxi.com/anthropic"},
}


def provider_env(provider: str, model: str) -> dict:
    if provider == "anthropic":
        return {}
    env = dict(PROVIDERS[provider])
    env["ANTHROPIC_AUTH_TOKEN"] = os.environ["MINIMAX_API_KEY"]
    for k in ("ANTHROPIC_MODEL", "ANTHROPIC_SMALL_FAST_MODEL", "ANTHROPIC_DEFAULT_OPUS_MODEL",
              "ANTHROPIC_DEFAULT_SONNET_MODEL", "ANTHROPIC_DEFAULT_HAIKU_MODEL", "CLAUDE_CODE_SUBAGENT_MODEL"):
        env[k] = model
    env["CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC"] = "1"
    return env


def run_one(case: dict, out: Path, model: str, timeout: int, provider: str = "anthropic") -> dict:
    d = out / "cases" / case["id"]
    if (d / "result.json").exists():
        return json.loads((d / "result.json").read_text("utf-8"))
    d.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("ZEBRA_CASE") and k != "ANTHROPIC_API_KEY"}
    env.update(provider_env(provider, model))
    cmd = [shutil.which("claude", path=str(Path.home() / ".local/bin")) or "claude", "-p", prompt_for(case),
           "--model", model, "--plugin-dir", str(ROOT),
           "--settings", json.dumps({"enabledPlugins": {"zebra-mod@zebra-mod": False}}),
           "--allowedTools", "mcp__zebra-mod__*", "Skill",
           "--disallowedTools", "WebSearch", "WebFetch", "Bash", "Agent",
           "--permission-mode", "dontAsk", "--max-turns", "40",
           "--no-session-persistence", "--output-format", "json"]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, cwd=d, env=env, capture_output=True, text=True, timeout=timeout)
        raw, err, code = p.stdout, p.stderr, p.returncode
    except subprocess.TimeoutExpired as e:
        raw, err, code = (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or ""), "timeout", -9
    wall = time.time() - t0
    (d / "transcript.json").write_text(raw or "", "utf-8")
    if err:
        (d / "stderr.txt").write_text(err, "utf-8")
    final, cost, turns, tools = "", None, None, []
    try:
        msgs = json.loads(raw)
        for m in msgs:
            if m.get("type") == "result":
                final, cost, turns = m.get("result") or "", m.get("total_cost_usd"), m.get("num_turns")
            if m.get("type") == "assistant":
                for c in m["message"].get("content", []):
                    if c.get("type") == "tool_use":
                        tools.append(c["name"])
    except Exception:
        pass
    res = {"id": case["id"], "dataset": case["dataset"], "differential": parse_final(final), "cost_usd": cost,
           "turns": turns, "wall_s": round(wall, 1), "exit": code, "tools": tools, "model": model}
    if res["differential"] is not None:  # failures stay retryable
        (d / "result.json").write_text(json.dumps(res, ensure_ascii=False), "utf-8")
    else:
        (d / "failed.json").write_text(json.dumps(res, ensure_ascii=False), "utf-8")
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--ids", nargs="*")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--model", default="claude-opus-5-5")
    ap.add_argument("--timeout", type=int, default=1200)
    ap.add_argument("--provider", choices=["anthropic", *PROVIDERS], default="anthropic")
    a = ap.parse_args()
    cases = [json.loads(l) for l in open(a.cases, encoding="utf-8")]
    if a.only:
        cases = [c for c in cases if c["dataset"] in a.only]
    if a.ids:
        cases = [c for c in cases if c["id"] in a.ids]
    if a.limit:
        cases = cases[: a.limit]
    a.out.mkdir(parents=True, exist_ok=True)
    meta = a.out / "run_meta.json"
    if not meta.exists():
        pj = json.loads((ROOT / ".claude-plugin/plugin.json").read_text("utf-8"))
        head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        cv = subprocess.run([str(Path.home() / ".local/bin/claude"), "--version"], capture_output=True, text=True).stdout.strip()
        meta.write_text(json.dumps({"zebra_mod": pj["version"], "commit": head, "claude_code": cv, "model": a.model, "provider": a.provider,
                                    "prompt": PROMPT, "started": time.strftime("%Y-%m-%dT%H:%M:%S%z")}, indent=1), "utf-8")
    done = cost = fail = 0
    t0 = time.time()
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(run_one, c, a.out, a.model, a.timeout, a.provider): c for c in cases}
        for f in as_completed(futs):
            r = f.result()
            done += 1
            cost += r.get("cost_usd") or 0
            fail += r["differential"] is None
            print(f"[{done}/{len(cases)}] {r['id']} turns={r['turns']} ${r.get('cost_usd') or 0:.2f} "
                  f"{r['wall_s']}s parsed={r['differential'] is not None} | total ${cost:.2f} fail={fail} "
                  f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
