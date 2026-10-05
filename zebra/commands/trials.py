"""zebra trials: clinical trials from ClinicalTrials.gov (API v2) by condition, keyword, country and status.

Each trial carries its eligibility criteria, who may enrol, who to write to and
any reason to distrust the record (see zebra.sources.ctgov). Criteria are cut to
ctgov.ELIGIBILITY_BUDGET characters unless `--full-eligibility`, which sends the whole text and warns
when the answer passes FULL_TEXT_BUDGET characters, because the mod cuts a tool
result at 60,000. ChiCTR, where Chinese trials are also registered, is not
searched here and never was; the note saying so is printed for every query and
raised to a warning when the country asked for is China, Hong Kong or Taiwan.
"""

from __future__ import annotations

import argparse
import json
from typing import Any, Dict, List, Optional

from zebra.core import Outcome, UsageError
from zebra.sources import ctgov

FULL_TEXT_BUDGET = 50_000  # characters; past this --full-eligibility warns and asks for a smaller --limit
CHICTR_NOTE = ("Chinese trials are also registered in ChiCTR (chictr.org.cn), which is not searched here: "
               "a trial registered only in ChiCTR will not appear in this answer")
CHINESE_REGISTRIES = ("China", "Hong Kong", "Taiwan")


def _contact_str(c: Dict[str, Any]) -> str:
    name = c.get("name") or "?"
    role = f" ({(c.get('role') or '').lower()})" if c.get("role") else ""
    rest = " · ".join(x for x in (c.get("phone"), c.get("email")) if x)
    return f"{name}{role}" + (f" · {rest}" if rest else "")


def _site_str(x: Dict[str, Any]) -> str:
    s = f"{x.get('facility') or '?'} ({x.get('city') or '?'}, {x.get('status') or '?'})"
    cs = x.get("contacts") or []
    if cs:
        s += " [" + "; ".join(_contact_str(c) for c in cs) + "]"
    return s


def _line(s: Dict[str, Any], country: Optional[str] = None) -> str:
    phase = "/".join(p.replace("PHASE", "P") for p in s.get("phases") or []) or (s.get("study_type") or "-").lower()
    iv = ", ".join(f"{i['name']} ({(i.get('type') or '').lower()})" for i in (s.get("interventions") or [])[:4]) or "-"
    age = f"{s.get('min_age') or '?'}–{s.get('max_age') or 'no max'}"
    ctry = ", ".join(s.get("countries") or []) or "no sites listed"
    if len(s.get("countries") or []) > 6:
        ctry = ", ".join(s["countries"][:6]) + f" +{len(s['countries']) - 6}"
    line = (f"{s['nct_id']} · {phase} · {s.get('status')} · {s.get('title')}\n"
            f"    {iv} · n={s.get('enrollment') or '?'} · age {age} · start {s.get('start_date') or '?'}"
            f" · primary completion {s.get('primary_completion_date') or '?'}\n"
            f"    {s.get('n_sites')} sites: {ctry} · sponsor {s.get('sponsor') or '-'} · {s.get('url')}")
    for flag in s.get("status_flags") or []:
        text = flag
        if flag.startswith("site status conflicts"):
            text += f": a site still reads recruiting, but the overall status ({s.get('status')}) is the one to trust"
        line += f"\n    ! {text}"
    who: List[str] = []
    if s.get("sex"):
        who.append(f"sex {str(s['sex']).lower()}")
    if s.get("healthy_volunteers") is not None:
        who.append("healthy volunteers accepted" if s["healthy_volunteers"] else "patients only")
    if who:
        line += "\n    " + " · ".join(who)
    if s.get("eligibility"):
        tail = ""
        if s.get("eligibility_truncated"):
            tail = f" [cut from {s.get('eligibility_chars')} characters; --full-eligibility for all of it]"
        line += f"\n    eligibility: {s['eligibility']}{tail}"
    elif "eligibility" in s:
        line += "\n    eligibility: none in the response for this trial"
    if s.get("contacts"):
        line += "\n    contact: " + "; ".join(_contact_str(c) for c in s["contacts"])
    if country is not None and s.get("sites") is not None:
        sites = "; ".join(_site_str(x) for x in s["sites"])
        more = f" +{s['sites_in_country'] - len(s['sites'])}" if s["sites_in_country"] > len(s["sites"]) else ""
        line += f"\n    in {country}: {s['sites_in_country']} site(s): {sites}{more}"
    return line


def _trials(args: argparse.Namespace) -> Outcome:
    condition = " ".join(args.condition).strip()
    if not condition:
        raise UsageError("give a condition, e.g. zebra trials \"Dravet syndrome\"")
    if args.limit < 1:
        raise UsageError("--limit must be at least 1")
    status = args.status.upper()
    if status not in ctgov.STATUSES:
        raise UsageError(f"--status must be one of {', '.join(ctgov.STATUSES)}")
    term = (args.term or "").strip() or None
    country = (args.country or "").strip() or None
    full = bool(getattr(args, "full_eligibility", False))
    out = ctgov.search(condition, term=term, country=country, status=status, limit=args.limit,
                       full_eligibility=full)
    r = out.result
    head = f"ClinicalTrials.gov: {condition}" + (f" + {term}" if term else "")
    head += f" · status {status}" + (f" · sites in {r['country']}" if r["country"] else "")
    head += f" — {r['total']} trials, showing {r['returned']}"
    lines: List[str] = [head]
    lines += [_line(s, r["country"]) for s in r["studies"]]
    if not r["studies"]:
        if r.get("total_any_status"):
            lines.append(f"no {status} trials; {r['total_any_status']} in any status (--status ANY)")
        else:
            lines.append("no trials")
    lines.append("note: " + CHICTR_NOTE)
    if r["country"] in CHINESE_REGISTRIES:
        out.warnings.append(CHICTR_NOTE + f"; the {r['country']} sites listed here are only the ones "
                            "whose trial was also registered on ClinicalTrials.gov")
    if r["country"] == "China":
        lines.append("note: ClinicalTrials.gov lists Hong Kong and Taiwan as separate countries")
    out.text = "\n".join(lines)
    if full:
        rendered = len(out.text)
        payload = len(json.dumps(r, ensure_ascii=False))
        worst = max(rendered, payload)
        if worst > FULL_TEXT_BUDGET:
            out.warnings.append(
                f"--full-eligibility produced {rendered} characters of text and {payload} of JSON for "
                f"{r['returned']} trials, over the {FULL_TEXT_BUDGET}-character guide; a tool result is cut at "
                f"60,000, so lower --limit (about {max(1, int(r['returned'] * FULL_TEXT_BUDGET / worst))} trials "
                "fits) or drop --full-eligibility")
    out.query = {"condition": condition, "term": term, "country": country, "status": status, "limit": args.limit,
                 "full_eligibility": full}
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("trials", help="clinical trials (ClinicalTrials.gov v2): by condition, keyword, country, status")
    p.add_argument("condition", nargs="+", help='condition, e.g. "Dravet syndrome"')
    p.add_argument("--term", help="extra keyword: gene, drug, modality")
    p.add_argument("--country", help="country with a trial site, as ClinicalTrials.gov names it (China, United States, ...)")
    p.add_argument("--status", default="RECRUITING", type=str.upper,
                   help="RECRUITING (default), NOT_YET_RECRUITING, ACTIVE_NOT_RECRUITING, COMPLETED, ... or ANY")
    p.add_argument("--limit", type=int, default=20)
    p.add_argument("--full-eligibility", action="store_true",
                   help=f"send the whole eligibility criteria instead of the first {ctgov.ELIGIBILITY_BUDGET} "
                        "characters; use a small --limit")
    p.set_defaults(func=_trials)
