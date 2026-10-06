"""zebra trials: clinical trials from ClinicalTrials.gov (API v2) by condition, keyword, country and status.

Each trial carries its eligibility criteria, who may enrol, who to write to and
any reason to distrust the record (see zebra.sources.ctgov). Criteria are cut to
ctgov.ELIGIBILITY_BUDGET characters unless `--full-eligibility`, which sends the whole text and warns
when the answer passes FULL_TEXT_BUDGET characters, because the mod cuts a tool
result at 60,000.

Relevance (CP1-8): ClinicalTrials.gov expands a condition through MeSH, so a
search for spinal muscular atrophy in China returned androgen-receptor prostate
and breast cancer trials. Trials that name the condition nowhere (conditions,
title, keywords, MeSH terms) are set aside in `result.filtered` with a warning
naming them; `--keep-unrelated` turns the check off.

ChiCTR (chictr.org.cn), where most Chinese investigator-initiated trials are
registered, cannot be searched by a script: checked 2026-10-06, its search page
(searchproj.html) and record pages (showproj.html) answer every scripted request
with HTTP 405 and the Alibaba Cloud WAF block page ("your request has been
blocked"), also with a browser User-Agent, the site's own acw_tc cookie and a
Referer, while the home page loads the vendor's anti-bot script; only a real
browser gets through. WHO ICTRP (trialsearch.who.int) mirrors ChiCTR and its
search form can be posted by a script (676 records for "spinal muscular atrophy",
including ChiCTR2600132926 registered 2026-09-20), but its robots.txt is
"User-agent: * / Disallow: /", so zebra does not query it. The note saying so is
printed for every query and raised to a warning when the country asked for is
China, Hong Kong or Taiwan.
"""

from __future__ import annotations

import argparse
import json
from typing import Any, Dict, List, Optional

from zebra.core import Outcome, UsageError
from zebra.sources import ctgov

FULL_TEXT_BUDGET = 50_000  # characters; past this --full-eligibility warns and asks for a smaller --limit
CHICTR_NOTE = ("Chinese trials are also registered in ChiCTR (chictr.org.cn), which is not searched here: its search "
               "and record pages answer scripted requests with HTTP 405 from the site's Alibaba Cloud WAF (checked "
               "2026-10-06; only a real browser gets through). A trial registered only in ChiCTR will not appear in "
               "this answer — search it by hand at https://www.chictr.org.cn/searchproj.html, or in the WHO ICTRP "
               "portal (https://trialsearch.who.int/), which mirrors ChiCTR records but whose robots.txt disallows "
               "automated access, so zebra does not query it either")
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
    keep = bool(getattr(args, "keep_unrelated", False))
    out = ctgov.search(condition, term=term, country=country, status=status, limit=args.limit,
                       full_eligibility=full, filter_relevance=not keep)
    r = out.result
    head = f"ClinicalTrials.gov: {condition}" + (f" + {term}" if term else "")
    head += f" · status {status}" + (f" · sites in {r['country']}" if r["country"] else "")
    head += f" — {r['total']} trials" + (" (before the relevance check)" if r.get("filtered") else "") \
        + f", showing {r['returned']}"
    lines: List[str] = [head]
    lines += [_line(s, r["country"]) for s in r["studies"]]
    if not r["studies"]:
        if r.get("total_any_status"):
            lines.append(f"no {status} trials; {r['total_any_status']} in any status (--status ANY)")
        else:
            lines.append("no trials")
    if r.get("filtered"):
        lines.append(f"set aside as unrelated ({len(r['filtered'])}; they name '{condition}' nowhere in their "
                     "conditions, title, keywords or MeSH terms; --keep-unrelated shows them):")
        for f in r["filtered"][:10]:
            lines.append(f"    {f['nct_id']} · {f['title']} · conditions: {', '.join(f['conditions'][:3]) or '-'}")
    elif r.get("relevance_check", "").startswith("not applied"):
        lines.append(f"note: relevance {r['relevance_check']}")
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
                 "full_eligibility": full, "keep_unrelated": keep}
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("trials", help="clinical trials (ClinicalTrials.gov v2): by condition, keyword, country, status")
    p.add_argument("condition", nargs="+", help='condition, e.g. "Dravet syndrome"')
    p.add_argument("--term", help="extra keyword: gene, drug, modality")
    p.add_argument("--country", help="country with a trial site, as ClinicalTrials.gov names it (China, United States, ...)")
    p.add_argument("--status", default="RECRUITING", type=str.upper,
                   help="RECRUITING (default), NOT_YET_RECRUITING, ACTIVE_NOT_RECRUITING, COMPLETED, ... or ANY")
    p.add_argument("--limit", type=int, default=20)
    p.add_argument("--keep-unrelated", action="store_true",
                   help="do not set aside trials that name the condition nowhere (CT.gov's MeSH expansion)")
    p.add_argument("--full-eligibility", action="store_true",
                   help=f"send the whole eligibility criteria instead of the first {ctgov.ELIGIBILITY_BUDGET} "
                        "characters; use a small --limit")
    p.set_defaults(func=_trials)
