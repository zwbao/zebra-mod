"""zebra trials: clinical trials from ClinicalTrials.gov (API v2) by condition, keyword, country and status."""

from __future__ import annotations

import argparse
from typing import Any, Dict, List, Optional

from zebra.core import Outcome, UsageError
from zebra.sources import ctgov


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
    if country is not None and s.get("sites") is not None:
        sites = "; ".join(f"{x.get('facility') or '?'} ({x.get('city') or '?'}, {x.get('status') or '?'})" for x in s["sites"])
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
    out = ctgov.search(condition, term=term, country=country, status=status, limit=args.limit)
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
    if r["country"] in ("China",):
        lines.append("note: ClinicalTrials.gov lists Hong Kong and Taiwan as separate countries; Chinese trials are "
                     "also registered in ChiCTR (not searched here)")
    out.text = "\n".join(lines)
    out.query = {"condition": condition, "term": term, "country": country, "status": status, "limit": args.limit}
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("trials", help="clinical trials (ClinicalTrials.gov v2): by condition, keyword, country, status")
    p.add_argument("condition", nargs="+", help='condition, e.g. "Dravet syndrome"')
    p.add_argument("--term", help="extra keyword: gene, drug, modality")
    p.add_argument("--country", help="country with a trial site, as ClinicalTrials.gov names it (China, United States, ...)")
    p.add_argument("--status", default="RECRUITING", type=str.upper,
                   help="RECRUITING (default), NOT_YET_RECRUITING, ACTIVE_NOT_RECRUITING, COMPLETED, ... or ANY")
    p.add_argument("--limit", type=int, default=20)
    p.set_defaults(func=_trials)
