"""What changed in the world since a case was last checked (CP1-12).

A diagnostic odyssey lasts years; the answer to "is this VUS pathogenic yet?",
"is there a trial now?" or "what was published?" changes underneath a case
that nobody touches. `recheck` asks the same questions again and compares:

- each recorded sequence variant: ClinVar classification and review stars,
  gnomAD maximum group frequency (the variant card, without the slow
  Chinese-cohort and MaveDB lookups);
- each gene in the case: ClinGen gene-disease validity rows;
- each leading or considered hypothesis: recruiting trials (ClinicalTrials.gov)
  and papers first published since the last check (Europe PMC).

The answers are kept in `<case>/evidence/recheck.json` (this run and the one
before it), so the next run can say what moved. Nothing is concluded here: a
ClinVar change is a reason to re-run the ACMG reading, a new trial a reason
to read its eligibility criteria. Every query carries biology only and is
refused if it would contain one of the case's protected identifiers.
"""

from __future__ import annotations

import json
import os
import re
from concurrent.futures import ThreadPoolExecutor
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from zebra.core import Outcome

SNAPSHOT = "recheck.json"
MAX_VARIANTS = 6
MAX_HYPOTHESES = 5
MAX_GENES = 8
FIRST_LOOKBACK_DAYS = 365


def _today() -> date:
    return datetime.now(timezone.utc).date()


def _snapshot_path(case_dir: str) -> Path:
    return Path(case_dir).expanduser() / "evidence" / SNAPSHOT


def load_snapshot(case_dir: str) -> Optional[Dict[str, Any]]:
    p = _snapshot_path(case_dir)
    if not p.exists():
        return None
    try:
        data = json.loads(p.read_text("utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) and isinstance(data.get("items"), dict) else None


def save_snapshot(case_dir: str, current: Dict[str, Any], previous: Optional[Dict[str, Any]]) -> None:
    p = _snapshot_path(case_dir)
    p.parent.mkdir(parents=True, exist_ok=True)
    body = dict(current)
    if previous:
        body["previous"] = {k: v for k, v in previous.items() if k != "previous"}
    tmp = p.with_name(f".{p.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(body, ensure_ascii=False, indent=1, allow_nan=False), "utf-8")
    os.replace(tmp, p)


# ------------------------------------------------------------------ what to ask

def _variant_query(v: Dict[str, Any]) -> Optional[str]:
    """The sequence-variant notation a card can read, or None (CNVs, repeats and copy numbers are not rechecked)."""
    if (v.get("kind") or "small") != "small":
        return None
    for key in ("vcf", "hgvs_g", "hgvs_c"):
        val = v.get(key)
        if isinstance(val, str) and val.strip():
            return val.strip()
    return None


def plan(data: Dict[str, Any]) -> Dict[str, List[Dict[str, Any]]]:
    """What a recheck of this case would ask, before anything is sent."""
    variants = []
    for v in data.get("variants") or []:
        q = _variant_query(v)
        if q:
            variants.append({"key": f"variant:{v.get('id')}", "id": v.get("id"), "query": q,
                             "assembly": v.get("assembly") or "GRCh38", "gene": v.get("gene")})
    genes = sorted({str(v.get("gene")).strip().upper() for v in data.get("variants") or []
                    if isinstance(v.get("gene"), str) and re.fullmatch(r"[A-Za-z0-9-]{2,20}", v["gene"].strip())})
    hyps = [h for h in data.get("hypotheses") or [] if h.get("status") in ("leading", "considered", "confirmed")
            and isinstance(h.get("disease"), str) and h["disease"].strip()]
    order = {"confirmed": 0, "leading": 1, "considered": 2}
    hyps.sort(key=lambda h: order.get(h.get("status"), 9))
    return {
        "variants": variants[:MAX_VARIANTS],
        "genes": [{"key": f"gene:{g}", "symbol": g} for g in genes[:MAX_GENES]],
        "hypotheses": [{"key": f"disease:{h.get('id')}", "id": h.get("id"), "disease": h["disease"].strip()}
                       for h in hyps[:MAX_HYPOTHESES]],
        "left_out": {
            "variants": max(0, len(variants) - MAX_VARIANTS),
            "genes": max(0, len(genes) - MAX_GENES),
            "hypotheses": max(0, len(hyps) - MAX_HYPOTHESES),
        },
    }


# ------------------------------------------------------------------ asking

def _ask_variant(item: Dict[str, Any]) -> Tuple[Dict[str, Any], Outcome]:
    from zebra.sources import variant

    out = variant.card(item["query"], assembly=item["assembly"], gene=item.get("gene"), contracts=False)
    r = out.result or {}
    clin = r.get("clinvar") or {}
    pop = r.get("acmg_inputs") or {}
    # the card turns a failed source into a warning and an empty field: that is "not checked",
    # never "no ClinVar record" or "absent from gnomAD" (a transient error must not read as a change)
    clinvar_failed = any(w.startswith("ClinVar") for w in out.warnings)
    gnomad_failed = any(w.startswith("gnomAD") or "gnomAD unavailable" in w for w in out.warnings)
    state = {
        "clinvar": {k: clin.get(k) for k in ("vcv", "classification", "stars", "last_evaluated", "url")} if clin else None,
        "clinvar_checked": not clinvar_failed,
        "grpmax_af": pop.get("grpmax_af"),
        "gnomad_ac": pop.get("gnomad_ac"),
        # a frequency from VEP's fallback is another dataset, not comparable with gnomAD's own answer
        "freq_checked": not gnomad_failed and not str(pop.get("frequency_source") or "").startswith("VEP"),
    }
    return state, out


def _ask_gene(item: Dict[str, Any]) -> Tuple[Dict[str, Any], Outcome]:
    from zebra.sources import clingen

    out = clingen.validity(symbol=item["symbol"])
    rows = sorted(({"disease": r.get("disease"), "mondo": r.get("mondo"), "classification": r.get("classification"),
                    "moi": r.get("moi"), "date": r.get("date"), "url": r.get("url")} for r in out.result or []),
                  key=lambda r: (str(r.get("mondo")), str(r.get("disease")), str(r.get("moi"))))
    return {"validity": rows}, out


def _ask_trials(item: Dict[str, Any]) -> Tuple[Dict[str, Any], Outcome]:
    from zebra.sources import ctgov

    out = ctgov.search(item["disease"], status="RECRUITING", limit=100)
    r = out.result or {}
    studies = r.get("studies") or []
    total = r.get("total")
    seen = len(studies) + len(r.get("filtered") or [])
    # a list ranked by relevance and cut at a page is a sliding window: only a complete list can say
    # which trials are new or no longer recruiting
    complete = isinstance(total, int) and total <= seen
    return {"recruiting": sorted({s.get("nct_id") for s in studies if s.get("nct_id")}),
            "titles": {s.get("nct_id"): s.get("title") for s in studies if s.get("nct_id")},
            "total": total, "complete": complete}, out


def _ask_papers(item: Dict[str, Any], since: date, until: date) -> Tuple[Dict[str, Any], Outcome]:
    """Papers first published after `since` (exclusive) up to `until`: windows never overlap or leave a gap."""
    from zebra.sources import europepmc

    start = since + timedelta(days=1)
    if start > until:  # checked again the same day: an empty window, nothing to ask
        return {"since": since.isoformat(), "checked_on": until.isoformat(), "count": 0, "latest": []}, Outcome({})
    q = f'"{item["disease"]}" AND FIRST_PDATE:[{start.isoformat()} TO {until.isoformat()}]'
    out = europepmc.search(q, limit=8, sort="date")
    r = out.result or {}
    hits = [{"pmid": h.get("pmid"), "title": h.get("title"), "year": h.get("year"), "url": h.get("url")}
            for h in r.get("hits") or []]
    return {"since": since.isoformat(), "checked_on": until.isoformat(), "count": r.get("hitCount"), "latest": hits}, out


def _guarded(identifiers: List[str]) -> Callable[[str], Optional[str]]:
    """A check that refuses any outgoing query holding one of the case's protected identifiers."""
    from zebra.report_export import identifier_hits

    def check(text: str) -> Optional[str]:
        hit, _ = identifier_hits(text, identifiers)
        return hit
    return check


def run(case_dir: str, data: Dict[str, Any], workers: int = 4) -> Tuple[Outcome, Dict[str, Any], Optional[Dict[str, Any]]]:
    """(the answer, the snapshot to save, the snapshot it was compared with)."""
    previous = load_snapshot(case_dir)
    today = _today()
    since = None
    if previous and previous.get("checked_on"):
        try:
            since = date.fromisoformat(str(previous["checked_on"]))
        except ValueError:
            since = None
    first_run = since is None
    since = since or (today - timedelta(days=FIRST_LOOKBACK_DAYS))
    todo = plan(data)
    unsafe = _guarded(list((data.get("privacy") or {}).get("identifiers") or []))

    before = (previous or {}).get("items") or {}

    def paper_since(key: str) -> date:
        """The day this question was last answered, or the first-run lookback."""
        old = before.get(key) or {}
        try:
            return date.fromisoformat(str(old.get("checked_on")))
        except ValueError:
            return since

    jobs: List[Tuple[str, str, Callable[[], Tuple[Dict[str, Any], Outcome]]]] = []
    refused: List[str] = []
    for v in todo["variants"]:
        if unsafe(v["query"]) or (v.get("gene") and unsafe(v["gene"])):
            refused.append(v["key"])
            continue
        jobs.append((v["key"], f"variant {v['id']} ({v['query']})", lambda v=v: _ask_variant(v)))
    for g in todo["genes"]:
        jobs.append((g["key"], f"ClinGen validity for {g['symbol']}", lambda g=g: _ask_gene(g)))
    for h in todo["hypotheses"]:
        if unsafe(h["disease"]):
            refused.append(h["key"])
            continue
        jobs.append((f"trials:{h['id']}", f"recruiting trials for {h['disease']}", lambda h=h: _ask_trials(h)))
        ps = paper_since(f"papers:{h['id']}")
        jobs.append((f"papers:{h['id']}", f"papers on {h['disease']} since {ps.isoformat()}",
                     lambda h=h, ps=ps: _ask_papers(h, ps, today)))

    results: Dict[str, Dict[str, Any]] = {}
    sources: List[Dict[str, Any]] = []
    warnings: List[str] = []
    labels = {key: label for key, label, _ in jobs}

    def one(job: Tuple[str, str, Callable[[], Tuple[Dict[str, Any], Outcome]]]) -> Tuple[str, Any, Any]:
        key, label, fn = job
        try:
            state, out = fn()
            return key, state, out
        except Exception as err:  # noqa: BLE001 - one failed question costs only its own line
            return key, None, f"{type(err).__name__}: {str(err)[:200]}"

    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        for key, state, out in pool.map(one, jobs):
            if state is None:
                warnings.append(f"not rechecked: {labels[key]} — {out}")
                continue
            results[key] = state
            sources.extend(out.sources)
            warnings.extend(f"{labels[key]}: {w}" for w in out.warnings)
    for key in refused:
        warnings.append(f"not rechecked: {key} — its query would carry a protected identifier of this case")

    changes = compare(before, results, labels, first_run)
    # The snapshot keeps the last known answer of every question: one failed source must not erase a
    # baseline (or the next recheck would compare against nothing and report no change).
    items: Dict[str, Any] = dict(before)
    for key, state in results.items():
        old = before.get(key) or {}
        if key.startswith("variant:"):
            if not state.get("clinvar_checked", True) and old:
                state = {**state, "clinvar": old.get("clinvar"), "clinvar_checked": old.get("clinvar_checked", True)}
            if not state.get("freq_checked", True) and old:
                state = {**state, "grpmax_af": old.get("grpmax_af"), "gnomad_ac": old.get("gnomad_ac"),
                         "freq_checked": old.get("freq_checked", True)}
        items[key] = state
    current = {"checked_on": today.isoformat(), "checked_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
               "items": items}
    unchecked = sorted(set(labels) - set(results)) + refused
    result = {
        "since": None if first_run else since.isoformat(),
        "checked_on": today.isoformat(),
        "first_check": first_run,
        "changes": changes,
        "unchanged": sorted(k for k in results if k in before and not any(c["key"] == k for c in changes)),
        "not_checked": unchecked,
        "left_out": todo["left_out"],
        "note": ("first recheck: this is the baseline the next one is compared with; papers are counted over the last "
                 f"{FIRST_LOOKBACK_DAYS} days") if first_run else
                "a change is a reason to look again (re-run the ACMG reading, read a trial's eligibility), not a conclusion",
    }
    return Outcome(result, sources=sources, warnings=warnings), current, previous


# ------------------------------------------------------------------ comparing

def _clin_text(c: Optional[Dict[str, Any]]) -> str:
    if not c:
        return "no ClinVar record"
    stars = c.get("stars")
    return f"{c.get('classification') or 'unclassified'}" + (f" ({stars} of 4 stars)" if stars is not None else "")


def compare(before: Dict[str, Any], now: Dict[str, Any], labels: Dict[str, str], first_run: bool) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for key, state in now.items():
        old = before.get(key)
        label = labels.get(key, key)
        if key.startswith("variant:"):
            if old is None:
                continue
            a, b = old.get("clinvar"), state.get("clinvar")
            clin_known = state.get("clinvar_checked", True) and old.get("clinvar_checked", True)
            if clin_known and ((a or {}).get("classification") != (b or {}).get("classification")
                               or (a or {}).get("stars") != (b or {}).get("stars")):
                out.append({"key": key, "what": "ClinVar", "about": label, "before": _clin_text(a), "after": _clin_text(b),
                            "url": (b or a or {}).get("url"), "act": "re-run the ACMG reading for this variant"})
            fa, fb = old.get("grpmax_af"), state.get("grpmax_af")
            freq_known = state.get("freq_checked", True) and old.get("freq_checked", True)
            if freq_known and ((fa is None) != (fb is None) or (fa and fb and abs(fb - fa) / max(fa, fb) > 0.2)):
                out.append({"key": key, "what": "gnomAD maximum group frequency", "about": label,
                            "before": fa, "after": fb, "act": "recheck PM2 / BS1 / BA1"})
        elif key.startswith("gene:"):
            if old is None:
                continue
            a = {(r.get("mondo") or r.get("disease"), r.get("moi")): r for r in old.get("validity") or []}
            b = {(r.get("mondo") or r.get("disease"), r.get("moi")): r for r in state.get("validity") or []}
            for k in sorted(set(a) | set(b), key=str):
                ra, rb = a.get(k), b.get(k)
                if (ra or {}).get("classification") != (rb or {}).get("classification"):
                    out.append({"key": key, "what": "ClinGen gene-disease validity", "about": label,
                                "disease": (rb or ra or {}).get("disease"),
                                "before": (ra or {}).get("classification") or "not curated",
                                "after": (rb or {}).get("classification") or "no longer listed",
                                "url": (rb or ra or {}).get("url"), "act": "a P/LP reading for this disease depends on it"})
        elif key.startswith("trials:"):
            if old is None:
                continue
            if not (state.get("complete", True) and old.get("complete", True)):
                if state.get("total") != old.get("total"):
                    out.append({"key": key, "what": "recruiting trials (count)", "about": label,
                                "before": old.get("total"), "after": state.get("total"),
                                "act": "the list is longer than one page, so which trials changed is not computed: search the trials"})
                continue
            new = sorted(set(state.get("recruiting") or []) - set(old.get("recruiting") or []))
            gone = sorted(set(old.get("recruiting") or []) - set(state.get("recruiting") or []))
            if new:
                out.append({"key": key, "what": "new recruiting trials", "about": label,
                            "trials": [{"nct_id": n, "title": (state.get("titles") or {}).get(n),
                                        "url": f"https://clinicaltrials.gov/study/{n}"} for n in new],
                            "act": "read each trial's eligibility against this patient"})
            if gone:
                out.append({"key": key, "what": "no longer recruiting", "about": label, "trials": gone,
                            "act": "update any lead that relied on them"})
        elif key.startswith("papers:"):
            if old is not None and (state.get("count") or 0) > 0:
                out.append({"key": key, "what": "new papers", "about": label, "count": state.get("count"),
                            "latest": state.get("latest"), "act": "screen the latest for this genotype and phenotype"})
    return out


def summary_line(changes: List[Dict[str, Any]]) -> str:
    if not changes:
        return "recheck: nothing changed"
    parts = []
    for c in changes[:4]:
        if c["what"] == "ClinVar":
            parts.append(f"ClinVar {c['before']} → {c['after']}")
        elif c["what"] == "new recruiting trials":
            parts.append(f"{len(c['trials'])} new trial(s)")
        elif c["what"] == "new papers":
            parts.append(f"{c['count']} new paper(s)")
        else:
            parts.append(c["what"])
    more = f" (+{len(changes) - 4} more)" if len(changes) > 4 else ""
    return "recheck: " + "; ".join(parts) + more
