"""Case workspace: one folder per patient (or family, or cohort question).

    <case>/case.json              structured profile (schema zebra.case/1)
    <case>/evidence/ledger.jsonl  every source a zebra command used, E1, E2, ...
    <case>/records/               the person's own files; never uploaded
    <case>/reports/               what gets written for clinicians and family

Everything stays on this machine. `privacy.identifiers` lists strings (name,
date of birth, record numbers) that the mod refuses to send anywhere.
"""

from __future__ import annotations

import json
import os
import re
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

from zebra.http import now_iso

SCHEMA = "zebra.case/1"
HPO_RE = re.compile(r"^HP:\d{7}$")
ROLES = ("family", "patient", "clinician", "researcher")
PHENO_STATUS = ("present", "excluded")
HYP_STATUS = ("leading", "considered", "excluded", "confirmed")
LEAD_KINDS = ("approved", "trial", "repurposing", "n-of-1", "supportive", "other")


class CaseError(Exception):
    pass


def _path(case_dir: str) -> Path:
    return Path(os.path.expanduser(case_dir)).resolve()


def case_file(case_dir: str) -> Path:
    return _path(case_dir) / "case.json"


def init(case_dir: str, title: str = "", role: str = "family", language: str = "zh") -> Dict[str, Any]:
    root = _path(case_dir)
    if (root / "case.json").exists():
        raise CaseError(f"{root} already holds a case")
    if role not in ROLES:
        raise CaseError(f"role must be one of {', '.join(ROLES)}")
    for sub in ("evidence", "records", "reports"):
        (root / sub).mkdir(parents=True, exist_ok=True)
    stamp = now_iso()
    data: Dict[str, Any] = {
        "schema": SCHEMA,
        "id": root.name,
        "title": title or root.name,
        "role": role,
        "language": language,
        "created_at": stamp,
        "updated_at": stamp,
        "proband": {"sex": "unknown", "age": None, "ancestry": None, "consanguinity": None},
        "phenotypes": [],
        "variants": [],
        "family": {"members": [], "notes": ""},
        "hypotheses": [],
        "therapy_leads": [],
        "questions": [],
        "timeline": [],
        "privacy": {"identifiers": []},
    }
    _write(root, data)
    (root / "evidence" / "ledger.jsonl").touch()
    gitignore = root / ".gitignore"
    if not gitignore.exists():
        gitignore.write_text("# patient data: keep out of any repository\n*\n", "utf-8")
    return data


def load(case_dir: str) -> Dict[str, Any]:
    path = case_file(case_dir)
    if not path.exists():
        raise CaseError(f"no case at {path.parent} (run: zebra case init {case_dir})")
    try:
        data = json.loads(path.read_text("utf-8"))
    except ValueError as err:
        raise CaseError(f"{path} is not valid JSON: {err}") from None
    if data.get("schema") != SCHEMA:
        raise CaseError(f"{path}: unknown schema {data.get('schema')!r}")
    return data


def _write(root: Path, data: Dict[str, Any]) -> None:
    data["updated_at"] = now_iso()
    tmp = root / "case.json.tmp"
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", "utf-8")
    os.replace(tmp, root / "case.json")


@contextmanager
def _locked(path: Path) -> Iterator[None]:
    path.parent.mkdir(parents=True, exist_ok=True)
    lock = open(str(path) + ".lock", "a+")
    try:
        try:
            import fcntl

            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        except ImportError:  # Windows: best effort
            pass
        yield
    finally:
        lock.close()


@contextmanager
def editing(case_dir: str) -> Iterator[Dict[str, Any]]:
    """Load, let the caller change, write back under a lock."""
    root = _path(case_dir)
    with _locked(root / "case.json"):
        data = load(case_dir)
        yield data
        _write(root, data)


def _next_id(items: List[Dict[str, Any]], prefix: str) -> str:
    taken = {str(i.get("id")) for i in items}
    n = 1
    while f"{prefix}{n}" in taken:
        n += 1
    return f"{prefix}{n}"


def add_phenotype(case_dir: str, hpo_id: str, label: str, status: str = "present", onset: Optional[str] = None,
                  source: Optional[str] = None, note: Optional[str] = None) -> Dict[str, Any]:
    if not HPO_RE.match(hpo_id):
        raise CaseError(f"{hpo_id!r} is not an HPO id (HP:0000000)")
    if status not in PHENO_STATUS:
        raise CaseError(f"status must be one of {', '.join(PHENO_STATUS)}")
    entry = {"id": hpo_id, "label": label, "status": status, "onset": onset, "source": source, "note": note}
    with editing(case_dir) as data:
        data["phenotypes"] = [p for p in data["phenotypes"] if p["id"] != hpo_id] + [entry]
    return entry


def add_variant(case_dir: str, **fields: Any) -> Dict[str, Any]:
    if not any(fields.get(k) for k in ("hgvs_c", "hgvs_g", "vcf", "hgvs_p", "description")):
        raise CaseError("a variant needs at least one of hgvs_c, hgvs_g, vcf (chr-pos-ref-alt), hgvs_p or description")
    with editing(case_dir) as data:
        entry = {"id": _next_id(data["variants"], "v")}
        entry.update({k: v for k, v in fields.items() if v is not None})
        data["variants"].append(entry)
    return entry


def add_hypothesis(case_dir: str, disease: str, status: str = "considered", ids: Optional[Dict[str, str]] = None,
                   support: Optional[List[str]] = None, against: Optional[List[str]] = None,
                   note: Optional[str] = None) -> Dict[str, Any]:
    if status not in HYP_STATUS:
        raise CaseError(f"status must be one of {', '.join(HYP_STATUS)}")
    with editing(case_dir) as data:
        existing = next((h for h in data["hypotheses"] if h["disease"].lower() == disease.lower()), None)
        if existing:
            existing.update({"status": status})
            if ids:
                existing.setdefault("ids", {}).update(ids)
            if support:
                existing["support"] = sorted(set(existing.get("support", [])) | set(support))
            if against:
                existing["against"] = sorted(set(existing.get("against", [])) | set(against))
            if note:
                existing["note"] = note
            return existing
        entry = {"id": _next_id(data["hypotheses"], "h"), "disease": disease, "status": status, "ids": ids or {},
                 "support": support or [], "against": against or [], "note": note}
        data["hypotheses"].append(entry)
    return entry


def add_lead(case_dir: str, name: str, kind: str, status: Optional[str] = None,
             evidence: Optional[List[str]] = None, note: Optional[str] = None) -> Dict[str, Any]:
    if kind not in LEAD_KINDS:
        raise CaseError(f"kind must be one of {', '.join(LEAD_KINDS)}")
    with editing(case_dir) as data:
        entry = {"id": _next_id(data["therapy_leads"], "t"), "name": name, "kind": kind, "status": status,
                 "evidence": evidence or [], "note": note}
        data["therapy_leads"].append(entry)
    return entry


def set_acmg(case_dir: str, variant_id: str, classification: str, points: Optional[int], codes: List[str],
             note: Optional[str] = None) -> Dict[str, Any]:
    """Store zebra's research-grade ACMG reading on a recorded variant."""
    with editing(case_dir) as data:
        v = next((x for x in data["variants"] if str(x.get("id")) == variant_id), None)
        if v is None:
            raise CaseError(f"no variant {variant_id} in the case")
        v["zebra_acmg"] = {"classification": classification, "points": points, "codes": codes, "note": note,
                           "at": now_iso(), "grade": "research"}
        return v


def add_question(case_dir: str, text: str) -> List[str]:
    with editing(case_dir) as data:
        if text not in data["questions"]:
            data["questions"].append(text)
        return list(data["questions"])


def set_profile(case_dir: str, **fields: Any) -> Dict[str, Any]:
    """Title, role, language and proband basics (sex, age, ancestry, consanguinity)."""
    with editing(case_dir) as data:
        for key in ("title", "role", "language"):
            value = fields.get(key)
            if value is None:
                continue
            if key == "role" and value not in ROLES:
                raise CaseError(f"role must be one of {', '.join(ROLES)}")
            data[key] = value
        for key in ("sex", "age", "ancestry", "consanguinity"):
            if fields.get(key) is not None:
                data["proband"][key] = fields[key]
        return {"title": data["title"], "role": data["role"], "language": data["language"], "proband": dict(data["proband"])}


def set_identifiers(case_dir: str, identifiers: List[str]) -> List[str]:
    clean = sorted({s.strip() for s in identifiers if s and len(s.strip()) >= 2})
    with editing(case_dir) as data:
        data["privacy"]["identifiers"] = clean
    return clean


def remove(case_dir: str, kind: str, item_id: str) -> bool:
    key = {"phenotype": "phenotypes", "variant": "variants", "hypothesis": "hypotheses", "lead": "therapy_leads"}.get(kind)
    if not key:
        raise CaseError("kind must be phenotype, variant, hypothesis or lead")
    with editing(case_dir) as data:
        before = len(data[key])
        data[key] = [i for i in data[key] if str(i.get("id")) != item_id]
        return len(data[key]) < before


def append_ledger(case_dir: str, command: str, query: Dict[str, Any], sources: List[Dict[str, Any]]) -> List[str]:
    """Append one ledger row per source; returns their evidence ids (E1, E2, ...)."""
    root = _path(case_dir)
    if not (root / "case.json").exists():
        raise CaseError(f"no case at {root}")
    ledger = root / "evidence" / "ledger.jsonl"
    ids: List[str] = []
    with _locked(ledger):
        count = 0
        if ledger.exists():
            with open(ledger, "r", encoding="utf-8") as fh:
                count = sum(1 for line in fh if line.strip())
        with open(ledger, "a", encoding="utf-8") as fh:
            for src in sources:
                count += 1
                eid = f"E{count}"
                ids.append(eid)
                fh.write(json.dumps({"eid": eid, "command": command, "query": query, **src}, ensure_ascii=False) + "\n")
    return ids


def read_ledger(case_dir: str) -> List[Dict[str, Any]]:
    ledger = _path(case_dir) / "evidence" / "ledger.jsonl"
    if not ledger.exists():
        return []
    rows = []
    with open(ledger, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                try:
                    rows.append(json.loads(line))
                except ValueError:
                    continue
    return rows


def summary(case_dir: str) -> Dict[str, Any]:
    """What the mod's case board draws."""
    data = load(case_dir)
    ledger = read_ledger(case_dir)
    order = {s: i for i, s in enumerate(("confirmed", "leading", "considered", "excluded"))}
    hyps = sorted(data["hypotheses"], key=lambda h: order.get(h.get("status"), 9))
    return {
        "path": str(_path(case_dir)),
        "id": data["id"],
        "title": data["title"],
        "role": data["role"],
        "language": data.get("language"),
        "updated_at": data["updated_at"],
        "phenotypes": [
            {"id": p["id"], "label": p.get("label"), "status": p.get("status", "present")} for p in data["phenotypes"]
        ],
        "variants": [
            {
                "id": v["id"],
                "gene": v.get("gene"),
                "label": v.get("hgvs_c") or v.get("hgvs_g") or v.get("vcf") or v.get("hgvs_p") or v.get("description"),
                "zygosity": v.get("zygosity"),
                "classification": (v.get("zebra_acmg") or {}).get("classification") or v.get("classification_lab"),
            }
            for v in data["variants"]
        ],
        "hypotheses": [
            {"id": h["id"], "disease": h["disease"], "status": h["status"], "support": len(h.get("support", [])),
             "against": len(h.get("against", []))}
            for h in hyps
        ],
        "therapy_leads": [{"id": t["id"], "name": t["name"], "kind": t["kind"], "status": t.get("status")} for t in data["therapy_leads"]],
        "questions": data["questions"],
        "evidence_count": len(ledger),
        "identifiers": len(data["privacy"].get("identifiers", [])),
    }
