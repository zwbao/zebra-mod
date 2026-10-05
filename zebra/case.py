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
from typing import Any, Dict, Iterator, List, Optional, Set

from zebra.http import now_iso

SCHEMA = "zebra.case/1"
HPO_RE = re.compile(r"^HP:\d{7}$")
ROLES = ("family", "patient", "clinician", "researcher")
PHENO_STATUS = ("present", "excluded")
HYP_STATUS = ("leading", "considered", "excluded", "confirmed")
LEAD_KINDS = ("approved", "trial", "repurposing", "n-of-1", "supportive", "other")
EVIDENCE_RE = re.compile(r"^E\d+$")
_EID_IN_TEXT = re.compile(r'"eid"\s*:\s*"E(\d+)"')

# Variant kinds. "small" (SNV/indel, the only kind schema zebra.case/1 could
# hold) stays the default, so a case.json written before these existed reads
# unchanged; the report forms Chinese families usually hold get their own kinds.
VARIANT_KINDS = ("small", "cnv", "exon_cnv", "copy_number", "repeat_expansion")
CNV_TYPES = ("gain", "loss", "duplication", "deletion", "amplification", "unknown")

# Identifier shapes, so an id that cannot exist is refused at the door. Values
# are checked for shape only: existence is a source lookup, not a regex.
ID_SHAPES = {
    "HP": (re.compile(r"^\d{7}$"), "HP:0001250"),
    "OMIM": (re.compile(r"^\d{6}(\.\d{4})?$"), "OMIM:607208"),
    "ORPHA": (re.compile(r"^\d{1,7}$"), "ORPHA:33069"),
    "MONDO": (re.compile(r"^\d{7}$"), "MONDO:0100135"),
    "MEDGEN": (re.compile(r"^(C\d{6,7}|CN\d{6})$"), "MEDGEN:C1843367"),
    "DOID": (re.compile(r"^\d+$"), "DOID:0050434"),
    "MESH": (re.compile(r"^[CD]\d{6}$"), "MESH:D004831"),
    "NCT": (re.compile(r"^\d{8}$"), "NCT04006210"),
    "PMID": (re.compile(r"^\d{1,9}$"), "PMID:28919360"),
    "PMC": (re.compile(r"^\d+$"), "PMC:5760072"),
    "HGNC": (re.compile(r"^\d+$"), "HGNC:10585"),
    "CA": (re.compile(r"^\d+$"), "CA:116077"),
    "CLINVAR": (re.compile(r"^(VCV|RCV|SCV)\d{9}(\.\d+)?$"), "CLINVAR:VCV000038596"),
    "DBSNP": (re.compile(r"^rs\d+$", re.I), "DBSNP:rs113993960"),
    "ENSEMBL": (re.compile(r"^ENS[A-Z]*\d{6,}(\.\d+)?$"), "ENSEMBL:ENSG00000198947"),
}
# HGVS: a versioned reference sequence (or a bare gene symbol for c./p.) and a
# sequence type. Accepts the exon-boundary form NM_004006.3:c.6439-?_7309+?del.
HGVS_RE = re.compile(
    r"^(?:(?:N[CGMRPTW]_\d+(?:\.\d+)?|ENS[TGP]\d+(?:\.\d+)?|LRG_\d+(?:t\d+|p\d+)?|[A-Za-z0-9][A-Za-z0-9_.-]{0,20})"
    r"(?:\([A-Za-z0-9_.-]+\))?):[cgmnopr]\.\S+$"
)


class CaseError(Exception):
    pass


def check_id(pair: str) -> str:
    """Validate a `PREFIX:VALUE` identifier's shape; returns the normalised form."""
    text = str(pair).strip()
    if ":" not in text:
        raise CaseError(f"{text!r} is not an identifier: it takes PREFIX:VALUE (ORPHA:33069, OMIM:607208, "
                        "MONDO:0100135, NCT04006210, PMID:28919360)")
    prefix, value = text.split(":", 1)
    prefix, value = prefix.strip().upper(), value.strip()
    if prefix == "NCT" and value == "":  # NCT04006210 written without a colon
        raise CaseError("NCT ids are written NCT04006210")
    shape = ID_SHAPES.get(prefix)
    if shape is None:
        raise CaseError(f"unknown identifier prefix {prefix!r} in {text!r}; zebra validates "
                        + ", ".join(sorted(ID_SHAPES)))
    if not shape[0].match(value):
        raise CaseError(f"{text!r} is not a well-formed {prefix} id (expected like {shape[1]})")
    return f"{prefix}:{value}"


def check_hgvs(text: str, field: str = "hgvs") -> str:
    value = str(text).strip()
    if not HGVS_RE.match(value):
        raise CaseError(f"{field}: {value!r} is not HGVS (expected like NM_000492.4:c.1521_1523del, "
                        "NC_000007.14:g.117559590_117559592del or NP_000483.3:p.Phe508del)")
    return value


def str_list(value: Any, field: str) -> List[str]:
    """A list of non-empty strings. A bare string is the commonest shape error and is refused by name."""
    if value is None:
        return []
    if isinstance(value, str):
        raise CaseError(f"{field} must be a list of strings, not the string {value!r} "
                        f'(write ["{value}"]; a string would be stored one character per item)')
    if not isinstance(value, (list, tuple)):
        raise CaseError(f"{field} must be a list of strings, not {type(value).__name__}")
    out = []
    for item in value:
        if not isinstance(item, str) or not item.strip():
            raise CaseError(f"{field}: every item must be a non-empty string, got {item!r}")
        out.append(item.strip())
    return out


def check_evidence_ids(value: Any, field: str, known: Optional[Set[str]] = None) -> List[str]:
    """Ledger ids (E1, E7). When `known` is given, every id must already be in the ledger."""
    ids = str_list(value, field)
    for eid in ids:
        if not EVIDENCE_RE.match(eid):
            raise CaseError(f"{field}: {eid!r} is not an evidence id (ledger ids look like E7; "
                            "run `zebra case ledger` to see them)")
        if known is not None and eid not in known:
            raise CaseError(f"{field}: {eid} is not in this case's evidence ledger"
                            + (f" (it holds {min(known, key=_eid_num)}..{max(known, key=_eid_num)})" if known
                               else " (the ledger is empty: run a zebra command with --case first)"))
    return ids


def _eid_num(eid: str) -> int:
    try:
        return int(str(eid)[1:])
    except (ValueError, TypeError):
        return 0


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


_LIST_KEYS = ("phenotypes", "variants", "hypotheses", "therapy_leads", "questions", "timeline")
_DICT_KEYS = {"proband": {"sex": "unknown", "age": None, "ancestry": None, "consanguinity": None},
              "family": {"members": [], "notes": ""},
              "privacy": {"identifiers": []}}


def _validate(data: Any, where: str) -> Dict[str, Any]:
    """Structural check of a case.json, which a person may have edited by hand.

    Missing containers are filled with their defaults (an older or partial file
    still loads); a container of the wrong type is an error naming the key,
    rather than a KeyError or TypeError deep inside `summary`.
    """
    if not isinstance(data, dict):
        raise CaseError(f"{where}: the case must be a JSON object, not {type(data).__name__}")
    for key in _LIST_KEYS:
        value = data.get(key)
        if value is None:
            data[key] = []
        elif not isinstance(value, list):
            raise CaseError(f"{where}: {key!r} must be a list, not {type(value).__name__}")
    for key, default in _DICT_KEYS.items():
        value = data.get(key)
        if value is None:
            data[key] = dict(default)
        elif not isinstance(value, dict):
            raise CaseError(f"{where}: {key!r} must be an object, not {type(value).__name__}")
        else:
            for sub, sub_default in default.items():
                value.setdefault(sub, sub_default)
    if not isinstance(data["privacy"].get("identifiers"), list):
        raise CaseError(f"{where}: privacy.identifiers must be a list")
    for key, default in (("title", ""), ("role", "family"), ("id", ""), ("language", "zh"), ("updated_at", "")):
        if not isinstance(data.get(key), str):
            data[key] = str(data.get(key) or default)
    for key in _LIST_KEYS:
        if key == "questions":
            data[key] = [str(q) for q in data[key]]
        elif key != "timeline":
            bad = [i for i in data[key] if not isinstance(i, dict)]
            if bad:
                raise CaseError(f"{where}: every item in {key!r} must be an object, got {bad[0]!r}")
    return data


def load(case_dir: str) -> Dict[str, Any]:
    path = case_file(case_dir)
    if not path.exists():
        raise CaseError(f"no case at {path.parent} (run: zebra case init {case_dir})")
    try:
        data = json.loads(path.read_text("utf-8"))
    except ValueError as err:
        raise CaseError(f"{path} is not valid JSON: {err}") from None
    schema = data.get("schema") if isinstance(data, dict) else None
    if schema != SCHEMA:
        raise CaseError(f"{path}: unknown schema {schema!r}")
    return _validate(data, str(path))


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


# Each mutator comes in two halves: `apply_*` changes a loaded case dict and is
# what `case apply` calls (all ops under one lock, one write), and the public
# function wraps it in `editing()` for single-op callers.

def apply_phenotype(data: Dict[str, Any], hpo_id: str, label: str, status: str = "present",
                    onset: Optional[str] = None, source: Optional[str] = None,
                    note: Optional[str] = None) -> Dict[str, Any]:
    if not HPO_RE.match(str(hpo_id)):
        raise CaseError(f"{hpo_id!r} is not an HPO id (HP:0000000)")
    if status not in PHENO_STATUS:
        raise CaseError(f"status must be one of {', '.join(PHENO_STATUS)}")
    entry = {"id": hpo_id, "label": label, "status": status, "onset": onset, "source": source, "note": note}
    data["phenotypes"] = [p for p in data["phenotypes"] if str(p.get("id")) != hpo_id] + [entry]
    return entry


def add_phenotype(case_dir: str, hpo_id: str, label: str, status: str = "present", onset: Optional[str] = None,
                  source: Optional[str] = None, note: Optional[str] = None) -> Dict[str, Any]:
    with editing(case_dir) as data:
        return apply_phenotype(data, hpo_id, label, status=status, onset=onset, source=source, note=note)


# Per kind: the fields that identify the finding (at least one is required).
_VARIANT_IDENTITY = {
    "small": ("hgvs_c", "hgvs_g", "vcf", "hgvs_p", "description"),
    "cnv": ("region", "iscn", "description"),
    "exon_cnv": ("exons", "hgvs_c", "description"),
    "copy_number": ("copy_number", "description"),
    "repeat_expansion": ("repeat_count", "description"),
}
VARIANT_FIELDS = ("gene", "hgvs_c", "hgvs_g", "hgvs_p", "vcf", "assembly", "zygosity", "inheritance",
                  "classification_lab", "source", "description", "kind", "region", "iscn", "cnv_type",
                  "copy_number", "exons", "genes", "motif", "repeat_count", "method", "note")


def check_variant(**fields: Any) -> Dict[str, Any]:
    """Validate a variant record of any kind, without writing. Returns the normalised fields."""
    kind = fields.get("kind") or "small"
    if kind not in VARIANT_KINDS:
        raise CaseError(f"variant kind must be one of {', '.join(VARIANT_KINDS)}, got {kind!r}")
    fields["kind"] = kind
    needed = _VARIANT_IDENTITY[kind]
    if not any(fields.get(k) not in (None, "", [], {}) for k in needed):
        raise CaseError(f"a {kind} variant needs at least one of " + ", ".join(needed))
    if fields.get("cnv_type") and fields["cnv_type"] not in CNV_TYPES:
        raise CaseError(f"cnv_type must be one of {', '.join(CNV_TYPES)}")
    for key in ("hgvs_c", "hgvs_g", "hgvs_p"):
        if fields.get(key):
            fields[key] = check_hgvs(fields[key], key)
    if fields.get("copy_number") is not None:
        try:
            fields["copy_number"] = int(fields["copy_number"])
        except (TypeError, ValueError):
            raise CaseError(f"copy_number must be a whole number, got {fields['copy_number']!r}") from None
        if fields["copy_number"] < 0:
            raise CaseError("copy_number cannot be negative")
    if fields.get("repeat_count") is not None and not isinstance(fields["repeat_count"], (int, float, str)):
        raise CaseError("repeat_count must be a number or the range as the report writes it")
    if fields.get("genes") is not None:
        fields["genes"] = str_list(fields["genes"], "genes")
    if fields.get("region") is not None and not isinstance(fields["region"], str):
        raise CaseError("region must be a string like 15:23000000-28500000")
    return fields


def apply_variant(data: Dict[str, Any], **fields: Any) -> Dict[str, Any]:
    fields = check_variant(**fields)
    entry = {"id": _next_id(data["variants"], "v")}
    entry.update({k: v for k, v in fields.items() if v is not None})
    data["variants"].append(entry)
    return entry


def add_variant(case_dir: str, **fields: Any) -> Dict[str, Any]:
    with editing(case_dir) as data:
        return apply_variant(data, **fields)


def apply_hypothesis(data: Dict[str, Any], disease: str, status: Optional[str] = None,
                     ids: Optional[Dict[str, str]] = None, support: Optional[List[str]] = None,
                     against: Optional[List[str]] = None, note: Optional[str] = None) -> Dict[str, Any]:
    """Add or update a hypothesis. `status=None` keeps the status an existing hypothesis already has."""
    if status is not None and status not in HYP_STATUS:
        raise CaseError(f"status must be one of {', '.join(HYP_STATUS)}")
    if not isinstance(disease, str) or not disease.strip():
        raise CaseError(f"disease must be a non-empty string, got {disease!r}")
    disease = disease.strip()
    existing = next((h for h in data["hypotheses"] if str(h.get("disease", "")).lower() == disease.lower()), None)
    if existing:
        if status is not None:  # an update without a status never demotes a leading/confirmed hypothesis
            existing["status"] = status
        existing.setdefault("status", "considered")
        if ids:
            existing.setdefault("ids", {}).update(ids)
        if support:
            existing["support"] = sorted(set(existing.get("support") or []) | set(support), key=_eid_num)
        if against:
            existing["against"] = sorted(set(existing.get("against") or []) | set(against), key=_eid_num)
        if note:
            existing["note"] = note
        return existing
    entry = {"id": _next_id(data["hypotheses"], "h"), "disease": disease, "status": status or "considered",
             "ids": ids or {}, "support": support or [], "against": against or [], "note": note}
    data["hypotheses"].append(entry)
    return entry


def add_hypothesis(case_dir: str, disease: str, status: Optional[str] = None, ids: Optional[Dict[str, str]] = None,
                   support: Optional[List[str]] = None, against: Optional[List[str]] = None,
                   note: Optional[str] = None) -> Dict[str, Any]:
    with editing(case_dir) as data:
        return apply_hypothesis(data, disease, status=status, ids=ids, support=support, against=against, note=note)


def apply_lead(data: Dict[str, Any], name: str, kind: str, status: Optional[str] = None,
               evidence: Optional[List[str]] = None, note: Optional[str] = None,
               ids: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    if kind not in LEAD_KINDS:
        raise CaseError(f"kind must be one of {', '.join(LEAD_KINDS)}")
    if not isinstance(name, str) or not name.strip():
        raise CaseError(f"a therapy lead needs a name (a non-empty string), got {name!r}")
    entry = {"id": _next_id(data["therapy_leads"], "t"), "name": name.strip(), "kind": kind, "status": status,
             "evidence": evidence or [], "note": note}
    if ids:
        entry["ids"] = ids
    data["therapy_leads"].append(entry)
    return entry


def add_lead(case_dir: str, name: str, kind: str, status: Optional[str] = None,
             evidence: Optional[List[str]] = None, note: Optional[str] = None,
             ids: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    with editing(case_dir) as data:
        return apply_lead(data, name, kind, status=status, evidence=evidence, note=note, ids=ids)


def apply_acmg(data: Dict[str, Any], variant_id: str, classification: str, points: Optional[int],
               codes: List[str], note: Optional[str] = None) -> Dict[str, Any]:
    v = next((x for x in data["variants"] if str(x.get("id")) == variant_id), None)
    if v is None:
        have = ", ".join(str(x.get("id")) for x in data["variants"]) or "none recorded"
        raise CaseError(f"no variant {variant_id} in the case (recorded: {have})")
    v["zebra_acmg"] = {"classification": classification, "points": points, "codes": codes, "note": note,
                       "at": now_iso(), "grade": "research"}
    return v


def set_acmg(case_dir: str, variant_id: str, classification: str, points: Optional[int], codes: List[str],
             note: Optional[str] = None) -> Dict[str, Any]:
    """Store zebra's research-grade ACMG reading on a recorded variant."""
    with editing(case_dir) as data:
        return apply_acmg(data, variant_id, classification, points, codes, note=note)


def apply_question(data: Dict[str, Any], text: str) -> List[str]:
    if not isinstance(text, str) or not text.strip():
        raise CaseError(f"a question must be a non-empty string, got {text!r}")
    if text not in data["questions"]:
        data["questions"].append(text)
    return list(data["questions"])


def add_question(case_dir: str, text: str) -> List[str]:
    with editing(case_dir) as data:
        return apply_question(data, text)


PROFILE_FIELDS = ("title", "role", "language", "sex", "age", "ancestry", "consanguinity")


def apply_profile(data: Dict[str, Any], **fields: Any) -> Dict[str, Any]:
    """Title, role, language and proband basics (sex, age, ancestry, consanguinity)."""
    for key in ("title", "role", "language"):
        value = fields.get(key)
        if value is None:
            continue
        if not isinstance(value, str) or not value.strip():
            raise CaseError(f"profile.{key} must be a non-empty string, got {value!r}")
        if key == "role" and value not in ROLES:
            raise CaseError(f"role must be one of {', '.join(ROLES)}")
        data[key] = value
    for key in ("sex", "age", "ancestry", "consanguinity"):
        if fields.get(key) is not None:
            data["proband"][key] = fields[key]
    return {"title": data["title"], "role": data["role"], "language": data["language"],
            "proband": dict(data["proband"])}


def set_profile(case_dir: str, **fields: Any) -> Dict[str, Any]:
    with editing(case_dir) as data:
        return apply_profile(data, **fields)


def set_identifiers(case_dir: str, identifiers: List[str]) -> List[str]:
    clean = sorted({s.strip() for s in identifiers if s and len(s.strip()) >= 2})
    with editing(case_dir) as data:
        data["privacy"]["identifiers"] = clean
    return clean


REMOVE_KINDS = {"phenotype": "phenotypes", "variant": "variants", "hypothesis": "hypotheses",
                "lead": "therapy_leads"}


def apply_remove(data: Dict[str, Any], kind: str, item_id: str) -> bool:
    key = REMOVE_KINDS.get(kind)
    if not key:
        raise CaseError("kind must be phenotype, variant, hypothesis or lead")
    if not isinstance(item_id, str) or not item_id.strip():
        raise CaseError(f"remove needs the item's id as a string (v1, h2, t1, or an HPO id), got {item_id!r}")
    before = len(data[key])
    data[key] = [i for i in data[key] if str(i.get("id")) != item_id]
    return len(data[key]) < before


def remove(case_dir: str, kind: str, item_id: str) -> bool:
    with editing(case_dir) as data:
        return apply_remove(data, kind, item_id)


def append_ledger(case_dir: str, command: str, query: Dict[str, Any], sources: List[Dict[str, Any]]) -> List[str]:
    """Append one ledger row per source; returns their evidence ids (E1, E2, ...)."""
    root = _path(case_dir)
    if not (root / "case.json").exists():
        raise CaseError(f"no case at {root}")
    ledger = root / "evidence" / "ledger.jsonl"
    ids: List[str] = []
    with _locked(ledger):
        # Number from the highest id ever issued, not from the line count: a
        # deleted line must not let an id be re-used, and a row whose eid
        # cannot be read still counts.
        highest = 0
        unterminated = False
        if ledger.exists():
            with open(ledger, "r", encoding="utf-8") as fh:
                last = ""
                for line in fh:
                    last = line
                    text = line.strip()
                    if not text:
                        continue
                    try:
                        highest = max(highest, _eid_num(json.loads(text).get("eid", "")))
                    except ValueError:
                        # a truncated row: recover its id from the raw text
                        found = _EID_IN_TEXT.search(text)
                        if found:
                            highest = max(highest, int(found.group(1)))
                unterminated = bool(last) and not last.endswith("\n")
        with open(ledger, "a", encoding="utf-8") as fh:
            if unterminated:
                # a crash mid-write left a partial line; close it so the next
                # row is not concatenated onto it and silently lost
                fh.write("\n")
            for src in sources:
                highest += 1
                eid = f"E{highest}"
                ids.append(eid)
                fh.write(json.dumps({"eid": eid, "command": command, "query": query, **src}, ensure_ascii=False) + "\n")
    return ids


def ledger_ids(case_dir: str) -> Set[str]:
    """Every evidence id in the ledger: what a hypothesis's support/against may cite."""
    return {str(r.get("eid")) for r in read_ledger(case_dir) if r.get("eid")}


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

    def label(v: Dict[str, Any]) -> Optional[str]:
        """How the finding reads on one line, whatever kind it is."""
        kind = v.get("kind") or "small"
        if kind == "cnv":
            size = v.get("region") or v.get("iscn")
            parts = [str(x) for x in (v.get("cnv_type"), size) if x]
            if v.get("copy_number") is not None:
                parts.append(f"CN{v['copy_number']}")
            return " ".join(parts) or v.get("description")
        if kind == "exon_cnv":
            parts = [x for x in (v.get("gene"), (f"exon {v['exons']}" if v.get("exons") else None),
                                 v.get("cnv_type"), v.get("hgvs_c")) if x]
            return " ".join(str(p) for p in parts) or v.get("description")
        if kind == "copy_number":
            return " ".join(str(x) for x in (v.get("gene"), f"copy number {v.get('copy_number')}") if x)
        if kind == "repeat_expansion":
            return " ".join(str(x) for x in (v.get("gene"), v.get("motif"),
                                             f"{v.get('repeat_count')} repeats") if x)
        return v.get("hgvs_c") or v.get("hgvs_g") or v.get("vcf") or v.get("hgvs_p") or v.get("description")
    return {
        "path": str(_path(case_dir)),
        "id": data["id"],
        "title": data["title"],
        "role": data["role"],
        "language": data.get("language"),
        "updated_at": data["updated_at"],
        "phenotypes": [
            {"id": p.get("id"), "label": p.get("label"), "status": p.get("status", "present")} for p in data["phenotypes"]
        ],
        "variants": [
            {
                "id": v.get("id"),
                "kind": v.get("kind") or "small",
                "gene": v.get("gene"),
                "label": label(v),
                "zygosity": v.get("zygosity"),
                "classification": (v.get("zebra_acmg") or {}).get("classification") or v.get("classification_lab"),
            }
            for v in data["variants"]
        ],
        "hypotheses": [
            {"id": h.get("id"), "disease": str(h.get("disease", "")), "status": h.get("status", "considered"),
             "support": len(h.get("support") or []), "against": len(h.get("against") or [])}
            for h in hyps
        ],
        "therapy_leads": [{"id": t.get("id"), "name": t.get("name"), "kind": t.get("kind"), "status": t.get("status")}
                          for t in data["therapy_leads"]],
        "questions": data["questions"],
        "evidence_count": len(ledger),
        "identifiers": len(data["privacy"].get("identifiers", [])),
    }
