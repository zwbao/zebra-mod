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
# HGVS shape: an optional reference sequence (or a gene symbol), then a sequence
# type and a description. The reference is optional because that is how reports
# and families write it ("p.Phe508del"); which transcript was used is then not
# recorded, and `hgvs_reference` says so, so the caller can warn instead of
# refusing the commonest form. Accepts NM_004006.3:c.6439-?_7309+?del.
HGVS_REF = (r"(?:N[CGMRPTW]_\d+(?:\.\d+)?|ENS[TGP]\d+(?:\.\d+)?|LRG_\d+(?:t\d+|p\d+)?"
            r"|[A-Za-z0-9][A-Za-z0-9_.-]{0,20})(?:\([A-Za-z0-9_.-]+\))?")
# the description must at least start the way HGVS descriptions do: a position
# (c.1521del, g.117559590A>G, p.Phe508del) or the *-/? forms reports use
HGVS_DESC = r"(?:[*-]?\d|\(|\[|[A-Z][a-z]{2}\d|[A-Z]\d|=|\?)\S*"
HGVS_RE = re.compile(r"^(?:(?P<ref>" + HGVS_REF + r"):)?(?P<kind>[cgmnopr])\.(?P<desc>" + HGVS_DESC + r")$")


class CaseError(Exception):
    pass


_BARE_ID_RE = re.compile(r"^(NCT|PMC|rs)(\d+)$", re.I)


def check_id(pair: str) -> str:
    """Validate an identifier's shape; returns the normalised `PREFIX:VALUE` form.

    The forms people actually paste are accepted: `NCT04006210` and `rs113993960`
    carry their prefix without a colon, which is how every registry prints them.
    """
    text = str(pair).strip()
    bare = _BARE_ID_RE.match(text)
    if bare:
        prefix = bare.group(1).upper()
        text = f"{'DBSNP' if prefix == 'RS' else prefix}:{bare.group(0) if prefix == 'RS' else bare.group(2)}"
    if ":" not in text:
        raise CaseError(f"{text!r} is not an identifier: it takes PREFIX:VALUE (ORPHA:33069, OMIM:607208, "
                        "MONDO:0100135, PMID:28919360) or a registry id that carries its own prefix "
                        "(NCT04006210, rs113993960)")
    prefix, value = text.split(":", 1)
    prefix, value = prefix.strip().upper(), value.strip()
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
                        "c.1521_1523del, NC_000007.14:g.117559590_117559592del or p.Phe508del)")
    return value


def hgvs_reference(text: str) -> Optional[str]:
    """The reference sequence an HGVS string names, or None when it names none."""
    got = HGVS_RE.match(str(text).strip())
    return got.group("ref") if got else None


def variant_notes(fields: Dict[str, Any]) -> List[str]:
    """What is worth saying about a recorded variant without refusing it."""
    notes = []
    for key in ("hgvs_c", "hgvs_g", "hgvs_p"):
        value = fields.get(key)
        if value and not hgvs_reference(value):
            notes.append(f"{key} {value!r} names no reference sequence: the transcript or genome the laboratory "
                         "used is not recorded with it, so the position cannot be checked")
    return notes


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
    # the nested objects summary() and the report reach into
    for v in data["variants"]:
        if v.get("zebra_acmg") is not None and not isinstance(v["zebra_acmg"], dict):
            raise CaseError(f"{where}: variant {v.get('id')!r}: zebra_acmg must be an object "
                            f"(classification/points/codes), not {type(v['zebra_acmg']).__name__}")
    for h in data["hypotheses"]:
        for field in ("support", "against"):
            if h.get(field) is not None and not isinstance(h[field], list):
                raise CaseError(f"{where}: hypothesis {h.get('id')!r}: {field!r} must be a list of evidence ids "
                                f"(E1, E7), not {type(h[field]).__name__}")
        if h.get("ids") is not None and not isinstance(h["ids"], dict):
            raise CaseError(f"{where}: hypothesis {h.get('id')!r}: 'ids' must be an object of PREFIX: VALUE")
    for t in data["therapy_leads"]:
        if t.get("evidence") is not None and not isinstance(t["evidence"], list):
            raise CaseError(f"{where}: lead {t.get('id')!r}: 'evidence' must be a list of evidence ids")
    if data.get("_issued_ids") is not None and not isinstance(data["_issued_ids"], dict):
        data["_issued_ids"] = {}
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
    tmp = root / f"case.json.{os.getpid()}.tmp"
    try:
        tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", "utf-8")
        os.replace(tmp, root / "case.json")
    except OSError as err:
        raise CaseError(f"could not write {root / 'case.json'}: {err.strerror or err}. The case was not changed; "
                        "check the folder's permissions and that the disk is not full.") from None


@contextmanager
def _locked(path: Path) -> Iterator[None]:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        lock = open(str(path) + ".lock", "a+")
    except OSError as err:
        raise CaseError(f"could not take the lock on {path}: {err.strerror or err}. Nothing was changed; check the "
                        "folder's permissions.") from None
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
    if not (root / "case.json").exists():  # never create a folder for a mistyped path
        raise CaseError(f"no case at {root} (run: zebra case init {case_dir})")
    with _locked(root / "case.json"):
        data = load(case_dir)
        yield data
        _write(root, data)


def _next_id(data: Dict[str, Any], key: str, prefix: str) -> str:
    """The next id for `data[key]`, numbered from the highest ever issued.

    A deleted item must not let its id be re-used: `v1` quoted from an earlier
    `case summary` would otherwise attach an ACMG reading to a different gene.
    """
    issued = data.setdefault("_issued_ids", {})
    highest = 0
    for item in data.get(key) or []:
        text = str(item.get("id") or "")
        if text.startswith(prefix) and text[len(prefix):].isdigit():
            highest = max(highest, int(text[len(prefix):]))
    try:
        highest = max(highest, int(issued.get(prefix, 0)))
    except (TypeError, ValueError):
        pass
    issued[prefix] = highest + 1
    return f"{prefix}{highest + 1}"


# Each mutator comes in two halves: `apply_*` changes a loaded case dict and is
# what `case apply` calls (all ops under one lock, one write), and the public
# function wraps it in `editing()` for single-op callers.

def apply_phenotype(data: Dict[str, Any], hpo_id: str, label: str, status: Optional[str] = None,
                    onset: Optional[str] = None, source: Optional[str] = None,
                    note: Optional[str] = None) -> Dict[str, Any]:
    """Add or update one phenotype.

    Fields that are not given are kept from the existing entry: re-recording a
    term must never turn a clinician's `excluded` into `present`, nor drop the
    note that said why it was excluded.
    """
    if not HPO_RE.match(str(hpo_id)):
        raise CaseError(f"{hpo_id!r} is not an HPO id (HP:0000000)")
    if status is not None and status not in PHENO_STATUS:
        raise CaseError(f"status must be one of {', '.join(PHENO_STATUS)}")
    existing = next((p for p in data["phenotypes"] if str(p.get("id")) == hpo_id), None) or {}
    entry = {"id": hpo_id, "label": label or existing.get("label"),
             "status": status or existing.get("status") or "present",
             "onset": onset if onset is not None else existing.get("onset"),
             "source": source if source is not None else existing.get("source"),
             "note": note if note is not None else existing.get("note")}
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
    for key in ("region", "iscn", "exons", "motif", "method", "note", "source", "description", "gene"):
        if fields.get(key) is not None and not isinstance(fields[key], str):
            raise CaseError(f"{key} must be a string, not {type(fields[key]).__name__}")
    if fields.get("region"):
        got = re.match(r"^(?:chr)?([0-9]{1,2}|X|Y|MT)[:\s](\d+)[-_](\d+)$", fields["region"].strip(), re.I)
        if not got:
            raise CaseError(f"region {fields['region']!r} must look like 15:23000000-28500000")
        if int(got.group(3)) < int(got.group(2)):
            raise CaseError(f"region {fields['region']!r} ends before it starts")
    return fields


def apply_variant(data: Dict[str, Any], **fields: Any) -> Dict[str, Any]:
    fields = check_variant(**fields)
    entry = {"id": _next_id(data, "variants", "v")}
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
    entry = {"id": _next_id(data, "hypotheses", "h"), "disease": disease, "status": status or "considered",
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
    entry = {"id": _next_id(data, "therapy_leads", "t"), "name": name.strip(), "kind": kind, "status": status,
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


def _clean_identifiers(values: List[str]) -> List[str]:
    return sorted({s.strip() for s in values if isinstance(s, str) and len(s.strip()) >= 2})


def set_identifiers(case_dir: str, identifiers: List[str]) -> List[str]:
    """Replace the protected identifiers."""
    clean = _clean_identifiers(identifiers)
    with editing(case_dir) as data:
        data["privacy"]["identifiers"] = clean
    return clean


def add_identifiers(case_dir: str, identifiers: List[str]) -> List[str]:
    """Add protected identifiers. Read and write under one lock, so parallel adds never lose one."""
    with editing(case_dir) as data:
        merged = _clean_identifiers(list(data["privacy"].get("identifiers") or []) + list(identifiers))
        data["privacy"]["identifiers"] = merged
    return merged


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
    mark = root / "evidence" / ".high_water"
    ids: List[str] = []
    with _locked(ledger):
        # Number from the highest id ever issued, not from the line count: a
        # deleted line must not let an id be re-used, and a row whose eid
        # cannot be read still counts.
        highest = 0
        if mark.exists():  # survives deleting the last row, which the file alone cannot
            try:
                highest = max(highest, int(mark.read_text("utf-8").strip() or 0))
            except (ValueError, OSError):
                pass
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
                fh.write(json.dumps({"eid": eid, "command": command, "query": query, **src},
                                     ensure_ascii=False, allow_nan=False) + "\n")
        try:
            mark.write_text(str(highest) + "\n", "utf-8")
        except OSError:
            pass  # the ids in the file are still monotonic; only a deletion could repeat one
    return ids


def ledger_ids(case_dir: str) -> Set[str]:
    """Every evidence id in the ledger: what a hypothesis's support/against may cite."""
    return {str(r.get("eid")) for r in read_ledger(case_dir) if isinstance(r, dict) and r.get("eid")}


def read_ledger(case_dir: str) -> List[Dict[str, Any]]:
    """Every readable ledger row. A line that is not a JSON object is skipped, not raised."""
    ledger = _path(case_dir) / "evidence" / "ledger.jsonl"
    if not ledger.exists():
        return []
    rows = []
    with open(ledger, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                row = json.loads(line)
            except ValueError:
                continue
            if isinstance(row, dict):  # `null`, a list or a number is not a row
                rows.append(row)
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
             "support": len(h.get("support") or []), "against": len(h.get("against") or []),
             "ids": dict(h.get("ids") or {})}
            for h in hyps
        ],
        "therapy_leads": [{"id": t.get("id"), "name": t.get("name"), "kind": t.get("kind"), "status": t.get("status")}
                          for t in data["therapy_leads"]],
        "questions": data["questions"],
        "evidence_count": len(ledger),
        "identifiers": len(data["privacy"].get("identifiers", [])),
    }
