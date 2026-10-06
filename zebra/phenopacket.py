"""GA4GH Phenopacket v2: read one into zebra's terms, and write a zebra case as one.

`read(path_or_dict)` returns what zebra needs from a phenopacket:

    {"id", "schema_version", "subject": {"id", "sex", "age"},
     "present":  [{"id", "label", "onset", "description"}],   # HPO terms observed
     "excluded": [{"id", "label", "onset", "description"}],   # HPO terms explicitly absent
     "diseases": [{"id", "label", "excluded", "from", "status", "onset"}],
     "genes":    [{"symbol", "id", "status"}],
     "variants": [{"kind", "gene", "gene_id", "hgvs_c", "hgvs_g", "hgvs_p", "vcf", "assembly",
                   "zygosity", "acmg", "status", "description", "structural_type"}],
     "references": ["PMID:..."], "warnings": [...]}

`diseases[].from` says where a disease came from: "interpretation" (a
diagnosis, with its `progressStatus` as `status`) or "diseases" (the
phenopacket's disease list, whose `excluded` flag is kept). Genes and variants
come from the genomic interpretations, each with its `interpretationStatus`.
Version 1 phenopackets are read too (`negated` features, top-level `genes` and
`variants`). Malformed input raises `PhenopacketError` (a ValueError) naming
what is wrong; a feature that is not an HPO term is skipped with a warning, not
silently.

`from_case(case_dict)` writes a zebra case (schema zebra.case/1) as a Phenopacket
v2 document. What leaves the case is deliberately narrow, because a phenopacket
is made to be shared: phenotypes (with `excluded` and onset), sex and age,
diagnoses and hypotheses that carry ontology ids, and the recorded variants as
structured fields. No free text the person typed is copied: not the case title,
notes, record file names, family, questions, a variant's description or method,
an onset that is not an age or an onset class; labels come from the HPO release
(or are the term id) rather than from the case. The phenopacket id is random
and the subject is "proband". Finally the whole document is checked against the
case's protected identifiers (and Chinese resident-ID shapes) with the same
check reports use (zebra.report_export.identifier_hits); a hit refuses the
export with PhenopacketError instead of writing it.

Onset: zebra stores it as text ("6 months", "出生时"). It becomes a Phenopacket
TimeElement only when it parses without guessing — an ISO 8601 duration, a
number with a unit, or one of the HPO onset classes below (ids checked against
the HPO 2026-09-01 release); anything else is kept word for word in the
feature's `description` as "onset: <text>", never converted to a made-up age.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from zebra import __version__

MAX_BYTES = 20 * 1024 * 1024  # a single phenopacket is kilobytes; refuse anything absurd
HPO_RE = re.compile(r"^HP[:_](\d{7})$")
ISO_DURATION = re.compile(r"^P(?=\d)(?:(\d+)Y)?(?:(\d+)M)?(?:(\d+)W)?(?:(\d+)D)?$")
SCHEMA_VERSION = "2.0"

# HPO onset classes (HP:0003674 Onset and below), each id verified in hp.json 2026-09-01
ONSET_CLASSES: Dict[str, Tuple[str, str]] = {
    "congenital": ("HP:0003577", "Congenital onset"),
    "at birth": ("HP:0003577", "Congenital onset"),
    "birth": ("HP:0003577", "Congenital onset"),
    "出生时": ("HP:0003577", "Congenital onset"),
    "出生": ("HP:0003577", "Congenital onset"),
    "先天": ("HP:0003577", "Congenital onset"),
    "先天性": ("HP:0003577", "Congenital onset"),
    "neonatal": ("HP:0003623", "Neonatal onset"),
    "新生儿期": ("HP:0003623", "Neonatal onset"),
    "新生儿": ("HP:0003623", "Neonatal onset"),
    "infantile": ("HP:0003593", "Infantile onset"),
    "infancy": ("HP:0003593", "Infantile onset"),
    "婴儿期": ("HP:0003593", "Infantile onset"),
    "childhood": ("HP:0011463", "Childhood onset"),
    "儿童期": ("HP:0011463", "Childhood onset"),
    "juvenile": ("HP:0003621", "Juvenile onset"),
    "青少年期": ("HP:0003621", "Juvenile onset"),
    "adult": ("HP:0003581", "Adult onset"),
    "成年期": ("HP:0003581", "Adult onset"),
    "成年": ("HP:0003581", "Adult onset"),
    "antenatal": ("HP:0030674", "Antenatal onset"),
    "prenatal": ("HP:0030674", "Antenatal onset"),
    "产前": ("HP:0030674", "Antenatal onset"),
    "胎儿期": ("HP:0030674", "Antenatal onset"),
}
_UNITS = {
    "y": "Y", "yr": "Y", "yrs": "Y", "year": "Y", "years": "Y", "岁": "Y", "周岁": "Y",
    "m": "M", "mo": "M", "mos": "M", "month": "M", "months": "M", "个月": "M", "月龄": "M",
    "w": "W", "wk": "W", "wks": "W", "week": "W", "weeks": "W", "周": "W",
    "d": "D", "day": "D", "days": "D", "天": "D",
}
_AGE_TEXT = re.compile(r"^(\d+(?:\.\d+)?)\s*(" + "|".join(sorted(map(re.escape, _UNITS), key=len, reverse=True))
                       + r")(?:\s*(?:old|of age|大))?$", re.I)

# Genotype Ontology / Sequence Ontology classes used below (verified in OLS4, 2026-10-06)
ZYGOSITY = {
    "heterozygous": ("GENO:0000135", "heterozygous"), "het": ("GENO:0000135", "heterozygous"),
    "杂合": ("GENO:0000135", "heterozygous"),
    "homozygous": ("GENO:0000136", "homozygous"), "hom": ("GENO:0000136", "homozygous"),
    "纯合": ("GENO:0000136", "homozygous"),
    "hemizygous": ("GENO:0000134", "hemizygous"), "hemi": ("GENO:0000134", "hemizygous"),
    "半合子": ("GENO:0000134", "hemizygous"),
    "compound heterozygous": ("GENO:0000402", "compound heterozygous"),
    "compound_heterozygous": ("GENO:0000402", "compound heterozygous"),
    "复合杂合": ("GENO:0000402", "compound heterozygous"),
}
_GENO_LABEL = {v[0]: v[1] for v in ZYGOSITY.values()}
_GENO_LABEL["GENO:0000137"] = "unspecified zygosity"
SO_LOSS = ("SO:0001743", "copy_number_loss")
SO_GAIN = ("SO:0001742", "copy_number_gain")
SO_CNV = ("SO:0001019", "copy_number_variation")
SO_REPEAT = ("SO:0002162", "short_tandem_repeat_expansion")
ACMG = {
    "pathogenic": "PATHOGENIC", "likely pathogenic": "LIKELY_PATHOGENIC",
    "uncertain significance": "UNCERTAIN_SIGNIFICANCE", "vus": "UNCERTAIN_SIGNIFICANCE",
    "variant of uncertain significance": "UNCERTAIN_SIGNIFICANCE",
    "likely benign": "LIKELY_BENIGN", "benign": "BENIGN",
    "p": "PATHOGENIC", "lp": "LIKELY_PATHOGENIC", "lb": "LIKELY_BENIGN", "b": "BENIGN",
    "p/lp": "LIKELY_PATHOGENIC", "lb/b": "LIKELY_BENIGN",
    "致病": "PATHOGENIC", "致病性": "PATHOGENIC", "致病变异": "PATHOGENIC",
    "可能致病": "LIKELY_PATHOGENIC", "可能致病性": "LIKELY_PATHOGENIC", "疑似致病": "LIKELY_PATHOGENIC",
    "意义不明": "UNCERTAIN_SIGNIFICANCE", "意义未明": "UNCERTAIN_SIGNIFICANCE",
    "临床意义未明": "UNCERTAIN_SIGNIFICANCE", "临床意义不明": "UNCERTAIN_SIGNIFICANCE",
    "可能良性": "LIKELY_BENIGN", "良性": "BENIGN",
}
# MONDO's root class, used as the diagnosis of an interpretation that carries variants but no
# diagnosis yet: Phenopacket v2 requires Diagnosis.disease, and "disease" is all that is known.
# read() does not report it as a diagnosis. (MONDO:0000001 "disease", checked in OLS4 2026-10-06.)
UNKNOWN_DISEASE = {"id": "MONDO:0000001", "label": "disease"}
DISEASE_PREFIXES = ("OMIM", "ORPHA", "MONDO")
_ID_VALUE = {"OMIM": re.compile(r"^\d{6}$"), "ORPHA": re.compile(r"^\d{1,7}$"), "MONDO": re.compile(r"^\d{7}$")}
MAX_DEPTH = 64
RESOURCES = {
    "HP": {"id": "hp", "name": "human phenotype ontology", "url": "http://purl.obolibrary.org/obo/hp.owl",
           "namespacePrefix": "HP", "iriPrefix": "http://purl.obolibrary.org/obo/HP_"},
    "GENO": {"id": "geno", "name": "Genotype Ontology", "url": "http://purl.obolibrary.org/obo/geno.owl",
             "namespacePrefix": "GENO", "iriPrefix": "http://purl.obolibrary.org/obo/GENO_"},
    "SO": {"id": "so", "name": "Sequence types and features ontology", "url": "http://purl.obolibrary.org/obo/so.owl",
           "namespacePrefix": "SO", "iriPrefix": "http://purl.obolibrary.org/obo/SO_"},
    "OMIM": {"id": "omim", "name": "An Online Catalog of Human Genes and Genetic Disorders",
             "url": "https://www.omim.org", "namespacePrefix": "OMIM", "iriPrefix": "https://www.omim.org/entry/"},
    "ORPHA": {"id": "orpha", "name": "Orphanet", "url": "https://www.orpha.net",
              "namespacePrefix": "ORPHA", "iriPrefix": "https://www.orpha.net/en/disease/detail/"},
    "MONDO": {"id": "mondo", "name": "Mondo Disease Ontology", "url": "http://purl.obolibrary.org/obo/mondo.owl",
              "namespacePrefix": "MONDO", "iriPrefix": "http://purl.obolibrary.org/obo/MONDO_"},
    "hgnc.symbol": {"id": "hgnc.symbol", "name": "HGNC gene symbol", "url": "https://www.genenames.org",
                    "namespacePrefix": "hgnc.symbol", "iriPrefix": "https://identifiers.org/hgnc.symbol:"},
}


class PhenopacketError(ValueError):
    """The input is not a phenopacket zebra can read (the message says why)."""


# ---------------------------------------------------------------- read


def _load(source: Union[str, Path, Dict[str, Any]]) -> Dict[str, Any]:
    if isinstance(source, dict):
        return source
    if not isinstance(source, (str, Path)):
        raise PhenopacketError(f"a phenopacket is a file path or a parsed JSON object, not {type(source).__name__}")
    path = Path(source)
    try:
        size = path.stat().st_size
    except OSError as err:
        raise PhenopacketError(f"cannot read {path}: {err.strerror or err}") from None
    if size > MAX_BYTES:
        raise PhenopacketError(f"{path} is {size:,} bytes; a phenopacket is at most {MAX_BYTES:,} (is this a cohort dump?)")
    try:
        text = path.read_text("utf-8-sig")
    except UnicodeDecodeError as err:
        raise PhenopacketError(f"{path} is not UTF-8 text ({err.reason})") from None
    except OSError as err:  # a directory, a permission problem
        raise PhenopacketError(f"cannot read {path}: {err.strerror or err}") from None
    try:
        data = json.loads(text)
    except ValueError as err:
        raise PhenopacketError(f"{path} is not JSON ({err})") from None
    except RecursionError:
        raise PhenopacketError(f"{path} is nested too deeply to be a phenopacket") from None
    if not isinstance(data, dict):
        raise PhenopacketError(f"{path}: a phenopacket is a JSON object, not {type(data).__name__}")
    return data


def _obj(value: Any) -> Dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _list(value: Any) -> List[Any]:
    return value if isinstance(value, list) else []


def _term(cls: Any) -> Tuple[Optional[str], Optional[str]]:
    c = _obj(cls)
    tid, label = c.get("id"), c.get("label")
    return (tid.strip() if isinstance(tid, str) and tid.strip() else None,
            label.strip() if isinstance(label, str) and label.strip() else None)


def render_time(element: Any, _depth: int = 0) -> Optional[str]:
    """A Phenopacket TimeElement (v2) or v1 onset/age as compact text; None when absent or unreadable."""
    if _depth > 4:
        return None
    if isinstance(element, str):
        return element.strip() or None
    t = _obj(element)
    if not t:
        return None
    if isinstance(t.get("iso8601duration"), str):  # bare Age
        return t["iso8601duration"]
    age = _obj(t.get("age"))
    if isinstance(age.get("iso8601duration"), str):
        return age["iso8601duration"]
    if isinstance(t.get("age"), str):  # v1 Age as a plain string
        return t["age"]
    rng = _obj(t.get("ageRange"))
    if rng:
        start = render_time(rng.get("start"), _depth + 1) or "?"
        end = render_time(rng.get("end"), _depth + 1) or "?"
        return f"{start}-{end}"
    ga = _obj(t.get("gestationalAge"))
    if isinstance(ga.get("weeks"), int) and not isinstance(ga.get("weeks"), bool):
        days = ga.get("days") if isinstance(ga.get("days"), int) and not isinstance(ga.get("days"), bool) else 0
        return f"gestational {ga['weeks']}w{days}d"
    tid, label = _term(t.get("ontologyClass") or (t if "id" in t else None))
    if tid:
        return f"{tid} {label}" if label else tid
    if isinstance(t.get("timestamp"), str):
        return t["timestamp"]
    iv = _obj(t.get("interval"))
    if iv:
        return f"{iv.get('start', '?')}/{iv.get('end', '?')}"
    return None


def _hpo_id(raw: Optional[str]) -> Optional[str]:
    m = HPO_RE.match(raw or "")
    return f"HP:{m.group(1)}" if m else None


def _expressions(desc: Dict[str, Any]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for e in _list(desc.get("expressions")):
        e = _obj(e)
        syntax, value = e.get("syntax"), e.get("value")
        if isinstance(syntax, str) and isinstance(value, str) and syntax.lower() in ("hgvs.c", "hgvs.g", "hgvs.p", "hgvs"):
            key = {"hgvs.c": "hgvs_c", "hgvs.g": "hgvs_g", "hgvs.p": "hgvs_p"}.get(syntax.lower())
            if key is None:  # bare "hgvs": decide by the sequence type
                key = "hgvs_g" if ":g." in value else "hgvs_p" if ":p." in value else "hgvs_c"
            out.setdefault(key, value)
    return out


def _variant_from_descriptor(desc: Dict[str, Any], status: Optional[str], acmg: Optional[str]) -> Dict[str, Any]:
    gene = _obj(desc.get("geneContext"))
    vcf = _obj(desc.get("vcfRecord"))
    zyg_id, zyg_label = _term(desc.get("allelicState"))
    stype_id, stype_label = _term(desc.get("structuralType"))
    row: Dict[str, Any] = {
        "kind": "cnv" if stype_id else "small",
        "gene": gene.get("symbol"), "gene_id": gene.get("valueId"),
        **_expressions(desc),
        "vcf": (f"{vcf.get('chrom')}-{vcf.get('pos')}-{vcf.get('ref')}-{vcf.get('alt')}"
                if all(vcf.get(k) not in (None, "") for k in ("chrom", "pos", "ref", "alt")) else None),
        "assembly": vcf.get("genomeAssembly"),
        "zygosity": zyg_label or _GENO_LABEL.get(zyg_id or ""),
        "acmg": acmg, "status": status,
        "description": desc.get("description") or desc.get("label"),
        "structural_type": f"{stype_id} {stype_label}".strip() if stype_id else None,
    }
    if stype_id == SO_REPEAT[0]:
        row["kind"] = "repeat_expansion"
    return {k: v for k, v in row.items() if v not in (None, "")}


def read(source: Union[str, Path, Dict[str, Any]]) -> Dict[str, Any]:
    """Read a Phenopacket (v2, or v1) from a file path or a parsed object. See the module docstring."""
    data = _load(source)
    try:
        return _read(data)
    except RecursionError:
        raise PhenopacketError("the phenopacket is nested too deeply to read") from None


def _flag(value: Any) -> bool:
    """A protobuf-JSON boolean: true, or the string "true"; anything else (false, "false", NaN) is false."""
    return value is True or (isinstance(value, str) and value.strip().lower() == "true")


def _read(data: Dict[str, Any]) -> Dict[str, Any]:
    warnings: List[str] = []
    if "phenotypicFeatures" not in data and isinstance(data.get("members"), list) and "proband" not in data:
        raise PhenopacketError("this is a Cohort message: read each member's phenopacket on its own")
    if "phenotypicFeatures" not in data and isinstance(data.get("proband"), dict):
        warnings.append("this is a Family message: only the proband's phenopacket was read")
        data = data["proband"]
    if not any(k in data for k in ("phenotypicFeatures", "interpretations", "diseases", "subject", "metaData")):
        raise PhenopacketError("not a phenopacket: none of phenotypicFeatures, interpretations, diseases, subject, "
                               "metaData is present")
    present: List[Dict[str, Any]] = []
    excluded: List[Dict[str, Any]] = []
    seen: Dict[str, str] = {}
    for i, feat in enumerate(_list(data.get("phenotypicFeatures"))):
        feat = _obj(feat)
        raw, label = _term(feat.get("type"))
        hid = _hpo_id(raw)
        if not hid:
            warnings.append(f"phenotypicFeatures[{i}]: {raw or 'no type id'} is not an HPO term; skipped")
            continue
        is_excluded = _flag(feat.get("excluded")) or _flag(feat.get("negated"))
        status = "excluded" if is_excluded else "present"
        if hid in seen:
            if seen[hid] != status:
                warnings.append(f"{hid} is recorded both present and excluded; kept as {seen[hid]} (first entry)")
            continue
        seen[hid] = status
        row = {"id": hid, "label": label, "onset": render_time(feat.get("onset") or feat.get("classOfOnset")
                                                               or feat.get("ageOfOnset"))}
        if isinstance(feat.get("description"), str) and feat["description"].strip():
            row["description"] = feat["description"].strip()
        (excluded if is_excluded else present).append(row)

    diseases: List[Dict[str, Any]] = []
    genes: List[Dict[str, Any]] = []
    variants: List[Dict[str, Any]] = []

    def add_disease(cls: Any, *, excluded_flag: bool, origin: str, status: Optional[str], onset: Any = None) -> None:
        did, dlabel = _term(cls)
        if not did or did == UNKNOWN_DISEASE["id"]:
            return  # MONDO:0000001 is the placeholder for "no diagnosis yet", not a diagnosis
        for d in diseases:
            if d["id"] == did:
                if status and not d.get("status"):
                    d["status"] = status
                return
        diseases.append({"id": did, "label": dlabel, "excluded": excluded_flag, "from": origin, "status": status,
                         "onset": render_time(onset)})

    def add_gene(symbol: Any, gid: Any, status: Optional[str]) -> None:
        if not (isinstance(symbol, str) and symbol.strip()) and not (isinstance(gid, str) and gid.strip()):
            return
        sym = symbol.strip() if isinstance(symbol, str) else None
        if any(g["symbol"] == sym and (g["id"] == gid or not gid) for g in genes):
            return
        genes.append({"symbol": sym, "id": gid if isinstance(gid, str) else None, "status": status})

    for k, interp in enumerate(_list(data.get("interpretations"))):
        interp = _obj(interp)
        progress = interp.get("progressStatus") if isinstance(interp.get("progressStatus"), str) else None
        diag = _obj(interp.get("diagnosis"))
        add_disease(diag.get("disease"), excluded_flag=False, origin="interpretation", status=progress)
        for gi in _list(diag.get("genomicInterpretations")):
            gi = _obj(gi)
            status = gi.get("interpretationStatus") if isinstance(gi.get("interpretationStatus"), str) else None
            if gi.get("gene"):
                g = _obj(gi["gene"])
                add_gene(g.get("symbol"), g.get("valueId"), status)
            vi = _obj(gi.get("variantInterpretation"))
            desc = _obj(vi.get("variationDescriptor"))
            if desc:
                acmg = vi.get("acmgPathogenicityClassification")
                row = _variant_from_descriptor(desc, status, acmg if isinstance(acmg, str) else None)
                variants.append(row)
                gc = _obj(desc.get("geneContext"))
                add_gene(gc.get("symbol"), gc.get("valueId"), status)
            elif not gi.get("gene"):
                warnings.append(f"interpretations[{k}]: a genomic interpretation with neither a gene nor a variant")
    for d in _list(data.get("diseases")):
        d = _obj(d)
        add_disease(d.get("term"), excluded_flag=_flag(d.get("excluded")), origin="diseases", status=None,
                    onset=d.get("onset") or d.get("classOfOnset") or d.get("ageOfOnset"))
    # v1: top-level genes and variants
    for g in _list(data.get("genes")):
        g = _obj(g)
        add_gene(g.get("symbol"), g.get("id"), None)
    for v in _list(data.get("variants")):
        v = _obj(v)
        hgvs = _obj(v.get("hgvsAllele")).get("hgvs")
        vcf = _obj(v.get("vcfAllele"))
        _zid, zlabel = _term(v.get("zygosity"))
        row = {"kind": "small", "status": None}
        if isinstance(hgvs, str):
            row["hgvs_g" if ":g." in hgvs else "hgvs_c"] = hgvs
        if vcf.get("chr") and vcf.get("pos") is not None:
            row["vcf"] = f"{vcf.get('chr')}-{vcf.get('pos')}-{vcf.get('ref')}-{vcf.get('alt')}"
            row["assembly"] = vcf.get("genomeAssembly")
        if zlabel:
            row["zygosity"] = zlabel
        variants.append({k: x for k, x in row.items() if x is not None})

    subject = _obj(data.get("subject"))
    meta = _obj(data.get("metaData"))
    sid = subject.get("id") if isinstance(subject.get("id"), str) else None
    sex = subject.get("sex") if isinstance(subject.get("sex"), str) else None
    refs = [r.get("id") for r in _list(meta.get("externalReferences")) if isinstance(_obj(r).get("id"), str)]
    return {
        "id": data.get("id") if isinstance(data.get("id"), str) else None,
        "schema_version": meta.get("phenopacketSchemaVersion"),
        "subject": {"id": sid, "sex": sex,
                    "age": render_time(subject.get("timeAtLastEncounter") or subject.get("ageAtCollection"))},
        "present": present, "excluded": excluded, "diseases": diseases, "genes": genes, "variants": variants,
        "references": refs, "warnings": warnings,
    }


# ---------------------------------------------------------------- write


def parse_age(text: Any) -> Optional[str]:
    """An ISO 8601 duration for age text that states one unambiguously, else None.

    "P1Y6M" stays; "6 months", "3岁", "18 mo" convert; a fraction converts only when it
    is a whole number of the next unit down (1.5 years = P18M). Anything else is None:
    an age is never estimated.
    """
    if not isinstance(text, str):
        return None  # a bare number does not say its unit
    s = text.strip()
    if ISO_DURATION.match(s.upper()):
        return s.upper()
    m = _AGE_TEXT.match(s)
    if not m:
        return None
    value = float(m.group(1))
    unit = _UNITS.get(m.group(2).lower()) or _UNITS.get(m.group(2))
    if unit is None:
        return None
    if value.is_integer():
        return f"P{int(value)}{unit}"
    down = {"Y": ("M", 12), "W": ("D", 7)}.get(unit)
    if down and (value * down[1]).is_integer():
        return f"P{int(value * down[1])}{down[0]}"
    return None


def onset_element(text: Any, hpo_labels: Optional[Dict[str, str]] = None) -> Optional[Dict[str, Any]]:
    """A TimeElement for onset text, or None when it cannot be stated without guessing."""
    if not isinstance(text, str) or not text.strip():
        return None
    s = text.strip()
    hid = _hpo_id(s.split()[0]) if s.split() else None
    if hid:
        label = (hpo_labels or {}).get(hid) or next((lab for i, lab in ONSET_CLASSES.values() if i == hid), None)
        return {"ontologyClass": {"id": hid, "label": label}} if label else None
    key = s.lower()
    cls = ONSET_CLASSES.get(key) or ONSET_CLASSES.get(key.rstrip("期"))
    if cls:
        return {"ontologyClass": {"id": cls[0], "label": cls[1]}}
    iso = parse_age(s)
    if iso:
        return {"age": {"iso8601duration": iso}}
    return None


def _sex(value: Any) -> str:
    s = str(value or "").strip().lower()
    if s in ("male", "m", "男", "男性", "boy", "男孩"):
        return "MALE"
    if s in ("female", "f", "女", "女性", "girl", "女孩"):
        return "FEMALE"
    if s in ("other", "other_sex", "其他"):
        return "OTHER_SEX"
    return "UNKNOWN_SEX"


def _acmg(variant: Dict[str, Any]) -> str:
    for raw in (variant.get("classification_lab"), _obj(variant.get("zebra_acmg")).get("classification")):
        if isinstance(raw, str) and raw.strip():
            key = re.sub(r"\s+", " ", raw.strip().lower().replace("_", " "))
            if key in ACMG:
                return ACMG[key]
    return "NOT_PROVIDED"


def _disease_class(hyp: Dict[str, Any], names: Dict[str, str], warn: List[str]) -> Optional[Dict[str, str]]:
    ids = hyp.get("ids")
    pairs: List[Tuple[str, Any]] = []
    if isinstance(ids, dict):
        pairs = [(str(k).upper(), v) for k, v in ids.items()]
    elif isinstance(ids, list):
        pairs = [tuple(x.split(":", 1)) for x in ids if isinstance(x, str) and ":" in x]  # type: ignore[misc]
    by_prefix: Dict[str, str] = {}
    for prefix, value in pairs:
        prefix = str(prefix).strip().upper()
        if not isinstance(value, (str, int)) or isinstance(value, bool):
            continue
        value = str(value).split(":", 1)[-1].strip()
        if prefix in _ID_VALUE and _ID_VALUE[prefix].match(value):
            by_prefix.setdefault(prefix, value)
        elif prefix in _ID_VALUE:
            warn.append(f"hypothesis id {prefix}:{value} is not a well-formed {prefix} id; left out")
    for prefix in DISEASE_PREFIXES:
        if prefix in by_prefix:
            did = f"{prefix}:{by_prefix[prefix]}"
            # the curated name when the HPO release knows the disease; the case's own wording otherwise
            # (checked with the rest of the document against the case's protected identifiers)
            label = names.get(did) or str(hyp.get("disease") or did).strip() or did
            return {"id": did, "label": label}
    return None


_SAFE_VARIANT_TEXT = ("region", "iscn", "cnv_type", "exons", "copy_number", "motif", "repeat_count", "hgvs_c")


def _variant_descriptor(v: Dict[str, Any], used: set, warn: List[str]) -> Dict[str, Any]:
    from zebra.case import KIND_ALIASES

    raw_kind = str(v.get("kind") or "small").strip().lower()
    kind = KIND_ALIASES.get(raw_kind, raw_kind)
    if kind not in ("small", "cnv", "exon_cnv", "copy_number", "repeat_expansion"):
        warn.append(f"variant {v.get('id')!r}: kind {v.get('kind')!r} is not one zebra knows; exported as a small variant")
        kind = "small"
    desc: Dict[str, Any] = {"id": str(v.get("id") or "variant"), "moleculeContext": "genomic"}
    gene = v.get("gene")
    if isinstance(gene, str) and re.match(r"^[A-Za-z0-9][A-Za-z0-9._@-]{0,30}$", gene.strip()):
        desc["geneContext"] = {"valueId": f"hgnc.symbol:{gene.strip()}", "symbol": gene.strip()}
        used.add("hgnc.symbol")
    elif gene not in (None, ""):
        warn.append(f"variant {v.get('id')!r}: gene {gene!r} is not a gene symbol; left out")
    if v.get("description") or v.get("method") or v.get("note"):
        warn.append(f"variant {v.get('id')!r}: free-text description/method/note not exported")
    if kind == "small":
        exprs = [{"syntax": syntax, "value": v[key]} for key, syntax in
                 (("hgvs_c", "hgvs.c"), ("hgvs_g", "hgvs.g"), ("hgvs_p", "hgvs.p")) if isinstance(v.get(key), str)]
        if exprs:
            desc["expressions"] = exprs
        vcf = v.get("vcf")
        if isinstance(vcf, str):
            parts = re.split(r"[-:\s]+", vcf.strip())
            assembly = v.get("assembly")
            if not (len(parts) == 4 and parts[1].isdigit()):
                warn.append(f"variant {v.get('id')!r}: vcf {vcf!r} is not chrom-pos-ref-alt; left out")
            elif not (isinstance(assembly, str) and assembly.strip()):
                warn.append(f"variant {v.get('id')!r}: no genome assembly recorded, so its VCF position is left out "
                            "(the HGVS is kept)")
            else:
                desc["vcfRecord"] = {"genomeAssembly": assembly.strip(), "chrom": parts[0], "pos": int(parts[1]),
                                     "ref": parts[2], "alt": parts[3]}
    else:
        bits = [f"{k}: {v[k]}" for k in _SAFE_VARIANT_TEXT if v.get(k) not in (None, "", [])]
        desc["label"] = f"{kind}{' ' + gene.strip() if isinstance(gene, str) and 'geneContext' in desc else ''}"
        desc["description"] = "; ".join(bits) or kind
        if kind == "repeat_expansion":
            stype = SO_REPEAT
        else:
            ctype = str(v.get("cnv_type") or "").lower()
            if ctype in ("loss", "deletion"):
                stype = SO_LOSS
            elif ctype in ("gain", "duplication", "amplification"):
                stype = SO_GAIN
            else:  # a copy number alone does not say loss or gain without the normal count (X in males)
                stype = SO_CNV
        desc["structuralType"] = {"id": stype[0], "label": stype[1]}
        used.add("SO")
    zraw = str(v.get("zygosity") or "").strip().lower()
    zyg = ZYGOSITY.get(zraw)
    if zyg:
        desc["allelicState"] = {"id": zyg[0], "label": zyg[1]}
        used.add("GENO")
    elif zraw in ("unknown", "unspecified", "未知", "不详"):
        desc["allelicState"] = {"id": "GENO:0000137", "label": "unspecified zygosity"}
        used.add("GENO")
    elif zraw:
        warn.append(f"variant {v.get('id')!r}: zygosity {v.get('zygosity')!r} has no GENO class here (mosaic is not "
                    "a zygosity); left out")
    return desc


def _strings(obj: Any, out: List[str]) -> List[str]:
    if isinstance(obj, str):
        out.append(obj)
    elif isinstance(obj, (int, float)) and not isinstance(obj, bool):
        out.append(str(obj))  # a record number stored as a number (copy_number, repeat_count) is checked too
    elif isinstance(obj, dict):
        for k, v in obj.items():
            _strings(v, out)
    elif isinstance(obj, list):
        for v in obj:
            _strings(v, out)
    return out


def _privacy_check(packet: Dict[str, Any], identifiers: List[str]) -> None:
    """Refuse the export when any protected identifier (or a Chinese resident-ID shape) is in it."""
    text = "\n".join(_strings(packet, []))
    try:
        from zebra.report_export import identifier_hits
    except ImportError:  # the report layer is optional; fall back to a plain containment check
        identifier_hits = None
    if identifier_hits is not None:
        why, _warnings = identifier_hits(text, identifiers)
    else:
        why = next((f"protected identifier #{n} of the case" for n, ident in enumerate(identifiers, 1)
                    if isinstance(ident, str) and len(ident.strip()) >= 2 and ident.strip() in text), None)
    if why:
        raise PhenopacketError(f"not exported: the phenopacket would contain {why}. Remove it from the case's "
                               "labels or hypothesis names and export again")


def from_case(case: Dict[str, Any], *, created: Optional[str] = None, hpo_version: Optional[str] = None,
              warnings: Optional[List[str]] = None) -> Dict[str, Any]:
    """A zebra case (zebra.case/1, as loaded by zebra.case.load) as a Phenopacket v2 JSON object.

    `created` (ISO 8601) defaults to now; `hpo_version` to the local release when
    one is installed. Whatever could not be represented is appended to `warnings`.
    Raises PhenopacketError when the result would carry a protected identifier.
    """
    import uuid

    if not isinstance(case, dict):
        raise PhenopacketError(f"a zebra case is a JSON object, not {type(case).__name__}")
    warn = warnings if warnings is not None else []
    used = {"HP"}
    hpo_names: Dict[str, str] = {}
    disease_names: Dict[str, str] = {}
    try:
        from zebra import hpo_local

        idx = hpo_local.load()
        hpo_names, disease_names = idx.names, idx.disease_names
        if hpo_version is None:
            hpo_version = idx.version
    except Exception:  # noqa: BLE001 - optional: labels then come from the case, the version is "unknown"
        if hpo_version is None:
            hpo_version = "unknown"
    subject_id = "proband"
    proband = _obj(case.get("proband"))
    subject: Dict[str, Any] = {"id": subject_id, "sex": _sex(proband.get("sex"))}
    age = parse_age(proband.get("age"))
    if age:
        subject["timeAtLastEncounter"] = {"age": {"iso8601duration": age}}
    elif proband.get("age") not in (None, ""):
        warn.append("proband age is not an unambiguous age (a number with a unit); left out")

    features = []
    seen = set()
    for p in _list(case.get("phenotypes")):
        p = _obj(p)
        hid = _hpo_id(str(p.get("id") or ""))
        if not hid:
            warn.append(f"phenotype {p.get('id')!r} is not an HPO id; left out")
            continue
        status = p.get("status", "present")
        if status in (None, ""):
            status = "present"
        if status not in ("present", "excluded"):
            warn.append(f"phenotype {hid}: status {status!r} is neither present nor excluded; left out")
            continue
        if hid in seen:
            warn.append(f"phenotype {hid} is recorded twice; the first entry was exported")
            continue
        seen.add(hid)
        label = hpo_names.get(hid) or (p.get("label") if isinstance(p.get("label"), str) else None) or hid
        feat: Dict[str, Any] = {"type": {"id": hid, "label": label}}
        if status == "excluded":
            feat["excluded"] = True
        onset = p.get("onset")
        el = onset_element(onset, hpo_names)
        if el:
            feat["onset"] = el
        elif isinstance(onset, str) and onset.strip():
            warn.append(f"phenotype {hid}: onset is free text, not an age or an HPO onset class; not exported")
        features.append(feat)

    interpretations: List[Dict[str, Any]] = []
    diseases: List[Dict[str, Any]] = []
    attach_to: Optional[Dict[str, Any]] = None
    hyps = [_obj(h) for h in _list(case.get("hypotheses"))]
    for status in ("confirmed", "leading"):
        for h in hyps:
            if h.get("status") != status:
                continue
            cls = _disease_class(h, disease_names, warn)
            if cls is None:
                warn.append(f"a {status} hypothesis has no OMIM/ORPHA/MONDO id; not exported")
                continue
            used.add(cls["id"].split(":")[0])
            interp = {"id": f"interpretation-{len(interpretations) + 1}",
                      "progressStatus": "SOLVED" if status == "confirmed" else "IN_PROGRESS",
                      "diagnosis": {"disease": cls}}
            interpretations.append(interp)
            if status == "confirmed":
                diseases.append({"term": dict(cls)})
            if attach_to is None:
                attach_to = interp
    for h in hyps:
        if h.get("status") == "excluded":
            cls = _disease_class(h, disease_names, warn)
            if cls:
                used.add(cls["id"].split(":")[0])
                diseases.append({"term": cls, "excluded": True})

    genomic = []
    for v in _list(case.get("variants")):
        v = _obj(v)
        genomic.append({"subjectOrBiosampleId": subject_id, "interpretationStatus": "CANDIDATE",
                        "variantInterpretation": {"acmgPathogenicityClassification": _acmg(v),
                                                  "therapeuticActionability": "UNKNOWN_ACTIONABILITY",
                                                  "variationDescriptor": _variant_descriptor(v, used, warn)}})
    if genomic:
        if attach_to is None:
            used.add("MONDO")
            attach_to = {"id": "variants", "progressStatus": "IN_PROGRESS",
                         "summary": "variants recorded in the case; no diagnosis with an ontology id yet",
                         "diagnosis": {"disease": dict(UNKNOWN_DISEASE)}}
            interpretations.append(attach_to)
            warn.append("variants exported in an interpretation without a diagnosis (disease MONDO:0000001, "
                        "status IN_PROGRESS): the case has no confirmed or leading hypothesis with an id")
        attach_to["diagnosis"]["genomicInterpretations"] = genomic

    if created is None:
        from datetime import datetime, timezone

        created = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    resources = []
    for key in ("HP", "GENO", "SO", "OMIM", "ORPHA", "MONDO", "hgnc.symbol"):
        if key in used:
            res = dict(RESOURCES[key])
            res["version"] = hpo_version if key == "HP" else "unknown"
            resources.append(res)
    packet: Dict[str, Any] = {
        "id": "zebra-" + uuid.uuid4().hex[:16],  # random: a hash of the case id could be reversed to a name
        "subject": subject,
        "phenotypicFeatures": features,
    }
    if interpretations:
        packet["interpretations"] = interpretations
    if diseases:
        packet["diseases"] = diseases
    packet["metaData"] = {"created": created, "createdBy": f"zebra-mod {__version__}", "resources": resources,
                          "phenopacketSchemaVersion": SCHEMA_VERSION,
                          "externalReferences": []}
    identifiers = [i for i in _list(_obj(case.get("privacy")).get("identifiers")) if isinstance(i, str)]
    _privacy_check(packet, identifiers)
    return packet
