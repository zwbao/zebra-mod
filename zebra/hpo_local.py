"""Offline phenotype-driven ranking from the HPO release files.

Data (fetched once by `zebra hpo fetch`, ~80 MB, into ~/.cache/zebra-mod/hpo):
  hp.json                 the ontology (obographs JSON)
  phenotype.hpoa          disease -> phenotype annotations (OMIM, ORPHA, DECIPHER)
  genes_to_phenotype.txt  gene -> phenotype, with the disease each pair comes from

Two scoring methods share one index (`rank(..., method=)`):

"resnik" (zebra 0.1.0): information content IC(t) = -ln(share of diseases
annotated with t or a descendant); a query term's match to a disease is the IC of
their most informative common ancestor (Resnik); a disease's score is the mean over
the query terms of their best match (query -> disease best-match average, as the
Phenomizer does). Each excluded term the disease is annotated with (directly or
below) costs its IC times the annotation frequency, SUMMED over distinct excluded
terms.

"lr" (added in 0.2, after the phenopacket-store benchmark, docs/BENCHMARK.md): a
likelihood-ratio score in the spirit of LIRICAL (Robinson et al., AJHG 2020) -- not
LIRICAL itself: no genotype term, no onset or sex model, no pretest probabilities.
For each present term q, log LR(q | D) = max over D's annotations t of
[ IC(MICA(q, t)) + ln f(t) ]: P(q | D) is estimated as the annotation frequency f
times the share of diseases with the common ancestor that also have q, and P(q | not
D) as the share of all diseases that have q, so the ratio reduces to that
expression. Frequencies come from HPOA (case counts pooled across papers, the
Orphanet frequency classes, percentages); an annotation with no stated frequency
counts as `params["f_unknown"]`. Each excluded term e adds ln(1 - P(e | D)) with
P(e | D) the frequency of D's annotations at or below e, capped so one curated
feature cannot veto a disease outright. A disease's score is the sum (natural log).

Curated negative evidence (both methods): a disease whose annotation carries HPOA's
`NOT` qualifier, or the "Excluded" (HP:0040285) frequency class, for a term the
patient is said to HAVE (or an ancestor of it) is penalised -- unless the same
disease also annotates that term positively. A case count such as "0/2" is NOT
negative curation: it says two reported patients lacked the feature, so it is pooled
with the disease's other counts for that term (6/11 + 1/1 + 0/2 = 7/14) or, when it
is the only row, simply ignored (D-P1-3).

Excluded terms that contradict the query (the same term given as present, or an
ancestor of a present term: febrile seizure present, seizure excluded) are dropped
with a note rather than scored. Scores are reported raw and as a share of the
query's own maximum; `ties` says how many diseases share the top score.

These are transparent baselines that run offline, reproducibly, from public files;
they are not Exomiser or LIRICAL.
"""

from __future__ import annotations

import csv
import json
import math
import os
import pickle
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

from zebra.http import cache_dir

FILES = ("hp.json", "phenotype.hpoa", "genes_to_phenotype.txt")
RELEASE_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/{name}"
# official Chinese labels (HPO translation project); optional, enables Chinese search
ZH_FILE = "hp-zh.babelon.tsv"
ZH_URL = "https://raw.githubusercontent.com/obophenotype/hpo-translations/main/babelon/hp-zh.babelon.tsv"
_CJK = re.compile(r"[㐀-鿿]")
ROOT_PHENO = "HP:0000118"  # Phenotypic abnormality
EXCLUDED_CLASS = "HP:0040285"  # frequency class "Excluded" (0%): curated absence

METHODS = ("resnik", "lr")
DEFAULT_METHOD = "resnik"
# LR parameters. Chosen on the development half of the phenopacket-store benchmark only
# (tools/bench/tune.py); the held-out half is reported, never tuned on.
LR_PARAMS: Dict[str, float] = {
    "f_unknown": 1.0,     # frequency of an annotation that states none
    "f_floor": 0.02,      # smallest frequency used inside ln f (a 0/1 pooled to 0 cannot be -inf)
    "excluded_cap": 0.9,  # P(e | D) for an excluded term is capped here: ln(1 - 0.9) = -2.3 at most
    "excluded_min_f": 0.0,  # an excluded term costs nothing unless D has it at least this often
    "not_penalty": -2.3,  # ln LR for a present term the disease is curated NOT to have (ln 0.1)
    "unexplained": 0.0,   # ln LR for a query term the disease matches only at the root
}
# Resnik: the 0.1.0 method; `excluded_weight` scales its excluded-term penalty (1.0 = 0.1.0)
_PARAM_RANGE: Dict[str, Tuple[float, float]] = {
    "f_unknown": (0.001, 1.0), "f_floor": (0.0001, 1.0), "excluded_cap": (0.0, 0.999),
    "excluded_min_f": (0.0, 1.0), "not_penalty": (-50.0, 0.0), "unexplained": (-50.0, 0.0),
    "excluded_weight": (0.0, 10.0), "resnik_excluded_min_f": (0.0, 1.0), "tiebreak_frequency": (0.0, 1.0),
}
RESNIK_PARAMS: Dict[str, float] = {
    "excluded_weight": 1.0,          # scales the excluded-term penalty (1.0 = 0.1.0)
    "resnik_excluded_min_f": 0.0,    # only annotations at least this frequent are penalised
    "tiebreak_frequency": 0.0,       # 1: inside a tie, more frequent matched annotations first (CP1-16)
}

# Lay phrases families actually use, each mapped to the official HPO term(s) a clinician
# would record, best match first. This is a translation table, not a guess: every id below
# was checked against the HPO 2026-09-01 release's own labels, and
# tests/test_hpo_local.py::test_f36_lay_terms_point_at_real_ids_in_the_release re-checks
# that each still exists, while test_d_p1_4_live_* checks that the phrase ranks its
# intended term FIRST in the installed release (a wrong but real id passes the first
# check, not the second).
# Needed because the ontology's labels are clinical: "走路晚" (walked late) matches nothing
# literally, "智力低下" is closer by character bigrams to 面部肌张力低下 (facial hypotonia)
# than to its own official label 智力障碍, and "听力下降" (hearing loss) is closer to
# 视力下降 (reduced visual acuity) than to 听力受损.
# A phrase that appears inside a longer query ("孩子走路晚", "my son walked late") is
# matched too, unless a negation precedes it ("没有抽搐", "no seizures").
LAY_TERMS: Dict[str, List[str]] = {
    # --- Chinese: motor milestones
    "走路晚": ["HP:0031936", "HP:0001270"],      # Delayed ability to walk; Motor delay
    "学步晚": ["HP:0031936", "HP:0001270"],
    "还不会走路": ["HP:0031936", "HP:0001270"],
    "还不会走": ["HP:0031936", "HP:0001270"],
    "不会走路": ["HP:0002540", "HP:0031936"],    # Inability to walk; Delayed ability to walk (D-P2-7)
    "不会走": ["HP:0002540", "HP:0031936"],
    "翻身晚": ["HP:0032989", "HP:0001270"],      # Delayed ability to roll over
    "不会翻身": ["HP:0032989", "HP:0001270"],
    "坐得晚": ["HP:0025336", "HP:0001270"],      # Delayed ability to sit
    "不会坐": ["HP:0025336", "HP:0001270"],
    "坐不稳": ["HP:0025336", "HP:0001270"],
    "抬头晚": ["HP:0002421", "HP:0032988"],      # Poor head control; Persistent head lag
    "软趴趴": ["HP:0008947", "HP:0001252"],      # Floppy infant; Hypotonia
    "孩子软": ["HP:0001252", "HP:0008947"],      # Hypotonia; Floppy infant
    "身体软": ["HP:0001252", "HP:0008947"],
    "肌张力低": ["HP:0001252"],                   # Hypotonia (official zh label 肌张力减退)
    "肌张力低下": ["HP:0001252"],
    "肌张力偏低": ["HP:0001252"],
    "肌张力减低": ["HP:0001252"],
    "没力气": ["HP:0001324", "HP:0001252"],      # Muscle weakness; Hypotonia
    "肌无力": ["HP:0001324"],
    "走路不稳": ["HP:0002317", "HP:0001251"],    # Unsteady gait; Ataxia
    "走路摇晃": ["HP:0002317", "HP:0001251"],
    "走路摇摇晃晃": ["HP:0002317", "HP:0001251"],
    "退行": ["HP:0002376"],                       # Developmental regression
    "倒退": ["HP:0002376"],
    "发育倒退": ["HP:0002376"],
    "发育退行": ["HP:0002376"],
    # --- Chinese: seizures
    "抽风": ["HP:0001250"],                       # Seizure
    "抽搐": ["HP:0001250"],
    "惊厥": ["HP:0001250"],
    "发热惊厥": ["HP:0002373"],                   # Febrile seizure
    "高热惊厥": ["HP:0002373"],
    "发烧抽筋": ["HP:0002373"],
    "发烧抽搐": ["HP:0002373"],
    "发烧抽风": ["HP:0002373"],
    "发热抽搐": ["HP:0002373"],
    # --- Chinese: speech and cognition
    "不会说话": ["HP:0001344", "HP:0000750"],     # Absent speech; Delayed speech and language development
    "不说话": ["HP:0001344"],
    "说话晚": ["HP:0000750"],
    "说话不清楚": ["HP:0001350", "HP:0001260"],   # Slurred speech; Dysarthria
    "说话不清": ["HP:0001350", "HP:0001260"],
    "口齿不清": ["HP:0001350", "HP:0001260"],
    "智力低下": ["HP:0001249"],                   # Intellectual disability
    "智力落后": ["HP:0001249"],
    "智力障碍": ["HP:0001249"],
    "发育迟缓": ["HP:0001263"],                   # Global developmental delay
    "发育落后": ["HP:0001263"],
    # --- Chinese: growth, head, senses, feeding, labs
    "不长个": ["HP:0004322"],                     # Short stature
    "长不高": ["HP:0004322"],
    "个子矮": ["HP:0004322"],
    "头小": ["HP:0000252"],                       # Microcephaly
    "小头": ["HP:0000252"],
    "头围小": ["HP:0040195", "HP:0000252"],       # Decreased head circumference; Microcephaly
    "头围偏小": ["HP:0040195", "HP:0000252"],
    "头围减小": ["HP:0040195", "HP:0000252"],
    "头大": ["HP:0000256"],                       # Macrocephaly
    "头围大": ["HP:0040194", "HP:0000256"],       # Increased head circumference; Macrocephaly
    "头围偏大": ["HP:0040194", "HP:0000256"],
    "听不见": ["HP:0000365"],                     # Hearing impairment
    "耳朵听不见": ["HP:0000365"],
    "听力下降": ["HP:0000365"],
    "听力差": ["HP:0000365"],
    "听力不好": ["HP:0000365"],
    "听力减退": ["HP:0000365"],
    "听力损失": ["HP:0000365"],
    "耳聋": ["HP:0000365"],                       # (HP:0000365 carries the synonym "Deafness")
    "聋": ["HP:0000365"],
    "看不清": ["HP:0000505"],                     # Visual impairment
    "视力差": ["HP:0000505", "HP:0007663"],       # Visual impairment; Reduced visual acuity
    "喂奶困难": ["HP:0008872"],                   # Feeding difficulties in infancy
    "吃奶费劲": ["HP:0008872"],
    "吞咽困难": ["HP:0002015"],                   # Dysphagia
    "尿有怪味": ["HP:0012088"],                   # Abnormal urinary odor
    "尿味怪": ["HP:0012088"],
    "尿有异味": ["HP:0012088"],
    "尿液有异味": ["HP:0012088"],
    "CK高": ["HP:0003236"],                       # Elevated circulating creatine kinase activity
    "肌酸激酶高": ["HP:0003236"],
    "肌酸激酶升高": ["HP:0003236"],
    # --- English lay phrasing
    "floppy baby": ["HP:0008947", "HP:0001252"],
    "floppy infant": ["HP:0008947"],
    "late walker": ["HP:0031936", "HP:0001270"],
    "walking late": ["HP:0031936", "HP:0001270"],
    "walked late": ["HP:0031936", "HP:0001270"],
    "not walking": ["HP:0031936", "HP:0001270"],
    "cannot walk": ["HP:0002540", "HP:0031936"],
    "not talking": ["HP:0001344", "HP:0000750"],
    "no speech": ["HP:0001344"],
    "cannot talk": ["HP:0001344"],
    "speech delay": ["HP:0000750"],
    "fit": ["HP:0001250"],
    "fits": ["HP:0001250"],
    "fitting": ["HP:0001250"],
    "convulsion": ["HP:0001250"],
    "convulsions": ["HP:0001250"],
    "febrile convulsion": ["HP:0002373"],
    "febrile convulsions": ["HP:0002373"],
    "small head": ["HP:0000252"],
    "mental retardation": ["HP:0001249"],
    "slow learner": ["HP:0001328"],               # Specific learning disability (D-P2-7: was ID)
    "developmental delay": ["HP:0001263"],
    "weak muscles": ["HP:0001324", "HP:0001252"],
    "poor feeding": ["HP:0008872"],
    "hard of hearing": ["HP:0000365"],
    "hearing loss": ["HP:0000365"],
    "smelly urine": ["HP:0012088"],
    "high ck": ["HP:0003236"],
}
# a lay phrase inside a longer query is ignored when one of these precedes it
_NEGATION_ZH = ("没有", "没", "无", "否认", "未见", "未", "不是", "不")
_NEGATION_EN = ("no", "not", "without", "never", "denies", "denied", "absent", "negative")
_PURL = re.compile(r"^http://purl\.obolibrary\.org/obo/HP_(\d{7})$")
INDEX_VERSION = 6  # bumped: 0/n counts pooled, not curated NOT (D-P1-3); closure carries frequencies


class HpoDataMissing(Exception):
    pass


def data_dir() -> Path:
    return Path(os.environ.get("ZEBRA_HPO_DIR") or (cache_dir() / "hpo"))


def missing_files() -> List[str]:
    return [f for f in FILES if not (data_dir() / f).exists() or (data_dir() / f).stat().st_size == 0]


def _curie(iri: str) -> Optional[str]:
    m = _PURL.match(iri)
    return f"HP:{m.group(1)}" if m else (iri if iri.startswith("HP:") else None)


_CLASS_FREQ = {
    "HP:0040280": 1.0,    # Obligate (100%)
    "HP:0040281": 0.9,    # Very frequent (80-99%)
    "HP:0040282": 0.55,   # Frequent (30-79%)
    "HP:0040283": 0.17,   # Occasional (5-29%)
    "HP:0040284": 0.025,  # Very rare (1-4%)
    EXCLUDED_CLASS: 0.0,  # Excluded (0%): curated absence
}


def _freq_weight(raw: str) -> float:
    """phenotype.hpoa frequency -> a weight in [0, 1]; 1.0 when no frequency is stated."""
    raw = (raw or "").strip()
    if raw in _CLASS_FREQ:
        return _CLASS_FREQ[raw]
    m = re.match(r"^(\d+)/(\d+)$", raw)
    if m and int(m.group(2)) > 0:
        return min(1.0, int(m.group(1)) / int(m.group(2)))
    m = re.match(r"^([\d.]+)%$", raw)
    if m:
        return min(1.0, float(m.group(1)) / 100)
    return 1.0


def _encode(known: float, unknown: bool) -> float:
    """One float per (term, disease) in the closure: the highest stated frequency of the
    disease's annotations at or below the term, plus 2.0 when one of them states none."""
    return known + (2.0 if unknown else 0.0)


def freq_at(code: float, f_unknown: float) -> float:
    """Decode `_encode`: the frequency to use, with `f_unknown` for an annotation that states none."""
    if code >= 2.0:
        return max(code - 2.0, f_unknown)
    return code


@dataclass
class Index:
    names: Dict[str, str]
    synonyms: Dict[str, List[str]]
    parents: Dict[str, List[str]]
    obsolete: Dict[str, Optional[str]]  # id -> replaced_by
    alt: Dict[str, str]  # alt id -> primary
    disease_names: Dict[str, str]
    # disease -> term -> frequency weight (pooled case counts, class or percentage; 1.0 if none stated)
    disease_terms: Dict[str, Dict[str, float]]
    disease_genes: Dict[str, List[str]]
    ic: Dict[str, float] = field(default_factory=dict)
    # term -> {disease index: encoded frequency (see _encode)} over the annotation closure:
    # a disease is listed under every ancestor of every term it is annotated with
    term_diseases: Dict[str, Dict[int, float]] = field(default_factory=dict)
    disease_list: List[str] = field(default_factory=list)
    version: str = ""
    zh: Dict[str, str] = field(default_factory=dict)
    # Curated negative evidence, kept rather than discarded: disease -> term -> why.
    # Two channels in phenotype.hpoa say "this disease does NOT have this feature" --
    # the `NOT` qualifier, and the "Excluded" (HP:0040285, 0%) frequency class. A case
    # count "0/n" is not one of them (D-P1-3).
    disease_excluded: Dict[str, Dict[str, str]] = field(default_factory=dict)
    # disease -> terms annotated without any stated frequency
    freq_unknown: Dict[str, Set[str]] = field(default_factory=dict)
    _anc: Dict[str, Set[str]] = field(default_factory=dict, repr=False, compare=False)

    def ancestors(self, term: str) -> Set[str]:
        """The term and all its ancestors (memoised; treat the returned set as read-only)."""
        got = self._anc.get(term)
        if got is not None:
            return got
        seen: Set[str] = set()
        stack = [term]
        while stack:
            t = stack.pop()
            if t in seen:
                continue
            seen.add(t)
            stack.extend(self.parents.get(t, ()))
        self._anc[term] = seen
        return seen

    def primary(self, term: str) -> Tuple[Optional[str], Optional[str]]:
        """Resolve alt ids and obsolete ids; returns (id, note)."""
        if term in self.alt:
            return self.alt[term], f"{term} is an alternative id of {self.alt[term]}"
        if term in self.obsolete:
            repl = self.obsolete[term]
            if repl:
                return repl, f"{term} is obsolete; replaced by {repl}"
            return None, f"{term} is obsolete with no replacement"
        if term not in self.names:
            return None, f"{term} is not in this HPO release"
        return term, None


def _parse_ontology(path: Path):
    doc = json.loads(path.read_text("utf-8"))
    graph = doc["graphs"][0]
    names: Dict[str, str] = {}
    synonyms: Dict[str, List[str]] = {}
    obsolete: Dict[str, Optional[str]] = {}
    alt: Dict[str, str] = {}
    for node in graph.get("nodes", []):
        cid = _curie(node.get("id", ""))
        if not cid:
            continue
        meta = node.get("meta") or {}
        label = node.get("lbl") or ""
        props = {p.get("pred", "").rsplit("#", 1)[-1].rsplit("/", 1)[-1]: p.get("val") for p in meta.get("basicPropertyValues", []) or []}
        if meta.get("deprecated"):
            repl = props.get("IAO_0100001") or props.get("term_replaced_by")
            obsolete[cid] = _curie(repl) if repl else None
            continue
        names[cid] = label
        synonyms[cid] = [s.get("val", "") for s in meta.get("synonyms", []) or [] if s.get("val")]
        for p in meta.get("basicPropertyValues", []) or []:
            if p.get("pred", "").endswith("hasAlternativeId") and p.get("val"):
                alt[p["val"]] = cid
    parents: Dict[str, List[str]] = defaultdict(list)
    for edge in graph.get("edges", []):
        if edge.get("pred") not in ("is_a", "http://www.w3.org/2000/01/rdf-schema#subClassOf"):
            continue
        sub, obj = _curie(edge.get("sub", "")), _curie(edge.get("obj", ""))
        if sub and obj and sub in names and obj in names:
            parents[sub].append(obj)
    version = ""
    meta = graph.get("meta") or {}
    for p in meta.get("basicPropertyValues", []) or []:
        if p.get("pred", "").endswith("versionInfo"):
            version = p.get("val", "")
    if not version:
        version = str(meta.get("version", ""))
    return names, synonyms, dict(parents), obsolete, alt, version


def _parse_hpoa(path: Path):
    """disease names, positive annotations (with frequencies), curated NOT annotations, unknown-frequency terms.

    Frequencies of one disease-term pair from several rows are combined: case counts
    are pooled (6/11 + 1/1 + 0/2 = 7/14), otherwise the highest class or percentage is
    kept, and a stated frequency beats a row that states none. A pair whose only rows
    are zero counts ("0/2": two reported patients lacked it) is not annotated at all:
    it is neither a feature of the disease nor a curated absence (D-P1-3). Only the
    `NOT` qualifier and the "Excluded" class HP:0040285 are curated absence.
    """
    disease_names: Dict[str, str] = {}
    counts: Dict[Tuple[str, str], List[int]] = {}
    stated: Dict[Tuple[str, str], float] = {}
    unstated: Set[Tuple[str, str]] = set()
    disease_excluded: Dict[str, Dict[str, str]] = defaultdict(dict)
    with open(path, encoding="utf-8") as fh:
        rows = (line for line in fh if not line.startswith("#"))
        reader = csv.DictReader(rows, delimiter="\t")
        for row in reader:
            if row.get("aspect") != "P":
                continue
            did = row["database_id"]
            disease_names[did] = row.get("disease_name", did)
            term = row["hpo_id"]
            freq = (row.get("frequency") or "").strip()
            if (row.get("qualifier") or "").upper() == "NOT":
                disease_excluded[did][term] = "HPOA NOT qualifier"
                continue
            if freq == EXCLUDED_CLASS:
                disease_excluded[did][term] = f"HPOA frequency class {EXCLUDED_CLASS} (Excluded, 0%)"
                continue
            key = (did, term)
            m = re.match(r"^(\d+)/(\d+)$", freq)
            if m and int(m.group(2)) > 0:
                c = counts.setdefault(key, [0, 0])
                c[0] += min(int(m.group(1)), int(m.group(2)))
                c[1] += int(m.group(2))
            elif freq and not m:
                stated[key] = max(_freq_weight(freq), stated.get(key, 0.0))
            else:  # no frequency, or a count with a zero denominator
                unstated.add(key)
    disease_terms: Dict[str, Dict[str, float]] = defaultdict(dict)
    freq_unknown: Dict[str, Set[str]] = defaultdict(set)
    for key in set(counts) | set(stated) | unstated:
        did, term = key
        w = 0.0
        if key in counts and counts[key][0] > 0:
            w = counts[key][0] / counts[key][1]
        if stated.get(key, 0.0) > 0:  # a count and a class for the same pair: keep the higher estimate
            w = max(w, stated[key])
        if w > 0:
            disease_terms[did][term] = w
        elif key in unstated:  # curated with no frequency (zero counts elsewhere do not override it)
            disease_terms[did][term] = 1.0
            freq_unknown[did].add(term)
        # else: only zero counts ("0/2"): neither a feature nor a curated absence
    return (disease_names, dict(disease_terms), {d: dict(t) for d, t in disease_excluded.items()},
            {d: set(t) for d, t in freq_unknown.items()})


def _parse_g2p(path: Path) -> Dict[str, List[str]]:
    genes: Dict[str, Set[str]] = defaultdict(set)
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t")
        for row in reader:
            did = row.get("disease_id")
            sym = row.get("gene_symbol")
            if did and sym and sym.strip() not in ("-", ""):
                genes[did].add(sym.strip())
    return {d: sorted(g) for d, g in genes.items()}


def _parse_zh(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("predicate_id") == "rdfs:label" and row.get("translation_value") and row.get("subject_id", "").startswith("HP:"):
                out[row["subject_id"]] = row["translation_value"].strip()
    return out


def load(rebuild: bool = False) -> Index:
    missing = missing_files()
    if missing:
        raise HpoDataMissing(f"HPO files missing in {data_dir()}: {', '.join(missing)} (run: zebra hpo fetch)")
    d = data_dir()
    present = list(FILES) + ([ZH_FILE] if (d / ZH_FILE).exists() else [])
    stamp = "|".join(f"{f}:{(d / f).stat().st_size}:{int((d / f).stat().st_mtime)}" for f in present)
    cache = d / "index.pickle"
    if cache.exists() and not rebuild:
        try:
            with open(cache, "rb") as fh:
                saved = pickle.load(fh)
            if saved.get("stamp") == stamp and saved.get("v") == INDEX_VERSION:
                return saved["index"]
        except Exception:  # noqa: BLE001 - rebuild on any problem
            pass
    names, synonyms, parents, obsolete, alt, version = _parse_ontology(d / "hp.json")
    disease_names, disease_terms, disease_excluded, freq_unknown = _parse_hpoa(d / "phenotype.hpoa")
    disease_genes = _parse_g2p(d / "genes_to_phenotype.txt")
    idx = Index(names, synonyms, parents, obsolete, alt, disease_names, disease_terms, disease_genes,
                version=version, disease_excluded=disease_excluded, freq_unknown=freq_unknown)
    if (d / ZH_FILE).exists():
        idx.zh = _parse_zh(d / ZH_FILE)
    _build_ic(idx)
    try:
        with open(cache, "wb") as fh:
            pickle.dump({"stamp": stamp, "v": INDEX_VERSION, "index": idx}, fh, protocol=pickle.HIGHEST_PROTOCOL)
    except OSError:
        pass
    return idx


def _build_ic(idx: Index) -> None:
    idx.disease_list = sorted(idx.disease_terms)
    closure: Dict[str, Dict[int, float]] = defaultdict(dict)
    for i, did in enumerate(idx.disease_list):
        unknown = idx.freq_unknown.get(did, set())
        best: Dict[str, Tuple[float, bool]] = {}
        for term, w in idx.disease_terms[did].items():
            is_unknown = term in unknown
            if term not in idx.names:
                term = idx.alt.get(term, term)
            for a in idx.ancestors(term):
                k, u = best.get(a, (0.0, False))
                best[a] = (k if is_unknown else max(k, w), u or is_unknown)
        for a, (k, u) in best.items():
            closure[a][i] = _encode(k, u)
    n = len(idx.disease_list)
    idx.term_diseases = dict(closure)
    idx.ic = {t: -math.log(len(ds) / n) for t, ds in closure.items() if ds}
    idx._anc = {}  # not pickled warm: rebuilt lazily


def _norm(text: str) -> str:
    words = re.sub(r"[^a-z0-9]+", " ", text.lower()).split()
    out = []
    for w in words:
        if len(w) > 3 and w.endswith("s") and not w.endswith(("ss", "is", "us")):
            w = w[:-1]
        out.append(w)
    return " ".join(out)


_ZH_PUNCT = re.compile(r"[\s，。、；：,.;:()（）!！?？\"'“”‘’]+")


def _zh_key(text: str) -> str:
    return _ZH_PUNCT.sub("", text)


# What may surround a lay phrase inside a longer query without changing what it says: who it is
# about, how often or how much, sentence particles. Anything else ("舌头大": tongue; "耳聋家族史":
# family history; "无明显抽搐": negated; "听力下降不明显": negated; "眼睑抽搐": another finding)
# means the phrase is part of something else, and it is not used (R-A-P0-2).
_ZH_FILLER = sorted((
    "我家", "我们家", "我的", "我们的", "孩子", "宝宝", "小孩", "患儿", "儿子", "女儿", "闺女", "娃", "他", "她",
    "我", "现在", "最近", "近来", "一直", "总是", "经常", "老是", "还是", "还", "也", "都", "有点儿", "有点",
    "有些", "比较", "很", "特别", "非常", "出现", "有", "就", "会", "了", "的", "啊", "呀", "呢", "吧", "哦", "嘛",
    "明显", "厉害", "严重", "一些", "一点", "得",
), key=len, reverse=True)
_EN_FILLER = {
    "my", "our", "his", "her", "their", "the", "a", "an", "son", "daughter", "child", "kid", "baby", "boy", "girl",
    "he", "she", "they", "is", "was", "are", "were", "has", "had", "have", "been", "very", "quite", "bit", "little",
    "always", "often", "still", "also", "really", "so", "seems", "seem", "to", "be", "year", "month", "old", "aged",
    "age", "at", "since", "and", "now",
}


def _only_filler_zh(text: str) -> bool:
    rest = text
    while rest:
        for f in _ZH_FILLER:
            if rest.startswith(f):
                rest = rest[len(f):]
                break
        else:
            return False
    return True


def lay_lookup(text: str) -> Tuple[List[str], Optional[str]]:
    """The curated HPO ids for a lay phrase, best first, and the LAY_TERMS key that matched.

    An exact match wins. Otherwise the longest lay phrase found inside the query is used
    ("孩子走路晚" -> 走路晚; "my son walked late" -> walked late), but only when everything
    around it is filler (who, how often, particles: _ZH_FILLER/_EN_FILLER). A negation
    ("没有抽搐", "听力下降不明显", "no speech delay") or any other word ("舌头大", "耳聋家族史")
    therefore stops the match: a negated or different finding must never come back as present.
    Pure: needs no HPO files, so the online search can use it too.
    """
    raw = (text or "").strip()
    if not raw:
        return [], None
    if _CJK.search(raw):
        key = _zh_key(raw)
        if key in LAY_TERMS:
            return list(LAY_TERMS[key]), key
        upper = key.upper()
        for k in sorted((k for k in LAY_TERMS if _CJK.search(k) or any(c.isupper() for c in k)),
                        key=len, reverse=True):
            if len(k) < 2:
                continue
            pos = upper.find(k.upper())
            if pos < 0:
                continue
            before, after = key[:pos], key[pos + len(k):]
            if any(before.endswith(neg) or (len(neg) > 1 and neg in before) for neg in _NEGATION_ZH):
                return [], None  # negated: nothing, not a shorter phrase
            if _only_filler_zh(before) and _only_filler_zh(after):
                return list(LAY_TERMS[k]), k
        return [], None
    norm = _norm(raw)
    for cand in (norm, raw.lower()):
        if cand in LAY_TERMS:
            return list(LAY_TERMS[cand]), cand
    words = norm.split()
    best: Optional[Tuple[int, str, int]] = None
    for k in LAY_TERMS:
        if _CJK.search(k):
            continue
        kw = _norm(k).split()
        # a single word inside a sentence is too ambiguous ("fit" in "fit and well"), and a phrase
        # that itself starts with a negation ("no speech") is only taken as the whole query
        if len(kw) < 2 or kw[0] in _NEGATION_EN or kw[0] == "cannot":
            continue
        for i in range(len(words) - len(kw) + 1):
            if words[i:i + len(kw)] == kw and (best is None or len(kw) > best[2]):
                best = (i, k, len(kw))
                break
    if best is None:
        return [], None
    i, k, n = best
    rest = words[:i] + words[i + n:]
    if any(w in _NEGATION_EN for w in words[:i]) or not all(w in _EN_FILLER or w.isdigit() for w in rest):
        return [], None
    return list(LAY_TERMS[k]), k


def negation_warning(text: str) -> Optional[str]:
    """A warning when the phrase negates something (不伴抽搐, 肌张力不高, "no seizures"), unless it is a
    curated lay phrase as a whole (说话不清楚, 不会走路, "not walking"). HPO terms describe findings
    that are present: search can only return the finding itself, which must then be recorded as
    excluded, or not at all when the phrase means "normal"."""
    raw = (text or "").strip()
    ids, key = lay_lookup(raw)
    if ids and key is not None and (_zh_key(raw).upper() == key.upper() if _CJK.search(raw)
                                    else _norm(raw) == _norm(key) or raw.lower() == key.lower()):
        return None
    if _CJK.search(raw):
        hit = next((n for n in _NEGATION_ZH if n in _zh_key(raw)), None)
    else:
        hit = next((w for w in _norm(raw).split() if w in _NEGATION_EN or w in ("cannot", "denies")), None)
    if hit is None:
        return None
    return (f"the phrase contains a negation ({hit!r}): the terms below describe the finding itself. If the "
            "person does NOT have it, record it with status excluded; if the phrase means 'normal', record nothing")


def _lay_hits(idx: Index, text: str) -> List[Dict[str, str]]:
    """Curated lay-phrase matches for `text`, best first, dropping ids not in this release."""
    ids, key = lay_lookup(text)
    out: List[Dict[str, str]] = []
    for tid in ids:
        pid, _note = idx.primary(tid)
        if not pid or pid in {h["id"] for h in out}:
            continue
        out.append({"id": pid, "label": idx.names[pid], "matched": key or text.strip(), "matched_on": "lay phrase",
                    **({"label_zh": idx.zh[pid]} if pid in idx.zh else {})})
    return out


def _lay_is_exact(text: str, hits: List[Dict[str, str]]) -> bool:
    if not hits:
        return False
    key = hits[0]["matched"]
    if _CJK.search(text):
        return _zh_key(text).upper() == str(key).upper()
    return _norm(text) == _norm(str(key)) or text.strip().lower() == str(key).lower()


def search(idx: Index, text: str, limit: int = 10) -> List[Dict[str, str]]:
    """Offline label/synonym search.

    Ranking: a curated lay phrase (LAY_TERMS) first, then exact > prefix > whole words >
    all query words present > substring, after normalising case, punctuation and simple
    plurals. Within a rank the shortest matching wording wins, because a query matched
    inside a much longer label is usually the wrong term ("small head" inside
    "Microcephalic sperm head"), and a match on the term's own label beats the same match
    on a synonym. Non-phenotype terms (onset, inheritance, modifiers) are ranked last.
    `matched_on` says which channel matched: "lay phrase", "label" or "synonym".
    """
    if _CJK.search(text):
        return _search_zh(idx, text, limit)
    lay = _lay_hits(idx, text)
    q = _norm(text)
    if not q:
        return lay[:limit]
    q_words = set(q.split())
    scored = []
    for tid, label in idx.names.items():
        best = None
        for i, cand in enumerate([label] + idx.synonyms.get(tid, [])):
            c = _norm(cand)
            if c == q:
                rank = 0
            elif c.startswith(q + " ") or c.startswith(q):
                rank = 1
            elif re.search(rf"\b{re.escape(q)}\b", c):
                rank = 2
            elif q_words <= set(c.split()):
                rank = 3
            elif q in c:
                rank = 4
            else:
                continue
            # excess length first: how much longer the matching wording is than the query.
            # Then label (0) before synonym (1) at equal length.
            key = (rank, max(0, len(c) - len(q)), 0 if i == 0 else 1, len(c))
            if best is None or key < best[0]:
                best = (key, cand, i)
        if best:
            under_pheno = ROOT_PHENO in idx.ancestors(tid)
            rank, excess, is_syn, _clen = best[0]
            scored.append((rank, 0 if under_pheno else 1, excess, is_syn, len(label), tid, label, best[1], best[2]))
    scored.sort()
    rows = [{"id": s[5], "label": s[6], "matched": s[7], "matched_on": "label" if s[8] == 0 else "synonym",
             **({"label_zh": idx.zh[s[5]]} if s[5] in idx.zh else {})} for s in scored[:limit + len(lay)]]
    # a lay phrase found INSIDE the query yields to an official label or synonym that IS the query
    # ("motor developmental delay" is a synonym of Motor delay; its 'developmental delay' is a lay key)
    exact = [r for r, s in zip(rows, scored) if s[0] == 0] if lay and not _lay_is_exact(text, lay) else []
    return _merge(exact, lay, rows, limit=limit)


def _bigrams(text: str) -> Set[str]:
    return {text[i:i + 2] for i in range(len(text) - 1)} or {text}


# Direction words. A query that says "small/low/less" must not get a label that says
# "increased/large/high" first: 头围小 (small head circumference) put 头围增加 (Increased
# head circumference) first because the two share only the bigram 头围 and the tie was
# broken by length (D-P1-4). Markers are words, not single characters, except at the
# end of a phrase, because 小 and 大 also occur in 小脑 (cerebellum) and 大脑 (cerebrum).
_DOWN = ("减小", "减少", "降低", "下降", "减退", "低下", "缩小", "变小", "偏小", "过小", "不足", "缺乏",
         "减弱", "变短", "过短", "偏低", "过低", "变窄", "变薄")
_UP = ("增加", "增大", "增多", "升高", "增高", "上升", "偏大", "过大", "亢进", "过多", "增宽", "增厚",
       "肥大", "延长", "过长", "偏高", "过高", "增强")
_DOWN_END = ("小", "低", "少", "短", "窄", "薄", "差")
_UP_END = ("大", "高", "多", "长", "宽", "厚")


_POL_NEG = ("不", "未", "无", "没")


def _polarity(text: str) -> int:
    """-1 for 'less/smaller/lower', +1 for 'more/larger/higher', 0 when it says neither (or both).

    A direction word right after a negation says nothing about direction (肌张力不高 is "tone
    not raised", not "raised"), so it is ignored.
    """
    def said(words: Sequence[str], ends: Sequence[str]) -> bool:
        for w in words:
            i = text.find(w)
            while i >= 0:
                if not (i > 0 and text[i - 1] in _POL_NEG):
                    return True
                i = text.find(w, i + 1)
        return any(text.endswith(e) and not (len(text) > len(e) and text[-len(e) - 1] in _POL_NEG) for e in ends)

    down = said(_DOWN, _DOWN_END)
    up = said(_UP, _UP_END)
    return -1 if down and not up else 1 if up and not down else 0


def _search_zh(idx: Index, text: str, limit: int) -> List[Dict[str, str]]:
    """Chinese search over the official Chinese labels.

    Ranking: a curated lay phrase (LAY_TERMS, exact or inside the query) first, then an
    exact label, then any label containing the query, then character-bigram overlap.
    Containment is ordered by how much LONGER the label is than the query -- not by
    prefix-before-substring, which used to put 抽搐样不自主自我拥抱 (self hugging, 10
    characters) above 手足抽搐 (tetany, 4) for the query 抽搐. Within the bigram tier,
    overlap decides, then character overlap (头围减小 shares 头, 围 and 小 with 头围小;
    头围增加 shares two), then the closest length. A label whose direction contradicts
    the query's (增加 for 小, 降低 for 大) is moved behind every label that does not.
    """
    lay = _lay_hits(idx, text)
    if not idx.zh:
        return lay[:limit]
    q = _zh_key(text)
    if not q:
        return lay[:limit]
    qb = _bigrams(q)
    qc = set(q)
    qpol = _polarity(q)
    scored = []
    for tid, zh in idx.zh.items():
        if tid not in idx.names:
            continue
        z = _zh_key(zh)
        if z == q:
            key = (0, 0, 0.0, 0.0)
        elif q in z:
            # one containment tier ordered by excess length; a prefix match only breaks
            # ties between labels of the same length
            key = (1, len(z) - len(q), 0.0 if z.startswith(q) else 0.5, 0.0)
        else:
            zb = _bigrams(z)
            overlap = len(qb & zb) / len(qb | zb)
            if overlap < 0.2:
                continue
            chars = len(qc & set(z)) / len(qc | set(z))
            key = (2, 0, -overlap, -chars)
        opposed = 1 if qpol and _polarity(z) == -qpol else 0
        under_pheno = ROOT_PHENO in idx.ancestors(tid)
        scored.append((opposed, key[0], 0 if under_pheno else 1, key[1], key[2], key[3], abs(len(z) - len(q)),
                       len(z), tid, zh))
    scored.sort()
    rows = [{"id": s[8], "label": idx.names[s[8]], "label_zh": s[9], "matched": s[9], "matched_on": "label_zh"}
            for s in scored[:limit + len(lay)]]
    # an official label that IS the query beats a lay phrase found inside it (语言发育迟缓 is the
    # official label of Delayed speech and language development; 发育迟缓 inside it is a lay key)
    exact = [r for r, s in zip(rows, scored) if s[1] == 0 and s[0] == 0] if lay and not _lay_is_exact(text, lay) else []
    return _merge(exact, lay, rows, limit=limit)


def _merge(*groups: List[Dict[str, str]], limit: int = 10) -> List[Dict[str, str]]:
    out: List[Dict[str, str]] = []
    seen: Set[str] = set()
    for group in groups:
        for row in group:
            if row["id"] not in seen:
                seen.add(row["id"])
                out.append(row)
    return out[:limit]


def resolve_query(idx: Index, present: Sequence[str], excluded: Sequence[str] = ()
                  ) -> Tuple[List[str], List[str], List[str], List[Dict[str, str]]]:
    """(present ids, excluded ids, notes, contradictions) after resolving alt/obsolete ids.

    Non-phenotype terms (onset, inheritance, modifiers) are dropped from the present
    list. An excluded term that is also present, or is an ancestor of a present term
    (seizure excluded, febrile seizure present), contradicts the query: it is dropped
    and reported in `contradictions` (D-P2-1) rather than silently scored.
    """
    notes: List[str] = []
    query: List[str] = []
    for t in present:
        pid, note = idx.primary(t)
        if note:
            notes.append(note)
        if pid and pid not in query:
            query.append(pid)
    excl: List[str] = []
    for t in excluded:
        pid, note = idx.primary(t)
        if note:
            notes.append(note)
        if pid and pid not in excl:
            excl.append(pid)
    if not query:
        raise ValueError("no usable HPO terms in the query")
    under = [t for t in query if ROOT_PHENO not in idx.ancestors(t)]
    for t in under:
        notes.append(f"{t} ({idx.names.get(t)}) is not a phenotypic abnormality term (onset/inheritance/modifier); it is ignored for matching")
    query = [t for t in query if t not in under]
    if not query:
        raise ValueError("no phenotypic abnormality terms left in the query")
    contradictions: List[Dict[str, str]] = []
    kept: List[str] = []
    for e in excl:
        clash = next((q for q in query if q == e or e in idx.ancestors(q)), None)
        if clash is None:
            kept.append(e)
            continue
        why = ("given as both present and excluded" if clash == e else
               f"is an ancestor of the present term {clash} ({idx.names.get(clash)}): a patient with that "
               "finding has this one")
        contradictions.append({"excluded": e, "label": idx.names.get(e) or "", "present": clash, "why": why})
        notes.append(f"excluded {e} ({idx.names.get(e)}) {why}; the exclusion was dropped")
    return query, kept, notes, contradictions


def _not_contradictions(idx: Index, did: str, di: int, query: Sequence[str]) -> Dict[str, Tuple[str, str]]:
    """Curated NOT annotations of `did` that contradict a present query term: NOT-term -> (query term, why).

    A NOT on the patient's term, or on an ancestor of it, contradicts it -- unless the same
    disease also annotates that term, or the patient's term or something below it, positively:
    curations then disagree, and that must not cost the disease its own feature (D-P1-3).
    """
    out: Dict[str, Tuple[str, str]] = {}
    pos = idx.disease_terms.get(did) or {}
    for term, why in (idx.disease_excluded.get(did) or {}).items():
        pterm = idx.alt.get(term, term)
        if pterm in pos or term in pos:
            continue
        for q in query:
            if di in idx.term_diseases.get(q, {}):
                continue  # the disease has q (or a more specific form of it) positively
            if pterm == q or pterm in idx.ancestors(q):
                out[pterm] = (q, why)
                break
    return out


def _score_resnik(idx: Index, query: Sequence[str], excl: Sequence[str], keep: Sequence[int],
                  params: Optional[Dict[str, float]] = None):
    """Per disease: (score, raw, penalty, hits, best[qi], via[qi], contradicted)."""
    nq = len(query)
    best_ic: List[Dict[int, float]] = []
    best_via: List[Dict[int, str]] = []
    for q in query:
        ic_q: Dict[int, float] = {}
        via_q: Dict[int, str] = {}
        for a in sorted(idx.ancestors(q), key=lambda x: idx.ic.get(x, 0.0), reverse=True):
            ic = idx.ic.get(a, 0.0)
            if ic <= 0:
                break
            for di in idx.term_diseases.get(a, ()):
                if di not in ic_q:
                    ic_q[di] = ic
                    via_q[di] = a
        best_ic.append(ic_q)
        best_via.append(via_q)
    # excluded: a disease annotated with the excluded term or below it costs IC(e) x frequency
    # (an annotation that states no frequency counts as 1), summed over distinct terms (F29)
    weight = (params or RESNIK_PARAMS).get("excluded_weight", 1.0)
    min_f = (params or RESNIK_PARAMS).get("resnik_excluded_min_f", 0.0)
    excl_pen: Dict[int, Dict[str, float]] = defaultdict(dict)
    for e in excl:
        ic_e = idx.ic.get(e, 0.0)
        for di, code in idx.term_diseases.get(e, {}).items():
            f = freq_at(code, 1.0)
            if f >= min_f:  # recorded even at weight 0: `excluded_hits` flags the candidate either way
                excl_pen[di][e] = weight * ic_e * f / nq
    out = {}
    for di in keep:
        did = idx.disease_list[di]
        raw = sum(best_ic[qi].get(di, 0.0) for qi in range(nq)) / nq
        per_term = excl_pen.get(di, {})
        contradicted = {}
        not_pen = 0.0
        if did in idx.disease_excluded:
            for pterm, (_q, why) in _not_contradictions(idx, did, di, query).items():
                not_pen += idx.ic.get(pterm, 0.0) / nq
                contradicted[pterm] = why
        excl_part = sum(per_term.values())
        penalty = excl_part + not_pen
        out[di] = (raw - penalty, raw, penalty, sorted(per_term), [best_ic[qi].get(di, 0.0) for qi in range(nq)],
                   [best_via[qi].get(di) for qi in range(nq)], contradicted, excl_part, not_pen)
    return out


def _score_lr(idx: Index, query: Sequence[str], excl: Sequence[str], keep: Sequence[int], params: Dict[str, float]):
    """Per disease: (score, raw, penalty, hits, best[qi], via[qi], contradicted); score = sum of ln LR."""
    f_unknown = params["f_unknown"]
    f_floor = params["f_floor"]
    nq = len(query)
    best_v: List[Dict[int, float]] = []
    best_via: List[Dict[int, str]] = []
    log_floor = math.log(f_floor)
    for q in query:
        v_q: Dict[int, float] = {}
        via_q: Dict[int, str] = {}
        for a in sorted(idx.ancestors(q), key=lambda x: idx.ic.get(x, 0.0), reverse=True):
            ic = idx.ic.get(a, 0.0)
            if ic <= 0:
                break
            for di, code in idx.term_diseases.get(a, {}).items():
                cur = v_q.get(di)
                if cur is not None and cur >= ic:
                    continue  # ln f <= 0, so a less informative ancestor cannot beat it
                f = freq_at(code, f_unknown)
                val = ic + (math.log(f) if f > f_floor else log_floor)
                if cur is None or val > cur:
                    v_q[di] = val
                    via_q[di] = a
        best_v.append(v_q)
        best_via.append(via_q)
    cap = params["excluded_cap"]
    min_f = params["excluded_min_f"]
    excl_pen: Dict[int, Dict[str, float]] = defaultdict(dict)
    for e in excl:
        for di, code in idx.term_diseases.get(e, {}).items():
            f = freq_at(code, f_unknown)
            if f < min_f or f <= 0:
                continue
            pe = min(cap, f)
            excl_pen[di][e] = -math.log(1.0 - pe)
    unexplained = params["unexplained"]
    not_pen = params["not_penalty"]
    out = {}
    for di in keep:
        did = idx.disease_list[di]
        per_q = [best_v[qi].get(di, unexplained) for qi in range(nq)]
        contradicted = {}
        not_cost = 0.0
        if did in idx.disease_excluded:
            for pterm, (qterm, why) in _not_contradictions(idx, did, di, query).items():
                qi = query.index(qterm)
                if per_q[qi] > not_pen:
                    not_cost += per_q[qi] - not_pen
                    per_q[qi] = not_pen
                contradicted[pterm] = why
        raw = sum(per_q) + not_cost  # before the NOT correction, so penalty = excluded + NOT
        per_term = excl_pen.get(di, {})
        excl_part = sum(per_term.values())
        penalty = excl_part + not_cost
        out[di] = (raw - penalty, raw, penalty, sorted(per_term), per_q,
                   [best_via[qi].get(di) for qi in range(nq)], contradicted, excl_part, not_cost)
    return out


def rank(idx: Index, present: Sequence[str], excluded: Sequence[str] = (), top: int = 20,
         db: Sequence[str] = ("OMIM", "ORPHA"), method: Optional[str] = None,
         params: Optional[Dict[str, float]] = None) -> Dict[str, object]:
    method = method or DEFAULT_METHOD
    if method not in METHODS:
        raise ValueError(f"method must be one of {', '.join(METHODS)}")
    top = max(1, int(top))
    excluded = list(excluded or ())
    lr_params = dict(LR_PARAMS)
    resnik_params = dict(RESNIK_PARAMS)
    for k, v in (params or {}).items():
        try:
            v = float(v)
        except (TypeError, ValueError):
            raise ValueError(f"ranking parameter {k!r} must be a number, got {v!r}") from None
        lo, hi = _PARAM_RANGE.get(k, (None, None))
        if lo is None:
            raise ValueError(f"unknown ranking parameter {k!r}; known: {', '.join(sorted(LR_PARAMS) + sorted(RESNIK_PARAMS))}")
        if not (lo <= v <= hi):
            raise ValueError(f"ranking parameter {k!r} must be between {lo} and {hi}, got {v}")
        (lr_params if k in lr_params else resnik_params)[k] = v
    query, excl, notes, contradictions = resolve_query(idx, present, excluded)
    nq = len(query)
    keep = [i for i, d in enumerate(idx.disease_list) if d.split(":")[0] in db]
    # the best a query term can score: the IC of its most informative ANNOTATED ancestor (a term no
    # disease is annotated with has IC 0 while its ancestors still match, which made `relative` > 1)
    reach = [max((idx.ic.get(a, 0.0) for a in idx.ancestors(q)), default=0.0) for q in query]
    if method == "lr":
        scored = _score_lr(idx, query, excl, keep, lr_params)
        max_possible = sum(reach)
    else:
        scored = _score_resnik(idx, query, excl, keep, resnik_params)
        max_possible = sum(reach) / nq

    results = []
    spec_cache: Dict[int, float] = {}
    for di, (adj, raw, penalty, hits, per_q, via, contradicted, _ep, _np) in scored.items():
        # F35: a principled tie-break, so equal scores are not resolved by disease id.
        # (1) how many query terms matched the disease's own term exactly rather than
        # through a common ancestor; (2) the disease's annotation specificity, the mean IC
        # of the terms it is annotated with -- a tightly specified disorder is a more
        # informative explanation than one annotated only with generic features.
        exact_hits = sum(1 for qi in range(nq) if via[qi] == query[qi])
        results.append((adj, raw, penalty, di, hits, exact_hits, contradicted))

    def specificity(di: int) -> float:
        if di not in spec_cache:
            ann = [idx.ic.get(idx.alt.get(t, t), 0.0) for t in idx.disease_terms[idx.disease_list[di]]]
            spec_cache[di] = (sum(ann) / len(ann)) if ann else 0.0
        return spec_cache[di]

    freq_cache: Dict[int, float] = {}

    def matched_frequency(di: int) -> float:
        """Mean frequency of the disease's annotations behind each matched query term (stated
        frequencies; none stated counts as 1, as in the excluded penalty)."""
        if di not in freq_cache:
            via = scored[di][5]
            vals = [freq_at(idx.term_diseases.get(a, {}).get(di, 0.0), 1.0) if a else 0.0 for a in via]
            freq_cache[di] = sum(vals) / len(vals) if vals else 0.0
        return freq_cache[di]

    # the later keys are only needed where the first two tie; compute them lazily
    use_freq = method == "resnik" and resnik_params.get("tiebreak_frequency", 0.0) > 0
    results = _sort_lazy(results, specificity, idx, matched_frequency if use_freq else None)
    # how many diseases share the top adjusted score (1 = no tie); the caller needs to
    # know when the order it is reading was decided by a tie-break rather than the score
    top_score = round(results[0][0], 9) if results else None
    ties = sum(1 for r in results if round(r[0], 9) == top_score) if results else 0
    tie_counts: Dict[float, int] = defaultdict(int)
    for r in results:
        tie_counts[round(r[0], 9)] += 1

    out = []
    for adj, raw, penalty, di, hits, exact_hits, contradicted in results[:top]:
        did = idx.disease_list[di]
        _a, _r, _p, _h, per_q, via, _c, excl_part, not_part = scored[di]
        matched = []
        for qi, q in enumerate(query):
            mica = via[qi]
            row = {"query": q, "query_label": idx.names.get(q),
                   "via": mica, "via_label": idx.names.get(mica) if mica else None,
                   "exact": mica == q, ("ln_lr" if method == "lr" else "ic"): round(per_q[qi], 3)}
            matched.append(row)
        out.append({
            "disease": did, "name": idx.disease_names.get(did, did), "score": round(adj, 4),
            "relative": round(adj / max_possible, 3) if max_possible else None,
            "excluded_hits": hits, "penalty": round(penalty, 4),
            "excluded_penalty": round(excl_part, 4), "curated_not_penalty": round(not_part, 4),
            "contradicted_by_curation": [{"term": t, "label": idx.names.get(t), "source": why}
                                         for t, why in sorted(contradicted.items())],
            "exact_matches": exact_hits, "annotation_specificity": round(specificity(di), 3),
            "tied_at_this_score": tie_counts[round(adj, 9)],
            "genes": idx.disease_genes.get(did, []), "matches": matched,
        })
    # a gene scores as its best disease; ties are broken by its support across the next-best diseases
    gene_scores: Dict[str, Tuple[float, str]] = {}
    gene_support: Dict[str, float] = defaultdict(float)
    for rank_i, (adj, raw, penalty, di, hits, _ex, _contra) in enumerate(results):
        if rank_i >= 200 or (method == "resnik" and adj <= 0):
            break
        did = idx.disease_list[di]
        for g in idx.disease_genes.get(did, []):
            gene_support[g] += adj if method == "resnik" else 1.0
            if g not in gene_scores or gene_scores[g][0] < adj:
                gene_scores[g] = (adj, did)
    genes = sorted(gene_scores.items(), key=lambda kv: (kv[1][0], gene_support[kv[0]]), reverse=True)[:top]
    if method == "lr":
        method_text = ("likelihood-ratio score in the spirit of LIRICAL (not LIRICAL): per present term "
                       "ln LR = max over the disease's annotations of [IC(common ancestor) + ln frequency]; per "
                       "excluded term ln(1 - frequency of the disease's annotations at or below it, capped at "
                       f"{lr_params['excluded_cap']}); curated NOT on a present term costs {lr_params['not_penalty']}. "
                       "Score = sum of ln LR. Equal scores are broken by exactly matched query terms, then annotation "
                       "specificity, then disease id.")
    else:
        w = resnik_params.get("excluded_weight", 1.0)
        excl_text = ("Excluded terms are not scored (weight 0): a candidate annotated with an excluded feature is "
                     "flagged in `excluded_hits` instead, because on the phenopacket-store benchmark penalising them "
                     "cost accuracy at every number of excluded terms (docs/BENCHMARK.md). " if w == 0 else
                     f"Excluded terms are penalised by {w:g} x the SUM of their IC x annotation frequency over distinct "
                     "terms (bounded by sum(IC of the excluded terms) / number of query terms). ")
        method_text = ("Resnik best-match average (query→disease), IC from HPO disease annotations. " + excl_text +
                       "An HPOA `NOT` or \"Excluded\" "
                       "(HP:0040285) annotation that contradicts a present query term costs that term's IC "
                       "(a case count of 0/n is not a NOT). A negative score means the curated phenotype contradicts "
                       "the query. Equal scores are broken by the number of exactly "
                       "matched query terms, then by the disease's annotation specificity (mean IC of its annotated "
                       "terms), then by disease id for determinism — see `ties`.")
    return {
        "method": method_text,
        "method_id": method,
        "ties": ties,
        "tie_break": "exact query-term matches, then annotation specificity, then disease id",
        "hpo_version": idx.version,
        "query": [{"id": q, "label": idx.names.get(q), **({"label_zh": idx.zh[q]} if q in idx.zh else {})} for q in query],
        "excluded": [{"id": e, "label": idx.names.get(e)} for e in excl],
        "contradictions": contradictions,
        "diseases_scored": len(keep),
        "max_possible": round(max_possible, 4),
        "diseases": out,
        "genes": [{"gene": g, "score": round(s, 4), "support": round(gene_support[g], 3), "via": d,
                   "via_name": idx.disease_names.get(d)} for g, (s, d) in genes],
        "notes": notes,
    }


def _sort_lazy(rows, specificity, idx, frequency=None):
    """Sort by (score, exact matches, [matched frequency,] specificity, id) descending; the later keys
    are computed only inside ties."""
    rows = sorted(rows, key=lambda r: (round(r[0], 9), r[5]), reverse=True)
    out = []
    i = 0
    while i < len(rows):
        j = i + 1
        key = (round(rows[i][0], 9), rows[i][5])
        while j < len(rows) and (round(rows[j][0], 9), rows[j][5]) == key:
            j += 1
        block = rows[i:j]
        if len(block) > 1:
            if frequency is not None:
                block.sort(key=lambda r: (round(frequency(r[3]), 9), round(specificity(r[3]), 9),
                                          idx.disease_list[r[3]]), reverse=True)
            else:
                block.sort(key=lambda r: (round(specificity(r[3]), 9), idx.disease_list[r[3]]), reverse=True)
        out.extend(block)
        i = j
    return out
