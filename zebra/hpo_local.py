"""Offline phenotype-driven ranking from the HPO release files.

Data (fetched once by `zebra hpo fetch`, ~80 MB, into ~/.cache/zebra-mod/hpo):
  hp.json                 the ontology (obographs JSON)
  phenotype.hpoa          disease -> phenotype annotations (OMIM, ORPHA, DECIPHER)
  genes_to_phenotype.txt  gene -> phenotype, with the disease each pair comes from

Method: information content IC(t) = -ln(share of diseases annotated with t
or a descendant); a query term's match to a disease is the IC of their most
informative common ancestor (Resnik); a disease's score is the mean over the
query terms of their best match (query -> disease best-match average, as the
Phenomizer does). Each excluded term the disease is annotated with (directly or
below) costs its IC times the annotation frequency, and those costs are SUMMED
over distinct excluded terms. Diseases whose
curated annotation carries HPOA's `NOT` qualifier, or the "Excluded" (0%)
frequency class, for a term the patient is said to HAVE are penalised the same
way. Scores are reported raw and as a share of the query's own maximum (mean IC
of the query terms); `ties` says how many diseases share the top score.

This is a transparent baseline, not LIRICAL or Exomiser: no likelihood ratios,
no genotype term. It runs offline, reproducibly, from public files.
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
from typing import Dict, List, Optional, Sequence, Set, Tuple

from zebra.http import cache_dir

FILES = ("hp.json", "phenotype.hpoa", "genes_to_phenotype.txt")
RELEASE_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/{name}"
# official Chinese labels (HPO translation project); optional, enables Chinese search
ZH_FILE = "hp-zh.babelon.tsv"
ZH_URL = "https://raw.githubusercontent.com/obophenotype/hpo-translations/main/babelon/hp-zh.babelon.tsv"
_CJK = re.compile(r"[\u3400-\u9fff]")
ROOT_PHENO = "HP:0000118"  # Phenotypic abnormality

# Lay phrases families actually use, each mapped to the official HPO term(s) a clinician
# would record, best match first. This is a translation table, not a guess: every id below
# was checked against the HPO release's own labels, and
# tests/test_hpo_local.py::test_f36_lay_terms_point_at_real_ids re-checks them against the
# live release files so a mapping cannot silently rot across HPO versions.
# Needed because the ontology's labels are clinical: "走路晚" (walked late) matches nothing
# literally, and "智力低下" is closer by character bigrams to 面部肌张力低下 (facial
# hypotonia) than to its own official label 智力障碍.
LAY_TERMS: Dict[str, List[str]] = {
    # --- Chinese: motor milestones
    "走路晚": ["HP:0031936", "HP:0001270"],      # Delayed ability to walk; Motor delay
    "学步晚": ["HP:0031936", "HP:0001270"],
    "不会走路": ["HP:0031936", "HP:0001270"],
    "不会走": ["HP:0031936", "HP:0001270"],
    "翻身晚": ["HP:0032989", "HP:0001270"],      # Delayed ability to roll over
    "不会翻身": ["HP:0032989", "HP:0001270"],
    "抬头晚": ["HP:0002421", "HP:0032988"],      # Poor head control; Persistent head lag
    "软趴趴": ["HP:0008947", "HP:0001252"],      # Floppy infant; Hypotonia
    "没力气": ["HP:0001324", "HP:0001252"],      # Muscle weakness; Hypotonia
    "肌无力": ["HP:0001324"],
    # --- Chinese: seizures
    "抽风": ["HP:0001250"],                       # Seizure
    "抽搐": ["HP:0001250"],
    "惊厥": ["HP:0001250"],
    "发热惊厥": ["HP:0002373"],                   # Febrile seizure
    "高热惊厥": ["HP:0002373"],
    "发烧抽筋": ["HP:0002373"],
    "发烧抽搐": ["HP:0002373"],
    # --- Chinese: speech and cognition
    "不会说话": ["HP:0001344", "HP:0000750"],     # Absent speech; Delayed speech and language development
    "不说话": ["HP:0001344"],
    "说话晚": ["HP:0000750"],
    "智力低下": ["HP:0001249"],                   # Intellectual disability
    "智力落后": ["HP:0001249"],
    "智力障碍": ["HP:0001249"],
    "发育迟缓": ["HP:0001263"],                   # Global developmental delay
    "发育落后": ["HP:0001263"],
    # --- Chinese: growth, head, senses, feeding
    "不长个": ["HP:0004322"],                     # Short stature
    "长不高": ["HP:0004322"],
    "个子矮": ["HP:0004322"],
    "头小": ["HP:0000252"],                       # Microcephaly
    "听不见": ["HP:0000365"],                     # Hearing impairment
    "耳朵听不见": ["HP:0000365"],
    "喂奶困难": ["HP:0008872"],                   # Feeding difficulties in infancy
    "吃奶费劲": ["HP:0008872"],
    "吞咽困难": ["HP:0002015"],                   # Dysphagia
    # --- English lay phrasing
    "floppy baby": ["HP:0008947", "HP:0001252"],
    "floppy infant": ["HP:0008947"],
    "late walker": ["HP:0031936", "HP:0001270"],
    "walking late": ["HP:0031936", "HP:0001270"],
    "walked late": ["HP:0031936", "HP:0001270"],
    "not walking": ["HP:0031936", "HP:0001270"],
    "not talking": ["HP:0001344", "HP:0000750"],
    "no speech": ["HP:0001344"],
    "cannot talk": ["HP:0001344"],
    "speech delay": ["HP:0000750"],
    "fits": ["HP:0001250"],
    "convulsion": ["HP:0001250"],
    "convulsions": ["HP:0001250"],
    "febrile convulsion": ["HP:0002373"],
    "febrile convulsions": ["HP:0002373"],
    "small head": ["HP:0000252"],
    "mental retardation": ["HP:0001249"],
    "slow learner": ["HP:0001249"],
    "developmental delay": ["HP:0001263"],
    "weak muscles": ["HP:0001324", "HP:0001252"],
    "poor feeding": ["HP:0008872"],
    "hard of hearing": ["HP:0000365"],
}
_PURL = re.compile(r"^http://purl\.obolibrary\.org/obo/HP_(\d{7})$")
INDEX_VERSION = 5  # bumped: the index now carries disease_excluded (F29)


class HpoDataMissing(Exception):
    pass


def data_dir() -> Path:
    return Path(os.environ.get("ZEBRA_HPO_DIR") or (cache_dir() / "hpo"))


def missing_files() -> List[str]:
    return [f for f in FILES if not (data_dir() / f).exists() or (data_dir() / f).stat().st_size == 0]


def _curie(iri: str) -> Optional[str]:
    m = _PURL.match(iri)
    return f"HP:{m.group(1)}" if m else (iri if iri.startswith("HP:") else None)


def _freq_weight(raw: str) -> float:
    """phenotype.hpoa frequency → a weight in (0, 1]."""
    raw = (raw or "").strip()
    named = {
        "HP:0040280": 1.0, "HP:0040281": 0.9, "HP:0040282": 0.55, "HP:0040283": 0.17, "HP:0040284": 0.025,
        "HP:0040285": 0.0,
    }
    if raw in named:
        return named[raw]
    m = re.match(r"^(\d+)/(\d+)$", raw)
    if m and int(m.group(2)) > 0:
        return int(m.group(1)) / int(m.group(2))
    m = re.match(r"^([\d.]+)%$", raw)
    if m:
        return float(m.group(1)) / 100
    return 1.0


@dataclass
class Index:
    names: Dict[str, str]
    synonyms: Dict[str, List[str]]
    parents: Dict[str, List[str]]
    obsolete: Dict[str, Optional[str]]  # id -> replaced_by
    alt: Dict[str, str]  # alt id -> primary
    disease_names: Dict[str, str]
    disease_terms: Dict[str, Dict[str, float]]  # disease -> term -> frequency weight
    disease_genes: Dict[str, List[str]]
    ic: Dict[str, float] = field(default_factory=dict)
    term_diseases: Dict[str, Set[int]] = field(default_factory=dict)
    disease_list: List[str] = field(default_factory=list)
    version: str = ""
    zh: Dict[str, str] = field(default_factory=dict)
    # Curated negative evidence, kept rather than discarded: disease -> term -> why.
    # Two channels in phenotype.hpoa say "this disease does NOT have this feature" --
    # the `NOT` qualifier, and the "Excluded" (HP:0040285, 0%) frequency class. Both
    # are the strongest available evidence against a candidate and used to be dropped.
    disease_excluded: Dict[str, Dict[str, str]] = field(default_factory=dict)

    def ancestors(self, term: str) -> Set[str]:
        seen: Set[str] = set()
        stack = [term]
        while stack:
            t = stack.pop()
            if t in seen:
                continue
            seen.add(t)
            stack.extend(self.parents.get(t, ()))
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
    disease_names: Dict[str, str] = {}
    disease_terms: Dict[str, Dict[str, float]] = defaultdict(dict)
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
            # Both channels of curated negative evidence are kept in their own map:
            # the `NOT` qualifier states the disease does not have the feature, and
            # frequency class HP:0040285 ("Excluded", 0%) says the same numerically.
            if (row.get("qualifier") or "").upper() == "NOT":
                disease_excluded[did][term] = "HPOA NOT qualifier"
                continue
            w = _freq_weight(row.get("frequency", ""))
            if w <= 0:
                disease_excluded[did][term] = f"HPOA frequency class {row.get('frequency', '')!r} (0%, Excluded)"
                continue
            disease_terms[did][term] = max(w, disease_terms[did].get(term, 0.0))
    return disease_names, dict(disease_terms), {d: dict(t) for d, t in disease_excluded.items()}


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
    disease_names, disease_terms, disease_excluded = _parse_hpoa(d / "phenotype.hpoa")
    disease_genes = _parse_g2p(d / "genes_to_phenotype.txt")
    idx = Index(names, synonyms, parents, obsolete, alt, disease_names, disease_terms, disease_genes,
                version=version, disease_excluded=disease_excluded)
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
    anc_cache: Dict[str, Set[str]] = {}
    term_diseases: Dict[str, Set[int]] = defaultdict(set)
    for i, did in enumerate(idx.disease_list):
        closure: Set[str] = set()
        for term in idx.disease_terms[did]:
            if term not in idx.names:
                term = idx.alt.get(term, term)
            if term not in anc_cache:
                anc_cache[term] = idx.ancestors(term)
            closure |= anc_cache[term]
        for t in closure:
            term_diseases[t].add(i)
    n = len(idx.disease_list)
    idx.term_diseases = dict(term_diseases)
    idx.ic = {t: -math.log(len(ds) / n) for t, ds in term_diseases.items() if ds}


def _norm(text: str) -> str:
    words = re.sub(r"[^a-z0-9]+", " ", text.lower()).split()
    out = []
    for w in words:
        if len(w) > 3 and w.endswith("s") and not w.endswith(("ss", "is", "us")):
            w = w[:-1]
        out.append(w)
    return " ".join(out)


def _lay_hits(idx: Index, text: str) -> List[Dict[str, str]]:
    """Curated lay-phrase matches for `text`, best first, dropping ids not in this release."""
    key = _norm(text) if not _CJK.search(text) else re.sub(r"[\s，。、；：,.;:()（）]+", "", text)
    ids = LAY_TERMS.get(key)
    if ids is None and not _CJK.search(text):
        ids = LAY_TERMS.get(text.strip().lower())
    out: List[Dict[str, str]] = []
    for tid in ids or []:
        pid, _note = idx.primary(tid)
        if not pid or pid in {h["id"] for h in out}:
            continue
        out.append({"id": pid, "label": idx.names[pid], "matched": text.strip(), "matched_on": "lay phrase",
                    **({"label_zh": idx.zh[pid]} if pid in idx.zh else {})})
    return out


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
    out = _lay_hits(idx, text)
    seen = {h["id"] for h in out}
    q = _norm(text)
    if not q:
        return out[:limit]
    q_words = set(q.split())
    scored = []
    for tid, label in idx.names.items():
        if tid in seen:
            continue
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
    for s in scored[:max(0, limit - len(out))]:
        out.append({"id": s[5], "label": s[6], "matched": s[7],
                    "matched_on": "label" if s[8] == 0 else "synonym",
                    **({"label_zh": idx.zh[s[5]]} if s[5] in idx.zh else {})})
    return out[:limit]


def _bigrams(text: str) -> Set[str]:
    return {text[i:i + 2] for i in range(len(text) - 1)} or {text}


def _search_zh(idx: Index, text: str, limit: int) -> List[Dict[str, str]]:
    """Chinese search over the official Chinese labels.

    Ranking: a curated lay phrase (LAY_TERMS) first, then an exact label, then any label
    containing the query, then character-bigram overlap. Containment is ordered by how
    much LONGER the label is than the query -- not by prefix-before-substring, which
    used to put 抽搐样不自主自我拥抱 (self hugging, 10 characters) above 手足抽搐
    (tetany, 4) for the query 抽搐. Within the bigram tier, overlap decides and the
    closest length breaks ties.
    """
    out = _lay_hits(idx, text)
    seen = {h["id"] for h in out}
    if not idx.zh:
        return out[:limit]
    q = re.sub(r"[\s，。、；：,.;:()（）]+", "", text)
    if not q:
        return out[:limit]
    qb = _bigrams(q)
    scored = []
    for tid, zh in idx.zh.items():
        if tid not in idx.names or tid in seen:
            continue
        z = re.sub(r"[\s，。、；：,.;:()（）]+", "", zh)
        if z == q:
            key = (0, 0, 0.0)
        elif q in z:
            # one containment tier ordered by excess length; a prefix match only breaks
            # ties between labels of the same length
            key = (1, len(z) - len(q), 0.0 if z.startswith(q) else 0.5)
        else:
            overlap = len(qb & _bigrams(z)) / len(qb | _bigrams(z))
            if overlap < 0.2:
                continue
            key = (2, 0, -overlap)
        under_pheno = ROOT_PHENO in idx.ancestors(tid)
        scored.append((key[0], 0 if under_pheno else 1, key[1], key[2], abs(len(z) - len(q)), len(z), tid, zh))
    scored.sort()
    for s in scored[:max(0, limit - len(out))]:
        out.append({"id": s[6], "label": idx.names[s[6]], "label_zh": s[7], "matched": s[7],
                    "matched_on": "label_zh"})
    return out[:limit]


def rank(idx: Index, present: Sequence[str], excluded: Sequence[str] = (), top: int = 20,
         db: Sequence[str] = ("OMIM", "ORPHA")) -> Dict[str, object]:
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
        if pid:
            excl.append(pid)
    if not query:
        raise ValueError("no usable HPO terms in the query")
    under = [t for t in query if ROOT_PHENO not in idx.ancestors(t)]
    for t in under:
        notes.append(f"{t} ({idx.names.get(t)}) is not a phenotypic abnormality term (onset/inheritance/modifier); it is ignored for matching")
    query = [t for t in query if t not in under]
    if not query:
        raise ValueError("no phenotypic abnormality terms left in the query")

    n_dis = len(idx.disease_list)
    keep = [i for i, d in enumerate(idx.disease_list) if d.split(":")[0] in db]
    best = [[0.0] * len(query) for _ in range(n_dis)]
    best_term: List[List[Optional[str]]] = [[None] * len(query) for _ in range(n_dis)]
    for qi, q in enumerate(query):
        ancestors = sorted(idx.ancestors(q), key=lambda a: idx.ic.get(a, 0.0), reverse=True)
        assigned: Set[int] = set()
        for a in ancestors:
            ic = idx.ic.get(a, 0.0)
            if ic <= 0:
                break
            for di in idx.term_diseases.get(a, ()):
                if di not in assigned:
                    assigned.add(di)
                    best[di][qi] = ic
                    best_term[di][qi] = a
    max_possible = sum(idx.ic.get(q, 0.0) for q in query) / len(query)

    results = []
    for di in keep:
        did = idx.disease_list[di]
        score = sum(best[di]) / len(query)
        # The docstring has always described the penalty as a sum over excluded terms
        # ("Each excluded term ... costs its IC times the annotation frequency"), and the
        # sum is the right reading: each excluded finding the disease is annotated with is
        # an independent contradiction, so five of them must cost more than one. The code
        # took the maximum, discarding most of the force of the one lever a clinician has
        # for pushing a wrong candidate down. Fixed here rather than in the docstring.
        # No artificial ceiling is put on the sum: it is already bounded by
        # sum(IC(e) for e in excluded) / len(query), since each term contributes at most
        # its own IC at full annotation weight and is counted once. A candidate
        # contradicted by five explicitly-absent features should fall below everything,
        # and a negative adjusted score says exactly that.
        per_term: Dict[str, float] = {}
        if excl:
            for term, w in idx.disease_terms[did].items():
                anc = idx.ancestors(idx.alt.get(term, term))
                for e in excl:
                    if e in anc:  # disease has the excluded term or something below it
                        p = idx.ic.get(e, 0.0) * w / len(query)
                        per_term[e] = max(per_term.get(e, 0.0), p)
        # Curated negative evidence the other way round: the disease is annotated NOT to
        # have a term the patient is said to HAVE. The weight is the term's own IC on the
        # same scale, at full annotation weight, because a curator asserted it.
        contradicted: Dict[str, str] = {}
        for term, why in (idx.disease_excluded.get(did) or {}).items():
            pterm = idx.alt.get(term, term)
            for q in query:
                if pterm == q or pterm in idx.ancestors(q):
                    per_term[f"not:{pterm}"] = max(per_term.get(f"not:{pterm}", 0.0),
                                                   idx.ic.get(pterm, 0.0) / len(query))
                    contradicted[pterm] = why
        penalty = sum(per_term.values()) if per_term else 0.0
        hits = sorted(e for e in per_term if not e.startswith("not:"))
        # F35: a principled tie-break, so equal scores are not resolved by disease id.
        # (1) how many query terms matched the disease's own term exactly rather than
        # through a common ancestor; (2) the disease's annotation specificity, the mean IC
        # of the terms it is annotated with -- a tightly specified disorder is a more
        # informative explanation than one annotated only with generic features.
        exact_hits = sum(1 for qi in range(len(query)) if best_term[di][qi] == query[qi])
        ann = [idx.ic.get(idx.alt.get(t, t), 0.0) for t in idx.disease_terms[did]]
        specificity = (sum(ann) / len(ann)) if ann else 0.0
        results.append((score - penalty, score, penalty, di, hits, exact_hits, specificity, contradicted))
    results.sort(key=lambda r: (round(r[0], 9), r[5], round(r[6], 9), idx.disease_list[r[3]]), reverse=True)
    # how many diseases share the top adjusted score (1 = no tie); the caller needs to
    # know when the order it is reading was decided by a tie-break rather than the score
    top_score = round(results[0][0], 9) if results else None
    ties = sum(1 for r in results if round(r[0], 9) == top_score) if results else 0

    out = []
    for adj, raw, penalty, di, hits, exact_hits, specificity, contradicted in results[:top]:
        did = idx.disease_list[di]
        matched = []
        for qi, q in enumerate(query):
            mica = best_term[di][qi]
            matched.append({
                "query": q, "query_label": idx.names.get(q),
                "via": mica, "via_label": idx.names.get(mica) if mica else None,
                "exact": mica == q, "ic": round(best[di][qi], 3),
            })
        out.append({
            "disease": did, "name": idx.disease_names.get(did, did), "score": round(adj, 4),
            "relative": round(adj / max_possible, 3) if max_possible else None,
            "excluded_hits": hits, "penalty": round(penalty, 4),
            "contradicted_by_curation": [{"term": t, "label": idx.names.get(t), "source": why}
                                         for t, why in sorted(contradicted.items())],
            "exact_matches": exact_hits, "annotation_specificity": round(specificity, 3),
            "tied_at_this_score": sum(1 for r in results if round(r[0], 9) == round(adj, 9)),
            "genes": idx.disease_genes.get(did, []), "matches": matched,
        })
    # a gene scores as its best disease; ties are broken by its support across the next-best diseases
    gene_scores: Dict[str, Tuple[float, str]] = {}
    gene_support: Dict[str, float] = defaultdict(float)
    for rank_i, (adj, raw, penalty, di, hits, _ex, _spec, _contra) in enumerate(results):
        if adj <= 0 or rank_i >= 200:
            break
        did = idx.disease_list[di]
        for g in idx.disease_genes.get(did, []):
            gene_support[g] += adj
            if g not in gene_scores or gene_scores[g][0] < adj:
                gene_scores[g] = (adj, did)
    genes = sorted(gene_scores.items(), key=lambda kv: (kv[1][0], gene_support[kv[0]]), reverse=True)[:top]
    return {
        "method": ("Resnik best-match average (query→disease), IC from HPO disease annotations. Excluded terms are "
                   "penalised by the SUM of their IC x annotation frequency over distinct terms (bounded by "
                   "sum(IC of the excluded terms) / number of query terms); HPOA `NOT` and \"Excluded\" (0%) "
                   "annotations that contradict a present query term are penalised the same way. A negative "
                   "score means the curated phenotype contradicts the query. Equal scores are broken by the number of exactly "
                   "matched query terms, then by the disease's annotation specificity (mean IC of its annotated "
                   "terms), then by disease id for determinism — see `ties`."),
        "ties": ties,
        "tie_break": "exact query-term matches, then annotation specificity, then disease id",
        "hpo_version": idx.version,
        "query": [{"id": q, "label": idx.names.get(q), **({"label_zh": idx.zh[q]} if q in idx.zh else {})} for q in query],
        "excluded": [{"id": e, "label": idx.names.get(e)} for e in excl],
        "diseases_scored": len(keep),
        "max_possible": round(max_possible, 4),
        "diseases": out,
        "genes": [{"gene": g, "score": round(s, 4), "support": round(gene_support[g], 3), "via": d,
                   "via_name": idx.disease_names.get(d)} for g, (s, d) in genes],
        "notes": notes,
    }
