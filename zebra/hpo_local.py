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
below) costs its IC times the annotation frequency. Scores are reported raw and
as a share of the query's own maximum (mean IC of the query terms).

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
_PURL = re.compile(r"^http://purl\.obolibrary\.org/obo/HP_(\d{7})$")
INDEX_VERSION = 3


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
    with open(path, encoding="utf-8") as fh:
        rows = (line for line in fh if not line.startswith("#"))
        reader = csv.DictReader(rows, delimiter="\t")
        for row in reader:
            if row.get("aspect") != "P" or (row.get("qualifier") or "").upper() == "NOT":
                continue
            did = row["database_id"]
            disease_names[did] = row.get("disease_name", did)
            w = _freq_weight(row.get("frequency", ""))
            if w <= 0:
                continue
            term = row["hpo_id"]
            disease_terms[did][term] = max(w, disease_terms[did].get(term, 0.0))
    return disease_names, dict(disease_terms)


def _parse_g2p(path: Path) -> Dict[str, List[str]]:
    genes: Dict[str, Set[str]] = defaultdict(set)
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t")
        for row in reader:
            did = row.get("disease_id")
            sym = row.get("gene_symbol")
            if did and sym:
                genes[did].add(sym)
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
    disease_names, disease_terms = _parse_hpoa(d / "phenotype.hpoa")
    disease_genes = _parse_g2p(d / "genes_to_phenotype.txt")
    idx = Index(names, synonyms, parents, obsolete, alt, disease_names, disease_terms, disease_genes, version=version)
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


def search(idx: Index, text: str, limit: int = 10) -> List[Dict[str, str]]:
    """Offline label/synonym search: exact > prefix > whole words > substring, after normalising
    case, punctuation and simple plurals; ties go to the shorter (more general) wording."""
    if _CJK.search(text):
        return _search_zh(idx, text, limit)
    q = _norm(text)
    if not q:
        return []
    q_words = set(q.split())
    scored = []
    for tid, label in idx.names.items():
        best = None
        for cand in [label] + idx.synonyms.get(tid, []):
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
            key = (rank, len(c))
            if best is None or key < best[0]:
                best = (key, cand)
        if best:
            under_pheno = ROOT_PHENO in idx.ancestors(tid)
            scored.append((best[0][0], 0 if under_pheno else 1, best[0][1], len(label), tid, label, best[1]))
    scored.sort()
    return [{"id": s[4], "label": s[5], "matched": s[6], **({"label_zh": idx.zh[s[4]]} if s[4] in idx.zh else {})}
            for s in scored[:limit]]


def _bigrams(text: str) -> Set[str]:
    return {text[i:i + 2] for i in range(len(text) - 1)} or {text}


def _search_zh(idx: Index, text: str, limit: int) -> List[Dict[str, str]]:
    """Chinese search over the official Chinese labels: exact > prefix > contains > character-bigram overlap."""
    if not idx.zh:
        return []
    q = re.sub(r"[\s，。、；：,.;:()（）]+", "", text)
    qb = _bigrams(q)
    scored = []
    for tid, zh in idx.zh.items():
        if tid not in idx.names:
            continue
        z = re.sub(r"[\s，。、；：,.;:()（）]+", "", zh)
        if z == q:
            key = (0, 0.0)
        elif z.startswith(q):
            key = (1, 0.0)
        elif q in z:
            key = (2, 0.0)
        else:
            overlap = len(qb & _bigrams(z)) / len(qb | _bigrams(z))
            if overlap < 0.2:
                continue
            key = (3, -overlap)
        under_pheno = ROOT_PHENO in idx.ancestors(tid)
        scored.append((key[0], key[1], 0 if under_pheno else 1, len(z), tid, zh))
    scored.sort()
    return [{"id": s[4], "label": idx.names[s[4]], "label_zh": s[5], "matched": s[5]} for s in scored[:limit]]


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

    excl_closure = {e: idx.ancestors(e) for e in excl}
    results = []
    for di in keep:
        did = idx.disease_list[di]
        score = sum(best[di]) / len(query)
        penalty = 0.0
        hits = []
        if excl:
            for term, w in idx.disease_terms[did].items():
                anc = idx.ancestors(idx.alt.get(term, term))
                for e in excl:
                    if e in anc:  # disease has the excluded term or something below it
                        p = idx.ic.get(e, 0.0) * w / len(query)
                        penalty = max(penalty, p) if hits else p
                        hits.append(e)
        results.append((score - penalty, score, penalty, di, sorted(set(hits))))
    results.sort(key=lambda r: r[0], reverse=True)

    out = []
    for adj, raw, penalty, di, hits in results[:top]:
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
            "genes": idx.disease_genes.get(did, []), "matches": matched,
        })
    # a gene scores as its best disease; ties are broken by its support across the next-best diseases
    gene_scores: Dict[str, Tuple[float, str]] = {}
    gene_support: Dict[str, float] = defaultdict(float)
    for rank_i, (adj, raw, penalty, di, hits) in enumerate(results):
        if adj <= 0 or rank_i >= 200:
            break
        did = idx.disease_list[di]
        for g in idx.disease_genes.get(did, []):
            gene_support[g] += adj
            if g not in gene_scores or gene_scores[g][0] < adj:
                gene_scores[g] = (adj, did)
    genes = sorted(gene_scores.items(), key=lambda kv: (kv[1][0], gene_support[kv[0]]), reverse=True)[:top]
    return {
        "method": "Resnik best-match average (query→disease), IC from HPO disease annotations; excluded terms penalised",
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
