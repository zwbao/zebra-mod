"""zebra doctor: Python, data files, case, API reachability, keys and helper tools — in under 30 s."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import threading
import time
import uuid
from typing import Any, Callable, Dict, List, Optional, Tuple

from zebra import __version__
from zebra.core import Outcome
from zebra.http import SourceError, cache_dir, request

PROBE_TIMEOUT = 8.0  # per attempt; one retry (TLS resets are the usual failure and clear on retry)
DEADLINE = 25.0  # whole doctor run
_PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# name, url, method, body, accept, validator(text) -> detail or raises, note
_GET = "GET"


def _json_has(*path: str) -> Callable[[str], str]:
    def check(text: str) -> str:
        data: Any = json.loads(text)
        for p in path:
            data = data[int(p)] if isinstance(data, list) else data[p]
        return ""
    return check


def _ensembl_ping(text: str) -> str:
    if json.loads(text).get("ping") != 1:
        raise ValueError("ping != 1")
    return ""


def _spliceai(text: str) -> str:
    data = json.loads(text)
    if data.get("error") or data.get("inputError"):
        raise ValueError(str(data.get("error") or "inputError"))
    if not data.get("scores"):
        raise ValueError("no scores")
    return "Cloud Run backend of spliceailookup.broadinstitute.org (undocumented URL)"


def _html(text: str) -> str:
    if "<html" not in text.lower():
        raise ValueError("no HTML page")
    return "host reachable (homepage; the ranking API returns whole lists and is not exercised here)"


NETWORK: List[Tuple[str, str, str, Optional[Any], str, Callable[[str], str]]] = [
    ("HPO JAX API", "https://ontology.jax.org/api/hp/terms/HP:0001250", _GET, None, "application/json", _json_has("id")),
    ("Monarch", "https://api-v3.monarchinitiative.org/v3/api/entity/MONDO:0100135", _GET, None, "application/json",
     _json_has("id")),
    ("PubCaseFinder", "https://pubcasefinder.dbcls.jp/", _GET, None, "text/html", _html),
    ("Orphadata", "https://api.orphadata.com/rd-cross-referencing/orphacodes/33069?lang=en", _GET, None,
     "application/json", _json_has("data")),
    ("OLS", "https://www.ebi.ac.uk/ols4/api/ontologies/hp", _GET, None, "application/json", _json_has("ontologyId")),
    ("NCBI E-utilities",
     "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=SCN1A%5Bsym%5D+AND+human%5Borgn%5D&retmode=json&retmax=1",
     _GET, None, "application/json", _json_has("esearchresult", "idlist")),
    ("gnomAD", "https://gnomad.broadinstitute.org/api?query=%7Bmeta%7Bclinvar_release_date%7D%7D", _GET, None,
     "application/json", _json_has("data", "meta")),
    ("Ensembl REST GRCh38", "https://rest.ensembl.org/info/ping?content-type=application/json", _GET, None,
     "application/json", _ensembl_ping),
    ("Ensembl REST GRCh37", "https://grch37.rest.ensembl.org/info/ping?content-type=application/json", _GET, None,
     "application/json", _ensembl_ping),
    ("SpliceAI lookup",
     "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=38&variant=chr8-140300616-T-G&distance=50&mask=0",
     _GET, None, "application/json", _spliceai),
    ("Europe PMC",
     "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=SCN1A&format=json&pageSize=1&resultType=idlist",
     _GET, None, "application/json", _json_has("hitCount")),
    ("PubTator", "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/entity/autocomplete/?query=SCN1A&limit=1", _GET,
     None, "application/json", _json_has("0", "_id")),
    ("LitVar", "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/?query=SCN1A", _GET, None,
     "application/json", _json_has("0", "_id")),
    ("ClinicalTrials.gov", "https://clinicaltrials.gov/api/v2/version", _GET, None, "application/json",
     _json_has("apiVersion")),
    ("Open Targets", "https://api.platform.opentargets.org/api/v4/graphql", "POST",
     {"query": "{meta{apiVersion{x y z}}}"}, "application/json", _json_has("data", "meta", "apiVersion")),
    ("PanelApp", "https://panelapp.genomicsengland.co.uk/api/v1/genes/SCN1A/?page=1", _GET, None, "application/json",
     _json_has("count")),
    ("ClinGen", "https://search.clinicalgenome.org/api/genes/look/SCN1A", _GET, None, "application/json, text/html",
     _json_has("0", "hgnc")),
    ("UniProt", "https://rest.uniprot.org/uniprotkb/P35498?fields=accession&format=json", _GET, None,
     "application/json", _json_has("primaryAccession")),
]

# UniProt double-gzips its body for this client (Content-Encoding: gzip over an already gzipped payload);
# zebra.http._gunzip unwraps up to 3 gzip layers; ask for an uncompressed
# answer here so this probe measures the service, not the decoder.
_HEADERS = {"UniProt": {"Accept-Encoding": "identity"}}

KEYS = [
    ("AlphaGenome key", ("ALPHAGENOME_API_KEY",), "AlphaGenome variant-effect predictions (zebra-s2f)"),
    ("Evo 2 key", ("NVCF_RUN_KEY", "EVO2_API_KEY"), "Evo 2 scoring via NVIDIA NIM (zebra-s2f)"),
    ("NCBI API key", ("NCBI_API_KEY",), "10 instead of 3 requests/s to NCBI E-utilities"),
    # no source in this tree reads OMIM_API_KEY; it is listed so a user who set it
    # is told it is unused rather than assuming it is in effect.
    ("OMIM API key", ("OMIM_API_KEY",), "NOT USED by any zebra source: OMIM entries come from Orphanet "
                                        "cross-references and HPO annotations"),
]


def _check(name: str, ok: bool, detail: str) -> Dict[str, Any]:
    return {"name": name, "ok": bool(ok), "detail": detail}


def _probe(name: str, url: str, method: str, body: Optional[Any], accept: str,
           validate: Callable[[str], str]) -> Dict[str, Any]:
    t0 = time.monotonic()
    try:
        resp = request(url, source=name, method=method, body=body, accept=accept, timeout=PROBE_TIMEOUT,
                       retries=1, cache_ttl=0, headers=_HEADERS.get(name))
    except SourceError as err:
        took = time.monotonic() - t0
        what = f"HTTP {err.status}" if err.status else "unreachable"
        return _check(f"network: {name}", False, f"{what} after {took:.1f} s: {err.message[:120]}")
    took = time.monotonic() - t0
    try:
        extra = validate(resp.text)
    except Exception as err:  # noqa: BLE001 - any shape problem is a failed check
        return _check(f"network: {name}", False,
                      f"HTTP {resp.status} in {took:.1f} s but unexpected body ({type(err).__name__}: {str(err)[:80]})")
    return _check(f"network: {name}", True, f"HTTP {resp.status} in {took:.1f} s" + (f"; {extra}" if extra else ""))


def _tool(name: str, exe: str, args: List[str], missing_hint: str) -> Dict[str, Any]:
    path = shutil.which(exe)
    if not path:
        return _check(name, False, f"not on PATH — {missing_hint}")
    try:
        ran = subprocess.run([path] + args, capture_output=True, text=True, timeout=10)
        out = (ran.stdout or ran.stderr or "").strip().splitlines()
        version = out[0].strip() if out else f"exit {ran.returncode}"
        return _check(name, ran.returncode == 0, f"{version} ({path})")
    except subprocess.TimeoutExpired:
        return _check(name, False, f"{path} did not answer {' '.join(args)} within 10 s")
    except OSError as err:
        return _check(name, False, f"{path}: {err}")


def _python() -> Dict[str, Any]:
    v = sys.version_info
    return _check("Python ≥ 3.9", v >= (3, 9), f"{v.major}.{v.minor}.{v.micro} ({sys.executable})")


def _cache() -> Dict[str, Any]:
    root = cache_dir()
    try:
        root.mkdir(parents=True, exist_ok=True)
        probe = root / f".doctor-{uuid.uuid4().hex}"
        probe.write_text("ok", "utf-8")
        probe.unlink()
        return _check("cache dir writable", True, str(root))
    except OSError as err:
        return _check("cache dir writable", False, f"{root}: {err} (set ZEBRA_CACHE_DIR to a writable folder)")


def _hpo() -> Dict[str, Any]:
    from zebra import hpo_local

    d = hpo_local.data_dir()
    missing = hpo_local.missing_files()
    if missing:
        return _check("HPO local files", False, f"missing in {d}: {', '.join(missing)} — run: zebra hpo fetch")
    version, hp_release = None, None
    try:
        with open(d / "phenotype.hpoa", encoding="utf-8") as fh:
            for _ in range(10):
                line = fh.readline()
                if not line.startswith("#"):
                    break
                if line.startswith("#version:"):
                    version = line.split(":", 1)[1].strip()
                elif line.startswith("#hpo-version:"):
                    hp_release = line.split(":", 1)[1].strip().rsplit("/releases/", 1)[-1].split("/")[0]
    except OSError as err:
        return _check("HPO local files", False, f"{d}: {err}")
    size = sum((d / f).stat().st_size for f in hpo_local.FILES) / 1e6
    detail = f"annotations {version or '?'}, ontology {hp_release or '?'} ({size:.0f} MB in {d})"
    return _check("HPO local files", True, detail)


def _case(case_dir: Optional[str]) -> Dict[str, Any]:
    if not case_dir:
        return _check("case", True, "no active case (set ZEBRA_CASE, pass --case, or /zebra case <dir>)")
    from zebra import case as case_mod

    try:
        data = case_mod.load(case_dir)
    except (case_mod.CaseError, OSError) as err:
        return _check("case", False, str(err))
    n_ph = sum(1 for p in data.get("phenotypes", []) if p.get("status", "present") == "present")
    return _check("case", True, f"{case_mod.case_file(case_dir).parent}: {data.get('title')!r}, "
                                f"{n_ph} phenotypes, {len(data.get('variants', []))} variants")


def _keys() -> List[Dict[str, Any]]:
    out = []
    for name, envs, what in KEYS:
        present = [e for e in envs if os.environ.get(e, "").strip()]
        if present:
            out.append(_check(name, True, f"{present[0]} is set"))
        else:
            out.append(_check(name, False, f"{' / '.join(envs)} not set — optional: {what}"))
    return out


def run_checks(case_dir: Optional[str] = None, deadline: float = DEADLINE,
               network: Optional[List[Tuple[str, str, str, Optional[Any], str, Callable[[str], str]]]] = None
               ) -> Dict[str, Any]:
    """All checks; slow ones run on daemon threads so a hung host can never hold the process past the deadline."""
    t0 = time.monotonic()
    network = NETWORK if network is None else network
    jobs: List[Tuple[str, Callable[[], Dict[str, Any]]]] = []
    for item in network:
        jobs.append((f"network: {item[0]}", (lambda it=item: _probe(*it))))
    jobs.append(("s2f CLI", lambda: _tool("s2f CLI", "s2f", ["--version"],
                                          'optional: uv tool install "git+https://github.com/zwbao/s2f-penguin"')))
    jobs.append(("tabix", lambda: _tool("tabix", "tabix", ["--version"], "optional: htslib, needed by GPN-MSA in s2f")))
    jobs.append(("HPO local files", _hpo))
    results: Dict[str, Dict[str, Any]] = {}
    lock = threading.Lock()

    def runner(label: str, fn: Callable[[], Dict[str, Any]]) -> None:
        try:
            res = fn()
        except Exception as err:  # noqa: BLE001 - the doctor never raises
            res = _check(label, False, f"check crashed: {type(err).__name__}: {err}")
        with lock:
            results[label] = res

    threads = []
    for label, fn in jobs:
        th = threading.Thread(target=runner, args=(label, fn), daemon=True)
        th.start()
        threads.append(th)

    local: List[Dict[str, Any]] = []
    for fn in (_python, lambda: _check("zebra", True, f"zebra {__version__} ({_PKG_DIR})"), _cache,
               lambda: _case(case_dir)):
        try:
            local.append(fn())
        except Exception as err:  # noqa: BLE001 - the doctor never raises
            local.append(_check("local check", False, f"crashed: {type(err).__name__}: {err}"))
    keys = _keys()

    for th in threads:
        th.join(max(0.0, deadline - (time.monotonic() - t0)))
    with lock:
        done = dict(results)

    def got(label: str) -> Dict[str, Any]:
        return done.get(label) or _check(label, False, f"timed out (no answer within the {deadline:.0f} s budget)")

    checks = local[:3] + [got("HPO local files")] + local[3:] + keys + [got("s2f CLI"), got("tabix")]
    checks += [got(lbl) for lbl, _ in jobs if lbl.startswith("network: ")]
    seconds = round(time.monotonic() - t0, 1)
    n_ok = sum(1 for c in checks if c["ok"])
    return {"checks": checks, "summary": {"ok": n_ok, "failed": len(checks) - n_ok, "seconds": seconds}}


def _doctor(args: argparse.Namespace) -> Outcome:
    res = run_checks(case_dir=getattr(args, "case", None))
    lines = [f"{'✓' if c['ok'] else '✗'} {c['name']} — {c['detail']}" for c in res["checks"]]
    s = res["summary"]
    lines.append(f"{s['ok']} ok, {s['failed']} not ok, {s['seconds']} s")
    return Outcome(res, text="\n".join(lines), query={})


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("doctor", help="check Python, data files, case, API reachability, keys and helper tools")
    p.set_defaults(func=_doctor)
