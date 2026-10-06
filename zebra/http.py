"""HTTP client: urllib with retries, per-host pacing, an on-disk cache and provenance.

Every request made through here can be turned into a `Source` record
(database, record id, URL, retrieval time, whether it came from cache), which
is what the evidence ledger stores.
"""

from __future__ import annotations

import gzip
import hashlib
import http.client
import json
import math
import os
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import zlib
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Mapping, Optional

from zebra import __version__

USER_AGENT = f"zebra-mod/{__version__} (+https://github.com/zwbao/zebra-mod)"

# Minimum seconds between two requests to the same host. NCBI allows 3 rps
# without a key and 10 with one; the others are courtesy limits.
_HOST_INTERVAL = {
    "eutils.ncbi.nlm.nih.gov": 0.34,
    "www.ncbi.nlm.nih.gov": 0.34,
    "rest.ensembl.org": 0.07,
    "grch37.rest.ensembl.org": 0.07,
    "gnomad.broadinstitute.org": 6.0,  # gnomAD asks for ~10 requests/minute
    "panelapp.genomicsengland.co.uk": 0.2,
    "panelapp-aus.org": 0.2,
    "rest.genenames.org": 0.2,
    "rest.uniprot.org": 0.1,
    "alphafold.ebi.ac.uk": 0.1,
    "search.clinicalgenome.org": 0.5,
    "ftp.clinicalgenome.org": 0.5,
    "spliceai-38-xwkwwwxdwq-uc.a.run.app": 2.0,
    "spliceai-37-xwkwwwxdwq-uc.a.run.app": 2.0,
    "pangolin-38-xwkwwwxdwq-uc.a.run.app": 2.0,
    "pangolin-37-xwkwwwxdwq-uc.a.run.app": 2.0,
    "pubcasefinder.dbcls.jp": 6.0,  # 10 requests/minute published limit
    "clinicaltrials.gov": 0.2,
    "www.ebi.ac.uk": 0.1,
    "api.platform.opentargets.org": 0.1,
    "www.ema.europa.eu": 1.0,
    "api.orphadata.com": 0.2,
    "ftp.ncbi.nlm.nih.gov": 0.34,
    "api-v3.monarchinitiative.org": 0.1,
}
_last_call: Dict[str, float] = {}
_pace_lock = threading.Lock()

# Wall-clock budget for the whole process. The mod's tools run the CLI under a
# timeout; they pass that timeout minus a margin as ZEBRA_DEADLINE_MS, so a
# hung source becomes a named SourceError (which `core.attempt` turns into a
# warning) instead of the host killing the process and losing every answer.
_PROCESS_START = time.monotonic()  # import time, which for the CLI is within milliseconds of process start
_MIN_ATTEMPT_TIMEOUT = 0.2  # the floor an attempt's socket timeout is clamped to
READ_CHUNK = 1 << 16  # the body is read in chunks so the deadline is re-checked while it arrives


def deadline_seconds() -> Optional[float]:
    """Seconds left of ZEBRA_DEADLINE_MS, or None when no usable deadline is set."""
    raw = os.environ.get("ZEBRA_DEADLINE_MS")
    if not raw:
        return None
    try:
        budget = float(raw) / 1000.0
    except ValueError:
        return None
    if not math.isfinite(budget) or budget <= 0:  # nan/inf would disable every comparison below
        return None
    return budget - (time.monotonic() - _PROCESS_START)


def _read_body(resp: Any, source: str, url: str) -> bytes:
    """Read the response body, giving up if the deadline passes while it trickles in."""
    left = deadline_seconds()
    if left is None:
        return resp.read()
    parts = []
    while True:
        if deadline_seconds() <= 0:  # type: ignore[operator]
            raise SourceError(source, url, None,
                              f"deadline reached after {sum(len(p) for p in parts)} byte(s) of the response body: "
                              "the source answered but did not finish in the time available")
        # read1, not read: read() blocks until it has the whole chunk (or the
        # whole body, when Content-Length is set), so the deadline check
        # between chunks would never run against a slow trickle
        chunk = resp.read1(READ_CHUNK)
        if not chunk:
            return b"".join(parts)
        parts.append(chunk)


class SourceError(Exception):
    """A request to an upstream source failed after retries."""

    def __init__(self, source: str, url: str, status: Optional[int], message: str):
        super().__init__(f"{source}: {message} ({status or 'no status'}) <{url}>")
        self.source = source
        self.url = url
        self.status = status
        self.message = message


@dataclass
class Response:
    url: str  # credentials stripped: zebra.sources.public_url is applied where it is set
    status: int
    text: str
    retrieved_at: str
    cached: bool
    # the cache entry this answer was read from or written to, so a caller that finds
    # the body unusable can drop it (`evict`) instead of having it replayed for days
    cache_key: Optional[str] = None

    def json(self) -> Any:
        return json.loads(self.text)


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def cache_dir() -> Path:
    root = os.environ.get("ZEBRA_CACHE_DIR") or os.path.join(
        os.path.expanduser("~"), ".cache", "zebra-mod"
    )
    return Path(root)


def _cache_path(key: str) -> Path:
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
    return cache_dir() / "http" / digest[:2] / f"{digest}.json"


def _pace(host: str, max_wait: Optional[float] = None) -> Optional[float]:
    """Wait out the per-host courtesy interval. Returns the seconds waited.

    `max_wait` is the deadline's remaining budget: when the wait would exceed
    it, nothing is reserved and None is returned, so the caller can fail with a
    named error instead of sleeping past its budget.
    """
    interval = _HOST_INTERVAL.get(host, 0.0)
    if host.endswith("ncbi.nlm.nih.gov") and os.environ.get("NCBI_API_KEY"):
        interval = 0.11
    if interval <= 0:
        return 0.0
    with _pace_lock:
        now = time.monotonic()
        slot = max(now, _last_call.get(host, 0.0) + interval)
        wait = slot - now
        if max_wait is not None and wait > max_wait:
            return None  # not reserved: another caller may still use this slot
        _last_call[host] = slot
    if wait > 0:
        time.sleep(wait)
    return wait


def _wait_for_retry(delay: float) -> bool:
    """Sleep before the next attempt. False means the deadline leaves no room for one."""
    left = deadline_seconds()
    if left is None:
        time.sleep(delay)
        return True
    if left - delay <= _MIN_ATTEMPT_TIMEOUT:
        return False
    time.sleep(delay)
    return True


def _opener() -> urllib.request.OpenerDirector:
    # urllib honours http_proxy/https_proxy; a socks all_proxy is ignored by
    # ProxyHandler, which is what we want (urllib cannot speak socks).
    proxies = {
        k: v
        for k, v in urllib.request.getproxies().items()
        if k in ("http", "https", "no") and not str(v).startswith("socks")
    }
    return urllib.request.build_opener(urllib.request.ProxyHandler(proxies))


def request(
    url: str,
    *,
    source: str,
    params: Optional[Mapping[str, Any]] = None,
    method: str = "GET",
    body: Optional[Any] = None,
    headers: Optional[Mapping[str, str]] = None,
    timeout: float = 30.0,
    retries: int = 3,
    cache_ttl: float = 3 * 86400,
    accept: str = "application/json",
    ok_statuses: tuple = (200,),
    validate: Optional[Callable[[str], Optional[str]]] = None,
    refresh: bool = False,
    not_found_ttl: Optional[float] = None,
) -> Response:
    """Fetch `url` and return a Response, from cache when fresh.

    `body` is JSON-encoded unless it is already bytes or str. POST bodies are
    part of the cache key. `cache_ttl=0` bypasses the cache; so does
    ZEBRA_NO_CACHE=1.

    `validate(text)` returns why a body is unusable, or None when it is good: a
    body it rejects is never written to the cache and never returned, so a
    service answering an error, an empty document or an HTML page with HTTP 200
    cannot be served back for days (F10). `refresh=True` skips the cache read
    and still writes, which is how a caller recovers from a bad entry stored
    before this check existed. `not_found_ttl` caches an accepted
    "not found" status for its own, shorter time: whether a record exists
    changes far sooner than the record's content does, and the TTL has to be
    chosen before the status is known.
    """
    if params:
        clean = {k: v for k, v in params.items() if v is not None}
        query = urllib.parse.urlencode(clean, doseq=True)
        url = f"{url}{'&' if '?' in url else '?'}{query}"
    data: Optional[bytes] = None
    hdrs = {"User-Agent": USER_AGENT, "Accept": accept, "Accept-Encoding": "gzip"}
    if body is not None:
        if isinstance(body, bytes):
            data = body
        elif isinstance(body, str):
            data = body.encode("utf-8")
        else:
            data = json.dumps(body).encode("utf-8")
            hdrs["Content-Type"] = "application/json"
    if headers:
        hdrs.update(headers)

    use_cache = cache_ttl > 0 and os.environ.get("ZEBRA_NO_CACHE") != "1"
    key = f"{method} {url} {data.decode('utf-8', 'replace') if data else ''}"
    path = _cache_path(key)
    # Whether a record exists changes far sooner than its content: an accepted
    # "not found" lives a day unless the caller says otherwise (C-P1-4).
    if not_found_ttl is None:
        not_found_ttl = min(cache_ttl, NOT_FOUND_TTL)
    if use_cache and not refresh and path.exists():
        try:
            entry = json.loads(path.read_text("utf-8"))
            ttl = cache_ttl
            if entry["status"] != 200:
                ttl = min(ttl, not_found_ttl)
            if time.time() - entry["stored"] < ttl:
                bad = _body_problem(entry["text"], accept, validate) if entry["status"] == 200 else None
                if bad is None:
                    return Response(url, entry["status"], entry["text"], entry["retrieved_at"], True, key)
                # a bad body stored earlier: fetch again rather than serve it
        except (ValueError, KeyError, OSError):
            pass

    host = urllib.parse.urlsplit(url).hostname or ""
    opener = _opener()
    last_error: Optional[BaseException] = None
    status: Optional[int] = None
    for attempt in range(retries + 1):
        left = deadline_seconds()
        if left is not None and left <= 0:
            raise SourceError(
                source, url, status,
                f"deadline reached: gave up after {attempt} attempt(s) without an answer"
                + (f"; last error: {_short(str(last_error), 120)}" if last_error else ""),
            ) from None
        attempt_timeout = timeout if left is None else max(_MIN_ATTEMPT_TIMEOUT, min(timeout, left))
        waited = _pace(host, None if left is None else left)
        if waited is None:
            raise SourceError(source, url, status,
                              f"deadline reached while waiting out the {_HOST_INTERVAL.get(host, 0)} s courtesy "
                              f"interval for {host}: no request was sent") from None
        req = urllib.request.Request(url, data=data, headers=hdrs, method=method)
        try:
            with opener.open(req, timeout=attempt_timeout) as resp:
                raw = _gunzip(_read_body(resp, source, url))
                status = resp.status
                text = raw.decode("utf-8", "replace")
        except urllib.error.HTTPError as err:
            status = err.code
            try:
                raw = _gunzip(err.read())
                text = raw.decode("utf-8", "replace")
            except Exception:  # noqa: BLE001 - best effort error body
                text = ""
            if status in ok_statuses:
                pass
            elif status in (429, 500, 502, 503, 504) and attempt < retries:
                # 4xx other than 429 is the server's final answer: never retried,
                # whatever the body looks like
                retry_after = err.headers.get("Retry-After") if err.headers else None
                delay = float(retry_after) if retry_after and retry_after.isdigit() else 1.5 * (2**attempt)
                last_error = err
                if _wait_for_retry(min(delay, 20)):
                    continue
                raise SourceError(source, url, status,
                                  f"HTTP {status} and no time left to retry before the deadline: {_short(text) or str(err)}") from None
            else:
                raise SourceError(source, url, status, _short(text) or str(err)) from None
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError, http.client.HTTPException,
                EOFError, zlib.error) as err:
            last_error = err
            if attempt < retries and _wait_for_retry(1.5 * (2**attempt)):
                continue
            raise SourceError(source, url, None, f"network error: {err}") from None

        if status not in ok_statuses:
            raise SourceError(source, url, status, _short(text))
        bad = _body_problem(text, accept, validate) if status == 200 else None
        if bad is not None:
            # a body the caller cannot use: never cached, and only retried while
            # the server might still answer differently
            last_error = SourceError(source, url, status, bad)
            if attempt < retries and _wait_for_retry(1.5 * (2**attempt)):
                continue
            raise last_error
        retrieved = now_iso()
        if use_cache:
            try:
                path.parent.mkdir(parents=True, exist_ok=True)
                tmp = path.with_suffix(f".{os.getpid()}.{threading.get_ident()}.tmp")
                tmp.write_text(
                    json.dumps({"stored": time.time(), "status": status, "text": text, "retrieved_at": retrieved}),
                    "utf-8",
                )
                os.replace(tmp, path)
            except OSError:
                pass
        return Response(url, status, text, retrieved, False, key if use_cache else None)
    raise SourceError(source, url, status, f"gave up: {last_error}")


NOT_FOUND_TTL = 86400.0


def _body_problem(text: str, accept: str, validate: Optional[Callable[[str], Optional[str]]]) -> Optional[str]:
    """Why a 200 body cannot be used, or None. Checked before a body is cached and when it is read back."""
    if accept == "application/json":
        head = text.lstrip()[:1]
        if not head:
            return "empty body where JSON was expected"
        if head == "<":
            return "HTML page where JSON was expected"
    return validate(text) if validate is not None else None


def evict(resp: "Response") -> None:
    """Drop the cache entry behind `resp`: its body turned out to be unusable (an error document
    sent with HTTP 200), and it must not be served again for the rest of its time to live."""
    key = resp.cache_key or f"GET {resp.url} "
    try:
        _cache_path(key).unlink()
    except OSError:
        pass


def get_json(url: str, *, source: str, **kw: Any) -> Response:
    return request(url, source=source, **kw)


def post_json(url: str, payload: Any, *, source: str, **kw: Any) -> Response:
    return request(url, source=source, method="POST", body=payload, **kw)


def _gunzip(raw: bytes) -> bytes:
    """Undo gzip however many times it was applied (some servers double-encode).

    A body that claims to be gzip and is not decodable is returned as it came:
    a corrupt upstream body must become a SourceError upstream, never a bare
    zlib.error that no caller catches.
    """
    for _ in range(3):
        if raw[:2] != b"\x1f\x8b":
            break
        try:
            raw = gzip.decompress(raw)
        except (OSError, EOFError, zlib.error):
            break
    return raw


def _short(text: str, limit: int = 300) -> str:
    text = " ".join((text or "").split())
    return text[:limit]


def source_record(db: str, record: Optional[str], resp: Optional[Response] = None, *, url: Optional[str] = None, note: Optional[str] = None) -> Dict[str, Any]:
    """A provenance record for the evidence ledger."""
    rec: Dict[str, Any] = {"db": db, "record": record, "url": url or (resp.url if resp else None)}
    if resp is not None:
        rec["retrieved_at"] = resp.retrieved_at
        rec["cached"] = resp.cached
    else:
        rec["retrieved_at"] = now_iso()
    if note:
        rec["note"] = note
    return rec
