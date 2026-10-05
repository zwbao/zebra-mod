"""HTTP client: urllib with retries, per-host pacing, an on-disk cache and provenance.

Every request made through here can be turned into a `Source` record
(database, record id, URL, retrieval time, whether it came from cache), which
is what the evidence ledger stores.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from zebra import __version__

USER_AGENT = f"zebra-mod/{__version__} (+https://github.com/zwbao/zebra-mod)"

# Minimum seconds between two requests to the same host. NCBI allows 3 rps
# without a key and 10 with one; the others are courtesy limits.
_HOST_INTERVAL = {
    "eutils.ncbi.nlm.nih.gov": 0.34,
    "www.ncbi.nlm.nih.gov": 0.34,
    "rest.ensembl.org": 0.07,
    "grch37.rest.ensembl.org": 0.07,
    "gnomad.broadinstitute.org": 0.5,
    "pubcasefinder.dbcls.jp": 0.5,
    "clinicaltrials.gov": 0.2,
}
_last_call: Dict[str, float] = {}


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
    url: str
    status: int
    text: str
    retrieved_at: str
    cached: bool

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


def _pace(host: str) -> None:
    interval = _HOST_INTERVAL.get(host, 0.0)
    if host.endswith("ncbi.nlm.nih.gov") and os.environ.get("NCBI_API_KEY"):
        interval = 0.11
    if interval <= 0:
        return
    last = _last_call.get(host, 0.0)
    wait = last + interval - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_call[host] = time.monotonic()


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
) -> Response:
    """Fetch `url` and return a Response, from cache when fresh.

    `body` is JSON-encoded unless it is already bytes or str. POST bodies are
    part of the cache key. `cache_ttl=0` bypasses the cache; so does
    ZEBRA_NO_CACHE=1.
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
    if use_cache and path.exists():
        try:
            entry = json.loads(path.read_text("utf-8"))
            if time.time() - entry["stored"] < cache_ttl:
                return Response(url, entry["status"], entry["text"], entry["retrieved_at"], True)
        except (ValueError, KeyError, OSError):
            pass

    host = urllib.parse.urlsplit(url).hostname or ""
    opener = _opener()
    last_error: Optional[BaseException] = None
    status: Optional[int] = None
    for attempt in range(retries + 1):
        _pace(host)
        req = urllib.request.Request(url, data=data, headers=hdrs, method=method)
        try:
            with opener.open(req, timeout=timeout) as resp:
                raw = resp.read()
                if resp.headers.get("Content-Encoding") == "gzip":
                    raw = gzip.decompress(raw)
                status = resp.status
                text = raw.decode("utf-8", "replace")
        except urllib.error.HTTPError as err:
            status = err.code
            try:
                raw = err.read()
                if err.headers.get("Content-Encoding") == "gzip":
                    raw = gzip.decompress(raw)
                text = raw.decode("utf-8", "replace")
            except Exception:  # noqa: BLE001 - best effort error body
                text = ""
            if status in ok_statuses:
                pass
            elif status in (429, 500, 502, 503, 504) and attempt < retries:
                retry_after = err.headers.get("Retry-After") if err.headers else None
                delay = float(retry_after) if retry_after and retry_after.isdigit() else 1.5 * (2**attempt)
                time.sleep(min(delay, 20))
                last_error = err
                continue
            else:
                raise SourceError(source, url, status, _short(text) or str(err)) from None
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as err:
            last_error = err
            if attempt < retries:
                time.sleep(1.5 * (2**attempt))
                continue
            raise SourceError(source, url, None, f"network error: {err}") from None

        if status not in ok_statuses:
            raise SourceError(source, url, status, _short(text))
        retrieved = now_iso()
        if use_cache:
            try:
                path.parent.mkdir(parents=True, exist_ok=True)
                tmp = path.with_suffix(".tmp")
                tmp.write_text(
                    json.dumps({"stored": time.time(), "status": status, "text": text, "retrieved_at": retrieved}),
                    "utf-8",
                )
                os.replace(tmp, path)
            except OSError:
                pass
        return Response(url, status, text, retrieved, False)
    raise SourceError(source, url, status, f"gave up: {last_error}")


def get_json(url: str, *, source: str, **kw: Any) -> Response:
    return request(url, source=source, **kw)


def post_json(url: str, payload: Any, *, source: str, **kw: Any) -> Response:
    return request(url, source=source, method="POST", body=payload, **kw)


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
