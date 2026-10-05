"""Upstream sources: one module per database, each returning data with provenance.

This package root holds the three guards every source module shares.

`attempt` is `zebra.core.attempt` widened with `AttributeError` (F11). The
parsers call `.get` on whatever the service sent, so a list or a string where a
dict was expected raises `AttributeError`; `zebra.core.attempt` does not catch
it, the exception escapes `cli.main`, stdout stays empty and every *other*
source that did answer is lost. Source modules and the commands that drive them
import `attempt` from here.

`public_url` removes credential parameters from a URL. `zebra.http.request`
bakes `params` into `Response.url`, and `source_record` copies that URL into the
evidence ledger, the JSON the model reads and any report that cites sources
(E8). Nothing that carries a key may reach `source_record`, so source modules
call `record()` (defined here) instead of `zebra.http.source_record`.

`validated_json` rejects the bodies a 200 can still carry — empty, HTML, or
JSON that is really an error document — before a caller parses them (F10).
`zebra.http.request` writes its cache before any caller has seen the body, so
such an answer is otherwise served from disk for 7-30 days; raising here turns
it into a warning now, and `refetch=` re-requests once when the bad body came
from the cache.
"""

from __future__ import annotations

import urllib.parse
from typing import Any, Callable, Dict, List, Optional, TypeVar

from zebra.http import Response, SourceError
from zebra.http import source_record as _source_record

T = TypeVar("T")

# query parameters that must never appear in a recorded URL
CREDENTIAL_PARAMS = ("api_key", "apikey", "api-key", "key", "token", "access_token", "email", "password", "secret")


def attempt(label: str, fn: Callable[[], T], warnings: List[str]) -> Optional[T]:
    """Run one source call; on failure record a warning and return None.

    Same contract as `zebra.core.attempt`, plus `AttributeError`: an upstream
    shape change must cost one card section, not the whole command.
    """
    try:
        return fn()
    except SourceError as err:
        warnings.append(f"{label} unavailable: {err.message} (HTTP {err.status or '-'})")
    except (KeyError, ValueError, TypeError, IndexError, AttributeError) as err:
        warnings.append(f"{label}: unexpected response shape ({type(err).__name__}: {err})")
    return None


def public_url(url: Optional[str]) -> Optional[str]:
    """`url` with every credential parameter removed; the rest is left byte-for-byte."""
    if not url:
        return url
    split = urllib.parse.urlsplit(url)
    if not split.query:
        return url
    kept = [(k, v) for k, v in urllib.parse.parse_qsl(split.query, keep_blank_values=True)
            if k.lower() not in CREDENTIAL_PARAMS]
    if len(kept) == len(urllib.parse.parse_qsl(split.query, keep_blank_values=True)):
        return url
    return urllib.parse.urlunsplit((split.scheme, split.netloc, split.path,
                                    urllib.parse.urlencode(kept), split.fragment))


def record(db: str, rec_id: Optional[str], resp: Optional[Response] = None, *, url: Optional[str] = None,
           note: Optional[str] = None) -> Dict[str, Any]:
    """`zebra.http.source_record` with credentials stripped from the URL it stores."""
    out = _source_record(db, rec_id, resp, url=public_url(url), note=note)
    out["url"] = public_url(out.get("url"))
    return out


def validated_json(resp: Response, source: str, *, require: Any = None,
                   refetch: Optional[Callable[[], Response]] = None) -> Any:
    """Parse `resp` as JSON, rejecting an error body a 200 can still carry.

    Rejected: an empty body, a body that starts with `<`, a body that is not
    JSON, a JSON object with a top-level `error`/`errors`, and — when `require`
    names a key — an object that lacks it. When the bad body came from the cache
    and `refetch` is given, the request is made once more without the cache
    before giving up, so a transient upstream error is not replayed for days.
    """
    try:
        return _validate(resp, source, require)
    except SourceError:
        if not (resp.cached and refetch):
            raise
    fresh = refetch()
    return _validate(fresh, source, require)


def _validate(resp: Response, source: str, require: Any) -> Any:
    text = (resp.text or "").strip()
    if not text:
        raise SourceError(source, public_url(resp.url) or "", resp.status, "empty body where JSON was expected")
    if text[:1] == "<":
        raise SourceError(source, public_url(resp.url) or "", resp.status, "HTML page where JSON was expected")
    try:
        data = resp.json()
    except ValueError as err:
        raise SourceError(source, public_url(resp.url) or "", resp.status, f"body is not JSON ({err})") from None
    if isinstance(data, dict):
        err = data.get("error") or data.get("errors")
        if err:
            raise SourceError(source, public_url(resp.url) or "", resp.status, f"error body: {str(err)[:200]}")
        if require is not None and require not in data:
            raise SourceError(source, public_url(resp.url) or "", resp.status,
                              f"response has no {require!r} (keys: {', '.join(sorted(data)[:8]) or 'none'})")
    return data


def validated_text(resp: Response, source: str, *, must_contain: Optional[str] = None,
                   refetch: Optional[Callable[[], Response]] = None) -> str:
    """`resp.text` for a TSV/CSV/plain source, rejecting empty and HTML bodies (F10).

    `must_contain` is the header marker the parser needs; a maintenance page or a
    truncated download is rejected here instead of being cached for weeks.
    """
    try:
        return _validate_text(resp, source, must_contain)
    except SourceError:
        if not (resp.cached and refetch):
            raise
    return _validate_text(refetch(), source, must_contain)


def _validate_text(resp: Response, source: str, must_contain: Optional[str]) -> str:
    text = resp.text or ""
    if not text.strip():
        raise SourceError(source, public_url(resp.url) or "", resp.status, "empty body")
    if text.lstrip()[:1] == "<":
        raise SourceError(source, public_url(resp.url) or "", resp.status, "HTML page where tabular text was expected")
    if must_contain is not None and must_contain not in text:
        raise SourceError(source, public_url(resp.url) or "", resp.status,
                          f"body does not contain {must_contain!r}: not the expected table")
    return text
