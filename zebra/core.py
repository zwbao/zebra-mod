"""The shape every command returns, and helpers for gathering several sources.

A command returns an `Outcome`: the result, the sources it stands on and the
warnings it raised. A source that fails is never dropped silently: it becomes
a warning naming the source and the reason, and the rest of the answer stands.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, TypeVar

from zebra.http import SourceError

T = TypeVar("T")


@dataclass
class Outcome:
    result: Any
    sources: List[Dict[str, Any]] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    text: Optional[str] = None  # human rendering; JSON is the default for tools
    query: Dict[str, Any] = field(default_factory=dict)

    def add(self, other: "Outcome") -> Any:
        self.sources.extend(other.sources)
        self.warnings.extend(other.warnings)
        return other.result


def attempt(label: str, fn: Callable[[], T], warnings: List[str]) -> Optional[T]:
    """Run one source call; on failure record a warning and return None."""
    try:
        return fn()
    except SourceError as err:
        warnings.append(f"{label} unavailable: {err.message} (HTTP {err.status or '-'})")
    except (KeyError, ValueError, TypeError, IndexError) as err:
        warnings.append(f"{label}: unexpected response shape ({type(err).__name__}: {err})")
    return None


class UsageError(Exception):
    """Bad input from the caller (exit code 2, message shown as is)."""
