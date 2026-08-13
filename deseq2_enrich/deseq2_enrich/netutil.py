"""Retry and provenance helpers for the network-backed stages.

ORA, orthology and gene-set fetching all call keyless public services
(g:Profiler, Enrichr). A single transient blip should not lose a stage that
may already have cost the user a minute of compute, so those calls get one
cheap retry. Kept deliberately small: no new dependency, no backoff library.

Provenance lives here too because it answers the same question -- *what did
this run actually talk to, and when* -- which matters for citing a result.
Enrichr libraries update in place and ``GO_*_2026`` will drift, so the fetch
date is part of reproducing a figure.
"""
from __future__ import annotations

import time
from datetime import datetime, timezone
from functools import lru_cache
from typing import Callable, TypeVar

T = TypeVar("T")

RETRY_ATTEMPTS = 2
RETRY_BACKOFF_SECONDS = 1.5


def retry_call(
    fn: Callable[..., T],
    *args,
    attempts: int = RETRY_ATTEMPTS,
    backoff: float = RETRY_BACKOFF_SECONDS,
    **kwargs,
) -> T:
    """Call ``fn``, retrying once on failure before giving up.

    Catches ``Exception`` rather than a narrow tuple because the three
    services raise wildly different types for what is the same transient
    condition (urllib ``URLError``, requests' ``ConnectionError``, and
    g:Profiler surfacing an HTTP error as ``ValueError``). ``BaseException``
    is deliberately not caught, so Ctrl-C still interrupts immediately.

    The final failure is re-raised unchanged, so callers keep their existing
    error handling and messages.
    """
    last: BaseException | None = None
    for attempt in range(max(1, attempts)):
        try:
            return fn(*args, **kwargs)
        except Exception as exc:
            last = exc
            if attempt < attempts - 1:
                time.sleep(backoff * (attempt + 1))
    raise last  # type: ignore[misc]


@lru_cache(maxsize=1)
def package_versions() -> dict[str, str]:
    """Versions of the analysis backends, for the run manifest."""
    import importlib.metadata as md

    out: dict[str, str] = {}
    for dist, label in (("gseapy", "gseapy"), ("gprofiler-official", "gprofiler")):
        try:
            out[label] = md.version(dist)
        except Exception:
            out[label] = "unknown"
    return out


def utc_now_iso() -> str:
    """Timestamp for the run manifest (UTC, second resolution)."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()
