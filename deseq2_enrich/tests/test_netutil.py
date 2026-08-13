"""Phase 3: retry behaviour and run provenance."""
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from deseq2_enrich import genesets, netutil, ortho


def test_retry_call_succeeds_on_second_attempt():
    """A single transient blip must not lose the stage."""
    fn = MagicMock(side_effect=[ConnectionError("boom"), "ok"])
    with patch("deseq2_enrich.netutil.time.sleep"):  # no real backoff in tests
        assert netutil.retry_call(fn) == "ok"
    assert fn.call_count == 2


def test_retry_call_reraises_original_error_after_exhausting():
    fn = MagicMock(side_effect=ValueError("persistent"))
    with patch("deseq2_enrich.netutil.time.sleep"):
        with pytest.raises(ValueError, match="persistent"):
            netutil.retry_call(fn)
    assert fn.call_count == netutil.RETRY_ATTEMPTS


def test_retry_call_does_not_retry_on_success():
    fn = MagicMock(return_value=42)
    assert netutil.retry_call(fn) == 42
    assert fn.call_count == 1


def test_fetch_library_retries(monkeypatch):
    genesets.fetch_library.cache_clear()
    calls = {"n": 0}

    def flaky(name, organism):
        calls["n"] += 1
        if calls["n"] == 1:
            raise ConnectionError("enrichr blip")
        return {"TERM": ["A", "B"]}

    with patch("deseq2_enrich.genesets.gp.get_library", side_effect=flaky), \
         patch("deseq2_enrich.netutil.time.sleep"):
        assert genesets.fetch_library("X", "human") == {"TERM": ["A", "B"]}
    assert calls["n"] == 2
    genesets.fetch_library.cache_clear()


def test_orth_cached_retries():
    ortho.clear_cache()
    frame = pd.DataFrame({"incoming": ["G1"], "ortholog_name": ["HS1"]})
    client = MagicMock()
    client.orth.side_effect = [ConnectionError("gprofiler 503"), frame]

    with patch("deseq2_enrich.ortho._client", return_value=client), \
         patch("deseq2_enrich.netutil.time.sleep"):
        out = ortho._orth_cached(("G1",), "ggallus", "hsapiens")
    assert client.orth.call_count == 2
    assert out.equals(frame)
    ortho.clear_cache()


def test_package_versions_and_timestamp():
    versions = netutil.package_versions()
    assert set(versions) == {"gseapy", "gprofiler"}
    # Real versions, not placeholders, in a working environment.
    assert versions["gseapy"] != "unknown"
    stamp = netutil.utc_now_iso()
    assert stamp.endswith("+00:00") and "T" in stamp
