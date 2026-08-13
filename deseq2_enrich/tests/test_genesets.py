from unittest.mock import patch

from deseq2_enrich import config, genesets


def test_go_gsea_libraries_are_grouped_and_labeled():
    assert "GO_Biological_Process_2026" in config.GSEA_LIBRARIES
    assert "GO_Molecular_Function_2026" in config.GSEA_LIBRARIES
    assert "GO_Cellular_Component_2026" in config.GSEA_LIBRARIES
    # GO slim leads the group: it is the only GO option that fits the
    # deployed app's memory budget, so it should be the easiest to reach.
    assert config.GSEA_LIBRARY_GROUPS["Gene Ontology"] == [
        "GO_Slim_Chicken",
        "GO_Biological_Process_2026",
        "GO_Molecular_Function_2026",
        "GO_Cellular_Component_2026",
    ]


def test_combine_libraries_prefixes_go_terms():
    with patch(
        "deseq2_enrich.genesets.fetch_library",
        return_value={"positive regulation": ["A", "B", "C"]},
    ):
        combined = genesets.combine_libraries(["GO_Biological_Process_2026"])

    assert combined == {"GO:BP | positive regulation": ["A", "B", "C"]}


# --------------------------------------------------------------------------
# GO sub-collections (slim / namespace split)
# --------------------------------------------------------------------------
_FAKE_GMT = "\n".join([
    "GO:0000001\tbp term\tGENE1\tGENE2",
    "GO:0000002\tmf term\tGENE3\tGENE4",
    "GO:0000003\tslim term\tGENE5\tGENE6",
    "REAC:R-GGA-1\treac term\tGENE7\tGENE8",
]).encode()


def _patch_sources(monkeypatch=None):
    return (
        patch("deseq2_enrich.genesets._read_obo"),
        patch("deseq2_enrich.genesets.urllib.request.urlopen"),
    )


def test_go_namespace_split_selects_only_requested_branch():
    """GO:BP must not silently include MF/CC (the GMT has no namespace)."""
    genesets.fetch_chicken_gmt.cache_clear()
    genesets.fetch_go_namespaces.cache_clear()
    genesets.fetch_go_slim_ids.cache_clear()
    ns = {"GO:0000001": "biological_process", "GO:0000002": "molecular_function",
          "GO:0000003": "biological_process"}

    class Resp:
        def read(self): return _FAKE_GMT
        def __enter__(self): return self
        def __exit__(self, *a): return False

    with patch("deseq2_enrich.genesets.fetch_go_namespaces", return_value=ns), \
         patch("deseq2_enrich.genesets.urllib.request.urlopen", return_value=Resp()):
        sets = genesets.fetch_chicken_gmt(sources=("GO",), go_subsets=("GO:BP",))

    assert set(sets) == {"GO | bp term", "GO | slim term"}
    assert "GO | mf term" not in sets
    genesets.fetch_chicken_gmt.cache_clear()


def test_go_slim_selects_only_slim_terms():
    genesets.fetch_chicken_gmt.cache_clear()

    class Resp:
        def read(self): return _FAKE_GMT
        def __enter__(self): return self
        def __exit__(self, *a): return False

    with patch("deseq2_enrich.genesets.fetch_go_slim_ids",
               return_value=frozenset({"GO:0000003"})), \
         patch("deseq2_enrich.genesets.urllib.request.urlopen", return_value=Resp()):
        sets = genesets.fetch_chicken_gmt(sources=("GO",), go_subsets=("GO:SLIM",))

    assert set(sets) == {"GO | slim term"}
    genesets.fetch_chicken_gmt.cache_clear()


def test_obo_parsers_reject_a_redirect_stub():
    """A 167-byte redirect stub must raise, not silently yield zero terms."""
    genesets.fetch_go_slim_ids.cache_clear()
    with patch("deseq2_enrich.genesets._read_obo", return_value=["<html>moved</html>"]):
        try:
            genesets.fetch_go_slim_ids()
        except RuntimeError as exc:
            assert "zero terms" in str(exc)
        else:
            raise AssertionError("expected RuntimeError on an empty parse")
    genesets.fetch_go_slim_ids.cache_clear()
