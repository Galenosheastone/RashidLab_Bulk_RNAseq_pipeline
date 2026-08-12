from deseq2_enrich.ora import _flatten_intersection


def test_flat_gene_list():
    assert _flatten_intersection(["IL6", "IL1B", "IFNG"]) == ["IL6", "IL1B", "IFNG"]


def test_none_returns_empty():
    assert _flatten_intersection(None) == []


def test_empty_returns_empty():
    assert _flatten_intersection([]) == []


def test_scalar_returns_empty():
    assert _flatten_intersection("IL6") == []
