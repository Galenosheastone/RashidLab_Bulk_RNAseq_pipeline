import io

from deseq2_enrich import io as pio


def test_column_alias_detection():
    csv = (
        "gene\tsymbol\tbasemean\tlog2fc\tp\tp_adj\n"
        "ENSGALG001\tGENE1\t100\t2.5\t1e-10\t1e-9\n"
        "ENSGALG002\tGENE2\t50\t-1.8\t1e-8\t1e-7\n"
    )
    df, rpt = pio.load_deseq2(io.StringIO(csv))
    assert "gene_id" in df.columns
    assert "log2FoldChange" in df.columns
    assert "padj" in df.columns
    assert rpt.n_rows == 2


def test_missing_required_columns_reports_cleanly():
    csv = "gene_id\tlog2FoldChange\nG1\t1.0\n"
    df, rpt = pio.load_deseq2(io.StringIO(csv))
    assert rpt.missing_required
    assert any("padj" in m or "pvalue" in m for m in rpt.missing_required)
