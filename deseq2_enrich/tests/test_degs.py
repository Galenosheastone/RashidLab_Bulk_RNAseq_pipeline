from deseq2_enrich import degs


def test_universe_excludes_na_padj(toy_deseq2):
    d = degs.select_degs(toy_deseq2, padj_threshold=0.05, lfc_threshold=1.0)
    assert len(d.universe) == 14  # 15 tested rows, with one duplicate gene_id
    na_ids = set(toy_deseq2.loc[toy_deseq2["padj"].isna(), "gene_id"])
    assert not (set(d.universe) & na_ids)


def test_directional_split(toy_deseq2):
    d = degs.select_degs(toy_deseq2, padj_threshold=0.05, lfc_threshold=1.0)
    assert len(d.up) == 3
    assert len(d.down) == 3
    assert set(d.up).isdisjoint(set(d.down))


def test_all_sig_is_union(toy_deseq2):
    d = degs.select_degs(toy_deseq2, padj_threshold=0.05, lfc_threshold=1.0)
    assert set(d.all_sig) == set(d.up) | set(d.down)
