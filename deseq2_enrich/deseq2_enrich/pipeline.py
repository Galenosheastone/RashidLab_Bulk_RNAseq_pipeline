"""End-to-end orchestration shared by the CLI and the Streamlit app.

``run_contrast`` executes the whole flow for one DESeq2 table:

    load -> DEG selection -> ORA (native chicken, g:Profiler)
         -> ortholog map -> ranked list -> GSEA (MSigDB via gseapy)

Each stage is optional and independently catchable so a network hiccup in one
service does not lose the others. The returned bundle is plain data
(dataframes + the gseapy result object) so the app can cache and plot it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import pandas as pd

from . import config, io, degs, rank, ortho, ora, genesets, gsea


@dataclass
class ContrastResult:
    name: str
    df: pd.DataFrame
    report: object
    deg_sets: object
    ora: pd.DataFrame = field(default_factory=pd.DataFrame)
    gsea: Optional[object] = None          # gsea.GSEAResult
    ranking: Optional[pd.Series] = None
    errors: dict = field(default_factory=dict)


def run_contrast(
    path_or_buffer,
    contrast_name: str = "contrast_1",
    *,
    padj_threshold: float = config.PADJ_THRESHOLD,
    lfc_threshold: float = config.LFC_THRESHOLD,
    id_col: str = "gene_id",
    rank_metric: str = config.RANK_METRIC,
    ora_sources: Optional[list[str]] = None,
    ora_directions: tuple[str, ...] = ("up", "down"),
    gsea_libraries: Optional[list[str]] = None,
    custom_gmt: Optional[dict] = None,
    organism: str = config.ORGANISM,
    do_ora: bool = True,
    do_gsea: bool = True,
    gsea_permutations: int = config.GSEA_PERMUTATIONS,
) -> ContrastResult:
    df, report = io.load_deseq2(path_or_buffer, contrast_name=contrast_name)
    if report.missing_required:
        return ContrastResult(contrast_name, df, report, None,
                              errors={"load": report.missing_required})

    deg_sets = degs.select_degs(df, padj_threshold, lfc_threshold, id_col)
    result = ContrastResult(contrast_name, df, report, deg_sets)

    # --- ORA: native chicken via g:Profiler -------------------------------
    if do_ora:
        try:
            result.ora = ora.run_ora_directional(
                deg_sets,
                sources=ora_sources or config.ORA_DEFAULT_SOURCES,
                organism=organism,
                directions=ora_directions,
            )
        except (KeyError, ValueError, RuntimeError, AssertionError, TypeError,
                OSError, ConnectionError) as exc:
            import traceback

            result.errors["ora"] = f"{type(exc).__name__}: {exc}" or repr(exc)
            result.errors["ora_traceback"] = traceback.format_exc()

    # --- GSEA: ortholog -> ranked human symbols -> gseapy -----------------
    if do_gsea:
        try:
            gene_sets = {}
            if gsea_libraries:
                gene_sets.update(
                    genesets.combine_libraries(gsea_libraries, config.ORTHOLOG_TARGET)
                )
            if custom_gmt:
                gene_sets.update(custom_gmt)
            if not gene_sets:
                raise ValueError("No gene sets selected for GSEA.")

            mapped = ortho.attach_human_symbol(df, id_col=id_col, source=organism)
            ranking = rank.build_rank(mapped, metric=rank_metric,
                                      key_col="human_symbol")
            result.ranking = ranking
            result.gsea = gsea.run_prerank(
                ranking, gene_sets,
                min_size=config.GSEA_MIN_SIZE, max_size=config.GSEA_MAX_SIZE,
                permutations=gsea_permutations, seed=config.GSEA_SEED,
            )
        except (KeyError, ValueError, RuntimeError, AssertionError, TypeError,
                OSError, ConnectionError) as exc:
            import traceback

            result.errors["gsea"] = f"{type(exc).__name__}: {exc}" or repr(exc)
            result.errors["gsea_traceback"] = traceback.format_exc()

    return result
