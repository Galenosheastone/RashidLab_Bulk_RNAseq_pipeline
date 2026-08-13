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
    gsea_metadata: dict = field(default_factory=dict)
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
    strict_one_to_one: bool = False,
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

    if do_gsea:
        run_gsea_for_result(
            result,
            id_col=id_col,
            rank_metric=rank_metric,
            gsea_libraries=gsea_libraries,
            custom_gmt=custom_gmt,
            organism=organism,
            gsea_permutations=gsea_permutations,
            strict_one_to_one=strict_one_to_one,
        )

    return result


def run_gsea_for_result(
    result: ContrastResult,
    *,
    id_col: str = "gene_id",
    rank_metric: str = config.RANK_METRIC,
    gsea_libraries: Optional[list[str]] = None,
    custom_gmt: Optional[dict] = None,
    organism: str = config.ORGANISM,
    gsea_permutations: int = config.GSEA_PERMUTATIONS,
    strict_one_to_one: bool = False,
) -> ContrastResult:
    """Run or rerun only the GSEA stage on an existing contrast result.

    ``strict_one_to_one`` keeps only unambiguous chicken<->human ortholog
    pairs. It defaults to False to preserve existing behaviour for
    programmatic callers; the app exposes it as a toggle.
    """
    result.errors.pop("gsea", None)
    result.errors.pop("gsea_traceback", None)
    result.gsea = None
    result.ranking = None
    result.gsea_metadata = {}

    try:
        gene_sets = {}
        libraries = list(gsea_libraries or [])
        if libraries:
            gene_sets.update(
                genesets.combine_libraries(libraries, config.ORTHOLOG_TARGET)
            )
        if custom_gmt:
            gene_sets.update(custom_gmt)
        if not gene_sets:
            raise ValueError("No gene sets selected for GSEA.")

        mapped, ortho_report = ortho.attach_human_symbol(
            result.df,
            id_col=id_col,
            source=organism,
            strict_one_to_one=strict_one_to_one,
        )
        result.gsea_metadata["ortholog_report"] = ortho_report.as_text()
        result.gsea_metadata["mapping_rate"] = round(ortho_report.mapping_rate, 3)

        # A pre-ranked GSEA takes the ranked list as its universe. If most of
        # the chicken genes dropped out, the null is computed on a biased
        # remnant and the NES/FDR are not interpretable -- fail loudly rather
        # than hand back confident-looking numbers.
        if (ortho_report.n_query_mapped < config.ORTHO_MIN_MAPPED_GENES
                or ortho_report.mapping_rate < config.ORTHO_MIN_MAPPING_RATE):
            raise ValueError(
                "Ortholog mapping too sparse for reliable GSEA "
                f"({ortho_report.as_text()}). "
                "Check the ID column (prefer ENSGALG Ensembl IDs), or run "
                "native-chicken GSEA for GO/KEGG/Reactome/WP instead."
            )

        ranking = rank.build_rank(mapped, metric=rank_metric, key_col="human_symbol")
        result.ranking = ranking
        result.gsea = gsea.run_prerank(
            ranking, gene_sets,
            min_size=config.GSEA_MIN_SIZE, max_size=config.GSEA_MAX_SIZE,
            permutations=gsea_permutations, seed=config.GSEA_SEED,
        )
        # update(), not reassign: the ortholog report was recorded above and
        # must survive into the metadata the UI renders.
        result.gsea_metadata.update({
            "libraries": libraries,
            "custom_gmt_terms": len(custom_gmt or {}),
            "gene_sets_requested": len(gene_sets),
            "id_col": id_col,
            "rank_metric": rank_metric,
            "organism": organism,
            "permutations": gsea_permutations,
            "ranking_size": len(ranking),
            "strict_one_to_one": strict_one_to_one,
        })
    except (KeyError, ValueError, RuntimeError, AssertionError, TypeError,
            OSError, ConnectionError) as exc:
        import traceback

        result.errors["gsea"] = f"{type(exc).__name__}: {exc}" or repr(exc)
        result.errors["gsea_traceback"] = traceback.format_exc()

    return result
