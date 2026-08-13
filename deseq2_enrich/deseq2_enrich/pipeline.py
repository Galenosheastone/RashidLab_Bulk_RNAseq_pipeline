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
    gsea_mode: str = "ortholog",
    chicken_gmt_keying: str = config.CHICKEN_GMT_KEYING,
    min_size: int = config.GSEA_MIN_SIZE,
    max_size: int = config.GSEA_MAX_SIZE,
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
            gsea_mode=gsea_mode,
            chicken_gmt_keying=chicken_gmt_keying,
            min_size=min_size,
            max_size=max_size,
        )

    return result


def _partition_libraries(
    libraries: list[str], gsea_mode: str
) -> tuple[set[str], list[str]]:
    """Split selected collections into (native chicken sources, ortholog libs).

    In ``auto`` a collection goes native whenever chicken annotations exist for
    it (GO/Reactome/WikiPathways); Hallmark and Oncogenic have no chicken
    equivalent and keep the ortholog route. KEGG has no native GMT either --
    g:Profiler does not redistribute it -- so it also falls to orthologs.
    """
    if gsea_mode == "ortholog":
        return set(), list(libraries)
    if gsea_mode == "native_chicken":
        native = {
            config.LIBRARY_NATIVE_SOURCE[lib]
            for lib in libraries
            if lib in config.LIBRARY_NATIVE_SOURCE
        }
        # Nothing selected that maps natively -> fall back to the defaults so
        # the user still gets the chicken collections they asked for.
        return (native or set(config.NATIVE_GSEA_DEFAULT_SOURCES)), []

    native, human_only = set(), []
    for lib in libraries:
        source = config.LIBRARY_NATIVE_SOURCE.get(lib)
        if source is not None:
            native.add(source)
        else:
            human_only.append(lib)
    return native, human_only


def _native_key_column(df: pd.DataFrame, keying: str) -> str:
    """Pick the dataframe column matching the chicken GMT's ID space.

    ``ensg`` keying needs an Ensembl column; if the export has none we fall
    back to symbols rather than silently ranking on an ID space the GMT does
    not use, which would score nothing at all.
    """
    if keying == "ensg":
        if "ensembl_id" in df.columns and df["ensembl_id"].notna().any():
            return "ensembl_id"
        if "gene_name" in df.columns and df["gene_name"].notna().any():
            return "gene_name"
        raise ValueError(
            "Native chicken GSEA with Ensembl keying needs an 'ensembl_id' or "
            "'gene_name' column; this table has neither."
        )
    for col in ("gene_name", "gene_id"):
        if col in df.columns and df[col].notna().any():
            return col
    raise ValueError(
        "Native chicken GSEA with symbol keying needs a 'gene_name' or "
        "'gene_id' column carrying chicken gene symbols."
    )


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
    gsea_mode: str = "ortholog",
    chicken_gmt_keying: str = config.CHICKEN_GMT_KEYING,
    min_size: int = config.GSEA_MIN_SIZE,
    max_size: int = config.GSEA_MAX_SIZE,
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
        if gsea_mode not in config.GSEA_MODES:
            raise ValueError(
                f"Unknown gsea_mode {gsea_mode!r}; expected one of {config.GSEA_MODES}"
            )
        libraries = list(gsea_libraries or [])
        native_sources, ortholog_libs = _partition_libraries(libraries, gsea_mode)

        routes: dict[str, gsea.GSEAResult] = {}
        tables: list[pd.DataFrame] = []
        requested: dict[str, int] = {}

        # --- Native chicken route: no ortholog step, no dropout bias -------
        if native_sources or (gsea_mode == "native_chicken" and custom_gmt):
            native_sets = {}
            if native_sources:
                native_sets.update(
                    genesets.fetch_chicken_gmt(
                        keying=chicken_gmt_keying, sources=tuple(sorted(native_sources))
                    )
                )
            if custom_gmt:
                native_sets.update(custom_gmt)
            if native_sets:
                key_col = _native_key_column(result.df, chicken_gmt_keying)
                native_rank = rank.build_rank(
                    result.df, metric=rank_metric, key_col=key_col
                )
                requested["native"] = len(native_sets)
                routes["native"] = gsea.run_prerank(
                    native_rank, native_sets,
                    min_size=min_size, max_size=max_size,
                    permutations=gsea_permutations, seed=config.GSEA_SEED,
                )
                result.gsea_metadata["native_key_col"] = key_col
                result.gsea_metadata["native_sources"] = sorted(native_sources)
                result.gsea_metadata["native_ranking_size"] = len(native_rank)

        # --- Ortholog route: human-only collections (Hallmark / Oncogenic) --
        if ortholog_libs or (gsea_mode == "ortholog" and custom_gmt):
            ortho_sets = {}
            if ortholog_libs:
                ortho_sets.update(
                    genesets.combine_libraries(ortholog_libs, config.ORTHOLOG_TARGET)
                )
            if custom_gmt and gsea_mode == "ortholog":
                ortho_sets.update(custom_gmt)
            if ortho_sets:
                mapped, ortho_report = ortho.attach_human_symbol(
                    result.df,
                    id_col=id_col,
                    source=organism,
                    strict_one_to_one=strict_one_to_one,
                )
                result.gsea_metadata["ortholog_report"] = ortho_report.as_text()
                result.gsea_metadata["mapping_rate"] = round(
                    ortho_report.mapping_rate, 3
                )

                # A pre-ranked GSEA takes the ranked list as its universe. If
                # most of the chicken genes dropped out, the null is computed
                # on a biased remnant and the NES/FDR are not interpretable --
                # fail loudly rather than hand back confident-looking numbers.
                if (ortho_report.n_query_mapped < config.ORTHO_MIN_MAPPED_GENES
                        or ortho_report.mapping_rate < config.ORTHO_MIN_MAPPING_RATE):
                    raise ValueError(
                        "Ortholog mapping too sparse for reliable GSEA "
                        f"({ortho_report.as_text()}). "
                        "Check the ID column (prefer ENSGALG Ensembl IDs), or "
                        "use native-chicken GSEA for GO/Reactome/WP instead."
                    )

                ortho_rank = rank.build_rank(
                    mapped, metric=rank_metric, key_col="human_symbol"
                )
                requested["ortholog"] = len(ortho_sets)
                routes["ortholog"] = gsea.run_prerank(
                    ortho_rank, ortho_sets,
                    min_size=min_size, max_size=max_size,
                    permutations=gsea_permutations, seed=config.GSEA_SEED,
                )
                result.gsea_metadata["ortholog_ranking_size"] = len(ortho_rank)

        if not routes:
            raise ValueError("No gene sets selected for GSEA.")

        # FDR is computed within each prerank and is NOT pooled across routes:
        # the two runs have different nulls and different gene universes, so a
        # combined q-value would be meaningless. The route is carried in the
        # table so the distinction stays visible downstream.
        for route_name, route_result in routes.items():
            tbl = route_result.table.copy()
            tbl["gsea_route"] = route_name
            tables.append(tbl)

        primary = routes.get("native") or routes["ortholog"]
        result.gsea = gsea.GSEAResult(
            table=pd.concat(tables, ignore_index=True) if len(tables) > 1
            else tables[0],
            raw=primary.raw,
            ranking=primary.ranking,
            running_terms=tuple(
                t for r in routes.values() for t in (r.running_terms or ())
            ),
            routes=routes if len(routes) > 1 else {},
        )
        result.ranking = primary.ranking

        # update(), not reassign: the ortholog report was recorded above and
        # must survive into the metadata the UI renders.
        result.gsea_metadata.update({
            "libraries": libraries,
            "custom_gmt_terms": len(custom_gmt or {}),
            "id_col": id_col,
            "rank_metric": rank_metric,
            "organism": organism,
            "permutations": gsea_permutations,
            "ranking_size": len(result.ranking),
            "strict_one_to_one": strict_one_to_one,
            # gseapy drops any set whose intersection with the ranked list
            # falls outside [min_size, max_size], silently. Reporting requested
            # vs scored is the only way to see how much was actually tested.
            "gene_sets_requested": sum(requested.values()),
            "gene_sets_scored": int(len(result.gsea.table)),
            "gene_sets_requested_by_route": dict(requested),
            "gene_sets_scored_by_route": {
                name: int(len(r.table)) for name, r in routes.items()
            },
            "min_size": min_size,
            "max_size": max_size,
            "gsea_mode": gsea_mode,
            "gsea_routes": sorted(routes),
            "chicken_gmt_keying": chicken_gmt_keying,
            "fdr_per_route": len(routes) > 1,
        })
    except (KeyError, ValueError, RuntimeError, AssertionError, TypeError,
            OSError, ConnectionError) as exc:
        import traceback

        result.errors["gsea"] = f"{type(exc).__name__}: {exc}" or repr(exc)
        result.errors["gsea_traceback"] = traceback.format_exc()

    return result
