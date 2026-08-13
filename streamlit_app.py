"""Streamlit front-end for the DESeq2 -> ORA/GSEA pipeline.

Design notes
------------
* Enrichment is computed once, on an explicit "Run" click, and stored in
  ``st.session_state`` so tweaking a plot control (top-N, term selection) never
  re-hits the network or re-runs GSEA.
* Network calls (g:Profiler, Enrichr libraries) are memoised at module level in
  the package, so repeated runs in a session are cheap.
* Uploaded tables are held in memory only; nothing is written server-side.
"""
from __future__ import annotations

import gc
import io as _io
import os
import shutil
import sys
import json
import zipfile
from datetime import datetime

import pandas as pd
import streamlit as st

APP_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_PARENT = os.path.join(APP_DIR, "deseq2_enrich")
if not os.path.isdir(os.path.join(PACKAGE_PARENT, "deseq2_enrich")):
    PACKAGE_PARENT = os.path.dirname(APP_DIR)


def _purge_project_bytecode(root: str) -> None:
    """Avoid stale .pyc files after Streamlit Cloud redeploys."""
    sys.dont_write_bytecode = True
    for dirpath, dirnames, _ in os.walk(root):
        if "__pycache__" not in dirnames:
            continue
        cache_dir = os.path.join(dirpath, "__pycache__")
        shutil.rmtree(cache_dir, ignore_errors=True)
        dirnames.remove("__pycache__")


_purge_project_bytecode(APP_DIR)
sys.path.insert(0, PACKAGE_PARENT)

from deseq2_enrich import config, plots, genesets, ora, ortho  # noqa: E402
from deseq2_enrich.pipeline import run_contrast, run_gsea_for_result, ContrastResult  # noqa: E402

st.set_page_config(page_title="DESeq2 Enrichment", layout="wide",
                   page_icon="🧬")

SAMPLE_PATH = os.path.join(
    PACKAGE_PARENT,
    "sample_data", "DESeq2_Sacral_vs_Sacralized_Caudal_all.tsv",
)


def _stored_params(p: dict) -> dict:
    return {
        k: v for k, v in p.items()
        if k not in ("uploads", "custom_gmt_file", "run")
    }


def _clear_package_caches() -> None:
    for module in (genesets, ora, ortho):
        clear = getattr(module, "clear_cache", None)
        if clear is not None:
            clear()


def _clear_streamlit_caches() -> None:
    for attr in ("cache_data", "cache_resource"):
        cache = getattr(st, attr, None)
        clear = getattr(cache, "clear", None)
        if clear is not None:
            clear()


def _clear_analysis_memory(clear_results: bool = False) -> None:
    if clear_results:
        st.session_state.pop("results", None)
        st.session_state.pop("params", None)
    _clear_package_caches()
    _clear_streamlit_caches()
    gc.collect()


def _cache_status() -> str:
    parts = []
    for label, module in (
        ("Gene sets", genesets),
        ("ORA", ora),
        ("Orthologs", ortho),
    ):
        info = getattr(module, "cache_info", lambda: None)()
        if info is not None:
            parts.append(f"{label}: {info.currsize}/{info.maxsize}")
    return " | ".join(parts) if parts else "No package caches active"


# --------------------------------------------------------------------------
# Sidebar — inputs and parameters
# --------------------------------------------------------------------------
def sidebar() -> dict:
    st.sidebar.title("🧬 DESeq2 Enrichment")
    st.sidebar.caption(
        "Native chicken throughout: ORA via g:Profiler, GSEA on chicken "
        "GO / Reactome / WikiPathways. Orthologs only for human-only collections."
    )

    st.sidebar.subheader("1 · Data")
    use_sample = st.sidebar.toggle("Use bundled sample data", value=True,
                                   help="The synthetic Sacral vs Caudal DESeq2 table.")
    uploads = st.sidebar.file_uploader(
        "Or upload DESeq2 result table(s)", type=["tsv", "csv", "txt"],
        accept_multiple_files=True,
        help="One file per contrast. Standard DESeq2 columns are auto-detected.",
    )

    st.sidebar.subheader("2 · DEG thresholds")
    padj = st.sidebar.number_input("padj <", value=config.PADJ_THRESHOLD,
                                   min_value=0.0, max_value=1.0, step=0.01, format="%.3f")
    lfc = st.sidebar.number_input("|log2FC| >", value=config.LFC_THRESHOLD,
                                  min_value=0.0, step=0.25)
    id_col = st.sidebar.selectbox("Gene ID column", ["gene_id", "entrez_id", "gene_name"],
                                  help="Identifier passed to g:Profiler.")

    st.sidebar.subheader("3 · GSEA ranking")
    rank_metric = st.sidebar.selectbox(
        "Ranking metric", config.RANK_METRICS, index=0,
        format_func=lambda m: {"stat": "Wald stat (recommended)",
                               "signed_logp": "sign(FC)·-log10 p",
                               "log2fc": "log2 fold change"}[m],
    )

    st.sidebar.subheader("4 · ORA sources (native chicken)")
    ora_sources = st.sidebar.multiselect(
        "g:Profiler sources", list(config.ORA_SOURCES.keys()),
        default=config.ORA_DEFAULT_SOURCES,
        format_func=lambda s: config.ORA_SOURCES[s],
    )
    ora_dirs = st.sidebar.multiselect("ORA directions", ["up", "down", "all"],
                                      default=["up", "down"])

    st.sidebar.subheader("5 · Compute")
    quick = st.sidebar.toggle("Quick mode (fewer permutations)", value=True,
                              help="1000 → 100 permutations for a fast preview.")
    do_ora = st.sidebar.checkbox("Run ORA", value=True)
    do_gsea = st.sidebar.checkbox("Run GSEA", value=True)
    try:
        run = st.sidebar.button("▶ Run enrichment", type="primary", width="stretch")
    except TypeError:
        run = st.sidebar.button("▶ Run enrichment", type="primary",
                                use_container_width=True)

    st.sidebar.caption("Uploaded data is processed in memory and not stored.")
    with st.sidebar.expander("Memory & cache", expanded=False):
        st.caption(_cache_status())
        if st.button("Clear cached analysis data", key="clear_cached_analysis"):
            _clear_analysis_memory(clear_results=True)
            st.success("Cleared cached results and network caches.")

    return dict(
        use_sample=use_sample, uploads=uploads, padj=padj, lfc=lfc, id_col=id_col,
        rank_metric=rank_metric, ora_sources=ora_sources, ora_dirs=tuple(ora_dirs),
        gsea_libs=st.session_state.get("gsea_libs", config.GSEA_DEFAULT_LIBRARIES),
        custom_gmt_file=None,
        permutations=100 if quick else config.GSEA_PERMUTATIONS,
        do_ora=do_ora, do_gsea=do_gsea, run=run,
    )


# --------------------------------------------------------------------------
# Run + cache in session_state
# --------------------------------------------------------------------------
def execute(p: dict) -> list[ContrastResult]:
    inputs: list[tuple[str, object]] = []
    if p["uploads"]:
        for f in p["uploads"]:
            name = os.path.splitext(f.name)[0]
            inputs.append((name, f))
    elif p["use_sample"]:
        inputs.append(("Sacral_vs_Sacralized_Caudal", SAMPLE_PATH))
    if not inputs:
        st.warning("Upload a DESeq2 table or enable the sample data.")
        return []

    custom = None
    if p["custom_gmt_file"] is not None:
        lines = p["custom_gmt_file"].getvalue().decode("utf-8", "replace").splitlines()
        custom = genesets.load_gmt(lines)

    results = []
    prog = st.progress(0.0, text="Starting…")
    for i, (name, src) in enumerate(inputs):
        prog.progress(i / len(inputs), text=f"Enriching {name}…")
        res = run_contrast(
            src, contrast_name=name,
            padj_threshold=p["padj"], lfc_threshold=p["lfc"], id_col=p["id_col"],
            rank_metric=p["rank_metric"], ora_sources=p["ora_sources"],
            ora_directions=p["ora_dirs"], gsea_libraries=p["gsea_libs"],
            custom_gmt=custom, do_ora=p["do_ora"], do_gsea=p["do_gsea"],
            gsea_permutations=p["permutations"],
            # Chicken-first: native gene sets wherever they exist, orthologs
            # only for the human-only collections.
            gsea_mode=p.get("gsea_mode", config.GSEA_DEFAULT_MODE),
            strict_one_to_one=p.get("strict_one_to_one", True),
        )
        results.append(res)
    prog.progress(1.0, text="Done")
    prog.empty()
    return results


def _gene_columns(df: pd.DataFrame) -> list[str]:
    preferred = ["gene_id", "gene_name", "entrez_id"]
    return [c for c in preferred if c in df.columns]


def _parse_gene_text(text: str) -> list[str]:
    text = (text or "").replace(",", "\n").replace(";", "\n")
    return [g.strip() for g in text.splitlines() if g.strip()]


def _gene_label(row: pd.Series) -> str:
    parts = []
    for col in ("gene_name", "gene_id", "entrez_id"):
        if col in row.index and pd.notna(row[col]):
            value = str(row[col])
            if value and value not in parts:
                parts.append(value)
    return " / ".join(parts) if parts else str(row.name)


def _top_gene_ids(df: pd.DataFrame, preset: str, n: int) -> list[str]:
    d = df[df["padj"].notna()].copy()
    if d.empty:
        return []
    if preset == "Top up":
        d = d[d["log2FoldChange"] > 0].sort_values(
            ["padj", "log2FoldChange"], ascending=[True, False]
        )
    elif preset == "Top down":
        d = d[d["log2FoldChange"] < 0].sort_values(
            ["padj", "log2FoldChange"], ascending=[True, True]
        )
    elif preset == "Most significant":
        d = d.sort_values("padj")
    else:
        return []
    id_col = "gene_name" if "gene_name" in d.columns else "gene_id"
    return d[id_col].dropna().astype(str).head(n).tolist()


def _parse_ora_genes(value) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        return [str(v) for v in value if str(v)]
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    for ch in "[](){}'\"":
        text = text.replace(ch, "")
    sep = ";" if ";" in text else "," if "," in text else "|"
    return [part.strip() for part in text.split(sep) if part.strip()]


def _ora_label(row: pd.Series) -> str:
    term = row.get("term_name", row.get("term_id", "term"))
    source = row.get("source", "ORA")
    direction = row.get("direction", "all")
    return f"{term} ({source}, {direction})"


def _plotly(fig, key: str):
    try:
        st.plotly_chart(fig, width="stretch", key=key)
    except TypeError:
        st.plotly_chart(fig, use_container_width=True, key=key)


def _dataframe(df: pd.DataFrame, *, key: str | None = None,
               height: int | None = None):
    kwargs = {"key": key} if key else {}
    if height is not None:
        kwargs["height"] = height
    try:
        st.dataframe(df, width="stretch", **kwargs)
    except TypeError:
        st.dataframe(df, use_container_width=True, **kwargs)


def _gsea_library_options() -> list[str]:
    ordered = []
    for names in config.GSEA_LIBRARY_GROUPS.values():
        ordered.extend(name for name in names if name in config.GSEA_LIBRARIES)
    ordered.extend(name for name in config.GSEA_LIBRARIES if name not in ordered)
    return ordered


def _format_gsea_library(name: str) -> str:
    label = config.GSEA_LIBRARIES.get(name, name)
    for group, names in config.GSEA_LIBRARY_GROUPS.items():
        if name in names:
            return f"{group}: {label}"
    return label


def _parse_custom_gmt_upload(upload) -> tuple[dict | None, str | None]:
    if upload is None:
        return None, None
    lines = upload.getvalue().decode("utf-8", "replace").splitlines()
    custom = genesets.load_gmt(lines)
    return custom, upload.name


def _gsea_term_label(row: pd.Series) -> str:
    collection = str(row.get("collection", "GSEA"))
    term = str(row.get("term_short", row.get("term", "")))
    nes = row.get("NES", float("nan"))
    fdr = row.get("fdr", float("nan"))
    return f"{collection}: {term}  |  NES {nes:.2f}, FDR {fdr:.2g}"


def _gsea_terms_with_running_data(gsea_result, terms: list[str]) -> list[str]:
    retained = set(getattr(gsea_result, "running_terms", ()) or ())
    raw_results = getattr(getattr(gsea_result, "raw", None), "results", {}) or {}
    available = []
    for term in terms:
        payload = raw_results.get(term, {})
        if term in retained or ("RES" in payload and "hits" in payload):
            available.append(term)
    return available


def _render_gsea_controls(res: ContrastResult, params: dict) -> None:
    st.markdown("**Gene-set collections**")
    default_libs = res.gsea_metadata.get(
        "libraries",
        params.get("gsea_libs", config.GSEA_DEFAULT_LIBRARIES),
    )
    selected_libs = st.multiselect(
        "Collections for GSEA",
        _gsea_library_options(),
        default=[lib for lib in default_libs if lib in config.GSEA_LIBRARIES],
        format_func=_format_gsea_library,
        key=f"gsea_libs_{res.name}",
        help=(
            "Under the default Auto route, GO / Reactome / WikiPathways selections "
            "are served from native chicken gene sets — the Enrichr library name is "
            "just how you pick the collection. Hallmark and Oncogenic Signatures "
            "have no chicken equivalent and run through ortholog mapping instead."
        ),
    )
    custom_gmt_file = st.file_uploader(
        "Custom GMT gene sets",
        type=["gmt"],
        key=f"gsea_custom_gmt_{res.name}",
        help="Optional pathway/module file. Terms are combined with selected collections.",
    )
    slow = [lib for lib in selected_libs if lib in config.GSEA_SLOW_LIBRARIES]
    if slow:
        perms = params.get("permutations", config.GSEA_PERMUTATIONS)
        st.error(
            f"**{', '.join(config.GSEA_LIBRARIES.get(l, l) for l in slow)}** "
            "selected. Native chicken GO is ~18,700 terms and scores ~5,900 sets "
            "on a typical table. Measured peak memory: **~4.2 GB**, against the "
            "~1 GB limit on Streamlit Community Cloud — the app will most likely "
            "be killed or thrash rather than finish. Runtime is ~6 min at 1000 "
            f"permutations even with enough RAM (currently {perms}).\n\n"
            "Reactome + WikiPathways score 570 sets at ~830 MB in under a minute. "
            "Run GO locally, or deselect it here.",
            icon="🛑",
        )

    mode_labels = {
        "auto": "Auto — native chicken where possible (recommended)",
        "native_chicken": "Native chicken only",
        "ortholog": "Human orthologs only (legacy)",
    }
    mode = st.radio(
        "GSEA route",
        list(mode_labels),
        format_func=lambda m: mode_labels[m],
        index=list(mode_labels).index(
            res.gsea_metadata.get("gsea_mode", config.GSEA_DEFAULT_MODE)
        ),
        key=f"gsea_mode_{res.name}",
        help=(
            "Native chicken runs GSEA directly on chicken GO/Reactome/WikiPathways "
            "annotations — no ortholog step, so no genes are dropped for lacking a "
            "human counterpart. Hallmark and Oncogenic Signatures are human-only "
            "collections and always use orthologs. KEGG has no downloadable chicken "
            "gene sets; it remains available in the ORA tab."
        ),
    )
    strict = st.checkbox(
        "Strict 1:1 orthologs (cleaner, fewer genes)",
        value=res.gsea_metadata.get("strict_one_to_one", True),
        key=f"gsea_strict_{res.name}",
        help=(
            "Only unambiguous chicken↔human pairs. Prevents one chicken "
            "measurement appearing as several correlated human paralog ranks. "
            "Applies to the ortholog route only."
        ),
    )
    s1, s2 = st.columns(2)
    with s1:
        min_size = st.number_input(
            "Min gene-set size", min_value=2, max_value=100,
            value=int(res.gsea_metadata.get("min_size", config.GSEA_MIN_SIZE)),
            step=1, key=f"gsea_min_size_{res.name}",
            help=(
                "Sets are filtered AFTER intersecting with your ranked list, so "
                "a high floor silently discards valid small pathways. Sets "
                "outside these bounds are dropped without warning by gseapy."
            ),
        )
    with s2:
        max_size = st.number_input(
            "Max gene-set size", min_value=50, max_value=5000,
            value=int(res.gsea_metadata.get("max_size", config.GSEA_MAX_SIZE)),
            step=50, key=f"gsea_max_size_{res.name}",
        )
    c1, c2 = st.columns([1, 3])
    with c1:
        update = st.button("Update GSEA", type="primary",
                           key=f"gsea_update_{res.name}")
    with c2:
        if selected_libs:
            labels = [config.GSEA_LIBRARIES.get(lib, lib) for lib in selected_libs]
            st.caption("Selected: " + ", ".join(labels))
        else:
            st.caption("Selected: custom GMT only")
    if not update:
        return

    custom_gmt, custom_name = _parse_custom_gmt_upload(custom_gmt_file)
    # Native mode falls back to the default chicken sources when nothing
    # selected maps natively, so it does not need an explicit selection.
    if not selected_libs and not custom_gmt and mode != "native_chicken":
        st.warning("Select at least one collection or upload a GMT file.")
        return

    with st.spinner("Updating GSEA for this contrast..."):
        run_gsea_for_result(
            res,
            id_col=params.get("id_col", "gene_id"),
            rank_metric=params.get("rank_metric", config.RANK_METRIC),
            gsea_libraries=selected_libs,
            custom_gmt=custom_gmt,
            organism=config.ORGANISM,
            gsea_permutations=params.get("permutations", config.GSEA_PERMUTATIONS),
            gsea_mode=mode,
            strict_one_to_one=strict,
            min_size=int(min_size),
            max_size=int(max_size),
        )
    gc.collect()
    st.session_state["gsea_libs"] = selected_libs
    params["gsea_libs"] = selected_libs
    params["gsea_custom_gmt_name"] = custom_name
    params["do_gsea"] = True
    st.session_state["params"] = params
    if "gsea" in res.errors:
        st.error(f"GSEA failed: {res.errors['gsea']}")
    else:
        st.success("GSEA updated for this contrast.")


# --------------------------------------------------------------------------
# Tabs
# --------------------------------------------------------------------------
def tab_qc(res: ContrastResult, p: dict):
    st.subheader(f"Quality control — {res.name}")
    c1, c2, c3, c4 = st.columns(4)
    counts = res.deg_sets.counts
    c1.metric("Tested genes", res.report.n_tested)
    c2.metric("Up", counts["up"])
    c3.metric("Down", counts["down"])
    c4.metric("NA padj (dropped)", res.report.n_na_padj)
    with st.expander("Coverage & load report", expanded=False):
        st.text(res.report.as_text())
    left, right = st.columns(2)
    with left:
        _plotly(plots.de_count_bar(res.df, p["padj"], p["lfc"]),
                key=f"qc_de_count_{res.name}")
    with right:
        _plotly(plots.ma_plot(res.df, p["padj"], p["lfc"]),
                key=f"qc_ma_{res.name}")


def tab_genes(res: ContrastResult, p: dict):
    st.subheader(f"Gene explorer — {res.name}")
    d = plots.de_plot_table(res.df, p["padj"], p["lfc"])
    gene_cols = _gene_columns(d)
    if not gene_cols:
        st.info("No gene identifier columns were detected.")
        return

    filters = st.columns([1.5, 1, 1, 1])
    with filters[0]:
        query = st.text_input("Search genes", key=f"gene_query_{res.name}")
    with filters[1]:
        statuses = st.multiselect(
            "Status", ["Up", "Down", "n.s."],
            default=["Up", "Down", "n.s."], key=f"gene_status_{res.name}"
        )
    with filters[2]:
        max_padj = st.number_input(
            "Max padj", min_value=0.0, max_value=1.0, value=1.0,
            step=0.01, format="%.3f", key=f"gene_padj_{res.name}"
        )
    with filters[3]:
        min_abs_lfc = st.number_input(
            "Min |log2FC|", min_value=0.0, value=0.0, step=0.25,
            key=f"gene_lfc_{res.name}"
        )

    filtered = d.copy()
    if query:
        q = query.strip().lower()
        mask = pd.Series(False, index=filtered.index)
        for col in gene_cols:
            mask |= filtered[col].astype(str).str.lower().str.contains(q, na=False)
        filtered = filtered[mask]
    if statuses:
        filtered = filtered[filtered["de_status"].isin(statuses)]
    filtered = filtered[
        (filtered["padj"].isna() | (filtered["padj"] <= max_padj))
        & (filtered["log2FoldChange"].abs() >= min_abs_lfc)
    ]

    sort_choice = st.selectbox(
        "Sort genes",
        ["padj ascending", "|log2FC| descending", "baseMean descending"],
        key=f"gene_sort_{res.name}",
    )
    if sort_choice == "|log2FC| descending":
        filtered = filtered.assign(abs_lfc=filtered["log2FoldChange"].abs()).sort_values(
            "abs_lfc", ascending=False
        )
    elif sort_choice == "baseMean descending" and "baseMean" in filtered.columns:
        filtered = filtered.sort_values("baseMean", ascending=False)
    else:
        filtered = filtered.sort_values("padj", na_position="last")

    c1, c2, c3 = st.columns(3)
    c1.metric("Filtered genes", f"{len(filtered):,}")
    c2.metric("Up", int((filtered["de_status"] == "Up").sum()))
    c3.metric("Down", int((filtered["de_status"] == "Down").sum()))

    show_cols = [
        c for c in [
            "gene_id", "gene_name", "entrez_id", "de_status", "baseMean",
            "log2FoldChange", "stat", "pvalue", "padj", "neg_log10_padj",
        ] if c in filtered.columns
    ]
    _dataframe(filtered[show_cols], key=f"gene_table_{res.name}", height=360)
    st.download_button(
        "⬇ Filtered genes (CSV)",
        filtered[show_cols].to_csv(index=False),
        file_name=f"{res.name}_filtered_genes.csv",
        mime="text/csv",
    )

    if filtered.empty:
        return
    options = filtered.head(2500).index.tolist()
    selected_idx = st.selectbox(
        "Gene detail", options, format_func=lambda i: _gene_label(filtered.loc[i]),
        key=f"gene_detail_{res.name}"
    )
    row = filtered.loc[selected_idx]
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Status", str(row["de_status"]))
    m2.metric("log2FC", f"{row['log2FoldChange']:.3g}")
    m3.metric("padj", f"{row['padj']:.3g}" if pd.notna(row["padj"]) else "NA")
    if "baseMean" in row.index and pd.notna(row["baseMean"]):
        m4.metric("baseMean", f"{row['baseMean']:.3g}")
    else:
        m4.metric("baseMean", "NA")
    with st.expander("Full selected-gene row", expanded=False):
        detail = row.astype(str).to_frame("value")
        _dataframe(detail, key=f"gene_detail_row_{res.name}")
    if st.button("Highlight this gene in volcano", key=f"gene_highlight_{res.name}"):
        st.session_state[f"highlight_genes_{res.name}"] = [_gene_label(row).split(" / ")[0]]


def tab_de_plots(res: ContrastResult, p: dict):
    st.subheader(f"Differential expression plots — {res.name}")
    label_cols = _gene_columns(res.df) or ["gene_id"]
    controls = st.columns([1, 1, 1, 1])
    with controls[0]:
        label_col = st.selectbox("Volcano label", label_cols,
                                 key=f"volcano_label_{res.name}")
    with controls[1]:
        label_top_n = st.slider("Top labels", 0, 40, 12,
                                key=f"volcano_label_n_{res.name}")
    with controls[2]:
        preset = st.selectbox(
            "Highlight preset", ["Custom", "Top up", "Top down", "Most significant"],
            key=f"volcano_preset_{res.name}",
        )
    with controls[3]:
        preset_n = st.slider("Highlight N", 1, 30, 8,
                             key=f"volcano_preset_n_{res.name}")

    highlight_key = f"highlight_genes_{res.name}"
    saved_highlights = st.session_state.get(highlight_key, [])
    if preset == "Custom":
        highlight_text = st.text_area(
            "Highlight genes", value="\n".join(saved_highlights), height=80,
            key=f"volcano_highlight_text_{res.name}",
        )
        highlight_genes = _parse_gene_text(highlight_text)
        st.session_state[highlight_key] = highlight_genes
    else:
        highlight_genes = _top_gene_ids(res.df, preset, preset_n)

    left, right = st.columns([1.4, 1])
    with left:
        _plotly(
            plots.volcano(
                res.df, p["padj"], p["lfc"], label_top_n=label_top_n,
                label_col=label_col, highlight_genes=highlight_genes,
            ),
            key=f"de_volcano_{res.name}",
        )
    with right:
        _plotly(plots.de_count_bar(res.df, p["padj"], p["lfc"]),
                key=f"de_count_{res.name}")
        _plotly(plots.ma_plot(res.df, p["padj"], p["lfc"]),
                key=f"de_ma_{res.name}")

    d1, d2 = st.columns(2)
    with d1:
        _plotly(plots.de_histogram(res.df, "log2FoldChange", p["padj"], p["lfc"]),
                key=f"de_hist_lfc_{res.name}")
        _plotly(plots.de_histogram(res.df, "padj", p["padj"], p["lfc"]),
                key=f"de_hist_padj_{res.name}")
    with d2:
        if "baseMean" in res.df.columns:
            _plotly(plots.de_histogram(res.df, "baseMean", p["padj"], p["lfc"]),
                    key=f"de_hist_basemean_{res.name}")
        _plotly(plots.pvalue_adjustment_scatter(res.df),
                key=f"de_pvalue_scatter_{res.name}")

    _plotly(plots.threshold_sensitivity(res.df),
            key=f"de_threshold_sensitivity_{res.name}")


def tab_ora(res: ContrastResult):
    st.subheader(f"Over-representation (ORA) — {res.name}")
    if "ora" in res.errors:
        st.error(f"ORA failed: {res.errors['ora']}")
        return
    if res.ora is None or res.ora.empty:
        st.info("No ORA results. Run ORA with at least one source selected.")
        return
    ora = res.ora.copy()
    controls = st.columns([1, 1, 1, 1])
    with controls[0]:
        sources = sorted(ora["source"].dropna().astype(str).unique()) if "source" in ora else []
        source_sel = st.multiselect("Source", sources, default=sources,
                                    key=f"ora_source_filter_{res.name}")
    with controls[1]:
        directions = sorted(ora["direction"].dropna().astype(str).unique()) if "direction" in ora else []
        direction_sel = st.multiselect("Direction", directions, default=directions,
                                       key=f"ora_dir_filter_{res.name}")
    with controls[2]:
        max_p = st.number_input("Max ORA p", min_value=0.0, max_value=1.0,
                                value=1.0, step=0.01, format="%.3f",
                                key=f"ora_p_filter_{res.name}")
    with controls[3]:
        term_query = st.text_input("Search terms", key=f"ora_search_{res.name}")

    if source_sel and "source" in ora.columns:
        ora = ora[ora["source"].astype(str).isin(source_sel)]
    if direction_sel and "direction" in ora.columns:
        ora = ora[ora["direction"].astype(str).isin(direction_sel)]
    if "p_value" in ora.columns:
        ora = ora[pd.to_numeric(ora["p_value"], errors="coerce") <= max_p]
    if term_query:
        q = term_query.strip().lower()
        mask = pd.Series(False, index=ora.index)
        for col in ("term_name", "term_id", "source"):
            if col in ora.columns:
                mask |= ora[col].astype(str).str.lower().str.contains(q, na=False)
        ora = ora[mask]

    if ora.empty:
        st.info("No ORA terms match the current filters.")
        return

    top_n = st.slider("Top terms", 5, 60, 20, key=f"ora_n_{res.name}")
    xmetric = st.radio("Dotplot X axis", ["gene_ratio", "recall", "neg_log10_p"],
                       horizontal=True, key=f"ora_x_{res.name}")

    v_dot, v_bar, v_bubble, v_heat, v_net, v_table = st.tabs(
        ["Dotplot", "Bar", "Source Bubble", "Gene Heatmap", "Gene Network", "Table"]
    )
    with v_dot:
        _plotly(plots.ora_dotplot(ora, top_n, xmetric),
                key=f"ora_dot_{res.name}")
    with v_bar:
        _plotly(plots.ora_barplot(ora, top_n),
                key=f"ora_bar_{res.name}")
    with v_bubble:
        _plotly(plots.ora_source_bubble(ora, top_n),
                key=f"ora_bubble_{res.name}")
    with v_heat:
        _plotly(plots.ora_gene_heatmap(ora, min(top_n, 25)),
                key=f"ora_heatmap_{res.name}")
    with v_net:
        _plotly(plots.ora_gene_network(ora, min(top_n, 18)),
                key=f"ora_network_{res.name}")
    with v_table:
        _dataframe(ora, key=f"ora_table_{res.name}", height=320)

    term_options = ora.sort_values("p_value").index.tolist() if "p_value" in ora.columns else ora.index.tolist()
    selected_term = st.selectbox(
        "Term detail", term_options, format_func=lambda i: _ora_label(ora.loc[i]),
        key=f"ora_term_detail_{res.name}",
    )
    term_row = ora.loc[selected_term]
    t1, t2, t3, t4 = st.columns(4)
    t1.metric("Source", str(term_row.get("source", "ORA")))
    t2.metric("Direction", str(term_row.get("direction", "all")))
    t3.metric("p-value", f"{term_row.get('p_value', float('nan')):.3g}")
    t4.metric("Genes", str(term_row.get("intersection_size", "NA")))
    genes = _parse_ora_genes(term_row.get("genes", []))
    if genes:
        _dataframe(pd.DataFrame({"genes": genes}),
                   key=f"ora_term_genes_{res.name}", height=180)
        st.download_button(
            "⬇ Term genes (CSV)",
            pd.DataFrame({"genes": genes}).to_csv(index=False),
            file_name=f"{res.name}_{term_row.get('term_id', 'ORA_term')}_genes.csv",
            mime="text/csv",
        )
    st.download_button("⬇ Filtered ORA table (CSV)", ora.to_csv(index=False),
                       file_name=f"{res.name}_ORA_filtered.csv", mime="text/csv")


def _render_ortholog_report(metadata: dict) -> None:
    """Surface how much of the chicken list actually reached human GSEA.

    Pre-ranked GSEA treats the ranked list as its universe, so the mapping
    rate is part of interpreting the result -- it belongs next to the table,
    not buried in a log.
    """
    report = metadata.get("ortholog_report")
    if not report:
        return
    rate = metadata.get("mapping_rate")
    if rate is not None and rate < config.ORTHO_MIN_MAPPING_RATE:
        st.warning(report, icon="⚠️")
    else:
        st.caption(f"🧬 {report}")


def tab_gsea(res: ContrastResult):
    st.subheader(f"Gene-set enrichment (GSEA) — {res.name}")
    params = st.session_state.get("params", {})
    _render_gsea_controls(res, params)
    st.divider()

    if params.get("permutations", 0) < 1000:
        st.warning(
            f"GSEA was run with {params.get('permutations')} permutations "
            "(Quick mode). Q-values are noisy. Rerun with full 1000 "
            "permutations before reporting.",
            icon="⚠️",
        )
    if "gsea" in res.errors:
        st.error(f"GSEA failed: {res.errors['gsea']}")
        # A sparse-mapping abort is the most common failure; show the report
        # that explains it rather than making the user read the traceback.
        _render_ortholog_report(res.gsea_metadata or {})
        tb = res.errors.get("gsea_traceback")
        if tb:
            with st.expander("Show error details", expanded=False):
                st.code(tb, language="text")
        return
    if res.gsea is None or res.gsea.table.empty:
        st.info("No GSEA results yet. Select collections above and click Update GSEA.")
        return
    table = res.gsea.table

    with st.expander("Reading this table", expanded=False):
        st.markdown(
            """
            `NES` is the normalized enrichment score. Positive values mean the
            gene set is concentrated near the positive end of the ranked gene
            list; negative values mean it is concentrated near the negative end.
            `FDR q-value` is the multiple-testing adjusted GSEA significance
            estimate. Leading-edge genes are the genes that drive the selected
            enrichment signal.
            """
        )

    metadata = res.gsea_metadata or {}
    fdr_default = 0.25
    controls = st.columns([1, 1, 1.2, 1.4])
    with controls[0]:
        fdr_cutoff = st.slider(
            "FDR cutoff", 0.0, 1.0, fdr_default, 0.01,
            key=f"gsea_fdr_{res.name}",
        )
    with controls[1]:
        direction_choice = st.radio(
            "Direction", ["Both", "Positive NES", "Negative NES"],
            horizontal=False, key=f"gsea_direction_{res.name}",
        )
    with controls[2]:
        collections = sorted(table["collection"].dropna().astype(str).unique())
        collection_sel = st.multiselect(
            "Collections", collections, default=collections,
            key=f"gsea_collection_filter_{res.name}",
        )
    with controls[3]:
        term_query = st.text_input("Search terms", key=f"gsea_search_{res.name}")

    filtered = table.copy()
    if direction_choice == "Positive NES":
        filtered = filtered[filtered["NES"] >= 0]
    elif direction_choice == "Negative NES":
        filtered = filtered[filtered["NES"] < 0]
    if collection_sel:
        filtered = filtered[filtered["collection"].astype(str).isin(collection_sel)]
    if term_query:
        q = term_query.strip().lower()
        term_cols = [c for c in ("term", "term_short", "collection") if c in filtered]
        mask = pd.Series(False, index=filtered.index)
        for col in term_cols:
            mask |= filtered[col].astype(str).str.lower().str.contains(q, na=False)
        filtered = filtered[mask]

    sig = table[pd.to_numeric(table["fdr"], errors="coerce") <= fdr_cutoff]
    sig_pos = int((sig["NES"] >= 0).sum())
    sig_neg = int((sig["NES"] < 0).sum())
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Ranked genes used", f"{metadata.get('ranking_size', len(res.gsea.ranking)):,}")
    requested = int(metadata.get("gene_sets_requested", 0) or 0)
    dropped = max(requested - len(table), 0)
    m2.metric(
        "Gene sets scored", f"{len(table):,}",
        delta=f"-{dropped:,} dropped" if dropped else None,
        delta_color="off",
        help=(
            f"{requested:,} sets requested; {dropped:,} fell outside the "
            f"size bounds ({metadata.get('min_size')}–{metadata.get('max_size')} "
            "genes) after intersecting with the ranked list and were not tested."
        ) if requested else None,
    )
    m3.metric(f"FDR <= {fdr_cutoff:.2f}", f"{len(sig):,}")
    m4.metric("Significant direction", f"{sig_pos} pos / {sig_neg} neg")

    routes = metadata.get("gsea_routes") or []
    if "native" in routes:
        st.caption(
            "🐔 Native chicken GSEA on "
            f"{', '.join(metadata.get('native_sources', []))} — no ortholog "
            f"step, {metadata.get('native_ranking_size', 0):,} chicken genes ranked."
        )
    if metadata.get("fdr_per_route"):
        st.info(
            "This run used two independent preranks (native chicken and human "
            "ortholog). **FDR is computed within each route, not pooled** — "
            "q-values are comparable inside a route, not across the two. The "
            "`gsea_route` column marks which is which.",
            icon="ℹ️",
        )
    _render_ortholog_report(metadata)

    selected_labels = [
        config.GSEA_LIBRARIES.get(lib, lib)
        for lib in metadata.get("libraries", [])
    ]
    if metadata.get("custom_gmt_terms"):
        selected_labels.append(f"Custom GMT ({metadata['custom_gmt_terms']} terms)")
    if selected_labels:
        st.caption("Computed with: " + ", ".join(selected_labels))

    if filtered.empty:
        st.info("No GSEA terms match the current filters.")
        _dataframe(table, key=f"gsea_table_{res.name}", height=320)
        return

    top_n = st.slider("Top terms by |NES|", 5, 40, 20, key=f"gsea_n_{res.name}")
    _plotly(plots.gsea_bar(filtered, top_n, fdr_alpha=fdr_cutoff),
            key=f"gsea_bar_{res.name}")

    st.markdown("**Running-enrichment plot**")
    term_rows = filtered.reset_index(drop=True)
    running_terms = _gsea_terms_with_running_data(
        res.gsea, term_rows["term"].astype(str).tolist()
    )
    if running_terms:
        running_rows = term_rows[term_rows["term"].astype(str).isin(running_terms)]
        term = st.selectbox(
            "Term",
            running_rows["term"].tolist(),
            format_func=lambda t: _gsea_term_label(
                running_rows.loc[running_rows["term"].eq(t)].iloc[0]
            ),
            key=f"gsea_term_{res.name}",
        )
        trimmed = len(term_rows) - len(running_rows)
        if trimmed > 0:
            st.caption(
                f"Running curves are retained for the top "
                f"{config.GSEA_MAX_RUNNING_PLOT_TERMS} terms to limit memory use."
            )
        _plotly(plots.gsea_running_plot(res.gsea, term),
                key=f"gsea_running_{res.name}")
    else:
        term = term_rows.iloc[0]["term"]
        st.info("Running plot data was not retained for the filtered terms.")

    st.markdown("**Leading-edge genes**")
    selected_raw = res.gsea.raw.results.get(term, {})
    leading = [
        g for g in str(selected_raw.get("lead_genes", "")).split(";") if g
    ]
    if leading:
        st.caption(f"Top leading-edge genes for selected term: {', '.join(leading[:20])}")
    default_terms = filtered.head(5)["term"].tolist()
    sel = st.multiselect("Terms for heatmap", filtered["term"].tolist(),
                         default=default_terms, key=f"gsea_le_{res.name}")
    if sel:
        _plotly(plots.leading_edge_heatmap(res.gsea, sel),
                key=f"gsea_leading_edge_{res.name}")

    st.markdown("**Enrichment map**")
    jac = st.slider("Jaccard edge threshold", 0.05, 0.6, 0.25, 0.05,
                    key=f"gsea_jac_{res.name}")
    net_terms = filtered.head(40)["term"].tolist()
    term_map = {}
    for t in net_terms:
        payload = res.gsea.raw.results.get(t, {})
        genes = [g for g in str(payload.get("lead_genes", "")).split(";") if g]
        if genes:
            term_map[t] = genes
    scores = dict(zip(table["term"], table["NES"]))
    _plotly(plots.enrichment_network(term_map, scores, jac),
            key=f"gsea_network_{res.name}")

    _dataframe(filtered, key=f"gsea_table_{res.name}", height=320)
    st.download_button("⬇ Filtered GSEA table (CSV)", filtered.to_csv(index=False),
                       file_name=f"{res.name}_GSEA.csv", mime="text/csv")
    if leading:
        st.download_button(
            "⬇ Leading-edge genes (CSV)",
            pd.DataFrame({"gene": leading}).to_csv(index=False),
            file_name=f"{res.name}_GSEA_leading_edge.csv",
            mime="text/csv",
        )


def tab_compare(results: list[ContrastResult]):
    st.subheader("Compare contrasts")
    if len(results) < 2:
        st.info("Upload ≥2 contrasts to compare shared vs unique signal.")
        return
    what = st.radio("Compare", ["Up-regulated DEGs", "Down-regulated DEGs",
                                "Significant ORA terms"], horizontal=True)
    named = {}
    for r in results:
        if what == "Up-regulated DEGs":
            named[r.name] = r.deg_sets.up
        elif what == "Down-regulated DEGs":
            named[r.name] = r.deg_sets.down
        else:
            if r.ora is not None and not r.ora.empty:
                named[r.name] = r.ora["term_id"].tolist()
    named = {k: v for k, v in named.items() if v}
    if len(named) < 2:
        st.info("Not enough non-empty sets to compare.")
        return
    st.pyplot(plots.upset_comparison(named))


def tab_export(results: list[ContrastResult], p: dict):
    st.subheader("Export bundle")
    st.caption("Zip of result tables, interactive HTML figures, and a run manifest.")
    if st.button("📦 Build export bundle"):
        buf = _io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            manifest = {
                "generated": datetime.utcnow().isoformat() + "Z",
                "parameters": {k: v for k, v in p.items()
                               if k not in ("uploads", "custom_gmt_file", "run")},
                "contrasts": [r.name for r in results],
                # Which route produced these numbers, how many sets were
                # actually scored, and (if orthologs were used at all) the
                # mapping rate. Without this the exported tables cannot be
                # interpreted or reproduced.
                "gsea_metadata": {r.name: (r.gsea_metadata or {}) for r in results},
            }
            z.writestr("run_manifest.json", json.dumps(manifest, indent=2, default=str))
            for r in results:
                if r.ora is not None and not r.ora.empty:
                    z.writestr(f"{r.name}/ORA.csv", r.ora.to_csv(index=False))
                    z.writestr(f"{r.name}/ORA_dotplot.html",
                               plots.ora_dotplot(r.ora).to_html(include_plotlyjs="cdn"))
                if r.gsea is not None and not r.gsea.table.empty:
                    z.writestr(f"{r.name}/GSEA.csv", r.gsea.table.to_csv(index=False))
                    z.writestr(f"{r.name}/GSEA_bar.html",
                               plots.gsea_bar(r.gsea.table).to_html(include_plotlyjs="cdn"))
                z.writestr(f"{r.name}/volcano.html",
                           plots.volcano(r.df, p["padj"], p["lfc"]).to_html(include_plotlyjs="cdn"))
        st.download_button("⬇ Download bundle (.zip)", buf.getvalue(),
                           file_name="enrichment_bundle.zip", mime="application/zip")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    p = sidebar()
    if p["run"]:
        st.session_state.pop("results", None)
        st.session_state.pop("params", None)
        gc.collect()
        with st.spinner("Running enrichment…"):
            st.session_state["results"] = execute(p)
        st.session_state["params"] = _stored_params(p)
        gc.collect()

    results = st.session_state.get("results")
    if not results:
        st.title("DESeq2 → ORA / GSEA — native chicken")
        st.markdown(
            "Upload one or more **DESeq2 result tables** (or use the bundled sample), "
            "set thresholds in the sidebar, and click **Run enrichment**. Both "
            "enrichment steps run on *Gallus gallus* annotations by default — no "
            "human ortholog mapping is involved unless you ask for it.\n\n"
            "* **ORA** runs on native chicken annotations via g:Profiler "
            "(GO / KEGG / Reactome / WikiPathways) with your tested genes as background.\n"
            "* **GSEA** ranks by the DESeq2 Wald statistic and scores against "
            "**native chicken gene sets** (GO / Reactome / WikiPathways). No genes "
            "are dropped for lacking a human counterpart, so the ranked list stays "
            "your real gene universe.\n"
            "* **Hallmark** and **Oncogenic Signatures** are human-only collections "
            "with no chicken equivalent. Select either in the GSEA tab and it runs as "
            "a separate ortholog-mapped prerank, with its own FDR and a visible "
            "mapping rate. KEGG has no downloadable chicken gene sets and is "
            "available in the ORA tab only.\n"
            "* Add a **custom .gmt** in the GSEA tab to score your own curated pathway modules."
        )
        return

    params = st.session_state.get("params", _stored_params(p))
    names = [r.name for r in results]
    chosen = st.selectbox("Contrast", names) if len(names) > 1 else names[0]
    res = next(r for r in results if r.name == chosen)

    if res.errors.get("load"):
        st.error(f"Could not load {res.name}: missing columns {res.errors['load']}")
        return

    t_qc, t_genes, t_de, t_ora, t_gsea, t_cmp, t_exp = st.tabs(
        ["Upload & QC", "Genes", "DE Plots", "ORA", "GSEA", "Compare", "Export"]
    )
    with t_qc:
        tab_qc(res, params)
    with t_genes:
        tab_genes(res, params)
    with t_de:
        tab_de_plots(res, params)
    with t_ora:
        tab_ora(res)
    with t_gsea:
        tab_gsea(res)
    with t_cmp:
        tab_compare(results)
    with t_exp:
        tab_export(results, params)


if __name__ == "__main__":
    main()
