"""All figures for the pipeline and the app.

Interactive figures (volcano, MA, dotplots, GSEA curve, heatmap, network) are
Plotly ``go.Figure`` objects so they work natively in Streamlit and export to
HTML/PNG/SVG/PDF via kaleido. The UpSet comparison is matplotlib (via
``upsetplot``) and returns a ``matplotlib.figure.Figure``.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from . import config


# --------------------------------------------------------------------------
# DE-level plots
# --------------------------------------------------------------------------
def de_plot_table(
    df: pd.DataFrame,
    padj_threshold: float = config.PADJ_THRESHOLD,
    lfc_threshold: float = config.LFC_THRESHOLD,
) -> pd.DataFrame:
    """Return a plotting copy with common DE display columns."""
    d = df.copy()
    d["neg_log10_padj"] = np.nan
    mask = d["padj"].notna()
    d.loc[mask, "neg_log10_padj"] = -np.log10(
        d.loc[mask, "padj"].clip(lower=config.PVALUE_FLOOR)
    )
    sig = (
        d["padj"].notna()
        & (d["padj"] < padj_threshold)
        & (d["log2FoldChange"].abs() > lfc_threshold)
    )
    d["de_status"] = np.where(
        sig & (d["log2FoldChange"] > 0), "Up",
        np.where(sig & (d["log2FoldChange"] < 0), "Down", "n.s."),
    )
    return d


def _highlight_mask(df: pd.DataFrame, genes: list[str] | None) -> pd.Series:
    if not genes:
        return pd.Series(False, index=df.index)
    wanted = {str(g).strip().lower() for g in genes if str(g).strip()}
    mask = pd.Series(False, index=df.index)
    for col in ("gene_id", "gene_name", "entrez_id"):
        if col in df.columns:
            mask |= df[col].astype(str).str.lower().isin(wanted)
    return mask


def volcano(
    df: pd.DataFrame,
    padj_threshold: float = config.PADJ_THRESHOLD,
    lfc_threshold: float = config.LFC_THRESHOLD,
    label_top_n: int = 12,
    label_col: str = "gene_name",
    highlight_genes: list[str] | None = None,
) -> go.Figure:
    """Volcano plot coloured by significance/direction with top-gene labels."""
    d = de_plot_table(df, padj_threshold, lfc_threshold)
    d = d[d["padj"].notna()]
    color_map = {"Up": config.COLOR_UP, "Down": config.COLOR_DOWN, "n.s.": config.COLOR_NS}
    hover_col = label_col if label_col in d.columns else "gene_id"

    fig = go.Figure()
    for cls in ("n.s.", "Down", "Up"):
        sub = d[d["de_status"] == cls]
        fig.add_trace(go.Scattergl(
            x=sub["log2FoldChange"], y=sub["neg_log10_padj"],
            mode="markers", name=cls,
            marker=dict(size=5, color=color_map[cls], opacity=0.6 if cls == "n.s." else 0.85),
            text=sub[hover_col] if hover_col in sub.columns else None,
            hovertemplate="%{text}<br>log2FC=%{x:.2f}<br>-log10 padj=%{y:.2f}<extra></extra>",
        ))

    highlight = d[_highlight_mask(d, highlight_genes)]
    if not highlight.empty:
        fig.add_trace(go.Scatter(
            x=highlight["log2FoldChange"], y=highlight["neg_log10_padj"],
            mode="markers+text", name="Highlighted",
            marker=dict(size=12, color="#F1C40F", symbol="star",
                        line=dict(width=1.2, color="#111")),
            text=highlight[hover_col] if hover_col in highlight.columns else None,
            textposition="top center",
            hovertemplate="%{text}<br>log2FC=%{x:.2f}<br>-log10 padj=%{y:.2f}<extra></extra>",
        ))

    # Threshold guides.
    fig.add_vline(x=lfc_threshold, line=dict(dash="dot", color="#888", width=1))
    fig.add_vline(x=-lfc_threshold, line=dict(dash="dot", color="#888", width=1))
    fig.add_hline(y=-np.log10(padj_threshold), line=dict(dash="dot", color="#888", width=1))

    # Label the most extreme significant genes (by combined rank).
    if label_col in d.columns and label_top_n > 0:
        sig_d = d[d["de_status"] != "n.s."].copy()
        sig_d["score"] = sig_d["neg_log10_padj"] * sig_d["log2FoldChange"].abs()
        for _, row in sig_d.nlargest(label_top_n, "score").iterrows():
            fig.add_annotation(
                x=row["log2FoldChange"], y=row["neg_log10_padj"],
                text=str(row[label_col]), showarrow=False,
                font=dict(size=10), yshift=8,
            )

    fig.update_layout(
        template="simple_white",
        xaxis_title="log2 fold change",
        yaxis_title="-log10 adjusted p",
        legend_title="",
        height=520,
        title=f"Volcano  (padj<{padj_threshold}, |log2FC|>{lfc_threshold})",
    )
    return fig


def ma_plot(df: pd.DataFrame,
            padj_threshold: float = config.PADJ_THRESHOLD,
            lfc_threshold: float = config.LFC_THRESHOLD) -> go.Figure:
    """MA plot: log2FC vs mean expression."""
    d = df.copy()
    d = d[d["baseMean"] > 0]
    d["log10_mean"] = np.log10(d["baseMean"])
    sig = (d["padj"] < padj_threshold) & (d["log2FoldChange"].abs() > lfc_threshold)
    d["cls"] = np.where(sig & (d["log2FoldChange"] > 0), "Up",
                        np.where(sig & (d["log2FoldChange"] < 0), "Down", "n.s."))
    cmap = {"Up": config.COLOR_UP, "Down": config.COLOR_DOWN, "n.s.": config.COLOR_NS}
    fig = go.Figure()
    for cls in ("n.s.", "Down", "Up"):
        sub = d[d["cls"] == cls]
        fig.add_trace(go.Scattergl(
            x=sub["log10_mean"], y=sub["log2FoldChange"], mode="markers", name=cls,
            marker=dict(size=4, color=cmap[cls], opacity=0.55 if cls == "n.s." else 0.85),
        ))
    fig.add_hline(y=0, line=dict(color="#555", width=1))
    fig.update_layout(template="simple_white", xaxis_title="log10 mean expression",
                      yaxis_title="log2 fold change", height=460, title="MA plot")
    return fig


def de_count_bar(
    df: pd.DataFrame,
    padj_threshold: float = config.PADJ_THRESHOLD,
    lfc_threshold: float = config.LFC_THRESHOLD,
) -> go.Figure:
    d = de_plot_table(df, padj_threshold, lfc_threshold)
    counts = d["de_status"].value_counts().reindex(["Up", "Down", "n.s."], fill_value=0)
    fig = go.Figure(go.Bar(
        x=counts.index,
        y=counts.values,
        marker=dict(color=[config.COLOR_UP, config.COLOR_DOWN, config.COLOR_NS]),
        text=counts.values,
        textposition="outside",
        hovertemplate="%{x}<br>genes=%{y}<extra></extra>",
    ))
    fig.update_layout(template="simple_white", height=320,
                      title="DE gene counts", xaxis_title="", yaxis_title="Genes")
    return fig


def de_histogram(
    df: pd.DataFrame,
    column: str,
    padj_threshold: float = config.PADJ_THRESHOLD,
    lfc_threshold: float = config.LFC_THRESHOLD,
    bins: int = 60,
) -> go.Figure:
    d = de_plot_table(df, padj_threshold, lfc_threshold)
    d = d[pd.to_numeric(d[column], errors="coerce").notna()]
    title_map = {
        "padj": "Adjusted p-value distribution",
        "pvalue": "Raw p-value distribution",
        "log2FoldChange": "log2 fold-change distribution",
        "baseMean": "Mean expression distribution",
    }
    fig = go.Figure()
    for cls in ("n.s.", "Down", "Up"):
        sub = d[d["de_status"] == cls]
        x = pd.to_numeric(sub[column], errors="coerce")
        if column == "baseMean":
            x = np.log10(x[x > 0])
        fig.add_trace(go.Histogram(
            x=x, nbinsx=bins, name=cls,
            marker_color={"Up": config.COLOR_UP, "Down": config.COLOR_DOWN,
                          "n.s.": config.COLOR_NS}[cls],
            opacity=0.72,
        ))
    x_title = "log10 baseMean" if column == "baseMean" else column
    fig.update_layout(template="simple_white", barmode="overlay", height=360,
                      title=title_map.get(column, f"{column} distribution"),
                      xaxis_title=x_title, yaxis_title="Genes")
    return fig


def pvalue_adjustment_scatter(df: pd.DataFrame) -> go.Figure:
    d = df.copy()
    d = d[d["pvalue"].notna() & d["padj"].notna()]
    fig = go.Figure(go.Scattergl(
        x=d["pvalue"], y=d["padj"], mode="markers",
        marker=dict(size=4, color=config.COLOR_ACCENT, opacity=0.45),
        text=d["gene_name"] if "gene_name" in d.columns else d.get("gene_id"),
        hovertemplate="%{text}<br>p=%{x:.2e}<br>padj=%{y:.2e}<extra></extra>",
    ))
    fig.add_trace(go.Scatter(
        x=[0, 1], y=[0, 1], mode="lines", showlegend=False,
        line=dict(color="#888", dash="dot", width=1),
        hoverinfo="skip",
    ))
    fig.update_layout(template="simple_white", height=420,
                      title="Raw p-value vs adjusted p-value",
                      xaxis_title="pvalue", yaxis_title="padj")
    return fig


def threshold_sensitivity(
    df: pd.DataFrame,
    padj_values: list[float] | None = None,
    lfc_values: list[float] | None = None,
) -> go.Figure:
    d = df[df["padj"].notna()].copy()
    padj_values = padj_values or [0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2]
    lfc_values = lfc_values or [0, 0.5, 1.0, 1.5, 2.0]
    z = []
    text = []
    for padj in padj_values:
        row = []
        text_row = []
        for lfc in lfc_values:
            up = ((d["padj"] < padj)
                  & (d["log2FoldChange"] > lfc)).sum()
            down = ((d["padj"] < padj)
                    & (d["log2FoldChange"] < -lfc)).sum()
            total = int(up + down)
            row.append(total)
            text_row.append(f"total={total}<br>up={int(up)}<br>down={int(down)}")
        z.append(row)
        text.append(text_row)

    fig = go.Figure(go.Heatmap(
        z=z,
        x=[str(v) for v in lfc_values],
        y=[str(v) for v in padj_values],
        text=text,
        colorscale="Viridis",
        colorbar=dict(title="DEGs"),
        hovertemplate="padj < %{y}<br>|log2FC| > %{x}<br>%{text}<extra></extra>",
    ))
    fig.update_layout(template="simple_white", height=420,
                      title="Threshold sensitivity",
                      xaxis_title="|log2 fold change| threshold",
                      yaxis_title="padj threshold")
    return fig


# --------------------------------------------------------------------------
# ORA plots
# --------------------------------------------------------------------------
def _flatten_column_index(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with simple string column names."""
    d = df.copy()
    if isinstance(d.columns, pd.MultiIndex):
        flat = []
        for col in d.columns:
            parts = [str(part) for part in col
                     if pd.notna(part) and str(part) not in ("", "nan")]
            flat.append(parts[0] if parts else "")
        d.columns = flat
    else:
        d.columns = [str(col) for col in d.columns]
    return d


def _first_column(df: pd.DataFrame, candidates: list[str],
                  default=None) -> pd.Series | None:
    for col in candidates:
        if col not in df.columns:
            continue
        try:
            value = df.loc[:, col]
        except KeyError:
            continue
        if isinstance(value, pd.DataFrame):
            if value.shape[1] == 0:
                continue
            value = value.iloc[:, 0]
        return value
    if default is None:
        return None
    return pd.Series(default, index=df.index)


def _normalise_ora_plot_df(ora_df: pd.DataFrame,
                           x: str) -> tuple[pd.DataFrame | None, str | None]:
    d = _flatten_column_index(ora_df)
    term_name = _first_column(d, ["term_name", "name", "term_id", "native"])
    if term_name is None:
        return None, "ORA terms are missing labels; rerun ORA"

    out = d.copy()
    out["term_name"] = term_name.astype(str)
    out["source"] = _first_column(out, ["source"], "ORA").fillna("ORA").astype(str)
    out["direction"] = (
        _first_column(out, ["direction"], "all").fillna("all").astype(str)
    )

    p_value = _first_column(out, ["p_value"])
    if "neg_log10_p" not in out.columns and p_value is not None:
        p = pd.to_numeric(p_value, errors="coerce")
        out["neg_log10_p"] = -np.log10(p.clip(lower=config.PVALUE_FLOOR))
    if "gene_ratio" not in out.columns:
        inter = _first_column(out, ["intersection_size"])
        query = _first_column(out, ["query_size"])
        if inter is not None and query is not None:
            out["gene_ratio"] = (
                pd.to_numeric(inter, errors="coerce")
                / pd.to_numeric(query, errors="coerce").replace(0, np.nan)
            )
    if "recall" not in out.columns:
        inter = _first_column(out, ["intersection_size"])
        term_size = _first_column(out, ["term_size"])
        if inter is not None and term_size is not None:
            out["recall"] = (
                pd.to_numeric(inter, errors="coerce")
                / pd.to_numeric(term_size, errors="coerce").replace(0, np.nan)
            )

    required = [x, "p_value", "neg_log10_p", "intersection_size", "term_size"]
    missing = [col for col in required if col not in out.columns]
    if missing:
        return None, "ORA table is missing: " + ", ".join(missing)

    for col in required:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out.dropna(subset=[x, "p_value", "neg_log10_p", "intersection_size"])
    if out.empty:
        return None, "No plottable ORA terms"
    return out, None


def ora_dotplot(ora_df: pd.DataFrame, top_n: int = 15,
                x: str = "gene_ratio") -> go.Figure:
    """Dotplot of top ORA terms per direction.

    Dot size = intersection size, colour = -log10 p. Terms are grouped by
    direction (up above, down below) for readability.
    """
    if ora_df is None or ora_df.empty:
        return _empty_fig("No ORA terms")
    ora_df, error = _normalise_ora_plot_df(ora_df, x)
    if error:
        return _empty_fig(error)

    frames = []
    for direction, grp in ora_df.groupby("direction"):
        frames.append(grp.nsmallest(top_n, "p_value"))
    d = pd.concat(frames, ignore_index=True) if frames else ora_df
    d = d.sort_values(["direction", x])
    label = (d["term_name"].astype(str).str.slice(0, 45)
             + " (" + d["source"].astype(str) + ")")
    size_max = d["intersection_size"].max()
    if pd.isna(size_max) or size_max <= 0:
        size_max = 1

    fig = go.Figure(go.Scatter(
        x=d[x], y=label, mode="markers",
        marker=dict(
            size=8 + 22 * (d["intersection_size"] / size_max),
            color=d["neg_log10_p"], colorscale="Viridis",
            colorbar=dict(title="-log10 p"), line=dict(width=0.5, color="#333"),
        ),
        text=d.apply(lambda r: f"{r['term_name']}<br>{r['source']} · {r['direction']}"
                                f"<br>genes={r['intersection_size']}/{r['term_size']}"
                                f"<br>p={r['p_value']:.2e}", axis=1),
        hovertemplate="%{text}<extra></extra>",
    ))
    fig.update_layout(template="simple_white",
                      xaxis_title=x.replace("_", " "),
                      height=max(360, 24 * len(d)),
                      title=f"ORA — top {top_n} terms per direction",
                      margin=dict(l=10, r=10))
    return fig


def _top_ora_terms(ora_df: pd.DataFrame, top_n: int,
                   direction: str | None = None,
                   source: str | None = None) -> pd.DataFrame:
    d, error = _normalise_ora_plot_df(ora_df, "neg_log10_p")
    if error or d is None:
        return pd.DataFrame()
    if direction and direction != "all directions":
        d = d[d["direction"] == direction]
    if source and source != "all sources":
        d = d[d["source"] == source]
    if d.empty:
        return d
    return d.nsmallest(top_n, "p_value")


def ora_barplot(ora_df: pd.DataFrame, top_n: int = 20,
                direction: str | None = None,
                source: str | None = None) -> go.Figure:
    d = _top_ora_terms(ora_df, top_n, direction, source)
    if d.empty:
        return _empty_fig("No ORA terms for this filter")
    d = d.sort_values("neg_log10_p")
    label = d["term_name"].astype(str).str.slice(0, 58)
    colors = np.where(d["direction"].eq("up"), config.COLOR_UP,
                      np.where(d["direction"].eq("down"), config.COLOR_DOWN,
                               config.COLOR_ACCENT))
    fig = go.Figure(go.Bar(
        x=d["neg_log10_p"], y=label, orientation="h",
        marker=dict(color=colors),
        text=d["source"],
        hovertemplate="%{y}<br>-log10 p=%{x:.2f}<br>%{text}<extra></extra>",
    ))
    fig.update_layout(template="simple_white", height=max(360, 24 * len(d)),
                      title=f"ORA top {top_n} terms",
                      xaxis_title="-log10 p", yaxis_title="")
    return fig


def ora_source_bubble(ora_df: pd.DataFrame, top_n: int = 40) -> go.Figure:
    d = _top_ora_terms(ora_df, top_n)
    if d.empty:
        return _empty_fig("No ORA terms")
    size_max = d["intersection_size"].max()
    if pd.isna(size_max) or size_max <= 0:
        size_max = 1
    fig = go.Figure(go.Scatter(
        x=d["gene_ratio"], y=d["source"], mode="markers",
        marker=dict(
            size=8 + 26 * (d["intersection_size"] / size_max),
            color=d["neg_log10_p"], colorscale="Viridis",
            colorbar=dict(title="-log10 p"),
            line=dict(width=0.5, color="#333"),
        ),
        text=d.apply(lambda r: f"{r['term_name']}<br>{r['direction']}"
                                f"<br>genes={r['intersection_size']}/{r['term_size']}"
                                f"<br>p={r['p_value']:.2e}", axis=1),
        hovertemplate="%{text}<extra></extra>",
    ))
    fig.update_layout(template="simple_white", height=420,
                      title=f"ORA source bubble plot (top {top_n})",
                      xaxis_title="gene ratio", yaxis_title="source")
    return fig


def _parse_gene_list(value) -> list[str]:
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


def _ora_term_gene_frame(ora_df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    d = _top_ora_terms(ora_df, top_n)
    if d.empty or "genes" not in d.columns:
        return pd.DataFrame()
    rows = []
    for _, row in d.iterrows():
        for gene in _parse_gene_list(row["genes"]):
            rows.append({
                "term": str(row["term_name"])[:60],
                "gene": gene,
                "source": row["source"],
                "direction": row["direction"],
                "neg_log10_p": row["neg_log10_p"],
            })
    return pd.DataFrame(rows)


def ora_gene_heatmap(ora_df: pd.DataFrame, top_n: int = 15,
                     max_genes: int = 45) -> go.Figure:
    edges = _ora_term_gene_frame(ora_df, top_n)
    if edges.empty:
        return _empty_fig("No ORA gene intersections available")
    genes = edges["gene"].value_counts().head(max_genes).index.tolist()
    terms = edges["term"].drop_duplicates().tolist()
    z = []
    for term in terms:
        term_genes = set(edges.loc[edges["term"] == term, "gene"])
        z.append([1 if gene in term_genes else 0 for gene in genes])
    fig = go.Figure(go.Heatmap(
        z=z, x=genes, y=terms, colorscale=[[0, "#f2f2f2"], [1, config.COLOR_ACCENT]],
        showscale=False,
        hovertemplate="%{y}<br>%{x}<extra></extra>",
    ))
    fig.update_layout(template="simple_white", height=max(360, 24 * len(terms)),
                      title=f"ORA term-gene heatmap (top {top_n})",
                      xaxis=dict(tickangle=45), xaxis_title="Genes",
                      yaxis_title="Terms")
    return fig


def ora_gene_network(ora_df: pd.DataFrame, top_n: int = 12,
                     max_genes: int = 40) -> go.Figure:
    edges = _ora_term_gene_frame(ora_df, top_n)
    if edges.empty:
        return _empty_fig("No ORA gene intersections available")
    keep_genes = set(edges["gene"].value_counts().head(max_genes).index)
    edges = edges[edges["gene"].isin(keep_genes)]
    if edges.empty:
        return _empty_fig("No ORA gene intersections available")

    import networkx as nx
    graph = nx.Graph()
    for _, row in edges.iterrows():
        term = "term:" + row["term"]
        gene = "gene:" + row["gene"]
        graph.add_node(term, label=row["term"], kind="term", score=row["neg_log10_p"])
        graph.add_node(gene, label=row["gene"], kind="gene", score=0)
        graph.add_edge(term, gene)
    pos = nx.spring_layout(graph, seed=config.GSEA_SEED, k=0.75)
    edge_x, edge_y = [], []
    for a, b in graph.edges():
        edge_x += [pos[a][0], pos[b][0], None]
        edge_y += [pos[a][1], pos[b][1], None]
    edge_trace = go.Scatter(x=edge_x, y=edge_y, mode="lines",
                            line=dict(width=0.7, color="#cccccc"),
                            hoverinfo="none")
    term_nodes = [n for n, data in graph.nodes(data=True) if data["kind"] == "term"]
    gene_nodes = [n for n, data in graph.nodes(data=True) if data["kind"] == "gene"]
    traces = [edge_trace]
    for nodes, name, color, size in (
        (term_nodes, "Terms", config.COLOR_ACCENT, 16),
        (gene_nodes, "Genes", "#6C7A89", 9),
    ):
        traces.append(go.Scatter(
            x=[pos[n][0] for n in nodes],
            y=[pos[n][1] for n in nodes],
            mode="markers+text",
            name=name,
            marker=dict(size=size, color=color, line=dict(width=0.5, color="#222")),
            text=[graph.nodes[n]["label"][:24] for n in nodes],
            textposition="top center",
            hovertemplate="%{text}<extra></extra>",
        ))
    fig = go.Figure(traces)
    fig.update_layout(template="simple_white", height=620,
                      title=f"ORA term-gene network (top {top_n})",
                      xaxis=dict(visible=False), yaxis=dict(visible=False))
    return fig


# --------------------------------------------------------------------------
# GSEA plots
# --------------------------------------------------------------------------
def gsea_bar(gsea_df: pd.DataFrame, top_n: int = 20,
             fdr_alpha: float = config.SIG_ALPHA) -> go.Figure:
    """Horizontal NES bars for the top |NES| terms, coloured by direction."""
    if gsea_df is None or gsea_df.empty:
        return _empty_fig("No GSEA terms")
    d = gsea_df.copy()
    d["absnes"] = d["NES"].abs()
    d = d.nlargest(top_n, "absnes").sort_values("NES")
    label = d["term_short"] if "term_short" in d.columns else d["term"]
    colors = np.where(d["NES"] >= 0, config.COLOR_UP, config.COLOR_DOWN)
    opacity = np.where(d["fdr"] < fdr_alpha, 1.0, 0.4)
    fig = go.Figure(go.Bar(
        x=d["NES"], y=label.astype(str).str.slice(0, 50), orientation="h",
        marker=dict(color=colors, opacity=opacity),
        text=[f"FDR={v:.2g}" for v in d["fdr"]],
        hovertemplate="%{y}<br>NES=%{x:.2f}<br>%{text}<extra></extra>",
    ))
    fig.add_vline(x=0, line=dict(color="#555", width=1))
    fig.update_layout(template="simple_white", xaxis_title="NES",
                      height=max(360, 26 * len(d)),
                      title=f"GSEA — top {top_n} by |NES| (faded = FDR≥{fdr_alpha})")
    return fig


def gsea_running_plot(gsea_result, term: str) -> go.Figure:
    """Classic GSEA running-enrichment plot for one term (3 stacked panels)."""
    res = gsea_result.raw.results[term]
    RES = np.asarray(res["RES"], dtype=float)
    hits = np.asarray(res["hits"], dtype=int)
    ranking = gsea_result.ranking.values
    x = np.arange(len(RES))

    from plotly.subplots import make_subplots
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True,
                        row_heights=[0.6, 0.12, 0.28], vertical_spacing=0.02)
    # running ES
    fig.add_trace(go.Scatter(x=x, y=RES, mode="lines",
                             line=dict(color=config.COLOR_ACCENT, width=2),
                             name="Running ES"), row=1, col=1)
    peak = int(np.argmax(np.abs(RES)))
    fig.add_trace(go.Scatter(x=[peak], y=[RES[peak]], mode="markers",
                             marker=dict(color="black", size=7), showlegend=False),
                  row=1, col=1)
    # hit ticks
    fig.add_trace(go.Scatter(x=hits, y=np.zeros_like(hits), mode="markers",
                             marker=dict(symbol="line-ns-open", size=10,
                                         color="#444"), showlegend=False),
                  row=2, col=1)
    # ranking metric
    fig.add_trace(go.Scatter(x=x, y=ranking, mode="lines",
                             line=dict(color="#999", width=1), fill="tozeroy",
                             showlegend=False), row=3, col=1)
    nes = res["nes"]
    fdr = res["fdr"]
    fig.update_layout(
        template="simple_white", height=520,
        title=f"{term}<br><sup>NES={nes:.2f}  FDR={fdr:.2g}</sup>",
    )
    fig.update_yaxes(title_text="Enrichment score", row=1, col=1)
    fig.update_yaxes(showticklabels=False, row=2, col=1)
    fig.update_yaxes(title_text="Ranked metric", row=3, col=1)
    fig.update_xaxes(title_text="Rank in ordered gene list", row=3, col=1)
    return fig


def leading_edge_heatmap(gsea_result, terms: list[str], max_genes: int = 40) -> go.Figure:
    """Heatmap of leading-edge genes (union across terms) × terms.

    Cell value = the gene's ranking metric; shows which genes drive each set.
    """
    ranking = gsea_result.ranking
    gene_union: list[str] = []
    per_term: dict[str, set] = {}
    for term in terms:
        lead = gsea_result.raw.results[term]["lead_genes"].split(";")
        lead = [g for g in lead if g]
        per_term[term] = set(lead)
        for g in lead:
            if g not in gene_union:
                gene_union.append(g)
    # keep the most extreme genes if too many
    gene_union = sorted(gene_union, key=lambda g: -abs(ranking.get(g, 0)))[:max_genes]
    z = []
    for term in terms:
        row = [ranking.get(g, np.nan) if g in per_term[term] else np.nan
               for g in gene_union]
        z.append(row)
    fig = go.Figure(go.Heatmap(
        z=z, x=gene_union, y=[t[:40] for t in terms],
        colorscale=config.CONTINUOUS_SCALE, zmid=0,
        colorbar=dict(title="rank metric"),
        hovertemplate="%{y}<br>%{x}<br>metric=%{z:.2f}<extra></extra>",
    ))
    fig.update_layout(template="simple_white",
                      height=max(300, 40 * len(terms)),
                      title="Leading-edge genes per term",
                      xaxis=dict(tickangle=45))
    return fig


# --------------------------------------------------------------------------
# Network + comparison
# --------------------------------------------------------------------------
def enrichment_network(term_gene_map: dict[str, list[str]],
                       scores: dict[str, float] | None = None,
                       jaccard_min: float = 0.25,
                       max_terms: int = 40) -> go.Figure:
    """Term-term similarity network (edges = Jaccard of member genes).

    ``term_gene_map`` : {term: [genes]}. ``scores`` optionally colours nodes
    (e.g. NES or -log10 p).
    """
    import networkx as nx
    terms = list(term_gene_map.keys())[:max_terms]
    sets = {t: set(term_gene_map[t]) for t in terms}
    G = nx.Graph()
    for t in terms:
        G.add_node(t)
    for a, b in combinations(terms, 2):
        inter = len(sets[a] & sets[b])
        if inter == 0:
            continue
        union = len(sets[a] | sets[b])
        j = inter / union if union else 0
        if j >= jaccard_min:
            G.add_edge(a, b, weight=j)
    if G.number_of_nodes() == 0:
        return _empty_fig("No terms to network")
    pos = nx.spring_layout(G, seed=config.GSEA_SEED, k=0.6)

    edge_x, edge_y = [], []
    for a, b in G.edges():
        edge_x += [pos[a][0], pos[b][0], None]
        edge_y += [pos[a][1], pos[b][1], None]
    edge_trace = go.Scatter(x=edge_x, y=edge_y, mode="lines",
                            line=dict(width=0.8, color="#bbb"), hoverinfo="none")
    node_x = [pos[t][0] for t in G.nodes()]
    node_y = [pos[t][1] for t in G.nodes()]
    node_color = ([scores.get(t, 0) for t in G.nodes()] if scores
                  else [len(sets[t]) for t in G.nodes()])
    node_trace = go.Scatter(
        x=node_x, y=node_y, mode="markers",
        marker=dict(size=[8 + 0.5 * len(sets[t]) for t in G.nodes()],
                    color=node_color, colorscale="RdBu_r" if scores else "Viridis",
                    colorbar=dict(title="NES" if scores else "size"),
                    line=dict(width=0.5, color="#333")),
        text=[t[:45] for t in G.nodes()], hovertemplate="%{text}<extra></extra>",
    )
    fig = go.Figure([edge_trace, node_trace])
    fig.update_layout(template="simple_white", showlegend=False, height=560,
                      title=f"Enrichment map (Jaccard ≥ {jaccard_min})",
                      xaxis=dict(visible=False), yaxis=dict(visible=False))
    return fig


def upset_comparison(named_sets: dict[str, list[str]], max_intersections: int = 25):
    """Self-contained UpSet plot of overlaps across contrasts/directions.

    Hand-rolled with matplotlib (no ``upsetplot`` dependency) so it is robust
    across pandas/matplotlib versions. ``named_sets`` : {label: [ids]}.
    Returns a ``matplotlib.figure.Figure``.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    contents = {k: set(v) for k, v in named_sets.items() if v}
    labels = list(contents.keys())
    if len(contents) < 2:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, "Need ≥2 non-empty sets for UpSet",
                ha="center", va="center")
        ax.axis("off")
        return fig

    # Assign each element the exact combination (membership pattern) it belongs to.
    universe = set().union(*contents.values())
    from collections import Counter
    patterns: Counter = Counter()
    for elem in universe:
        member = tuple(lbl for lbl in labels if elem in contents[lbl])
        if member:
            patterns[member] += 1
    combos = sorted(patterns.items(), key=lambda kv: -kv[1])[:max_intersections]

    n_sets = len(labels)
    n_combos = len(combos)
    fig = plt.figure(figsize=(max(7, 0.55 * n_combos + 2), 2.2 + 0.45 * n_sets))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, max(1, n_sets)], hspace=0.05)
    ax_bar = fig.add_subplot(gs[0])
    ax_mat = fig.add_subplot(gs[1], sharex=ax_bar)

    xs = np.arange(n_combos)
    heights = [c for _, c in combos]
    ax_bar.bar(xs, heights, color=config.COLOR_ACCENT, width=0.6)
    for x, h in zip(xs, heights):
        ax_bar.text(x, h, str(h), ha="center", va="bottom", fontsize=8)
    ax_bar.set_ylabel("intersection size")
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.tick_params(labelbottom=False)

    # Dot matrix.
    for yi, lbl in enumerate(labels):
        y = n_sets - 1 - yi
        for xi, (member, _) in enumerate(combos):
            on = lbl in member
            ax_mat.scatter(xi, y, s=90,
                           color=config.COLOR_ACCENT if on else "#dddddd",
                           zorder=3)
        # connect dots in a combo
    for xi, (member, _) in enumerate(combos):
        ys = [n_sets - 1 - labels.index(m) for m in member]
        if len(ys) > 1:
            ax_mat.plot([xi, xi], [min(ys), max(ys)], color=config.COLOR_ACCENT,
                        lw=2, zorder=2)
    ax_mat.set_yticks(range(n_sets))
    ax_mat.set_yticklabels(labels[::-1])
    ax_mat.set_xticks([])
    ax_mat.set_ylim(-0.5, n_sets - 0.5)
    ax_mat.spines[["top", "right", "bottom"]].set_visible(False)
    fig.suptitle("Shared vs unique across sets", y=0.98)
    fig.tight_layout()
    return fig


def _empty_fig(msg: str) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=msg, showarrow=False, font=dict(size=14))
    fig.update_layout(template="simple_white", height=300,
                      xaxis=dict(visible=False), yaxis=dict(visible=False))
    return fig
