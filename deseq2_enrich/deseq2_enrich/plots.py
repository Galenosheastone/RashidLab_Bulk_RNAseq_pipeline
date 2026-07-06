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
def volcano(
    df: pd.DataFrame,
    padj_threshold: float = config.PADJ_THRESHOLD,
    lfc_threshold: float = config.LFC_THRESHOLD,
    label_top_n: int = 12,
    label_col: str = "gene_name",
) -> go.Figure:
    """Volcano plot coloured by significance/direction with top-gene labels."""
    d = df.copy()
    d = d[d["padj"].notna()]
    d["neg_log10_padj"] = -np.log10(d["padj"].clip(lower=config.PVALUE_FLOOR))
    sig = (d["padj"] < padj_threshold) & (d["log2FoldChange"].abs() > lfc_threshold)
    d["cls"] = np.where(
        sig & (d["log2FoldChange"] > 0), "Up",
        np.where(sig & (d["log2FoldChange"] < 0), "Down", "n.s."),
    )
    color_map = {"Up": config.COLOR_UP, "Down": config.COLOR_DOWN, "n.s.": config.COLOR_NS}

    fig = go.Figure()
    for cls in ("n.s.", "Down", "Up"):
        sub = d[d["cls"] == cls]
        fig.add_trace(go.Scattergl(
            x=sub["log2FoldChange"], y=sub["neg_log10_padj"],
            mode="markers", name=cls,
            marker=dict(size=5, color=color_map[cls], opacity=0.6 if cls == "n.s." else 0.85),
            text=sub[label_col] if label_col in sub.columns else None,
            hovertemplate="%{text}<br>log2FC=%{x:.2f}<br>-log10 padj=%{y:.2f}<extra></extra>",
        ))

    # Threshold guides.
    fig.add_vline(x=lfc_threshold, line=dict(dash="dot", color="#888", width=1))
    fig.add_vline(x=-lfc_threshold, line=dict(dash="dot", color="#888", width=1))
    fig.add_hline(y=-np.log10(padj_threshold), line=dict(dash="dot", color="#888", width=1))

    # Label the most extreme significant genes (by combined rank).
    if label_col in d.columns and label_top_n > 0:
        sig_d = d[d["cls"] != "n.s."].copy()
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


# --------------------------------------------------------------------------
# ORA plots
# --------------------------------------------------------------------------
def ora_dotplot(ora_df: pd.DataFrame, top_n: int = 15,
                x: str = "gene_ratio") -> go.Figure:
    """Dotplot of top ORA terms per direction.

    Dot size = intersection size, colour = -log10 p. Terms are grouped by
    direction (up above, down below) for readability.
    """
    if ora_df is None or ora_df.empty:
        return _empty_fig("No ORA terms")
    ora_df = ora_df.copy()
    if "term_name" not in ora_df.columns:
        if "name" in ora_df.columns:
            ora_df["term_name"] = ora_df["name"]
        elif "term_id" in ora_df.columns:
            ora_df["term_name"] = ora_df["term_id"]
        elif "native" in ora_df.columns:
            ora_df["term_name"] = ora_df["native"]
        else:
            return _empty_fig("ORA terms are missing labels; rerun ORA")
    if "direction" not in ora_df.columns:
        ora_df["direction"] = "all"
    if "source" not in ora_df.columns:
        ora_df["source"] = "ORA"
    missing = [c for c in (x, "p_value", "neg_log10_p", "intersection_size", "term_size")
               if c not in ora_df.columns]
    if missing:
        return _empty_fig("ORA table is missing: " + ", ".join(missing))
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
