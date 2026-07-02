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

import io as _io
import os
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
sys.path.insert(0, PACKAGE_PARENT)

from deseq2_enrich import config, plots, genesets  # noqa: E402
from deseq2_enrich.pipeline import run_contrast, ContrastResult  # noqa: E402

st.set_page_config(page_title="DESeq2 Enrichment", layout="wide",
                   page_icon="🧬")

SAMPLE_PATH = os.path.join(
    PACKAGE_PARENT,
    "sample_data", "DESeq2_Sacral_vs_Sacralized_Caudal_all.tsv",
)


# --------------------------------------------------------------------------
# Sidebar — inputs and parameters
# --------------------------------------------------------------------------
def sidebar() -> dict:
    st.sidebar.title("🧬 DESeq2 Enrichment")
    st.sidebar.caption("ORA (g:Profiler, native chicken) + GSEA (MSigDB via orthologs)")

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

    st.sidebar.subheader("5 · GSEA gene sets (human orthologs)")
    gsea_libs = st.sidebar.multiselect(
        "MSigDB / pathway libraries", list(config.GSEA_LIBRARIES.keys()),
        default=config.GSEA_DEFAULT_LIBRARIES,
        format_func=lambda s: config.GSEA_LIBRARIES[s],
    )
    custom_gmt_file = st.sidebar.file_uploader(
        "Custom .gmt (optional)", type=["gmt"],
        help="Your curated modules (e.g. cGAS-STING, necroptosis, osteoclast).",
    )

    st.sidebar.subheader("6 · Compute")
    quick = st.sidebar.toggle("Quick mode (fewer permutations)", value=True,
                              help="1000 → 100 permutations for a fast preview.")
    do_ora = st.sidebar.checkbox("Run ORA", value=True)
    do_gsea = st.sidebar.checkbox("Run GSEA", value=True)
    run = st.sidebar.button("▶ Run enrichment", type="primary", use_container_width=True)

    st.sidebar.caption("Uploaded data is processed in memory and not stored.")

    return dict(
        use_sample=use_sample, uploads=uploads, padj=padj, lfc=lfc, id_col=id_col,
        rank_metric=rank_metric, ora_sources=ora_sources, ora_dirs=tuple(ora_dirs),
        gsea_libs=gsea_libs, custom_gmt_file=custom_gmt_file,
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
        )
        results.append(res)
    prog.progress(1.0, text="Done")
    prog.empty()
    return results


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
        st.plotly_chart(plots.volcano(res.df, p["padj"], p["lfc"]),
                        use_container_width=True)
    with right:
        st.plotly_chart(plots.ma_plot(res.df, p["padj"], p["lfc"]),
                        use_container_width=True)


def tab_ora(res: ContrastResult):
    st.subheader(f"Over-representation (ORA) — {res.name}")
    if "ora" in res.errors:
        st.error(f"ORA failed: {res.errors['ora']}")
        return
    if res.ora is None or res.ora.empty:
        st.info("No ORA results. Run ORA with at least one source selected.")
        return
    top_n = st.slider("Top terms per direction", 5, 40, 15, key=f"ora_n_{res.name}")
    xmetric = st.radio("X axis", ["gene_ratio", "recall", "neg_log10_p"],
                       horizontal=True, key=f"ora_x_{res.name}")
    st.plotly_chart(plots.ora_dotplot(res.ora, top_n, xmetric),
                    use_container_width=True)
    st.dataframe(res.ora, use_container_width=True, height=320)
    st.download_button("⬇ ORA table (CSV)", res.ora.to_csv(index=False),
                       file_name=f"{res.name}_ORA.csv", mime="text/csv")


def tab_gsea(res: ContrastResult):
    st.subheader(f"Gene-set enrichment (GSEA) — {res.name}")
    if "gsea" in res.errors:
        st.error(f"GSEA failed: {res.errors['gsea']}")
        return
    if res.gsea is None or res.gsea.table.empty:
        st.info("No GSEA results. Select at least one library or upload a GMT.")
        return
    table = res.gsea.table
    top_n = st.slider("Top terms by |NES|", 5, 40, 20, key=f"gsea_n_{res.name}")
    st.plotly_chart(plots.gsea_bar(table, top_n), use_container_width=True)

    st.markdown("**Running-enrichment plot**")
    term = st.selectbox("Term", table["term"].tolist(), key=f"gsea_term_{res.name}")
    st.plotly_chart(plots.gsea_running_plot(res.gsea, term),
                    use_container_width=True)

    st.markdown("**Leading-edge genes**")
    default_terms = table.head(5)["term"].tolist()
    sel = st.multiselect("Terms for heatmap", table["term"].tolist(),
                         default=default_terms, key=f"gsea_le_{res.name}")
    if sel:
        st.plotly_chart(plots.leading_edge_heatmap(res.gsea, sel),
                        use_container_width=True)

    st.markdown("**Enrichment map**")
    jac = st.slider("Jaccard edge threshold", 0.05, 0.6, 0.25, 0.05,
                    key=f"gsea_jac_{res.name}")
    net_terms = table.head(40)["term"].tolist()
    term_map = {t: res.gsea.raw.results[t]["lead_genes"].split(";") for t in net_terms}
    scores = dict(zip(table["term"], table["NES"]))
    st.plotly_chart(plots.enrichment_network(term_map, scores, jac),
                    use_container_width=True)

    st.dataframe(table, use_container_width=True, height=320)
    st.download_button("⬇ GSEA table (CSV)", table.to_csv(index=False),
                       file_name=f"{res.name}_GSEA.csv", mime="text/csv")


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
        with st.spinner("Running enrichment…"):
            st.session_state["results"] = execute(p)
        st.session_state["params"] = p

    results = st.session_state.get("results")
    if not results:
        st.title("DESeq2 → ORA / GSEA")
        st.markdown(
            "Upload one or more **DESeq2 result tables** (or use the bundled sample), "
            "set thresholds in the sidebar, and click **Run enrichment**.\n\n"
            "* **ORA** runs on native chicken annotations via g:Profiler "
            "(GO / KEGG / Reactome / WikiPathways) with your tested genes as background.\n"
            "* **GSEA** ranks by the DESeq2 Wald statistic, maps chicken → human "
            "orthologs, and scores against MSigDB / Reactome collections.\n"
            "* Add a **custom .gmt** to score your own curated pathway modules."
        )
        return

    params = st.session_state.get("params", p)
    names = [r.name for r in results]
    chosen = st.selectbox("Contrast", names) if len(names) > 1 else names[0]
    res = next(r for r in results if r.name == chosen)

    if res.errors.get("load"):
        st.error(f"Could not load {res.name}: missing columns {res.errors['load']}")
        return

    t_qc, t_ora, t_gsea, t_cmp, t_exp = st.tabs(
        ["Upload & QC", "ORA", "GSEA", "Compare", "Export"]
    )
    with t_qc:
        tab_qc(res, params)
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
