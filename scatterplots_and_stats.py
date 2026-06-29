import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr, false_discovery_control

# ── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(page_title="Library Comparisons", layout="wide")

# ── Data paths ───────────────────────────────────────────────────────────────
LIBRARY_PATHS = {
    "MinLib": "/lustre/scratch126/casm/teams/team215/projects/Genmap/library2_analysis/n_guides_analysis/minlib_guides_only/results/GI_calls.txt",
    "Library 1": "/lustre/scratch126/casm/teams/team215/projects/Genmap/library1_results/GI_calls.txt",
    "Library 2": "/lustre/scratch126/casm/teams/team215/projects/Genmap/library2_analysis/results/GI_calls.txt",
}

METRICS = ["dLFC", "HSA", "lfc_combo"]


# ── Cached data loader ──────────────────────────────────────────────────────
@st.cache_data
def load_library(name: str) -> pd.DataFrame:
    return pd.read_csv(LIBRARY_PATHS[name], sep="\t", usecols=["gene_combo", "lfc_combo", "dLFC", "HSA"])


# ── Helper functions ─────────────────────────────────────────────────────────
def common_gene_combos(df1: pd.DataFrame, df2: pd.DataFrame) -> set:
    """Return the set of gene_combo values present in both dataframes."""
    return set(df1["gene_combo"]).intersection(df2["gene_combo"])


def correlation_stats(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    label1: str = "df1",
    label2: str = "df2",
    metrics: list[str] | None = None,
    p_adjust_method: str = "bh",
) -> pd.DataFrame:
    if metrics is None:
        metrics = METRICS
    common = common_gene_combos(df1, df2)
    d1 = df1[df1["gene_combo"].isin(common)].set_index("gene_combo")
    d2 = df2[df2["gene_combo"].isin(common)].set_index("gene_combo")
    d1, d2 = d1.align(d2, join="inner")

    rows = []
    for metric in metrics:
        r_spear, p_spear = spearmanr(d1[metric], d2[metric])
        r_pear, p_pear = pearsonr(d1[metric], d2[metric])
        rows.append(
            {
                "metric": metric,
                "n_pairs": len(d1),
                "spearman_r": r_spear,
                "spearman_p": p_spear,
                "pearson_r": r_pear,
                "pearson_p": p_pear,
            }
        )
    result = pd.DataFrame(rows)
    result["spearman_p_adj"] = false_discovery_control(result["spearman_p"].values, method=p_adjust_method)
    result["pearson_p_adj"] = false_discovery_control(result["pearson_p"].values, method=p_adjust_method)
    result.attrs["comparison"] = f"{label1} vs {label2}"
    return result


def plot_metric_scatter(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    label1: str = "df1",
    label2: str = "df2",
    metrics: list[str] | None = None,
    alpha: float = 0.3,
    s: float = 8,
) -> plt.Figure:
    if metrics is None:
        metrics = METRICS
    common = common_gene_combos(df1, df2)
    d1 = df1[df1["gene_combo"].isin(common)].set_index("gene_combo")
    d2 = df2[df2["gene_combo"].isin(common)].set_index("gene_combo")
    d1, d2 = d1.align(d2, join="inner")

    fig, axes = plt.subplots(1, len(metrics), figsize=(6 * len(metrics), 5))
    if len(metrics) == 1:
        axes = [axes]

    for ax, metric in zip(axes, metrics):
        x = d1[metric].values
        y = d2[metric].values

        ax.scatter(x, y, alpha=alpha, s=s, edgecolors="none", rasterized=True)

        lo = min(x.min(), y.min())
        hi = max(x.max(), y.max())
        ax.plot([lo, hi], [lo, hi], ls="--", c="grey", lw=0.8)

        r_spear, p_spear = spearmanr(x, y)
        r_pear, p_pear = pearsonr(x, y)
        ax.annotate(
            f"Spearman r = {r_spear:.3f} (p = {p_spear:.2e})\n"
            f"Pearson  r = {r_pear:.3f} (p = {p_pear:.2e})",
            xy=(0.05, 0.95),
            xycoords="axes fraction",
            va="top",
            ha="left",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="grey", alpha=0.8),
        )
        ax.set_xlabel(label1)
        ax.set_ylabel(label2)
        ax.set_title(metric)

    fig.suptitle(f"{label1} vs {label2}", y=1.02, fontsize=13)
    fig.tight_layout()
    return fig


# ── Sidebar controls ────────────────────────────────────────────────────────
st.sidebar.header("Comparison Settings")
library_names = list(LIBRARY_PATHS.keys())

lib_a = st.sidebar.selectbox("Library A", library_names, index=0)
lib_b = st.sidebar.selectbox("Library B", library_names, index=1)
metric = st.sidebar.selectbox("Metric", METRICS)

if lib_a == lib_b:
    st.error("Please select two different libraries to compare.")
    st.stop()

# ── Load data ────────────────────────────────────────────────────────────────
df_a = load_library(lib_a)
df_b = load_library(lib_b)

# ── Main panel ───────────────────────────────────────────────────────────────
st.title(f"{lib_a} vs {lib_b}")

common = common_gene_combos(df_a, df_b)
st.metric("Common gene combos", f"{len(common):,}")

# Scatter plot
fig = plot_metric_scatter(df_a, df_b, label1=lib_a, label2=lib_b, metrics=[metric])
st.pyplot(fig)

# Correlation stats table
stats = correlation_stats(df_a, df_b, label1=lib_a, label2=lib_b, metrics=[metric])
st.subheader("Correlation Statistics")
st.dataframe(stats, use_container_width=True)
