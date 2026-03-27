#!/usr/bin/env python3
"""
ATB AMR Explorer — Streamlit app
Interactive exploration of AMRFinderPlus annotations on the AllTheBacteria tree.

Data loaded from GitHub Release assets.
"""

import gzip
import io
import requests
import pandas as pd
import streamlit as st
import plotly.express as px

# ── Release assets ──────────────────────────────────────────────────────────────
RELEASE_BASE = (
    "https://github.com/AMR-genomics-hackathon-2026/atb_tree"
    "/releases/latest/download"
)

DATASETS = {
    "Summary metadata (AMR class + GTDB taxonomy)": {
        "url": f"{RELEASE_BASE}/taxonium_metadata.tsv.gz",
        "desc": "One row per sample. Columns: species, family, genus, AMR class summary, gene count, AMR presence.",
        "sep": "\t",
        "mode": "summary",
    },
    "Binary gene matrix (top 500 genes + classes)": {
        "url": f"{RELEASE_BASE}/taxonium_metadata_binary.tsv.gz",
        "desc": "One row per sample. Binary 0/1 columns for top 500 genes, 100 drug classes, 200 subclasses.",
        "sep": "\t",
        "mode": "binary",
    },
    "Drug class matrix — R script output": {
        "url": f"{RELEASE_BASE}/AMRFP_res_for_atb_tree_Class.csv.gz",
        "desc": "Binary presence-absence by drug class, produced by AMRFP_to_treemeta.R (Rebecca Gladstone).",
        "sep": ",",
        "mode": "binary",
    },
}

TAXONIUM_URL = (
    "https://taxonium.org/atb?addMetadataUrl="
    f"{RELEASE_BASE}/taxonium_metadata.tsv.gz"
)
TAXONIUM_BINARY_URL = (
    "https://taxonium.org/atb?addMetadataUrl="
    f"{RELEASE_BASE}/taxonium_metadata_binary.tsv.gz"
)

st.set_page_config(
    page_title="ATB AMR Explorer",
    page_icon="🧬",
    layout="wide",
)


# ── Data loading ────────────────────────────────────────────────────────────────

@st.cache_data(show_spinner="Downloading data from GitHub Release…")
def load_data(url: str, sep: str) -> pd.DataFrame:
    # Accept-Encoding: identity prevents requests from auto-decompressing,
    # so we always receive raw gzip bytes regardless of server headers
    resp = requests.get(url, allow_redirects=True, timeout=300,
                        headers={"Accept-Encoding": "identity"})
    resp.raise_for_status()
    content = gzip.decompress(resp.content)
    return pd.read_csv(io.BytesIO(content), sep=sep, dtype=str, low_memory=False)


def prepare_summary(df: pd.DataFrame) -> pd.DataFrame:
    df["amr_gene_count"] = pd.to_numeric(df.get("amr_gene_count", 0), errors="coerce").fillna(0).astype(int)
    df["amr_present"] = df.get("amr_present", pd.Series("False", index=df.index)).map(
        {"True": True, "False": False, True: True, False: False}
    )
    return df


# ── Header ──────────────────────────────────────────────────────────────────────

st.title("🧬 AllTheBacteria — AMR Annotation Explorer")
st.markdown(
    "Explore AMRFinderPlus annotations across ~2.4 M bacterial genomes, "
    "filtered using [AMRrules](https://github.com/AMRverse/AMRrules) to exclude "
    "intrinsic wildtype-susceptible genes."
)

col1, col2, col3 = st.columns(3)
col1.link_button("▶ View in Taxonium (summary)", TAXONIUM_URL, use_container_width=True)
col2.link_button("▶ View in Taxonium (gene matrix)", TAXONIUM_BINARY_URL, use_container_width=True)
col3.link_button("GitHub Repo", "https://github.com/AMR-genomics-hackathon-2026/atb_tree", use_container_width=True)

st.divider()

# ── Dataset selector ─────────────────────────────────────────────────────────────

st.subheader("Select dataset")
selected_name = st.radio(
    "Choose which metadata file to explore:",
    options=list(DATASETS.keys()),
    index=0,
    horizontal=False,
)
meta = DATASETS[selected_name]
st.caption(f"**Source:** [{meta['url'].split('/')[-1]}]({meta['url']})  \n{meta['desc']}")

st.divider()

# ── Load data ───────────────────────────────────────────────────────────────────

try:
    df_raw = load_data(meta["url"], meta["sep"])
except Exception as e:
    st.error(f"Could not load dataset: {e}")
    st.stop()

is_summary = meta["mode"] == "summary"

if is_summary:
    df = prepare_summary(df_raw.copy())
else:
    df = df_raw.copy()

# ── Sidebar filters ─────────────────────────────────────────────────────────────

st.sidebar.header("Filters")

# Species filter (available in both modes if Species column exists)
if "Species" in df.columns:
    all_species = sorted(df["Species"].dropna().unique())
    selected_species = st.sidebar.multiselect("Species", all_species, placeholder="All species")
else:
    selected_species = []

if is_summary:
    # Drug class filter
    all_classes = sorted(
        {c.strip() for row in df.get("amr_class_summary", pd.Series(dtype=str)).dropna()
         for c in row.split(";") if c.strip()}
    )
    selected_classes = st.sidebar.multiselect("Drug class (any)", all_classes, placeholder="All drug classes")

    # AMR presence filter
    amr_filter = st.sidebar.radio("AMR status", ["All", "AMR positive", "AMR negative"])
else:
    selected_classes = []
    amr_filter = "All"

    # Binary matrix: filter by presence of specific columns
    binary_cols = [c for c in df.columns if c.startswith(("gene_", "class_", "subclass_"))
                   or (c not in ("Sample", "sample_id", "Species", "Genus", "Family", "Name")
                       and df[c].isin(["0", "1", 0, 1]).all())]
    if binary_cols:
        selected_cols = st.sidebar.multiselect(
            "Show only samples positive for…",
            options=sorted(binary_cols),
            placeholder="Any gene/class",
        )
    else:
        selected_cols = []

# Sample search
sample_col = "sample_id" if "sample_id" in df.columns else ("Sample" if "Sample" in df.columns else "Name")
sample_search = st.sidebar.text_input("Search sample ID", placeholder="e.g. ERR123456")

# ── Apply filters ────────────────────────────────────────────────────────────────

filtered = df.copy()

if selected_species and "Species" in filtered.columns:
    filtered = filtered[filtered["Species"].isin(selected_species)]

if is_summary:
    if selected_classes:
        pattern = "|".join(selected_classes)
        filtered = filtered[filtered["amr_class_summary"].str.contains(pattern, na=False)]
    if amr_filter == "AMR positive":
        filtered = filtered[filtered["amr_present"] == True]
    elif amr_filter == "AMR negative":
        filtered = filtered[filtered["amr_present"] == False]
else:
    if selected_cols:
        for col in selected_cols:
            if col in filtered.columns:
                filtered = filtered[filtered[col].astype(str).isin(["1", "1.0"])]

if sample_search and sample_col in filtered.columns:
    filtered = filtered[filtered[sample_col].str.contains(sample_search, na=False)]

# ── Summary metrics ─────────────────────────────────────────────────────────────

total = len(filtered)

if is_summary:
    amr_pos = int(filtered["amr_present"].sum())
    pct = amr_pos / total * 100 if total > 0 else 0
    n_species = filtered["Species"].nunique() if "Species" in filtered.columns else "—"

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Samples", f"{total:,}")
    m2.metric("AMR positive", f"{amr_pos:,}")
    m3.metric("AMR prevalence", f"{pct:.1f}%")
    m4.metric("Species", f"{n_species:,}")
else:
    n_species = filtered["Species"].nunique() if "Species" in filtered.columns else "—"
    n_cols = len([c for c in filtered.columns if c not in ("Sample", "sample_id", "Species", "Genus", "Family", "Name")])

    m1, m2, m3 = st.columns(3)
    m1.metric("Samples", f"{total:,}")
    m2.metric("Species", f"{n_species:,}")
    m3.metric("AMR columns", f"{n_cols:,}")

st.divider()

# ── Charts (summary mode) ────────────────────────────────────────────────────────

if is_summary:
    chart_col1, chart_col2 = st.columns(2)

    with chart_col1:
        st.subheader("AMR prevalence by species (top 20)")
        species_stats = (
            filtered.groupby("Species")
            .agg(total=("sample_id", "count"), amr_pos=("amr_present", "sum"))
            .assign(pct=lambda x: x["amr_pos"] / x["total"] * 100)
            .sort_values("total", ascending=False)
            .head(20)
            .reset_index()
        )
        fig1 = px.bar(
            species_stats, x="pct", y="Species", orientation="h",
            color="pct", color_continuous_scale="Reds",
            labels={"pct": "% AMR positive", "Species": ""},
            hover_data={"total": True, "amr_pos": True},
        )
        fig1.update_layout(height=500, coloraxis_showscale=False, yaxis={"autorange": "reversed"})
        st.plotly_chart(fig1, use_container_width=True)

    with chart_col2:
        st.subheader("Drug class distribution (top 20)")
        class_counts = (
            filtered["amr_class_summary"].dropna()
            .str.split(";").explode().str.strip()
            .loc[lambda s: s != ""]
            .value_counts().head(20).reset_index()
        )
        class_counts.columns = ["Drug class", "Samples"]
        fig2 = px.bar(
            class_counts, x="Samples", y="Drug class", orientation="h",
            color="Samples", color_continuous_scale="Blues",
            labels={"Samples": "Samples with gene", "Drug class": ""},
        )
        fig2.update_layout(height=500, coloraxis_showscale=False, yaxis={"autorange": "reversed"})
        st.plotly_chart(fig2, use_container_width=True)

    st.subheader("AMR gene count distribution")
    fig3 = px.histogram(
        filtered[filtered["amr_present"] == True],
        x="amr_gene_count", nbins=50,
        labels={"amr_gene_count": "Number of acquired AMR genes", "count": "Samples"},
        color_discrete_sequence=["#e07b39"],
    )
    fig3.update_layout(height=300)
    st.plotly_chart(fig3, use_container_width=True)

else:
    # Binary mode: heatmap of top prevalent genes across top species
    st.subheader("Gene/class prevalence across species (top 20 × top 20)")
    binary_data_cols = [c for c in filtered.columns
                        if c not in ("Sample", "sample_id", "Species", "Genus", "Family", "Name")][:500]
    if "Species" in filtered.columns and binary_data_cols:
        heat_df = filtered.copy()
        heat_df[binary_data_cols] = heat_df[binary_data_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        top_genes = heat_df[binary_data_cols].sum().sort_values(ascending=False).head(20).index.tolist()
        top_species = heat_df["Species"].value_counts().head(20).index.tolist()
        heat = (
            heat_df[heat_df["Species"].isin(top_species)]
            .groupby("Species")[top_genes].mean() * 100
        )
        fig4 = px.imshow(
            heat, aspect="auto", color_continuous_scale="Blues",
            labels={"color": "% positive"},
        )
        fig4.update_layout(height=500)
        st.plotly_chart(fig4, use_container_width=True)

st.divider()

# ── Data table ───────────────────────────────────────────────────────────────────

st.subheader(f"Sample table ({total:,} rows)")

if is_summary:
    display_cols = ["sample_id", "Species", "Genus", "Family",
                    "amr_present", "amr_gene_count", "amr_class_summary", "amr_genes"]
    display_cols = [c for c in display_cols if c in filtered.columns]
else:
    # Show taxonomy + first 20 binary columns
    meta_cols = [c for c in ("Sample", "sample_id", "Species", "Genus", "Family") if c in filtered.columns]
    val_cols = [c for c in filtered.columns if c not in meta_cols][:20]
    display_cols = meta_cols + val_cols

st.dataframe(filtered[display_cols].reset_index(drop=True), use_container_width=True, height=400)

csv_bytes = filtered[display_cols].to_csv(index=False).encode()
st.download_button(
    label="Download filtered table (CSV)",
    data=csv_bytes,
    file_name="atb_amr_filtered.csv",
    mime="text/csv",
)

st.caption(
    "Data: [AllTheBacteria](https://allthebacteria.readthedocs.io/) · "
    "[AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) · "
    "[AMRrules](https://github.com/AMRverse/AMRrules) · "
    "[GTDB taxonomy](https://gtdb.ecogenomic.org/)"
)
