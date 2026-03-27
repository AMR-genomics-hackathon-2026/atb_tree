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
SUMMARY_URL = f"{RELEASE_BASE}/taxonium_metadata.tsv.gz"
TAXONIUM_URL = (
    "https://taxonium.org/atb?addMetadataUrl="
    f"{RELEASE_BASE}/taxonium_metadata.tsv.gz"
)
TAXONIUM_BINARY_URL = (
    "https://taxonium.org/atb?addMetadataUrl="
    f"{RELEASE_BASE}/taxonium_metadata_binary.tsv.gz"
)

DOWNLOADS = {
    "taxonium_metadata.tsv.gz": "Summary metadata (AMR class + GTDB taxonomy)",
    "taxonium_metadata_binary.tsv.gz": "Binary gene matrix (for Taxonium)",
    "AMRFP_res_for_atb_tree_Class.csv.gz": "Drug class matrix (R script output)",
    "meta.tar.gz": "Full sample metadata (Rebecca Gladstone)",
}

st.set_page_config(page_title="ATB AMR Explorer", page_icon="🧬", layout="wide")


# ── Data loading ────────────────────────────────────────────────────────────────

@st.cache_data(show_spinner="Downloading metadata (~13 MB)…")
def load_summary() -> pd.DataFrame:
    resp = requests.get(SUMMARY_URL, allow_redirects=True, timeout=300,
                        headers={"Accept-Encoding": "identity"})
    resp.raise_for_status()
    content = gzip.decompress(resp.content)
    df = pd.read_csv(io.BytesIO(content), sep="\t", dtype=str, low_memory=False)
    df["amr_gene_count"] = pd.to_numeric(df["amr_gene_count"], errors="coerce").fillna(0).astype(int)
    df["amr_present"] = df["amr_present"].map({"True": True, "False": False})
    return df


# ── Header ──────────────────────────────────────────────────────────────────────

st.title("🧬 AllTheBacteria — AMR Annotation Explorer")
st.markdown(
    "Explore AMRFinderPlus annotations across ~2.4 M bacterial genomes, "
    "filtered using [AMRrules](https://github.com/AMRverse/AMRrules) to exclude "
    "intrinsic wildtype-susceptible genes."
)

c1, c2, c3 = st.columns(3)
c1.link_button("▶ View tree (summary)", TAXONIUM_URL, use_container_width=True)
c2.link_button("▶ View tree (gene matrix)", TAXONIUM_BINARY_URL, use_container_width=True)
c3.link_button("GitHub Repo", "https://github.com/AMR-genomics-hackathon-2026/atb_tree", use_container_width=True)

st.divider()

# ── Load data ───────────────────────────────────────────────────────────────────

try:
    df = load_summary()
except Exception as e:
    st.error(f"Could not load metadata: {e}")
    st.stop()

# ── View selector ────────────────────────────────────────────────────────────────

view = st.radio("View", ["Overview", "By gene", "By drug class", "Species deep-dive", "Downloads"],
                horizontal=True)

st.divider()

# ── Sidebar filters ──────────────────────────────────────────────────────────────

st.sidebar.header("Filters")

all_species = sorted(df["Species"].dropna().unique())
selected_species = st.sidebar.multiselect("Species", all_species, placeholder="All species")

all_classes = sorted(
    {c.strip() for row in df["amr_class_summary"].dropna()
     for c in row.split(";") if c.strip()}
)
selected_classes = st.sidebar.multiselect("Drug class", all_classes, placeholder="All classes")

amr_filter = st.sidebar.radio("AMR status", ["All", "AMR positive", "AMR negative"])
sample_search = st.sidebar.text_input("Search sample ID", placeholder="e.g. ERR123456")

# Apply filters
filtered = df.copy()
if selected_species:
    filtered = filtered[filtered["Species"].isin(selected_species)]
if selected_classes:
    pattern = "|".join(selected_classes)
    filtered = filtered[filtered["amr_class_summary"].str.contains(pattern, na=False)]
if amr_filter == "AMR positive":
    filtered = filtered[filtered["amr_present"] == True]
elif amr_filter == "AMR negative":
    filtered = filtered[filtered["amr_present"] == False]
if sample_search:
    filtered = filtered[filtered["sample_id"].str.contains(sample_search, na=False)]

total = len(filtered)
amr_pos = int(filtered["amr_present"].sum())
pct = amr_pos / total * 100 if total > 0 else 0

# ── Metrics ──────────────────────────────────────────────────────────────────────

m1, m2, m3, m4 = st.columns(4)
m1.metric("Samples", f"{total:,}")
m2.metric("AMR positive", f"{amr_pos:,}")
m3.metric("AMR prevalence", f"{pct:.1f}%")
m4.metric("Species", f"{filtered['Species'].nunique():,}")

st.divider()

# ── Views ────────────────────────────────────────────────────────────────────────

if view == "Overview":
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Top species by AMR-positive samples")
        sp = (
            filtered.groupby("Species")
            .agg(total=("sample_id", "count"), amr_pos=("amr_present", "sum"))
            .assign(pct=lambda x: (x["amr_pos"] / x["total"] * 100).round(1))
            .sort_values("amr_pos", ascending=False).head(20).reset_index()
        )
        fig = px.bar(sp, x="amr_pos", y="Species", orientation="h",
                     color="pct", color_continuous_scale="Reds",
                     hover_data={"total": True, "pct": True},
                     labels={"amr_pos": "AMR-positive samples", "Species": "",
                             "pct": "Prevalence %", "total": "Total samples"})
        fig.update_layout(height=500, coloraxis_showscale=True,
                          yaxis={"autorange": "reversed"})
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Top drug classes")
        cc = (
            filtered["amr_class_summary"].dropna()
            .str.split(";").explode().str.strip()
            .loc[lambda s: s != ""].value_counts().head(20).reset_index()
        )
        cc.columns = ["Drug class", "Samples"]
        fig2 = px.bar(cc, x="Samples", y="Drug class", orientation="h",
                      color="Samples", color_continuous_scale="Blues",
                      labels={"Drug class": ""})
        fig2.update_layout(height=500, coloraxis_showscale=False,
                           yaxis={"autorange": "reversed"})
        st.plotly_chart(fig2, use_container_width=True)

    st.subheader("AMR gene count distribution")
    fig3 = px.histogram(
        filtered[filtered["amr_present"] == True],
        x="amr_gene_count", nbins=50,
        labels={"amr_gene_count": "Acquired AMR genes per sample"},
        color_discrete_sequence=["#e07b39"],
    )
    fig3.update_layout(height=280)
    st.plotly_chart(fig3, use_container_width=True)

elif view == "By gene":
    st.subheader("Gene prevalence")

    # Explode amr_genes column
    gene_df = (
        filtered[["sample_id", "Species", "amr_genes"]].dropna(subset=["amr_genes"])
        .assign(gene=filtered["amr_genes"].str.split(";"))
        .explode("gene")
    )
    gene_df["gene"] = gene_df["gene"].str.strip()
    gene_df = gene_df[gene_df["gene"] != ""]

    # Gene filter
    all_genes = sorted(gene_df["gene"].unique())
    selected_genes = st.multiselect("Filter by gene", all_genes, placeholder="All genes")
    if selected_genes:
        gene_df = gene_df[gene_df["gene"].isin(selected_genes)]

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Top 20 genes by sample count**")
        top_genes = gene_df["gene"].value_counts().head(20).reset_index()
        top_genes.columns = ["Gene", "Samples"]
        fig_g = px.bar(top_genes, x="Samples", y="Gene", orientation="h",
                       color="Samples", color_continuous_scale="Greens",
                       labels={"Gene": ""})
        fig_g.update_layout(height=500, coloraxis_showscale=False,
                            yaxis={"autorange": "reversed"})
        st.plotly_chart(fig_g, use_container_width=True)

    with col2:
        st.markdown("**Gene prevalence across top 15 species**")
        top_sp = filtered["Species"].value_counts().head(15).index.tolist()
        top_gn = gene_df["gene"].value_counts().head(15).index.tolist()
        heat_data = (
            gene_df[gene_df["Species"].isin(top_sp) & gene_df["gene"].isin(top_gn)]
            .groupby(["Species", "gene"])["sample_id"].count()
            .unstack(fill_value=0)
        )
        if not heat_data.empty:
            # Normalise by species total
            sp_totals = filtered[filtered["Species"].isin(top_sp)]["Species"].value_counts()
            heat_pct = heat_data.div(sp_totals, axis=0).fillna(0) * 100
            fig_h = px.imshow(heat_pct, aspect="auto", color_continuous_scale="Blues",
                              labels={"color": "% samples"})
            fig_h.update_layout(height=500)
            st.plotly_chart(fig_h, use_container_width=True)

elif view == "By drug class":
    st.subheader("Drug class breakdown")

    class_df = (
        filtered[["sample_id", "Species", "amr_class_summary"]].dropna(subset=["amr_class_summary"])
        .assign(drug_class=filtered["amr_class_summary"].str.split(";"))
        .explode("drug_class")
    )
    class_df["drug_class"] = class_df["drug_class"].str.strip()
    class_df = class_df[class_df["drug_class"] != ""]

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Samples per drug class**")
        dc = class_df["drug_class"].value_counts().head(20).reset_index()
        dc.columns = ["Drug class", "Samples"]
        fig_dc = px.bar(dc, x="Samples", y="Drug class", orientation="h",
                        color="Samples", color_continuous_scale="Oranges",
                        labels={"Drug class": ""})
        fig_dc.update_layout(height=500, coloraxis_showscale=False,
                             yaxis={"autorange": "reversed"})
        st.plotly_chart(fig_dc, use_container_width=True)

    with col2:
        st.markdown("**Drug class prevalence across top 15 species**")
        top_sp = filtered["Species"].value_counts().head(15).index.tolist()
        top_dc = class_df["drug_class"].value_counts().head(10).index.tolist()
        heat_data = (
            class_df[class_df["Species"].isin(top_sp) & class_df["drug_class"].isin(top_dc)]
            .groupby(["Species", "drug_class"])["sample_id"].count()
            .unstack(fill_value=0)
        )
        if not heat_data.empty:
            sp_totals = filtered[filtered["Species"].isin(top_sp)]["Species"].value_counts()
            heat_pct = heat_data.div(sp_totals, axis=0).fillna(0) * 100
            fig_h2 = px.imshow(heat_pct, aspect="auto", color_continuous_scale="Oranges",
                               labels={"color": "% samples"})
            fig_h2.update_layout(height=500)
            st.plotly_chart(fig_h2, use_container_width=True)

elif view == "Species deep-dive":
    st.subheader("Species deep-dive")

    species_list = sorted(df["Species"].dropna().unique())
    selected_sp = st.selectbox("Select a species", species_list,
                               index=species_list.index("Escherichia coli") if "Escherichia coli" in species_list else 0)

    sp_df = df[df["Species"] == selected_sp].copy()
    sp_amr = sp_df[sp_df["amr_present"] == True]
    sp_total = len(sp_df)
    sp_pos = len(sp_amr)
    sp_pct = sp_pos / sp_total * 100 if sp_total > 0 else 0
    avg_genes = sp_amr["amr_gene_count"].mean() if sp_pos > 0 else 0

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Total samples", f"{sp_total:,}")
    m2.metric("AMR positive", f"{sp_pos:,}")
    m3.metric("AMR prevalence", f"{sp_pct:.1f}%")
    m4.metric("Avg genes (positive)", f"{avg_genes:.1f}")

    st.divider()
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Top AMR genes**")
        gene_counts = (
            sp_amr["amr_genes"].dropna()
            .str.split(";").explode().str.strip()
            .loc[lambda s: s != ""].value_counts().head(20).reset_index()
        )
        gene_counts.columns = ["Gene", "Samples"]
        gene_counts["% of species"] = (gene_counts["Samples"] / sp_total * 100).round(1)
        fig_sg = px.bar(gene_counts, x="% of species", y="Gene", orientation="h",
                        color="% of species", color_continuous_scale="Greens",
                        hover_data={"Samples": True},
                        labels={"Gene": "", "% of species": "% of all samples"})
        fig_sg.update_layout(height=500, coloraxis_showscale=False,
                             yaxis={"autorange": "reversed"})
        st.plotly_chart(fig_sg, use_container_width=True)

    with col2:
        st.markdown("**Drug class breakdown**")
        class_counts = (
            sp_amr["amr_class_summary"].dropna()
            .str.split(";").explode().str.strip()
            .loc[lambda s: s != ""].value_counts().head(20).reset_index()
        )
        class_counts.columns = ["Drug class", "Samples"]
        class_counts["% of species"] = (class_counts["Samples"] / sp_total * 100).round(1)
        fig_sc = px.bar(class_counts, x="% of species", y="Drug class", orientation="h",
                        color="% of species", color_continuous_scale="Reds",
                        hover_data={"Samples": True},
                        labels={"Drug class": "", "% of species": "% of all samples"})
        fig_sc.update_layout(height=500, coloraxis_showscale=False,
                             yaxis={"autorange": "reversed"})
        st.plotly_chart(fig_sc, use_container_width=True)

    col3, col4 = st.columns(2)

    with col3:
        # Pre-aggregate gene count distribution server-side
        st.markdown("**AMR gene burden distribution**")
        burden_pos = sp_amr["amr_gene_count"].value_counts().sort_index().reset_index()
        burden_pos.columns = ["Gene count", "Samples"]
        burden_pos["Status"] = "AMR positive"
        burden_neg = pd.DataFrame({"Gene count": [0],
                                   "Samples": [sp_total - sp_pos],
                                   "Status": ["AMR negative"]})
        burden = pd.concat([burden_neg, burden_pos], ignore_index=True)
        fig_dist = px.bar(burden, x="Gene count", y="Samples", color="Status",
                          color_discrete_map={"AMR positive": "#c0392b",
                                             "AMR negative": "#aaaaaa"},
                          labels={"Gene count": "Number of acquired AMR genes"},
                          barmode="overlay")
        fig_dist.update_layout(height=300)
        st.plotly_chart(fig_dist, use_container_width=True)

    with col4:
        # Pre-aggregate MDR profile server-side
        st.markdown("**Multi-drug resistance (MDR) profile**")
        n_classes = sp_amr["amr_class_summary"].str.count(";").add(1).fillna(0).astype(int)
        mdr_counts = n_classes.value_counts().sort_index().reset_index()
        mdr_counts.columns = ["Drug classes", "Samples"]
        mdr_counts["Category"] = mdr_counts["Drug classes"].apply(
            lambda n: "Single class" if n == 1 else
                      ("2 classes" if n == 2 else "MDR (3+ classes)")
        )
        fig_mdr = px.bar(mdr_counts, x="Drug classes", y="Samples", color="Category",
                         color_discrete_map={"Single class": "#f39c12",
                                            "2 classes": "#e67e22",
                                            "MDR (3+ classes)": "#c0392b"},
                         labels={"Drug classes": "Resistance drug classes"})
        fig_mdr.update_layout(height=300)
        st.plotly_chart(fig_mdr, use_container_width=True)

elif view == "Downloads":
    st.subheader("Download files")
    for fname, desc in DOWNLOADS.items():
        url = f"{RELEASE_BASE}/{fname}"
        st.markdown(f"**[{fname}]({url})** — {desc}")

st.divider()

# ── Sample table ─────────────────────────────────────────────────────────────────

TABLE_LIMIT = 2000
st.subheader(f"Sample table ({total:,} rows — showing first {min(total, TABLE_LIMIT):,})")
display_cols = [c for c in ["sample_id", "Species", "Genus", "Family",
                             "amr_present", "amr_gene_count",
                             "amr_class_summary", "amr_genes"] if c in filtered.columns]
st.dataframe(filtered[display_cols].head(TABLE_LIMIT).reset_index(drop=True),
             use_container_width=True, height=380)

if total > TABLE_LIMIT:
    st.caption(f"Table capped at {TABLE_LIMIT:,} rows. Use download for full dataset.")

csv_bytes = filtered[display_cols].to_csv(index=False).encode()
st.download_button(f"Download filtered table ({total:,} rows, CSV)",
                   csv_bytes, "atb_amr_filtered.csv", "text/csv")

st.caption(
    "Data: [AllTheBacteria](https://allthebacteria.readthedocs.io/) · "
    "[AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) · "
    "[AMRrules](https://github.com/AMRverse/AMRrules) · "
    "[GTDB taxonomy](https://gtdb.ecogenomic.org/)"
)
