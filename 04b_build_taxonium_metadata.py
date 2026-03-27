#!/usr/bin/env python3
"""
04b_build_taxonium_metadata.py
Rebuild Taxonium metadata using GTDB taxonomy (much better than ENA for bacteria)
and produce two complementary output files:

  A) taxonium_metadata.tsv.gz  — one row per sample, summary columns
       sample_id, species, genus, family, order, phylum,
       cluster (GTDB subspecies cluster), is_representative,
       amr_gene_count, amr_genes, amr_class_summary, amr_present, …

  B) taxonium_metadata_binary.tsv.gz — one row per sample, binary AMR matrix
       sample_id, species, cluster, <gene1>, <gene2>, … (0/1 per top-N gene)
       Mirrors the R script approach but with intrinsic genes filtered out.
       Top 500 genes, top 100 drug classes, top 200 subclasses.

Inputs:
  data/raw/atb_tree/ATB_tree_taxonomic_assignments_GTDB.csv
  data/processed/amr_per_sample.parquet

Outputs:
  data/processed/taxonium_metadata.tsv.gz
  data/processed/taxonium_metadata_binary.tsv.gz
"""

from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np

GTDB_FILE  = Path(__file__).parent.parent / "data" / "raw" / "atb_tree" / "ATB_tree_taxonomic_assignments_GTDB.csv"
AMR_TABLE  = Path(__file__).parent.parent / "data" / "processed" / "amr_per_sample.parquet"
STATUS_FILE = Path(__file__).parent.parent / "data" / "raw" / "amrfp" / "AMRFP_status.tsv.gz"
AMRFP_RAW  = Path(__file__).parent.parent / "data" / "raw" / "amrfp" / "AMRFP_results.tsv.gz"

OUT_SUMMARY = Path(__file__).parent.parent / "data" / "processed" / "taxonium_metadata.tsv.gz"
OUT_BINARY  = Path(__file__).parent.parent / "data" / "processed" / "taxonium_metadata_binary.tsv.gz"

# Top-N settings for binary matrix (mirrors R script, extended)
TOP_GENES    = 500
TOP_CLASSES  = 100
TOP_SUBCLASS = 200


def load_gtdb(path: Path) -> pd.DataFrame:
    print("Loading GTDB taxonomy …")
    df = pd.read_csv(path, dtype=str, low_memory=False)
    df = df.rename(columns={"Sample": "sample_id"})
    # Derive a clean 'subspecies' field from Species + Cluster
    # Species already includes genus, e.g. 'Salmonella enterica'
    # Cluster is the GTDB cluster number — useful for subspecies-level coloring
    print(f"  {len(df):,} samples, {df['Species'].nunique():,} unique species")
    return df


def load_amr_summary(path: Path) -> pd.DataFrame:
    print("Loading AMR per-sample summary …")
    df = pd.read_parquet(path)
    print(f"  {len(df):,} samples with acquired AMR genes")
    return df


def build_summary(gtdb: pd.DataFrame, amr: pd.DataFrame, status: pd.DataFrame) -> pd.DataFrame:
    """Merge GTDB + AMR + status into one-row-per-sample summary."""
    print("Building summary metadata …")
    merged = gtdb.merge(amr, on="sample_id", how="left")

    if status is not None:
        merged = merged.merge(status, on="sample_id", how="left")
        merged["amrfp_status"] = merged["amrfp_status"].fillna("not_run")

    merged["amr_gene_count"]       = merged["amr_gene_count"].fillna(0).astype(int)
    merged["amr_genes"]            = merged["amr_genes"].fillna("")
    merged["amr_class_summary"]    = merged["amr_class_summary"].fillna("")
    merged["amr_subclass_summary"] = merged.get("amr_subclass_summary",
                                                  pd.Series("", index=merged.index)).fillna("")
    merged["amr_present"]          = merged["amr_present"].fillna(False).astype(bool)

    amr_pos = merged["amr_present"].sum()
    print(f"  {amr_pos:,} / {len(merged):,} samples have acquired AMR genes ({amr_pos/len(merged)*100:.1f}%)")
    return merged


def build_binary_matrix(gtdb: pd.DataFrame, amrfp_raw_path: Path,
                         intrinsic_path: Path | None = None) -> pd.DataFrame:
    """
    Build wide binary presence-absence matrix for top AMR genes/classes,
    with intrinsic genes already filtered (Scope=plus rows only from raw file).
    """
    print("Building binary AMR matrix …")
    print("  Loading raw AMRFinderPlus results (plus-scope only) …")

    cols_needed = ["Name", "Gene symbol", "Class", "Subclass", "Scope", "Element type"]
    raw = pd.read_csv(amrfp_raw_path, sep="\t", dtype=str, low_memory=False,
                      usecols=lambda c: c in cols_needed)

    # Keep only acquired AMR (Scope=plus, Element type=AMR)
    raw = raw[
        (raw["Scope"].str.lower() == "plus") &
        (raw["Element type"].str.upper() == "AMR")
    ].copy()
    print(f"  {len(raw):,} acquired AMR rows across {raw['Name'].nunique():,} samples")

    results = gtdb[["sample_id", "Family", "Genus", "Species"]].copy()

    for level, col, top_n, prefix in [
        ("Gene symbol", "Gene symbol", TOP_GENES,    "gene_"),
        ("Class",       "Class",       TOP_CLASSES,  "class_"),
        ("Subclass",    "Subclass",    TOP_SUBCLASS, "subclass_"),
    ]:
        print(f"  Building {level} matrix (top {top_n}) …")
        # Split on "/" separator before counting (mirrors R: strsplit(Target, "/"))
        sub = raw[["Name", col]].dropna().copy()
        sub = sub.assign(**{col: sub[col].str.split("/")}).explode(col)
        sub[col] = sub[col].str.strip()

        top_vals = (sub[col].dropna()
                            .value_counts()
                            .head(top_n)
                            .index.tolist())

        sub = sub[sub[col].isin(top_vals)].drop_duplicates()
        sub["_val"] = 1
        wide = sub.pivot_table(index="Name", columns=col, values="_val",
                               aggfunc="max", fill_value=0)
        wide.columns = [f"{prefix}{c.replace(' ', '_').replace('/', '_')}"
                        for c in wide.columns]
        wide = wide.reset_index().rename(columns={"Name": "sample_id"})
        results = results.merge(wide, on="sample_id", how="left")
        # Fill missing (samples with no hits in this category) with 0
        binary_cols = [c for c in results.columns if c.startswith(prefix)]
        results[binary_cols] = results[binary_cols].fillna(0).astype(np.int8)

    print(f"  Binary matrix: {len(results):,} rows, {len(results.columns):,} columns")
    return results


def main():
    # 1. Load GTDB taxonomy
    gtdb = load_gtdb(GTDB_FILE)

    # 2. Load AMR summary
    amr = load_amr_summary(AMR_TABLE)

    # 3. Load AMRFP run status
    status = None
    if STATUS_FILE.exists():
        st = pd.read_csv(STATUS_FILE, sep="\t", dtype=str)
        status = st.rename(columns={"sample": "sample_id", "status": "amrfp_status"})
        print(f"Loaded AMRFP status for {len(status):,} samples.")

    # 4. Build summary metadata
    summary = build_summary(gtdb, amr, status)

    # Reorder columns — sample_id first (Taxonium requirement)
    priority = ["sample_id", "Species", "Genus", "Family", "Order", "Phylum", "Domain",
                "Cluster", "is_representative",
                "amr_gene_count", "amr_genes", "amr_class_summary",
                "amr_subclass_summary", "amr_present", "amr_highest_identity",
                "amrfp_status"]
    cols = priority + [c for c in summary.columns if c not in priority]
    summary = summary[[c for c in cols if c in summary.columns]]

    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    print(f"\nSaved summary: {OUT_SUMMARY}  ({len(summary):,} rows, {len(summary.columns)} cols)")

    # 5. Build binary matrix
    if AMRFP_RAW.exists():
        binary = build_binary_matrix(gtdb, AMRFP_RAW)
        binary.to_csv(OUT_BINARY, sep="\t", index=False)
        print(f"Saved binary matrix: {OUT_BINARY}  ({len(binary):,} rows, {len(binary.columns):,} cols)")
    else:
        print(f"WARNING: {AMRFP_RAW} not found — skipping binary matrix.")

    print("\nDone.")
    print(f"\nSpecies breakdown (top 15):")
    print(summary["Species"].value_counts().head(15).to_string())


if __name__ == "__main__":
    main()
