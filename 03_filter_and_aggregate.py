#!/usr/bin/env python3
"""
03_filter_and_aggregate.py
Process raw AMRFinderPlus TSVs:
  1. Concatenate all per-sample TSVs (may be split by batch/species).
  2. Remove intrinsic WT-susceptible genes using the AMRrules lookup
     (data/processed/intrinsic_genes.parquet).
  3. Aggregate to one row per sample with summary columns:
       - amr_gene_count     : number of distinct acquired AMR genes
       - amr_genes          : semicolon-joined sorted gene symbols
       - amr_class_summary  : semicolon-joined sorted drug classes
       - amr_present        : boolean (true/false)
       - amr_highest_level  : highest resistance level (if column present)

Input:  data/raw/amrfp/*.tsv  (or .tsv.gz)
        data/processed/intrinsic_genes.parquet
Output: data/processed/amr_per_sample.parquet
        data/processed/amr_per_sample.tsv.gz   (Taxonium metadata input)
"""

import sys
import glob
from pathlib import Path

import pandas as pd

RAW_DIR        = Path(__file__).parent.parent / "data" / "raw" / "amrfp"
INTRINSIC_PATH = Path(__file__).parent.parent / "data" / "processed" / "intrinsic_genes.parquet"
OUT_PARQUET    = Path(__file__).parent.parent / "data" / "processed" / "amr_per_sample.parquet"
OUT_TSV        = Path(__file__).parent.parent / "data" / "processed" / "amr_per_sample.tsv.gz"

# AMRFinderPlus column names (standard output)
COL_SAMPLE   = "Name"                            # sample/genome identifier
COL_GENE     = "Gene symbol"
COL_CLASS    = "Class"                           # high-level drug class (BETA-LACTAM, etc.)
COL_SUBCLASS = "Subclass"                        # specific subclass (CEPHALOTHIN, etc.)
COL_METHOD   = "Method"
COL_SCOPE    = "Scope"                           # 'core' = intrinsic, 'plus' = acquired
COL_ELEMENT  = "Element type"                    # AMR, STRESS, VIRULENCE, …
COL_LEVEL    = "% Identity to reference sequence"
# Alias used in filter logic — no Taxonomic group col in this dataset
COL_SPECIES  = "Taxonomic group"

# We only want AMR-type hits
ELEMENT_TYPE_FILTER = "AMR"

# AMRFinderPlus 'Scope' column: 'core' = intrinsic, 'plus' = acquired
# Filtering on scope=='plus' is the most reliable intrinsic exclusion when available.
ACQUIRED_SCOPE = "plus"


def load_amrfp_files(raw_dir: Path) -> pd.DataFrame:
    patterns = [str(raw_dir / p) for p in ("*.tsv", "*.tsv.gz", "**/*.tsv", "**/*.tsv.gz")]
    paths = []
    for pat in patterns:
        paths.extend(glob.glob(pat, recursive=True))
    paths = list(set(paths))
    # Exclude status/QC files — they have a different schema (sample, status, comments)
    paths = [p for p in paths if "status" not in Path(p).name.lower()]

    if not paths:
        print(f"ERROR: No TSV files found in {raw_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {len(paths)} AMRFinderPlus file(s) …")
    dfs = []
    for p in sorted(paths):
        try:
            df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
            dfs.append(df)
        except Exception as e:
            print(f"  WARNING: could not read {p}: {e}", file=sys.stderr)

    combined = pd.concat(dfs, ignore_index=True)
    print(f"  Total rows (pre-filter): {len(combined):,}")
    return combined


def normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Handle minor column-name variations across AMRFinderPlus versions."""
    renames = {}
    for col in df.columns:
        low = col.lower()
        if low in ("name", "target identifier", "sample", "genome", "accession"):
            renames[col] = COL_SAMPLE
        elif low in ("gene symbol", "gene_symbol", "gene"):
            renames[col] = COL_GENE
        elif low in ("subclass",):
            renames[col] = COL_SUBCLASS
        elif low in ("class", "drug class", "drug_class"):
            renames[col] = COL_CLASS
        elif low in ("taxonomic group", "taxgroup", "taxonomy"):
            renames[col] = COL_SPECIES
        elif low in ("scope",):
            renames[col] = COL_SCOPE
        elif low in ("element type", "element_type"):
            renames[col] = COL_ELEMENT
    return df.rename(columns=renames)


def filter_intrinsic_by_scope(df: pd.DataFrame) -> pd.DataFrame:
    """If AMRFinderPlus 'Scope' column is present, keep only 'plus' (acquired)."""
    if COL_SCOPE in df.columns:
        before = len(df)
        df = df[df[COL_SCOPE].str.lower() == ACQUIRED_SCOPE].copy()
        removed = before - len(df)
        print(f"  Scope filter (core/intrinsic removed): {removed:,} rows")
    return df


def filter_intrinsic_by_amrrules(df: pd.DataFrame, intrinsic: pd.DataFrame) -> pd.DataFrame:
    """
    Remove gene–species pairs that appear in the AMRrules intrinsic lookup.
    Matching is done on gene_symbol (case-insensitive) and species (if available).
    """
    if intrinsic.empty:
        return df

    # Normalise gene symbols for matching
    df["_gene_lower"]    = df[COL_GENE].str.lower().str.strip()
    intrinsic["_gl"]     = intrinsic["gene_symbol"].str.lower().str.strip()

    if COL_SPECIES in df.columns:
        df["_species_lower"] = df[COL_SPECIES].str.lower().str.strip()
        intrinsic["_sl"]     = intrinsic["species"].str.lower().str.strip()

        merge = df.merge(
            intrinsic[["_gl", "_sl"]].drop_duplicates(),
            left_on=["_gene_lower", "_species_lower"],
            right_on=["_gl", "_sl"],
            how="left",
            indicator=True,
        )
        before = len(df)
        df = merge[merge["_merge"] == "left_only"].drop(columns=["_gl", "_sl", "_merge"])
        df = df.drop(columns=["_gene_lower", "_species_lower"], errors="ignore")
        print(f"  AMRrules intrinsic filter (species-aware): {before - len(df):,} rows removed")
    else:
        # Species column absent — filter on gene symbol only
        intrinsic_genes = set(intrinsic["_gl"])
        before = len(df)
        df = df[~df["_gene_lower"].isin(intrinsic_genes)].copy()
        df = df.drop(columns=["_gene_lower"], errors="ignore")
        print(f"  AMRrules intrinsic filter (gene-only): {before - len(df):,} rows removed")

    return df


def aggregate_per_sample(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse to one row per sample with AMR summary columns."""
    # Keep only AMR element type if column present
    if COL_ELEMENT in df.columns:
        df = df[df[COL_ELEMENT].str.upper() == ELEMENT_TYPE_FILTER].copy()

    groups = df.groupby(COL_SAMPLE, sort=False)

    agg = pd.DataFrame({
        "amr_gene_count":    groups[COL_GENE].nunique(),
        "amr_genes":         groups[COL_GENE].agg(
            lambda x: ";".join(sorted(x.dropna().unique()))
        ),
        "amr_class_summary": groups[COL_CLASS].agg(
            lambda x: ";".join(sorted(x.dropna().unique()))
        ) if COL_CLASS in df.columns else "",
        "amr_subclass_summary": groups[COL_SUBCLASS].agg(
            lambda x: ";".join(sorted(x.dropna().unique()))
        ) if COL_SUBCLASS in df.columns else "",
    })
    agg["amr_present"] = agg["amr_gene_count"] > 0

    if COL_LEVEL in df.columns:
        agg["amr_highest_identity"] = (
            df.groupby(COL_SAMPLE)[COL_LEVEL]
            .apply(lambda x: pd.to_numeric(x, errors="coerce").max())
            .round(1)
            .astype(str)
        )

    agg.index.name = "sample_id"
    agg = agg.reset_index()
    return agg


def main():
    # 1. Load raw AMRFinderPlus data
    raw = load_amrfp_files(RAW_DIR)
    raw = normalise_columns(raw)

    if COL_SAMPLE not in raw.columns:
        print(f"ERROR: Could not identify a sample/name column. Columns: {list(raw.columns)}", file=sys.stderr)
        sys.exit(1)

    # 2. Load AMRrules intrinsic lookup
    if INTRINSIC_PATH.exists():
        intrinsic = pd.read_parquet(INTRINSIC_PATH)
        print(f"Loaded {len(intrinsic):,} intrinsic gene–species pairs from AMRrules.")
    else:
        print("WARNING: intrinsic_genes.parquet not found. Run 01_fetch_amrrules.py first.")
        intrinsic = pd.DataFrame(columns=["species", "gene_symbol", "drug_class"])

    # 3. Filter
    raw = filter_intrinsic_by_scope(raw)
    raw = filter_intrinsic_by_amrrules(raw, intrinsic)
    print(f"  Rows after all filters: {len(raw):,}")

    # 4. Aggregate
    print("Aggregating per sample …")
    agg = aggregate_per_sample(raw)
    print(f"  Unique samples: {len(agg):,}")

    # 5. Save
    OUT_PARQUET.parent.mkdir(parents=True, exist_ok=True)
    agg.to_parquet(OUT_PARQUET, index=False)
    agg.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"Saved:\n  {OUT_PARQUET}\n  {OUT_TSV}")


if __name__ == "__main__":
    main()
