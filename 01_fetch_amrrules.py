#!/usr/bin/env python3
"""
01_fetch_amrrules.py
Download AMRrules TSV files from GitHub and build a lookup table of
intrinsic (wildtype-susceptible) gene–species pairs to exclude from
the AMRFinderPlus annotation.

Intrinsic genes are identified by:
  - `phenotype` column == "Susceptible" (or no resistance phenotype)  AND
  - `mutation_string` is empty (wildtype, not a specific mutation)
  - OR any rule explicitly flagged as intrinsic

Output: data/processed/intrinsic_genes.parquet
  Columns: species (str), gene_symbol (str), drug_class (str)
"""

import io
import re
import sys
import requests
import pandas as pd
from pathlib import Path

GITHUB_API = "https://api.github.com/repos/AMRverse/AMRrules/contents/rules"
RAW_BASE   = "https://raw.githubusercontent.com/AMRverse/AMRrules/main/rules/"

OUT_PATH = Path(__file__).parent.parent / "data" / "processed" / "intrinsic_genes.parquet"
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)


def list_rule_files() -> list[dict]:
    """Return list of {name, download_url} dicts for every TSV in /rules."""
    resp = requests.get(GITHUB_API, timeout=30)
    resp.raise_for_status()
    return [f for f in resp.json() if f["name"].endswith(".tsv")]


def parse_species_from_filename(name: str) -> str:
    """'Escherichia_coli.tsv' -> 'Escherichia coli'"""
    return re.sub(r"\.tsv$", "", name).replace("_", " ")


def fetch_rules(file_meta: dict, species: str) -> pd.DataFrame:
    url = file_meta.get("download_url") or RAW_BASE + file_meta["name"]
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    df = pd.read_csv(io.StringIO(resp.text), sep="\t", dtype=str)
    df["_species"] = species
    return df


def is_intrinsic(df: pd.DataFrame) -> pd.Series:
    """
    Mark a rule row as intrinsic (wildtype susceptible) using the actual
    AMRrules column schema:
      - 'gene context' == 'core'  (gene is intrinsic/chromosomal in this species)
      - OR 'phenotype' == 'wildtype'
      - AND 'clinical category' == 'S'  (susceptible — not conferring resistance)
    We keep rows where the gene is core AND susceptible, i.e. no need to show.
    """
    mask = pd.Series(False, index=df.index)

    # Primary signal: gene context = core (intrinsic chromosomal gene)
    if "gene context" in df.columns:
        mask |= df["gene context"].str.lower().str.strip() == "core"

    # Secondary: wildtype phenotype
    if "phenotype" in df.columns:
        mask |= df["phenotype"].str.lower().str.strip() == "wildtype"

    # Restrict to susceptible clinical category (S) — exclude intermediate/resistant rows
    if "clinical category" in df.columns:
        susceptible = df["clinical category"].str.upper().str.strip().isin(["S", ""])
        mask &= susceptible | df["clinical category"].isna()

    return mask


def main():
    print("Fetching AMRrules file list from GitHub …")
    files = list_rule_files()
    print(f"  Found {len(files)} rule files.")

    rows = []
    for f in files:
        species = parse_species_from_filename(f["name"])
        print(f"  Parsing {f['name']} ({species})")
        try:
            df = fetch_rules(f, species)
        except Exception as e:
            print(f"    WARNING: could not fetch {f['name']}: {e}", file=sys.stderr)
            continue

        intrinsic_mask = is_intrinsic(df)
        intrinsic_df   = df[intrinsic_mask].copy()

        # Normalise gene symbol column name (actual col is 'gene' in AMRrules)
        gene_col = next(
            (c for c in df.columns if c.lower() in ("gene", "gene_symbol", "gene_name")),
            None
        )
        # 'drug class' is the column in AMRrules (with a space)
        drug_col = next(
            (c for c in df.columns if c.lower() in ("drug class", "drug_class", "class", "subclass")),
            None
        )

        if gene_col is None:
            print(f"    WARNING: no gene column found in {f['name']}, skipping.", file=sys.stderr)
            continue

        # Use 'organism' column if present (overrides filename-derived species)
        org_col = next((c for c in df.columns if c.lower() == "organism"), None)

        for _, row in intrinsic_df.iterrows():
            # AMRrules organism field looks like "s__Escherichia coli" — strip prefix
            org_val = str(row[org_col]).lstrip("s__").strip() if org_col else species
            rows.append({
                "species":     org_val,
                "gene_symbol": str(row[gene_col]).strip(),
                "drug_class":  str(row[drug_col]).strip() if drug_col else "",
            })

    if not rows:
        print("No intrinsic rules found — check column names in AMRrules TSVs.")
        sys.exit(1)

    result = pd.DataFrame(rows).drop_duplicates()
    print(f"\nTotal intrinsic gene–species pairs: {len(result)}")
    result.to_parquet(OUT_PATH, index=False)
    print(f"Saved to {OUT_PATH}")


if __name__ == "__main__":
    main()
