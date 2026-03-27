#!/usr/bin/env python3
"""
04_fetch_atb_metadata.py
Build the final Taxonium metadata TSV by:

  1. Downloading the ATB assembly manifest (sample_accession, dataset, filter status)
     from OSF (guid=sdx5m, xz-compressed TSV).
  2. Batch-fetching species / taxonomy / country from the ENA portal API
     (up to 500 accessions per request). Results are cached in
     data/processed/ena_metadata_cache.parquet so interrupted runs resume.
  3. Merging with the AMRFinderPlus per-sample table and AMRFP status file.
  4. Writing data/processed/taxonium_metadata.tsv.gz

Output columns (first column = sample_id, used by Taxonium to match node names):
  sample_id, run_accession, assembly_accession, dataset, atb_filter,
  scientific_name, tax_id, country, collection_date,
  amr_gene_count, amr_genes, amr_class_summary, amr_subclass_summary,
  amr_present, amr_highest_identity, amrfp_status
"""

import io
import sys
import time
import lzma
from pathlib import Path

import requests
import pandas as pd
from tqdm import tqdm

# ── Input paths ────────────────────────────────────────────────────────────────
AMR_TABLE   = Path(__file__).parent.parent / "data" / "processed" / "amr_per_sample.parquet"
STATUS_FILE = Path(__file__).parent.parent / "data" / "raw" / "amrfp" / "AMRFP_status.tsv.gz"
ENA_CACHE       = Path(__file__).parent.parent / "data" / "processed" / "ena_metadata_cache.parquet"
MANIFEST_CACHE  = Path(__file__).parent.parent / "data" / "processed" / "atb_manifest.parquet"
OUT_PATH    = Path(__file__).parent.parent / "data" / "processed" / "taxonium_metadata.tsv.gz"
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

# ── ATB assembly manifest (OSF) ────────────────────────────────────────────────
ATB_MANIFEST_URL = "https://osf.io/download/sdx5m/"   # guid=sdx5m, xz TSV

# ── ENA portal API ─────────────────────────────────────────────────────────────
ENA_API      = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_FIELDS   = "sample_accession,scientific_name,tax_id,country,collection_date"
ENA_BATCH    = 200     # accessions per request (500 causes 400 URL-too-long errors)
ENA_SLEEP    = 0.25    # seconds between requests (be polite)
ENA_SAVE_EVERY = 50    # save cache every N batches (~10k samples)
CHUNK_SIZE   = 2 * 1024 * 1024


# ── Step 1: ATB assembly manifest ─────────────────────────────────────────────

def fetch_atb_manifest() -> pd.DataFrame:
    if MANIFEST_CACHE.exists():
        print(f"Loading ATB manifest from cache …")
        df = pd.read_parquet(MANIFEST_CACHE)
        print(f"  {len(df):,} rows (cached).")
        return df

    print(f"Downloading ATB assembly manifest (xz) …")
    print(f"  {ATB_MANIFEST_URL}")
    resp = requests.get(ATB_MANIFEST_URL, stream=True, timeout=300)
    resp.raise_for_status()
    total = int(resp.headers.get("content-length", 0))
    buf   = io.BytesIO()
    with tqdm(total=total, unit="B", unit_scale=True, desc="manifest", leave=False) as bar:
        for chunk in resp.iter_content(chunk_size=CHUNK_SIZE):
            buf.write(chunk)
            bar.update(len(chunk))
    buf.seek(0)
    with lzma.open(buf, "rt") as f:
        df = pd.read_csv(f, sep="\t", dtype=str, low_memory=False)
    print(f"  {len(df):,} rows. Columns: {list(df.columns)}")
    MANIFEST_CACHE.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(MANIFEST_CACHE, index=False)
    print(f"  Manifest cached to {MANIFEST_CACHE}")
    return df


def normalise_manifest(df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns to standard names."""
    renames = {}
    for col in df.columns:
        low = col.lower()
        if low == "sample_accession":
            renames[col] = "sample_id"
        elif low in ("filter", "atb_filter"):
            renames[col] = "atb_filter"
    df = df.rename(columns=renames)
    # Keep only useful columns
    keep = [c for c in ("sample_id", "run_accession", "assembly_accession",
                         "dataset", "atb_filter", "comments") if c in df.columns]
    return df[keep]


# ── Step 2: ENA species metadata ───────────────────────────────────────────────

def load_ena_cache() -> pd.DataFrame:
    if ENA_CACHE.exists():
        return pd.read_parquet(ENA_CACHE)
    return pd.DataFrame(columns=["sample_id", "scientific_name", "tax_id",
                                  "country", "collection_date"])


def save_ena_cache(df: pd.DataFrame) -> None:
    df.to_parquet(ENA_CACHE, index=False)


def fetch_ena_batch(accessions: list[str]) -> pd.DataFrame:
    """Query ENA search API via POST for a batch of sample accessions."""
    query = " OR ".join(f'sample_accession="{a}"' for a in accessions)
    payload = {
        "result": "sample",
        "query":  query,
        "fields": ENA_FIELDS,
        "format": "tsv",
        "limit":  len(accessions),
    }
    try:
        resp = requests.post(ENA_API, data=payload, timeout=60)
        resp.raise_for_status()
        text = resp.text.strip()
        if not text or len(text.splitlines()) < 2:
            return pd.DataFrame()
        df = pd.read_csv(io.StringIO(text), sep="\t", dtype=str)
        if "sample_accession" in df.columns:
            df = df.rename(columns={"sample_accession": "sample_id"})
        return df
    except Exception as e:
        print(f"    WARNING: ENA batch failed: {e}", file=sys.stderr)
        return pd.DataFrame()


def fetch_ena_metadata(all_accessions: list[str]) -> pd.DataFrame:
    """Fetch species/taxonomy for all accessions, using cache for already-fetched ones."""
    cache = load_ena_cache()
    already = set(cache["sample_id"].tolist())
    to_fetch = [a for a in all_accessions if a not in already]

    if not to_fetch:
        print(f"  ENA metadata: all {len(all_accessions):,} samples in cache.")
        return cache

    print(f"  ENA metadata: {len(already):,} cached, fetching {len(to_fetch):,} new …")
    new_rows = []
    batches  = [to_fetch[i:i+ENA_BATCH] for i in range(0, len(to_fetch), ENA_BATCH)]

    for batch in tqdm(batches, desc="ENA batches"):
        df = fetch_ena_batch(batch)
        if not df.empty:
            new_rows.append(df)
        time.sleep(ENA_SLEEP)

        # Save cache periodically
        if len(new_rows) % ENA_SAVE_EVERY == 0 and new_rows:
            interim = pd.concat([cache] + new_rows, ignore_index=True)
            save_ena_cache(interim)

    if new_rows:
        result = pd.concat([cache] + new_rows, ignore_index=True).drop_duplicates("sample_id")
    else:
        result = cache

    save_ena_cache(result)
    print(f"  ENA cache now has {len(result):,} samples.")
    return result


# ── Step 3: extract subspecies / lineage from scientific name ──────────────────

def add_subspecies(df: pd.DataFrame) -> pd.DataFrame:
    """
    Derive a 'subspecies' column from 'scientific_name'.
    E.g. 'Escherichia coli O157:H7' → subspecies='O157:H7'
         'Klebsiella pneumoniae subsp. pneumoniae' → subspecies='subsp. pneumoniae'
    """
    if "scientific_name" not in df.columns:
        df["subspecies"] = ""
        return df

    def _sub(name: str) -> str:
        parts = str(name).split()
        if len(parts) <= 2:
            return ""
        return " ".join(parts[2:])

    df["subspecies"] = df["scientific_name"].apply(_sub)
    # Shorten genus+species to just 'species' column
    def _sp(name: str) -> str:
        return " ".join(str(name).split()[:2])
    df["species"] = df["scientific_name"].apply(_sp)
    return df


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--amr-only", action="store_true",
                        help="Only fetch ENA metadata for samples with acquired AMR genes (faster).")
    args = parser.parse_args()

    # 1. ATB assembly manifest
    manifest = fetch_atb_manifest()
    manifest = normalise_manifest(manifest)

    accessions = manifest["sample_id"].dropna().unique().tolist()
    print(f"  {len(accessions):,} unique sample accessions in manifest.")

    # Optional: prioritise / restrict to AMR-positive samples
    if args.amr_only and AMR_TABLE.exists():
        amr_ids = set(pd.read_parquet(AMR_TABLE)["sample_id"].tolist())
        accessions = [a for a in accessions if a in amr_ids]
        print(f"  --amr-only: restricting ENA fetch to {len(accessions):,} AMR-positive samples.")

    # 2. ENA species metadata
    print("Fetching ENA species metadata …")
    ena = fetch_ena_metadata(accessions)
    ena = add_subspecies(ena)

    # 3. Merge manifest + ENA
    meta = manifest.merge(ena, on="sample_id", how="left")

    # 4. AMRFP run status
    if STATUS_FILE.exists():
        st = pd.read_csv(STATUS_FILE, sep="\t", dtype=str)
        st = st.rename(columns={"sample": "sample_id"})
        meta = meta.merge(st[["sample_id", "status"]].rename(columns={"status": "amrfp_status"}),
                          on="sample_id", how="left")
        meta["amrfp_status"] = meta["amrfp_status"].fillna("not_run")
        print(f"  AMRFinderPlus run status merged ({st['status'].value_counts().to_dict()}).")

    # 5. AMR per-sample table
    if AMR_TABLE.exists():
        amr = pd.read_parquet(AMR_TABLE)
        meta = meta.merge(amr, on="sample_id", how="left")
        meta["amr_gene_count"]       = meta["amr_gene_count"].fillna(0).astype(int)
        meta["amr_genes"]            = meta["amr_genes"].fillna("")
        meta["amr_class_summary"]    = meta["amr_class_summary"].fillna("")
        meta["amr_subclass_summary"] = meta.get("amr_subclass_summary",
                                                  pd.Series("", index=meta.index)).fillna("")
        meta["amr_present"]          = meta["amr_present"].fillna(False)
        print(f"  AMR data merged: {(meta['amr_present']==True).sum():,} / {len(meta):,} samples have acquired AMR genes.")
    else:
        print("WARNING: amr_per_sample.parquet not found. Run steps 02 and 03 first.")

    # 6. Save
    meta.to_csv(OUT_PATH, sep="\t", index=False)
    print(f"\nSaved Taxonium metadata: {OUT_PATH}")
    print(f"  Rows: {len(meta):,}   Columns: {list(meta.columns)}")


if __name__ == "__main__":
    main()
