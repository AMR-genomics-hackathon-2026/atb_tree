#!/usr/bin/env python3
"""
05_build_taxonium_jsonl.py
Merge the Taxonium metadata TSV into the ATB Taxonium JSONL tree to produce
an annotated tree ready for upload/hosting.

Two workflows are supported:

A) ONLINE TAXONIUM (upload metadata only)
   If you are using taxonium.org and just need a metadata TSV to upload
   alongside the ATB tree, the output of 04_fetch_atb_metadata.py is already
   sufficient — no JSONL reprocessing needed. Skip to step B.

B) SELF-HOSTED / JSONL REPROCESSING
   If you have downloaded the ATB Taxonium JSONL (e.g., from S3) and want to
   produce a new JSONL with AMR metadata embedded, this script handles that.

Usage:
    python 05_build_taxonium_jsonl.py \
        --input  <path_to_atb.jsonl.gz>  \
        --meta   data/processed/taxonium_metadata.tsv.gz \
        --config config/taxonium_config.json \
        --output data/processed/atb_amr_annotated.jsonl.gz

The output file can be:
  - Loaded directly into taxonium.org (drag-and-drop)
  - Hosted locally via the Taxonium backend Docker container
"""

import argparse
import gzip
import json
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def load_metadata(meta_path: str) -> dict[str, dict]:
    """Load metadata TSV into a dict keyed by sample_id."""
    df = pd.read_csv(meta_path, sep="\t", dtype=str, low_memory=False,
                     compression="infer")
    id_col = df.columns[0]
    meta   = {}
    for _, row in df.iterrows():
        key = str(row[id_col]).strip()
        meta[key] = {k: ("" if pd.isna(v) else str(v)) for k, v in row.items() if k != id_col}
    return meta


def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return json.load(f)


def open_jsonl(path: str):
    """Open a plain or gzipped JSONL file for reading line by line."""
    p = Path(path)
    if p.suffix == ".gz" or str(path).endswith(".jsonl.gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def annotate_jsonl(input_path: str, meta: dict, config: dict, output_path: str) -> None:
    out_path = Path(output_path)
    open_fn  = gzip.open if output_path.endswith(".gz") else open

    # Count lines for progress bar (requires a full pass — skip for huge files)
    n_lines = None
    try:
        with open_jsonl(input_path) as fh:
            n_lines = sum(1 for _ in fh)
    except Exception:
        pass

    matched = 0
    total   = 0

    with open_jsonl(input_path) as fin, open_fn(out_path, "wt", encoding="utf-8") as fout:
        for i, line in enumerate(tqdm(fin, total=n_lines, desc="Annotating nodes")):
            line = line.strip()
            if not line:
                continue

            obj = json.loads(line)

            # First line is tree-wide metadata — merge our config into it
            if i == 0:
                obj.setdefault("config", {})
                obj["config"].update(config)
                fout.write(json.dumps(obj) + "\n")
                continue

            # Remaining lines are nodes — look up metadata by name
            node_name = obj.get("name") or obj.get("node_id") or obj.get("id", "")
            if node_name and node_name in meta:
                obj.update(meta[node_name])
                matched += 1
            total += 1

            fout.write(json.dumps(obj) + "\n")

    print(f"Annotated {matched:,} / {total:,} nodes with metadata.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",  required=True, help="ATB Taxonium JSONL (plain or .gz)")
    parser.add_argument("--meta",   default="data/processed/taxonium_metadata.tsv.gz",
                        help="Metadata TSV produced by step 04")
    parser.add_argument("--config", default="config/taxonium_config.json",
                        help="Taxonium config JSON")
    parser.add_argument("--output", default="data/processed/atb_amr_annotated.jsonl.gz",
                        help="Output annotated JSONL (.gz)")
    args = parser.parse_args()

    print("Loading metadata …")
    meta = load_metadata(args.meta)
    print(f"  {len(meta):,} samples in metadata table.")

    print("Loading Taxonium config …")
    config = load_config(args.config)

    print("Annotating JSONL …")
    annotate_jsonl(args.input, meta, config, args.output)
    print(f"Done. Output: {args.output}")


if __name__ == "__main__":
    main()
