#!/usr/bin/env bash
# run_pipeline.sh — end-to-end pipeline runner
# Usage: bash pipeline/run_pipeline.sh [--skip-download]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# Optionally skip slow download steps if data already present
SKIP_DOWNLOAD="${1:-}"

log "=== ATB AMR Taxonium Annotation Pipeline ==="

# 0. Environment
if [ ! -d ".venv" ]; then
    log "Creating virtual environment …"
    python3 -m venv .venv
fi
source .venv/bin/activate
pip install -q -r requirements.txt

# 1. Fetch AMRrules intrinsic gene list
log "Step 1: Fetch AMRrules intrinsic gene definitions"
python pipeline/01_fetch_amrrules.py

# 2. Download AMRFinderPlus results from OSF
if [ "$SKIP_DOWNLOAD" != "--skip-download" ]; then
    log "Step 2: Download AMRFinderPlus results from OSF (this may take a while)"
    python pipeline/02_fetch_amrfp.py
else
    log "Step 2: Skipped (--skip-download)"
fi

# 3. Filter intrinsic genes and aggregate per sample
log "Step 3: Filter intrinsic genes & aggregate per sample"
python pipeline/03_filter_and_aggregate.py

# 4. Fetch ATB metadata (species, subspecies) and merge
log "Step 4: Fetch ATB metadata and merge with AMR results"
python pipeline/04_fetch_atb_metadata.py

log "=== Pipeline complete ==="
log ""
log "Output files:"
log "  data/processed/intrinsic_genes.parquet    — AMRrules intrinsic gene lookup"
log "  data/processed/amr_per_sample.parquet     — per-sample AMR summary"
log "  data/processed/taxonium_metadata.tsv.gz   — ready to upload to taxonium.org"
log ""
log "To visualise:"
log "  1. Go to https://taxonium.org/atb"
log "  2. Click 'Add metadata' and upload data/processed/taxonium_metadata.tsv.gz"
log "     (match column 'sample_id' to node names)"
log "  3. Colour by 'amr_class_summary', 'subspecies', etc."
log ""
log "To self-host with an annotated JSONL (optional):"
log "  python pipeline/05_build_taxonium_jsonl.py \\"
log "      --input <atb.jsonl.gz> \\"
log "      --output data/processed/atb_amr_annotated.jsonl.gz"
log "  docker compose up"
