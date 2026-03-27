# AllTheBacteria tree — AMR annotation

Annotates the [AllTheBacteria (ATB)](https://allthebacteria.readthedocs.io/) phylogenetic tree with [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) results, filtered using [AMRrules](https://github.com/AMRverse/AMRrules) to exclude intrinsic wildtype-susceptible genes.

## Explore the data

> **[🌐 Landing page + interactive charts](https://amr-genomics-hackathon-2026.github.io/atb_tree/)**

> **[🔬 Interactive Streamlit explorer](https://atb-amr-explorer.streamlit.app)** — filter by species, drug class, gene; download subsets

## View in Taxonium

> **[▶ Open ATB tree with AMR annotation](https://taxonium.org/atb?addMetadataUrl=https://github.com/AMR-genomics-hackathon-2026/atb_tree/releases/latest/download/taxonium_metadata.tsv.gz)**

Colour by:
- `Species` — GTDB species assignment
- `Genus` / `Family` — higher-level taxonomy
- `amr_class_summary` — drug classes with acquired resistance genes
- `amr_present` — any acquired AMR gene (true/false)
- `amr_gene_count` — number of distinct acquired genes

A binary presence-absence matrix version (one column per gene, top 500 genes) is also available:

> **[▶ Open with binary gene matrix](https://taxonium.org/atb?addMetadataUrl=https://github.com/AMR-genomics-hackathon-2026/atb_tree/releases/latest/download/taxonium_metadata_binary.tsv.gz)**

## Data sources

| Resource | URL |
|----------|-----|
| AMRFinderPlus results (AllTheBacteria) | https://osf.io/ck7st |
| ATB file list (latest) | https://osf.io/download/4yv85/ |
| AMRrules intrinsic gene rules | https://github.com/AMRverse/AMRrules/tree/main/rules |
| GTDB taxonomy + tree | bundled in `atb_tree.tar.bz2` |

## Intrinsic gene filtering

Acquired AMR genes are identified using a two-layer filter:

1. **AMRFinderPlus `Scope=core`** — genes flagged as chromosomally intrinsic by AMRFinderPlus for that taxon
2. **AMRrules species-aware lookup** — genes appearing in the AMRrules rules TSVs with `gene context=core`, `phenotype=wildtype`, or `clinical category=S`

Only `Scope=plus` (acquired) genes passing both filters are retained.

## Output files (GitHub Release)

| File | Size | Description |
|------|------|-------------|
| `taxonium_metadata.tsv.gz` | ~13 MB | Summary: one row/sample, AMR class summary + GTDB taxonomy |
| `taxonium_metadata_binary.tsv.gz` | ~13 MB | Binary matrix: top 500 genes, 100 classes, 200 subclasses as 0/1 columns |

## Reproduce the pipeline

```bash
git clone https://github.com/AMR-genomics-hackathon-2026/atb_tree
cd atb_tree

python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# Full pipeline (downloads ~700 MB of AMRFinderPlus data)
bash run_pipeline.sh

# Or step by step:
python 01_fetch_amrrules.py              # AMRrules intrinsic gene list
python 02_fetch_amrfp.py                 # Download AMRFinderPlus results from OSF
python 03_filter_and_aggregate.py        # Filter intrinsic genes, aggregate per sample
python 04b_build_taxonium_metadata.py    # Merge with GTDB taxonomy, build outputs
```

## Files in this repo

| File | Description |
|------|-------------|
| `01_fetch_amrrules.py` | Download AMRrules intrinsic gene definitions |
| `02_fetch_amrfp.py` | Download AMRFinderPlus results from OSF |
| `03_filter_and_aggregate.py` | Filter intrinsic genes; aggregate per sample |
| `04b_build_taxonium_metadata.py` | Build Taxonium metadata using GTDB taxonomy |
| `05_build_taxonium_jsonl.py` | (Optional) Embed metadata into ATB JSONL for self-hosting |
| `AMRFP_to_treemeta.R` | R script: binary presence-absence matrices; updated to split multi-class entries and reduce metadata to Family/Genus/Species |
| `ReferenceGeneHierarchy.txt` | AMRFinderPlus gene hierarchy reference |
| `atb_tree.tar.bz2` | ATB Newick tree + GTDB taxonomy assignments |
| `config/taxonium_config.json` | Taxonium display configuration |

## Citation

If you use this resource, please cite:
- AllTheBacteria: [Blackwell et al., 2024](https://www.biorxiv.org/content/10.1101/2024.03.08.584059)
- AMRFinderPlus: [Feldgarden et al., 2021](https://doi.org/10.1038/s41598-021-91456-0)
- AMRrules: [Mather et al., AMRverse](https://github.com/AMRverse/AMRrules)
- Taxonium: [Sanderson, 2022](https://elifesciences.org/articles/82392)
