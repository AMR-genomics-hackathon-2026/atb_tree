#!/usr/bin/env python3
"""
02_fetch_amrfp.py
Download AMRFinderPlus results from the AllTheBacteria OSF project.

OSF project: https://osf.io/7nwrx  (AMRFinderPlus v3.12.8, DB 2024-01-31)
Structure discovered:
  latest/
    AMRFP_results.tsv.gz   (691 MB, guid=ck7st)  ← full combined results
    AMRFP_status.tsv.gz    ( 11 MB, guid=t2dkf)  ← per-sample run QC
  Incremental_release_2024-08/
    AMRFP_results.tsv.gz   (149 MB, guid=rd9j3)  ← incremental batch
    AMRFP_status.tsv.gz    (  2 MB, guid=97qav)

By default, only the 'latest/' files are downloaded (results + status).
Use --release to also grab the incremental batch.

Outputs land in data/raw/amrfp/
"""

import argparse
import sys
import tarfile
from pathlib import Path

import requests
from tqdm import tqdm

OSF_DOWNLOAD = "https://osf.io/download/{guid}/"

# Known files — keyed by (release, filename): (guid, size_bytes)
OSF_FILES = {
    ("latest",      "AMRFP_results.tsv.gz"):  ("ck7st",  690_757_436),
    ("latest",      "AMRFP_status.tsv.gz"):   ("t2dkf",  11_100_000),
    ("incremental", "AMRFP_results.tsv.gz"):  ("rd9j3",  148_800_000),
    ("incremental", "AMRFP_status.tsv.gz"):   ("97qav",   2_200_000),
}

RAW_DIR    = Path(__file__).parent.parent / "data" / "raw" / "amrfp"
CHUNK_SIZE = 2 * 1024 * 1024  # 2 MB


def list_osf_files(api_url: str) -> list[dict]:
    """Recursively list all files in an OSF storage node (for --list-only)."""
    files = []
    url = api_url
    while url:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        for item in data.get("data", []):
            a = item["attributes"]
            if a["kind"] == "file":
                files.append({
                    "name":   a["name"],
                    "size":   a.get("size") or 0,
                    "guid":   a.get("guid") or "",
                    "url":    item["links"].get("download", ""),
                })
            elif a["kind"] == "folder":
                sub = item["relationships"]["files"]["links"]["related"]["href"]
                files.extend(list_osf_files(sub))
        url = data.get("links", {}).get("next")
    return files


def print_file_listing(files: list[dict]) -> None:
    total = sum(f["size"] for f in files)
    print(f"\n  {'Name':<45} {'Size':>10}  Guid")
    print("  " + "-" * 70)
    for f in files:
        print(f"  {f['name']:<45} {f['size']/1e6:>7.1f} MB  {f['guid']}")
    print("  " + "-" * 70)
    print(f"  {len(files)} file(s)   total: {total/1e9:.2f} GB")


def download_file(guid: str, dest: Path, expected_bytes: int = 0) -> None:
    if dest.exists():
        print(f"  Already present: {dest.name}")
        return
    url = OSF_DOWNLOAD.format(guid=guid)
    print(f"  Downloading {dest.name} ({expected_bytes/1e6:.0f} MB) …")
    resp = requests.get(url, stream=True, timeout=300)
    resp.raise_for_status()
    total = int(resp.headers.get("content-length", expected_bytes))
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as f, tqdm(
        total=total, unit="B", unit_scale=True, desc=dest.name, leave=False
    ) as bar:
        for chunk in resp.iter_content(chunk_size=CHUNK_SIZE):
            f.write(chunk)
            bar.update(len(chunk))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--list-only", action="store_true",
        help="Query the OSF API and list all files without downloading."
    )
    parser.add_argument(
        "--release", choices=["latest", "incremental", "all"], default="latest",
        help="Which release to download (default: latest)."
    )
    parser.add_argument(
        "--results-only", action="store_true",
        help="Skip AMRFP_status files; download only AMRFP_results."
    )
    args = parser.parse_args()

    if args.list_only:
        print("Querying OSF project 7nwrx …")
        OSF_API = "https://api.osf.io/v2/nodes/7nwrx/files/osfstorage/"
        try:
            files = list_osf_files(OSF_API)
        except requests.HTTPError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            sys.exit(1)
        print(f"Found {len(files)} file(s):")
        print_file_listing(files)
        return

    # Decide which releases to grab
    releases = {"latest", "incremental"} if args.release == "all" else {args.release}

    to_download = [
        (release, fname, guid, size)
        for (release, fname), (guid, size) in OSF_FILES.items()
        if release in releases
        and not (args.results_only and "status" in fname.lower())
    ]

    if not to_download:
        print("Nothing to download.")
        return

    print(f"Files to download ({args.release} release):")
    total = sum(s for _, _, _, s in to_download)
    for release, fname, guid, size in to_download:
        print(f"  [{release}] {fname:<35} {size/1e6:>7.1f} MB  (guid={guid})")
    print(f"  Total: {total/1e9:.2f} GB\n")

    for release, fname, guid, size in to_download:
        dest = RAW_DIR / fname
        # Disambiguate if downloading multiple releases
        if args.release == "all":
            dest = RAW_DIR / release / fname
        download_file(guid, dest, size)

    print(f"\nDone. Data is in {RAW_DIR}")
    print("Next: python pipeline/03_filter_and_aggregate.py")


if __name__ == "__main__":
    main()
