#!/usr/bin/env python3
"""
GEO Dataset Downloader
───────────────────────
Downloads curated GEO datasets produced by omics_agent.py.

For each GSE accession it fetches:
  1. Series matrix file  (GSE*_series_matrix.txt.gz)  — metadata + processed data
  2. RAW supplementary archive  (GSE*_RAW.tar)        — count matrices, peak files, etc.
  3. Individual supplementary files listed in the GEO FTP index

File layout on disk:
  <out_dir>/
    <accession>/
      <accession>_series_matrix.txt.gz
      <accession>_RAW.tar
      supplementary/
        *.gz  (individual count/peak files)
      download_manifest.json   ← what was downloaded, sizes, checksums

Usage:
    # Download all curated datasets from agent output
    python geo_downloader.py --csv omics_results/omics_datasets_psoriasis_curated.csv

    # Download specific accessions
    python geo_downloader.py --accessions GSE200637 GSE301804 GSE314390

    # Dry-run: show what would be downloaded without fetching
    python geo_downloader.py --csv omics_results/omics_datasets_psoriasis_curated.csv --dry-run

    # Skip RAW.tar (large) and only get series matrix + supplementary files
    python geo_downloader.py --csv ... --no-raw

    # Resume interrupted downloads (skips already-complete files)
    python geo_downloader.py --csv ... --resume

Requirements:
    pip install requests rich
"""

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import requests

try:
    from rich.console import Console
    from rich.progress import (Progress, SpinnerColumn, BarColumn,
                                DownloadColumn, TransferSpeedColumn,
                                TimeRemainingColumn, TextColumn)
    from rich.table import Table
    from rich.panel import Panel
    from rich.markup import escape
    HAS_RICH = True
    console = Console()
except ImportError:
    HAS_RICH = False
    console = None


GEO_FTP_BASE  = "https://ftp.ncbi.nlm.nih.gov/geo/series"
GEO_SOFT_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CHUNK_SIZE    = 1024 * 1024   # 1 MB read chunks
DL_WORKERS    = 3             # parallel dataset downloads (be polite to NCBI)
NCBI_DELAY    = 0.5


# ════════════════════════════════════════════════════════════
#  GEO FTP PATH HELPERS
# ════════════════════════════════════════════════════════════

def geo_ftp_stub(accession: str) -> str:
    """
    GEO FTP uses a stub directory based on the accession prefix.
    GSE200637  → GSE200nnn
    GSE1234567 → GSE1234nnn
    """
    prefix = accession[:-3] + "nnn"   # replace last 3 digits with 'nnn'
    return f"{GEO_FTP_BASE}/{prefix}/{accession}"


def geo_file_urls(accession: str) -> dict[str, str]:
    """
    Return the standard GEO file URLs for a given GSE accession.
    These URLs always exist (NCBI convention) even before we check
    if the files are actually present.
    """
    stub = geo_ftp_stub(accession)
    return {
        "series_matrix": f"{stub}/matrix/{accession}_series_matrix.txt.gz",
        "raw_tar":        f"{stub}/suppl/{accession}_RAW.tar",
        "suppl_index":    f"{stub}/suppl/",
    }


# ════════════════════════════════════════════════════════════
#  SUPPLEMENTARY FILE DISCOVERY
# ════════════════════════════════════════════════════════════

def list_supplementary_files(accession: str) -> list[dict]:
    """
    Fetch the FTP index page for the supplementary directory and
    parse out individual file names, sizes, and URLs.
    Returns list of {name, url, size_bytes}.
    """
    index_url = geo_file_urls(accession)["suppl_index"]
    try:
        r = requests.get(index_url, timeout=20,
                         headers={"User-Agent": "GEODownloader/1.0"})
        if r.status_code == 404:
            return []
        r.raise_for_status()
        html = r.text

        files = []
        # NCBI FTP index pages use <a href="filename"> links
        for match in re.finditer(r'href="([^"]+)".*?(\d+)\s*$',
                                  html, re.MULTILINE):
            name = match.group(1).strip("/")
            size = int(match.group(2))
            # Skip parent directory links and the RAW.tar (handled separately)
            if name in ("../", "/", "") or name.endswith("/"):
                continue
            if name.endswith("_RAW.tar"):
                continue
            files.append({
                "name":       name,
                "url":        f"{index_url}{name}",
                "size_bytes": size,
            })
        return files
    except Exception as e:
        log(f"Could not list supplementary files for {accession}: {e}", "warn")
        return []


# ════════════════════════════════════════════════════════════
#  DOWNLOADER
# ════════════════════════════════════════════════════════════

def log(msg: str, status: str = "info"):
    ts = datetime.now().strftime("%H:%M:%S")
    sym = {"info": "·", "success": "✓", "warn": "!", "error": "✗"}.get(status, "·")
    if HAS_RICH:
        col = {"info": "dim", "success": "green",
               "warn": "yellow", "error": "red"}.get(status, "white")
        console.print(f"[dim]{ts}[/dim]  [{col}]{sym}[/{col}]  {escape(str(msg))}")
    else:
        print(f"{ts}  {sym}  {msg}")


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(CHUNK_SIZE), b""):
            h.update(chunk)
    return h.hexdigest()


def download_file(url: str, dest: Path, resume: bool = True) -> dict:
    """
    Download a single file with optional resume support.
    Returns a result dict: {url, path, size_bytes, md5, status, error}.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    existing_size = dest.stat().st_size if dest.exists() else 0

    headers = {"User-Agent": "GEODownloader/1.0"}
    if resume and existing_size > 0:
        headers["Range"] = f"bytes={existing_size}-"

    try:
        r = requests.get(url, headers=headers, stream=True, timeout=60)

        # 416 = range not satisfiable → file already complete
        if r.status_code == 416:
            return {
                "url": url, "path": str(dest),
                "size_bytes": existing_size,
                "md5": md5_file(dest),
                "status": "already_complete", "error": ""
            }

        if r.status_code == 404:
            return {
                "url": url, "path": str(dest),
                "size_bytes": 0, "md5": "",
                "status": "not_found",
                "error": "404 — file not present on GEO FTP"
            }

        r.raise_for_status()

        total = int(r.headers.get("content-length", 0))
        mode  = "ab" if (resume and existing_size > 0) else "wb"

        with open(dest, mode) as f:
            downloaded = existing_size
            for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)

        final_size = dest.stat().st_size
        return {
            "url": url, "path": str(dest),
            "size_bytes": final_size,
            "md5": md5_file(dest),
            "status": "ok", "error": ""
        }

    except requests.exceptions.RequestException as e:
        return {
            "url": url, "path": str(dest),
            "size_bytes": existing_size,
            "md5": "", "status": "error", "error": str(e)
        }


def fmt_size(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


# ════════════════════════════════════════════════════════════
#  PER-ACCESSION DOWNLOAD ORCHESTRATION
# ════════════════════════════════════════════════════════════

def download_geo_dataset(accession: str, out_dir: Path,
                          include_raw: bool = True,
                          resume: bool = True,
                          dry_run: bool = False) -> dict:
    """
    Download all files for a single GSE accession.
    Returns a manifest dict describing what was (or would be) downloaded.
    """
    acc_dir   = out_dir / accession
    suppl_dir = acc_dir / "supplementary"
    acc_dir.mkdir(parents=True, exist_ok=True)
    suppl_dir.mkdir(parents=True, exist_ok=True)

    urls      = geo_file_urls(accession)
    suppl_files = list_supplementary_files(accession)
    time.sleep(NCBI_DELAY)

    # Build the full download plan
    plan = []

    plan.append({
        "label":  "series_matrix",
        "url":    urls["series_matrix"],
        "dest":   acc_dir / f"{accession}_series_matrix.txt.gz",
        "required": True,
    })

    if include_raw:
        plan.append({
            "label":  "raw_tar",
            "url":    urls["raw_tar"],
            "dest":   acc_dir / f"{accession}_RAW.tar",
            "required": False,   # not all series have a RAW.tar
        })

    for sf in suppl_files:
        name = sf["name"]
        # Skip the series matrix (already in plan) and avoid duplication
        if "series_matrix" in name.lower():
            continue
        plan.append({
            "label":  f"suppl/{name}",
            "url":    sf["url"],
            "dest":   suppl_dir / name,
            "required": False,
        })

    if dry_run:
        log(f"[DRY RUN] {accession}: {len(plan)} files planned")
        for item in plan:
            log(f"  {item['label']:40s}  {item['url']}", "info")
        return {
            "accession": accession,
            "status":    "dry_run",
            "files":     [{"label": p["label"], "url": p["url"]} for p in plan],
        }

    # ── Execute downloads ──
    log(f"Downloading {accession} ({len(plan)} files)...")
    manifest_files = []
    errors = []

    for item in plan:
        dest: Path = item["dest"]

        # Skip if already complete and resuming
        if resume and dest.exists() and dest.stat().st_size > 0:
            # Quick check: if it's a gz, make sure it's not truncated
            if str(dest).endswith(".gz"):
                try:
                    with gzip.open(dest, "rb") as gf:
                        gf.read(128)
                    log(f"  {item['label']}: already complete, skipping", "info")
                    manifest_files.append({
                        "label":      item["label"],
                        "url":        item["url"],
                        "path":       str(dest),
                        "size_bytes": dest.stat().st_size,
                        "status":     "skipped_complete",
                    })
                    continue
                except Exception:
                    log(f"  {item['label']}: existing file corrupt, re-downloading", "warn")
            else:
                log(f"  {item['label']}: already complete, skipping", "info")
                manifest_files.append({
                    "label":      item["label"],
                    "url":        item["url"],
                    "path":       str(dest),
                    "size_bytes": dest.stat().st_size,
                    "status":     "skipped_complete",
                })
                continue

        result = download_file(item["url"], dest, resume=resume)

        if result["status"] == "not_found" and not item["required"]:
            log(f"  {item['label']}: not present on FTP (optional, skipping)", "info")
        elif result["status"] in ("ok", "already_complete"):
            size_str = fmt_size(result["size_bytes"])
            log(f"  {item['label']}: {size_str}  [{result['status']}]", "success")
        else:
            log(f"  {item['label']}: {result['error']}", "warn")
            errors.append(result["error"])

        manifest_files.append({
            "label":      item["label"],
            "url":        result["url"],
            "path":       result["path"],
            "size_bytes": result["size_bytes"],
            "md5":        result.get("md5", ""),
            "status":     result["status"],
            "error":      result.get("error", ""),
        })

    total_bytes = sum(f["size_bytes"] for f in manifest_files)
    manifest = {
        "accession":    accession,
        "downloaded_at": datetime.now().isoformat(),
        "total_size":   fmt_size(total_bytes),
        "total_bytes":  total_bytes,
        "file_count":   len(manifest_files),
        "error_count":  len(errors),
        "files":        manifest_files,
    }

    manifest_path = acc_dir / "download_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    status = "complete" if not errors else "partial"
    log(f"{accession}: {status} — {len(manifest_files)} files, {fmt_size(total_bytes)}",
        "success" if status == "complete" else "warn")
    return manifest


# ════════════════════════════════════════════════════════════
#  BATCH RUNNER
# ════════════════════════════════════════════════════════════

def run_downloads(accessions: list[str], out_dir: Path,
                  include_raw: bool = True, resume: bool = True,
                  dry_run: bool = False, workers: int = DL_WORKERS):

    out_dir.mkdir(parents=True, exist_ok=True)
    log(f"Downloading {len(accessions)} GEO datasets → {out_dir}")
    if dry_run:
        log("DRY RUN — no files will be written", "warn")

    manifests = {}

    # Use parallel workers for independent dataset downloads
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(
                download_geo_dataset,
                acc, out_dir, include_raw, resume, dry_run
            ): acc
            for acc in accessions
        }
        for future in as_completed(futures):
            acc = futures[future]
            try:
                manifests[acc] = future.result()
            except Exception as e:
                log(f"{acc}: unexpected error — {e}", "error")
                manifests[acc] = {
                    "accession": acc, "status": "error", "error": str(e)
                }

    # ── Summary ──
    total_bytes = sum(
        m.get("total_bytes", 0) for m in manifests.values()
    )
    complete = sum(1 for m in manifests.values()
                   if m.get("error_count", 0) == 0 and m.get("status") != "error")
    errors   = len(accessions) - complete

    if HAS_RICH:
        t = Table(title="Download Summary", box=None,
                  header_style="bold", padding=(0, 2))
        t.add_column("Accession", style="cyan", no_wrap=True)
        t.add_column("Files", justify="right")
        t.add_column("Size")
        t.add_column("Status")
        for acc in accessions:
            m = manifests.get(acc, {})
            n_files   = m.get("file_count", 0)
            size_str  = m.get("total_size", "—")
            err_count = m.get("error_count", 0)
            status    = "✓ complete" if err_count == 0 else f"! {err_count} errors"
            color     = "green" if err_count == 0 else "yellow"
            t.add_row(acc, str(n_files), size_str,
                      f"[{color}]{status}[/{color}]")
        console.print()
        console.print(t)
        console.print(Panel(
            f"[green]✓[/green]  {complete}/{len(accessions)} complete\n"
            f"[dim]Total size: {fmt_size(total_bytes)}[/dim]\n"
            f"[dim]Output:     {out_dir}[/dim]",
            title="Done", border_style="dim"
        ))
    else:
        print(f"\nDownload complete: {complete}/{len(accessions)}")
        print(f"Total size: {fmt_size(total_bytes)}")
        print(f"Output: {out_dir}")

    # Write master manifest
    master = {
        "downloaded_at": datetime.now().isoformat(),
        "accessions":    accessions,
        "total_bytes":   total_bytes,
        "total_size":    fmt_size(total_bytes),
        "datasets":      manifests,
    }
    master_path = out_dir / "master_manifest.json"
    with open(master_path, "w") as f:
        json.dump(master, f, indent=2)
    log(f"Master manifest: {master_path}", "info")

    return master


# ════════════════════════════════════════════════════════════
#  CSV READER
# ════════════════════════════════════════════════════════════

def accessions_from_csv(csv_path: Path) -> list[str]:
    """Read GSE accessions from a curated CSV produced by omics_agent.py."""
    accs = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc = row.get("accession", "").strip()
            if acc.startswith("GSE"):
                accs.append(acc)
    return list(dict.fromkeys(accs))   # deduplicate preserving order


# ════════════════════════════════════════════════════════════
#  CLI
# ════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Download GEO datasets curated by omics_agent.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download all curated datasets from agent output
  python geo_downloader.py --csv omics_results/omics_datasets_psoriasis_curated.csv

  # Download specific accessions
  python geo_downloader.py --accessions GSE200637 GSE301804 GSE314390

  # Dry run — see what would be downloaded
  python geo_downloader.py --csv omics_results/omics_datasets_psoriasis_curated.csv --dry-run

  # Skip RAW.tar (large files) — only get series matrix and supplementary files
  python geo_downloader.py --csv omics_results/... --no-raw

  # Resume interrupted download
  python geo_downloader.py --csv omics_results/... --resume
        """
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--csv",         "-c", type=str,
                     help="Path to curated CSV from omics_agent.py")
    src.add_argument("--accessions",  "-a", nargs="+",
                     help="One or more GSE accession numbers")

    parser.add_argument("--out",      "-o", type=str, default="geo_downloads",
                        help="Output directory (default: geo_downloads)")
    parser.add_argument("--no-raw",   action="store_true",
                        help="Skip RAW.tar download (saves disk space)")
    parser.add_argument("--resume",   action="store_true", default=True,
                        help="Resume interrupted downloads (default: on)")
    parser.add_argument("--no-resume",action="store_true",
                        help="Force re-download even if file exists")
    parser.add_argument("--dry-run",  action="store_true",
                        help="Show download plan without fetching")
    parser.add_argument("--workers",  "-w", type=int, default=DL_WORKERS,
                        help=f"Parallel download workers (default: {DL_WORKERS})")
    args = parser.parse_args()

    # Resolve accessions
    if args.csv:
        csv_path = Path(args.csv)
        if not csv_path.exists():
            print(f"CSV not found: {csv_path}")
            sys.exit(1)
        accessions = accessions_from_csv(csv_path)
        if not accessions:
            print("No GSE accessions found in CSV.")
            sys.exit(1)
        log(f"Found {len(accessions)} GSE accessions in {csv_path.name}")
    else:
        accessions = [a.strip() for a in args.accessions if a.strip().startswith("GSE")]
        if not accessions:
            print("No valid GSE accessions provided (must start with 'GSE').")
            sys.exit(1)

    resume = args.resume and not args.no_resume

    try:
        run_downloads(
            accessions    = accessions,
            out_dir       = Path(args.out),
            include_raw   = not args.no_raw,
            resume        = resume,
            dry_run       = args.dry_run,
            workers       = args.workers,
        )
    except KeyboardInterrupt:
        print("\nInterrupted. Partial files are kept — re-run with --resume to continue.")
    except Exception as e:
        print(f"Error: {e}")
        raise


if __name__ == "__main__":
    main()
