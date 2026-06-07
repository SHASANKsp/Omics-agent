#!/usr/bin/env python3
"""
sra_downloader.py — Download SRA runs via Docker
──────────────────────────────────────────────────
Reads sra_manifest.json produced by geo_downloader.py.
For each SRR run: prefetch → fasterq-dump → gzip → update manifest status.
Runs one SRR at a time so peak disk usage = one sample's worth of FASTQ.

Usage:
    python processing/sra_downloader.py --manifest geo_downloads/sra_manifest.json
    python processing/sra_downloader.py --manifest ... --dry-run
    python processing/sra_downloader.py --manifest ... --srr SRR12345 SRR12346

Requirements:
    Docker running with the omics-pipeline image built.
    pip install requests rich pyyaml
"""

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import yaml

try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    from rich.markup import escape
    HAS_RICH = True
    console = Console()
except ImportError:
    HAS_RICH = False
    console = None

CONFIG_FILE = Path("pipeline_config.yaml")


# ════════════════════════════════════════════════════════════
#  LOGGING
# ════════════════════════════════════════════════════════════

def log(msg: str, status: str = "info"):
    ts  = datetime.now().strftime("%H:%M:%S")
    sym = {"info":"·","success":"✓","warn":"!","error":"✗"}.get(status,"·")
    if HAS_RICH:
        col = {"info":"dim","success":"green","warn":"yellow","error":"red"}.get(status,"white")
        console.print(f"[dim]{ts}[/dim]  [{col}]{sym}[/{col}]  {escape(str(msg))}")
    else:
        print(f"{ts}  {sym}  {msg}")


def fmt_size(mb: float) -> str:
    if mb < 1024: return f"{mb:.0f} MB"
    return f"{mb/1024:.1f} GB"


# ════════════════════════════════════════════════════════════
#  CONFIG
# ════════════════════════════════════════════════════════════

def load_config() -> dict:
    if not CONFIG_FILE.exists():
        log(f"Config not found: {CONFIG_FILE}", "error"); sys.exit(1)
    with open(CONFIG_FILE) as f:
        return yaml.safe_load(f)


# ════════════════════════════════════════════════════════════
#  DOCKER CHECK
# ════════════════════════════════════════════════════════════

def check_docker(cfg: dict) -> bool:
    image = cfg["docker"]["image"]
    r = subprocess.run(f"docker image inspect {image}",
                       shell=True, capture_output=True)
    if r.returncode != 0:
        log(f"Docker image '{image}' not found.", "error")
        log("Build it:  docker build -t omics-pipeline:1.0 processing/", "warn")
        return False
    log(f"Docker image '{image}' found", "success")
    return True


# ════════════════════════════════════════════════════════════
#  DOCKER RUNNER
# ════════════════════════════════════════════════════════════

def docker_run(cfg: dict, command: str,
               extra_mounts: list[str] = None,
               workdir: str = "/data") -> tuple[int, str, str]:
    """
    Run a shell command inside the omics Docker container.
    Returns (returncode, stdout, stderr).
    """
    image    = cfg["docker"]["image"]
    mounts   = []
    if extra_mounts:
        mounts.extend(extra_mounts)

    cmd = (
        f"docker run --rm "
        f"{' '.join(mounts)} "
        f"-w {workdir} "
        f"{image} "
        f"bash -c \"{command}\""
    )
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.returncode, result.stdout, result.stderr


# ════════════════════════════════════════════════════════════
#  PER-SRR DOWNLOAD
# ════════════════════════════════════════════════════════════

def download_srr(srr: str, gse: str, gsm: str, omics_type: str,
                 cfg: dict, fastq_tmp: Path,
                 dry_run: bool = False) -> dict:
    """
    Download one SRR accession:
      1. prefetch  → .sra cache file
      2. fasterq-dump → .fastq files
      3. pigz → .fastq.gz
      4. delete uncompressed .fastq and .sra cache

    FASTQ files land in: fastq_tmp/{gse}/{gsm}/{srr}/
    Returns status dict that gets written back to sra_manifest.json.
    """
    out_dir = fastq_tmp / gse / gsm / srr
    out_dir.mkdir(parents=True, exist_ok=True)

    # Paths as seen inside the container
    host_out   = str(out_dir.resolve())
    cont_out   = "/data/fastq"
    threads    = cfg["processing"]["threads"]
    sra_cache  = cfg["sra"]["cache_dir"]   # inside container

    if dry_run:
        log(f"  [DRY RUN] {srr} → {out_dir}")
        return {"srr": srr, "status": "dry_run", "fastq_dir": str(out_dir)}

    log(f"  {srr}: prefetch...")
    # ── prefetch ──
    prefetch_cmd = (
        f"prefetch {srr} "
        f"--output-directory {sra_cache} "
        f"--max-size 50g "
        f"-p"   # show progress
    )
    rc, out, err = docker_run(
        cfg, prefetch_cmd,
        extra_mounts=[f"-v {host_out}:{cont_out}"],
    )
    if rc != 0:
        log(f"  {srr}: prefetch failed\n{err[:300]}", "error")
        return {"srr": srr, "status": "prefetch_failed", "error": err[:500]}

    log(f"  {srr}: fasterq-dump...")
    # ── fasterq-dump ──
    fasterq_cmd = (
        f"fasterq-dump {sra_cache}/{srr}/{srr}.sra "
        f"--outdir {cont_out} "
        f"--threads {threads} "
        f"--split-3 "    # paired-end aware: _1.fastq _2.fastq; single-end: .fastq
        f"--skip-technical"
    )
    rc, out, err = docker_run(
        cfg, fasterq_cmd,
        extra_mounts=[f"-v {host_out}:{cont_out}"],
    )
    if rc != 0:
        log(f"  {srr}: fasterq-dump failed\n{err[:300]}", "error")
        return {"srr": srr, "status": "fasterq_failed", "error": err[:500]}

    log(f"  {srr}: compressing FASTQ...")
    # ── gzip ──
    gzip_cmd = f"pigz -p {threads} {cont_out}/*.fastq"
    rc, out, err = docker_run(
        cfg, gzip_cmd,
        extra_mounts=[f"-v {host_out}:{cont_out}"],
    )
    if rc != 0:
        # pigz may not be available, fall back to gzip
        gzip_cmd = f"gzip {cont_out}/*.fastq"
        rc, out, err = docker_run(
            cfg, gzip_cmd,
            extra_mounts=[f"-v {host_out}:{cont_out}"],
        )

    # ── delete SRA cache ──
    del_cache_cmd = f"rm -rf {sra_cache}/{srr}"
    docker_run(cfg, del_cache_cmd)

    # Collect output FASTQ files
    fastq_files = sorted(out_dir.glob("*.fastq.gz"))
    total_bytes = sum(f.stat().st_size for f in fastq_files)

    if not fastq_files:
        log(f"  {srr}: no FASTQ files produced", "error")
        return {"srr": srr, "status": "no_output", "error": "no fastq.gz files found"}

    log(f"  {srr}: done — {len(fastq_files)} FASTQ files, "
        f"{total_bytes/1e9:.1f} GB", "success")

    return {
        "srr":         srr,
        "gse":         gse,
        "gsm":         gsm,
        "omicsType":   omics_type,
        "status":      "downloaded",
        "fastq_dir":   str(out_dir),
        "fastq_files": [str(f) for f in fastq_files],
        "size_bytes":  total_bytes,
        "downloaded_at": datetime.now().isoformat(),
    }


# ════════════════════════════════════════════════════════════
#  MANIFEST UPDATER
# ════════════════════════════════════════════════════════════

def load_manifest(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def save_manifest(manifest: dict, path: Path):
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2)


# ════════════════════════════════════════════════════════════
#  MAIN RUNNER
# ════════════════════════════════════════════════════════════

def run(manifest_path: Path, cfg: dict,
        filter_srrs: list[str] = None,
        dry_run: bool = False):

    manifest  = load_manifest(manifest_path)
    runs      = manifest.get("runs", [])
    fastq_tmp = Path(cfg["paths"]["fastq_tmp"]).resolve()
    fastq_tmp.mkdir(parents=True, exist_ok=True)

    # Filter to specific SRRs if requested
    if filter_srrs:
        runs = [r for r in runs if r["srr"] in filter_srrs]
        log(f"Filtered to {len(runs)} specific SRR(s)")

    # Skip already-downloaded runs
    pending = [r for r in runs if r.get("status") not in ("downloaded","skipped")]
    skipped = len(runs) - len(pending)
    if skipped:
        log(f"Skipping {skipped} already-downloaded SRRs")

    total_gb = sum(r.get("size_mb",0) for r in pending) / 1024
    log(f"{len(pending)} SRR runs to download, est. {total_gb:.1f} GB")

    if HAS_RICH:
        t = Table(title="Download Plan", box=None, header_style="bold", padding=(0,2))
        t.add_column("SRR",      style="cyan")
        t.add_column("GSE")
        t.add_column("GSM")
        t.add_column("Omics")
        t.add_column("Est. size", justify="right")
        for r in pending[:20]:  # show first 20
            t.add_row(r["srr"], r.get("gse",""), r.get("gsm",""),
                      r.get("omicsType",""), fmt_size(r.get("size_mb",0)))
        if len(pending) > 20:
            t.add_row(f"... and {len(pending)-20} more", "", "", "", "")
        console.print(t)

    if not dry_run:
        confirm = input(f"\nProceed with downloading {len(pending)} SRR runs? [y/N] ").strip().lower()
        if confirm != "y":
            log("Aborted.", "warn"); return

    # Build a lookup for fast manifest updates
    run_by_srr = {r["srr"]: r for r in manifest["runs"]}

    # Process ONE SRR at a time to keep peak disk use low
    success = error = 0
    for i, run_entry in enumerate(pending):
        srr        = run_entry["srr"]
        gse        = run_entry.get("gse","unknown")
        gsm        = run_entry.get("gsm","unknown")
        omics_type = run_entry.get("omicsType","")

        log(f"[{i+1}/{len(pending)}] {srr}  ({gse} / {gsm})")

        result = download_srr(
            srr        = srr,
            gse        = gse,
            gsm        = gsm,
            omics_type = omics_type,
            cfg        = cfg,
            fastq_tmp  = fastq_tmp,
            dry_run    = dry_run,
        )

        # Update manifest entry in-place and save immediately
        # (so a crash doesn't lose progress)
        if srr in run_by_srr:
            run_by_srr[srr].update(result)
        save_manifest(manifest, manifest_path)

        if result["status"] == "downloaded":
            success += 1
        else:
            error += 1

    log(f"Downloads complete: {success} ok, {error} errors", "success" if not error else "warn")
    log(f"FASTQ files are in: {fastq_tmp}")
    log(f"Next step:  python processing/processor.py --manifest {manifest_path}")


# ════════════════════════════════════════════════════════════
#  CLI
# ════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Download SRA runs from sra_manifest.json via Docker"
    )
    parser.add_argument("--manifest", "-m", required=True,
                        help="Path to sra_manifest.json from geo_downloader.py")
    parser.add_argument("--srr", nargs="+",
                        help="Download only specific SRR accessions")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--config",  default="pipeline_config.yaml")
    args = parser.parse_args()

    global CONFIG_FILE
    CONFIG_FILE = Path(args.config)
    cfg = load_config()

    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        print(f"Manifest not found: {manifest_path}"); sys.exit(1)

    if not check_docker(cfg):
        sys.exit(1)

    try:
        run(manifest_path, cfg,
            filter_srrs = args.srr,
            dry_run     = args.dry_run)
    except KeyboardInterrupt:
        print("\nInterrupted. Re-run to resume — completed SRRs are marked in the manifest.")
    except Exception as e:
        print(f"Error: {e}"); raise


if __name__ == "__main__":
    main()
