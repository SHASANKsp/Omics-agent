#!/usr/bin/env python3
"""
reference_setup.py — One-time reference genome setup
──────────────────────────────────────────────────────
Downloads GRCh38 genome + GTF from Ensembl and builds:
  • STAR index       (bulk RNA-seq + STARsolo scRNA-seq)
  • HISAT2 index     (downloaded pre-built from AWS — no build needed)
  • Bismark index    (methylation)
  • STARsolo whitelist (10x Chromium v3 barcode list)

Everything is downloaded automatically — nothing needs to be done manually.
Decompression and index building run inside Docker (works on Windows).
16 GB RAM is fine — STAR adapts automatically, just takes longer (~2 hours).

Usage:
    python processing/reference_setup.py              # full setup
    python processing/reference_setup.py --check      # verify existing files
    python processing/reference_setup.py --skip-star
    python processing/reference_setup.py --skip-bismark
    python processing/reference_setup.py --genome-only

Requirements:
    Docker running with omics-pipeline:1.0 image built.
    pip install requests rich pyyaml
"""

import argparse
import json
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

import requests
import yaml

try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    from rich.progress import (
        Progress, BarColumn, DownloadColumn,
        TransferSpeedColumn, TimeRemainingColumn,
        TimeElapsedColumn, TextColumn, SpinnerColumn,
        TaskProgressColumn,
    )
    from rich.markup import escape
    from rich.live import Live
    from rich.text import Text
    HAS_RICH = True
    console = Console()
except ImportError:
    HAS_RICH = False
    console = None

CONFIG_FILE = Path("pipeline_config.yaml")
CHUNK_SIZE  = 2 * 1024 * 1024   # 2 MB chunks (finer progress granularity)


# ════════════════════════════════════════════════════════════
#  LOGGING
# ════════════════════════════════════════════════════════════

def log(msg: str, status: str = "info"):
    ts  = datetime.now().strftime("%H:%M:%S")
    sym = {"info":"·","success":"✓","warn":"!","error":"✗"}.get(status, "·")
    if HAS_RICH:
        col = {"info":"dim","success":"green","warn":"yellow","error":"red"}.get(status,"white")
        console.print(f"[dim]{ts}[/dim]  [{col}]{sym}[/{col}]  {escape(str(msg))}")
    else:
        print(f"{ts}  {sym}  {msg}")


# ════════════════════════════════════════════════════════════
#  CONFIG
# ════════════════════════════════════════════════════════════

def load_config() -> dict:
    if not CONFIG_FILE.exists():
        log(f"Config not found: {CONFIG_FILE}", "error")
        sys.exit(1)
    with open(CONFIG_FILE) as f:
        return yaml.safe_load(f)


# ════════════════════════════════════════════════════════════
#  DOWNLOAD WITH RICH PROGRESS BAR
# ════════════════════════════════════════════════════════════

def fmt_size(n: int) -> str:
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {u}"
        n /= 1024
    return f"{n:.1f} TB"


def download_file(url: str, dest: Path, label: str = "") -> bool:
    """
    Download a file with a live rich progress bar showing:
      filename | bar | downloaded size | speed | time remaining
    Supports resume — if dest exists and is partial, continues from where
    it left off (HTTP Range request).
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    existing = dest.stat().st_size if dest.exists() else 0
    headers  = {"User-Agent": "OmicsSetup/1.0"}
    if existing:
        headers["Range"] = f"bytes={existing}-"

    desc = label or dest.name

    try:
        r = requests.get(url, headers=headers, stream=True, timeout=60)

        if r.status_code == 416:
            # Range not satisfiable = file already complete
            log(f"{desc}: already complete ({fmt_size(existing)})", "success")
            return True
        if r.status_code == 404:
            log(f"{desc}: not found at URL", "error")
            return False
        r.raise_for_status()

        content_len = int(r.headers.get("content-length", 0))
        total       = content_len + existing   # total including already-downloaded
        mode        = "ab" if existing else "wb"

        if HAS_RICH:
            with Progress(
                TextColumn("[bold cyan]{task.description}"),
                BarColumn(bar_width=35),
                TaskProgressColumn(),
                DownloadColumn(),
                TransferSpeedColumn(),
                TimeRemainingColumn(),
                TimeElapsedColumn(),
                console=console,
                transient=False,   # keep bar visible after completion
            ) as progress:
                resume_note = f" (resuming from {fmt_size(existing)})" if existing else ""
                task = progress.add_task(
                    f"{desc}{resume_note}",
                    total=total,
                    completed=existing,
                )
                with open(dest, mode) as f:
                    for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                        if chunk:
                            f.write(chunk)
                            progress.advance(task, len(chunk))
        else:
            # Fallback: plain text progress every 5 seconds
            done     = existing
            last_log = time.time()
            with open(dest, mode) as f:
                for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                    if chunk:
                        f.write(chunk)
                        done += len(chunk)
                        if time.time() - last_log > 5:
                            pct = done / total * 100 if total else 0
                            print(f"  {desc}: {fmt_size(done)} / {fmt_size(total)} ({pct:.0f}%)")
                            last_log = time.time()

        final = dest.stat().st_size
        log(f"{desc}: {fmt_size(final)} — complete", "success")
        return True

    except requests.exceptions.ConnectionError as e:
        log(f"{desc}: connection error — {e}", "error")
        return False
    except Exception as e:
        log(f"{desc}: {e}", "error")
        return False


# ════════════════════════════════════════════════════════════
#  DOCKER RUNNER WITH LIVE SPINNER
# ════════════════════════════════════════════════════════════

def check_docker_image(cfg: dict) -> bool:
    image = cfg["docker"]["image"]
    r = subprocess.run(f"docker image inspect {image}",
                       shell=True, capture_output=True)
    if r.returncode != 0:
        log(f"Docker image '{image}' not found.", "error")
        log("Build it first:  docker build -t omics-pipeline:1.0 processing/", "warn")
        return False
    log(f"Docker image '{image}' found ✓", "success")
    return True


def win_docker_mount(path: Path) -> str:
    """Convert Windows path to Docker-compatible mount path."""
    p = str(path.resolve()).replace("\\", "/")
    if len(p) > 1 and p[1] == ":":
        p = "/" + p[0].lower() + p[2:]
    return p


def docker_run(cfg: dict, command: str, spinner_label: str = "") -> bool:
    """
    Run a command inside the omics Docker container.
    Shows a live spinner with elapsed time while the command runs.
    Streams Docker stdout/stderr to the terminal so you can see progress.
    """
    ref_dir   = Path(cfg["paths"]["references"]).resolve()
    image     = cfg["docker"]["image"]
    ref_mount = win_docker_mount(ref_dir)
    ref_dir.mkdir(parents=True, exist_ok=True)

    cmd = (
        f'docker run --rm '
        f'-v "{ref_mount}:/references" '
        f'{image} '
        f'bash -c "{command}"'
    )

    label = spinner_label or command[:60]

    if not HAS_RICH:
        print(f"  Running: {label}...")
        result = subprocess.run(cmd, shell=True, text=True)
        if result.returncode != 0:
            print(f"  ✗ Failed (exit {result.returncode})")
            return False
        return True

    # Rich: show spinner + elapsed time while Docker runs in a thread
    start      = time.time()
    result_box = [None]
    stop_event = threading.Event()

    def run_docker():
        # Stream output so STAR/Bismark progress lines appear in terminal
        result_box[0] = subprocess.run(cmd, shell=True, text=True)
        stop_event.set()

    thread = threading.Thread(target=run_docker, daemon=True)
    thread.start()

    with Progress(
        SpinnerColumn(),
        TextColumn("[bold cyan]{task.description}"),
        TimeElapsedColumn(),
        console=console,
        transient=True,   # spinner disappears when done, replaced by log line
    ) as progress:
        task = progress.add_task(label, total=None)
        while not stop_event.wait(timeout=1.0):
            elapsed = time.time() - start
            progress.update(task, description=f"{label}  [dim]({elapsed:.0f}s)[/dim]")

    thread.join()
    rc = result_box[0].returncode if result_box[0] else 1

    elapsed = time.time() - start
    if rc != 0:
        log(f"{label}: failed (exit {rc}) after {elapsed:.0f}s", "error")
        return False

    log(f"{label}: done in {elapsed:.0f}s", "success")
    return True


# ════════════════════════════════════════════════════════════
#  STEP 1 — GENOME + GTF
# ════════════════════════════════════════════════════════════

def download_genome(cfg: dict) -> bool:
    release   = cfg["reference"]["ensembl_release"]
    ref_dir   = Path(cfg["paths"]["references"])
    genome_fa = ref_dir / cfg["reference"]["genome_fasta"]
    gtf_file  = ref_dir / cfg["reference"]["gtf"]
    genome_gz = ref_dir / "GRCh38" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    gtf_gz    = ref_dir / "GRCh38" / f"Homo_sapiens.GRCh38.{release}.gtf.gz"

    genome_url = (
        f"https://ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/dna/"
        f"Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    )
    gtf_url = (
        f"https://ftp.ensembl.org/pub/release-{release}/gtf/homo_sapiens/"
        f"Homo_sapiens.GRCh38.{release}.gtf.gz"
    )

    # ── Genome FASTA ──
    if genome_fa.exists():
        log(f"Genome FASTA already present ({fmt_size(genome_fa.stat().st_size)})", "success")
    else:
        log("Downloading GRCh38 genome FASTA from Ensembl (~900 MB compressed)...")
        if not download_file(genome_url, genome_gz, "Genome FASTA (.fa.gz)"):
            return False
        log("Decompressing genome FASTA inside Docker (~5 min)...")
        if not docker_run(cfg,
                          f"gunzip -k /references/GRCh38/{genome_gz.name}",
                          "Decompressing genome FASTA"):
            return False
        log("Genome FASTA ready", "success")

    # ── GTF ──
    if gtf_file.exists():
        log(f"GTF already present ({fmt_size(gtf_file.stat().st_size)})", "success")
    else:
        log("Downloading GRCh38 GTF from Ensembl (~60 MB compressed)...")
        if not download_file(gtf_url, gtf_gz,
                             f"GTF annotation (.gtf.gz)"):
            return False
        log("Decompressing GTF inside Docker...")
        if not docker_run(cfg,
                          f"gunzip -k /references/GRCh38/{gtf_gz.name}",
                          "Decompressing GTF"):
            return False
        log("GTF ready", "success")

    return True


# ════════════════════════════════════════════════════════════
#  STEP 2 — STARsolo whitelist
# ════════════════════════════════════════════════════════════

def download_whitelist(cfg: dict) -> bool:
    """
    Download the 10x Chromium v3 cell barcode whitelist for STARsolo.
    Uses wget inside Docker to handle GitHub Pages redirects reliably.
    """
    ref_dir   = Path(cfg["paths"]["references"])
    whitelist = cfg["processing"]["starsolo"]["whitelist"]
    dest_txt  = ref_dir / whitelist

    if dest_txt.exists():
        log(f"Whitelist already present ({fmt_size(dest_txt.stat().st_size)})", "success")
        return True

    log("Downloading STARsolo 10x Chromium v3 barcode whitelist...")
    # Use wget inside Docker — handles GitHub Pages and S3 redirects transparently.
    # The teichlab scg_lib_structs repo is the canonical public mirror
    # (10x Genomics removed this file from their CellRanger GitHub repo in 2024).
    url = "https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz"
    ok = docker_run(cfg,
        f"wget -q --show-progress '{url}' -O /references/{whitelist}.gz "
        f"&& gunzip -f /references/{whitelist}.gz",
        "Downloading + decompressing whitelist"
    )
    if not ok:
        return False
    log(f"Whitelist ready: {dest_txt}", "success")
    return True


# ════════════════════════════════════════════════════════════
#  STEP 3 — STAR index
# ════════════════════════════════════════════════════════════

def build_star_index(cfg: dict) -> bool:
    """
    Build STAR index inside Docker.
    Uses all available RAM automatically (reads /proc/meminfo at runtime).
    With 16 GB RAM: takes ~1.5-2 hours. With 30 GB: ~1 hour.
    --genomeSAindexNbases 13 halves SA-index RAM with no accuracy loss.
    """
    ref_cfg = cfg["reference"]
    threads = cfg["processing"]["threads"]
    idx_dir = Path(cfg["paths"]["references"]) / ref_cfg["star_index"]

    if (idx_dir / "genomeParameters.txt").exists():
        log("STAR index already built ✓", "success")
        return True

    log("Building STAR index inside Docker...")
    log("  This takes 1.5–2 hours on 16 GB RAM — safe to leave running", "warn")

    # RAM budget for 16 GB machine: 13 GB to STAR, leaving ~3 GB headroom.
    # --genomeSAindexNbases 13: halves SA index size vs default 14
    # --genomeSAsparseD 3:      sparse suffix array — cuts peak build RAM from
    #                           ~15 GB to ~8 GB. Alignment uses slightly more
    #                           RAM per run but well within 16 GB.
    # --genomeChrBinNbits 18:   reduces chromosome bin table memory
    ram_bytes = 13 * 1024 * 1024 * 1024   # 13 GB in bytes

    cmd = (
        "mkdir -p /references/star_index && "
        "STAR --runMode genomeGenerate "
        "--genomeDir /references/star_index "
        f"--genomeFastaFiles /references/{ref_cfg['genome_fasta']} "
        f"--sjdbGTFfile /references/{ref_cfg['gtf']} "
        "--sjdbOverhang 100 "
        f"--runThreadN {threads} "
        f"--limitGenomeGenerateRAM {ram_bytes} "
        "--genomeSAindexNbases 13 "
        "--genomeSAsparseD 3 "
        "--genomeChrBinNbits 18"
    )
    return docker_run(cfg, cmd, "Building STAR index (1.5-2h on 16 GB RAM)")


# ════════════════════════════════════════════════════════════
#  STEP 4 — HISAT2 index (pre-built download from AWS)
# ════════════════════════════════════════════════════════════

def download_hisat2_index(cfg: dict) -> bool:
    """
    Download the pre-built HISAT2 GRCh38+transcript index from AWS.
    ~8 GB download, no build step, no RAM requirement.
    """
    ref_dir = Path(cfg["paths"]["references"])
    idx_dir = ref_dir / "hisat2_index"
    idx_dir.mkdir(parents=True, exist_ok=True)

    existing_ht2 = list(idx_dir.glob("*.ht2")) + list(idx_dir.glob("*.ht2l"))
    if existing_ht2:
        log(f"HISAT2 index already present ({len(existing_ht2)} files) ✓", "success")
        return True

    url      = "https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz"
    dest_tar = ref_dir / "grch38_tran.tar.gz"

    log("Downloading pre-built HISAT2 GRCh38 index from AWS (~8 GB)...")
    if not download_file(url, dest_tar, "HISAT2 GRCh38 index (.tar.gz)"):
        return False

    log("Extracting HISAT2 index inside Docker...")
    ok = docker_run(cfg,
        "tar -xzf /references/grch38_tran.tar.gz "
        "-C /references/hisat2_index --strip-components=1",
        "Extracting HISAT2 index"
    )
    if not ok:
        return False

    dest_tar.unlink(missing_ok=True)
    cfg["reference"]["hisat2_index"] = "hisat2_index/genome_tran"

    ht2 = list(idx_dir.glob("*.ht2")) + list(idx_dir.glob("*.ht2l"))
    log(f"HISAT2 index ready: {len(ht2)} index files", "success")
    return True


# ════════════════════════════════════════════════════════════
#  STEP 5 — Bismark index
# ════════════════════════════════════════════════════════════

def build_bismark_index(cfg: dict) -> bool:
    """
    Build Bismark bisulfite genome index inside Docker.
    Needs ~12 GB RAM, takes ~2 hours. Safe on 16 GB.
    """
    ref_dir    = Path(cfg["paths"]["references"])
    bismark_ok = ref_dir / "GRCh38" / "Bisulfite_Genome"

    if bismark_ok.exists():
        log("Bismark index already built ✓", "success")
        return True

    log("Building Bismark index inside Docker...")
    log("  This takes ~2 hours and uses ~12 GB RAM — safe on 16 GB", "warn")

    return docker_run(cfg,
        "bismark_genome_preparation --parallel 4 /references/GRCh38",
        "Building Bismark index (~2h)"
    )


# ════════════════════════════════════════════════════════════
#  CHECK
# ════════════════════════════════════════════════════════════

def check_indices(cfg: dict) -> dict:
    ref_dir = Path(cfg["paths"]["references"])
    ref_cfg = cfg["reference"]

    checks = {
        "genome_fasta":  ref_dir / ref_cfg["genome_fasta"],
        "gtf":           ref_dir / ref_cfg["gtf"],
        "star_index":    ref_dir / ref_cfg["star_index"] / "genomeParameters.txt",
        "hisat2_index":  ref_dir / "hisat2_index" / "genome_tran.1.ht2",
        "bismark_index": ref_dir / "GRCh38" / "Bisulfite_Genome",
        "whitelist":     ref_dir / cfg["processing"]["starsolo"]["whitelist"],
    }

    status = {}
    if HAS_RICH:
        t = Table(title="Reference Index Status", box=None,
                  header_style="bold", padding=(0, 2))
        t.add_column("Component", style="cyan")
        t.add_column("Path")
        t.add_column("Size", justify="right")
        t.add_column("Status")

    for name, path in checks.items():
        exists = path.exists()
        status[name] = exists
        size_str = fmt_size(path.stat().st_size) if exists and path.is_file() else "—"
        if HAS_RICH:
            c = "green" if exists else "red"
            t.add_row(name, str(path), size_str,
                      f"[{c}]{'✓ present' if exists else '✗ missing'}[/{c}]")
        else:
            print(f"  {'✓' if exists else '✗'}  {name:20s}  {path}")

    if HAS_RICH:
        console.print(t)

    return status


# ════════════════════════════════════════════════════════════
#  references.json
# ════════════════════════════════════════════════════════════

def write_references_json(cfg: dict):
    ref_dir = Path(cfg["paths"]["references"]).resolve()
    ref_cfg = cfg["reference"]

    refs = {
        "genome_fasta":    str(ref_dir / ref_cfg["genome_fasta"]),
        "gtf":             str(ref_dir / ref_cfg["gtf"]),
        "star_index":      str(ref_dir / ref_cfg["star_index"]),
        "hisat2_index":    str(ref_dir / ref_cfg["hisat2_index"]),
        "bismark_genome":  str(ref_dir / "GRCh38"),
        "whitelist":       str(ref_dir / cfg["processing"]["starsolo"]["whitelist"]),
        "built_at":        datetime.now().isoformat(),
        "assembly":        ref_cfg["assembly"],
        "ensembl_release": ref_cfg["ensembl_release"],
    }

    dest = ref_dir / "references.json"
    with open(dest, "w") as f:
        json.dump(refs, f, indent=2)
    log(f"References manifest: {dest}", "success")


# ════════════════════════════════════════════════════════════
#  DISK + TIME ESTIMATE
# ════════════════════════════════════════════════════════════

def print_plan():
    if HAS_RICH:
        t = Table(title="What will be downloaded / built", box=None,
                  header_style="bold", padding=(0, 2))
        t.add_column("Step")
        t.add_column("Component")
        t.add_column("Size",    justify="right")
        t.add_column("Est. time (16 GB RAM)")
        t.add_column("Method")
        rows = [
            ("1", "Genome FASTA (GRCh38)",    "~3 GB",  "~10 min",       "Download → decompress in Docker"),
            ("1", "GTF annotation (Ensembl)",  "~1 GB",  "~3 min",        "Download → decompress in Docker"),
            ("2", "STARsolo whitelist",         "~1 MB",  "seconds",       "Download"),
            ("3", "STAR genome index",          "~30 GB", "~1.5-2 hours",  "Built in Docker (uses avail. RAM)"),
            ("4", "HISAT2 genome index",        "~8 GB",  "~15 min",       "Download from AWS (pre-built)"),
            ("5", "Bismark genome index",       "~6 GB",  "~2 hours",      "Built in Docker"),
            ("",  "TOTAL disk needed",          "~48 GB", "~4-5 hours",    ""),
        ]
        for r in rows:
            t.add_row(*r)
        console.print(t)
        console.print()
    else:
        print("\nSetup plan:")
        print("  Step 1  Genome FASTA + GTF   ~4 GB    ~13 min   (download)")
        print("  Step 2  STARsolo whitelist    ~1 MB    seconds   (download)")
        print("  Step 3  STAR index            ~30 GB   ~2 hours  (build in Docker)")
        print("  Step 4  HISAT2 index          ~8 GB    ~15 min   (download, pre-built)")
        print("  Step 5  Bismark index         ~6 GB    ~2 hours  (build in Docker)")
        print("  TOTAL                         ~48 GB   ~4-5h\n")


# ════════════════════════════════════════════════════════════
#  MAIN
# ════════════════════════════════════════════════════════════

def main():
    global CONFIG_FILE
    parser = argparse.ArgumentParser(
        description="One-time reference genome setup — downloads and builds everything automatically",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python processing/reference_setup.py              # full setup (recommended)
  python processing/reference_setup.py --check      # verify what exists
  python processing/reference_setup.py --genome-only  # just download genome + GTF
  python processing/reference_setup.py --skip-star    # skip STAR build
  python processing/reference_setup.py --skip-bismark # skip Bismark build
        """
    )
    parser.add_argument("--check",         action="store_true",
                        help="Check which components exist and exit")
    parser.add_argument("--skip-genome",   action="store_true")
    parser.add_argument("--skip-star",     action="store_true")
    parser.add_argument("--skip-hisat2",   action="store_true")
    parser.add_argument("--skip-bismark",  action="store_true")
    parser.add_argument("--genome-only",   action="store_true",
                        help="Download genome + GTF only (no index building)")
    parser.add_argument("--config", "-c",  default="pipeline_config.yaml")
    args = parser.parse_args()

    CONFIG_FILE = Path(args.config)
    cfg = load_config()

    if HAS_RICH:
        console.print(Panel.fit(
            "[bold]Omics Pipeline — Reference Setup[/bold]\n"
            f"[dim]Genome:[/dim] GRCh38 · Ensembl release {cfg['reference']['ensembl_release']}\n"
            f"[dim]Output:[/dim] {Path(cfg['paths']['references']).resolve()}\n"
            f"[dim]RAM:   [/dim] Everything runs inside Docker — 16 GB is fine",
            border_style="green"
        ))

    # ── Check mode ──
    if args.check:
        check_indices(cfg)
        sys.exit(0)

    print_plan()

    if not check_docker_image(cfg):
        sys.exit(1)

    ok = True

    # ── Step 1: Genome + GTF ──
    if not args.skip_genome:
        log("━━ Step 1/5 — Genome + GTF ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        ok = download_genome(cfg)
        if not ok:
            log("Genome download failed. Fix the error above and re-run.", "error")
            sys.exit(1)
    else:
        log("Step 1/5 — Genome skipped (--skip-genome)", "warn")

    if args.genome_only:
        log("--genome-only: stopping here.", "success")
        write_references_json(cfg)
        sys.exit(0)

    # ── Step 2: Whitelist ──
    log("━━ Step 2/5 — STARsolo whitelist ━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    download_whitelist(cfg)

    # ── Step 3: STAR index ──
    if not args.skip_star:
        log("━━ Step 3/5 — STAR index ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        if not build_star_index(cfg):
            log("STAR index build failed.", "error")
            ok = False
    else:
        log("Step 3/5 — STAR index skipped (--skip-star)", "warn")

    # ── Step 4: HISAT2 index ──
    if not args.skip_hisat2:
        log("━━ Step 4/5 — HISAT2 index ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        if not download_hisat2_index(cfg):
            log("HISAT2 index download failed.", "error")
            ok = False
    else:
        log("Step 4/5 — HISAT2 index skipped (--skip-hisat2)", "warn")

    # ── Step 5: Bismark index ──
    if not args.skip_bismark:
        log("━━ Step 5/5 — Bismark index ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        if not build_bismark_index(cfg):
            log("Bismark index build failed.", "error")
            ok = False
    else:
        log("Step 5/5 — Bismark index skipped (--skip-bismark)", "warn")

    # ── Write manifest and final check ──
    write_references_json(cfg)

    log("━━ Final verification ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    status = check_indices(cfg)
    all_ok = all(status.values())

    if HAS_RICH:
        if all_ok:
            console.print(Panel(
                "[green bold]✓  All reference files are ready[/green bold]\n"
                "[dim]You can now run: python data_curation/omics_agent.py[/dim]",
                border_style="green"
            ))
        else:
            missing = [k for k, v in status.items() if not v]
            console.print(Panel(
                f"[yellow]! Setup finished with {len(missing)} missing component(s)[/yellow]\n"
                f"[dim]Missing: {', '.join(missing)}[/dim]\n"
                "[dim]Re-run with the appropriate --skip flags to retry just what failed.[/dim]",
                border_style="yellow"
            ))
    else:
        print("\n✓ Setup complete." if all_ok else "\n! Some components missing — re-run to retry.")


if __name__ == "__main__":
    main()
