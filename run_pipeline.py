#!/usr/bin/env python3
"""
run_pipeline.py — End-to-end omics pipeline orchestrator
──────────────────────────────────────────────────────────
Chains all pipeline steps in order:
  1. omics_agent.py       — disease curation → curated CSV
  2. geo_downloader.py    — GEO resolution → sra_manifest.json
  3. sra_downloader.py    — FASTQ download from SRA
  4. processor.py         — alignment + counting + cleanup

Each step saves its state so the pipeline can be resumed
at any point after a failure or interruption.

Usage:
    # Full pipeline
    python run_pipeline.py --disease "psoriasis" --model llama3.2

    # Resume from a specific step (skip earlier steps)
    python run_pipeline.py --disease "psoriasis" --from-step geo
    python run_pipeline.py --disease "psoriasis" --from-step sra
    python run_pipeline.py --disease "psoriasis" --from-step process

    # Dry run — show what would happen without downloading/processing
    python run_pipeline.py --disease "psoriasis" --model llama3.2 --dry-run

    # Skip reference setup check (if already done)
    python run_pipeline.py --disease "psoriasis" --skip-ref-check

Steps:  curate → geo → sra → process

Requirements:
    pip install requests rich pyyaml
    Docker running with omics-pipeline:1.0 image built
    reference_setup.py already run
"""

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import yaml

try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.markup import escape
    HAS_RICH = True
    console = Console()
except ImportError:
    HAS_RICH = False
    console = None

CONFIG_FILE  = Path("pipeline_config.yaml")
STEP_ORDER   = ["curate", "geo", "sra", "process"]


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


def section(title: str):
    bar = "═" * 60
    if HAS_RICH:
        console.print(f"\n[bold green]{bar}[/bold green]")
        console.print(f"[bold]  {title}[/bold]")
        console.print(f"[bold green]{bar}[/bold green]\n")
    else:
        print(f"\n{bar}\n  {title}\n{bar}\n")


# ════════════════════════════════════════════════════════════
#  CONFIG
# ════════════════════════════════════════════════════════════

def load_config() -> dict:
    if not CONFIG_FILE.exists():
        log(f"pipeline_config.yaml not found.", "error"); sys.exit(1)
    with open(CONFIG_FILE) as f:
        return yaml.safe_load(f)


# ════════════════════════════════════════════════════════════
#  STATE FILE  (tracks what has completed, enables resume)
# ════════════════════════════════════════════════════════════

def state_path(disease: str, cfg: dict) -> Path:
    safe = disease.lower().replace(" ","_").replace("/","_")
    return Path(cfg["paths"]["curation_out"]) / f".pipeline_state_{safe}.json"


def load_state(disease: str, cfg: dict) -> dict:
    sp = state_path(disease, cfg)
    if sp.exists():
        with open(sp) as f:
            return json.load(f)
    return {"disease": disease, "steps": {}}


def save_state(state: dict, disease: str, cfg: dict):
    sp = state_path(disease, cfg)
    sp.parent.mkdir(parents=True, exist_ok=True)
    with open(sp, "w") as f:
        json.dump(state, f, indent=2)


def mark_done(state: dict, step: str, meta: dict, disease: str, cfg: dict):
    state["steps"][step] = {
        "status":      "done",
        "completed_at": datetime.now().isoformat(),
        **meta,
    }
    save_state(state, disease, cfg)


# ════════════════════════════════════════════════════════════
#  SUBPROCESS RUNNER
# ════════════════════════════════════════════════════════════

def run_step(cmd: list[str], step_name: str) -> bool:
    """Run a pipeline step as a subprocess, streaming output to terminal."""
    log(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        log(f"Step '{step_name}' failed (exit {result.returncode})", "error")
        return False
    return True


# ════════════════════════════════════════════════════════════
#  REFERENCE CHECK
# ════════════════════════════════════════════════════════════

def check_references(cfg: dict) -> bool:
    ref_dir  = Path(cfg["paths"]["references"])
    ref_json = ref_dir / "references.json"
    if not ref_json.exists():
        log("references.json not found — run reference_setup.py first", "error")
        log("  python processing/reference_setup.py", "warn")
        return False
    log("References found", "success")
    return True


def check_docker(cfg: dict) -> bool:
    image  = cfg["docker"]["image"]
    result = subprocess.run(f"docker image inspect {image}",
                            shell=True, capture_output=True)
    if result.returncode != 0:
        log(f"Docker image '{image}' not found", "error")
        log("  docker build -t omics-pipeline:1.0 processing/", "warn")
        return False
    log(f"Docker image '{image}' found", "success")
    return True


# ════════════════════════════════════════════════════════════
#  PIPELINE STEPS
# ════════════════════════════════════════════════════════════

def step_curate(disease: str, model: str, cfg: dict,
                state: dict, dry_run: bool) -> tuple[bool, dict]:
    """Step 1: Run omics_agent.py to produce curated CSV."""
    section("Step 1/4 — Disease curation")
    safe    = disease.lower().replace(" ","_")
    out_dir = cfg["paths"]["curation_out"]
    csv_out = Path(out_dir) / f"omics_datasets_{safe}_curated.csv"

    # If curated CSV already exists and step is marked done, skip
    if state["steps"].get("curate",{}).get("status") == "done" and csv_out.exists():
        log(f"Curation already done: {csv_out}", "success")
        return True, {"curated_csv": str(csv_out)}

    cmd = [
        sys.executable, "data_curation/omics_agent.py",
        "--disease", disease,
        "--model",   model,
        "--out",     out_dir,
    ]
    if dry_run:
        log(f"[DRY RUN] Would run: {' '.join(cmd)}")
        return True, {"curated_csv": str(csv_out)}

    ok = run_step(cmd, "curate")
    if not ok or not csv_out.exists():
        log("Curation failed or produced no CSV", "error")
        return False, {}

    log(f"Curated CSV: {csv_out}", "success")
    return True, {"curated_csv": str(csv_out)}


def step_geo(disease: str, cfg: dict, state: dict,
             curated_csv: str, dry_run: bool) -> tuple[bool, dict]:
    """Step 2: Run geo_downloader.py to resolve GSM→SRR."""
    section("Step 2/4 — GEO resolution")
    safe     = disease.lower().replace(" ","_")
    out_dir  = cfg["paths"]["geo_out"]
    manifest = Path(out_dir) / "sra_manifest.json"

    if state["steps"].get("geo",{}).get("status") == "done" and manifest.exists():
        log(f"GEO resolution already done: {manifest}", "success")
        return True, {"manifest": str(manifest)}

    cmd = [
        sys.executable, "data_curation/geo_downloader.py",
        "--csv",     curated_csv,
        "--out",     out_dir,
        "--disease", disease,
    ]
    if dry_run:
        cmd.append("--dry-run")

    ok = run_step(cmd, "geo")
    if not ok:
        return False, {}

    return True, {"manifest": str(manifest)}


def step_sra(disease: str, cfg: dict, state: dict,
             manifest: str, dry_run: bool) -> tuple[bool, dict]:
    """Step 3: Run sra_downloader.py to download FASTQ files."""
    section("Step 3/4 — SRA download")

    if state["steps"].get("sra",{}).get("status") == "done":
        log("SRA download already marked done", "success")
        return True, {}

    # Check if there are any SRA-linked datasets in the manifest
    with open(manifest) as f:
        mf = json.load(f)
    if not mf.get("runs"):
        log("No SRA runs in manifest — all datasets are GEO-direct", "success")
        return True, {"note": "no_sra_runs"}

    total_gb = mf.get("total_size_gb", 0)
    log(f"Manifest has {len(mf['runs'])} SRR runs, estimated {total_gb:.1f} GB")

    cmd = [
        sys.executable, "processing/sra_downloader.py",
        "--manifest", manifest,
    ]
    if dry_run:
        cmd.append("--dry-run")

    ok = run_step(cmd, "sra")
    return ok, {}


def step_process(disease: str, cfg: dict, state: dict,
                 manifest: str, dry_run: bool) -> tuple[bool, dict]:
    """Step 4: Run processor.py to align, count, and clean up."""
    section("Step 4/4 — Alignment + counting")

    if state["steps"].get("process",{}).get("status") == "done":
        log("Processing already done", "success")
        return True, {}

    cmd = [
        sys.executable, "processing/processor.py",
        "--manifest", manifest,
    ]
    if dry_run:
        cmd.append("--dry-run")

    ok = run_step(cmd, "process")
    return ok, {}


# ════════════════════════════════════════════════════════════
#  MAIN ORCHESTRATOR
# ════════════════════════════════════════════════════════════

def run_pipeline(disease: str, model: str, cfg: dict,
                 from_step: str = "curate",
                 dry_run: bool = False,
                 skip_ref_check: bool = False):

    if HAS_RICH:
        console.print(Panel.fit(
            f"[bold]Omics Pipeline[/bold]\n"
            f"[dim]Disease:[/dim] {escape(disease)}\n"
            f"[dim]Model:[/dim]   {model}\n"
            f"[dim]Start:[/dim]   step {from_step}"
            + (" [yellow](DRY RUN)[/yellow]" if dry_run else ""),
            border_style="green"
        ))

    # Pre-flight checks
    if not skip_ref_check:
        if not check_references(cfg): sys.exit(1)
    if not check_docker(cfg): sys.exit(1)

    state = load_state(disease, cfg)
    start_idx = STEP_ORDER.index(from_step) if from_step in STEP_ORDER else 0

    # Carry forward paths from state when resuming
    curated_csv = state["steps"].get("curate", {}).get("curated_csv", "")
    manifest    = state["steps"].get("geo",    {}).get("manifest", "")

    t_start = datetime.now()

    # ── Step 1: Curate ──
    if start_idx <= STEP_ORDER.index("curate"):
        ok, meta = step_curate(disease, model, cfg, state, dry_run)
        if not ok: sys.exit(1)
        mark_done(state, "curate", meta, disease, cfg)
        curated_csv = meta.get("curated_csv", curated_csv)

    # ── Step 2: GEO resolution ──
    if start_idx <= STEP_ORDER.index("geo"):
        if not curated_csv or not Path(curated_csv).exists():
            log("No curated CSV found. Run from step 'curate' first.", "error")
            sys.exit(1)
        ok, meta = step_geo(disease, cfg, state, curated_csv, dry_run)
        if not ok: sys.exit(1)
        mark_done(state, "geo", meta, disease, cfg)
        manifest = meta.get("manifest", manifest)

    # ── Step 3: SRA download ──
    if start_idx <= STEP_ORDER.index("sra"):
        if not manifest or not Path(manifest).exists():
            log("No sra_manifest.json found. Run from step 'geo' first.", "error")
            sys.exit(1)
        ok, meta = step_sra(disease, cfg, state, manifest, dry_run)
        if not ok: sys.exit(1)
        mark_done(state, "sra", meta, disease, cfg)

    # ── Step 4: Process ──
    if start_idx <= STEP_ORDER.index("process"):
        if not manifest or not Path(manifest).exists():
            log("No sra_manifest.json found. Run from step 'sra' first.", "error")
            sys.exit(1)
        ok, meta = step_process(disease, cfg, state, manifest, dry_run)
        if not ok: sys.exit(1)
        mark_done(state, "process", meta, disease, cfg)

    elapsed = (datetime.now() - t_start).total_seconds()
    mins    = int(elapsed // 60)
    secs    = int(elapsed % 60)

    processed_dir = Path(cfg["paths"]["processed_out"])
    if HAS_RICH:
        console.print(Panel(
            f"[green]✓  Pipeline complete[/green]  ({mins}m {secs}s)\n"
            f"[dim]Disease:   {disease}[/dim]\n"
            f"[dim]Output:    {processed_dir}[/dim]\n"
            f"[dim]Curated:   {curated_csv}[/dim]",
            title="Done", border_style="green"
        ))
    else:
        print(f"\nPipeline complete ({mins}m {secs}s)")
        print(f"Output: {processed_dir}")


# ════════════════════════════════════════════════════════════
#  CLI
# ════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="End-to-end omics pipeline: curation → download → process",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Steps:
  curate   — run omics_agent.py (LLM disease curation)
  geo      — run geo_downloader.py (GEO → SRR manifest)
  sra      — run sra_downloader.py (FASTQ download)
  process  — run processor.py (align + count + cleanup)

Examples:
  # Full pipeline
  python run_pipeline.py --disease "psoriasis" --model llama3.2

  # Resume from GEO step (curation already done)
  python run_pipeline.py --disease "psoriasis" --from-step geo

  # Dry run
  python run_pipeline.py --disease "psoriasis" --model llama3.2 --dry-run
        """
    )
    parser.add_argument("--disease",       "-d", required=True)
    parser.add_argument("--model",         "-m", default="mistral",
                        help="Ollama model for curation (default: mistral)")
    parser.add_argument("--from-step",     "-s",
                        choices=STEP_ORDER, default="curate",
                        help="Start from this step (default: curate)")
    parser.add_argument("--dry-run",       action="store_true")
    parser.add_argument("--skip-ref-check",action="store_true",
                        help="Skip reference genome check (if already confirmed)")
    parser.add_argument("--config",        default="pipeline_config.yaml")
    args = parser.parse_args()

    global CONFIG_FILE
    CONFIG_FILE = Path(args.config)
    cfg = load_config()

    try:
        run_pipeline(
            disease        = args.disease,
            model          = args.model,
            cfg            = cfg,
            from_step      = args.from_step,
            dry_run        = args.dry_run,
            skip_ref_check = args.skip_ref_check,
        )
    except KeyboardInterrupt:
        print("\nInterrupted. Re-run with --from-step to resume.")
    except Exception as e:
        print(f"Error: {e}"); raise


if __name__ == "__main__":
    main()
