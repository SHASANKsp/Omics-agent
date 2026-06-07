#!/usr/bin/env python3
"""
processor.py — Align, count, and clean up per SRR run
──────────────────────────────────────────────────────
Reads sra_manifest.json (updated by sra_downloader.py with fastq_files paths).
For each downloaded SRR, dispatches to the correct tool:

  RNA-seq (bulk)   → fastp → STAR → featureCounts → delete BAM + FASTQ
  scRNA-seq        → fastp → STARsolo → delete BAM + FASTQ
  ATAC-seq         → fastp → HISAT2 → MACS3 → delete BAM + FASTQ
  ChIP-seq         → fastp → HISAT2 → MACS3 → delete BAM + FASTQ
  Methylation      → fastp → Bismark → delete BAM + FASTQ

After all SRRs in a GSE are processed, merges per-sample count files
into a single GSE-level count matrix.

Usage:
    python processing/processor.py --manifest geo_downloads/sra_manifest.json
    python processing/processor.py --manifest ... --gse GSE314390  # one dataset only
    python processing/processor.py --manifest ... --dry-run

Requirements:
    Docker running with the omics-pipeline image built.
    reference_setup.py must have been run first.
    pip install pyyaml rich
"""

import argparse
import json
import subprocess
import sys
from collections import defaultdict
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

CONFIG_FILE = Path("pipeline_config.yaml")

# Map omicsType strings → pipeline branch
OMICS_MAP = {
    "rna-seq":       "rnaseq",
    "rna_seq":       "rnaseq",
    "bulk rna-seq":  "rnaseq",
    "scrna-seq":     "scrnaseq",
    "scrna_seq":     "scrnaseq",
    "single-cell":   "scrnaseq",
    "atac-seq":      "atacseq",
    "atac_seq":      "atacseq",
    "chip-seq":      "chipseq",
    "chip_seq":      "chipseq",
    "methylation":   "methylation",
    "bisulfite":     "methylation",
}

def resolve_pipeline(omics_type: str) -> str:
    return OMICS_MAP.get(omics_type.lower().strip(), "rnaseq")


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


# ════════════════════════════════════════════════════════════
#  CONFIG + REFERENCES
# ════════════════════════════════════════════════════════════

def load_config() -> dict:
    if not CONFIG_FILE.exists():
        log(f"Config not found: {CONFIG_FILE}", "error"); sys.exit(1)
    with open(CONFIG_FILE) as f:
        return yaml.safe_load(f)


def load_references(cfg: dict) -> dict:
    ref_json = Path(cfg["paths"]["references"]) / "references.json"
    if not ref_json.exists():
        log("references.json not found. Run reference_setup.py first.", "error")
        sys.exit(1)
    with open(ref_json) as f:
        return json.load(f)


# ════════════════════════════════════════════════════════════
#  DOCKER RUNNER
# ════════════════════════════════════════════════════════════

def docker_run(cfg: dict, command: str,
               mounts: list[str] = None) -> tuple[int, str, str]:
    image = cfg["docker"]["image"]
    mount_str = " ".join(mounts or [])
    cmd = f"docker run --rm {mount_str} {image} bash -c \"{command}\""
    log(f"    CMD: {command[:100]}...")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        log(f"    STDERR: {r.stderr[-500:]}", "warn")
    return r.returncode, r.stdout, r.stderr


def make_mounts(cfg: dict, refs: dict,
                fastq_dir: Path, out_dir: Path) -> list[str]:
    """Build Docker -v mount arguments."""
    ref_dir = Path(cfg["paths"]["references"]).resolve()
    return [
        f"-v {ref_dir}:/references:ro",
        f"-v {fastq_dir.resolve()}:/fastq:ro",
        f"-v {out_dir.resolve()}:/output",
    ]


# ════════════════════════════════════════════════════════════
#  FASTP QC (all pipelines)
# ════════════════════════════════════════════════════════════

def run_fastp(srr: str, fastq_dir: Path, out_dir: Path,
              cfg: dict, refs: dict) -> tuple[bool, list[str]]:
    """
    QC + adapter trimming with fastp.
    Auto-detects paired-end (two _1/_2 files) vs single-end.
    Returns (success, list of trimmed fastq paths inside container).
    """
    threads = cfg["processing"]["threads"]
    fq_cfg  = cfg["processing"]["fastp"]

    fastq_files = sorted(fastq_dir.glob("*.fastq.gz"))
    if not fastq_files:
        log(f"    {srr}: no FASTQ files found in {fastq_dir}", "error")
        return False, []

    paired = any("_1.fastq.gz" in f.name for f in fastq_files)

    if paired:
        r1 = f"/fastq/{srr}_1.fastq.gz"
        r2 = f"/fastq/{srr}_2.fastq.gz"
        out1 = f"/output/trimmed/{srr}_1_trimmed.fastq.gz"
        out2 = f"/output/trimmed/{srr}_2_trimmed.fastq.gz"
        input_args  = f"-i {r1} -I {r2}"
        output_args = f"-o {out1} -O {out2}"
        trimmed = [out1, out2]
    else:
        r1 = f"/fastq/{srr}.fastq.gz"
        out1 = f"/output/trimmed/{srr}_trimmed.fastq.gz"
        input_args  = f"-i {r1}"
        output_args = f"-o {out1}"
        trimmed = [out1]

    (out_dir / "trimmed").mkdir(parents=True, exist_ok=True)
    (out_dir / "qc").mkdir(parents=True, exist_ok=True)

    cmd = (
        f"fastp "
        f"{input_args} {output_args} "
        f"--thread {threads} "
        f"--length_required {fq_cfg['min_length']} "
        f"--qualified_quality_phred {fq_cfg['qualified_quality']} "
        f"--unqualified_percent_limit {fq_cfg['unqualified_percent_limit']} "
        f"--json /output/qc/{srr}_fastp.json "
        f"--html /output/qc/{srr}_fastp.html"
    )

    mounts = make_mounts(cfg, refs, fastq_dir, out_dir)
    rc, _, _ = docker_run(cfg, cmd, mounts)
    return rc == 0, trimmed


# ════════════════════════════════════════════════════════════
#  PIPELINE BRANCHES
# ════════════════════════════════════════════════════════════

def run_rnaseq(srr: str, fastq_dir: Path, out_dir: Path,
               cfg: dict, refs: dict, trimmed_fq: list[str]) -> bool:
    """STAR alignment → featureCounts → delete BAM."""
    threads  = cfg["processing"]["threads"]
    star_idx = "/references/star_index"
    gtf      = "/references/" + refs["gtf"].split("/references/")[-1] \
               if "/references/" in refs["gtf"] else f"/references/GRCh38/{Path(refs['gtf']).name}"
    bam_out  = f"/output/star/{srr}"
    star_extra = cfg["processing"]["star"]["extra_args"]
    mounts   = make_mounts(cfg, refs, fastq_dir, out_dir)
    (out_dir / "star" / srr).mkdir(parents=True, exist_ok=True)
    (out_dir / "counts").mkdir(parents=True, exist_ok=True)

    # ── STAR alignment ──
    fq_str = " ".join(trimmed_fq)
    cmd = (
        f"STAR --runThreadN {threads} "
        f"--genomeDir {star_idx} "
        f"--readFilesIn {fq_str} "
        f"--readFilesCommand zcat "
        f"--outFileNamePrefix {bam_out}/ "
        f"--quantMode GeneCounts "
        f"{star_extra}"
    )
    rc, _, _ = docker_run(cfg, cmd, mounts)
    if rc != 0:
        log(f"    {srr}: STAR failed", "error"); return False

    # ── featureCounts ──
    bam_path = f"{bam_out}/Aligned.sortedByCoord.out.bam"
    fc_extra = cfg["processing"]["featurecounts"]["extra_args"]
    fc_feat  = cfg["processing"]["featurecounts"]["feature_type"]
    fc_attr  = cfg["processing"]["featurecounts"]["attribute"]
    is_paired = len(trimmed_fq) == 2

    cmd = (
        f"featureCounts "
        f"-T {threads} "
        f"-a /references/GRCh38/{Path(refs['gtf']).name} "
        f"-t {fc_feat} -g {fc_attr} "
        f"{fc_extra if is_paired else fc_extra.replace('-p --countReadPairs','')} "
        f"-o /output/counts/{srr}_counts.txt "
        f"{bam_path}"
    )
    rc, _, _ = docker_run(cfg, cmd, mounts)
    if rc != 0:
        log(f"    {srr}: featureCounts failed", "error"); return False

    # ── delete BAM (default) ──
    if not cfg["processing"]["keep_bam"]:
        bam_local = out_dir / "star" / srr / "Aligned.sortedByCoord.out.bam"
        if bam_local.exists():
            bam_local.unlink()
            log(f"    {srr}: BAM deleted")

    return True


def run_scrnaseq(srr: str, fastq_dir: Path, out_dir: Path,
                 cfg: dict, refs: dict, trimmed_fq: list[str]) -> bool:
    """STARsolo for 10x Chromium scRNA-seq."""
    threads   = cfg["processing"]["threads"]
    star_idx  = "/references/star_index"
    sc_cfg    = cfg["processing"]["starsolo"]
    whitelist = f"/references/{sc_cfg['whitelist']}"
    cb_len    = sc_cfg["cell_barcode_length"]
    umi_len   = sc_cfg["umi_length"]
    mounts    = make_mounts(cfg, refs, fastq_dir, out_dir)
    (out_dir / "starsolo" / srr).mkdir(parents=True, exist_ok=True)
    bam_out   = f"/output/starsolo/{srr}"

    # STARsolo expects R1 (barcode+UMI) then R2 (cDNA)
    # _1.fastq.gz = R1 (barcode), _2.fastq.gz = R2 (cDNA)
    if len(trimmed_fq) == 2:
        r1, r2 = trimmed_fq[0], trimmed_fq[1]
    else:
        r1 = r2 = trimmed_fq[0]

    cmd = (
        f"STAR --runThreadN {threads} "
        f"--genomeDir {star_idx} "
        f"--readFilesIn {r2} {r1} "
        f"--readFilesCommand zcat "
        f"--soloType CB_UMI_Simple "
        f"--soloCBwhitelist {whitelist} "
        f"--soloCBstart 1 --soloCBlen {cb_len} "
        f"--soloUMIstart {cb_len+1} --soloUMIlen {umi_len} "
        f"--outSAMtype BAM SortedByCoordinate "
        f"--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM "
        f"--outFileNamePrefix {bam_out}/ "
        f"--soloOutDir {bam_out}/Solo.out"
    )
    rc, _, _ = docker_run(cfg, cmd, mounts)
    if rc != 0:
        log(f"    {srr}: STARsolo failed", "error"); return False

    if not cfg["processing"]["keep_bam"]:
        bam_local = out_dir / "starsolo" / srr / "Aligned.sortedByCoordinate.out.bam"
        if bam_local.exists():
            bam_local.unlink()
            log(f"    {srr}: BAM deleted")

    return True


def run_atacseq(srr: str, fastq_dir: Path, out_dir: Path,
                cfg: dict, refs: dict, trimmed_fq: list[str],
                is_chip: bool = False) -> bool:
    """HISAT2 alignment → samtools sort/index → MACS3 peak calling → delete BAM."""
    threads    = cfg["processing"]["threads"]
    hisat2_idx = "/references/" + refs["hisat2_index"].split("/references/")[-1] \
                 if "/references/" in refs["hisat2_index"] \
                 else f"/references/hisat2_index/genome"
    h2_extra   = cfg["processing"]["hisat2"]["extra_args"]
    macs_cfg   = cfg["processing"]["macs3"]
    mounts     = make_mounts(cfg, refs, fastq_dir, out_dir)

    (out_dir / "hisat2").mkdir(parents=True, exist_ok=True)
    (out_dir / "peaks" / srr).mkdir(parents=True, exist_ok=True)

    bam_path = f"/output/hisat2/{srr}.sorted.bam"
    paired   = len(trimmed_fq) == 2

    # ── HISAT2 alignment ──
    fq_args = f"-1 {trimmed_fq[0]} -2 {trimmed_fq[1]}" if paired else f"-U {trimmed_fq[0]}"
    cmd = (
        f"hisat2 -p {threads} {h2_extra} "
        f"-x {hisat2_idx} {fq_args} "
        f"| samtools sort -@ {threads} -o {bam_path} && "
        f"samtools index {bam_path}"
    )
    rc, _, _ = docker_run(cfg, cmd, mounts)
    if rc != 0:
        log(f"    {srr}: HISAT2 failed", "error"); return False

    # ── MACS3 peak calling ──
    extra = macs_cfg["atac_extra"] if not is_chip else macs_cfg["chip_extra"]
    fmt   = "BAMPE" if paired else "BAM"
    cmd = (
        f"macs3 callpeak "
        f"-t {bam_path} "
        f"-f {fmt} "
        f"-g {macs_cfg['genome_size']} "
        f"--outdir /output/peaks/{srr} "
        f"-n {srr} "
        f"--qvalue {macs_cfg['fdr']} "
        f"{extra}"
    )
    rc, _, _ = docker_run(cfg, cmd, mounts)
    if rc != 0:
        log(f"    {srr}: MACS3 failed", "error"); return False

    if not cfg["processing"]["keep_bam"]:
        bam_local = out_dir / "hisat2" / f"{srr}.sorted.bam"
        bai_local = out_dir / "hisat2" / f"{srr}.sorted.bam.bai"
        for f in [bam_local, bai_local]:
            if f.exists(): f.unlink()
        log(f"    {srr}: BAM deleted")

    return True


def run_methylation(srr: str, fastq_dir: Path, out_dir: Path,
                    cfg: dict, refs: dict, trimmed_fq: list[str]) -> bool:
    """Bismark bisulfite alignment → methylation extraction."""
    threads    = cfg["processing"]["threads"]
    bismark_g  = "/references/GRCh38"   # directory with Bisulfite_Genome/ subdir
    mounts     = make_mounts(cfg, refs, fastq_dir, out_dir)
    paired     = len(trimmed_fq) == 2
    (out_dir / "bismark" / srr).mkdir(parents=True, exist_ok=True)

    fq_args = f"-1 {trimmed_fq[0]} -2 {trimmed_fq[1]}" if paired else trimmed_fq[0]
    cmd = (
        f"bismark --genome {bismark_g} "
        f"--parallel {max(1, threads//4)} "
        f"{fq_args} "
        f"-o /output/bismark/{srr}"
    )
    rc, _, _ = docker_run(cfg, cmd, mounts)
    if rc != 0:
        log(f"    {srr}: Bismark failed", "error"); return False

    # Extract methylation calls
    bam_glob_cmd = f"ls /output/bismark/{srr}/*.bam | head -1"
    rc2, bam_path, _ = docker_run(cfg, bam_glob_cmd, mounts)
    bam_path = bam_path.strip()
    if bam_path:
        cmd = (
            f"bismark_methylation_extractor "
            f"{'--paired-end' if paired else '--single-end'} "
            f"--genome_folder {bismark_g} "
            f"--CX_context "
            f"--bedGraph "
            f"-o /output/bismark/{srr} "
            f"{bam_path}"
        )
        docker_run(cfg, cmd, mounts)

    if not cfg["processing"]["keep_bam"]:
        bismark_dir = out_dir / "bismark" / srr
        for bam in bismark_dir.glob("*.bam"):
            bam.unlink()
        log(f"    {srr}: BAM deleted")

    return True


# ════════════════════════════════════════════════════════════
#  DISPATCH
# ════════════════════════════════════════════════════════════

def process_srr(run_entry: dict, cfg: dict, refs: dict,
                processed_dir: Path, fastq_tmp: Path,
                dry_run: bool = False) -> dict:
    """
    Process one SRR:
      1. fastp QC
      2. dispatch to correct alignment pipeline
      3. delete FASTQ (unless keep_fastq=True)
    """
    srr        = run_entry["srr"]
    gse        = run_entry.get("gse","unknown")
    gsm        = run_entry.get("gsm","unknown")
    omics_type = run_entry.get("omicsType","RNA-seq")
    fastq_dir  = Path(run_entry.get("fastq_dir", fastq_tmp / gse / gsm / srr))
    pipeline   = resolve_pipeline(omics_type)
    out_dir    = processed_dir / gse / srr
    out_dir.mkdir(parents=True, exist_ok=True)

    log(f"  {srr}: pipeline={pipeline}, omics={omics_type}")

    if dry_run:
        return {"srr": srr, "status": "dry_run"}

    # ── fastp ──
    ok, trimmed_fq = run_fastp(srr, fastq_dir, out_dir, cfg, refs)
    if not ok:
        return {"srr": srr, "status": "fastp_failed"}

    # ── alignment + counting ──
    if pipeline == "rnaseq":
        ok = run_rnaseq(srr, fastq_dir, out_dir, cfg, refs, trimmed_fq)
    elif pipeline == "scrnaseq":
        ok = run_scrnaseq(srr, fastq_dir, out_dir, cfg, refs, trimmed_fq)
    elif pipeline == "atacseq":
        ok = run_atacseq(srr, fastq_dir, out_dir, cfg, refs, trimmed_fq, is_chip=False)
    elif pipeline == "chipseq":
        ok = run_atacseq(srr, fastq_dir, out_dir, cfg, refs, trimmed_fq, is_chip=True)
    elif pipeline == "methylation":
        ok = run_methylation(srr, fastq_dir, out_dir, cfg, refs, trimmed_fq)
    else:
        log(f"    {srr}: unknown pipeline '{pipeline}'", "warn")
        ok = False

    # ── delete trimmed FASTQ ──
    if ok and not cfg["processing"].get("keep_fastq", False):
        trimmed_dir = out_dir / "trimmed"
        for f in trimmed_dir.glob("*.fastq.gz"):
            f.unlink()
        log(f"    {srr}: trimmed FASTQ deleted")

    # ── delete raw FASTQ ──
    if ok and not cfg["processing"].get("keep_fastq", False):
        for f in fastq_dir.glob("*.fastq.gz"):
            f.unlink()
        log(f"    {srr}: raw FASTQ deleted")

    status = "processed" if ok else "processing_failed"
    log(f"  {srr}: {status}", "success" if ok else "error")

    return {
        "srr":           srr,
        "status":        status,
        "pipeline":      pipeline,
        "output_dir":    str(out_dir),
        "processed_at":  datetime.now().isoformat(),
    }


# ════════════════════════════════════════════════════════════
#  COUNT MATRIX MERGER (RNA-seq)
# ════════════════════════════════════════════════════════════

def merge_count_matrices(gse: str, processed_dir: Path):
    """
    Merge individual featureCounts output files into one matrix per GSE.
    Output: processed/{gse}/{gse}_count_matrix.tsv
    """
    count_files = sorted((processed_dir / gse).rglob("*_counts.txt"))
    if not count_files:
        return

    import csv as csv_mod
    matrices = {}
    gene_ids = None

    for cf in count_files:
        srr = cf.stem.replace("_counts","")
        with open(cf) as f:
            reader = csv_mod.reader(f, delimiter="\t")
            rows = [r for r in reader if not r[0].startswith("#")]
        # featureCounts: col0=geneid, col5=count
        if not gene_ids:
            gene_ids = [r[0] for r in rows[1:]]
        matrices[srr] = {r[0]: r[-1] for r in rows[1:]}

    if not gene_ids:
        return

    out_path = processed_dir / gse / f"{gse}_count_matrix.tsv"
    srrs = list(matrices.keys())
    with open(out_path, "w") as f:
        f.write("gene_id\t" + "\t".join(srrs) + "\n")
        for gid in gene_ids:
            row = [gid] + [matrices[s].get(gid,"0") for s in srrs]
            f.write("\t".join(row) + "\n")

    log(f"Count matrix: {out_path}  ({len(gene_ids)} genes × {len(srrs)} samples)", "success")


# ════════════════════════════════════════════════════════════
#  MAIN RUNNER
# ════════════════════════════════════════════════════════════

def run(manifest_path: Path, cfg: dict,
        filter_gse: str = None, dry_run: bool = False):

    with open(manifest_path) as f:
        manifest = json.load(f)

    refs          = load_references(cfg)
    fastq_tmp     = Path(cfg["paths"]["fastq_tmp"]).resolve()
    processed_dir = Path(cfg["paths"]["processed_out"]).resolve()
    processed_dir.mkdir(parents=True, exist_ok=True)

    runs = manifest.get("runs", [])

    # Only process SRRs that have been downloaded
    ready = [r for r in runs
             if r.get("status") == "downloaded"
             and r.get("fastq_files")]
    if filter_gse:
        ready = [r for r in ready if r.get("gse") == filter_gse]

    log(f"{len(ready)} SRR runs ready to process")

    # Group by GSE so we can merge after all samples in a dataset are done
    by_gse: dict = defaultdict(list)
    for r in ready:
        by_gse[r.get("gse","unknown")].append(r)

    run_by_srr = {r["srr"]: r for r in manifest["runs"]}

    for gse, gse_runs in by_gse.items():
        log(f"\n── {gse} ({len(gse_runs)} runs) ──")
        all_ok = True
        for run_entry in gse_runs:
            result = process_srr(run_entry, cfg, refs,
                                  processed_dir, fastq_tmp, dry_run)
            if run_entry["srr"] in run_by_srr:
                run_by_srr[run_entry["srr"]].update(result)
            with open(manifest_path, "w") as f:
                json.dump(manifest, f, indent=2)
            if result["status"] != "processed":
                all_ok = False

        if all_ok and not dry_run:
            log(f"Merging count matrices for {gse}...")
            merge_count_matrices(gse, processed_dir)

    log("Processing complete.", "success")
    log(f"Output: {processed_dir}")


# ════════════════════════════════════════════════════════════
#  CLI
# ════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Align and count reads for each downloaded SRR run"
    )
    parser.add_argument("--manifest", "-m", required=True)
    parser.add_argument("--gse",  help="Process only this GSE accession")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--config",  default="pipeline_config.yaml")
    args = parser.parse_args()

    global CONFIG_FILE
    CONFIG_FILE = Path(args.config)
    cfg = load_config()

    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        print(f"Manifest not found: {manifest_path}"); sys.exit(1)

    try:
        run(manifest_path, cfg,
            filter_gse = args.gse,
            dry_run    = args.dry_run)
    except KeyboardInterrupt:
        print("\nInterrupted. Re-run to resume — completed SRRs are marked in the manifest.")
    except Exception as e:
        print(f"Error: {e}"); raise


if __name__ == "__main__":
    main()
