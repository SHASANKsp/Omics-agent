# Omics-pipelines

An end-to-end, disease-driven pipeline that goes from a plain disease name to
processed omics results (gene count matrices, peaks, and methylation calls).
Given a disease such as `psoriasis`, it uses a **local LLM to reason about which
omics assays matter**, queries **real public databases** (NCBI GEO, NCBI SRA,
EBI PRIDE) for matching datasets, scores them for relevance, resolves them down
to individual sequencing runs, downloads the raw FASTQ, and runs the appropriate
alignment/quantification workflow — all inside a single reproducible Docker
image.

Two ideas run through the whole design:

- **The LLM reasons, it never invents data.** Accession IDs always come from
  live NCBI/EBI API responses, never from the model. The model is used only to
  pick omics types, generate search terms, and score relevance. This keeps the
  output verifiable and free of hallucinated GEO/SRA IDs.
- **Disk is treated as the scarce resource.** Runs are processed one SRR at a
  time and intermediate BAM/FASTQ files are deleted immediately after counting,
  so peak disk usage stays close to a single sample's footprint rather than the
  whole cohort's.

---

## Pipeline at a glance

```
                        pipeline_config.yaml  (single source of truth)
                                     │
   disease name ──▶ [1] omics_agent.py ──▶ curated CSV of GEO/SRA datasets
                                     │
                          [2] geo_downloader.py
                          parse series matrix, classify each dataset:
                            • GEO-direct  → download processed files from GEO FTP
                            • SRA-linked  → resolve GSM → SRX → SRR, write manifest
                                     │
                          [3] sra_downloader.py
                          prefetch → fasterq-dump → gzip  (one SRR at a time)
                                     │
                          [4] processor.py
                          fastp QC → aligner by assay type → counts/peaks/meth
                          delete BAM + FASTQ, merge per-sample counts per GSE
                                     │
                                     ▼
                     processed/<GSE>/<GSE>_count_matrix.tsv  (+ peaks, methylation, QC)
```

`run_pipeline.py` chains all four steps and can resume from any point after an
interruption.

---

## Repository layout

```
Omics-pipelines/
├── pipeline_config.yaml        # single source of truth for all settings
├── run_pipeline.py             # orchestrator: curate → geo → sra → process (with resume)
│
├── data_curation/
│   ├── omics_agent.py          # LLM-guided dataset discovery + relevance scoring (Ollama + NCBI/EBI APIs)
│   └── geo_downloader.py       # GEO series-matrix parser, GSM→SRR resolver, SRA manifest builder
│
├── processing/
│   ├── Dockerfile              # reproducible image with all bioinformatics tools
│   ├── reference_setup.py      # one-time GRCh38 reference + index setup
│   ├── sra_downloader.py       # prefetch → fasterq-dump → gzip, per SRR
│   └── processor.py            # fastp + aligners (STAR/STARsolo/HISAT2/Bismark) + counting/peaks
│
├── references/                 # genome, GTF, and prebuilt indexes (large; not versioned)
│   ├── references.json         # machine-written manifest of resolved reference paths
│   ├── GRCh38/                 # genome FASTA, GTF, Bismark bisulfite genome
│   ├── star_index/             # STAR / STARsolo index
│   └── hisat2_index/           # HISAT2 (grch38_tran) index
│
└── omics_results/              # curation outputs (example run: psoriasis)
    ├── omics_datasets_psoriasis_curated.csv   # relevant datasets (score ≥ 50)
    ├── omics_datasets_psoriasis_all.csv       # full audit trail incl. dropped datasets
    └── omics_report_psoriasis.json            # strategy + scores + full report
```

---

## The four stages in detail

### 1. Dataset curation — `data_curation/omics_agent.py`

Talks to a **local Ollama model** (default `mistral`; `llama3.2` and others work)
purely for reasoning:

1. **Strategy.** The model proposes a disease category, primary tissue,
   mechanistic focus, and 3–5 relevant omics types, each with up to six diverse
   NCBI search terms. A compact fallback prompt is used automatically if the
   model truncates its response.
2. **Live fetch.** For each omics type the tool queries the real APIs —
   **NCBI GEO** (`esearch`/`esummary` on the `gds` database), **NCBI SRA**
   (resolved to study-level `SRP` accessions to avoid per-run duplicates), and
   **EBI PRIDE** for proteomics. Human-only filtering (`"Homo sapiens"[Organism]`)
   is enforced at the query level, not left to the model. NCBI calls are
   rate-limited (0.8 s spacing) with automatic 429/5xx retry; summary fetches
   run in a small thread pool.
3. **Relevance scoring.** A fast text pre-filter drops records whose title and
   summary never mention the disease or a known synonym (no LLM call needed).
   Surviving candidates are scored 0–100 by the model in parallel batches, with
   hard drops for wrong-disease or non-human records.

Outputs three files per disease: a curated CSV (score ≥ 50), an audit CSV of
**all** fetched datasets with drop reasons, and a full JSON report. Nothing is
silently discarded.

```bash
python data_curation/omics_agent.py --disease "psoriasis" --model llama3.2
python data_curation/omics_agent.py --disease "ALS" --model mistral --max-results 20
python data_curation/omics_agent.py --list-models
```

### 2. GEO resolution — `data_curation/geo_downloader.py`

For each `GSE` accession in the curated CSV it downloads and parses the GEO
**series matrix**, then classifies the dataset:

- **GEO-direct** (samples have supplementary files on the GEO FTP) — the
  supplementary files (`_RAW.tar`, processed matrices, etc.) are downloaded
  directly, with resume support and MD5 checks.
- **SRA-linked** (samples point only to raw reads) — each `GSM` is resolved to
  its `SRR` runs via NCBI `eLink` (`GSM → GEO UID → SRA → SRX → SRR`), with FASTQ
  size estimates and a warning for very large (>500 GB) cohorts.

SRA-linked results are written to `sra_manifest.json`, the contract consumed by
the next two stages.

```bash
python data_curation/geo_downloader.py --csv omics_results/omics_datasets_psoriasis_curated.csv
python data_curation/geo_downloader.py --accessions GSE200637 GSE301804
python data_curation/geo_downloader.py --csv ... --dry-run
```

### 3. SRA download — `processing/sra_downloader.py`

Reads the manifest and, for each pending `SRR`, runs `prefetch → fasterq-dump →
pigz` inside the Docker container. Runs are downloaded **one at a time** to cap
peak disk usage, the `.sra` cache is deleted after compression, and the manifest
is updated after every run so an interrupted download resumes cleanly.

```bash
python processing/sra_downloader.py --manifest geo_downloads/sra_manifest.json
python processing/sra_downloader.py --manifest ... --srr SRR12345 SRR12346
python processing/sra_downloader.py --manifest ... --dry-run
```

### 4. Processing — `processing/processor.py`

Each downloaded run goes through `fastp` QC/trimming and is then dispatched by
its omics type:

| Assay type            | Workflow                                             | Main output                       |
|-----------------------|------------------------------------------------------|-----------------------------------|
| Bulk RNA-seq          | `fastp → STAR → featureCounts`                       | per-GSE gene count matrix         |
| scRNA-seq             | `fastp → STARsolo` (10x Chromium v3)                 | cell × gene matrices (`Solo.out`) |
| ATAC-seq              | `fastp → HISAT2 → MACS3` (`--nomodel --shift ...`)   | peak calls                        |
| ChIP-seq              | `fastp → HISAT2 → MACS3`                              | peak calls                        |
| Methylation / bisulfite | `fastp → Bismark → methylation extractor`          | methylation calls / bedGraph      |

Paired- vs single-end is auto-detected. BAM and FASTQ are deleted after counting
by default (`keep_bam` / `keep_fastq` in the config). Once every run in a `GSE`
is processed, the per-sample `featureCounts` outputs are merged into a single
`processed/<GSE>/<GSE>_count_matrix.tsv`.

```bash
python processing/processor.py --manifest geo_downloads/sra_manifest.json
python processing/processor.py --manifest ... --gse GSE314390   # one dataset only
python processing/processor.py --manifest ... --dry-run
```

---

## Requirements

- **Docker** (running) — all bioinformatics tools run inside the image, so
  nothing except Python and Docker needs to be installed on the host. This is
  what makes the pipeline work identically on Windows, macOS, and Linux.
- **Ollama** running locally (`ollama serve`) with at least one model pulled
  (`ollama pull mistral` or `ollama pull llama3.2`) — used only by stage 1.
- **Python 3.10+** on the host with: `pip install requests rich pyyaml`
- **Disk**: references alone are ~48 GB; raw data adds more per cohort (though
  it is deleted as it is processed).
- **RAM**: 16 GB is enough — index builds run inside Docker and adapt to
  available memory (STAR build takes ~1.5–2 h on 16 GB).

---

## Setup

**1. Build the Docker image** (SRA Toolkit, fastp, STAR/STARsolo, HISAT2,
featureCounts, MACS3, samtools, Bismark, bowtie2):

```bash
docker build -t omics-pipeline:1.0 processing/
```

**2. Build the reference genome and indexes** (one-time, ~4–5 hours, mostly
unattended). Everything is downloaded and built automatically — GRCh38 genome +
GTF from Ensembl (release 111), the STAR index, the prebuilt HISAT2 `grch38_tran`
index from AWS, the Bismark bisulfite genome, and the 10x Chromium v3 barcode
whitelist:

```bash
python processing/reference_setup.py            # full setup
python processing/reference_setup.py --check    # verify what already exists
python processing/reference_setup.py --genome-only
```

This writes `references/references.json`, which the processing stage reads to
locate each index.

---

## Running the full pipeline

`run_pipeline.py` chains all four stages and tracks progress in a per-disease
state file (`.pipeline_state_<disease>.json`), so any interrupted run can resume:

```bash
# Full pipeline
python run_pipeline.py --disease "psoriasis" --model llama3.2

# Resume from a specific step (earlier steps skipped if already done)
python run_pipeline.py --disease "psoriasis" --from-step geo
python run_pipeline.py --disease "psoriasis" --from-step sra
python run_pipeline.py --disease "psoriasis" --from-step process

# Preview without downloading/processing
python run_pipeline.py --disease "psoriasis" --model llama3.2 --dry-run

# Skip the reference check if you know it's already set up
python run_pipeline.py --disease "psoriasis" --skip-ref-check
```

Before running, the orchestrator checks that `references.json` and the Docker
image are present, and tells you the exact command to fix either if not.

---

## Configuration

All paths, thread counts, tool parameters, and cleanup policy live in
**`pipeline_config.yaml`** — every script reads from it, so it is the one file
to edit. Highlights:

- `docker.image` — image reference used at runtime (`omics-pipeline:1.0`).
- `paths.*` — where each stage reads and writes.
- `reference.*` — Ensembl release, assembly, and index paths (filled in by
  `reference_setup.py`).
- `processing.threads` — CPU threads for STAR, HISAT2, featureCounts.
- `processing.keep_bam` / `keep_fastq` — set `true` to retain intermediates
  (defaults `false` to save disk).
- Per-tool parameter blocks: `fastp`, `star`, `starsolo`, `hisat2`, `macs3`,
  `featurecounts`, `bismark`.
- `sra.*` — parallel prefetch jobs, cache dir, `fasterq-dump` threads.

---

## Example run (included)

`omics_results/` contains a completed curation run for **psoriasis**: a curated
dataset list, the full audit CSV with drop reasons, and the JSON report — useful
as a reference for the CSV/JSON schema the downstream stages expect.

---

## Third-party tools and licensing

The Docker image bundles established open-source bioinformatics tools —
SRA Toolkit (NCBI), fastp, STAR/STARsolo, HISAT2, Subread/featureCounts, MACS3,
samtools, bowtie2, and Bismark — plus Ollama and its models for the reasoning
step. Each carries its own license (several tool binaries and models have terms
that differ for academic vs. commercial use). Before using this pipeline in a
commercial setting, verify the license of every bundled tool and of the specific
Ollama model you run, and confirm that redistribution of the reference indexes
(e.g. the prebuilt HISAT2 index and the 10x barcode whitelist) is permitted for
your use case.
