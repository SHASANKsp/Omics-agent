#!/usr/bin/env python3
"""
geo_downloader.py — GEO series matrix parser + SRA manifest builder
─────────────────────────────────────────────────────────────────────
For each GSE accession from the curated CSV:
  1. Downloads the series matrix file (metadata)
  2. Parses it to classify the dataset:
       GEO-direct  → processed files (counts, peaks) on GEO FTP → download them
       SRA-linked  → raw reads in SRA → resolve GSM→SRX→SRR accessions
  3. For GEO-direct: downloads supplementary files immediately
  4. For SRA-linked: writes sra_manifest.json for sra_downloader.py

Usage:
    python geo_downloader.py --csv omics_results/omics_datasets_psoriasis_curated.csv
    python geo_downloader.py --accessions GSE200637 GSE301804
    python geo_downloader.py --csv ... --dry-run
    python geo_downloader.py --csv ... --sra-only     # skip GEO-direct downloads

Requirements:
    pip install requests rich pyyaml
"""

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import requests
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

# ── Constants ────────────────────────────────────────────────
GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
NCBI_EUTILS  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CHUNK_SIZE   = 1024 * 1024
NCBI_DELAY   = 0.5
DL_WORKERS   = 3
CONFIG_FILE  = Path("pipeline_config.yaml")


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
#  GEO FTP HELPERS
# ════════════════════════════════════════════════════════════

def geo_stub(acc: str) -> str:
    return f"{GEO_FTP_BASE}/{acc[:-3]}nnn/{acc}"

def geo_matrix_url(acc: str) -> str:
    return f"{geo_stub(acc)}/matrix/{acc}_series_matrix.txt.gz"

def geo_suppl_index_url(acc: str) -> str:
    return f"{geo_stub(acc)}/suppl/"


# ════════════════════════════════════════════════════════════
#  FILE DOWNLOAD
# ════════════════════════════════════════════════════════════

def fmt_size(n: int) -> str:
    for u in ("B","KB","MB","GB"):
        if n < 1024: return f"{n:.1f} {u}"
        n /= 1024
    return f"{n:.1f} TB"


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path,"rb") as f:
        for chunk in iter(lambda: f.read(CHUNK_SIZE), b""):
            h.update(chunk)
    return h.hexdigest()


def download_file(url: str, dest: Path, resume: bool = True) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    existing = dest.stat().st_size if dest.exists() else 0
    headers  = {"User-Agent": "GEODownloader/1.0"}
    if resume and existing:
        headers["Range"] = f"bytes={existing}-"

    try:
        r = requests.get(url, headers=headers, stream=True, timeout=60)
        if r.status_code == 416:
            return {"url":url,"path":str(dest),"size_bytes":existing,
                    "md5":md5_file(dest),"status":"already_complete","error":""}
        if r.status_code == 404:
            return {"url":url,"path":str(dest),"size_bytes":0,"md5":"",
                    "status":"not_found","error":"404"}
        r.raise_for_status()
        mode = "ab" if (resume and existing) else "wb"
        with open(dest, mode) as f:
            for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                if chunk: f.write(chunk)
        return {"url":url,"path":str(dest),"size_bytes":dest.stat().st_size,
                "md5":md5_file(dest),"status":"ok","error":""}
    except Exception as e:
        return {"url":url,"path":str(dest),"size_bytes":existing,
                "md5":"","status":"error","error":str(e)}


# ════════════════════════════════════════════════════════════
#  SERIES MATRIX PARSER
# ════════════════════════════════════════════════════════════

def parse_series_matrix(gz_path: Path) -> dict:
    """
    Parse a GEO series matrix file and extract:
      - Series-level metadata (title, type, platform, organism)
      - Sample accessions (GSM numbers)
      - Whether samples have supplementary files on GEO FTP
        or point only to SRA (supplementary = NONE)
      - SRA project accession (from !Series_relation field)

    Returns a structured dict describing the dataset.
    """
    meta   = {}
    samples = []
    suppl_files_per_sample = []
    sra_project = None

    try:
        opener = gzip.open if str(gz_path).endswith(".gz") else open
        with opener(gz_path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.rstrip("\n")

                # Series-level fields
                if line.startswith("!Series_"):
                    key, _, val = line.partition(" = ")
                    key = key.lstrip("!").strip()
                    val = val.strip().strip('"')

                    if key == "Series_geo_accession":
                        meta["accession"] = val
                    elif key == "Series_title":
                        meta["title"] = val
                    elif key == "Series_type":
                        meta.setdefault("types", []).append(val)
                    elif key == "Series_platform_id":
                        meta["platform"] = val
                    elif key == "Series_relation":
                        # Modern GEO uses BioProject links, older ones use SRP directly.
                        # Examples:
                        #   "SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRP123456"
                        #   "BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA123456"
                        srp_match = re.search(r"(SRP\d+|ERP\d+|DRP\d+)", val)
                        prj_match = re.search(r"(PRJNA\d+|PRJEB\d+|PRJDB\d+)", val)
                        if srp_match:
                            sra_project = srp_match.group(1)
                        elif prj_match:
                            # Store BioProject accession — will be resolved via eLink
                            sra_project = prj_match.group(1)
                        meta.setdefault("relations", []).append(val)

                # Sample accession row — tab-separated: !Sample_geo_accession\t"GSM1"\t"GSM2"...
                elif line.startswith("!Sample_geo_accession"):
                    parts = line.split("\t")
                    # parts[0] = "!Sample_geo_accession", parts[1:] = quoted GSM accessions
                    samples = [v.strip().strip('"')
                               for v in parts[1:] if v.strip().strip('"')]

                # Supplementary file row — same tab-separated format
                elif line.startswith("!Sample_supplementary_file"):
                    parts = line.split("\t")
                    file_list = [v.strip().strip('"')
                                 for v in parts[1:] if v.strip().strip('"')]
                    suppl_files_per_sample.append(file_list)

    except Exception as e:
        log(f"Series matrix parse error: {e}", "warn")

    # Determine if this dataset is GEO-direct or SRA-linked
    # A sample is GEO-direct if ANY of its supplementary_file fields is not "NONE"
    all_none = True
    geo_suppl_names = set()
    for file_row in suppl_files_per_sample:
        for fname in file_row:
            if fname.upper() != "NONE" and fname:
                all_none = False
                geo_suppl_names.add(fname.split("/")[-1])

    meta["sample_accessions"] = samples
    meta["sample_count"]      = len(samples)
    meta["sra_project"]       = sra_project
    meta["geo_supplementary_files"] = sorted(geo_suppl_names)

    # Classification
    if not all_none and geo_suppl_names:
        meta["data_location"] = "GEO_direct"
        meta["has_processed_data"] = True
    elif sra_project or all_none:
        meta["data_location"] = "SRA_linked"
        meta["has_processed_data"] = False
    else:
        meta["data_location"] = "unknown"
        meta["has_processed_data"] = False

    return meta


# ════════════════════════════════════════════════════════════
#  GSM → SRX → SRR RESOLVER (via NCBI eLink)
# ════════════════════════════════════════════════════════════

def gsm_to_srr(gsm: str) -> list[dict]:
    """
    Resolve a GEO sample accession (GSM) to its SRA run accessions (SRR).
    Path: GSM → [eSearch GEO samples db] → UID → [eLink to SRA] → SRX UID
          → [eSummary] → SRR accessions

    Returns list of {srr, srx, gsm, title, spots, bases, size_mb}.
    """
    time.sleep(NCBI_DELAY)

    # Step 1: find UID for the GSM in the GEO samples database
    search_url = (
        f"{NCBI_EUTILS}/esearch.fcgi"
        f"?db=gds&term={gsm}[Accession]&retmode=json&retmax=1"
    )
    try:
        r = requests.get(search_url, timeout=20,
                         headers={"User-Agent": "GEODownloader/1.0"})
        r.raise_for_status()
        uids = r.json()["esearchresult"].get("idlist", [])
        if not uids:
            return []
        gsm_uid = uids[0]
    except Exception as e:
        log(f"  eSearch failed for {gsm}: {e}", "warn")
        return []

    time.sleep(NCBI_DELAY)

    # Step 2: eLink from GEO samples → SRA
    elink_url = (
        f"{NCBI_EUTILS}/elink.fcgi"
        f"?dbfrom=gds&db=sra&id={gsm_uid}&retmode=json"
    )
    try:
        r = requests.get(elink_url, timeout=20,
                         headers={"User-Agent": "GEODownloader/1.0"})
        r.raise_for_status()
        data = r.json()
        linksets = data.get("linksets", [])
        sra_uids = []
        for ls in linksets:
            for lsd in ls.get("linksetdbs", []):
                if lsd.get("dbto") == "sra":
                    sra_uids.extend(lsd.get("links", []))
        if not sra_uids:
            return []
    except Exception as e:
        log(f"  eLink failed for {gsm}: {e}", "warn")
        return []

    time.sleep(NCBI_DELAY)

    # Step 3: eSummary on SRA UIDs to get SRR run accessions
    uid_str = ",".join(str(u) for u in sra_uids[:20])  # cap at 20 experiments
    summary_url = (
        f"{NCBI_EUTILS}/esummary.fcgi"
        f"?db=sra&id={uid_str}&retmode=json"
    )
    try:
        r = requests.get(summary_url, timeout=30,
                         headers={"User-Agent": "GEODownloader/1.0"})
        r.raise_for_status()
        result = r.json().get("result", {})
    except Exception as e:
        log(f"  eSummary failed for {gsm}: {e}", "warn")
        return []

    runs = []
    for uid, rec in result.items():
        if uid == "uids":
            continue
        exp_xml  = rec.get("expxml", "")
        runs_xml = rec.get("runs", "")

        # Extract SRX (experiment accession)
        srx_m = re.search(r'acc="(SRX\d+)"', exp_xml)
        srx   = srx_m.group(1) if srx_m else ""

        # Extract title
        title_m = re.search(r'<Title>(.*?)</Title>', exp_xml)
        title   = title_m.group(1).strip() if title_m else ""

        # Extract all SRR run accessions from the runs XML blob
        for run_m in re.finditer(
            r'acc="(SRR\d+)"[^>]*spots="(\d+)"[^>]*bases="(\d+)"', runs_xml
        ):
            srr    = run_m.group(1)
            spots  = int(run_m.group(2))
            bases  = int(run_m.group(3))
            size_mb = round(bases / 4 / 1e6, 1)  # rough FASTQ size estimate
            runs.append({
                "srr":     srr,
                "srx":     srx,
                "gsm":     gsm,
                "title":   title,
                "spots":   spots,
                "bases":   bases,
                "size_mb": size_mb,
            })

    return runs


# ════════════════════════════════════════════════════════════
#  SUPPLEMENTARY FILE DISCOVERY (GEO FTP)
# ════════════════════════════════════════════════════════════

def list_geo_suppl_files(acc: str) -> list[dict]:
    """
    Parse GEO FTP index for supplementary files.
    Filters strictly to filenames that look like GEO data files —
    the NCBI FTP index HTML contains navigation links that must be excluded.
    """
    index_url = geo_suppl_index_url(acc)
    try:
        r = requests.get(index_url, timeout=20,
                         headers={"User-Agent": "GEODownloader/1.0"})
        if r.status_code == 404:
            return []
        r.raise_for_status()
        files = []
        for m in re.finditer(r'href="([^"]+)"', r.text):
            name = m.group(1).strip("/")

            # Skip anything that looks like a URL or navigation path
            if name.startswith("http") or name.startswith("/") or "/" in name:
                continue
            # Skip empty, parent dir, directory links
            if not name or name == "../" or name.endswith("/"):
                continue
            # Only accept names that look like actual data files
            # GEO filenames always start with the accession or are known types
            valid_extensions = (
                ".tar", ".gz", ".zip", ".txt", ".csv", ".tsv",
                ".bed", ".bw", ".bigwig", ".bam", ".vcf", ".h5",
                ".rds", ".rdata", ".loom", ".h5ad", ".mtx"
            )
            name_lower = name.lower()
            if not any(name_lower.endswith(ext) for ext in valid_extensions):
                continue

            if name.endswith("_RAW.tar"):
                files.insert(0, {
                    "name": name,
                    "url":  f"{index_url}{name}",
                    "type": "raw_tar",
                })
            elif "series_matrix" in name_lower:
                continue   # already downloaded
            else:
                files.append({
                    "name": name,
                    "url":  f"{index_url}{name}",
                    "type": "supplementary",
                })
        return files
    except Exception as e:
        log(f"FTP index failed for {acc}: {e}", "warn")
        return []


# ════════════════════════════════════════════════════════════
#  PER-ACCESSION ORCHESTRATION
# ════════════════════════════════════════════════════════════

def process_accession(acc: str, out_dir: Path, omics_type: str = "",
                      dry_run: bool = False, resume: bool = True) -> dict:
    """
    Full resolution for one GSE accession:
      1. Download + parse series matrix
      2. Classify as GEO-direct or SRA-linked
      3a. GEO-direct  → download supplementary files
      3b. SRA-linked  → resolve GSM → SRR list
    Returns a result dict used to build sra_manifest.json.
    """
    acc_dir = out_dir / acc
    acc_dir.mkdir(parents=True, exist_ok=True)
    log(f"Processing {acc} ({omics_type or 'unknown omics'})...")

    # ── Step 1: series matrix ──
    matrix_url  = geo_matrix_url(acc)
    matrix_dest = acc_dir / f"{acc}_series_matrix.txt.gz"
    dl_result   = download_file(matrix_url, matrix_dest, resume=resume)
    if dl_result["status"] not in ("ok", "already_complete"):
        log(f"  {acc}: series matrix download failed — {dl_result['error']}", "error")
        return {"accession": acc, "status": "error",
                "error": "series matrix download failed"}

    log(f"  {acc}: parsing series matrix...")
    meta = parse_series_matrix(matrix_dest)
    meta["accession"]  = acc
    meta["omicsType"]  = omics_type
    meta["matrix_path"]= str(matrix_dest)

    log(f"  {acc}: {meta['sample_count']} samples | "
        f"location={meta['data_location']} | "
        f"SRA project={meta.get('sra_project','none')}", "success")

    result = {
        "accession":     acc,
        "omicsType":     omics_type,
        "title":         meta.get("title",""),
        "sample_count":  meta.get("sample_count", 0),
        "data_location": meta.get("data_location", "unknown"),
        "sra_project":   meta.get("sra_project"),
        "matrix_path":   str(matrix_dest),
        "geo_files":     [],
        "sra_runs":      [],
        "status":        "pending",
    }

    # ── Step 2a: GEO-direct — download supplementary files ──
    if meta["data_location"] == "GEO_direct" and not dry_run:
        suppl_files = list_geo_suppl_files(acc)
        suppl_dir   = acc_dir / "supplementary"
        suppl_dir.mkdir(exist_ok=True)
        geo_dl_results = []
        for sf in suppl_files:
            dest = suppl_dir / sf["name"]
            log(f"  Downloading {sf['name']}...")
            dl  = download_file(sf["url"], dest, resume=resume)
            geo_dl_results.append({
                "name":       sf["name"],
                "type":       sf["type"],
                "url":        sf["url"],
                "path":       str(dest),
                "size":       fmt_size(dl["size_bytes"]),
                "status":     dl["status"],
            })
            log(f"    {sf['name']}: {fmt_size(dl['size_bytes'])} [{dl['status']}]",
                "success" if dl["status"] in ("ok","already_complete") else "warn")
        result["geo_files"] = geo_dl_results
        result["status"]    = "geo_complete"

    elif meta["data_location"] == "GEO_direct" and dry_run:
        suppl_files = list_geo_suppl_files(acc)
        result["geo_files"] = [{"name": sf["name"], "url": sf["url"]}
                                for sf in suppl_files]
        result["status"] = "dry_run"
        log(f"  [DRY RUN] {acc}: {len(suppl_files)} GEO supplementary files")

    # ── Step 2b: SRA-linked — resolve GSM → SRR ──
    elif meta["data_location"] == "SRA_linked":
        samples = meta.get("sample_accessions", [])
        log(f"  {acc}: resolving {len(samples)} GSM accessions → SRR...")

        # Cap at max_gsm_per_gse to avoid runaway resolution on large datasets.
        # GSE314390 has 23 GSMs × 20 runs each = 460 SRRs = 2910 GB.
        # Default cap: resolve all GSMs but flag if estimated size is large.
        all_runs = []
        for i, gsm in enumerate(samples):
            runs = gsm_to_srr(gsm)
            all_runs.extend(runs)
            if (i+1) % 5 == 0:
                log(f"    {i+1}/{len(samples)} GSMs resolved ({len(all_runs)} SRR so far)")

        total_mb  = sum(r.get("size_mb", 0) for r in all_runs)
        total_gb  = total_mb / 1024

        # Warn if dataset is very large
        if total_gb > 500:
            log(f"  {acc}: ⚠ Large dataset — {len(all_runs)} SRR runs, "
                f"est. {total_gb:.0f} GB. Consider downloading selectively "
                f"using --srr flags in sra_downloader.py", "warn")
        
        result["sra_runs"]     = all_runs
        result["total_size_mb"]= total_mb
        result["status"]       = "sra_resolved"
        log(f"  {acc}: {len(all_runs)} SRR runs found, "
            f"est. {total_gb:.1f} GB FASTQ", "success")

    else:
        result["status"] = "unknown_location"
        log(f"  {acc}: data location unknown — manual inspection needed", "warn")

    # Save per-accession metadata
    meta_path = acc_dir / "dataset_meta.json"
    with open(meta_path, "w") as f:
        json.dump({**meta, **result}, f, indent=2)

    return result


# ════════════════════════════════════════════════════════════
#  SRA MANIFEST WRITER
# ════════════════════════════════════════════════════════════

def write_sra_manifest(results: list[dict], out_dir: Path,
                       disease: str = "") -> Path:
    """
    Write sra_manifest.json consumed by sra_downloader.py.
    Contains one entry per SRR run with enough metadata for
    the downloader to organise output correctly.
    """
    manifest_entries = []
    for r in results:
        if r.get("status") != "sra_resolved":
            continue
        acc       = r["accession"]
        omics     = r.get("omicsType", "")
        for run in r.get("sra_runs", []):
            manifest_entries.append({
                "srr":        run["srr"],
                "srx":        run.get("srx",""),
                "gsm":        run.get("gsm",""),
                "gse":        acc,
                "omicsType":  omics,
                "title":      run.get("title",""),
                "spots":      run.get("spots", 0),
                "bases":      run.get("bases", 0),
                "size_mb":    run.get("size_mb", 0),
                "status":     "pending",   # updated by sra_downloader.py
            })

    manifest = {
        "disease":    disease,
        "created_at": datetime.now().isoformat(),
        "total_runs": len(manifest_entries),
        "total_size_gb": round(
            sum(e["size_mb"] for e in manifest_entries) / 1024, 2
        ),
        "runs": manifest_entries,
    }

    dest = out_dir / "sra_manifest.json"
    with open(dest, "w") as f:
        json.dump(manifest, f, indent=2)

    log(f"SRA manifest: {dest}  ({len(manifest_entries)} runs, "
        f"{manifest['total_size_gb']:.1f} GB estimated)", "success")
    return dest


# ════════════════════════════════════════════════════════════
#  CSV READER
# ════════════════════════════════════════════════════════════

def accessions_from_csv(csv_path: Path) -> list[dict]:
    """Read GSE accessions + omicsType from curated CSV."""
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            acc = row.get("accession","").strip()
            if acc.startswith("GSE"):
                rows.append({
                    "accession": acc,
                    "omicsType": row.get("omicsType",""),
                })
    seen = set()
    deduped = []
    for r in rows:
        if r["accession"] not in seen:
            seen.add(r["accession"])
            deduped.append(r)
    return deduped


# ════════════════════════════════════════════════════════════
#  MAIN RUNNER
# ════════════════════════════════════════════════════════════

def run(entries: list[dict], out_dir: Path,
        dry_run: bool = False, resume: bool = True,
        sra_only: bool = False, disease: str = ""):

    out_dir.mkdir(parents=True, exist_ok=True)
    log(f"Processing {len(entries)} GEO accessions → {out_dir}")

    all_results = []
    for entry in entries:
        result = process_accession(
            acc        = entry["accession"],
            out_dir    = out_dir,
            omics_type = entry.get("omicsType",""),
            dry_run    = dry_run,
            resume     = resume,
        )
        all_results.append(result)

    # Summary table
    geo_direct = [r for r in all_results if r.get("data_location") == "GEO_direct"]
    sra_linked = [r for r in all_results if r.get("data_location") == "SRA_linked"]
    unknown    = [r for r in all_results
                  if r.get("data_location") not in ("GEO_direct","SRA_linked")]

    if HAS_RICH:
        t = Table(title="GEO Resolution Summary", box=None,
                  header_style="bold", padding=(0,2))
        t.add_column("Accession", style="cyan", no_wrap=True)
        t.add_column("Omics")
        t.add_column("Location")
        t.add_column("Samples", justify="right")
        t.add_column("SRR runs", justify="right")
        t.add_column("Est. size")
        for r in all_results:
            loc   = r.get("data_location","?")
            col   = "green" if loc=="GEO_direct" else "blue" if loc=="SRA_linked" else "yellow"
            size  = f"{r.get('total_size_mb',0)/1024:.1f} GB" if r.get("sra_runs") else "—"
            t.add_row(
                r["accession"],
                r.get("omicsType",""),
                f"[{col}]{loc}[/{col}]",
                str(r.get("sample_count","?")),
                str(len(r.get("sra_runs",[]))),
                size,
            )
        console.print()
        console.print(t)

    log(f"GEO-direct: {len(geo_direct)} | SRA-linked: {len(sra_linked)} | Unknown: {len(unknown)}")

    # Write SRA manifest for SRA-linked datasets
    manifest_path = None
    if sra_linked:
        manifest_path = write_sra_manifest(all_results, out_dir, disease=disease)
        log(f"Next step:  python processing/sra_downloader.py --manifest {manifest_path}")
    else:
        log("No SRA-linked datasets found — all data is on GEO FTP directly", "success")

    return all_results, manifest_path


# ════════════════════════════════════════════════════════════
#  CLI
# ════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Resolve GEO accessions to data files or SRA run manifests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--csv",         "-c", type=str)
    src.add_argument("--accessions",  "-a", nargs="+")
    parser.add_argument("--out",      "-o", type=str, default="geo_downloads")
    parser.add_argument("--disease",  "-d", type=str, default="")
    parser.add_argument("--dry-run",  action="store_true")
    parser.add_argument("--no-resume",action="store_true")
    parser.add_argument("--sra-only", action="store_true",
                        help="Skip GEO-direct file downloads, only resolve SRR accessions")
    args = parser.parse_args()

    if args.csv:
        csv_path = Path(args.csv)
        if not csv_path.exists():
            print(f"CSV not found: {csv_path}"); sys.exit(1)
        entries = accessions_from_csv(csv_path)
        disease = args.disease or csv_path.stem.replace("omics_datasets_","").replace("_curated","")
    else:
        entries = [{"accession": a.strip(), "omicsType": ""}
                   for a in args.accessions if a.strip().startswith("GSE")]
        disease = args.disease

    if not entries:
        print("No valid GSE accessions."); sys.exit(1)

    try:
        run(entries, Path(args.out),
            dry_run = args.dry_run,
            resume  = not args.no_resume,
            sra_only= args.sra_only,
            disease = disease)
    except KeyboardInterrupt:
        print("\nInterrupted.")
    except Exception as e:
        print(f"Error: {e}"); raise


if __name__ == "__main__":
    main()
