#!/usr/bin/env python3
"""
geo_fetch.py — Download data from GEO for a given GSE ID
────────────────────────────────────────────────────────────
A standalone, self-contained downloader. The only third-party
dependency is `requests` (pip install requests). No SRA Toolkit,
no Docker, no config file.

For each GSE accession it can fetch four kinds of data:

  1. Series matrix          — GSExxx_series_matrix.txt.gz (sample metadata)
  2. Series supplementary   — files on the GEO FTP for the series
                              (e.g. GSExxx_RAW.tar, processed count matrices)
  3. Per-sample (GSM) suppl — individual per-sample supplementary files
  4. Raw reads (FASTQ)      — resolved to SRA runs and downloaded from the
                              EBI ENA mirror over HTTPS (with md5 checks)

Raw FASTQ can be very large (often 10s–100s of GB), so it is downloaded only
when you pass --raw (or --all), and you are shown the total size first.

──────────────────────────────────────────────────────────────
Usage
──────────────────────────────────────────────────────────────
  # Metadata + all supplementary files (the light default)
  python geo_fetch.py GSE200637

  # Everything, including raw FASTQ (asks for size confirmation)
  python geo_fetch.py GSE200637 --all

  # Just see what's available and how big it is — download nothing
  python geo_fetch.py GSE200637 --all --list

  # Pick exactly what you want
  python geo_fetch.py GSE200637 --matrix --raw
  python geo_fetch.py GSE200637 --series-suppl --sample-suppl

  # Several series at once, custom output dir, skip the raw-size prompt
  python geo_fetch.py GSE200637 GSE301804 -o ./downloads --all --yes

Options
  -o, --outdir DIR    Output directory (default: ./geo_downloads)
  --matrix            Download the series matrix (metadata)
  --series-suppl      Download series-level supplementary files
  --sample-suppl      Download per-sample (GSM) supplementary files
  --raw               Download raw FASTQ reads (heavy; via ENA)
  --all               Shorthand for all four of the above
  --skip-raw-tar      Skip the series-level GSExxx_RAW.tar (redundant with per-sample files)
  --list              Dry run: list what would be downloaded, then exit
  --max-runs N        Cap the number of raw FASTQ runs (default: no cap)
  --no-md5            Skip md5 verification of downloaded FASTQ files
  --yes               Don't prompt before downloading large raw data
  --retries N         HTTP retry attempts per request (default: 4)

If no content flag (--matrix/--series-suppl/--sample-suppl/--raw/--all) is
given, the default is: matrix + series-suppl + sample-suppl (no raw reads).
"""

import argparse
import gzip
import hashlib
import re
import sys
import time
from datetime import datetime
from pathlib import Path

import requests

# ── Endpoints ────────────────────────────────────────────────
GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo"
NCBI_EUTILS  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA_PORTAL   = "https://www.ebi.ac.uk/ena/portal/api/filereport"

CHUNK_SIZE   = 1024 * 1024        # 1 MB streaming chunks
NCBI_DELAY   = 0.4                # polite spacing between NCBI calls (sec)
UA           = {"User-Agent": "geo_fetch/1.0 (standalone GEO downloader)"}
TIMEOUT      = 60


# ════════════════════════════════════════════════════════════
#  LOGGING (stdlib only — no rich dependency)
# ════════════════════════════════════════════════════════════

def log(msg: str, status: str = "info"):
    ts  = datetime.now().strftime("%H:%M:%S")
    sym = {"info": "·", "ok": "✓", "warn": "!", "err": "✗"}.get(status, "·")
    print(f"{ts}  {sym}  {msg}")


def human(n: float) -> str:
    """Human-readable byte size."""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"


# ════════════════════════════════════════════════════════════
#  HTTP HELPERS
# ════════════════════════════════════════════════════════════

def http_get(url: str, retries: int = 4, stream: bool = False,
             headers: dict = None) -> requests.Response:
    """GET with retry/backoff on timeouts and 5xx/429 responses."""
    hdrs = dict(UA)
    if headers:
        hdrs.update(headers)
    last_exc = None
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, headers=hdrs, stream=stream, timeout=TIMEOUT)
            if r.status_code in (429, 500, 502, 503, 504):
                wait = min(30, 2 ** attempt)
                log(f"  HTTP {r.status_code} — retry {attempt}/{retries} in {wait}s", "warn")
                time.sleep(wait)
                continue
            return r
        except requests.RequestException as e:
            last_exc = e
            wait = min(30, 2 ** attempt)
            log(f"  request error ({e}) — retry {attempt}/{retries} in {wait}s", "warn")
            time.sleep(wait)
    if last_exc:
        raise last_exc
    raise RuntimeError(f"GET failed after {retries} attempts: {url}")


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(CHUNK_SIZE), b""):
            h.update(chunk)
    return h.hexdigest()


def download_file(url: str, dest: Path, retries: int = 4,
                  resume: bool = True, expected_md5: str = "") -> dict:
    """
    Stream a file to disk with resume (HTTP Range) support.
    Returns {status, size_bytes, path, error}.
    status ∈ {ok, already_complete, not_found, md5_mismatch, error}
    """
    dest.parent.mkdir(parents=True, exist_ok=True)

    # If a complete, md5-matching file already exists, skip re-download.
    if dest.exists() and expected_md5:
        if md5_of(dest) == expected_md5:
            return {"status": "already_complete", "size_bytes": dest.stat().st_size,
                    "path": str(dest), "error": ""}

    existing = dest.stat().st_size if dest.exists() else 0
    headers  = {}
    if resume and existing:
        headers["Range"] = f"bytes={existing}-"

    try:
        r = http_get(url, retries=retries, stream=True, headers=headers)
        if r.status_code == 416:   # range beyond EOF → already complete
            r.close()
            return {"status": "already_complete", "size_bytes": existing,
                    "path": str(dest), "error": ""}
        if r.status_code == 404:
            r.close()
            return {"status": "not_found", "size_bytes": 0,
                    "path": str(dest), "error": "404"}
        r.raise_for_status()

        mode = "ab" if (resume and existing and r.status_code == 206) else "wb"
        with open(dest, mode) as f:
            for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                if chunk:
                    f.write(chunk)

        size = dest.stat().st_size
        if expected_md5:
            got = md5_of(dest)
            if got != expected_md5:
                return {"status": "md5_mismatch", "size_bytes": size,
                        "path": str(dest), "error": f"expected {expected_md5}, got {got}"}
        return {"status": "ok", "size_bytes": size, "path": str(dest), "error": ""}

    except Exception as e:
        return {"status": "error", "size_bytes": existing,
                "path": str(dest), "error": str(e)}


# ════════════════════════════════════════════════════════════
#  GEO FTP LAYOUT
# ════════════════════════════════════════════════════════════

def acc_stub(acc: str) -> str:
    """GSE200637 → GSE200nnn ; GSM1234567 → GSM1234nnn."""
    return acc[:-3] + "nnn"


def series_dir(gse: str) -> str:
    return f"{GEO_FTP_BASE}/series/{acc_stub(gse)}/{gse}"


def sample_suppl_dir(gsm: str) -> str:
    return f"{GEO_FTP_BASE}/samples/{acc_stub(gsm)}/{gsm}/suppl/"


def list_ftp_dir(url: str, retries: int = 4) -> list[str]:
    """
    Parse an NCBI FTP HTTPS directory index and return real file names
    (navigation links, parent dirs and query-string sort links excluded).
    """
    r = http_get(url, retries=retries)
    if r.status_code == 404:
        return []
    r.raise_for_status()
    names = []
    for m in re.finditer(r'href="([^"]+)"', r.text):
        name = m.group(1)
        # Skip absolute URLs, parent/sort links, and subdirectories.
        if name.startswith(("http", "/", "?")) or "/" in name.strip("/"):
            continue
        name = name.strip("/")
        if not name or name in ("..",):
            continue
        names.append(name)
    return sorted(set(names))


# ════════════════════════════════════════════════════════════
#  SERIES MATRIX
# ════════════════════════════════════════════════════════════

def download_series_matrix(gse: str, outdir: Path, retries: int = 4,
                           list_only: bool = False) -> list[Path]:
    """
    Download every series matrix file for the GSE (there may be more than one
    when a series spans multiple platforms). Returns local paths.
    """
    matrix_url_dir = f"{series_dir(gse)}/matrix/"
    files = [f for f in list_ftp_dir(matrix_url_dir, retries)
             if f.endswith("series_matrix.txt.gz")]
    if not files:
        log(f"  {gse}: no series matrix found on GEO FTP", "warn")
        return []

    dest_dir = outdir / gse / "matrix"
    paths = []
    for fname in files:
        dest = dest_dir / fname
        if list_only:
            log(f"  [list] series matrix: {fname}")
            paths.append(dest)
            continue
        res = download_file(matrix_url_dir + fname, dest, retries=retries)
        if res["status"] in ("ok", "already_complete"):
            log(f"  {gse}: series matrix {fname} ({human(res['size_bytes'])})", "ok")
            paths.append(dest)
        else:
            log(f"  {gse}: series matrix {fname} failed — {res['error']}", "err")
    return paths


def _parse_matrix_lines(lines) -> dict:
    """Shared parser: extract samples, per-sample suppl URLs, study accessions."""
    meta = {"samples": [], "sample_suppl_urls": [], "study_accessions": [],
            "title": ""}
    study = set()
    suppl = set()

    for line in lines:
        line = line.rstrip("\n")

        if line.startswith("!Series_title"):
            meta["title"] = line.partition(" = ")[2].strip().strip('"')

        elif line.startswith("!Series_relation"):
            val = line.partition(" = ")[2]
            study.update(re.findall(r"(?:SRP|ERP|DRP)\d+", val))
            study.update(re.findall(r"(?:PRJNA|PRJEB|PRJDB)\d+", val))

        elif line.startswith("!Sample_geo_accession"):
            meta["samples"] = [v.strip().strip('"')
                               for v in line.split("\t")[1:] if v.strip().strip('"')]

        elif line.startswith("!Sample_supplementary_file"):
            for v in line.split("\t")[1:]:
                u = v.strip().strip('"')
                if u and u.upper() != "NONE":
                    suppl.add(u)

    meta["study_accessions"] = sorted(study)
    meta["sample_suppl_urls"] = sorted(suppl)
    return meta


def parse_series_matrix(path: Path) -> dict:
    """Parse a series matrix .txt.gz from disk."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        return _parse_matrix_lines(f)


def load_matrix_meta_in_memory(gse: str, retries: int = 4) -> dict:
    """
    Fetch the first series matrix into memory and parse it — used by --list so
    we can report per-sample files and raw runs without writing to disk.
    """
    matrix_url_dir = f"{series_dir(gse)}/matrix/"
    files = [f for f in list_ftp_dir(matrix_url_dir, retries)
             if f.endswith("series_matrix.txt.gz")]
    if not files:
        return {"samples": [], "sample_suppl_urls": [], "study_accessions": [], "title": ""}
    r = http_get(matrix_url_dir + files[0], retries=retries)
    r.raise_for_status()
    try:
        text = gzip.decompress(r.content).decode("utf-8", errors="replace")
    except OSError:
        text = r.content.decode("utf-8", errors="replace")
    return _parse_matrix_lines(text.splitlines())


# ════════════════════════════════════════════════════════════
#  SUPPLEMENTARY FILES
# ════════════════════════════════════════════════════════════

def download_series_suppl(gse: str, outdir: Path, retries: int = 4,
                          list_only: bool = False,
                          skip_raw_tar: bool = False) -> list[dict]:
    """Download series-level supplementary files from GEO FTP suppl/ dir."""
    suppl_url = f"{series_dir(gse)}/suppl/"
    files = [f for f in list_ftp_dir(suppl_url, retries)
             if not f.endswith("series_matrix.txt.gz")]
    if skip_raw_tar:
        skipped = [f for f in files if f.upper().endswith("_RAW.TAR")]
        for f in skipped:
            log(f"  {gse}: skipping {f} (--skip-raw-tar)", "warn")
        files = [f for f in files if not f.upper().endswith("_RAW.TAR")]
    if not files:
        log(f"  {gse}: no series-level supplementary files", "warn")
        return []

    dest_dir = outdir / gse / "suppl"
    results = []
    for fname in files:
        dest = dest_dir / fname
        if list_only:
            log(f"  [list] series suppl: {fname}")
            results.append({"name": fname, "status": "listed"})
            continue
        log(f"  {gse}: downloading {fname} ...")
        res = download_file(suppl_url + fname, dest, retries=retries)
        log(f"    {fname}: {human(res['size_bytes'])} [{res['status']}]",
            "ok" if res["status"] in ("ok", "already_complete") else "warn")
        results.append({"name": fname, **res})
    return results


def _url_to_https(u: str) -> str:
    """Convert an ftp:// GEO/NCBI URL to https:// (NCBI serves both)."""
    if u.startswith("ftp://"):
        return "https://" + u[len("ftp://"):]
    return u


def download_sample_suppl(meta: dict, gse: str, outdir: Path,
                          retries: int = 4, list_only: bool = False) -> list[dict]:
    """
    Download per-sample (GSM) supplementary files. Uses the URLs listed in the
    series matrix; SRA experiment links are skipped (those are raw reads,
    handled separately by --raw).
    """
    urls = [u for u in meta.get("sample_suppl_urls", [])
            if "/sra/" not in u.lower() and "trace.ncbi" not in u.lower()]
    if not urls:
        log(f"  {gse}: no per-sample supplementary file URLs in matrix", "warn")
        return []

    dest_dir = outdir / gse / "samples"
    results = []
    for u in urls:
        u_https = _url_to_https(u)
        # Group each file under its GSM accession when discoverable in the URL.
        gsm_m = re.search(r"(GSM\d+)", u_https)
        sub = gsm_m.group(1) if gsm_m else ""
        fname = u_https.rstrip("/").split("/")[-1]
        dest = dest_dir / sub / fname if sub else dest_dir / fname
        if list_only:
            log(f"  [list] sample suppl: {sub + '/' if sub else ''}{fname}")
            results.append({"name": fname, "gsm": sub, "status": "listed"})
            continue
        res = download_file(u_https, dest, retries=retries)
        log(f"  {gse}: {sub + '/' if sub else ''}{fname} "
            f"{human(res['size_bytes'])} [{res['status']}]",
            "ok" if res["status"] in ("ok", "already_complete") else "warn")
        results.append({"name": fname, "gsm": sub, **res})
    return results


# ════════════════════════════════════════════════════════════
#  RAW READS  (GSE → SRA study → ENA FASTQ)
# ════════════════════════════════════════════════════════════

def ena_runs(accession: str, retries: int = 4) -> list[dict]:
    """
    Query the ENA portal for all sequencing runs under a study/run/sample
    accession. Returns run dicts with fastq URLs, byte sizes and md5 sums.
    """
    fields = ("run_accession,sample_accession,experiment_accession,"
              "fastq_ftp,fastq_bytes,fastq_md5,library_strategy,read_count")
    url = (f"{ENA_PORTAL}?accession={accession}&result=read_run"
           f"&fields={fields}&format=tsv&limit=0")
    r = http_get(url, retries=retries)
    if r.status_code != 200 or not r.text.strip():
        return []
    return parse_ena_tsv(r.text)


def parse_ena_tsv(text: str) -> list[dict]:
    """
    Parse an ENA filereport TSV into per-file run dicts. A run with paired-end
    reads lists two ';'-separated fastq_ftp entries, expanded into two dicts.
    Kept separate from the HTTP call so it can be unit-tested offline.
    """
    lines = text.rstrip("\n").split("\n")
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    runs = []
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t")))
        ftp = row.get("fastq_ftp", "")
        if not ftp:
            continue
        files = ftp.split(";")
        sizes = (row.get("fastq_bytes", "") or "").split(";")
        md5s  = (row.get("fastq_md5", "") or "").split(";")
        for i, fpath in enumerate(files):
            if not fpath:
                continue
            runs.append({
                "run":      row.get("run_accession", ""),
                "sample":   row.get("sample_accession", ""),
                "strategy": row.get("library_strategy", ""),
                "url":      "https://" + fpath,
                "bytes":    int(sizes[i]) if i < len(sizes) and sizes[i].isdigit() else 0,
                "md5":      md5s[i] if i < len(md5s) else "",
            })
    return runs


def ncbi_gsm_to_srr(gsm: str, retries: int = 4) -> list[str]:
    """Fallback: resolve a GSM to SRR run accessions via NCBI eutils."""
    time.sleep(NCBI_DELAY)
    try:
        r = http_get(f"{NCBI_EUTILS}/esearch.fcgi?db=gds&term={gsm}[Accession]"
                     f"&retmode=json&retmax=1", retries=retries)
        uids = r.json()["esearchresult"].get("idlist", [])
        if not uids:
            return []
        time.sleep(NCBI_DELAY)
        r = http_get(f"{NCBI_EUTILS}/elink.fcgi?dbfrom=gds&db=sra&id={uids[0]}"
                     f"&retmode=json", retries=retries)
        sra_uids = []
        for ls in r.json().get("linksets", []):
            for lsd in ls.get("linksetdbs", []):
                if lsd.get("dbto") == "sra":
                    sra_uids.extend(lsd.get("links", []))
        if not sra_uids:
            return []
        time.sleep(NCBI_DELAY)
        ids = ",".join(str(u) for u in sra_uids[:50])
        r = http_get(f"{NCBI_EUTILS}/esummary.fcgi?db=sra&id={ids}&retmode=json",
                     retries=retries)
        result = r.json().get("result", {})
        srrs = set()
        for uid, rec in result.items():
            if uid == "uids":
                continue
            srrs.update(re.findall(r'acc="(SRR\d+|ERR\d+|DRR\d+)"',
                                   rec.get("runs", "")))
        return sorted(srrs)
    except Exception as e:
        log(f"  GSM→SRR resolution failed for {gsm}: {e}", "warn")
        return []


def resolve_raw_runs(meta: dict, gse: str, retries: int = 4) -> list[dict]:
    """
    Build the list of FASTQ files to download for a series.
    Primary path: study accession from the series matrix → ENA.
    Fallback:     resolve each GSM → SRR via NCBI, then ENA per run.
    """
    runs: list[dict] = []
    seen = set()

    studies = meta.get("study_accessions", [])
    if studies:
        for acc in studies:
            log(f"  {gse}: querying ENA for study {acc} ...")
            for run in ena_runs(acc, retries):
                if run["url"] not in seen:
                    seen.add(run["url"])
                    runs.append(run)
    if not runs and meta.get("samples"):
        log(f"  {gse}: no study accession in matrix — resolving GSMs via NCBI "
            f"(this is slower) ...", "warn")
        srrs = set()
        for gsm in meta["samples"]:
            srrs.update(ncbi_gsm_to_srr(gsm, retries))
        for srr in sorted(srrs):
            for run in ena_runs(srr, retries):
                if run["url"] not in seen:
                    seen.add(run["url"])
                    runs.append(run)
    return runs


def download_raw_reads(runs: list[dict], gse: str, outdir: Path,
                       retries: int = 4, verify_md5: bool = True) -> list[dict]:
    """Download each resolved FASTQ file into outdir/GSE/fastq/<run>/."""
    dest_root = outdir / gse / "fastq"
    results = []
    for i, run in enumerate(runs, 1):
        fname = run["url"].split("/")[-1]
        dest  = dest_root / (run["run"] or "unknown") / fname
        log(f"  [{i}/{len(runs)}] {run['run']} {fname} "
            f"({human(run['bytes']) if run['bytes'] else '?'}) ...")
        res = download_file(run["url"], dest, retries=retries,
                            expected_md5=run["md5"] if verify_md5 else "")
        status = res["status"]
        log(f"    {fname}: {human(res['size_bytes'])} [{status}]",
            "ok" if status in ("ok", "already_complete") else
            "warn" if status == "md5_mismatch" else "err")
        results.append({"run": run["run"], "name": fname, **res})
    return results


# ════════════════════════════════════════════════════════════
#  PER-GSE ORCHESTRATION
# ════════════════════════════════════════════════════════════

def fetch_gse(gse: str, outdir: Path, want: dict, args) -> dict:
    print()
    log(f"══ {gse} ══", "info")
    summary = {"accession": gse, "title": "", "matrix": [], "series_suppl": [],
               "sample_suppl": [], "raw": []}

    # Matrix is downloaded whenever it's wanted OR needed to find per-sample
    # supplementary URLs / SRA study accessions.
    need_meta   = want["sample_suppl"] or want["raw"]
    need_matrix = want["matrix"] or need_meta
    meta = {"samples": [], "sample_suppl_urls": [], "study_accessions": [], "title": ""}
    matrix_paths = []

    if need_matrix:
        matrix_paths = download_series_matrix(gse, outdir, args.retries, args.list)

    if args.list:
        # Dry run: parse the matrix in memory (no disk write) so we can still
        # report per-sample files and raw runs.
        if need_meta:
            try:
                meta = load_matrix_meta_in_memory(gse, args.retries)
            except Exception as e:
                log(f"  {gse}: could not parse matrix for listing — {e}", "warn")
    elif matrix_paths:
        # Parse the first matrix file (samples/relations are consistent across
        # per-platform matrices for a series).
        meta = parse_series_matrix(matrix_paths[0])
        log(f"  {gse}: {len(meta['samples'])} samples | "
            f"study={','.join(meta['study_accessions']) or 'none'}", "ok")

    summary["title"] = meta.get("title", "")
    if want["matrix"]:
        summary["matrix"] = [str(p) for p in matrix_paths]

    if want["series_suppl"]:
        summary["series_suppl"] = download_series_suppl(
            gse, outdir, args.retries, args.list, skip_raw_tar=args.skip_raw_tar)

    if want["sample_suppl"]:
        summary["sample_suppl"] = download_sample_suppl(
            meta, gse, outdir, args.retries, args.list)

    if want["raw"]:
        runs = resolve_raw_runs(meta, gse, args.retries)
        if args.max_runs and len(runs) > args.max_runs:
            log(f"  {gse}: capping raw runs {len(runs)} → {args.max_runs} "
                f"(--max-runs)", "warn")
            runs = runs[:args.max_runs]
        total = sum(r["bytes"] for r in runs)
        log(f"  {gse}: {len(runs)} FASTQ file(s), est. {human(total)}",
            "ok" if runs else "warn")

        if args.list:
            for r in runs:
                log(f"  [list] fastq: {r['run']}  {r['url'].split('/')[-1]}  "
                    f"({human(r['bytes']) if r['bytes'] else '?'})")
        elif runs:
            if not args.yes:
                ans = input(f"    Download {len(runs)} FASTQ files "
                            f"(~{human(total)}) for {gse}? [y/N] ").strip().lower()
                if ans != "y":
                    log(f"  {gse}: raw download skipped by user", "warn")
                    return summary
            summary["raw"] = download_raw_reads(
                runs, gse, outdir, args.retries, verify_md5=not args.no_md5)

    return summary


# ════════════════════════════════════════════════════════════
#  CLI
# ════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(
        description="Download data from GEO for one or more GSE IDs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Usage")[1] if "Usage" in __doc__ else "")
    p.add_argument("gse", nargs="+", help="One or more GSE accessions (e.g. GSE200637)")
    p.add_argument("-o", "--outdir", default="geo_downloads",
                   help="Output directory (default: ./geo_downloads)")
    p.add_argument("--matrix",       action="store_true", help="Download series matrix (metadata)")
    p.add_argument("--series-suppl", action="store_true", help="Download series-level supplementary files")
    p.add_argument("--sample-suppl", action="store_true", help="Download per-sample (GSM) supplementary files")
    p.add_argument("--raw",          action="store_true", help="Download raw FASTQ reads via ENA")
    p.add_argument("--all",          action="store_true", help="Fetch matrix + series-suppl + sample-suppl + raw")
    p.add_argument("--list",         action="store_true", help="Dry run: list what would be downloaded, then exit")
    p.add_argument("--skip-raw-tar", action="store_true",
                   help="Skip the series-level GSExxx_RAW.tar bundle (redundant with per-sample files)")
    p.add_argument("--max-runs",     type=int, default=0, help="Cap number of raw FASTQ runs (0 = no cap)")
    p.add_argument("--no-md5",       action="store_true", help="Skip md5 verification of FASTQ files")
    p.add_argument("--yes",          action="store_true", help="Don't prompt before large raw downloads")
    p.add_argument("--retries",      type=int, default=4, help="HTTP retry attempts per request (default: 4)")
    args = p.parse_args()

    # Resolve which content types are wanted.
    if args.all:
        want = {"matrix": True, "series_suppl": True, "sample_suppl": True, "raw": True}
    elif any([args.matrix, args.series_suppl, args.sample_suppl, args.raw]):
        want = {"matrix": args.matrix, "series_suppl": args.series_suppl,
                "sample_suppl": args.sample_suppl, "raw": args.raw}
    else:
        # Default: everything except raw reads (which are heavy).
        want = {"matrix": True, "series_suppl": True, "sample_suppl": True, "raw": False}

    # Validate accessions.
    gses = []
    for g in args.gse:
        g = g.strip().upper()
        if re.fullmatch(r"GSE\d+", g):
            gses.append(g)
        else:
            log(f"Ignoring invalid GSE accession: {g}", "warn")
    if not gses:
        log("No valid GSE accessions given.", "err")
        sys.exit(1)

    outdir = Path(args.outdir)
    log(f"Output dir: {outdir.resolve()}")
    log(f"Fetching: " + ", ".join(k for k, v in want.items() if v)
        + ("   [LIST ONLY]" if args.list else ""))

    all_summaries = []
    try:
        for gse in gses:
            all_summaries.append(fetch_gse(gse, outdir, want, args))
    except KeyboardInterrupt:
        print()
        log("Interrupted. Re-run to resume — completed files are skipped.", "warn")
        sys.exit(130)

    # Final recap.
    print()
    for s in all_summaries:
        n_matrix = len(s["matrix"])
        n_ss     = len([x for x in s["series_suppl"] if x.get("status") in ("ok", "already_complete", "listed")])
        n_gsm    = len([x for x in s["sample_suppl"] if x.get("status") in ("ok", "already_complete", "listed")])
        n_raw    = len([x for x in s["raw"] if x.get("status") in ("ok", "already_complete")])
        log(f"{s['accession']}: matrix={n_matrix} series-suppl={n_ss} "
            f"sample-suppl={n_gsm} fastq={n_raw}", "ok")
    log("Done.", "ok")


if __name__ == "__main__":
    main()
