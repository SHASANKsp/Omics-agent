#!/usr/bin/env python3
"""
Omics Dataset Research Agent — Ollama Edition
───────────────────────────────────────────────
AI-powered dataset curation for disease omics research.
Runs fully locally via Ollama. No API key needed.

Recommended models (install via ollama pull):
  ollama pull mistral          # fast, good JSON
  ollama pull llama3.1         # best quality
  ollama pull llama3.2         # lighter
  ollama pull gemma3           # good alternative
  ollama pull qwen2.5          # strong on structured output

Usage:
    python omics_agent.py
    python omics_agent.py --disease "Systemic lupus erythematosus"
    python omics_agent.py --disease "ALS" --model llama3.1 --out results/
    python omics_agent.py --list-models

Requirements:
    pip install requests rich
    # Ollama must be running: https://ollama.com
"""

import argparse
import json
import os
import re
import sys
import time
import csv
from datetime import datetime
from pathlib import Path

import requests

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


OLLAMA_BASE = "http://localhost:11434"

# ─── PROMPTS ────────────────────────────────────────────────────────────────
# Extra explicit for local models — they need more nudging to stay on JSON

DISEASE_ANALYSIS_PROMPT = """You are an expert computational biologist and omics data scientist.

A researcher wants to study the disease: "{disease}"

Your task: analyze the disease and output a JSON object describing the best omics research strategy.

RULES:
- Output ONLY raw JSON. No explanation before or after. No markdown. No triple backticks.
- Start your response with {{ and end with }}

JSON structure to fill in:
{{
  "category": "<one of: oncology, neurodegenerative, metabolic, autoimmune, infectious, cardiovascular, rare_genetic, other>",
  "primaryTissue": "<most relevant tissue or cell type>",
  "mechanisticFocus": "<one sentence describing the core mechanism to study>",
  "researchGoals": ["<goal 1>", "<goal 2>"],
  "omicsTypes": [
    {{
      "name": "<omics name, e.g. RNA-seq>",
      "priority": "<primary or secondary>",
      "reason": "<one sentence: why this omics is needed for this disease>",
      "database": "<database name, e.g. NCBI SRA / GEO>",
      "dbCode": "<short code, e.g. SRA, GEO, TCGA, PRIDE, MetaboLights>",
      "searchTerms": ["<ncbi search term 1>", "<term 2>", "<term 3>"]
    }}
  ]
}}

Omics selection guide:
- RNA-seq: ALWAYS include as primary (bulk transcriptomics)
- scRNA-seq: add for cancer, neurodegeneration, autoimmune — any heterogeneous disease
- ATAC-seq: add for cancer or epigenetic/chromatin diseases
- ChIP-seq: add when transcription factor dysregulation is central
- WGS or WES: add for genetic diseases or cancer somatic mutations
- Proteomics (PRIDE): add when protein-level validation matters
- DNA methylation arrays (GEO): add for epigenetic diseases
- Metabolomics (MetaboLights): add for metabolic diseases
- Select 3 to 5 omics types total, most relevant ones only

Now output the JSON for: "{disease}"
"""

DATASET_CURATION_PROMPT = """You are a senior bioinformatics database curator with expertise in omics data repositories.

Disease to research: "{disease}"
Omics types needed: {omics_names}
Primary tissue/cells: {tissue}

Your task: generate a curated list of realistic research datasets.

RULES:
- Output ONLY a raw JSON array. No explanation. No markdown. No backticks.
- Start your response with [ and end with ]
- Generate exactly 10 to 14 dataset entries

Each dataset entry must follow this exact structure:
{{
  "id": "<realistic accession ID>",
  "title": "<descriptive dataset title>",
  "database": "<GEO, SRA, TCGA, ENCODE, PRIDE, or MetaboLights>",
  "omicsType": "<omics type from the list above>",
  "organism": "<Homo sapiens or Mus musculus>",
  "sampleCount": <integer number of samples>,
  "condition": "<e.g. SLE patients vs healthy controls>",
  "tissue": "<tissue or cell type>",
  "platform": "<sequencing or array platform>",
  "year": <year as integer, between 2018 and 2024>,
  "pmid": "<PubMed ID as string>",
  "tags": ["<tag1>", "<tag2>"],
  "downloadUrl": "<full URL to access the dataset>",
  "relevanceScore": <integer 1 to 100>,
  "notes": "<one sentence: why this dataset is useful for {disease} research>"
}}

Accession ID format rules:
- GEO: GSE followed by 6 digits, e.g. GSE198432
- SRA: SRP followed by 6 digits, e.g. SRP298765
- TCGA: TCGA-XX format, e.g. TCGA-LUAD
- PRIDE: PXD followed by 6 digits, e.g. PXD034521
- MetaboLights: MTBLS followed by 4 digits, e.g. MTBLS1234

Tag options: paired-end, single-end, bulk, single-cell, FFPE, fresh-frozen, longitudinal, treatment-naive, multi-omics, public, time-series

Spread datasets across all omics types: {omics_names}
Vary sample sizes between 8 and 500.
relevanceScore should reflect how directly this dataset addresses {disease}.

Now output the JSON array of datasets for: "{disease}"
"""


# ─── OLLAMA CLIENT ──────────────────────────────────────────────────────────

class OllamaClient:
    def __init__(self, model: str, base_url: str = OLLAMA_BASE):
        self.model = model
        self.base_url = base_url.rstrip("/")

    def check_connection(self) -> bool:
        try:
            r = requests.get(f"{self.base_url}/api/tags", timeout=5)
            return r.status_code == 200
        except requests.ConnectionError:
            return False

    def list_models(self) -> list[str]:
        try:
            r = requests.get(f"{self.base_url}/api/tags", timeout=5)
            r.raise_for_status()
            return [m["name"] for m in r.json().get("models", [])]
        except Exception:
            return []

    def generate(self, prompt: str, temperature: float = 0.1) -> str:
        """Call Ollama generate endpoint (non-streaming)."""
        payload = {
            "model": self.model,
            "prompt": prompt,
            "stream": False,
            "options": {
                "temperature": temperature,
                "top_p": 0.9,
                "num_predict": 4096,
            }
        }
        r = requests.post(
            f"{self.base_url}/api/generate",
            json=payload,
            timeout=300
        )
        r.raise_for_status()
        return r.json()["response"].strip()


# ─── JSON EXTRACTION ─────────────────────────────────────────────────────────
# Local models sometimes wrap output in markdown or add preamble text.
# These helpers robustly extract JSON from whatever the model returns.

def extract_json_object(text: str) -> dict:
    """Extract the first complete JSON object from model output."""
    # Try direct parse first
    text = text.strip()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

    # Strip markdown fences
    fenced = re.sub(r"```(?:json)?(.*?)```", r"\1", text, flags=re.DOTALL)
    try:
        return json.loads(fenced.strip())
    except json.JSONDecodeError:
        pass

    # Find first { ... } block
    start = text.find("{")
    if start == -1:
        raise ValueError("No JSON object found in model response")

    depth = 0
    for i, ch in enumerate(text[start:], start):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                try:
                    return json.loads(text[start:i+1])
                except json.JSONDecodeError:
                    break

    raise ValueError(f"Could not parse JSON object from response:\n{text[:500]}")


def extract_json_array(text: str) -> list:
    """Extract the first complete JSON array from model output."""
    text = text.strip()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

    # Strip markdown fences
    fenced = re.sub(r"```(?:json)?(.*?)```", r"\1", text, flags=re.DOTALL)
    try:
        return json.loads(fenced.strip())
    except json.JSONDecodeError:
        pass

    # Find first [ ... ] block
    start = text.find("[")
    if start == -1:
        raise ValueError("No JSON array found in model response")

    depth = 0
    for i, ch in enumerate(text[start:], start):
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
            if depth == 0:
                try:
                    return json.loads(text[start:i+1])
                except json.JSONDecodeError:
                    break

    raise ValueError(f"Could not parse JSON array from response:\n{text[:500]}")


# ─── AGENT ──────────────────────────────────────────────────────────────────

class OmicsAgent:
    def __init__(self, model: str):
        self.client = OllamaClient(model)
        self.model = model

    def analyze_disease(self, disease: str) -> dict:
        prompt = DISEASE_ANALYSIS_PROMPT.format(disease=disease)
        raw = self.client.generate(prompt, temperature=0.1)
        return extract_json_object(raw)

    def curate_datasets(self, disease: str, analysis: dict) -> list:
        omics_names = ", ".join(o["name"] for o in analysis["omicsTypes"])
        prompt = DATASET_CURATION_PROMPT.format(
            disease=disease,
            omics_names=omics_names,
            tissue=analysis.get("primaryTissue", "relevant tissue")
        )
        raw = self.client.generate(prompt, temperature=0.2)
        return extract_json_array(raw)


# ─── DISPLAY ────────────────────────────────────────────────────────────────

def print_step(msg: str, status: str = "info"):
    ts = datetime.now().strftime("%H:%M:%S")
    symbols = {"info": "·", "success": "✓", "warn": "!", "error": "✗"}
    colors  = {"info": "dim", "success": "green", "warn": "yellow", "error": "red"}
    sym = symbols.get(status, "·")
    if HAS_RICH:
        col = colors.get(status, "white")
        console.print(f"[dim]{ts}[/dim]  [{col}]{sym}[/{col}]  {escape(str(msg))}")
    else:
        print(f"{ts}  {sym}  {msg}")


def print_omics_strategy(analysis: dict):
    if HAS_RICH:
        table = Table(title="Omics Strategy", show_header=True,
                      header_style="bold", box=None, padding=(0, 2))
        table.add_column("Omics Type", style="bold cyan", no_wrap=True)
        table.add_column("Priority", no_wrap=True)
        table.add_column("Database")
        table.add_column("Reason", max_width=52)
        for o in analysis.get("omicsTypes", []):
            p = o.get("priority", "secondary")
            table.add_row(
                o.get("name", ""),
                f"[{'green' if p == 'primary' else 'blue'}]{p}[/{'green' if p == 'primary' else 'blue'}]",
                o.get("database", ""),
                o.get("reason", "")
            )
        console.print()
        console.print(table)
        console.print()
    else:
        print("\n=== Omics Strategy ===")
        for o in analysis.get("omicsTypes", []):
            print(f"  [{o.get('priority','?')}] {o.get('name','')} → {o.get('database','')}")
            print(f"    {o.get('reason','')}")
        print()


def print_datasets(datasets: list):
    db_icons = {
        "GEO": "🔵", "SRA": "🟢", "TCGA": "🟠",
        "ENCODE": "🟣", "PRIDE": "🟡", "MetaboLights": "🟤"
    }
    by_db: dict = {}
    for d in datasets:
        by_db.setdefault(d.get("database", "Other"), []).append(d)

    if not HAS_RICH:
        print("\n=== Curated Datasets ===")
        for db, items in by_db.items():
            print(f"\n── {db} ──")
            for d in sorted(items, key=lambda x: -x.get("relevanceScore", 0)):
                print(f"  {d.get('id','')} | {d.get('omicsType','')} | {d.get('sampleCount','')} samples | {d.get('year','')}")
                print(f"  {d.get('title','')}")
                print(f"  Score: {d.get('relevanceScore','')}% | {d.get('downloadUrl','')}")
        return

    for db, items in by_db.items():
        icon = db_icons.get(db, "⚪")
        items_sorted = sorted(items, key=lambda x: -x.get("relevanceScore", 0))
        table = Table(
            title=f"{icon}  {db}  ({len(items)} datasets)",
            show_header=True, header_style="bold", box=None, padding=(0, 1)
        )
        table.add_column("Accession", style="cyan", no_wrap=True)
        table.add_column("Omics", no_wrap=True)
        table.add_column("Title", max_width=40)
        table.add_column("Samples", justify="right")
        table.add_column("Year", justify="right")
        table.add_column("Match", justify="right")
        for d in items_sorted:
            score = d.get("relevanceScore", 0)
            sc = "green" if score >= 85 else "yellow" if score >= 70 else "dim"
            table.add_row(
                d.get("id", ""),
                d.get("omicsType", ""),
                d.get("title", ""),
                str(d.get("sampleCount", "")),
                str(d.get("year", "")),
                f"[{sc}]{score}%[/{sc}]"
            )
        console.print(table)
        console.print()


# ─── SAVE ───────────────────────────────────────────────────────────────────

def save_csv(datasets: list, disease: str, out_dir: Path) -> Path:
    safe = re.sub(r"[^\w\s-]", "", disease).strip().replace(" ", "_")
    fname = out_dir / f"omics_datasets_{safe}.csv"
    fields = [
        "id", "title", "database", "omicsType", "organism",
        "sampleCount", "condition", "tissue", "platform",
        "year", "pmid", "tags", "downloadUrl", "relevanceScore", "notes"
    ]
    with open(fname, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for d in sorted(datasets, key=lambda x: -x.get("relevanceScore", 0)):
            row = {**d, "tags": "|".join(d.get("tags", []))}
            w.writerow(row)
    return fname


def save_json(report: dict, disease: str, out_dir: Path) -> Path:
    safe = re.sub(r"[^\w\s-]", "", disease).strip().replace(" ", "_")
    fname = out_dir / f"omics_report_{safe}.json"
    with open(fname, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    return fname


# ─── RUN ────────────────────────────────────────────────────────────────────

def run(disease: str, model: str, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    agent = OmicsAgent(model)

    # Check Ollama is reachable
    if not agent.client.check_connection():
        print_step(
            "Cannot reach Ollama at http://localhost:11434\n"
            "  → Make sure Ollama is running: ollama serve\n"
            "  → Or download it from https://ollama.com",
            "error"
        )
        sys.exit(1)

    # Check model is available
    available = agent.client.list_models()
    model_base = model.split(":")[0]
    if available and not any(model_base in m for m in available):
        print_step(f'Model "{model}" not found locally.', "warn")
        print_step(f"Available: {', '.join(available) or 'none'}", "warn")
        print_step(f"Pull it with:  ollama pull {model}", "warn")
        sys.exit(1)

    if HAS_RICH:
        console.print(Panel.fit(
            f"[bold]Omics Dataset Research Agent[/bold]  [dim](Ollama · {model})[/dim]\n"
            f"[dim]Disease:[/dim] {escape(disease)}",
            border_style="green"
        ))
    else:
        print(f"\n{'='*62}")
        print(f"  Omics Dataset Research Agent  (Ollama · {model})")
        print(f"  Disease: {disease}")
        print(f"{'='*62}\n")

    # ── Phase 1: Disease analysis ──
    print_step(f'Analyzing disease: "{disease}"')
    t0 = time.time()
    try:
        analysis = agent.analyze_disease(disease)
    except (ValueError, json.JSONDecodeError) as e:
        print_step(f"JSON parse error in disease analysis: {e}", "error")
        print_step("Try a different model: --model llama3.1 or --model mistral", "warn")
        sys.exit(1)

    print_step(
        f"Category: {analysis.get('category','?')} | "
        f"Tissue: {analysis.get('primaryTissue','?')}  ({time.time()-t0:.1f}s)",
        "success"
    )
    print_step(f"Focus: {analysis.get('mechanisticFocus','')}")
    print_step(
        f"Omics: {', '.join(o['name'] for o in analysis.get('omicsTypes', []))}",
        "success"
    )
    print_omics_strategy(analysis)

    # ── Phase 2: Dataset curation ──
    print_step("Curating datasets from SRA, GEO, TCGA, ENCODE, PRIDE, MetaboLights...")
    t0 = time.time()
    try:
        datasets = agent.curate_datasets(disease, analysis)
    except (ValueError, json.JSONDecodeError) as e:
        print_step(f"JSON parse error in dataset curation: {e}", "error")
        print_step("Try: --model llama3.1 or --model mistral", "warn")
        sys.exit(1)

    dbs = len(set(d.get("database", "") for d in datasets))
    total_samples = sum(d.get("sampleCount", 0) for d in datasets)
    print_step(
        f"Curated {len(datasets)} datasets across {dbs} databases, "
        f"{total_samples:,} total samples  ({time.time()-t0:.1f}s)",
        "success"
    )

    print_datasets(datasets)

    # ── Save ──
    report = {
        "disease": disease,
        "model": model,
        "timestamp": datetime.now().isoformat(),
        "analysis": analysis,
        "datasets": datasets
    }
    csv_path  = save_csv(datasets, disease, out_dir)
    json_path = save_json(report, disease, out_dir)

    if HAS_RICH:
        console.print(Panel(
            f"[green]✓[/green]  CSV:  [cyan]{csv_path}[/cyan]\n"
            f"[green]✓[/green]  JSON: [cyan]{json_path}[/cyan]",
            title="Saved", border_style="dim"
        ))
    else:
        print(f"\nSaved:\n  {csv_path}\n  {json_path}\n")

    return report


# ─── CLI ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Omics Dataset Research Agent — powered by local Ollama models",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python omics_agent.py
  python omics_agent.py --disease "Systemic lupus erythematosus"
  python omics_agent.py --disease "ALS" --model llama3.1
  python omics_agent.py --disease "Breast cancer" --out ~/research/datasets/
  python omics_agent.py --list-models

Recommended models (pull once with ollama pull <name>):
  mistral        fast, reliable JSON output
  llama3.1       best overall quality
  llama3.2       lighter, good for low-RAM systems
  gemma3         strong alternative
  qwen2.5        excellent structured output
        """
    )
    parser.add_argument("--disease",  "-d", type=str, help="Disease name to research")
    parser.add_argument("--model",    "-m", type=str, default="mistral",
                        help="Ollama model to use (default: mistral)")
    parser.add_argument("--out",      "-o", type=str, default="omics_results",
                        help="Output directory (default: omics_results)")
    parser.add_argument("--ollama-url", type=str, default=OLLAMA_BASE,
                        help=f"Ollama base URL (default: {OLLAMA_BASE})")
    parser.add_argument("--list-models", action="store_true",
                        help="List locally available Ollama models and exit")
    args = parser.parse_args()

    # Optionally override Ollama URL
    global OLLAMA_BASE
    OLLAMA_BASE = args.ollama_url

    if args.list_models:
        client = OllamaClient("", base_url=args.ollama_url)
        if not client.check_connection():
            print("Cannot reach Ollama. Is it running?  ollama serve")
            sys.exit(1)
        models = client.list_models()
        if models:
            print("Locally available models:")
            for m in models:
                print(f"  {m}")
        else:
            print("No models installed. Try:  ollama pull mistral")
        sys.exit(0)

    disease = args.disease
    if not disease:
        if HAS_RICH:
            disease = console.input("[bold green]Disease name:[/bold green] ").strip()
        else:
            disease = input("Disease name: ").strip()
    if not disease:
        print("No disease entered. Exiting.")
        sys.exit(1)

    try:
        run(disease, args.model, Path(args.out))
    except KeyboardInterrupt:
        print("\nInterrupted.")
    except requests.ConnectionError:
        print_step("Lost connection to Ollama. Is it still running?", "error")
        sys.exit(1)
    except Exception as e:
        print_step(f"Unexpected error: {e}", "error")
        raise


if __name__ == "__main__":
    main()
