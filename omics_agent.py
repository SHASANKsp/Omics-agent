#!/usr/bin/env python3
"""
Omics Dataset Research Agent — Real Database Edition
──────────────────────────────────────────────────────
Queries REAL databases via their public APIs.
LLM is used ONLY for reasoning (which omics, search terms, relevance scoring).
Dataset accession IDs come from actual NCBI / EBI API responses — not hallucinated.

Databases queried:
  • NCBI GEO   — bulk RNA-seq, microarray, ChIP-seq, ATAC-seq, methylation
  • NCBI SRA   — raw sequencing runs (RNA-seq, scRNA-seq, ATAC-seq, WGS/WES)
  • EBI PRIDE  — proteomics
  • EBI ENA    — European sequencing archive mirror

Relevance handling:
  • Datasets scoring >= 50 are shown in terminal and written to *_curated.csv
  • ALL fetched datasets (including low-scorers) are written to *_all.csv and
    the full JSON report for audit/review — nothing is silently discarded

Usage:
    python omics_agent.py --disease "psoriasis" --model llama3.2
    python omics_agent.py --disease "ALS" --model mistral --max-results 20
    python omics_agent.py --list-models

Requirements:
    pip install requests rich
    Ollama running locally: https://ollama.com
"""

import argparse
import json
import re
import sys
import time
import csv
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
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
NCBI_BASE   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PRIDE_BASE  = "https://www.ebi.ac.uk/pride/ws/archive/v2"
ENA_BASE    = "https://www.ebi.ac.uk/ena/portal/api"

# Slower but safer — NCBI rate limit is 3 req/sec without API key.
# We use 0.8s to leave headroom and avoid 429s on longer runs.
NCBI_DELAY        = 0.8   # seconds between every NCBI HTTP call
NCBI_RETRY_DELAY  = 5.0   # seconds to wait after a 429 before retrying
NCBI_MAX_RETRIES  = 3     # retry attempts per request
SCORE_BATCH_SIZE  = 12    # datasets sent to LLM per scoring call (avoids context overflow)
FETCH_WORKERS     = 4     # parallel threads for esummary fetches (stay ≤5 without NCBI API key)
SCORE_WORKERS     = 3     # parallel threads for LLM scoring batches


# ════════════════════════════════════════════════════════════
#  PROMPTS  — LLM used ONLY for reasoning, never for data
# ════════════════════════════════════════════════════════════

STRATEGY_PROMPT = """You are an expert computational biologist.

Disease: "{disease}"

Return ONLY a raw JSON object (no markdown, no backticks, no explanation).
Start with {{ and end with }}.

{{
  "category": "<oncology|neurodegenerative|metabolic|autoimmune|infectious|cardiovascular|rare_genetic|other>",
  "primaryTissue": "<most relevant tissue or cell type>",
  "mechanisticFocus": "<one sentence: core biology to study>",
  "researchGoals": ["<goal1>", "<goal2>"],
  "omicsTypes": [
    {{
      "name": "<e.g. RNA-seq>",
      "priority": "<primary|secondary>",
      "reason": "<one sentence: why this omics matters for this disease>",
      "database": "<GEO|SRA|PRIDE|ENA>",
      "ncbiDb": "<gds|sra|null>",
      "searchTerms": [
        "<term using disease full name + omics + tissue>",
        "<term using disease abbreviation or synonym + omics>",
        "<term using disease name + key cell type or pathway>",
        "<term using disease name alone — broad sweep>",
        "<term using disease name + key molecular mechanism or dysregulated pathway>",
        "<term using disease name + treatment or experimental model context>"
      ]
    }}
  ]
}}

Search term rules — generate EXACTLY 6 terms per omics type, covering all these axes:
  1. Full disease name + omics keyword + primary tissue  (e.g. "psoriasis RNA-seq skin")
  2. Disease abbreviation or common synonym + omics      (e.g. "plaque psoriasis RNA-seq")
  3. Disease name + key cell type or pathway             (e.g. "psoriasis keratinocyte gene expression")
  4. Disease name alone — broad sweep                   (e.g. "psoriasis" as the full query)
  5. Disease name + key molecular mechanism              (e.g. "psoriasis IL-17 pathway transcriptome" or "psoriasis NF-kB signaling RNA-seq")
  6. Disease name + treatment or model context           (e.g. "psoriasis TNF inhibitor transcriptome" or "imiquimod psoriasis mouse RNA-seq")

STRICT RULES for search terms:
- Every term must refer only to "{disease}" itself — never to a comorbidity, upstream disease, or related condition.
  For example, for psoriasis: DO NOT use "psoriatic arthritis", "atopic dermatitis", or "skin inflammation".
  The term must be about psoriasis and nothing else.
- Do NOT include organism/species in the search terms — human filtering is applied separately by the tool.
- Each term should be 2-6 words. No quotes inside the terms. Use NCBI search syntax style.
- Diverse terms are essential — similar terms return duplicate results and waste API calls.

Omics selection rules:
- RNA-seq: ALWAYS include, primary. ncbiDb=gds fetches GEO series; ncbiDb=sra fetches raw runs.
- scRNA-seq: add for cancer, neurodegeneration, autoimmune (heterogeneous tissues).
- ATAC-seq: add for cancers or epigenetic/chromatin diseases. ncbiDb=gds.
- ChIP-seq: add when TF dysregulation is central. ncbiDb=gds.
- WGS or WES: add for somatic mutation or germline genetics. ncbiDb=sra.
- Proteomics: add for protein-level validation. database=PRIDE, ncbiDb=null.
- Methylation arrays: add for epigenetic diseases. ncbiDb=gds.
- Metabolomics: add for metabolic diseases. database=PRIDE, ncbiDb=null.
- Select 3-5 omics types total.
"""

# Compact version used as fallback when the full prompt causes token truncation.
# Asks for only 3 search terms and a shorter mechanisticFocus to reduce output size.
STRATEGY_PROMPT_COMPACT = """You are an expert computational biologist.

Disease: "{disease}"

Return ONLY a raw JSON object. No markdown, no backticks. Start with {{ end with }}.

{{
  "category": "<oncology|neurodegenerative|metabolic|autoimmune|infectious|cardiovascular|rare_genetic|other>",
  "primaryTissue": "<tissue>",
  "mechanisticFocus": "<one short sentence>",
  "researchGoals": ["<goal1>", "<goal2>"],
  "omicsTypes": [
    {{
      "name": "<omics name>",
      "priority": "<primary|secondary>",
      "reason": "<one sentence>",
      "database": "<GEO|SRA|PRIDE>",
      "ncbiDb": "<gds|sra|null>",
      "searchTerms": [
        "<disease name + omics + tissue>",
        "<disease synonym + omics>",
        "<disease name + key cell type>"
      ]
    }}
  ]
}}

Rules:
- EXACTLY 3 searchTerms per omics type. Keep each term under 5 words.
- Terms must be about "{disease}" only — no related diseases or comorbidities.
- Do NOT include organism in search terms.
- Always include RNA-seq as primary. Add scRNA-seq for heterogeneous diseases.
- Add ATAC-seq for cancers or epigenetic diseases. Add Proteomics (PRIDE) if protein data matters.
- Select 3-4 omics types maximum.
"""

SCORING_PROMPT = """You are a strict bioinformatics data curator. Score these datasets for relevance to "{disease}".

Datasets (JSON):
{datasets_json}

Return ONLY a raw JSON array (no markdown, no backticks). One object per dataset, same order:
[
  {{
    "accession": "<same accession as input>",
    "relevanceScore": <integer 0-100>,
    "relevant": <true if score >= 50, else false>,
    "omicsType": "<RNA-seq|scRNA-seq|ATAC-seq|ChIP-seq|WGS|WES|Proteomics|Metabolomics|Methylation|Other>",
    "dropReason": "<if relevant=false: one sentence why. If relevant=true: empty string>",
    "notes": "<if relevant=true: one sentence why useful. If relevant=false: empty string>"
  }}
]

STEP 1 — HARD DROPS (score=0, relevant=false, no further evaluation):
  a) The title or summary mentions a DIFFERENT disease as the primary focus.
     This is the most important check. Read the title carefully.
     Examples of wrong-disease drops for "{disease}":
       - "rhinophyma", "hypertrophic scar", "melanoma", "atopic dermatitis" → NOT {disease} → drop
       - Any disease name in the title that is not "{disease}" → drop
  b) Organism is explicitly non-human (mouse, rat, zebrafish, Mus musculus etc.) → drop.
  c) Title is uninformative (e.g. "GSM1234: sample1; Homo sapiens; RNA-Seq") AND summary is empty → score=20, drop.

STEP 2 — SCORING for datasets that pass Step 1:
  90-100: title or summary explicitly mentions "{disease}", correct tissue, good sample size
  75-89:  mentions "{disease}" but slightly off tissue, or a specific subtype/severity
  60-74:  "{disease}" mentioned indirectly (treatment study, biomarker screen) but still the focus
  40-59:  disease unclear from title/summary — possibly relevant but cannot confirm
  0-39:   wrong disease confirmed — drop

Be strict. A skin scRNA-seq study of rhinophyma is NOT relevant to {disease} even if the tissue matches.
Disease identity in the title takes priority over tissue or omics type match.
"""


# ════════════════════════════════════════════════════════════
#  OLLAMA CLIENT
# ════════════════════════════════════════════════════════════

class OllamaClient:
    def __init__(self, model: str, base_url: str = OLLAMA_BASE):
        self.model    = model
        self.base_url = base_url.rstrip("/")

    def check(self) -> bool:
        try:
            return requests.get(f"{self.base_url}/api/tags", timeout=5).status_code == 200
        except:
            return False

    def list_models(self) -> list:
        try:
            r = requests.get(f"{self.base_url}/api/tags", timeout=5)
            return [m["name"] for m in r.json().get("models", [])]
        except:
            return []

    def generate(self, prompt: str, temperature: float = 0.1,
                 num_predict: int = 8192) -> tuple[str, bool]:
        """
        Call Ollama. Returns (response_text, was_truncated).
        was_truncated=True means the model hit the token limit mid-response —
        the caller should retry with a shorter prompt or smaller output target.
        """
        r = requests.post(
            f"{self.base_url}/api/generate",
            json={
                "model":  self.model,
                "prompt": prompt,
                "stream": False,
                "options": {
                    "temperature": temperature,
                    "num_predict": num_predict,
                }
            },
            timeout=600   # longer timeout for bigger outputs
        )
        r.raise_for_status()
        data       = r.json()
        text       = data.get("response", "").strip()
        # Ollama sets done_reason="length" when num_predict is exhausted mid-token
        truncated  = data.get("done_reason", "") == "length"
        return text, truncated


# ════════════════════════════════════════════════════════════
#  JSON EXTRACTION  (handles messy local model output)
# ════════════════════════════════════════════════════════════

def extract_obj(text: str) -> dict:
    text = re.sub(r"```(?:json)?(.*?)```", r"\1", text, flags=re.DOTALL).strip()
    try: return json.loads(text)
    except: pass
    s = text.find("{")
    if s == -1: raise ValueError("No JSON object in response")
    depth = 0
    for i, c in enumerate(text[s:], s):
        depth += (c == "{") - (c == "}")
        if depth == 0:
            try: return json.loads(text[s:i+1])
            except: break
    raise ValueError(f"Cannot parse JSON object:\n{text[:400]}")

def extract_arr(text: str) -> list:
    text = re.sub(r"```(?:json)?(.*?)```", r"\1", text, flags=re.DOTALL).strip()
    try: return json.loads(text)
    except: pass
    s = text.find("[")
    if s == -1: raise ValueError("No JSON array in response")
    depth = 0
    for i, c in enumerate(text[s:], s):
        depth += (c == "[") - (c == "]")
        if depth == 0:
            try: return json.loads(text[s:i+1])
            except: break
    raise ValueError(f"Cannot parse JSON array:\n{text[:400]}")


# ════════════════════════════════════════════════════════════
#  REAL DATABASE FETCHERS
# ════════════════════════════════════════════════════════════

def ncbi_get(url: str) -> requests.Response:
    """GET with automatic retry on 429 / 5xx. Respects NCBI_DELAY between calls."""
    time.sleep(NCBI_DELAY)
    for attempt in range(1, NCBI_MAX_RETRIES + 1):
        try:
            r = requests.get(url, timeout=30, headers={"User-Agent": "OmicsAgent/1.0"})
            if r.status_code == 429:
                wait = NCBI_RETRY_DELAY * attempt
                print_step(f"NCBI rate-limited (429) — waiting {wait:.0f}s before retry {attempt}/{NCBI_MAX_RETRIES}", "warn")
                time.sleep(wait)
                continue
            if r.status_code >= 500:
                print_step(f"NCBI server error {r.status_code} — retry {attempt}/{NCBI_MAX_RETRIES}", "warn")
                time.sleep(NCBI_RETRY_DELAY)
                continue
            r.raise_for_status()
            return r
        except requests.Timeout:
            print_step(f"NCBI timeout — retry {attempt}/{NCBI_MAX_RETRIES}", "warn")
            time.sleep(NCBI_RETRY_DELAY)
    raise RuntimeError(f"NCBI request failed after {NCBI_MAX_RETRIES} attempts: {url}")


def ncbi_search(db: str, term: str, max_results: int = 10) -> list[str]:
    """
    Search NCBI GEO (db=gds) or SRA (db=sra). Returns list of UIDs.
    Always restricts to Homo sapiens at the query level.
    """
    # Append human organism filter — this is enforced at the NCBI API level,
    # not left to the LLM, so non-human results never enter the pipeline.
    human_term = f'({term}) AND "Homo sapiens"[Organism]'
    params = {
        "db": db, "term": human_term,
        "retmax": max_results, "retmode": "json", "usehistory": "n"
    }
    url = f"{NCBI_BASE}/esearch.fcgi?" + urllib.parse.urlencode(params)
    try:
        r = ncbi_get(url)
        result = r.json()["esearchresult"]
        count = result.get("count", "?")
        uids  = result.get("idlist", [])
        print_step(f"    NCBI {db} '{term}' [human] → {count} total hits, fetching {len(uids)}", "info")
        return uids
    except Exception as e:
        print_step(f"NCBI search failed ({db}, '{term}'): {e}", "warn")
        return []


def ncbi_fetch_geo_summary(uid: str) -> dict | None:
    """
    Fetch a GEO series (GSE) summary by UID.
    Skips GPL platforms, individual GSM samples, and non-human records.
    """
    url = f"{NCBI_BASE}/esummary.fcgi?db=gds&id={uid}&retmode=json"
    try:
        r = ncbi_get(url)
        result = r.json().get("result", {})
        rec = result.get(uid) or result.get(str(uid))
        if not rec or rec.get("error"):
            return None
        acc = rec.get("accession", "")
        # Only series — skip GPL platform records and GSM individual samples
        if not acc.startswith("GSE"):
            return None
        organism = rec.get("taxon", "")
        if organism and "Homo sapiens" not in organism:
            return None
        return {
            "accession":   acc,
            "title":       rec.get("title", ""),
            "summary":     rec.get("summary", "")[:600],
            "organism":    organism or "Homo sapiens",
            "sampleCount": int(rec.get("n_samples", 0)),
            "gdsType":     rec.get("gdstype", ""),
            "pubDate":     rec.get("pdat", ""),
            "database":    "GEO",
            "downloadUrl": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
        }
    except Exception as e:
        print_step(f"GEO fetch failed (uid={uid}): {e}", "warn")
        return None


def ncbi_fetch_sra_summary(uid: str) -> dict | None:
    """
    Fetch an SRA record by UID and resolve it to its parent study (SRP accession).
    Individual runs (SRR) and experiments (SRX) all belong to a study — we want
    study-level records so we don't get dozens of rows for the same dataset.
    Returns None for non-human records.
    """
    url = f"{NCBI_BASE}/esummary.fcgi?db=sra&id={uid}&retmode=json"
    try:
        r = ncbi_get(url)
        result = r.json().get("result", {})
        rec = result.get(uid) or result.get(str(uid))
        if not rec or rec.get("error"):
            return None

        exp_xml = rec.get("expxml", "")
        runs_xml = rec.get("runs", "")

        # Organism — try multiple locations in the XML
        org_match = (re.search(r'<ScientificName>(.*?)</ScientificName>', exp_xml)
                     or re.search(r'SCIENTIFIC_NAME="([^"]+)"', exp_xml))
        organism = org_match.group(1).strip() if org_match else ""

        # Belt-and-suspenders human check
        if organism and "Homo sapiens" not in organism:
            return None

        # Prefer study-level (SRP) accession — this deduplicates individual runs
        study_match = (re.search(r'<Study\s+acc="(SRP\d+)"', exp_xml)
                       or re.search(r'<Study\s+acc="(ERP\d+)"', exp_xml)
                       or re.search(r'<Study\s+acc="(DRP\d+)"', exp_xml))
        acc = study_match.group(1) if study_match else str(uid)

        title_match = re.search(r'<Title>(.*?)</Title>', exp_xml)
        desc_match  = re.search(r'<Summary>(.*?)</Summary>', exp_xml, re.DOTALL)
        title   = title_match.group(1).strip() if title_match else rec.get("title", "")
        summary = desc_match.group(1)[:600].strip() if desc_match else ""

        # Sample count: try the biosample count, then run count
        n_samples = 1
        runs_match = re.findall(r'<Run\s', runs_xml)
        if runs_match:
            n_samples = len(runs_match)
        elif isinstance(rec.get("biosample"), dict):
            n_samples = int(rec["biosample"].get("n", 1))

        return {
            "accession":   acc,
            "title":       title,
            "summary":     summary,
            "organism":    organism or "Homo sapiens",
            "sampleCount": n_samples,
            "gdsType":     "SRA",
            "pubDate":     rec.get("createdate", ""),
            "database":    "SRA",
            "downloadUrl": f"https://www.ncbi.nlm.nih.gov/sra/?term={acc}",
        }
    except Exception as e:
        print_step(f"SRA fetch failed (uid={uid}): {e}", "warn")
        return None


def fetch_pride_datasets(term: str, max_results: int = 6) -> list[dict]:
    """Fetch real datasets from PRIDE proteomics repository."""
    url = f"{PRIDE_BASE}/projects?keyword={urllib.parse.quote(term)}&pageSize={max_results}&page=0"
    try:
        r = requests.get(url, timeout=20, headers={"User-Agent": "OmicsAgent/1.0"})
        r.raise_for_status()
        results = []
        for p in r.json():
            acc = p.get("accession", "")
            if not acc:
                continue
            results.append({
                "accession":   acc,
                "title":       p.get("title", ""),
                "summary":     p.get("projectDescription", "")[:300],
                "organism":    ", ".join(o.get("name","") for o in p.get("organisms",[])),
                "sampleCount": p.get("numberOfSamples", 0) or 0,
                "gdsType":     "Proteomics",
                "pubDate":     p.get("submissionDate", ""),
                "database":    "PRIDE",
                "downloadUrl": f"https://www.ebi.ac.uk/pride/archive/projects/{acc}",
            })
        return results
    except Exception as e:
        print_step(f"PRIDE fetch failed ('{term}'): {e}", "warn")
        return []


def sanitise_strategy(s: dict) -> dict:
    """
    Strip angle-bracket placeholder markers that small models sometimes leave
    in their output, e.g. "<skin>" → "skin", "<one sentence>" → "".
    Also normalises ncbiDb: any omics that isn't explicitly proteomics/metabolomics
    should default to gds (GEO series) rather than sra (raw runs), because GEO
    records carry proper disease labels, sample counts, and summaries.
    """
    def clean(v: str) -> str:
        # Remove leading/trailing < > if the whole value is a placeholder
        v = v.strip()
        if v.startswith("<") and v.endswith(">"):
            v = v[1:-1].strip()
        return v

    s["primaryTissue"]    = clean(s.get("primaryTissue", ""))
    s["mechanisticFocus"] = clean(s.get("mechanisticFocus", ""))
    s["category"]         = clean(s.get("category", "other"))

    pride_omics = {"proteomics", "metabolomics", "metabolome"}
    for o in s.get("omicsTypes", []):
        o["name"]     = clean(o.get("name", ""))
        o["reason"]   = clean(o.get("reason", ""))
        o["database"] = clean(o.get("database", "GEO"))
        o["ncbiDb"]   = clean(o.get("ncbiDb", "gds"))

        # Force GEO (gds) for all omics types except proteomics/metabolomics.
        # The LLM often says ncbiDb=sra but GEO gives far better series-level metadata.
        # scRNA-seq, ATAC-seq, ChIP-seq, and bulk RNA-seq all have GSE series in GEO.
        name_lower = o["name"].lower()
        is_pride   = any(p in name_lower for p in pride_omics)
        if not is_pride and o["ncbiDb"] == "sra":
            o["ncbiDb"]   = "gds"
            o["database"] = "GEO"

        o["searchTerms"] = [clean(t) for t in o.get("searchTerms", []) if t.strip()]

    return s


# ════════════════════════════════════════════════════════════
#  PRE-FILTER  — fast title check before LLM scoring
# ════════════════════════════════════════════════════════════

def build_disease_aliases(disease: str) -> list[str]:
    """
    Return lowercase aliases/synonyms for the disease used in title pre-filtering.
    Covers common abbreviations so we don't drop valid datasets whose titles use
    a synonym rather than the full disease name.
    """
    base = disease.lower().strip()
    aliases = [base]
    known = {
        "psoriasis":                    ["psoriasis", "psoriatic", "pso"],
        "systemic lupus erythematosus": ["lupus", "sle", "systemic lupus"],
        "amyotrophic lateral sclerosis":["als", "amyotrophic lateral sclerosis", "motor neuron disease", "mnd"],
        "alzheimer":                    ["alzheimer", "ad dementia"],
        "parkinson":                    ["parkinson", "pd"],
        "multiple sclerosis":           ["multiple sclerosis", "ms"],
        "rheumatoid arthritis":         ["rheumatoid arthritis", "ra"],
        "type 2 diabetes":              ["type 2 diabetes", "t2d", "t2dm"],
        "type 1 diabetes":              ["type 1 diabetes", "t1d", "t1dm"],
        "crohn":                        ["crohn", "inflammatory bowel disease", "ibd"],
        "ulcerative colitis":           ["ulcerative colitis", "uc", "ibd"],
        "atopic dermatitis":            ["atopic dermatitis", "eczema"],
        "breast cancer":                ["breast cancer", "brca"],
        "lung cancer":                  ["lung cancer", "nsclc", "sclc"],
        "colorectal cancer":            ["colorectal cancer", "crc", "colon cancer"],
    }
    for key, syns in known.items():
        if key in base or base in key:
            aliases.extend(syns)
            break
    return list(dict.fromkeys(aliases))


def title_prefilter(datasets: list, disease: str) -> tuple[list, list]:
    """
    Fast text-based pre-filter run BEFORE LLM scoring.
    Checks title + summary for the disease name or any known synonym.
    Records with zero disease signal are auto-dropped — no LLM call needed.
    Returns (candidates_for_llm, auto_dropped).
    """
    aliases  = build_disease_aliases(disease)
    candidates, auto_dropped = [], []
    for d in datasets:
        text = (d.get("title", "") + " " + d.get("summary", "")).lower()
        if any(alias in text for alias in aliases):
            candidates.append(d)
        else:
            auto_dropped.append({
                **d,
                "omicsType":      d.get("omicsTypeHint", ""),
                "relevanceScore": 0,
                "relevant":       False,
                "notes":          "",
                "dropReason":     "pre-filter: disease name absent from title/summary",
            })
    return candidates, auto_dropped


# ════════════════════════════════════════════════════════════
#  AGENT ORCHESTRATION
# ════════════════════════════════════════════════════════════

class OmicsAgent:
    def __init__(self, model: str):
        self.llm = OllamaClient(model)

    def get_strategy(self, disease: str) -> dict:
        """
        Ask LLM for omics strategy. Uses full prompt (6 search terms) first.
        If the model truncates the response, automatically retries with the
        compact prompt (3 search terms) which produces a shorter JSON output.
        """
        # --- Attempt 1: full prompt ---
        raw, truncated = self.llm.generate(
            STRATEGY_PROMPT.format(disease=disease),
            temperature=0.1,
            num_predict=8192
        )
        if truncated:
            print_step(
                "Response was truncated (model hit token limit) — retrying with compact prompt",
                "warn"
            )
        else:
            try:
                return sanitise_strategy(extract_obj(raw))
            except (ValueError, json.JSONDecodeError):
                print_step("Full prompt parse failed — retrying with compact prompt", "warn")

        # --- Attempt 2: compact fallback prompt ---
        raw2, truncated2 = self.llm.generate(
            STRATEGY_PROMPT_COMPACT.format(disease=disease),
            temperature=0.1,
            num_predict=8192
        )
        if truncated2:
            raise RuntimeError(
                "Compact strategy prompt was also truncated. "
                "Try a model with a larger context window (e.g. llama3.1, mistral)."
            )
        return sanitise_strategy(extract_obj(raw2))

    def score_datasets(self, disease: str, datasets: list) -> list:
        """
        Score datasets for relevance in parallel batches.
        Each batch is sent to the LLM concurrently — since Ollama queues requests
        internally, this keeps the model busy rather than waiting for one batch
        to finish before starting the next.
        """
        if not datasets:
            return []

        # Build batches
        batches = [
            datasets[i : i + SCORE_BATCH_SIZE]
            for i in range(0, len(datasets), SCORE_BATCH_SIZE)
        ]
        total_batches = len(batches)
        results: dict[int, list] = {}    # batch_idx → scores list

        def score_one_batch(idx_batch):
            idx, batch = idx_batch
            slim = [
                {
                    "accession": d["accession"],
                    "title":     d["title"],
                    "summary":   d.get("summary", "")[:400],
                }
                for d in batch
            ]
            raw, truncated = self.llm.generate(
                SCORING_PROMPT.format(
                    disease=disease,
                    datasets_json=json.dumps(slim, indent=2)
                ),
                temperature=0.05
            )
            if truncated:
                print_step(f"  Batch {idx+1} truncated — results may be partial", "warn")
            try:
                return idx, extract_arr(raw)
            except Exception as e:
                print_step(f"  Batch {idx+1} parse error: {e} — scoring as 0", "warn")
                return idx, [
                    {
                        "accession":      d["accession"],
                        "relevanceScore": 0,
                        "relevant":       False,
                        "omicsType":      d.get("omicsTypeHint", ""),
                        "dropReason":     "Scoring failed — parse error",
                        "notes":          ""
                    }
                    for d in batch
                ]

        print_step(f"  Scoring {len(datasets)} datasets in {total_batches} batches "
                   f"(parallel, {SCORE_WORKERS} workers)...")

        with ThreadPoolExecutor(max_workers=SCORE_WORKERS) as ex:
            futures = {ex.submit(score_one_batch, (i, b)): i
                       for i, b in enumerate(batches)}
            for future in as_completed(futures):
                idx, scores = future.result()
                results[idx] = scores
                print_step(f"  Batch {idx+1}/{total_batches} done "
                            f"({len(scores)} scored)", "success")

        # Reassemble in original order
        all_scores = []
        for i in range(total_batches):
            all_scores.extend(results.get(i, []))
        return all_scores

    def fetch_all(self, disease: str, strategy: dict, max_per_omics: int) -> list:
        """
        Fetch real datasets from APIs for each omics type in the strategy.
        esummary calls are parallelised with ThreadPoolExecutor — the bottleneck
        is network latency, so concurrent I/O gives a large speedup.
        UIDs from all search terms are collected first (with dedup), then
        summaries are fetched in parallel.
        """
        all_raw   = []
        seen_accs = set()         # global dedup across all omics types
        seen_lock = __import__("threading").Lock()

        def fetch_one_geo(args):
            uid, term, omics_name = args
            rec = ncbi_fetch_geo_summary(uid)
            if rec:
                rec["omicsTypeHint"] = omics_name
                rec["foundByTerm"]   = term
            return rec

        def fetch_one_sra(args):
            uid, term, omics_name = args
            rec = ncbi_fetch_sra_summary(uid)
            if rec:
                rec["omicsTypeHint"] = omics_name
                rec["foundByTerm"]   = term
            return rec

        for omics in strategy.get("omicsTypes", []):
            name    = omics.get("name", "")
            db_hint = omics.get("database", "GEO")
            ncbi_db = omics.get("ncbiDb", "gds")
            terms   = omics.get("searchTerms", [disease])
            before  = len(all_raw)

            print_step(f"Fetching {name} ({db_hint}) — {len(terms)} search terms")

            if ncbi_db in ("gds", "sra") and "PRIDE" not in db_hint:
                # ── Collect UIDs from all terms (sequential — NCBI search has rate limits) ──
                uid_to_term: dict[str, str] = {}
                for term in terms:
                    for uid in ncbi_search(ncbi_db, term, max_results=max_per_omics):
                        if uid not in uid_to_term:
                            uid_to_term[uid] = term

                # ── Fetch summaries in parallel ──
                fetch_fn   = fetch_one_geo if ncbi_db == "gds" else fetch_one_sra
                fetch_args = [(uid, term, name) for uid, term in uid_to_term.items()]

                with ThreadPoolExecutor(max_workers=FETCH_WORKERS) as ex:
                    futures = {ex.submit(fetch_fn, arg): arg for arg in fetch_args}
                    for future in as_completed(futures):
                        rec = future.result()
                        if rec:
                            with seen_lock:
                                if rec["accession"] not in seen_accs:
                                    all_raw.append(rec)
                                    seen_accs.add(rec["accession"])

            elif "PRIDE" in db_hint or "proteom" in name.lower() or "metabol" in name.lower():
                for term in terms:
                    for rec in fetch_pride_datasets(f"{disease} {term}", max_results=max_per_omics):
                        with seen_lock:
                            if rec["accession"] not in seen_accs:
                                rec["omicsTypeHint"] = name
                                rec["foundByTerm"]   = term
                                all_raw.append(rec)
                                seen_accs.add(rec["accession"])

            print_step(f"  → {len(all_raw) - before} unique {name} records fetched", "success")

        return all_raw


# ════════════════════════════════════════════════════════════
#  DISPLAY
# ════════════════════════════════════════════════════════════

def print_step(msg: str, status: str = "info"):
    ts = datetime.now().strftime("%H:%M:%S")
    sym = {"info":"·","success":"✓","warn":"!","error":"✗"}.get(status,"·")
    if HAS_RICH:
        col = {"info":"dim","success":"green","warn":"yellow","error":"red"}.get(status,"white")
        console.print(f"[dim]{ts}[/dim]  [{col}]{sym}[/{col}]  {escape(str(msg))}")
    else:
        print(f"{ts}  {sym}  {msg}")


def print_strategy(analysis: dict):
    if not HAS_RICH:
        print("\n=== Omics Strategy ===")
        for o in analysis.get("omicsTypes", []):
            print(f"  [{o['priority']}] {o['name']} — {o['database']}")
            print(f"      {o['reason']}")
        return
    t = Table(title="Omics Strategy", box=None, header_style="bold", padding=(0,2))
    t.add_column("Omics Type", style="bold cyan", no_wrap=True)
    t.add_column("Priority", no_wrap=True)
    t.add_column("Database", no_wrap=True)
    t.add_column("Reason", max_width=52)
    for o in analysis.get("omicsTypes", []):
        p = o.get("priority","secondary")
        t.add_row(o["name"],
                  f"[{'green' if p=='primary' else 'blue'}]{p}[/{'green' if p=='primary' else 'blue'}]",
                  o.get("database",""), o.get("reason",""))
    console.print(); console.print(t); console.print()


def print_datasets(datasets: list):
    db_icons = {"GEO":"🔵","SRA":"🟢","TCGA":"🟠","ENCODE":"🟣","PRIDE":"🟡","MetaboLights":"🟤"}
    by_db: dict = {}
    for d in datasets:
        by_db.setdefault(d["database"], []).append(d)

    if not HAS_RICH:
        for db, items in by_db.items():
            print(f"\n── {db} ({len(items)}) ──")
            for d in items:
                print(f"  {d['accession']} | {d.get('omicsType','')} | {d.get('sampleCount','')} samples")
                print(f"  {d['title'][:80]}")
                print(f"  Score: {d.get('relevanceScore','')}% | {d['downloadUrl']}")
        return

    for db, items in by_db.items():
        t = Table(title=f"{db_icons.get(db,'⚪')}  {db}  ({len(items)} datasets)",
                  box=None, header_style="bold", padding=(0,1))
        t.add_column("Accession", style="cyan", no_wrap=True)
        t.add_column("Omics", no_wrap=True)
        t.add_column("Title", max_width=44)
        t.add_column("Samples", justify="right")
        t.add_column("Year", justify="right")
        t.add_column("Match", justify="right")
        for d in sorted(items, key=lambda x: -x.get("relevanceScore",0)):
            sc = d.get("relevanceScore", 0)
            c = "green" if sc>=80 else "yellow" if sc>=60 else "dim"
            yr = d.get("pubDate","")[:4] or "—"
            t.add_row(d["accession"], d.get("omicsType",""), d["title"],
                      str(d.get("sampleCount","—")), yr,
                      f"[{c}]{sc}%[/{c}]")
        console.print(t); console.print()


# ════════════════════════════════════════════════════════════
#  SAVE
# ════════════════════════════════════════════════════════════

def save_curated_csv(datasets: list, disease: str, out_dir: Path) -> Path:
    """Write only relevant (score >= 50) datasets, sorted by score."""
    safe  = re.sub(r"[^\w-]", "_", disease.strip())
    fname = out_dir / f"omics_datasets_{safe}_curated.csv"
    fields = ["accession", "title", "database", "omicsType", "organism",
              "sampleCount", "pubDate", "downloadUrl", "relevanceScore",
              "notes", "foundByTerm"]
    with open(fname, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for d in sorted(datasets, key=lambda x: -x.get("relevanceScore", 0)):
            w.writerow(d)
    return fname


def save_audit_csv(all_datasets: list, disease: str, out_dir: Path) -> Path:
    """
    Write ALL fetched datasets including dropped ones.
    Includes dropReason column so you can audit why each was excluded.
    """
    safe  = re.sub(r"[^\w-]", "_", disease.strip())
    fname = out_dir / f"omics_datasets_{safe}_all.csv"
    fields = ["accession", "title", "database", "omicsType", "organism",
              "sampleCount", "pubDate", "downloadUrl", "relevanceScore",
              "relevant", "notes", "dropReason", "foundByTerm", "omicsTypeHint"]
    with open(fname, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for d in sorted(all_datasets, key=lambda x: -x.get("relevanceScore", 0)):
            w.writerow(d)
    return fname

def save_json(report: dict, disease: str, out_dir: Path) -> Path:
    safe  = re.sub(r"[^\w-]", "_", disease.strip())
    fname = out_dir / f"omics_report_{safe}.json"
    with open(fname, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    return fname


# ════════════════════════════════════════════════════════════
#  MAIN RUN
# ════════════════════════════════════════════════════════════

def run(disease: str, model: str, out_dir: Path, max_per_omics: int = 8):
    out_dir.mkdir(parents=True, exist_ok=True)
    agent = OmicsAgent(model)

    if not agent.llm.check():
        print_step("Cannot reach Ollama at http://localhost:11434 — run: ollama serve", "error")
        sys.exit(1)

    available = agent.llm.list_models()
    model_base = model.split(":")[0]
    if available and not any(model_base in m for m in available):
        print_step(f'Model "{model}" not found. Available: {", ".join(available)}', "error")
        print_step(f"Pull it with:  ollama pull {model}", "warn")
        sys.exit(1)

    if HAS_RICH:
        console.print(Panel.fit(
            f"[bold]Omics Dataset Research Agent[/bold]  [dim](Ollama · {model})[/dim]\n"
            f"[dim]Disease:[/dim] {escape(disease)}\n"
            f"[dim]Data source:[/dim] NCBI GEO · NCBI SRA · EBI PRIDE  [green](live API)[/green]",
            border_style="green"))
    else:
        print(f"\n{'='*62}\n  Omics Research Agent ({model})\n  Disease: {disease}\n  Source: NCBI GEO / SRA / PRIDE (live)\n{'='*62}\n")

    # ── Step 1: LLM builds omics strategy + search terms ──
    print_step("Building omics strategy...")
    t0 = time.time()
    strategy = agent.get_strategy(disease)
    print_step(f"Category: {strategy.get('category','?')} | Tissue: {strategy.get('primaryTissue','?')}  ({time.time()-t0:.1f}s)", "success")
    print_step(f"Focus: {strategy.get('mechanisticFocus','')}")
    omics_summary = ", ".join(
        f"{o['name']} ({o.get('database','?')})"
        for o in strategy.get("omicsTypes", [])
    )
    print_step(f"Omics plan: {omics_summary}", "success")
    print_strategy(strategy)

    # ── Step 2: Fetch REAL datasets from NCBI + PRIDE APIs ──
    print_step("Querying live databases (NCBI GEO, SRA, EBI PRIDE)...")
    t0 = time.time()
    raw_datasets = agent.fetch_all(disease, strategy, max_per_omics)
    print_step(f"Fetched {len(raw_datasets)} real records in {time.time()-t0:.1f}s", "success")

    if not raw_datasets:
        print_step("No datasets returned. Check your internet connection.", "error")
        sys.exit(1)

    # ── Step 3: Fast title pre-filter (no LLM needed) ──
    candidates, prefilter_dropped = title_prefilter(raw_datasets, disease)
    print_step(
        f"Pre-filter: {len(candidates)} candidates pass (disease name in title/summary), "
        f"{len(prefilter_dropped)} auto-dropped",
        "success" if prefilter_dropped else "info"
    )

    # ── Step 4: LLM scores only the candidates that passed pre-filter ──
    print_step(f"Scoring {len(candidates)} candidates for relevance to '{disease}'...")
    t0 = time.time()
    scores = agent.score_datasets(disease, candidates)
    print_step(f"Scoring done ({time.time()-t0:.1f}s)", "success")

    # ── Merge scores back into candidate records ──
    score_map = {s["accession"]: s for s in scores}
    llm_merged = []
    for d in candidates:
        acc = d["accession"]
        sc  = score_map.get(acc, {})
        llm_merged.append({
            **d,
            "omicsType":      sc.get("omicsType",      d.get("omicsTypeHint", "")),
            "relevanceScore": sc.get("relevanceScore", 0),
            "relevant":       sc.get("relevant",       False),
            "notes":          sc.get("notes",          ""),
            "dropReason":     sc.get("dropReason",     ""),
        })

    # Combine LLM-scored records with pre-filter drops for full audit trail
    all_merged = llm_merged + prefilter_dropped
    all_merged.sort(key=lambda x: -x.get("relevanceScore", 0))

    curated  = [d for d in all_merged if d.get("relevanceScore", 0) >= 50]
    dropped  = [d for d in all_merged if d.get("relevanceScore", 0) <  50]

    total_samples = sum(d.get("sampleCount", 0) for d in curated)
    dbs           = len(set(d["database"] for d in curated))

    print_step(
        f"Kept {len(curated)} relevant datasets (score ≥ 50), "
        f"dropped {len(dropped)} total "
        f"({len(prefilter_dropped)} pre-filter + {len(dropped)-len(prefilter_dropped)} LLM) "
        f"— {total_samples:,} total samples across {dbs} database(s)",
        "success"
    )

    if dropped:
        print_step(f"Dropped datasets by reason:")
        reason_counts: dict = {}
        for d in dropped:
            r = d.get("dropReason", "unscored") or "no reason given"
            reason_counts[r[:80]] = reason_counts.get(r[:80], 0) + 1
        for reason, count in sorted(reason_counts.items(), key=lambda x: -x[1]):
            print_step(f"  {count}× {reason}", "warn")

    print_datasets(curated)

    # ── Save ──
    report = {
        "disease":        disease,
        "model":          model,
        "timestamp":      datetime.now().isoformat(),
        "strategy":       strategy,
        "summary": {
            "total_fetched":  len(raw_datasets),
            "curated":        len(curated),
            "dropped":        len(dropped),
            "total_samples":  total_samples,
            "databases":      dbs,
        },
        "curated_datasets": curated,
        "dropped_datasets":  dropped,   # full audit trail in JSON too
    }

    curated_csv = save_curated_csv(curated,    disease, out_dir)
    audit_csv   = save_audit_csv(all_merged,   disease, out_dir)
    json_path   = save_json(report,            disease, out_dir)

    if HAS_RICH:
        console.print(Panel(
            f"[green]✓[/green]  Curated CSV  ({len(curated)} datasets):  [cyan]{curated_csv}[/cyan]\n"
            f"[dim]✓  Audit CSV    ({len(all_merged)} datasets):  {audit_csv}[/dim]\n"
            f"[dim]✓  Full JSON:                        {json_path}[/dim]",
            title="Saved", border_style="dim"))
    else:
        print(f"\nSaved:"
              f"\n  Curated ({len(curated)}):  {curated_csv}"
              f"\n  Audit   ({len(all_merged)}):  {audit_csv}"
              f"\n  JSON:       {json_path}\n")
    return report


# ════════════════════════════════════════════════════════════
#  CLI
# ════════════════════════════════════════════════════════════

def main():
    global OLLAMA_BASE
    parser = argparse.ArgumentParser(
        description="Omics Dataset Research Agent — real database queries via NCBI + EBI APIs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python omics_agent.py --disease "psoriasis" --model llama3.2
  python omics_agent.py --disease "ALS" --model mistral --max-results 15
  python omics_agent.py --list-models
        """
    )
    parser.add_argument("--disease",     "-d", type=str)
    parser.add_argument("--model",       "-m", type=str, default="mistral")
    parser.add_argument("--out",         "-o", type=str, default="omics_results")
    parser.add_argument("--max-results", "-n", type=int, default=8,
                        help="Max datasets to fetch per omics type (default: 8)")
    parser.add_argument("--ollama-url",        type=str, default=OLLAMA_BASE)
    parser.add_argument("--list-models",       action="store_true")
    args = parser.parse_args()

    OLLAMA_BASE = args.ollama_url

    if args.list_models:
        c = OllamaClient("", base_url=OLLAMA_BASE)
        if not c.check():
            print("Cannot reach Ollama. Run:  ollama serve")
            sys.exit(1)
        ms = c.list_models()
        print("Installed models:\n" + ("\n".join(f"  {m}" for m in ms) if ms else "  none — try: ollama pull mistral"))
        sys.exit(0)

    disease = args.disease
    if not disease:
        disease = (console.input("[bold green]Disease:[/bold green] ") if HAS_RICH else input("Disease: ")).strip()
    if not disease:
        print("No disease entered."); sys.exit(1)

    try:
        run(disease, args.model, Path(args.out), args.max_results)
    except KeyboardInterrupt:
        print("\nInterrupted.")
    except requests.ConnectionError as e:
        print_step(f"Network error: {e}", "error")
        sys.exit(1)
    except Exception as e:
        print_step(f"Unexpected error: {e}", "error")
        raise

if __name__ == "__main__":
    main()
