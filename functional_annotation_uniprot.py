#!/usr/bin/env python3
"""
functional_annotation_uniprot.py

Annotate protein list with GO Biological Process (BP) and GO Molecular Function (MF)
using UniProt REST API. Includes robust normalization of the long table so plotting
functions always receive 'ontology' and 'go_term' columns.

Usage:
    python3 functional_annotation_uniprot.py -i "my_list.xlsx"

Outputs (in OUTDIR):
 - annotations_full.csv
 - annotations_long.csv
 - annotations_full.xlsx   (if xlsxwriter/openpyxl available)
 - go_annotations_long.csv
 - go_annotations_full.csv
 - go_summary.json
 - publication-grade plots (PNG + optional PDF/SVG)

Dependencies:
    pip install requests pandas matplotlib tqdm openpyxl xlsxwriter
"""

from __future__ import annotations
import argparse
import os
import time
import json
import re
from typing import List, Dict, Any, Optional, Tuple
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Optional progress bar
try:
    from tqdm import tqdm
except Exception:
    def tqdm(x, **_k): return x

# -------------------------
# Config
# -------------------------
UNIPROT_BASE = "https://rest.uniprot.org"
SLEEP_BETWEEN_REQUESTS = 0.3
CACHE_UNIPROT = "uniprot_cache.json"
REQUEST_TIMEOUT = 30
MAX_RETRIES = 3

# Publication-style rcParams
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans"],
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 300,
    "figure.facecolor": "white",
    "axes.grid": False,
})

# -------------------------
# Utility functions
# -------------------------
def load_json_cache(path: str) -> Dict[str, Any]:
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as fh:
                return json.load(fh)
        except Exception:
            return {}
    return {}

def save_json_cache(obj: Dict[str, Any], path: str):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, ensure_ascii=False)

def request_with_retries(url: str, params: dict = None, headers: dict = None, timeout: int = REQUEST_TIMEOUT) -> Tuple[int, str]:
    headers = headers or {}
    attempt = 0
    last_exc = None
    while attempt < MAX_RETRIES:
        try:
            r = requests.get(url, params=params, headers=headers, timeout=timeout)
            return r.status_code, r.text
        except Exception as e:
            last_exc = e
            attempt += 1
            time.sleep(1.5 ** attempt)
    if last_exc:
        raise last_exc
    raise RuntimeError(f"Failed request to {url}")

# -------------------------
# Input reading (heuristic)
# -------------------------
ID_HINTS = ["uniprot", "accession", "id", "protein", "gene", "symbol", "entry", "tair"]

def read_input_ids(input_path: str, sheet: Optional[str | int] = None) -> Tuple[List[str], pd.DataFrame, str]:
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"File not found: {input_path}")

    if input_path.lower().endswith((".xls", ".xlsx")):
        raw = pd.read_excel(input_path, sheet_name=sheet, engine="openpyxl")
        if isinstance(raw, dict):
            chosen_sheet = sheet if sheet is not None and sheet in raw else list(raw.keys())[0]
            df = raw[chosen_sheet]
            if sheet is None:
                print(f"Multiple sheets found — using first: '{chosen_sheet}'")
        else:
            df = raw
    elif input_path.lower().endswith(".csv"):
        df = pd.read_csv(input_path)
    elif input_path.lower().endswith((".tsv", ".txt")):
        df = pd.read_csv(input_path, sep="\t")
    else:
        raise ValueError("Unsupported file type. Use xlsx/csv/tsv/txt.")

    if df is None or df.shape[0] == 0:
        raise ValueError("Input file or sheet is empty.")

    cols = df.columns.tolist()
    chosen_col = None
    for c in cols:
        lc = str(c).lower()
        if any(h in lc for h in ID_HINTS):
            chosen_col = c
            break

    if chosen_col is None:
        # heuristic: find column with accession-like strings
        for c in cols:
            sample = df[c].astype(str).dropna().head(20).tolist()
            if any(re.match(r"^[A-NR-Z][0-9A-Z]{5}$", s) or re.match(r"^[A-Z0-9]{6,10}$", s) for s in sample):
                chosen_col = c
                print(f"Heuristically picked '{chosen_col}' as ID column.")
                break

    if chosen_col is None:
        chosen_col = cols[0]
        print(f"Defaulting to first column '{chosen_col}' as ID column (no heuristic match).")

    ids = df[chosen_col].astype(str).str.strip().replace({"": None, "nan": None, "NaN": None}).dropna().unique().tolist()
    print(f"Using column '{chosen_col}' with {len(ids)} unique IDs (showing up to 10): {ids[:10]}")
    return ids, df, chosen_col

# -------------------------
# UniProt fetching & parsing
# -------------------------
def fetch_uniprot_json(accession: str, cache: Dict[str, Any]) -> Dict[str, Any]:
    if accession in cache:
        return cache[accession]
    url = f"{UNIPROT_BASE}/uniprotkb/{accession}.json"
    try:
        status, text = request_with_retries(url)
        if status == 200:
            data = json.loads(text)
            out = {"found": True, "raw": data}
        elif status == 404:
            out = {"found": False}
        else:
            try:
                out = {"found": True, "raw": json.loads(text)}
            except Exception:
                out = {"found": False}
    except Exception as e:
        out = {"found": False, "error": str(e)}
    cache[accession] = out
    time.sleep(SLEEP_BETWEEN_REQUESTS)
    return out

def parse_uniprot_for_annotations(entry_json: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract basic annotation fields and GO BP/MF lists (readable string with optional GO ID).
    """
    out = {
        "accession": None,
        "primary_gene": None,
        "protein_name": None,
        "go_bp": [],
        "go_mf": [],
        "sequence_length": None,
        "sequence": None,
    }
    if not entry_json or not entry_json.get("found", False):
        return out

    data = entry_json.get("raw", {})

    # accession and names
    out["accession"] = data.get("primaryAccession") or out["accession"]
    prot = data.get("proteinDescription", {})
    try:
        rec = prot.get("recommendedName", {})
        pname = rec.get("fullName", {}).get("value")
        if not pname:
            anames = prot.get("alternativeNames", [])
            if anames:
                pname = anames[0].get("fullName", {}).get("value")
        out["protein_name"] = pname
    except Exception:
        out["protein_name"] = None

    gene_obj = data.get("genes", [])
    if gene_obj:
        try:
            out["primary_gene"] = gene_obj[0].get("geneName", {}).get("value") or gene_obj[0].get("primaryName", {}).get("value")
        except Exception:
            out["primary_gene"] = None

    seqobj = data.get("sequence", {})
    if seqobj:
        out["sequence_length"] = seqobj.get("length")
        out["sequence"] = seqobj.get("value")

    # collect GO from dbReferences / uniProtKBCrossReferences
    xrefs = data.get("uniProtKBCrossReferences", []) or data.get("dbReferences", []) or []
    for xr in xrefs:
        dbtype = (xr.get("database") or xr.get("type") or "").upper()
        if dbtype == "GO":
            props = xr.get("properties", []) or []
            term_val = None
            go_id = None
            for p in props:
                key = str(p.get("key", "")).lower()
                if key in {"term", "name"}:
                    term_val = p.get("value")
                if key in {"id", "go id", "accession"}:
                    go_id = p.get("value")
            if term_val:
                m = re.match(r"([PFC]):\s*(.+)", term_val)
                if m:
                    code, name = m.group(1), m.group(2)
                    name_full = f"{name} ({go_id})" if go_id else name
                    if code == "P":
                        out["go_bp"].append(name_full)
                    elif code == "F":
                        out["go_mf"].append(name_full)
                else:
                    name_full = f"{term_val} ({go_id})" if go_id else term_val
                    out["go_bp"].append(name_full)

    # Deduplicate
    out["go_bp"] = sorted(list(dict.fromkeys([g for g in out["go_bp"] if g])))
    out["go_mf"] = sorted(list(dict.fromkeys([g for g in out["go_mf"] if g])))

    return out

# -------------------------
# Annotation pipeline
# -------------------------
def annotate_list(ids: List[str], cache_file: str = CACHE_UNIPROT, pause: float = SLEEP_BETWEEN_REQUESTS) -> pd.DataFrame:
    cache = load_json_cache(cache_file)
    results = []
    for raw_id in tqdm(ids, desc="Annotating", unit="id"):
        identifier = str(raw_id).strip()
        if not identifier:
            continue
        if identifier in cache and isinstance(cache[identifier], dict) and cache[identifier].get("parsed"):
            parsed = cache[identifier]["parsed"]
            results.append({"input_id": identifier, **parsed})
            continue
        # fetch (if looks like accession use directly, else still attempt)
        entry = fetch_uniprot_json(identifier, cache)
        parsed = parse_uniprot_for_annotations(entry)
        cache[identifier] = {"raw_fetch": entry, "parsed": parsed}
        results.append({"input_id": identifier, **parsed})
        time.sleep(pause)
    save_json_cache(cache, cache_file)
    df = pd.DataFrame(results)
    # Ensure columns
    for col in ["accession","primary_gene","protein_name","go_bp","go_mf","sequence_length","sequence"]:
        if col not in df.columns:
            df[col] = None
    return df

def explode_long_table(df_full: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in df_full.iterrows():
        input_id = r["input_id"]
        for go in r.get("go_bp") or []:
            rows.append({"input_id": input_id, "type": "GO_BP", "value": go})
        for go in r.get("go_mf") or []:
            rows.append({"input_id": input_id, "type": "GO_MF", "value": go})
    return pd.DataFrame(rows)

# -------------------------
# Long-table normalization
# -------------------------
def normalize_long_table(df_long: Optional[pd.DataFrame]) -> pd.DataFrame:
    cols = ["input_id", "ontology", "go_term"]
    if df_long is None:
        return pd.DataFrame(columns=cols)
    df = df_long.copy()
    # Accept old keys: type->ontology, value->go_term
    if "type" in df.columns and "value" in df.columns:
        df = df.rename(columns={"type":"ontology","value":"go_term"})
    elif "type" in df.columns and "go_term" in df.columns:
        df = df.rename(columns={"type":"ontology"})
    elif "ontology" in df.columns and "value" in df.columns:
        df = df.rename(columns={"value":"go_term"})
    # If only 'value' - assume it's go_term and set ontology None
    if "value" in df.columns and "go_term" not in df.columns:
        df = df.rename(columns={"value":"go_term"})
        if "ontology" not in df.columns:
            df["ontology"] = None
    for c in cols:
        if c not in df.columns:
            df[c] = None
    # Standardize and drop rows lacking go_term
    df["go_term"] = df["go_term"].astype(str).replace({"None": None})
    df["ontology"] = df["ontology"].astype(str).replace({"None": None})
    df = df[~(df["go_term"].isna() | (df["go_term"] == "None"))].copy()
    df["go_term"] = df["go_term"].astype(str).str.strip()
    df["ontology"] = df["ontology"].astype(str).str.strip()
    return df[cols]

# -------------------------
# Plotting (publication grade)
# -------------------------
def wrap_label(text: str, width: int = 70) -> str:
    import textwrap
    return "\n".join(textwrap.wrap(str(text), width=width))

def save_figs(fig, base: str, fmt: str):
    if fmt in ("png","all"):
        fig.savefig(base + ".png", dpi=300, bbox_inches="tight")
    if fmt in ("pdf","all"):
        fig.savefig(base + ".pdf", bbox_inches="tight")
    if fmt in ("svg","all"):
        fig.savefig(base + ".svg", bbox_inches="tight")

def plot_single_go_enhanced(df_long: pd.DataFrame, outdir: str, ontology: str, top_n: int = 20, fmt: str = "all"):
    df_sub = df_long[df_long["ontology"] == ontology]
    if df_sub.empty:
        print(f"No GO {ontology} mappings to plot.")
        return
    counts = df_sub["go_term"].value_counts().head(top_n)
    labels_raw = counts.index.tolist()
    values = counts.values
    labels = [wrap_label(l, width=70) for l in labels_raw]
    total_proteins = df_long["input_id"].nunique()
    cum_counts = np.cumsum(values)
    cum_percent = 100.0 * cum_counts / total_proteins if total_proteins>0 else np.zeros_like(cum_counts)
    n = len(labels)
    height = max(3.2, 0.35*n)
    fig = plt.figure(figsize=(12, height))
    gs = fig.add_gridspec(1,3, width_ratios=[3,0.05,1], wspace=0.12)
    ax = fig.add_subplot(gs[0,0])
    ax_in = fig.add_subplot(gs[0,2])
    y = np.arange(n)
    cmap = plt.get_cmap("viridis" if ontology=="BP" else "plasma")
    colors = cmap(np.linspace(0.15,0.85,n))
    ax.barh(y, values[::-1], color=colors[::-1], edgecolor="black", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([labels[i] for i in range(n-1,-1,-1)])
    ax.invert_yaxis()
    ax.set_xlabel("Number of proteins annotated")
    ax.set_title("Top GO Biological Processes" if ontology=="BP" else "Top GO Molecular Functions")
    ax.xaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.6)
    for i,v in enumerate(values[::-1]):
        ax.text(v + max(1, max(values)*0.01), i, f"{v}", va="center", fontsize=9)
    ax_in.plot(range(1,n+1), cum_percent[::-1], marker="o")
    ax_in.set_ylim(0,100)
    ax_in.set_xlim(1,n)
    ax_in.set_xlabel("Rank")
    ax_in.set_ylabel("Cumulative % proteins")
    plt.tight_layout()
    base = os.path.join(outdir, f"top_go_{ontology.lower()}_top{n}")
    save_figs(fig, base, fmt)
    plt.close(fig)
    counts_df = pd.DataFrame({"go_term": labels_raw, "count": values, "cum_count": cum_counts, "cum_percent": cum_percent})
    counts_df.to_csv(base + "_counts.csv", index=False)
    print(f"Saved GO {ontology} plot: {base} (+ counts CSV)")

def plot_combined_go(df_long: pd.DataFrame, outdir: str, top_n: int = 15, fmt: str = "all"):
    df_bp = df_long[df_long["ontology"]=="BP"]
    df_mf = df_long[df_long["ontology"]=="MF"]
    counts_bp = df_bp["go_term"].value_counts().head(top_n)
    counts_mf = df_mf["go_term"].value_counts().head(top_n)
    n_bp = counts_bp.shape[0]
    n_mf = counts_mf.shape[0]
    if n_bp==0 and n_mf==0:
        print("No BP/MF terms for combined figure.")
        return
    height = max(4, 0.35*max(n_bp,n_mf,5))
    fig, axes = plt.subplots(1,2, figsize=(14,height))
    if n_bp>0:
        labels_bp = [wrap_label(l,60) for l in counts_bp.index.tolist()]
        vals_bp = counts_bp.values
        ybp = np.arange(len(labels_bp))
        axes[0].barh(ybp, vals_bp[::-1], color=plt.get_cmap("viridis")(np.linspace(0.2,0.8,len(labels_bp)))[::-1], edgecolor="black", linewidth=0.5)
        axes[0].set_yticks(ybp)
        axes[0].set_yticklabels([labels_bp[i] for i in range(len(labels_bp)-1,-1,-1)])
        axes[0].invert_yaxis()
        axes[0].set_title("GO Biological Processes")
        axes[0].set_xlabel("Number of proteins")
        for i,v in enumerate(vals_bp[::-1]):
            axes[0].text(v + max(1, max(vals_bp)*0.01), i, f"{v}", va="center", fontsize=9)
        axes[0].xaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.6)
    else:
        axes[0].axis("off")
    if n_mf>0:
        labels_mf = [wrap_label(l,60) for l in counts_mf.index.tolist()]
        vals_mf = counts_mf.values
        ymf = np.arange(len(labels_mf))
        axes[1].barh(ymf, vals_mf[::-1], color=plt.get_cmap("plasma")(np.linspace(0.2,0.8,len(labels_mf)))[::-1], edgecolor="black", linewidth=0.5)
        axes[1].set_yticks(ymf)
        axes[1].set_yticklabels([labels_mf[i] for i in range(len(labels_mf)-1,-1,-1)])
        axes[1].invert_yaxis()
        axes[1].set_title("GO Molecular Functions")
        axes[1].set_xlabel("Number of proteins")
        for i,v in enumerate(vals_mf[::-1]):
            axes[1].text(v + max(1, max(vals_mf)*0.01), i, f"{v}", va="center", fontsize=9)
        axes[1].xaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.6)
    else:
        axes[1].axis("off")
    plt.tight_layout()
    base = os.path.join(outdir, f"combined_go_bp_mf_top{top_n}")
    save_figs(fig, base, fmt)
    plt.close(fig)
    print(f"Saved combined GO BP/MF figure: {base}")

# -------------------------
# Main CLI
# -------------------------
def main():
    p = argparse.ArgumentParser(description="Annotate protein list via UniProt and plot GO BP/MF.")
    p.add_argument("-i","--input", required=True, help="Input file (xlsx/csv/tsv).")
    p.add_argument("-s","--sheet", default=None, help="Excel sheet name/index (optional).")
    p.add_argument("-o","--outdir", default="annotation_results", help="Output directory.")
    p.add_argument("--top", type=int, default=20, help="Top N terms to plot.")
    p.add_argument("--format", choices=["png","pdf","svg","all"], default="all", help="Figure output format(s).")
    p.add_argument("--no-plot", action="store_true", help="Do annotation only; skip plotting.")
    args = p.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    ids, df_input, chosen_col = read_input_ids(args.input, sheet=args.sheet)
    df_annotations = annotate_list(ids, cache_file=CACHE_UNIPROT, pause=SLEEP_BETWEEN_REQUESTS)

    # Save wide results
    full_csv = os.path.join(outdir, "annotations_full.csv")
    full_xlsx = os.path.join(outdir, "annotations_full.xlsx")
    df_annotations.to_csv(full_csv, index=False)
    try:
        df_annotations.to_excel(full_xlsx, index=False)
    except Exception:
        pass
    print(f"Saved: {full_csv} (and Excel if supported)")

    # Create long table and normalize
    df_long = explode_long_table(df_annotations)
    # save original long as well
    long_csv = os.path.join(outdir, "annotations_long_raw.csv")
    df_long.to_csv(long_csv, index=False)
    # normalize to expected columns for plotting
    df_long_norm = normalize_long_table(df_long)
    go_long_csv = os.path.join(outdir, "go_annotations_long.csv")
    df_long_norm.to_csv(go_long_csv, index=False)

    # Save full go annotations (wide)
    go_full_csv = os.path.join(outdir, "go_annotations_full.csv")
    df_annotations[["input_id","accession","primary_gene","protein_name","go_bp","go_mf"]].to_csv(go_full_csv, index=False)

    # Summary
    summary = {
        "n_input_ids": len(ids),
        "n_annotated": int(df_annotations['accession'].notna().sum()),
        "n_with_go_terms": int(df_long_norm['input_id'].nunique()),
        "n_go_mappings": int(df_long_norm.shape[0])
    }
    with open(os.path.join(outdir, "go_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print("Saved summary and annotation files.")

    # Plotting
    if not args.no_plot and not df_long_norm.empty:
        fmt = args.format
        plot_single_go_enhanced(df_long_norm, outdir, ontology="BP", top_n=args.top, fmt=fmt)
        plot_single_go_enhanced(df_long_norm, outdir, ontology="MF", top_n=args.top, fmt=fmt)
        plot_combined_go(df_long_norm, outdir, top_n=min(args.top,15), fmt=fmt)
    else:
        if df_long_norm.empty:
            print("No GO annotations available to plot.")
        else:
            print("Plotting skipped (--no-plot).")

    print("Done. Outputs in:", outdir)
    print("Summary:", json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
