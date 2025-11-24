#!/usr/bin/env python3
"""
Modified: relaxed interactors now also include proteins with frequency difference >= RELAXED_FREQ_DIFF_THRESHOLD

Original script: presence/absence curation (3 vs 3) with UniProt gene name mapping.

Usage: same as before. Outputs under ./curation_outputs/
"""
import os
from pathlib import Path
import time
import requests
import pandas as pd
import numpy as np
import io
import json

# Optional plotting/statistics libraries (import only if used)
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import importlib

# Optional venn plotting
_spec = importlib.util.find_spec("matplotlib_venn")
if _spec is not None:
    venn2 = importlib.import_module("matplotlib_venn").venn2
    HAVE_VENN = True
else:
    venn2 = None
    HAVE_VENN = False

# ---------------------------
# USER: Update these paths / settings
# ---------------------------
sample_files = [
    "18_GFP1_proteinList.xlsx",
    "18_GFP2_proteinList.xlsx",
    "18_GFP3_proteinList.xlsx",
]

control_files = [
    "Mito_GFP1_proteinList.xlsx",
    "Mito_GFP2_proteinList.xlsx",
    "Mito_GFP3_proteinList.xlsx",
]

OUTDIR = Path("curation_outputs")
OUTDIR.mkdir(exist_ok=True, parents=True)

# ---------------------------
# PARAMETERS / THRESHOLDS
# ---------------------------
MIN_UNIQUE_PEPTIDES = 1            # minimal peptide evidence to keep during curation
STRONG_PEPTIDE_THRESHOLD = 1       # mark as strong evidence if >= this
STRICT_SAMPLE_REPLICATES = 2       # for strict rule: >=2 of 3 sample reps
RELAXED_CONTROL_MAX = 1            # for relaxed rule: <=1 of 3 control reps
RELAXED_FREQ_DIFF_THRESHOLD = 0.33  # NEW: include proteins with freq_diff >= this into relaxed set
ALPHA = 0.05                       # significance threshold if running stats

# Control behaviour
RUN_STATS = True   # set True to run Fisher + BH correction
RUN_PLOTS = False   # plotting optional; default off for curation-focused run

# UniProt mapping config (gene names only)
UNIPROT_BATCH_SIZE = 200
UNIPROT_SLEEP = 0.5
UNIPROT_ENDPOINT = "https://rest.uniprot.org/uniprotkb/search"

# ---------------------------
# Helper: detect column names
# ---------------------------
def detect_id_col(df):
    candidates = [c for c in df.columns if any(k in c.lower() for k in
        ["accession", "master", "protein id", "protein ids", "uniprot", "entry"]) ]
    if candidates:
        for pref in ["accession", "master", "protein id", "protein ids", "uniprot", "entry", "protein"]:
            for c in candidates:
                if pref in c.lower():
                    return c
        return candidates[0]
    for c in df.columns:
        if pd.api.types.is_string_dtype(df[c]) or df[c].dtype == object:
            return c
    return df.columns[0]

def detect_unique_pep_col(df):
    candidates = [c for c in df.columns if any(k in c.lower() for k in
        ["unique peptide", "unique peptides", "uniquepep", "# unique", "peptide count", "psm", "peptides"]) ]
    if candidates:
        for pref in ["unique peptide", "unique peptides", "# unique", "uniquepep", "peptides", "psm"]:
            for c in candidates:
                if pref in c.lower():
                    return c
        return candidates[0]
    numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if numeric_cols:
        best = min(numeric_cols, key=lambda c: abs(df[c].dropna().median() - 2))
        return best
    return None

def detect_contaminant_col(df):
    candidates = [c for c in df.columns if "contaminant" in c.lower() or c.lower() in
                  ("iscontaminant", "is_contaminant", "contaminant_flag", "reverse", "decoy")]
    if candidates:
        for pref in ["contaminant", "iscontaminant", "is_contaminant"]:
            for c in candidates:
                if pref in c.lower():
                    return c
        return candidates[0]
    return None

# ---------------------------
# Contaminant marking (fallback)
# ---------------------------
def mark_contaminant_heuristic(row, accession_col, desc_cols):
    acc = str(row.get(accession_col, "")).upper()
    if acc.startswith("CON__") or acc.startswith("REV__") or acc.startswith("DECOY"):
        return True
    for dcol in desc_cols:
        txt = str(row.get(dcol, "")).lower()
        if "contaminant" in txt or "reverse" in txt or "decoy" in txt:
            return True
    return False

# ---------------------------
# Read sheet helper
# ---------------------------
def read_first_sheet_maybe(path):
    path = Path(path)
    if path.suffix.lower() in [".xls", ".xlsx"]:
        xls = pd.read_excel(path, sheet_name=None)
        first = list(xls.keys())[0]
        return xls[first]
    else:
        return pd.read_csv(path)

# ---------------------------
# Clean single replicate
# ---------------------------
def clean_replicate(path, outdir, sample_name=None, force_id_col=None, force_pep_col=None, force_contaminant_col=None):
    df = read_first_sheet_maybe(path)
    id_col = force_id_col if force_id_col else detect_id_col(df)
    pep_col = force_pep_col if force_pep_col else detect_unique_pep_col(df)
    contaminant_col = force_contaminant_col if force_contaminant_col else detect_contaminant_col(df)
    desc_candidates = [c for c in df.columns if "desc" in c.lower() or "protein" in c.lower() or "gene" in c.lower()]

    df = df.copy()
    df['__accession__'] = df[id_col].astype(str).str.strip()
    if pep_col:
        df['__unique_pep__'] = pd.to_numeric(df[pep_col], errors='coerce').fillna(0).astype(int)
    else:
        df['__unique_pep__'] = 1

    if contaminant_col:
        def parse_bool(x):
            if pd.isna(x):
                return False
            if isinstance(x, bool):
                return bool(x)
            s = str(x).strip().lower()
            return s in ("true", "1", "yes", "y", "t")
        df['__is_contaminant__'] = df[contaminant_col].apply(parse_bool)
    else:
        df['__is_contaminant__'] = df.apply(lambda r: mark_contaminant_heuristic(r, '__accession__', desc_candidates), axis=1)

    df_clean = df[~df['__is_contaminant__']].copy()
    df_clean = df_clean[df_clean['__unique_pep__'] >= MIN_UNIQUE_PEPTIDES].copy()
    df_clean = df_clean.sort_values('__unique_pep__', ascending=False).drop_duplicates('__accession__', keep='first')
    df_clean['__strong_evidence__'] = df_clean['__unique_pep__'] >= STRONG_PEPTIDE_THRESHOLD

    out = df_clean[['__accession__', '__unique_pep__', '__strong_evidence__']].copy()
    for c in ["Gene", "gene", "Protein", "protein", "Description", "description"]:
        if c in df.columns:
            out[c] = df_clean[c].astype(str)
            break

    sample_name = sample_name or Path(path).stem
    out_path = outdir / f"{sample_name}_cleaned.csv"
    out.to_csv(out_path, index=False)
    print(f"Cleaned {path} -> {out_path} (rows kept: {len(out)})  | contaminant_col: {contaminant_col}")
    return out_path, id_col, pep_col, contaminant_col

# ---------------------------
# Merge presence/absence
# ---------------------------
def build_presence_absence(cleaned_paths, sample_files_map, control_files_map):
    dfs = {}
    for p in cleaned_paths:
        df = pd.read_csv(p)
        if '__accession__' not in df.columns:
            raise ValueError(f"Cleaned file {p} missing __accession__ column")
        dfs[p] = df.set_index('__accession__')
    all_accessions = sorted(set().union(*[set(d.index) for d in dfs.values()]))
    merged = pd.DataFrame(index=all_accessions)
    def colname_from_path(p):
        return Path(p).stem
    for p in sample_files_map:
        cname = colname_from_path(p)
        df = dfs[p]
        merged[cname + "_present"] = 0
        merged[cname + "_unique_pep"] = 0
        merged.loc[df.index, cname + "_present"] = 1
        merged.loc[df.index, cname + "_unique_pep"] = df['__unique_pep__']
    for p in control_files_map:
        cname = colname_from_path(p)
        df = dfs[p]
        merged[cname + "_present"] = 0
        merged[cname + "_unique_pep"] = 0
        merged.loc[df.index, cname + "_present"] = 1
        merged.loc[df.index, cname + "_unique_pep"] = df['__unique_pep__']
    for p, df in dfs.items():
        if 'Gene' in df.columns:
            merged['Gene'] = df['Gene']
            break
    merged.reset_index(inplace=True)
    merged = merged.rename(columns={'index': 'Accession'})
    return merged

# ---------------------------
# Robust UniProt: fetch gene names with local fallback and per-accession JSON lookup
# ---------------------------
UNIPROT_PER_ACCESSION_ENDPOINT = "https://rest.uniprot.org/uniprotkb/"

def extract_gene_from_cleaned_files(cleaned_paths):
    mapping = {}
    for p in cleaned_paths:
        try:
            df = pd.read_csv(p, dtype=str)
        except Exception:
            continue
        acc_col = None
        for c in ['__accession__','Accession','accession']:
            if c in df.columns:
                acc_col = c
                break
        if acc_col is None:
            for c in df.columns:
                if 'acc' in c.lower() or 'entry' in c.lower():
                    acc_col = c
                    break
        gene_col = None
        for c in ['Gene','gene','Gene names','Gene names  (primary )','Description','description','Protein names']:
            if c in df.columns:
                gene_col = c
                break
        if acc_col is None:
            continue
        for _, row in df.iterrows():
            acc = str(row.get(acc_col, "")).strip()
            if not acc:
                continue
            if gene_col and pd.notna(row.get(gene_col, "")) and str(row.get(gene_col)).strip() != "":
                if acc not in mapping or not mapping[acc]:
                    mapping[acc] = str(row.get(gene_col)).strip()
    return mapping

def lookup_uniprot_gene_per_accession(accession, timeout=30):
    url = UNIPROT_PER_ACCESSION_ENDPOINT + f"{accession}.json"
    attempt = 0
    while attempt < 4:
        attempt += 1
        try:
            resp = requests.get(url, timeout=timeout)
            if resp.status_code == 200:
                try:
                    data = resp.json()
                except Exception:
                    return None
                genes = data.get('genes', []) or []
                if genes:
                    for g in genes:
                        gname = None
                        gn = g.get('geneName') or g.get('gene_name') or {}
                        if isinstance(gn, dict):
                            gname = gn.get('value') or gn.get('primary') or gn.get('name')
                        if not gname:
                            syn = g.get('synonyms') or []
                            if syn and isinstance(syn, list):
                                if syn:
                                    gname = syn[0] if isinstance(syn[0], str) else (syn[0].get('value') if isinstance(syn[0], dict) else None)
                        if gname:
                            return str(gname)
                return None
            elif resp.status_code in (429, 503):
                time.sleep(1 + attempt * 1.5)
            elif resp.status_code == 404:
                return None
            else:
                time.sleep(0.5 + attempt * 0.5)
        except Exception:
            time.sleep(0.5 + attempt * 0.5)
    return None

def fetch_uniprot_gene_names_robust(accessions, cleaned_paths=None, sleep=0.4):
    mapping_local = {}
    if cleaned_paths:
        mapping_local = extract_gene_from_cleaned_files(cleaned_paths)
    rows = []
    for acc in accessions:
        gene = None
        if mapping_local and acc in mapping_local and mapping_local[acc]:
            gene = mapping_local[acc]
        else:
            gene = lookup_uniprot_gene_per_accession(acc)
            time.sleep(sleep)
        rows.append({"Accession": acc, "gene_names": gene})
    df = pd.DataFrame(rows).set_index("Accession")
    return df

# ---------------------------
# Fisher and BH (optional)
# ---------------------------
def run_fisher_and_adjust(merged_df, sample_cols_present, control_cols_present):
    a_counts = merged_df[sample_cols_present].sum(axis=1).astype(int)
    b_counts = merged_df[control_cols_present].sum(axis=1).astype(int)
    pvals = []
    odds = []
    for a, b in zip(a_counts, b_counts):
        table = np.array([[a, b], [3 - a, 3 - b]])
        try:
            odr, p = fisher_exact(table, alternative='greater')
        except TypeError:
            odr, p = fisher_exact(table)
        pvals.append(p)
        odds.append(odr)
    merged_df['sample_present_count'] = a_counts
    merged_df['control_present_count'] = b_counts
    merged_df['freq_sample'] = a_counts / 3.0
    merged_df['freq_control'] = b_counts / 3.0
    merged_df['freq_diff'] = merged_df['freq_sample'] - merged_df['freq_control']
    merged_df['fisher_pvalue'] = pvals
    merged_df['odds_ratio'] = odds
    _, p_adj, _, _ = multipletests(merged_df['fisher_pvalue'].values, method='fdr_bh')
    merged_df['fisher_p_adj'] = p_adj
    return merged_df

# ---------------------------
# Classify candidates (frequency-based)
# ---------------------------
def classify_candidates(df, sample_prefixes=None):
    """
    Classification rules with added behavior:
    - strict_interactor: same as before (>= STRICT_SAMPLE_REPLICATES in sample AND 0 in control)
    - relaxed_candidate: (original rule) OR (freq_diff >= RELAXED_FREQ_DIFF_THRESHOLD), but never includes strict interactors.
    """
    df = df.copy()
    present_cols = [c for c in df.columns if c.endswith("_present")]
    sample_cols = [c for c in present_cols if any(pref in c for pref in sample_prefixes)] if sample_prefixes else present_cols[:3]
    control_cols = [c for c in present_cols if not any(pref in c for pref in sample_prefixes)] if sample_prefixes else present_cols[3:]
    if 'sample_present_count' not in df.columns:
        df['sample_present_count'] = df[sample_cols].sum(axis=1).astype(int)
    if 'control_present_count' not in df.columns:
        df['control_present_count'] = df[control_cols].sum(axis=1).astype(int)
    df['freq_sample'] = df['sample_present_count'] / 3.0
    df['freq_control'] = df['control_present_count'] / 3.0
    df['freq_diff'] = df['freq_sample'] - df['freq_control']

    # strict: >= STRICT_SAMPLE_REPLICATES in sample AND 0 in control
    df['strict_interactor'] = ((df['sample_present_count'] >= STRICT_SAMPLE_REPLICATES) & (df['control_present_count'] == 0))

    # original relaxed rule
    relaxed_orig = ((df['sample_present_count'] >= STRICT_SAMPLE_REPLICATES) & (df['control_present_count'] <= RELAXED_CONTROL_MAX))

    # new: include proteins with frequency difference >= RELAXED_FREQ_DIFF_THRESHOLD
    relaxed_by_freq = df['freq_diff'] >= RELAXED_FREQ_DIFF_THRESHOLD

    # combine: either original relaxed OR by-freq, but exclude strict_interactor
    df['relaxed_candidate'] = ( (relaxed_orig | relaxed_by_freq) & (~df['strict_interactor']) )

    pep_cols = [c for c in df.columns if c.endswith("_unique_pep")]
    if sample_prefixes:
        sample_pep_cols = [c for c in pep_cols if any(pref in c for pref in sample_prefixes)]
    else:
        sample_pep_cols = pep_cols[:3] if len(pep_cols) >= 3 else pep_cols
    df['max_sample_unique_pep'] = df[sample_pep_cols].max(axis=1)
    df['low_confidence_1peptide'] = df['max_sample_unique_pep'] < STRONG_PEPTIDE_THRESHOLD

    return df

# ---------------------------
# (Optional) Plot helpers - unchanged
# ---------------------------
def plot_heatmap_binary(df, sample_cols_present, control_cols_present, top_n=50, outpath=None):
    top = df.sort_values('freq_diff', ascending=False).head(top_n)
    cols = sample_cols_present + control_cols_present
    hm = top.set_index('Accession')[cols]
    plt.figure(figsize=(10, max(5, 0.15 * len(hm))))
    sns.heatmap(hm, cmap="Blues", cbar=False, linewidths=0.5, linecolor='gray')
    plt.title("Presence(1)/Absence(0) - top {} proteins by freq_diff".format(top_n))
    plt.tight_layout()
    if outpath:
        plt.savefig(outpath, dpi=300)
    plt.close()

def plot_freqdiff_bar(df, top_n=30, outpath=None):
    top = df.sort_values('freq_diff', ascending=False).head(top_n)
    plt.figure(figsize=(10,6))
    sns.barplot(x='freq_diff', y='Accession', data=top)
    plt.xlabel("freq_sample - freq_control")
    plt.ylabel("Accession")
    plt.title("Top {} proteins by freq_diff".format(top_n))
    plt.tight_layout()
    if outpath:
        plt.savefig(outpath, dpi=300)
    plt.close()

def plot_peptide_hist(df, sample_files, control_files, outpath=None):
    pep_cols = [c for c in df.columns if c.endswith("_unique_pep")]
    sample_stems = [Path(x).stem for x in sample_files]
    control_stems = [Path(x).stem for x in control_files]
    sample_cols = [c for c in pep_cols if any(s in c for s in sample_stems)]
    control_cols = [c for c in pep_cols if any(s in c for s in control_stems)]
    if not sample_cols:
        sample_cols = pep_cols[:3]
        control_cols = pep_cols[3:]
    sample_vals = df[sample_cols].values.flatten()
    control_vals = df[control_cols].values.flatten()
    sample_vals = sample_vals[~np.isnan(sample_vals)]
    control_vals = control_vals[~np.isnan(control_vals)]
    plt.figure(figsize=(10,4))
    plt.subplot(1,2,1)
    sns.histplot(sample_vals, discrete=True)
    plt.title("Sample unique peptide counts (all reps)")
    plt.subplot(1,2,2)
    sns.histplot(control_vals, discrete=True)
    plt.title("Control unique peptide counts (all reps)")
    plt.tight_layout()
    if outpath:
        plt.savefig(outpath, dpi=300)
    plt.close()

def plot_volcano_like(df, outpath=None):
    df_plot = df.copy()
    df_plot['neglog10_padj'] = -np.log10(df_plot.get('fisher_p_adj', 1.0).replace(0, 1e-300))
    plt.figure(figsize=(8,6))
    sns.scatterplot(x='freq_diff', y='neglog10_padj', data=df_plot, hue='strict_interactor', legend='brief', palette={True:'red', False:'gray'})
    plt.axvline(0, color='black', linestyle='--', linewidth=0.7)
    plt.xlabel("freq_sample - freq_control")
    plt.ylabel("-log10(BH FDR)")
    plt.title("Presence/absence enrichment (volcano-like)")
    plt.tight_layout()
    if outpath:
        plt.savefig(outpath, dpi=300)
    plt.close()

# ---------------------------
# Main pipeline
# ---------------------------
def main(sample_files, control_files):
    cleaned_paths = []
    contaminant_columns_used = {}
    detected_id_cols = {}
    detected_pep_cols = {}
    print("Cleaning replicates (using explicit contaminant column when available)...")
    for p in sample_files + control_files:
        p = Path(p)
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p}")
        stem = p.stem
        cleaned_path, idcol, pepcol, contcol = clean_replicate(p, OUTDIR, sample_name=stem)
        cleaned_paths.append(str(cleaned_path))
        contaminant_columns_used[str(cleaned_path)] = contcol
        detected_id_cols[str(cleaned_path)] = idcol
        detected_pep_cols[str(cleaned_path)] = pepcol

    sample_cleaned = cleaned_paths[:len(sample_files)]
    control_cleaned = cleaned_paths[len(sample_files):]

    print("Building presence/absence matrix...")
    merged = build_presence_absence(cleaned_paths, sample_cleaned, control_cleaned)
    merged_path = OUTDIR / "presence_absence_matrix.csv"
    merged.to_csv(merged_path, index=False)
    print("Saved presence/absence matrix:", merged_path)

    merged_cols = list(merged.columns)
    sample_stems = [Path(p).stem for p in sample_files]
    control_stems = [Path(p).stem for p in control_files]

    sample_cols_present = [c for c in merged_cols if c.endswith("_present") and any(st in c for st in sample_stems)]
    control_cols_present = [c for c in merged_cols if c.endswith("_present") and any(st in c for st in control_stems)]

    if not sample_cols_present:
        sample_cols_present = [Path(p).stem + "_present" for p in sample_cleaned]
    if not control_cols_present:
        control_cols_present = [Path(p).stem + "_present" for p in control_cleaned]

    print("Using sample present columns:", sample_cols_present)
    print("Using control present columns:", control_cols_present)

    merged['sample_present_count'] = merged[sample_cols_present].sum(axis=1).astype(int)
    merged['control_present_count'] = merged[control_cols_present].sum(axis=1).astype(int)
    merged['freq_sample'] = merged['sample_present_count'] / 3.0
    merged['freq_control'] = merged['control_present_count'] / 3.0
    merged['freq_diff'] = merged['freq_sample'] - merged['freq_control']

    merged = classify_candidates(merged, sample_prefixes=[Path(p).stem for p in sample_files])

    if RUN_STATS:
        print("Running Fisher's exact tests and BH correction...")
        merged = run_fisher_and_adjust(merged, sample_cols_present, control_cols_present)
    else:
        merged['fisher_pvalue'] = np.nan
        merged['odds_ratio'] = np.nan
        merged['fisher_p_adj'] = np.nan

    # ---------------------------
    # Fetch gene names only and merge
    # ---------------------------
    print("Mapping gene names via UniProt (accessions -> gene names). This requires internet access.")
    accessions = list(merged['Accession'].astype(str).unique())
    accessions = list(dict.fromkeys(accessions))
    gene_df = None
    try:
        # use robust fetcher from this file
        gene_df = fetch_uniprot_gene_names_robust(accessions, cleaned_paths=cleaned_paths)
        gene_df = gene_df.reset_index().rename(columns={"index":"Accession"}) if isinstance(gene_df.index, pd.Index) else gene_df.reset_index()
        if "Accession" in gene_df.columns:
            merged = merged.merge(gene_df[['Accession','gene_names']], on='Accession', how='left')
        else:
            merged['gene_names'] = np.nan
    except Exception as e:
        print("UniProt gene mapping failed:", e)
        merged['gene_names'] = np.nan

    final_path = OUTDIR / "final_candidates_annotated.csv"
    merged.to_csv(final_path, index=False)
    print("Saved final annotated table:", final_path)

    strict_df = merged[merged['strict_interactor']].copy()
    relaxed_df = merged[merged['relaxed_candidate']].copy()
    strict_path = OUTDIR / "strict_interactors.csv"
    relaxed_path = OUTDIR / "relaxed_interactors.csv"
    strict_df.to_csv(strict_path, index=False)
    relaxed_df.to_csv(relaxed_path, index=False)
    print("Saved strict and relaxed lists:", strict_path, relaxed_path)

    if RUN_PLOTS:
        print("Generating plots...")
        plot_heatmap_binary(merged, sample_cols_present, control_cols_present, top_n=50, outpath=OUTDIR / "heatmap_top50.png")
        plot_freqdiff_bar(merged, top_n=30, outpath=OUTDIR / "freqdiff_top30.png")
        plot_peptide_hist(merged, sample_files, control_files, outpath=OUTDIR / "peptide_histograms.png")
        if RUN_STATS:
            plot_volcano_like(merged, outpath=OUTDIR / "presence_volcano_like.png")
        try:
            sample_union = set(merged.loc[merged[sample_cols_present].sum(axis=1) > 0, 'Accession'])
            control_union = set(merged.loc[merged[control_cols_present].sum(axis=1) > 0, 'Accession'])
            if HAVE_VENN:
                plt.figure(figsize=(6,6))
                venn2([sample_union, control_union], set_labels=("Sample (18_GFP)", "Control (Mito_GFP)"))
                plt.title("Sample vs Control union overlap")
                plt.tight_layout()
                plt.savefig(OUTDIR / "venn_sample_control.png", dpi=300)
                plt.close()
            else:
                with open(OUTDIR / "venn_numbers.txt", "w") as fh:
                    fh.write(f"Sample union count: {len(sample_union)}\n")
                    fh.write(f"Control union count: {len(control_union)}\n")
                    fh.write(f"Overlap count: {len(sample_union & control_union)}\n")
        except Exception as e:
            print("Venn plot failed or skipped:", e)

    print("Contaminant columns used (per cleaned file):")
    for k,v in contaminant_columns_used.items():
        print(" ", k, "->", v)
    print("Done. All outputs in:", OUTDIR.resolve())

if __name__ == "__main__":
    main(sample_files, control_files)
