#!/usr/bin/env python3
"""
final_plots_top10_with_list_UPDATED.py

- Auto-detects Excel header row
- Generates a ranked protein list
- Optionally generates plots (toggleable)
- Ranked list includes TAIR_best_hit column (if present)
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ================= CONFIG =================
INPUT = "mitochondrial_interactors FINAL.xlsx"
SHEET_NAME = 0
OUTDIR = Path("Neha di Mass spec")
OUTDIR.mkdir(exist_ok=True)

TOP_N = 10
FIG_DPI = 300
FONT = "DejaVu Sans"

MAKE_PLOTS = True   # <<< MASTER SWITCH
# ==========================================

# ---------- COLUMN CANDIDATES ----------
EFFECT_CANDIDATES = ["freq_diff", "log2fc", "logfc", "effect", "delta"]
ADJP_CANDIDATES   = ["padj", "adj_p", "fdr", "qvalue"]
P_CANDIDATES      = ["p", "pval", "p_value"]
ID_CANDIDATES     = ["accession", "protein", "id", "gene"]
NAME_CANDIDATES   = ["protein_name", "gene_name", "description", "name"]
TAIR_CANDIDATES   = ["tair_best_hit", "tair_hit"]

ALL_CANDIDATES = (
    EFFECT_CANDIDATES +
    ADJP_CANDIDATES +
    P_CANDIDATES +
    ID_CANDIDATES +
    NAME_CANDIDATES +
    TAIR_CANDIDATES
)

sns.set(style="whitegrid", rc={"font.family": FONT})


# ---------- HELPERS ----------
def detect_header_row(excel_path, sheet=0, max_rows=10):
    """Detect which row contains the real header."""
    tmp = pd.read_excel(
        excel_path,
        sheet_name=sheet,
        header=None,
        nrows=max_rows,
        engine="openpyxl"
    )

    best_row = 0
    best_score = -1

    for i in range(len(tmp)):
        row = tmp.iloc[i].astype(str).str.lower()
        score = sum(
            any(key in cell for key in ALL_CANDIDATES)
            for cell in row
        )
        if score > best_score:
            best_score = score
            best_row = i

    return best_row


def find_col(columns, candidates):
    for cand in candidates:
        for c in columns:
            if cand.lower() == str(c).lower():
                return c
    for cand in candidates:
        for c in columns:
            if cand.lower() in str(c).lower():
                return c
    return None


def save_fig(fig, name):
    fig.savefig(OUTDIR / f"{name}.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


# ================= MAIN =================
def main():

    # ---- header detection ----
    header_row = detect_header_row(INPUT, SHEET_NAME)
    print("Using detected header row:", header_row)

    df = pd.read_excel(
        INPUT,
        sheet_name=SHEET_NAME,
        header=header_row,
        engine="openpyxl"
    )

    print("Loaded dataframe shape:", df.shape)
    print("Columns:")
    for c in df.columns:
        print(" -", c)

    # ---- detect columns ----
    effect_col = find_col(df.columns, EFFECT_CANDIDATES)
    adjp_col   = find_col(df.columns, ADJP_CANDIDATES)
    p_col      = find_col(df.columns, P_CANDIDATES)
    id_col     = find_col(df.columns, ID_CANDIDATES)
    name_col   = find_col(df.columns, NAME_CANDIDATES)
    tair_col   = find_col(df.columns, TAIR_CANDIDATES)

    print("\nDetected mappings:")
    print(" effect_col ->", effect_col)
    print(" adjp_col   ->", adjp_col)
    print(" p_col      ->", p_col)
    print(" id_col     ->", id_col)
    print(" name_col   ->", name_col)
    print(" tair_col   ->", tair_col)

    # ---- validation (STRICT) ----
    if effect_col is None:
        raise SystemExit("ERROR: No effect column detected — cannot rank.")
    if id_col is None:
        raise SystemExit("ERROR: No identifier column detected — cannot rank.")

    if adjp_col is None and p_col is not None:
        adjp_col = p_col
        print("NOTE: Using raw p-values as adj_p.")

    if adjp_col is None:
        raise SystemExit("ERROR: No p-value or adjusted p-value column detected.")

    # ---- numeric conversion ----
    df["_effect"] = pd.to_numeric(df[effect_col], errors="coerce")
    df["_adj_p"]  = pd.to_numeric(df[adjp_col], errors="coerce")

    df["_minuslog10p"] = df["_adj_p"].apply(
        lambda x: -np.log10(x) if pd.notna(x) and x > 0 else 0.0
    )
    df["_abs_effect"] = df["_effect"].abs().fillna(0.0)

    # ---- ranking ----
    df["rank_score"] = df["_abs_effect"] * (df["_minuslog10p"] + 1e-12)

    # ---- identifiers ----
    df["_id_raw"] = df[id_col].astype(str)
    df["_id"] = df["_id_raw"].str.strip().str.upper()

    # ---- labels ----
    if name_col:
        df["_label"] = df.apply(
            lambda r: f"{r['_id_raw']} — {r[name_col]}"
            if pd.notna(r[name_col]) else r["_id_raw"],
            axis=1
        )
    else:
        df["_label"] = df["_id_raw"]

    # ---- ranked table ----
    ranked = df.sort_values("rank_score", ascending=False).reset_index(drop=True)

    out_cols = [
        "_id_raw", "_id", "_label",
        "_effect", "_adj_p", "_abs_effect", "_minuslog10p", "rank_score"
    ]

    if tair_col and tair_col in ranked.columns:
        out_cols.insert(3, tair_col)
    else:
        ranked["TAIR_best_hit"] = pd.NA
        out_cols.insert(3, "TAIR_best_hit")

    ranked_out = ranked[out_cols].rename(columns={
        "_id_raw": "identifier_raw",
        "_id": "identifier",
        "_label": "label",
        "_effect": "effect",
        "_adj_p": "adj_p",
        "_abs_effect": "abs_effect",
        "_minuslog10p": "minuslog10p"
    })

    ranked_out.to_csv(OUTDIR / "mito_interactors_tair_hits.csv", index=False)
    print("\nSaved ranked list →", OUTDIR / "mito_interactors_tair_hits.csv")

    # ---- plotting (optional) ----
    if not MAKE_PLOTS:
        print("MAKE_PLOTS=False → skipping all plots.")
        return

    topN = ranked.head(TOP_N)

    fig, ax = plt.subplots(figsize=(8, 0.6 * len(topN) + 1))
    ax.barh(topN["_label"], topN["rank_score"])
    ax.invert_yaxis()
    ax.set_xlabel("Ranking score")
    save_fig(fig, "barplot_topN")

    print("Plots generated.")


# ============================================
if __name__ == "__main__":
    main()
