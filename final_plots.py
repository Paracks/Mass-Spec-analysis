#!/usr/bin/env python3
"""
final_plots_top10_with_list.py (labels with protein names)

- Same as your script, but adds protein names (when available) alongside accession numbers
  in barplot, volcano labels, and heatmap y-ticks.
All other behavior preserved.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler

# ---------- CONFIG ----------
INPUT = "mitochondrial_interactors FINAL.xlsx"   # path to the Excel file you uploaded
SHEET_NAME = 0    # use sheet index or name; 0 uses the first sheet
OUTDIR = Path("ms_plots_top10_plus_list")
OUTDIR.mkdir(exist_ok=True)
TOP_N = 10                    # top N by rank to always include
FIG_DPI = 300
FONT = "DejaVu Sans"          # publication-friendly font
# ----------------------------

# ---------- USER-PROVIDED ACCESSION LIST ----------
REQUESTED_ACCESSIONS = [
"A0A178UPX1", "A0A1P8BEB9", "A0A5S9XX21", "P82873", "A0A8S2ACP2", "A0A654EQ41", "A0A8S1ZU39", "A0A8S2A788",
]
REQUESTED_ACCESSIONS = set(a.strip().upper() for a in REQUESTED_ACCESSIONS if a and str(a).strip())

# set aesthetics
sns.set(style="whitegrid", rc={
    "font.family": FONT,
    "axes.grid": True,
    "grid.color": "#e9e9e9",
    "axes.facecolor": "white",
    "figure.facecolor": "white"
})

# candidate lists (unchanged)
EFFECT_CANDIDATES = [
    "freq_diff", "frequency_difference", "frequency diff", "freq.diff", "difference",
    "freqDiff", "freq_difference", "log2fc", "logfc", "log2_fold_change", "fold_change",
    "fc", "effect", "delta"
]
ADJP_CANDIDATES = [
    "fisher_p_adj", "fisher_pvalue_adj", "p_adj", "adj_p", "adj.p", "padj", "p.adjust",
    "fdr", "q_value", "qvalue", "adj.p.value", "adjp"
]
P_CANDIDATES = [
    "fisher_pvalue", "p_value", "p.value", "pval", "p"
]
ID_CANDIDATES = [
    "accession", "Accession", "protein", "Protein", "protein_name", "Protein Name", "id", "ID", "gene", "Gene"
]

# NEW: candidate name/description columns to find protein names
NAME_CANDIDATES = [
    "protein_name", "protein name", "name", "gene", "gene_name", "gene name",
    "description", "entry name", "entry_name", "protein description", "recommended name"
]

def detect_header_row(excel_path, sheet=0, max_rows=10):
    tmp = pd.read_excel(excel_path, sheet_name=sheet, header=None, nrows=max_rows, engine='openpyxl')
    candidates = set([c.lower() for c in EFFECT_CANDIDATES + ADJP_CANDIDATES + P_CANDIDATES + ID_CANDIDATES + NAME_CANDIDATES])
    best_row = 0
    best_matches = 0
    for i in range(len(tmp)):
        row = tmp.iloc[i].astype(str).str.lower().tolist()
        matches = sum(1 for cell in row if any(cand in cell for cand in candidates))
        nonempty = sum(1 for cell in row if cell.strip() != "nan" and cell.strip() != "")
        score = matches + 0.25 * nonempty
        if score > best_matches:
            best_matches = score
            best_row = i
    return best_row

def find_col(df_cols, candidates):
    cols = list(df_cols)
    low2orig = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in low2orig:
            return low2orig[cand.lower()]
    for cand in candidates:
        for c in cols:
            if cand.lower() in c.lower():
                return c
    return None

def save_fig(fig, name):
    path = OUTDIR / (name + ".png")
    fig.savefig(path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", path)

def _place_nonoverlapping_texts(ax, x, y, labels, fontsize=7, xoffset=0.0, ystep=0.03, max_iter=50, color="k"):
    texts = []
    renderer = ax.figure.canvas.get_renderer()
    placed_bboxes = []
    yrange = ax.get_ylim()
    y_span = abs(yrange[1] - yrange[0]) if yrange[1] != yrange[0] else 1.0
    ystep_data = ystep * y_span
    for xi, yi, lab in zip(x, y, labels):
        yi_loc = yi + xoffset
        yi_try = yi_loc
        txt = ax.text(xi, yi_try, lab, fontsize=fontsize, ha='center', va='bottom', color=color, zorder=6)
        iter_count = 0
        while iter_count < max_iter:
            txt_bbox = txt.get_window_extent(renderer=renderer)
            overlap = any(txt_bbox.overlaps(pb) for pb in placed_bboxes)
            if not overlap:
                placed_bboxes.append(txt_bbox)
                texts.append(txt)
                break
            yi_try = yi_try + ystep_data
            txt.set_position((xi, yi_try))
            iter_count += 1
        if iter_count >= max_iter:
            txt.set_alpha(0.7)
            placed_bboxes.append(txt.get_window_extent(renderer=renderer))
            texts.append(txt)
    return texts

def main():
    # detect header row
    print("Detecting header row in", INPUT)
    try:
        header_row = detect_header_row(INPUT, sheet=SHEET_NAME, max_rows=10)
    except Exception as e:
        print("Header detection failed:", e)
        header_row = 0

    print("Using header row index:", header_row)
    try:
        df = pd.read_excel(INPUT, sheet_name=SHEET_NAME, header=header_row, engine='openpyxl')
    except Exception as e:
        print("Failed to read Excel with header row. Attempting fallback read (header=None). Error:", e)
        df = pd.read_excel(INPUT, sheet_name=SHEET_NAME, header=None, engine='openpyxl')

    print("Loaded data. Shape:", df.shape)
    print("Columns detected:")
    for c in df.columns:
        print(" -", c)

    # detect columns
    effect_col = find_col(df.columns, EFFECT_CANDIDATES)
    adjp_col = find_col(df.columns, ADJP_CANDIDATES)
    p_col = find_col(df.columns, P_CANDIDATES)
    id_col = find_col(df.columns, ID_CANDIDATES)
    name_col = find_col(df.columns, NAME_CANDIDATES)  # NEW: detect name/description column

    print("\nAuto-detected columns:")
    print(" effect_col ->", effect_col)
    print(" adjp_col   ->", adjp_col)
    print(" p_col      ->", p_col)
    print(" id_col     ->", id_col)
    print(" name_col   ->", name_col)

    if adjp_col is None and p_col is not None:
        print("No adjusted-p detected; falling back to raw p column:", p_col)
        adjp_col = p_col

    # compute effect if missing (heuristic)
    if effect_col is None:
        freq_like = [c for c in df.columns if "freq" in str(c).lower() or "present" in str(c).lower() or "count" in str(c).lower()]
        print("No direct effect column. Looking for frequency-like columns:", freq_like)
        if len(freq_like) >= 2:
            a, b = freq_like[0], freq_like[1]
            print(f"Computing effect as {a} - {b}")
            df["_effect_computed"] = pd.to_numeric(df[a], errors='coerce') - pd.to_numeric(df[b], errors='coerce')
            effect_col = "_effect_computed"

    # numeric conversion and ranking
    df["_effect"] = pd.to_numeric(df[effect_col], errors='coerce') if (effect_col is not None and effect_col in df.columns) else np.nan
    df["_adj_p"] = pd.to_numeric(df[adjp_col], errors='coerce') if (adjp_col is not None and adjp_col in df.columns) else np.nan
    df["_minuslog10p"] = df["_adj_p"].apply(lambda x: -np.log10(x) if pd.notna(x) and x>0 else 0.0)
    df["_abs_effect"] = df["_effect"].abs().fillna(0.0)
    df["rank_score"] = df["_abs_effect"] * (df["_minuslog10p"] + 1e-12)

    # identifier fallback
    if id_col is not None and id_col in df.columns:
        df["_id_raw"] = df[id_col].astype(str)
    else:
        firstcol = df.columns[0]
        df["_id_raw"] = df[firstcol].astype(str)

    # normalized identifier column for matching (uppercase, strip)
    df["_id"] = df["_id_raw"].astype(str).str.strip().str.upper()

    # NEW: build readable label combining accession and name (if name_col found)
    if name_col is not None and name_col in df.columns:
        # coerce to string and trim long names
        def make_label(acc_raw, name_val, max_len=60):
            name_str = "" if pd.isna(name_val) else str(name_val).strip()
            label = f"{str(acc_raw).strip()} — {name_str}" if name_str else str(acc_raw).strip()
            if len(label) > max_len:
                return label[:max_len-3] + "..."
            return label
        df["_label"] = df.apply(lambda r: make_label(r["_id_raw"], r[name_col]), axis=1)
    else:
        df["_label"] = df["_id_raw"].astype(str)

    # normalized uppercase label for plotting lookups (we keep _id for matching)
    # (no change to matching logic — REQUESTED_ACCESSIONS still matched to _id)
    # find requested accessions present in the dataset
    present_requested = sorted(list(REQUESTED_ACCESSIONS.intersection(set(df["_id"].dropna().unique()))))
    missing_requested = sorted(list(REQUESTED_ACCESSIONS.difference(set(df["_id"].dropna().unique()))))

    # --- NEW: ensure these two requested accessions (if present) are placed at the end ---
    # the special accessions to push to bottom (uppercase)
    specials_to_bottom = ["A0A8S1ZU39", "A0A8S2A788"]
    # preserve original order among other requested accessions
    reordered_requested = [x for x in present_requested if x not in specials_to_bottom]
    # append the specials in the specified order only if they are present
    for s in specials_to_bottom:
        if s in present_requested:
            reordered_requested.append(s)
    # replace present_requested with reordered list for downstream plotting/order logic
    present_requested = reordered_requested
    # -------------------------------------------------------------------------------

    print(f"Requested accessions present: {len(present_requested)}; missing: {len(missing_requested)}")

    # save ranked CSV
    ranked = df.sort_values("rank_score", ascending=False).reset_index(drop=True)
    ranked_out = ranked[["_id_raw", "_id", "_label", "_effect", "_adj_p", "_abs_effect", "_minuslog10p", "rank_score"]].rename(
        columns={"_id_raw":"identifier_raw", "_id":"identifier", "_label":"label", "_effect":"effect", "_adj_p":"adj_p", "_abs_effect":"abs_effect", "_minuslog10p":"minuslog10p"})
    ranked_out.to_csv(OUTDIR/"ranked_proteins.csv", index=False)
    print("Saved ranked table to", OUTDIR/"ranked_proteins.csv")

    # determine top-N by rank
    top_by_rank = ranked.head(TOP_N)
    top_by_rank_ids = list(top_by_rank["_id"])
    top_by_rank_labels = list(top_by_rank["_label"])

    # build top_set = union(top_by_rank, requested present)
    top_set_ids = list(dict.fromkeys(top_by_rank_ids + present_requested))  # preserve order: top-ranked first, then requested
    print(f"Final number of proteins included in per-protein plots: {len(top_set_ids)} (TOP_N={TOP_N} + requested_present={len(present_requested)})")
    # --- FORCE A0A1P8BEB9 (PAM18-3) TO TOP OF PLOTS ---
    FORCE_TOP_ID = "A0A1P8BEB9"
    if FORCE_TOP_ID in top_set_ids:
        top_set_ids = [FORCE_TOP_ID] + [x for x in top_set_ids if x != FORCE_TOP_ID]

# ----------------------------------------------------------------


    # subset dataframe for plotting (preserve ranking order)
    subset_df = ranked[ranked["_id"].isin(top_set_ids)].copy()
    subset_df["_order"] = subset_df["_id"].apply(lambda x: top_set_ids.index(x) if x in top_set_ids else len(top_set_ids))
    subset_df = subset_df.sort_values("_order").reset_index(drop=True)

    # use labels for plotting display but preserve accession-based ids for logic
    subset_df["_display_label"] = subset_df["_label"]

    # -------------------------
    # Barplot for top_set (use _display_label)
    # -------------------------
    if not subset_df.empty and (subset_df["rank_score"] > 0).any():
        fig, ax = plt.subplots(figsize=(8, 0.7*len(subset_df) + 1.5))
        colors = []
        for xid in subset_df["_id"]:
            if xid in top_by_rank_ids:
                colors.append("#2b8cbe")   # top10 color (blue)
            elif xid in present_requested:
                colors.append("#fdae61")   # requested accession color (orange)
            else:
                colors.append("#9aa5b1")
        ax.barh(subset_df["_display_label"], subset_df["rank_score"], color=colors, edgecolor="k", linewidth=0.4)
        ax.invert_yaxis()
        ax.set_xlabel("Ranking score (|effect| × -log10(adj p))", fontsize=9)
        from matplotlib.patches import Patch
        legend_handles = [
            Patch(facecolor="#2b8cbe", edgecolor="k", label=f"Top {TOP_N} by rank"),
            Patch(facecolor="#fdae61", edgecolor="k", label="Requested accession (present)"),
        ]
        for spine in ["top","right"]:
            ax.spines[spine].set_visible(False)
        plt.tight_layout()
        save_fig(fig, "barplot_topN_plus_requested")
    else:
        print("Skipping barplot: no non-zero rank scores found for subset.")

    # -------------------------
    # Volcano: highlight top10 and requested with better label placement (use labels)
    # -------------------------
    if df["_effect"].notna().any() and df["_adj_p"].notna().any():
        fig, ax = plt.subplots(figsize=(8,6))
        ax.scatter(df["_effect"], df["_minuslog10p"], s=18, alpha=0.35, color="#c7c9cc", linewidth=0)

        # plot requested (but not in top10) with orange
        req_only = [x for x in present_requested if x not in top_by_rank_ids]
        req_df = df[df["_id"].isin(req_only)]
        if not req_df.empty:
            ax.scatter(req_df["_effect"], req_df["_minuslog10p"], s=50, color="#fdae61", edgecolor="k", linewidth=0.5, zorder=4)

        # plot top10 with blue
        top_df = df[df["_id"].isin(top_by_rank_ids)]
        if not top_df.empty:
            ax.scatter(top_df["_effect"], top_df["_minuslog10p"], s=70, color="#2b8cbe", edgecolor="k", linewidth=0.6, zorder=5)

        # label top10 (larger) and requested (smaller) with collision avoidance
        top_x = top_df["_effect"].tolist() if not top_df.empty else []
        top_y = top_df["_minuslog10p"].tolist() if not top_df.empty else []
        top_labels = top_df.apply(lambda r: r["_label"], axis=1).tolist() if not top_df.empty else []

        req_x = req_df["_effect"].tolist() if not req_df.empty else []
        req_y = req_df["_minuslog10p"].tolist() if not req_df.empty else []
        req_labels = req_df.apply(lambda r: r["_label"], axis=1).tolist() if not req_df.empty else []

        if top_labels:
            _place_nonoverlapping_texts(ax, top_x, top_y, top_labels, fontsize=7, ystep=0.03, color="#08306b")
        if req_labels:
            _place_nonoverlapping_texts(ax, req_x, req_y, req_labels, fontsize=6, ystep=0.025, color="#7f3b08")

        ax.axhline(-np.log10(0.05), linestyle="--", linewidth=0.8, color="#444444")
        ax.set_xlabel("Effect (frequency difference or computed)", fontsize=9)
        ax.set_ylabel("-log10(adjusted p-value)", fontsize=9)
        ax.set_title("Volcano plot — top10 & requested accessions highlighted", fontsize=10, fontweight="bold")
        ax.tick_params(labelsize=8)
        for spine in ["top","right"]:
            ax.spines[spine].set_visible(False)
        plt.tight_layout()
        save_fig(fig, "volcano_top10_and_requested")
    else:
        print("Skipping volcano: effect and/or p-values missing.")

    # -------------------------
    # Feature heatmap for subset (top + requested) — reduced features (use labels)
    # -------------------------
    numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c]) and c not in ["_effect","_adj_p","_minuslog10p","_abs_effect","rank_score"]]
    if not numeric_cols:
        numeric_cols = [c for c in df.columns if any(k in str(c).lower() for k in ["count","present","unique","odds","score","intensity","abundance"]) and pd.api.types.is_numeric_dtype(pd.to_numeric(df[c], errors='coerce'))]

    heatmap_cols = []
    if "_effect" in df.columns:
        heatmap_cols.append("_effect")
    priority_names = ["sample_present_count", "present_count", "max_sample_unique_pep", "odds_ratio", "unique_pep", "unique_peptides", "count"]
    extra_col = None
    for pn in priority_names:
        for c in numeric_cols:
            if pn in str(c).lower():
                extra_col = c
                break
        if extra_col:
            break
    if extra_col is None:
        for c in numeric_cols:
            if c != "_effect":
                extra_col = c
                break
    if extra_col is not None:
        heatmap_cols.append(extra_col)

    heatmap_cols = [c for c in heatmap_cols if c in df.columns]
    if len(heatmap_cols) == 0:
        print("Skipping feature heatmap: no relevant numeric columns found.")
    else:
        # Ensure subset_df is ordered by rank_score desc for heatmap display
        # If rank_score is all zero (uninformative), fallback to ordering by abs(effect)
        if subset_df["rank_score"].notna().any() and (subset_df["rank_score"].abs().sum() > 0):
            subset_for_heat_df = subset_df.sort_values("rank_score", ascending=False).copy()
        else:
            subset_for_heat_df = subset_df.sort_values("_abs_effect", ascending=False).copy()

        subset_for_heat = subset_for_heat_df[heatmap_cols].replace([np.inf, -np.inf], np.nan).fillna(0)
        if subset_for_heat.shape[0] >= 1 and subset_for_heat.shape[1] >= 1:
            if subset_for_heat.shape[1] == 1:
                colvals = subset_for_heat.iloc[:,0].astype(float).values
                if np.nanstd(colvals) == 0:
                    z = np.zeros_like(colvals).reshape(-1,1)
                else:
                    z = ((colvals - np.nanmean(colvals)) / (np.nanstd(colvals))).reshape(-1,1)
            else:
                scaler = StandardScaler()
                z = scaler.fit_transform(subset_for_heat.values)

            display_labels = subset_for_heat_df["_display_label"].tolist()

            fig, ax = plt.subplots(figsize=(max(4, 0.8*subset_for_heat.shape[1]), max(2.5, 0.6*subset_for_heat.shape[0])))
            sns.heatmap(z, cmap="vlag", yticklabels=display_labels, xticklabels=[c.replace("_"," ") for c in heatmap_cols],
                        cbar_kws={"label":"z-score"}, ax=ax)
            ax.set_title(f"Feature heatmap (top{TOP_N} + requested) — showing {len(heatmap_cols)} feature(s)", fontsize=10, fontweight="bold")
            ax.tick_params(axis='y', labelsize=7)
            ax.tick_params(axis='x', labelsize=8, rotation=45)
            plt.tight_layout()
            save_fig(fig, "feature_heatmap_topN_plus_requested_reduced")
        else:
            print("Not enough rows/columns for reduced feature heatmap.")

    # -------------------------
    # QC histograms (entire dataset)
    # -------------------------
    if df["_minuslog10p"].notna().any() and (df["_minuslog10p"]>0).any():
        fig, ax = plt.subplots(figsize=(6,3))
        ax.hist(df["_minuslog10p"].replace([np.inf, -np.inf], np.nan).dropna(), bins=40, edgecolor="k", linewidth=0.25)
        ax.set_xlabel("-log10(adj p)", fontsize=9)
        ax.set_title("Significance distribution", fontsize=9)
        ax.tick_params(labelsize=8)
        for spine in ["top","right"]:
            ax.spines[spine].set_visible(False)
        plt.tight_layout()
        save_fig(fig, "hist_minuslog10p")
    else:
        print("Skipping p-value histogram: no p-values available.")

    if df["_effect"].notna().any():
        fig, ax = plt.subplots(figsize=(6,3))
        ax.hist(df["_effect"].dropna(), bins=40, edgecolor="k", linewidth=0.25)
        ax.set_xlabel("Effect", fontsize=9)
        ax.set_title("Effect size distribution", fontsize=9)
        ax.tick_params(labelsize=8)
        for spine in ["top","right"]:
            ax.spines[spine].set_visible(False)
        plt.tight_layout()
        save_fig(fig, "hist_effect")
    else:
        print("Skipping effect histogram: no effect values found.")

    # write diagnostics
    with open(OUTDIR/"diagnostics.txt", "w") as fh:
        fh.write("Input file: " + str(INPUT) + "\n")
        fh.write("Used header_row index: " + str(header_row) + "\n\n")
        fh.write("Columns detected:\n")
        for c in df.columns:
            fh.write(" - " + str(c) + "\n")
        fh.write("\nAuto-mappings:\n")
        fh.write(f" effect_col -> {effect_col}\n")
        fh.write(f" adjp_col   -> {adjp_col}\n")
        fh.write(f" p_col      -> {p_col}\n")
        fh.write(f" id_col     -> {id_col}\n")
        fh.write(f" name_col   -> {name_col}\n")
        fh.write(f"\nRequested accessions (total): {len(REQUESTED_ACCESSIONS)}\n")
        fh.write(f"Requested accessions present in data: {len(present_requested)}\n")
        fh.write("Present (sample):\n")
        for a in present_requested[:200]:
            fh.write(" - " + a + "\n")
        fh.write("\nMissing (sample):\n")
        for a in missing_requested[:200]:
            fh.write(" - " + a + "\n")
        fh.write("\nSaved files:\n")
        for f in sorted(OUTDIR.iterdir()):
            fh.write(" - " + str(f) + "\n")

    print("Diagnostics saved to", OUTDIR/"diagnostics.txt")
    print("All done. Outputs are in", OUTDIR)

if __name__ == "__main__":
    main()
