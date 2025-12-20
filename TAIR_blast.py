#!/usr/bin/env python3
"""
uniprot_to_tair_blast_local.py

Reads a CSV with an accession column, fetches sequence from UniProt REST,
runs local blastp against TAIR10_proteins_db (built with makeblastdb),
adds column 'TAIR_best_hit' to CSV and saves output.

Requirements:
- Python 3.8+
- pip install pandas requests
- NCBI BLAST+ installed and 'blastp' & 'makeblastdb' available on PATH
- TAIR10_proteins_db created (see README above)
"""

import os
import subprocess
import tempfile
import time
import random
import pandas as pd
import requests

# ---------- USER CONFIG ----------
INPUT_CSV = "ranked_proteins.csv"
OUTPUT_CSV = "Mito_blast_ranked_proteins.csv"
ACCESSION_COLUMN_CANDIDATES = ["accession", "Accession", "Entry", "Acc", "UniProt", "identifier"]
TAIR_BLAST_DB = "TAIR10_proteins_db"   # the prefix you used with makeblastdb
UNIPROT_REST_TIMEOUT = 30
# BLAST parameters
BLASTP_MAX_TARGET_SEQS = 1
BLASTP_OUTFMT = "6 qseqid sseqid pident evalue bitscore stitle"
# polite delay between UniProt requests
MIN_DELAY = 0.5
MAX_DELAY = 1.5
# ---------- END CONFIG ----------


def find_accession_column(df):
    for c in ACCESSION_COLUMN_CANDIDATES:
        if c in df.columns:
            return c
    for c in df.columns:
        if "access" in c.lower():
            return c
    raise ValueError("Could not find accession column. Candidates tried: " + ", ".join(ACCESSION_COLUMN_CANDIDATES))


def fetch_uniprot_sequence(accession):
    """
    Query UniProt REST for fasta of accession.
    Returns sequence string or None.
    """
    urls = [
        f"https://rest.uniprot.org/uniprotkb/{accession}.fasta",
        f"https://rest.uniprot.org/uniprotkb/{accession}?format=fasta",
        f"https://www.uniprot.org/uniprot/{accession}.fasta"
    ]
    headers = {"User-Agent": "python-requests/automated-script"}
    for url in urls:
        try:
            r = requests.get(url, headers=headers, timeout=UNIPROT_REST_TIMEOUT)
            if r.status_code == 200 and r.text.startswith(">"):
                lines = [ln.strip() for ln in r.text.splitlines() if ln and not ln.startswith(">")]
                seq = "".join(lines)
                if seq:
                    return seq
        except requests.RequestException:
            pass
    return None


def run_blastp_get_best_hit(query_fasta_path, db_prefix):
    """
    Runs blastp with max_target_seqs=1 and returns the parsed best-hit line
    (tab-separated per BLAST outfmt fields). Returns None on failure / no hits.
    """
    cmd = [
        "blastp",
        "-query", query_fasta_path,
        "-db", db_prefix,
        "-max_target_seqs", str(BLASTP_MAX_TARGET_SEQS),
        "-outfmt", BLASTP_OUTFMT,
        "-evalue", "1e-5",
        "-num_threads", "1"
    ]
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=120)
    except Exception as e:
        print("blastp invocation failed:", e)
        return None

    out = cp.stdout.strip()
    if not out:
        return None
    # outfmt '6' returns one-line per hit; since max_target_seqs=1 there will be first line
    first_line = out.splitlines()[0]
    return first_line


def main():
    if not os.path.exists(INPUT_CSV):
        raise SystemExit(f"Input CSV not found: {INPUT_CSV}")
    # check BLAST DB exists (one of the index files)
    db_check = db_has_files = any(os.path.exists(TAIR_BLAST_DB + ext) for ext in [".phr", ".pin", ".psq", ".pal", ".psd", ".pnd"])
    if not db_check:
        # still allow proceeding but warn user
        print(f"WARNING: BLAST DB prefix '{TAIR_BLAST_DB}' not found (expected .pin/.psq etc). Make sure makeblastdb was run.")
    df = pd.read_csv(INPUT_CSV)
    try:
        accession_col = find_accession_column(df)
    except ValueError as e:
        raise SystemExit(str(e))

    df["TAIR_best_hit"] = pd.NA

    # temp dir to hold query FASTA files
    with tempfile.TemporaryDirectory() as tmpdir:
        for idx, row in df.iterrows():
            acc = str(row[accession_col]).strip()
            if not acc or acc.lower() in ["nan", "none"]:
                continue
            print(f"[{idx+1}/{len(df)}] Processing accession: {acc}")
            seq = fetch_uniprot_sequence(acc)
            if not seq:
                print(f"  -> Could not fetch UniProt sequence for {acc}; leaving blank.")
                df.at[idx, "TAIR_best_hit"] = pd.NA
                time.sleep(random.uniform(MIN_DELAY, MAX_DELAY))
                continue
            # write a small FASTA for this query
            query_path = os.path.join(tmpdir, f"{acc}.fasta")
            with open(query_path, "w") as fh:
                fh.write(f">{acc}\n")
                # wrap sequence lines at 80 chars (BLAST accepts long single line too)
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i+80] + "\n")

            # run blastp
            best_line = run_blastp_get_best_hit(query_path, TAIR_BLAST_DB)
            if best_line:
                # best_line contains qseqid sseqid pident evalue bitscore stitle
                df.at[idx, "TAIR_best_hit"] = best_line
                print("  -> Best hit:", best_line)
            else:
                df.at[idx, "TAIR_best_hit"] = pd.NA
                print("  -> No hit / parsing failed.")

            # small polite delay
            time.sleep(random.uniform(MIN_DELAY, MAX_DELAY))

    # save results
    df.to_csv(OUTPUT_CSV, index=False)
    print("Done. Output written to:", OUTPUT_CSV)


if __name__ == "__main__":
    main()
