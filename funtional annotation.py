#!/usr/bin/env python3
"""
functional_annotation_uniprot.py

Functional annotation of protein lists using UniProt REST API.

Produces:
 - annotations_full.csv   (one row per input ID with parsed annotation fields)
 - annotations_long.csv   (one row per (input ID x GO/domain/pathway) entry)
 - annotations_full.xlsx  (Excel workbook with 'full' and 'long' sheets)

Usage:
    python functional_annotation_uniprot.py -i "my_proteins.xlsx"

Dependencies:
    pip install requests pandas tqdm openpyxl
"""

from __future__ import annotations
import argparse
import os
import time
import json
import math
import re
from typing import List, Dict, Any, Optional
import requests
import pandas as pd
from tqdm import tqdm

# -------------------------
# Configuration / Defaults
# -------------------------
UNIPROT_BASE = "https://rest.uniprot.org"
SLEEP_BETWEEN_REQUESTS = 0.35  # seconds (be polite)
CACHE_DEFAULT = "uniprot_cache.json"
REQUEST_TIMEOUT = 30  # seconds
MAX_RETRIES = 3
# -------------------------

def is_uniprot_accession(s: str) -> bool:
    """Rudimentary check if string looks like a UniProt accession (SwissProt/TrEMBL style)."""
    if not isinstance(s, str):
        return False
    s = s.strip()
    # Typical patterns: 1 letter + 5 digits/letters (e.g. P12345) OR 6-10 alnum (TrEMBL) - use a permissive regex
    return bool(re.match(r"^[A-NR-Z][0-9A-Z]{5}$", s)) or bool(re.match(r"^[A-Z0-9]{6,10}$", s))

def read_input(input_path: str, sheet: Optional[str | int] = None, id_hints: List[str] = None) -> (List[str], pd.DataFrame, str):
    """Read excel/csv and try to detect ID column. Returns (ids_list, df, chosen_col)."""
    if id_hints is None:
        id_hints = ["uniprot", "accession", "id", "protein", "gene", "symbol", "entry"]

    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")

    if input_path.lower().endswith((".xls", ".xlsx")):
        raw = pd.read_excel(input_path, sheet_name=sheet, engine="openpyxl")
        if isinstance(raw, dict):
            if sheet is None:
                # pick first sheet
                sheet_name = list(raw.keys())[0]
                df = raw[sheet_name]
                print(f"Multiple sheets detected: using first sheet '{sheet_name}'")
            else:
                df = raw.get(sheet, raw[list(raw.keys())[0]])
        else:
            df = raw
    elif input_path.lower().endswith(".csv"):
        df = pd.read_csv(input_path)
    elif input_path.lower().endswith(".tsv") or input_path.lower().endswith(".txt"):
        df = pd.read_csv(input_path, sep="\t")
    else:
        raise ValueError("Unsupported file type. Use .xlsx/.xls/.csv/.tsv/.txt")

    if df is None or df.shape[0] == 0:
        raise ValueError("Input file or chosen sheet is empty")

    cols = df.columns.tolist()
    chosen = None
    for c in cols:
        lc = str(c).lower()
        if any(h in lc for h in id_hints):
            chosen = c
            break

    # Heuristic fallback: choose first column that looks like accessions
    if chosen is None:
        for c in cols:
            sample = df[c].astype(str).dropna().head(20).tolist()
            if any(is_uniprot_accession(s) for s in sample):
                chosen = c
                print(f"Heuristically selected column '{chosen}' as ID column")
                break

    if chosen is None:
        chosen = cols[0]
        print(f"Warning: couldn't detect ID column. Using first column '{chosen}'")

    ids = df[chosen].astype(str).str.strip().replace({"": None, "nan": None, "NaN": None}).dropna().unique().tolist()
    print(f"Using column '{chosen}' with {len(ids)} unique IDs (showing up to 10): {ids[:10]}")
    return ids, df, chosen

# -------------------------
# UniProt query functions
# -------------------------
def request_json_with_retries(url: str, params: dict = None, headers: dict = None, timeout: int = REQUEST_TIMEOUT) -> dict:
    """Request JSON with backoff retries."""
    headers = headers or {"Accept": "application/json"}
    attempt = 0
    while attempt < MAX_RETRIES:
        try:
            r = requests.get(url, params=params, headers=headers, timeout=timeout)
            if r.status_code == 200:
                return r.json()
            elif r.status_code == 404:
                return {"_status": "not_found", "http_status": 404}
            else:
                # transient error maybe; wait and retry
                attempt += 1
                time.sleep(1.2 ** attempt)
        except requests.exceptions.RequestException as e:
            attempt += 1
            time.sleep(1.2 ** attempt)
    # final attempt: raise
    r.raise_for_status()

def fetch_uniprot_entry(accession: str) -> dict:
    """
    Fetch UniProt entry JSON for a given accession using the UniProtKB endpoint.
    Endpoint: GET /uniprotkb/{accession}.json
    """
    url = f"{UNIPROT_BASE}/uniprotkb/{accession}.json"
    data = request_json_with_retries(url)
    # If not found, return a dict indicating not found
    if isinstance(data, dict) and data.get("_status") == "not_found":
        return {"uniprot_accession": accession, "found": False}
    # return the raw JSON
    return {"uniprot_accession": accession, "found": True, "raw": data}

def search_uniprot_by_query(q: str, size: int = 10) -> dict:
    """
    Fallback: search UniProt for a gene symbol or identifier using the search endpoint.
    Example query: 'gene_exact:GENE_NAME AND organism_id:3702' (3702 = Arabidopsis thaliana)
    We'll do a generic search by 'query' parameter.
    """
    url = f"{UNIPROT_BASE}/uniprotkb/search"
    params = {"query": q, "format": "json", "size": size}
    data = request_json_with_retries(url, params=params)
    return data

# -------------------------
# Parser: pull useful fields from UniProt JSON
# -------------------------
def parse_uniprot_json(entry_json: dict) -> dict:
    """Extract convenient fields: gene, protein name, go_terms, ec, domains, pathways, length, sequence."""
    out = {
        "accession": None,
        "primary_gene": None,
        "protein_name": None,
        "go_bp": [],
        "go_mf": [],
        "go_cc": [],
        "ec_numbers": [],
        "domains": [],  # tuples (domain_name, domain_id)
        "pathways": [],
        "sequence_length": None,
        "sequence": None,
    }
    if not entry_json or not entry_json.get("found", False):
        return out

    data = entry_json["raw"]

    # accession(s)
    acs = []
    # 'primaryAccession' key usually available
    if "primaryAccession" in data:
        out["accession"] = data["primaryAccession"]
    # names/protein recommended name
    prot = data.get("proteinDescription", {})
    # recommendedName -> fullName -> value
    try:
        rec = prot.get("recommendedName", {})
        pname = rec.get("fullName", {}).get("value")
        if not pname:
            # fallback to submittedName or alternativeNames
            anames = prot.get("alternativeNames", [])
            if anames:
                pname = anames[0].get("fullName", {}).get("value")
        out["protein_name"] = pname
    except Exception:
        out["protein_name"] = None

    # gene
    gene_obj = data.get("genes", [])
    if gene_obj:
        # pick the first gene's primary name 'value'
        try:
            out["primary_gene"] = gene_obj[0].get("geneName", {}).get("value") or gene_obj[0].get("primaryName", {}).get("value")
        except Exception:
            out["primary_gene"] = None

    # sequence length, sequence
    seqobj = data.get("sequence", {})
    if seqobj:
        out["sequence_length"] = seqobj.get("length")
        out["sequence"] = seqobj.get("value")

    # comments / cross-references: extract GO terms and EC numbers and pathways
    # GO annotations often under 'uniProtKBCrossReferences' or 'comments' sections or 'dbReferences'
    # More reliable: check 'uniProtKBCrossReferences' (dbReference with type 'GO' or 'InterPro', 'Pfam', 'KEGG', 'Reactome')
    xrefs = data.get("uniProtKBCrossReferences", []) or data.get("dbReferences", []) or []
    for xr in xrefs:
        dbtype = xr.get("database") or xr.get("type")
        if not dbtype:
            continue
        dbtype_up = dbtype.upper()
        if dbtype_up == "GO":
            # properties may include 'term' with format 'P:xxxxx' or 'F:' etc.
            # Some schemas: xr['properties'] list of {'key':..., 'value':...}
            props = xr.get("properties", [])
            term = None
            for p in props:
                if p.get("key", "").lower() in {"term", "name"}:
                    term = p.get("value")
            # term looks like 'P:biological process term name' or 'F:...
            if term:
                # split code and name
                m = re.match(r"([PFC]):\s*(.+)", term)
                if m:
                    k, name = m.group(1), m.group(2)
                    if k == "P":
                        out["go_bp"].append(name)
                    elif k == "F":
                        out["go_mf"].append(name)
                    elif k == "C":
                        out["go_cc"].append(name)
                else:
                    # fallback: just append entire term
                    out["go_bp"].append(term)
        elif dbtype_up in {"INTERPRO", "IPR"}:
            interpro_id = xr.get("id") or xr.get("accession")
            interpro_name = None
            # some xrefs include properties with 'entry name'
            props = xr.get("properties", [])
            for p in props:
                if p.get("key", "").lower() in {"name", "entry name", "description"}:
                    interpro_name = p.get("value")
            out["domains"].append((interpro_name or "", interpro_id or ""))
        elif dbtype_up in {"PFAM"}:
            pfam_id = xr.get("id")
            pfam_name = None
            for p in xr.get("properties", []):
                if p.get("key", "").lower() == "entry name":
                    pfam_name = p.get("value")
            out["domains"].append((pfam_name or "", pfam_id or ""))
        elif dbtype_up in {"KEGG", "REACTOME", "PATHWAY"} or dbtype_up in {"EC"}:
            # pathway or EC
            if dbtype_up == "EC":
                out["ec_numbers"].append(xr.get("id") or "")
            else:
                out["pathways"].append((dbtype_up, xr.get("id") or ""))
    # Additionally, scan 'comments' for catalytic activity (EC) or pathway info
    comments = data.get("comments", []) or []
    for c in comments:
        ctype = c.get("commentType", "").lower()
        if ctype == "pathway":
            # pathway entries in comments may include 'reaction' or 'pathway' details
            for e in c.get("pathway", []) if isinstance(c.get("pathway", []), list) else []:
                pname = e.get("name")
                pid = e.get("dbReference", {}).get("id")
                if pname or pid:
                    out["pathways"].append((pname, pid))
        elif ctype == "catalytic activity":
            # catalytic activity may include ecNumber
            cat = c.get("catalyticActivity", {})
            ec = cat.get("ecNumber")
            if ec:
                out["ec_numbers"].append(ec)
    # deduplicate lists
    out["go_bp"] = sorted(list(dict.fromkeys([x for x in out["go_bp"] if x])))
    out["go_mf"] = sorted(list(dict.fromkeys([x for x in out["go_mf"] if x])))
    out["go_cc"] = sorted(list(dict.fromkeys([x for x in out["go_cc"] if x])))
    out["ec_numbers"] = sorted(list(dict.fromkeys([x for x in out["ec_numbers"] if x])))
    out["domains"] = [d for d in out["domains"] if d[1] or d[0]]
    out["domains"] = list(dict.fromkeys(out["domains"]))
    out["pathways"] = list(dict.fromkeys(out["pathways"]))
    return out

# -------------------------
# Caching utilities
# -------------------------
def load_cache(cache_file: str) -> Dict[str, Any]:
    if os.path.exists(cache_file):
        try:
            with open(cache_file, "r", encoding="utf-8") as fh:
                return json.load(fh)
        except Exception:
            return {}
    return {}

def save_cache(cache: Dict[str, Any], cache_file: str):
    with open(cache_file, "w", encoding="utf-8") as fh:
        json.dump(cache, fh, indent=2, ensure_ascii=False)

# -------------------------
# Main annotation pipeline
# -------------------------
def annotate_list(ids: List[str], cache_file: str = CACHE_DEFAULT, pause: float = SLEEP_BETWEEN_REQUESTS, verbose: bool = True) -> pd.DataFrame:
    cache = load_cache(cache_file)
    results = []
    # iterate with progress bar
    for raw_id in tqdm(ids, desc="Annotating", unit="id"):
        identifier = str(raw_id).strip()
        if not identifier:
            continue

        # check cache
        if identifier in cache:
            entry = cache[identifier]
            # already parsed stored
            parsed = entry.get("parsed") if isinstance(entry, dict) else None
            if parsed:
                results.append({"input_id": identifier, **parsed})
                continue

        # determine whether it looks like a UniProt accession
        if is_uniprot_accession(identifier):
            # fetch directly
            try:
                entry_json = fetch_uniprot_entry(identifier)
            except Exception as e:
                print(f"Error fetching {identifier}: {e}")
                entry_json = {"uniprot_accession": identifier, "found": False}
        else:
            # fallback: search the UniProt API for a match (slower)
            # We'll try search: query = identifier
            try:
                search_res = search_uniprot_by_query(identifier, size=3)
                # look for first result and get its primaryAccession if possible
                entries = search_res.get("results", []) if isinstance(search_res, dict) else []
                if entries:
                    # each entry may have 'primaryAccession'
                    first = entries[0]
                    acc = first.get("primaryAccession") or first.get("primaryAccession")
                    if acc:
                        entry_json = fetch_uniprot_entry(acc)
                    else:
                        entry_json = {"uniprot_accession": identifier, "found": False}
                else:
                    entry_json = {"uniprot_accession": identifier, "found": False}
            except Exception as e:
                print(f"Error searching UniProt for {identifier}: {e}")
                entry_json = {"uniprot_accession": identifier, "found": False}

        # parse entry
        parsed = parse_uniprot_json(entry_json)
        # store into cache
        cache[identifier] = {"raw_fetch": entry_json, "parsed": parsed}
        # append with input id
        results.append({"input_id": identifier, **parsed})
        # be polite
        time.sleep(pause)

    # save cache
    save_cache(cache, cache_file)
    # build DataFrame
    df = pd.DataFrame(results)
    # Ensure columns exist
    for c in ["accession", "primary_gene", "protein_name", "go_bp", "go_mf", "go_cc", "ec_numbers", "domains", "pathways", "sequence_length", "sequence"]:
        if c not in df.columns:
            df[c] = None
    return df

def explode_long_table(df_full: pd.DataFrame) -> pd.DataFrame:
    """
    Convert the wide annotation DataFrame into a long tidy table:
    - one row per (input_id, annotation_type, annotation_value)
    Useful for downstream enrichment or filtering.
    """
    rows = []
    for _, r in df_full.iterrows():
        input_id = r["input_id"]
        # GO BP/MF/CC
        for go in r.get("go_bp") or []:
            rows.append({"input_id": input_id, "type": "GO_BP", "value": go})
        for go in r.get("go_mf") or []:
            rows.append({"input_id": input_id, "type": "GO_MF", "value": go})
        for go in r.get("go_cc") or []:
            rows.append({"input_id": input_id, "type": "GO_CC", "value": go})
        for ec in r.get("ec_numbers") or []:
            rows.append({"input_id": input_id, "type": "EC", "value": ec})
        for dom in r.get("domains") or []:
            dom_name, dom_id = dom if isinstance(dom, (list, tuple)) else (str(dom), "")
            rows.append({"input_id": input_id, "type": "DOMAIN", "value": f"{dom_name} ({dom_id})"})
        for p in r.get("pathways") or []:
            name, pid = p if isinstance(p, (list, tuple)) else (str(p), "")
            rows.append({"input_id": input_id, "type": "PATHWAY", "value": f"{name} ({pid})"})
    df_long = pd.DataFrame(rows)
    return df_long

# -------------------------
# CLI
# -------------------------
def main_cli():
    p = argparse.ArgumentParser(description="Functional annotation of protein list using UniProt REST API.")
    p.add_argument("-i", "--input", required=True, help="Input file (xlsx/csv/tsv).")
    p.add_argument("-s", "--sheet", default=None, help="Sheet name or index (for Excel).")
    p.add_argument("-o", "--outdir", default="annotation_results", help="Output directory.")
    p.add_argument("--cache", default=CACHE_DEFAULT, help="Cache JSON file path.")
    p.add_argument("--pause", type=float, default=SLEEP_BETWEEN_REQUESTS, help="Seconds to sleep between requests.")
    args = p.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    ids, df_input, chosen_col = read_input(args.input, sheet=args.sheet)
    df_annotations = annotate_list(ids, cache_file=args.cache, pause=args.pause)

    # Save wide results
    full_csv = os.path.join(outdir, "annotations_full.csv")
    full_xlsx = os.path.join(outdir, "annotations_full.xlsx")
    df_annotations.to_csv(full_csv, index=False)
    try:
        df_annotations.to_excel(full_xlsx, index=False)
    except Exception:
        pass
    print(f"Saved: {full_csv}  (and Excel if supported)")

    # Save long (tidy) table
    df_long = explode_long_table(df_annotations)
    long_csv = os.path.join(outdir, "annotations_long.csv")
    df_long.to_csv(long_csv, index=False)
    print(f"Saved: {long_csv}")

    # Simple summary
    summary = {
        "n_input_ids": len(ids),
        "n_annotated": int(df_annotations['accession'].notna().sum()),
        "n_with_go_terms": int(df_long[df_long['type'].str.startswith("GO")]['input_id'].nunique())
    }
    summary_file = os.path.join(outdir, "annotation_summary.json")
    with open(summary_file, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"Saved summary: {summary_file}")
    print("Done. Outputs in:", outdir)

if __name__ == "__main__":
    main_cli()
