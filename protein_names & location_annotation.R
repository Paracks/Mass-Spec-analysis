# Required packages
packages <- c("httr2", "readr", "readxl", "openxlsx", "dplyr", "stringr", "tibble")
missing <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(missing)) install.packages(missing)
library(httr2)
library(readr)
library(readxl)
library(openxlsx)
library(dplyr)
library(stringr)
library(tibble)

# ---------- USER SETTINGS ----------
input_path <- "relaxed_interactors.csv"
output_csv <- "Final_analysed_relaxed_interactors_subcellular_locations.csv"
output_xlsx <- "Final_analysed_relaxed_interactors_subcellular_locations.xlsx"
accession_col <- NA
# -----------------------------------

# Helper to read input
read_input_table <- function(path) {
  ext <- tools::file_ext(path) %>% tolower()
  if (ext %in% c("csv")) {
    df <- read_csv(path, show_col_types = FALSE, col_types = cols(.default = "c"))
  } else if (ext %in% c("tsv", "txt")) {
    df <- read_tsv(path, show_col_types = FALSE, col_types = cols(.default = "c"))
  } else if (ext %in% c("xls", "xlsx")) {
    df <- read_excel(path)
    df <- as.data.frame(df)
    df[] <- lapply(df, as.character)
  } else {
    stop("Unsupported file extension. Please export your .numbers file to CSV or XLSX and try again.")
  }
  return(df)
}

df <- read_input_table(input_path)

# Auto-detect accession column if needed
if (is.na(accession_col)) {
  possible_names <- c("uniprot", "uniprot_acc", "accession", "acc", "uniprot_accession", "uniprot_id", "Entry", "Accession")
  found <- intersect(tolower(names(df)), possible_names)
  if (length(found) == 0) {
    nm <- names(df)
    matches <- sapply(nm, function(x) any(str_detect(tolower(x), c("uniprot", "accession", "acc"))))
    if (any(matches)) {
      accession_col <- nm[which(matches)[1]]
      message("Auto-detected accession column: ", accession_col)
    } else {
      stop("Could not auto-detect accession column. Please set 'accession_col' to the column name containing UniProt accessions.")
    }
  } else {
    accession_col <- names(df)[which(tolower(names(df)) == found[1])[1]]
    message("Auto-detected accession column: ", accession_col)
  }
} else {
  if (!(accession_col %in% names(df))) stop("Provided accession_col not found in table column names.")
}

# Prepare accession list
accs <- df[[accession_col]] %>% as.character()
accs <- unique(accs[!is.na(accs) & accs != ""])
message("Unique accessions to query: ", length(accs))

# Query function (request protein_name + subcellular location)
query_uniprot_batch <- function(accession_vector, size = 500) {
  q <- paste0("accession:(", paste(accession_vector, collapse = " OR "), ")")
  base <- "https://rest.uniprot.org/uniprotkb/search"
  url <- paste0(base, "?query=", URLencode(q, reserved = TRUE),
                "&fields=accession,cc_subcellular_location,protein_name&format=tsv&size=", size)
  resp <- request(url) %>% req_perform()
  if (resp$status_code != 200L) {
    warning("UniProt request returned status ", resp$status_code)
    return(NULL)
  }
  text <- resp %>% resp_body_string()
  if (nchar(text) == 0) return(NULL)
  # force all columns to character
  df_res <- read_tsv(file = I(text), show_col_types = FALSE, col_types = cols(.default = "c"))
  return(df_res)
}

batch_size <- 100
batches <- split(accs, ceiling(seq_along(accs)/batch_size))

all_hits <- tibble()
for (i in seq_along(batches)) {
  batch <- batches[[i]]
  message(sprintf("Querying batch %d/%d ( %d accessions ) ...", i, length(batches), length(batch)))
  res <- tryCatch({
    query_uniprot_batch(batch, size = length(batch))
  }, error = function(e) {
    warning("Error in batch ", i, ": ", e$message)
    NULL
  })
  if (!is.null(res) && nrow(res) > 0) {
    all_hits <- bind_rows(all_hits, res)
  }
  Sys.sleep(1)
}

if (nrow(all_hits) == 0) {
  stop("No mapping results returned from UniProt. Check accession values and try again.")
}

message("Rows returned by UniProt (including possible repeated headers): ", nrow(all_hits))
message("Columns returned by UniProt (raw): ", paste(names(all_hits), collapse = ", "))

# Standardize names
names(all_hits) <- make.names(names(all_hits))

# Remove accidental header-rows that sometimes appear when concatenating TSVs:
# e.g., rows where a column that should be accession equals "Entry" or "accession"
acc_candidate_names <- c("accession", "Entry", "primaryAccession")
# find the actual accession-like column in returned data
acc_col_candidates <- intersect(names(all_hits), tolower(acc_candidate_names) %>% make.names())
if (length(acc_col_candidates) == 0) {
  # fallback: assume first column is accession
  uni_acc_colname <- names(all_hits)[1]
} else {
  # pick first candidate found
  uni_acc_colname <- acc_col_candidates[1]
}

# Drop rows where accession column equals header strings (case-insensitive)
header_rows <- which(tolower(all_hits[[uni_acc_colname]]) %in% c("entry", "accession", "primaryaccession"))
if (length(header_rows) > 0) {
  message("Dropping ", length(header_rows), " accidental header-like rows from UniProt results.")
  all_hits <- all_hits[-header_rows, , drop = FALSE]
}

# identify probable subcellular and protein name columns heuristically
possible_subcell_cols <- intersect(names(all_hits), c("cc_subcellular_location", "Subcellular.location.CC.", "Subcellular.location..CC.", "Subcellular.location.CC"))
possible_protein_cols <- intersect(names(all_hits), c("protein_name", "Protein.names", "Protein.names..full", "Protein.names.full"))

# heuristics if not found exactly
if (length(possible_subcell_cols) == 0) {
  cand <- names(all_hits)[which(str_detect(tolower(names(all_hits)), "subcell") | str_detect(tolower(names(all_hits)), "location"))]
  subcell_col <- if (length(cand)>0) cand[1] else names(all_hits)[2]
} else subcell_col <- possible_subcell_cols[1]

if (length(possible_protein_cols) == 0) {
  cand2 <- names(all_hits)[which(str_detect(tolower(names(all_hits)), "protein"))]
  protein_col <- if (length(cand2)>0) cand2[1] else NA_character_
} else protein_col <- possible_protein_cols[1]

message("Using accession column from UniProt results: ", uni_acc_colname)
message("Detected subcellular column: ", subcell_col)
message("Detected protein name column: ", protein_col)

# Keep only rows that have a non-empty accession in UniProt response
all_hits <- all_hits %>% filter(!is.na(.data[[uni_acc_colname]]) & .data[[uni_acc_colname]] != "")

# Build mapped data frame with vectorized extraction
# Safely extract columns (if protein_col missing, create NA column)
if (is.na(protein_col) || !(protein_col %in% names(all_hits))) {
  all_hits$protein_name_raw <- NA_character_
} else {
  all_hits$protein_name_raw <- as.character(all_hits[[protein_col]])
}
# subcell raw
all_hits$subcellular_location_raw <- as.character(all_hits[[subcell_col]])
# accession standardized
all_hits$uni_acc_std <- as.character(all_hits[[uni_acc_colname]])

# Vectorized cleanup functions using vapply
clean_protein_name <- function(x) {
  if (is.na(x) || x == "") return(NA_character_)
  # take first semicolon-separated segment
  segs <- str_split(x, ";")[[1]]
  first_seg <- segs[1] %>% str_trim()
  # remove parenthetical content
  first_seg <- str_replace_all(first_seg, "\\s*\\(.*?\\)", "") %>% str_trim()
  first_seg <- str_squish(first_seg)
  if (first_seg == "") NA_character_ else first_seg
}
clean_subcell <- function(x) {
  if (is.na(x) || x == "") return(NA_character_)
  s <- gsub("\\s+", " ", x)
  s <- str_squish(s)
  if (s == "") NA_character_ else s
}

mapped <- tibble(
  uni_acc = all_hits$uni_acc_std,
  protein_name_raw = all_hits$protein_name_raw,
  subcellular_location_raw = all_hits$subcellular_location_raw
) %>%
  mutate(
    protein_name = vapply(protein_name_raw, FUN = clean_protein_name, FUN.VALUE = character(1), USE.NAMES = FALSE),
    subcellular_location = vapply(subcellular_location_raw, FUN = clean_subcell, FUN.VALUE = character(1), USE.NAMES = FALSE)
  ) %>%
  select(uni_acc, protein_name, subcellular_location)

message("Unique UniProt accessions mapped: ", length(unique(mapped$uni_acc)))

# Merge back into original df
df[[accession_col]] <- as.character(df[[accession_col]])
mapped$uni_acc <- as.character(mapped$uni_acc)

df_annot <- df %>%
  left_join(mapped, by = setNames("uni_acc", accession_col))

# Report unmapped rows
missing_idx <- which(is.na(df_annot$subcellular_location) | df_annot$subcellular_location == "")
if (length(missing_idx) > 0) {
  message(length(missing_idx), " rows have no subcellular location.")
}
missing_proteins <- which(is.na(df_annot$protein_name) | df_annot$protein_name == "")
if (length(missing_proteins) > 0) {
  message(length(missing_proteins), " rows have no protein name returned from UniProt.")
}

# Write outputs
write_csv(df_annot, output_csv)
openxlsx::write.xlsx(df_annot, output_xlsx)

message("Done. Annotated files written to: ", output_csv, " and ", output_xlsx)
