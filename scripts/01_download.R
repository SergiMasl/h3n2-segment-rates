## 01_download.R ------------------------------------------------------------
## Scripted, resumable download of human seasonal H3N2 influenza A sequences
## and metadata from NCBI (GenBank via E-utilities).
##
## Strategy (3 stages), because "complete genome" is not a reliable Entrez
## filter for influenza -- most records are deposited as individual segments:
##
##   Stage A  For each year, search segment 4 (HA) records and harvest strain
##            names. HA is the best-sampled segment, so it gives the widest
##            candidate pool. Oversample: we need surplus because many strains
##            will fail the all-8-segments test in Stage B.
##   Stage B  For each candidate strain, query all 8 segments by strain name.
##            Keep only strains with a usable record for every segment.
##   Stage C  efetch GenBank flat files for the kept accessions, parse
##            metadata (strain, collection_date, country, segment, passage),
##            and write per-segment FASTA + a single metadata CSV.
##
## Everything is cached to data/raw/cache/*.rds, so a re-run resumes rather
## than re-hitting NCBI. Delete the cache dir to force a clean pull.
##
## Usage:  Rscript scripts/01_download.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(rentrez)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(readr)
})

## ---- Config ---------------------------------------------------------------

YEAR_MIN        <- 2000
YEAR_MAX        <- 2025
CANDIDATES_PER_YEAR <- 40    # oversample; Stage B will cut this down hard
TARGET_PER_YEAR     <- 6     # what 02_match_subsample.R will keep
MIN_SEG_LEN     <- 700       # drop obvious fragments (shortest segment, NS, ~890bp)
EMAIL           <- "YOUR_EMAIL@uncc.edu"   # <-- set this; NCBI requires it
CACHE_DIR       <- "data/raw/cache"
OUT_DIR         <- "data/raw"
SEED            <- 20261119

## NCBI API key raises the rate limit from 3 to 10 requests/sec.
## Get one at https://www.ncbi.nlm.nih.gov/account/settings/ and then either
## export NCBI_API_KEY=... in your shell, or paste it here.
API_KEY <- Sys.getenv("NCBI_API_KEY", unset = "")
if (nzchar(API_KEY)) {
  set_entrez_key(API_KEY)
  RATE_DELAY <- 0.11
  message("Using NCBI API key (10 req/sec).")
} else {
  RATE_DELAY <- 0.34
  message("No NCBI_API_KEY found -- throttling to 3 req/sec. This will be slow.")
}

SEGMENTS <- tibble::tribble(
  ~segment, ~gene, ~class,
  1L, "PB2", "internal",
  2L, "PB1", "internal",
  3L, "PA",  "internal",
  4L, "HA",  "surface",
  5L, "NP",  "internal",
  6L, "NA",  "surface",
  7L, "M",   "internal",
  8L, "NS",  "internal"
)

dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,   recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

## ---- Helpers --------------------------------------------------------------

#' Cache a value to RDS; return cached value if present.
cached <- function(key, expr) {
  path <- file.path(CACHE_DIR, paste0(key, ".rds"))
  if (file.exists(path)) {
    message("  [cache] ", key)
    return(readRDS(path))
  }
  val <- force(expr)
  saveRDS(val, path)
  val
}

#' esearch with retry + throttle.
safe_search <- function(term, retmax = 200, db = "nuccore", tries = 4) {
  for (i in seq_len(tries)) {
    res <- try(entrez_search(db = db, term = term, retmax = retmax,
                             use_history = FALSE), silent = TRUE)
    Sys.sleep(RATE_DELAY)
    if (!inherits(res, "try-error")) return(res)
    message("    retry ", i, " on search: ", substr(term, 1, 60))
    Sys.sleep(2 * i)
  }
  NULL
}

#' efetch with retry + throttle.
safe_fetch <- function(ids, rettype, retmode = "text", db = "nuccore", tries = 4) {
  for (i in seq_len(tries)) {
    res <- try(entrez_fetch(db = db, id = ids, rettype = rettype,
                            retmode = retmode), silent = TRUE)
    Sys.sleep(RATE_DELAY)
    if (!inherits(res, "try-error")) return(res)
    message("    retry ", i, " on fetch of ", length(ids), " ids")
    Sys.sleep(2 * i)
  }
  NULL
}

#' Pull the strain designation out of a flu record title.
#' Titles look like:
#'   "Influenza A virus (A/Texas/50/2012(H3N2)) segment 4 hemagglutinin ..."
extract_strain <- function(title) {
  m <- str_match(title, "\\(\\s*(A/[^()]+?)\\s*\\(H3N2\\)\\s*\\)")
  str_trim(m[, 2])
}

## ---- Stage A: harvest candidate strain names, year by year -----------------

message("\n=== Stage A: harvesting candidate strains by year ===")

candidates <- map_dfr(YEAR_MIN:YEAR_MAX, function(yr) {
  cached(sprintf("stageA_%d", yr), {
    message("  ", yr)
    term <- sprintf(
      paste0('"Influenza A virus"[Organism] AND H3N2[All Fields] ',
             'AND "segment 4"[Title] AND "complete cds"[Title] ',
             'AND "Homo sapiens"[Host] ',
             'AND ("%d/01/01"[Collection Date] : "%d/12/31"[Collection Date])'),
      yr, yr)

    sr <- safe_search(term, retmax = 500)
    if (is.null(sr) || sr$count == 0) {
      message("    no hits for ", yr)
      return(tibble(year = integer(), strain = character()))
    }

    ids <- sample(sr$ids, min(length(sr$ids), CANDIDATES_PER_YEAR * 3))
    summ <- entrez_summary(db = "nuccore", id = ids)
    Sys.sleep(RATE_DELAY)
    titles <- extract_from_esummary(summ, "title", simplify = TRUE)

    strains <- unique(na.omit(extract_strain(unlist(titles))))
    strains <- head(sample(strains), CANDIDATES_PER_YEAR)

    tibble(year = yr, strain = strains)
  })
})

message("Stage A: ", nrow(candidates), " candidate strains across ",
        n_distinct(candidates$year), " years.")

## ---- Stage B: keep only strains with all 8 segments ------------------------

message("\n=== Stage B: checking segment completeness per strain ===")

#' For one strain, find one accession per segment. Returns NULL if incomplete.
segments_for_strain <- function(strain) {
  hits <- map_dfr(SEGMENTS$segment, function(sg) {
    term <- sprintf('"%s"[All Fields] AND "segment %d"[Title] AND "Influenza A virus"[Organism]',
                    strain, sg)
    sr <- safe_search(term, retmax = 10)
    if (is.null(sr) || sr$count == 0) return(tibble())
    tibble(segment = sg, uid = sr$ids[1])
  })
  if (nrow(hits) < 8) return(NULL)
  hits$strain <- strain
  hits
}

complete <- cached("stageB_complete", {
  out <- vector("list", nrow(candidates))
  for (i in seq_len(nrow(candidates))) {
    s <- candidates$strain[i]
    if (i %% 10 == 0) message("  ", i, "/", nrow(candidates),
                              "  (kept so far: ", sum(!map_lgl(out, is.null)), ")")
    out[[i]] <- segments_for_strain(s)
  }
  bind_rows(out) |> left_join(candidates, by = "strain")
})

kept_strains <- n_distinct(complete$strain)
message("Stage B: ", kept_strains, " strains have all 8 segments.")
if (kept_strains < TARGET_PER_YEAR * (YEAR_MAX - YEAR_MIN + 1)) {
  warning("Fewer complete strains than the target subsample needs. ",
          "Raise CANDIDATES_PER_YEAR and re-run (delete the Stage B cache first).")
}

## ---- Stage C: fetch GenBank, parse metadata, write FASTA -------------------

message("\n=== Stage C: fetching records and parsing metadata ===")

#' Split a GenBank multi-record text blob into individual records.
split_gb <- function(txt) {
  recs <- str_split(txt, "\n//\n")[[1]]
  recs[str_detect(recs, "^LOCUS")]
}

#' Pull the qualifiers we need out of one GenBank flat record.
parse_gb <- function(rec) {
  acc <- str_match(rec, "VERSION\\s+(\\S+)")[, 2]
  defn <- str_match(rec, "DEFINITION\\s+((?:.|\\n)*?)\\nACCESSION")[, 2]
  defn <- str_squish(defn)

  qual <- function(tag) {
    m <- str_match(rec, paste0('/', tag, '="((?:[^"]|\\n)*)"'))
    if (is.na(m[, 2])) NA_character_ else str_squish(m[, 2])
  }

  seg <- suppressWarnings(as.integer(qual("segment")))
  if (is.na(seg)) {
    seg <- suppressWarnings(as.integer(str_match(defn, "segment (\\d)")[, 2]))
  }

  ## Sequence: everything after ORIGIN, digits and whitespace stripped.
  sq <- str_match(rec, "ORIGIN((?:.|\\n)*)$")[, 2]
  sq <- toupper(str_remove_all(sq %||% "", "[^acgtnACGTN]"))

  tibble(
    accession       = acc,
    definition      = defn,
    strain          = qual("strain") %||% extract_strain(defn),
    segment         = seg,
    collection_date = qual("collection_date"),
    country         = qual("geo_loc_name") %||% qual("country"),
    host            = qual("host"),
    note            = qual("note"),
    lab_host        = qual("lab_host"),
    seq             = sq,
    length          = nchar(sq)
  )
}

`%||%` <- function(a, b) if (is.null(a) || is.na(a) || !nzchar(a)) b else a

records <- cached("stageC_records", {
  uids <- unique(complete$uid)
  batches <- split(uids, ceiling(seq_along(uids) / 100))
  map_dfr(seq_along(batches), function(i) {
    message("  batch ", i, "/", length(batches))
    txt <- safe_fetch(batches[[i]], rettype = "gb")
    if (is.null(txt)) return(tibble())
    map_dfr(split_gb(txt), parse_gb)
  })
})

message("Stage C: parsed ", nrow(records), " records.")

## ---- Passage-history flag --------------------------------------------------
## Egg adaptation introduces artifactual HA mutations in antigenic sites, so
## flag these now; 02_match_subsample.R prefers clinical/original specimens.

records <- records |>
  mutate(
    passage_text = paste(coalesce(note, ""), coalesce(lab_host, "")),
    passage = case_when(
      str_detect(passage_text, regex("egg|E\\d|AM\\d", ignore_case = TRUE)) ~ "egg",
      str_detect(passage_text, regex("MDCK|SIAT|cell|C\\d", ignore_case = TRUE)) ~ "cell",
      str_detect(passage_text, regex("original|clinical|direct", ignore_case = TRUE)) ~ "clinical",
      TRUE ~ "unknown"
    )
  ) |>
  select(-passage_text)

## ---- Quality filters -------------------------------------------------------

before <- nrow(records)
records <- records |>
  filter(
    !is.na(segment), segment %in% 1:8,
    !is.na(collection_date),
    length >= MIN_SEG_LEN,
    str_count(seq, "N") / length < 0.01     # <1% ambiguous bases
  )
message("Quality filters: ", before, " -> ", nrow(records), " records.")

## ---- Write outputs ---------------------------------------------------------

metadata <- records |>
  select(accession, strain, segment, collection_date, country, host, passage, length) |>
  left_join(SEGMENTS, by = "segment") |>
  arrange(segment, strain)

write_csv(metadata, file.path(OUT_DIR, "metadata_raw.csv"))
message("Wrote ", file.path(OUT_DIR, "metadata_raw.csv"), " (", nrow(metadata), " rows)")

## Per-segment FASTA. Header = accession|strain|collection_date so downstream
## steps can parse tip dates without a second lookup.
for (sg in SEGMENTS$segment) {
  sub <- filter(records, segment == sg)
  if (nrow(sub) == 0) next
  gene <- SEGMENTS$gene[SEGMENTS$segment == sg]
  hdr <- sprintf(">%s|%s|%s", sub$accession,
                 str_replace_all(sub$strain, "[ /]", "_"),
                 sub$collection_date)
  lines <- as.vector(rbind(hdr, sub$seq))
  path <- file.path(OUT_DIR, sprintf("segment%d_%s_raw.fasta", sg, gene))
  writeLines(lines, path)
  message("Wrote ", path, " (", nrow(sub), " sequences)")
}

message("\nDone. Next: scripts/02_match_subsample.R")
