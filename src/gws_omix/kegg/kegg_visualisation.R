#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

gene_csv    <- args[[1]]
specie      <- args[[2]]
pathway_txt <- args[[3]]
kegg_dir    <- args[[4]]

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(pathview)
})

df <- read_csv(gene_csv, show_col_types = FALSE)

# Ensure first column is kegg_short if not already named
if (!("kegg_short" %in% colnames(df))) {
  colnames(df)[1] <- "kegg_short"
}

df <- df %>% filter(!is.na(kegg_short), kegg_short != "")

# Always produce an output file even if empty -> important for the Python task importer
write_empty_errors <- function() {
  write_csv(
    data.frame(pathway=character(), status=character(), message=character()),
    "list_pathway_error.csv"
  )
  # optional: full log file too
  write_csv(
    data.frame(pathway=character(), status=character(), message=character()),
    "list_pathway_log.csv"
  )
}

if (nrow(df) == 0) {
  # no genes
  write_csv(
    data.frame(pathway=character(), status="ERROR", message="gene_kegg.csv empty"),
    "list_pathway_error.csv"
  )
  write_csv(
    data.frame(pathway=character(), status="ERROR", message="gene_kegg.csv empty"),
    "list_pathway_log.csv"
  )
  quit(status=1)
}

# Build gene.data
# - If extra columns exist: they are fold-changes (multi-FC supported by pathview -> multi-state)
# - Else: dummy vector of 1s
if (ncol(df) > 1) {
  gene.data <- as.matrix(df[ , setdiff(colnames(df), "kegg_short"), drop=FALSE])
  rownames(gene.data) <- df$kegg_short
} else {
  vals <- rep(1, nrow(df))
  names(vals) <- df$kegg_short
  gene.data <- vals
}

# Read pathways list
pids <- readLines(pathway_txt, warn=FALSE)
pids <- trimws(pids)
pids <- pids[pids != ""]

if (length(pids) == 0) {
  write_csv(
    data.frame(pathway=character(), status="ERROR", message="No pathways provided"),
    "list_pathway_error.csv"
  )
  write_csv(
    data.frame(pathway=character(), status="ERROR", message="No pathways provided"),
    "list_pathway_log.csv"
  )
  quit(status=1)
}

log_rows <- data.frame(
  pathway=character(), status=character(), message=character(),
  stringsAsFactors=FALSE
)

for (pid in pids) {
  status <- "OK"
  msg <- "OK"

  tryCatch({
    # Run pathview
    pv <- pathview(
      gene.data   = gene.data,
      pathway.id  = pid,
      species     = specie,
      gene.idtype = "kegg",
      kegg.native = TRUE,
      same.layer  = FALSE,
      kegg.dir    = kegg_dir,
      out.suffix  = "pathview"
    )

    # Detect "PNG created but no gene colored"
    # If pv$plot.data.gene is all NA -> nothing colored
    if (!is.null(pv$plot.data.gene)) {
      mol <- pv$plot.data.gene
      if (is.data.frame(mol) || is.matrix(mol)) {
        # numeric columns only (in case of multi-state)
        mol_num <- suppressWarnings(as.matrix(mol))
        if (all(is.na(mol_num))) {
          status <- "NO_GENE_MAPPED"
          msg <- "PNG created but no gene colored (all mol.data NA)"
        }
      }
    }

  }, error=function(e) {
    status <<- "ERROR"
    msg <<- as.character(e$message)
  })

  log_rows <- rbind(
    log_rows,
    data.frame(pathway=pid, status=status, message=msg, stringsAsFactors=FALSE)
  )
}

# Always write the full log (optional but useful)
write_csv(log_rows, "list_pathway_log.csv")

# Write ONLY failures to list_pathway_error.csv (this is what you asked)
errors_only <- log_rows %>% filter(status != "OK")
write_csv(errors_only, "list_pathway_error.csv")

# If you want the task to "fail" when there are R errors (not NO_GENE_MAPPED),
# you can uncomment this:
# if (any(errors_only$status == "ERROR")) quit(status=1)

quit(status=0)
