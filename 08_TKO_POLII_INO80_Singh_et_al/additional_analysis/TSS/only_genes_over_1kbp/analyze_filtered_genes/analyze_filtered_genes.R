#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
infile <- args[1]

if (length(args) < 1 || is.na(infile)) {
  stop("Usage: Rscript audit_TSS_length_filter_gene_types.R computeMatrix2plot_mono.txt.gz")
}

read_matrix_rows <- function(path) {
  read.table(
    if (grepl("\\.gz$", path)) gzfile(path) else path,
    skip = 1,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
}

clean_gene_name <- function(x) {
  sub("_mRNA$", "", x)
}

read_bed_genes <- function(path) {
  if (!file.exists(path)) {
    warning("Missing BED file: ", path)
    return(character())
  }

  bed <- read.table(
    path,
    sep = "\t",
    header = FALSE,
    fill = TRUE,
    stringsAsFactors = FALSE
  )

  clean_gene_name(bed[, 4])
}

# ----------------------------
# Read input matrix
# ----------------------------

mat <- read_matrix_rows(infile)

gene_info <- data.frame(
  chr = mat[, 1],
  start = mat[, 2],
  end = mat[, 3],
  gene = mat[, 4],
  gene_clean = clean_gene_name(mat[, 4]),
  gene_length = abs(mat[, 3] - mat[, 2]),
  stringsAsFactors = FALSE
)

# ----------------------------
# Length filter
# ----------------------------

gene_info$filter_status <- ifelse(
  gene_info$gene_length > 1000,
  "kept_gt1000bp",
  "removed_le1000bp"
)

# ----------------------------
# Gene category files
# ----------------------------

group_files <- c(
  TES_convergent = "TES_regSorted_convergent.bed",
  TES_lonely = "TES_regSorted_lonely.bed",
  TES_tandem = "TES_regSorted_tandem.bed",
  TSS_divergent = "TSS_regSorted_divergent.bed",
  TSS_lonely = "TSS_regSorted_lonely.bed",
  TSS_tandem = "TSS_regSorted_tandem.bed"
)

group_genes <- lapply(
  group_files,
  read_bed_genes
)

# ----------------------------
# Assign gene types
# ----------------------------

gene_info$gene_type <- "unclassified"

for (grp in names(group_genes)) {

  hits <- gene_info$gene_clean %in% group_genes[[grp]]

  gene_info$gene_type[hits] <- ifelse(
    gene_info$gene_type[hits] == "unclassified",
    grp,
    paste(
      gene_info$gene_type[hits],
      grp,
      sep = ";"
    )
  )
}

# ----------------------------
# Split removed / kept genes
# ----------------------------

removed_genes <- gene_info[
  gene_info$filter_status == "removed_le1000bp",
]

kept_genes <- gene_info[
  gene_info$filter_status == "kept_gt1000bp",
]

# ----------------------------
# Order outputs
# ----------------------------

removed_genes <- removed_genes[
  order(
    removed_genes$gene_type,
    removed_genes$gene_clean,
    removed_genes$chr,
    removed_genes$start
  ),
]

kept_genes <- kept_genes[
  order(
    kept_genes$gene_type,
    kept_genes$gene_clean,
    kept_genes$chr,
    kept_genes$start
  ),
]

# ----------------------------
# Write detailed outputs
# ----------------------------

write.table(
  removed_genes,
  file = "genes_removed_by_TSS_gt1000bp_filter.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  kept_genes,
  file = "genes_kept_by_TSS_gt1000bp_filter.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------
# Summary table
# ----------------------------

summary_table <- as.data.frame(
  table(
    filter_status = gene_info$filter_status,
    gene_type = gene_info$gene_type
  )
)

summary_table <- summary_table[
  order(
    summary_table$filter_status,
    summary_table$gene_type
  ),
]

write.table(
  summary_table,
  file = "genes_removed_by_TSS_gt1000bp_filter_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------
# Console output
# ----------------------------

cat("Genes total:", nrow(gene_info), "\n")
cat("Genes kept (>1000 bp):", nrow(kept_genes), "\n")
cat("Genes removed (<=1000 bp):", nrow(removed_genes), "\n")

cat("\nRemoved genes by type:\n")
print(
  table(removed_genes$gene_type)
)

cat("\nKept genes by type:\n")
print(
  table(kept_genes$gene_type)
)

cat("\nSaved: genes_removed_by_TSS_gt1000bp_filter.tsv\n")
cat("Saved: genes_kept_by_TSS_gt1000bp_filter.tsv\n")
cat("Saved: genes_removed_by_TSS_gt1000bp_filter_summary.tsv\n")