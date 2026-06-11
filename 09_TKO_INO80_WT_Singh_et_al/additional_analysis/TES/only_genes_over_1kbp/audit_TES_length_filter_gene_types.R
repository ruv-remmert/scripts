#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
infile <- args[1]

if (length(args) < 1 || is.na(infile)) {
  stop("Usage: Rscript audit_TES_length_filter_gene_types.R computeMatrix_TES_matrix.gz")
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

mat <- read.table(
  if (grepl("\\.gz$", infile)) gzfile(infile) else infile,
  skip = 1,
  header = FALSE,
  sep = "\t",
  fill = TRUE,
  stringsAsFactors = FALSE
)

mat_gene <- clean_gene_name(mat[, 4])

# ----------------------------
# Gene length from GTF
# (TES BED has 1-bp coords, so length
#  must come from the GTF annotation)
# ----------------------------

gtf_file <- "/media/linuxmac/Storage2/scripts/08_tko_polII/data/genes_sacCer3.gtf"

gtf_genes <- read.table(
  gtf_file, sep = "\t", header = FALSE, comment.char = "#",
  stringsAsFactors = FALSE, quote = ""
)

gtf_genes <- gtf_genes[gtf_genes$V3 == "gene", ]

gtf_length <- data.frame(
  gene = gsub(".*gene_id \"([^\"]+)\".*", "\\1", gtf_genes$V9),
  length = gtf_genes$V5 - gtf_genes$V4,
  stringsAsFactors = FALSE
)

row.names(gtf_length) <- gtf_length$gene

gene_info <- data.frame(
  chr = mat[, 1],
  start = mat[, 2],
  end = mat[, 3],
  gene = mat[, 4],
  gene_clean = mat_gene,
  gene_length = ifelse(
    mat_gene %in% row.names(gtf_length),
    gtf_length[mat_gene, "length"],
    NA
  ),
  stringsAsFactors = FALSE
)

# ----------------------------
# Length filter
# ----------------------------

gene_info$filter_status <- ifelse(
  is.na(gene_info$gene_length),
  "unknown_length",
  ifelse(
    gene_info$gene_length > 1000,
    "kept_gt1000bp",
    "removed_le1000bp"
  )
)

# ----------------------------
# Gene category files
# ----------------------------

group_files <- c(
  TES_convergent = "../convergent/TES_regSorted_convergent.bed",
  TES_lonely    = "../lonely/TES_regSorted_lonely.bed",
  TES_tandem    = "../tandem/TES_regSorted_tandem.bed",
  TSS_divergent = "TSS_regSorted_divergent.bed",
  TSS_lonely    = "TSS_regSorted_lonely.bed",
  TSS_tandem    = "TSS_regSorted_tandem.bed"
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
  file = "genes_removed_by_TES_gt1000bp_filter.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  kept_genes,
  file = "genes_kept_by_TES_gt1000bp_filter.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------
# Summary table (full combined type)
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
  file = "genes_removed_by_TES_gt1000bp_filter_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------
# Per-category totals (TES_xxx and TSS_xxx aggregated separately)
# ----------------------------

group_prefixes <- list(
  TES = c("TES_convergent", "TES_lonely", "TES_tandem"),
  TSS = c("TSS_divergent", "TSS_lonely", "TSS_tandem")
)

summarise_category <- function(gene_info_df, labels) {

  gene_info_df$cat_col <- sapply(strsplit(as.character(gene_info_df$gene_type), ";"), function(parts) {
    matching <- parts[parts %in% labels]
    if (length(matching) == 1) matching[1] else NA_character_
  })

  kept    <- gene_info_df[gene_info_df$filter_status == "kept_gt1000bp", ]
  removed <- gene_info_df[gene_info_df$filter_status == "removed_le1000bp", ]

  total_counts    <- table(gene_info_df$cat_col)
  kept_counts     <- table(kept$cat_col)
  removed_counts  <- table(removed$cat_col)

  df <- data.frame(
    gene_type        = labels,
    total            = as.integer(total_counts[labels]),
    kept_gt1000bp    = as.integer(kept_counts[labels]),
    removed_le1000bp = as.integer(removed_counts[labels]),
    stringsAsFactors = FALSE
  )
  df[is.na(df)] <- 0L
  df
}

tes_summary <- summarise_category(gene_info, group_prefixes$TES)
tss_summary <- summarise_category(gene_info, group_prefixes$TSS)

write.table(
  tes_summary,
  file = "genes_removed_by_TES_gt1000bp_filter_TES_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tss_summary,
  file = "genes_removed_by_TES_gt1000bp_filter_TSS_summary.tsv",
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
cat("Genes with unknown length:", sum(gene_info$filter_status == "unknown_length"), "\n")

cat("\n--- TES gene type totals ---\n")
print(tes_summary)

cat("\n--- TSS gene type totals ---\n")
print(tss_summary)

cat("\nRemoved genes by full type:\n")
print(
  table(removed_genes$gene_type)
)

cat("\nKept genes by full type:\n")
print(
  table(kept_genes$gene_type)
)

cat("\nSaved: genes_removed_by_TES_gt1000bp_filter.tsv\n")
cat("Saved: genes_kept_by_TES_gt1000bp_filter.tsv\n")
cat("Saved: genes_removed_by_TES_gt1000bp_filter_summary.tsv\n")
cat("Saved: genes_removed_by_TES_gt1000bp_filter_TES_summary.tsv\n")
cat("Saved: genes_removed_by_TES_gt1000bp_filter_TSS_summary.tsv\n")