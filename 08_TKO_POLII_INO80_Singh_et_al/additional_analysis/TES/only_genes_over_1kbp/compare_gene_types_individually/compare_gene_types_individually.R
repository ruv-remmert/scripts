#!/usr/bin/env Rscript

library(RColorBrewer)
library(jsonlite)

######## DATA LOADING ###########
args <- commandArgs(TRUE)
infile <- args[1]

if (length(args) < 1 || is.na(infile)) {
  stop("Usage: Rscript compare_gene_types_individually.R <computeMatrix_matrix.gz>")
}

con <- if (grepl("\\.gz$", infile)) gzfile(infile, "rt") else file(infile, "r")
header <- readLines(con, n = 1, warn = FALSE)
close(con)

header <- header[1]
meta <- fromJSON(sub("^@", "", header))

upstream          <- meta$upstream[1]
downstream        <- meta$downstream[1]
binSize           <- meta$`bin size`[1]
sample_labels     <- meta$sample_labels
sample_boundaries <- meta$sample_boundaries

names.pre <- gsub("_monoNucs_profile", "", sample_labels)
names.pre <- gsub("_1\\.fastq\\.gz", "", names.pre)

mat <- read.table(
  if (grepl("\\.gz$", infile)) gzfile(infile) else infile,
  skip = 1,
  header = FALSE,
  sep = "\t",
  fill = TRUE,
  stringsAsFactors = FALSE
)

# --- Gene length filter (GTF) ---
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

mat_gene <- sub("_mRNA$", "", mat[, 4])
keep_length <- mat_gene %in% gtf_length$gene[gtf_length$length > 1000]

cat("Genes before filtering:", nrow(mat), "\n")
cat("Genes >1000 bp:", sum(keep_length), "\n")

# --- Gene type BED files ---
bed_dir <- "."

gene_types <- c(
  convergent = "TES_regSorted_convergent.bed",
  isolated    = "TES_regSorted_lonely.bed",
  tandem    = "TES_regSorted_tandem.bed"
)

type_genes <- lapply(gene_types, function(bed) {
  df <- read.table(
    file.path(bed_dir, bed),
    sep = "\t", header = FALSE, stringsAsFactors = FALSE
  )
  sub("_mRNA$", "", df[, 4])
})

# --- Compute average profiles per gene type ---
baseline_bins <- 1:50
signal_all <- mat[, 7:ncol(mat)]

compute_profiles <- function(rows, signal, sample_boundaries, sample_labels, names.pre) {
  sub_sig <- signal[rows, , drop = FALSE]
  profiles <- list()
  for (i in seq_along(sample_labels)) {
    col_start <- sample_boundaries[i] + 1
    col_end   <- sample_boundaries[i + 1]
    samp_mat <- sub_sig[, col_start:col_end, drop = FALSE]
    vals <- colMeans(samp_mat, na.rm = TRUE)
    baseline <- mean(vals[baseline_bins], na.rm = TRUE)
    profiles[[names.pre[i]]] <- vals / baseline
  }
  profiles
}

profiles_by_type <- list()

for (tp in names(gene_types)) {
  keep_type <- mat_gene %in% type_genes[[tp]]
  rows <- which(keep_length & keep_type)
  cat(tp, "genes >1000 bp:", length(rows), "\n")
  profiles_by_type[[tp]] <- compute_profiles(
    rows, signal_all, sample_boundaries, sample_labels, names.pre
  )
}

# --- Positions ---
start.plot <- -750
end.plot <- 750
pos <- seq(-upstream, downstream - binSize, binSize)

TES.idx <- which(pos == 0)
profile.index <- (TES.idx + start.plot / binSize):(TES.idx + end.plot / binSize)
pos.zoom <- pos[profile.index]

# --- Colors ---
type_colors <- c(convergent = "#D95F02", isolated = "#1B9E77", tandem = "#7570B3")
type_lty   <- c(convergent = 1, isolated = 1, tandem = 1)

# --- Define the 3 comparison panels ---
comparisons <- list(
  list(
    title    = "TKO",
    sample   = "_TKO_Rep_1",
    filename = "TES_compare_TKO_gene_types_gt1000bp"
  ),
  list(
    title    = "TKO - Pol II",
    sample   = "_TKO_Rpb1_Rapamycin_120_min",
    filename = "TES_compare_TKO_Rpb1_gene_types_gt1000bp"
  ),
  list(
    title    = "TKO - Pol II - INO80",
    sample   = "_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep_1",
    filename = "TES_compare_TKO_Rpb1_Ino80_gene_types_gt1000bp"
  )
)

for (comp in comparisons) {

  smp <- comp$sample

  smp_vals_by_type <- lapply(profiles_by_type, function(p) p[[smp]][profile.index])

  all_vals <- unlist(smp_vals_by_type)
  y_low  <- 0
  y_high <- max(all_vals)
  y_pad  <- (y_high - y_low) * 0.08
  y_high <- y_high + y_pad

  render_plot <- function(filename_base, draw_fun, width, height) {
    pdf(paste0(filename_base, ".pdf"), width = width, height = height)
    draw_fun()
    dev.off()
    cat("Saved:", paste0(filename_base, ".pdf"), "\n")

    png(paste0(filename_base, ".png"), width = width * 150, height = height * 150, res = 150)
    draw_fun()
    dev.off()
    cat("Saved:", paste0(filename_base, ".png"), "\n")
  }

  draw_comp <- function() {
    par(mar = c(11, 4.1, 3.1, 2.1), xaxs = "i", yaxs = "i")

    plot(pos.zoom, smp_vals_by_type[[1]],
         type = "n",
         xlab = NA,
         ylab = "Nucleosome Occupancy",
         main = paste0("TES Profile - ", comp$title, " (> 1000 bp)"),
         ylim = c(y_low, y_high),
         xlim = c(start.plot, end.plot),
         axes = TRUE)

    for (tp in names(gene_types)) {
      lines(pos.zoom, smp_vals_by_type[[tp]],
            col = type_colors[tp], lwd = 2.5, lty = type_lty[tp])
    }

    abline(v = 0, lty = 3, col = "grey70")

    usr <- par("usr")
    y_legend <- usr[3] - (usr[4] - usr[3]) * 0.12
    legend(start.plot, y_legend,
           bty = "n",
           legend = names(gene_types),
           col = type_colors,
           lwd = 2.5,
           lty = type_lty,
           xpd = TRUE)

    title(xlab = "distance from TES", line = 3)
  }

  render_plot(comp$filename, draw_comp, 7, 7)
}

cat("Done.\n")