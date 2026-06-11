#!/usr/bin/env Rscript

#### PACKAGE LOADING ####
library(RColorBrewer)
library(jsonlite)

######## DATA LOADING ###########
args <- commandArgs(TRUE)
infile <- args[1]

# --- Parse metadata from computeMatrix header ---
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

# clean names
names.pre <- gsub("_monoNucs_profile", "", sample_labels)
names.pre <- gsub("_1\\.fastq\\.gz", "", names.pre)

# --- Read computeMatrix matrix.gz ---
mat <- read.table(
  if (grepl("\\.gz$", infile)) gzfile(infile) else infile,
  skip = 1,
  header = FALSE,
  sep = "\t",
  fill = TRUE
)

# --------------------------------------------------
# Keep only genes >1000 bp
# (Gene length from GTF since TES bed has 1-bp coords)
# --------------------------------------------------

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
keep_genes <- mat_gene %in% gtf_length$gene[gtf_length$length > 1000]

cat("Genes before filtering:", nrow(mat), "\n")
cat("Genes >1000 bp:", sum(keep_genes), "\n")
cat("Genes removed:", sum(!keep_genes), "\n")

mat <- mat[keep_genes, , drop = FALSE]

signal <- mat[,7:ncol(mat)]

# --------------------------------------------------
# Average signal per sample
# --------------------------------------------------

baseline_bins <- 1:50

avg_profiles <- list()

for (i in seq_len(length(sample_labels))) {

  col_start <- sample_boundaries[i] + 1
  col_end   <- sample_boundaries[i + 1]

  samp_mat <- signal[, col_start:col_end, drop = FALSE]
  vals <- colMeans(samp_mat, na.rm = TRUE)

  baseline <- mean(vals[baseline_bins], na.rm = TRUE)

  avg_profiles[[names.pre[i]]] <- vals / baseline
}

# --- TES positions ---
pos <- seq(-upstream, downstream - binSize, binSize)

# --- Colors ---
dark2 <- c(
  brewer.pal(8,"Dark2"),
  brewer.pal(8,"Set1"),
  brewer.pal(8,"Set2")
)

# --- WT reference ---
wt_name <- "_Wild_type_Rep_1"
wt_vals <- avg_profiles[[wt_name]]

# --------------------------------------------------
# Combined TES plot
# --------------------------------------------------

all_vals <- unlist(avg_profiles)

y_pad  <- (max(all_vals)-min(all_vals))*0.1
y_low  <- min(all_vals)-y_pad
y_high <- max(all_vals)+y_pad

pdf("TES_profile_gt1000bp.pdf", width=14, height=7)

par(
  mar=c(6,6,5,16),
  cex.axis=1.3,
  cex.lab=1.5,
  cex.main=1.7
)

for(i in seq_along(names.pre)) {

  if(i==1){

    plot(
      pos,
      avg_profiles[[names.pre[i]]],
      type="l",
      xlab="Distance from TES (bp)",
      ylab="MNase fragment density",
      lwd=2.5,
      main="TES profile - all samples (genes >1000 bp)",
      col=dark2[i],
      ylim=c(y_low,y_high)
    )

    abline(v=0,lty=2,col="grey50")

  } else {

    lines(
      pos,
      avg_profiles[[names.pre[i]]],
      col=dark2[i],
      lwd=2.5
    )
  }
}

legend(
  downstream+50,
  y_high,
  bty="n",
  legend=names.pre,
  col=dark2[seq_along(names.pre)],
  lwd=2.5,
  cex=1.2,
  xpd=TRUE
)

dev.off()

cat("Saved: TES_profile_gt1000bp.pdf\n")

# --------------------------------------------------
# Individual plots vs WT
# --------------------------------------------------

for(i in seq_along(names.pre)) {

  if(names.pre[i] == wt_name) next

  curr_vals <- avg_profiles[[names.pre[i]]]

  y_pad  <- (max(c(curr_vals,wt_vals)) -
             min(c(curr_vals,wt_vals))) * 0.1
  y_low  <- min(c(curr_vals,wt_vals)) - y_pad
  y_high <- max(c(curr_vals,wt_vals)) + y_pad

  outfile <- paste0(
    "TES_profile_gt1000bp_",
    names.pre[i],
    "_vs_",
    wt_name,
    ".pdf"
  )

  pdf(outfile, width=10, height=6)

  par(
    mar=c(6,6,5,4),
    cex.axis=1.3,
    cex.lab=1.5,
    cex.main=1.7
  )

  plot(
    pos,
    curr_vals,
    type="l",
    xlab="Distance from TES (bp)",
    ylab="MNase fragment density",
    lwd=2.5,
    main=paste(
      "TES profile -",
      names.pre[i],
      "vs",
      wt_name,
      "(genes >1000 bp)"
    ),
    col=dark2[i],
    ylim=c(y_low,y_high)
  )

  lines(
    pos,
    wt_vals,
    col="grey40",
    lwd=2.5,
    lty=2
  )

  abline(v=0,lty=3,col="grey70")

  legend(
    "topright",
    legend=c(names.pre[i],wt_name),
    col=c(dark2[i],"grey40"),
    lwd=2.5,
    lty=c(1,2),
    bty="n",
    cex=1.2
  )

  dev.off()

  cat("Saved:", outfile, "\n")
}

# ---- Special combined plot: WT + TKO + TKO_Rpb1_Rapamycin + TKO_Rpb1_Ino80_Rapamycin ----
special_samples <- c(
  "_Wild_type_Rep_1",
  "_TKO_Rep_1",
  "_TKO_Rpb1_Rapamycin_120_min",
  "_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep_1"
)

special_colors <- c("#404040", "#A82649", "#E59B41", "#2E7D5B")
special_labels <- c("Wild type", "TKO", "TKO + Rpb1 depletion", "TKO + Rpb1 + Ino80 depletion")

missing <- special_samples[!special_samples %in% names.pre]
if (length(missing) > 0) {
  cat("WARNING: missing special samples:", paste(missing, collapse = ", "), "\n")
  cat("Available names:", paste(names.pre, collapse = ", "), "\n")
} else {

  special_vals <- lapply(special_samples, function(s) avg_profiles[[s]])

  all_special <- unlist(special_vals)
  y_pad  <- (max(all_special) - min(all_special)) * 0.1
  y_low  <- min(all_special) - y_pad
  y_high <- max(all_special) + y_pad

  pdf("TES_profile_rapamycin_summary_gt1000bp.pdf", width = 12, height = 6)
  par(mar = c(6, 6, 5, 16), cex.axis = 1.3, cex.lab = 1.5, cex.main = 1.7)

  for (i in seq_along(special_samples)) {
    if (i == 1) {
      plot(pos, special_vals[[i]],
           type = "l",
           xlab = "Distance from TES (bp)",
           ylab = "MNase fragment density",
           lwd  = 2.5,
           main = "TES profile - Rapamycin depletion summary (genes > 1000 bp)",
           col  = special_colors[i],
           ylim = c(y_low, y_high))
      abline(v = 0, lty = 2, col = "grey70")
    } else {
      lines(pos, special_vals[[i]], col = special_colors[i], lwd = 2.5)
    }
  }

  legend(downstream + 50, y_high,
         bty    = "n",
         legend = special_labels,
         col    = special_colors,
         lwd    = 2.5,
         cex    = 1.2,
         xpd    = TRUE)

  dev.off()
  cat("Saved: TES_profile_rapamycin_summary_gt1000bp.pdf\n")
}
