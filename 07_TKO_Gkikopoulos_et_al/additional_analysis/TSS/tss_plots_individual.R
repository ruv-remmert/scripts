#!/usr/bin/env Rscript

#### PACKAGE LOADING ####
library(RColorBrewer)
library(jsonlite)

######## DATA LOADING ###########
args <- commandArgs(TRUE)
infile <- args[1]

# --- Parse metadata from header line (handles .gz and plain text) ---
con <- if (grepl("\\.gz$", infile)) gzcon(file(infile, "rb")) else file(infile, "r")
header <- readLines(con, n = 1)
close(con)

meta <- fromJSON(sub("^@", "", header))

upstream          <- meta$upstream[1]
downstream        <- meta$downstream[1]
binSize           <- meta$`bin size`[1]
sample_labels     <- meta$sample_labels
sample_boundaries <- meta$sample_boundaries

# strip "_monoNucs_profile" and "_1.fastq.gz" suffixes to get clean sample names
names.pre <- gsub("_monoNucs_profile", "", sample_labels)
names.pre <- gsub("_1\\.fastq\\.gz", "", names.pre)

# --- Read signal data (skip header line, handles .gz and plain text) ---
mat <- read.table(
  if (grepl("\\.gz$", infile)) gzfile(infile) else infile,
  skip = 1, header = FALSE, sep = "\t", fill = TRUE
)
signal <- mat[, 7:ncol(mat)]

# --- Average signal per sample across all genes ---
avg_profiles <- list()
for (i in seq_len(length(sample_labels))) {
  col_start <- sample_boundaries[i] + 1
  col_end   <- sample_boundaries[i + 1]
  samp_mat  <- signal[, col_start:col_end, drop = FALSE]
  avg_profiles[[names.pre[i]]] <- colMeans(samp_mat, na.rm = TRUE)
}

# --- Position vector ---
pos <- seq(-upstream, downstream - binSize, binSize)

# --- Color palette ---
dark2 <- c(RColorBrewer::brewer.pal(8, "Dark2"),
           RColorBrewer::brewer.pal(8, "Set1"),
           RColorBrewer::brewer.pal(8, "Set2"))

# --- WT reference ---
wt_name <- "WT_Mnase"
wt_vals <- avg_profiles[[wt_name]]

# ---- Combined plot: all samples overlaid ----
all_vals <- unlist(avg_profiles)
y_pad  <- (max(all_vals) - min(all_vals)) * 0.1
y_low  <- min(all_vals) - y_pad
y_high <- max(all_vals) + y_pad

pdf("profile_monoNucs.pdf", width = 14, height = 7)
par(mar = c(6, 6, 5, 16), cex.axis = 1.3, cex.lab = 1.5, cex.main = 1.7)

for (i in seq_along(names.pre)) {
  if (i == 1) {
    plot(pos, avg_profiles[[names.pre[i]]],
         type = "l", xlab = "Distance from TSS (bp)",
         ylab = "MNase fragment density", lwd = 2.5,
         main = "TSS profile - all samples",
         col = dark2[i], ylim = c(y_low, y_high))
    abline(v = 0, lty = 2, col = "grey50")
  } else {
    lines(pos, avg_profiles[[names.pre[i]]], col = dark2[i], lwd = 2.5)
  }
}

legend(downstream + 50, y_high, bty = "n", legend = names.pre,
       col = dark2[seq_along(names.pre)], lwd = 2.5, cex = 1.2, xpd = TRUE)
dev.off()
cat("Saved: profile_monoNucs.pdf\n")

# ---- Individual plots: each mutant vs WT_Mnase ----
for (i in seq_along(names.pre)) {

  # skip WT itself
  if (names.pre[i] == wt_name) next

  curr_vals <- avg_profiles[[names.pre[i]]]

  y_pad  <- (max(c(curr_vals, wt_vals)) - min(c(curr_vals, wt_vals))) * 0.1
  y_low  <- min(c(curr_vals, wt_vals)) - y_pad
  y_high <- max(c(curr_vals, wt_vals)) + y_pad

  outfile <- paste0("profile_monoNucs_", names.pre[i], "_vs_", wt_name, ".pdf")
  pdf(outfile, width = 10, height = 6)
  par(mar = c(6, 6, 5, 4), cex.axis = 1.3, cex.lab = 1.5, cex.main = 1.7)

  plot(pos, curr_vals,
       type = "l",
       xlab = "Distance from TSS (bp)",
       ylab = "MNase fragment density",
       lwd  = 2.5,
       main = paste("TSS profile -", names.pre[i], "vs", wt_name),
       col  = dark2[i],
       ylim = c(y_low, y_high))

  lines(pos, wt_vals, col = "grey40", lwd = 2.5, lty = 2)
  abline(v = 0, lty = 3, col = "grey70")

  legend("topright",
         legend = c(names.pre[i], wt_name),
         col    = c(dark2[i], "grey40"),
         lwd    = 2.5,
         lty    = c(1, 2),
         bty    = "n",
         cex    = 1.2)

  dev.off()
  cat("Saved:", outfile, "\n")
}
