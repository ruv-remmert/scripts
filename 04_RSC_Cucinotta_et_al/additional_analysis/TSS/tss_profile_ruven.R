#!/usr/bin/env Rscript

#### PACKAGE LOADING ####
library(RColorBrewer)
library(jsonlite)

######## DATA LOADING ###########
args <- commandArgs(TRUE)
infile <- args[1]

# --- Parse metadata from header line ---
con <- file(infile, "r")
header <- readLines(con, n = 1)
close(con)

meta <- fromJSON(sub("^@", "", header))

upstream          <- meta$upstream[1]
downstream        <- meta$downstream[1]
binSize           <- meta$`bin size`[1]
sample_labels     <- meta$sample_labels
sample_boundaries <- meta$sample_boundaries

names.pre <- gsub("_monoNucs_profile", "", sample_labels)

# --- Read signal data (skip header line) ---
mat <- read.table(infile, skip = 1, header = FALSE, sep = "\t", fill = TRUE)
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

# --- Matched pairs: IAA treatment vs its corresponding control ---
pairs <- list(
  c(treatment = "MNase_IAA_Q_50U",   control = "MNase_crl_Q_50U"),
  c(treatment = "MNase_IAA_10m_50U", control = "MNase_crl_10m_50U"),
  c(treatment = "MNase_IAA_Q_5U",    control = "MNase_crl_Q_5U"),
  c(treatment = "MNase_IAA_10m_5U",  control = "MNase_crl_10m_5U")
)

# ---- Individual plots: each IAA sample vs its matched control ----
for (pair in pairs) {
  treatment <- pair["treatment"]
  control   <- pair["control"]

  trt_vals <- avg_profiles[[treatment]]
  crl_vals <- avg_profiles[[control]]

  y_pad  <- (max(c(trt_vals, crl_vals)) - min(c(trt_vals, crl_vals))) * 0.1
  y_low  <- min(c(trt_vals, crl_vals)) - y_pad
  y_high <- max(c(trt_vals, crl_vals)) + y_pad

  outfile <- paste0("profile_monoNucs_", treatment, "_vs_", control, ".pdf")
  pdf(outfile, width = 10, height = 6)
  par(mar = c(6, 6, 5, 4), cex.axis = 1.3, cex.lab = 1.5, cex.main = 1.7)

  plot(pos, trt_vals, type = "l",
       xlab = "Distance from TSS (bp)", ylab = "MNase fragment density",
       lwd = 2.5, main = paste("TSS profile -", treatment, "vs", control),
       col = "#2E7D5B", ylim = c(y_low, y_high))

  lines(pos, crl_vals, col = "grey40", lwd = 2.5, lty = 2)
  abline(v = 0, lty = 3, col = "grey70")

  legend("topright", legend = c(treatment, control),
         col = c("#2E7D5B", "grey40"), lwd = 2.5, lty = c(1, 2), bty = "n", cex = 1.2)

  dev.off()
  cat("Saved:", outfile, "\n")
}
