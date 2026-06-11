#!/usr/bin/env Rscript

library(jsonlite)

######## DATA LOADING ###########
args <- commandArgs(TRUE)
infile <- args[1]
out.dir <- args[2]
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

con <- if (grepl("\\.gz$", infile)) gzfile(infile, "rt") else file(infile, "r")
header <- readLines(con, n = 1, warn = FALSE)
close(con)

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
  fill = TRUE
)

cat("Total genes:", nrow(mat), "\n")

signal <- mat[, 7:ncol(mat)]

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

start.plot <- -750
end.plot <- 750
pos <- seq(-upstream, downstream - binSize, binSize)

TES.idx <- which(pos == 0)
profile.index <- (TES.idx + start.plot / binSize):(TES.idx + end.plot / binSize)
pos.zoom <- pos[profile.index]

fill_color       <- "grey80"
tko_color        <- "#CC0000"
ino80_tko_color  <- "#7B2D8E"
ino80_wt_color   <- "#2CA02C"

wt_name      <- "Wild_type_Rep_2"
tko_name     <- "TKO_Rep_1"
ino80_tko    <- "TKO_Ino80_Rapamycin_90_min_Rep_1"
ino80_wt     <- "WT_Ino80_Rapamycin_90_min_Rep_1"

wt_vals        <- avg_profiles[[wt_name]][profile.index]
tko_vals       <- avg_profiles[[tko_name]][profile.index]
ino80_tko_vals <- avg_profiles[[ino80_tko]][profile.index]
ino80_wt_vals  <- avg_profiles[[ino80_wt]][profile.index]

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

y_low_1  <- 0
y_high_1 <- max(c(wt_vals, tko_vals, ino80_tko_vals))
y_pad_1  <- (y_high_1 - y_low_1) * 0.08
y_high_1 <- y_high_1 + y_pad_1

draw_special <- function() {
  par(mar = c(11, 4.1, 3.1, 2.1), xaxs = "i", yaxs = "i")

  plot(pos.zoom, wt_vals,
       type = "n",
       xlab = NA,
       ylab = "Nucleosome Occupancy",
       main = "TES Profile (all genes)",
       ylim = c(y_low_1, y_high_1),
       xlim = c(start.plot, end.plot),
       axes = TRUE)

  polygon(c(pos.zoom, rev(pos.zoom)), c(wt_vals, rep(0, length(pos.zoom))),
          col = fill_color, border = NA)

  lines(pos.zoom, wt_vals, col = fill_color, lwd = 2)

  lines(pos.zoom, tko_vals, col = tko_color, lwd = 2.5)

  lines(pos.zoom, ino80_tko_vals, col = ino80_tko_color, lwd = 2.5)

  abline(v = 0, lty = 3, col = "grey70")

  usr <- par("usr")
  y_legend <- usr[3] - (usr[4] - usr[3]) * 0.12
  legend(start.plot, y_legend,
         bty = "n",
         legend = c("WT", "TKO", "TKO + Ino80 depletion"),
         col = c(fill_color, tko_color, ino80_tko_color),
         lwd = c(2, 2.5, 2.5),
         fill = c(fill_color, NA, NA),
         border = c(fill_color, NA, NA),
         xpd = TRUE)

  title(xlab = "distance from TES", line = 3)
}

render_plot(file.path(out.dir, "TES_special_comparison_all_genes"), draw_special, 7, 7)

y_low_2  <- 0
y_high_2 <- max(c(wt_vals, ino80_wt_vals))
y_pad_2  <- (y_high_2 - y_low_2) * 0.08
y_high_2 <- y_high_2 + y_pad_2

draw_wt_ino80 <- function() {
  par(mar = c(11, 4.1, 3.1, 2.1), xaxs = "i", yaxs = "i")

  plot(pos.zoom, wt_vals,
       type = "n",
       xlab = NA,
       ylab = "Nucleosome Occupancy",
       main = "TES Profile - WT + Ino80 depletion (all genes)",
       ylim = c(y_low_2, y_high_2),
       xlim = c(start.plot, end.plot),
       axes = TRUE)

  polygon(c(pos.zoom, rev(pos.zoom)), c(wt_vals, rep(0, length(pos.zoom))),
          col = fill_color, border = NA)

  lines(pos.zoom, wt_vals, col = fill_color, lwd = 2)

  lines(pos.zoom, ino80_wt_vals, col = ino80_wt_color, lwd = 2.5)

  abline(v = 0, lty = 3, col = "grey70")

  usr <- par("usr")
  y_legend <- usr[3] - (usr[4] - usr[3]) * 0.12
  legend(start.plot, y_legend,
         bty = "n",
         legend = c("WT", "WT + Ino80 depletion"),
         col = c(fill_color, ino80_wt_color),
         lwd = c(2, 2.5),
         fill = c(fill_color, NA),
         border = c(fill_color, NA),
         xpd = TRUE)

  title(xlab = "distance from TES", line = 3)
}

render_plot(file.path(out.dir, "TES_WT_vs_WT_Ino80_all_genes"), draw_wt_ino80, 7, 7)

y_low_3  <- 0
y_high_3 <- max(c(wt_vals, tko_vals, ino80_tko_vals, ino80_wt_vals))
y_pad_3  <- (y_high_3 - y_low_3) * 0.08
y_high_3 <- y_high_3 + y_pad_3

draw_combined <- function() {
  par(mar = c(11, 4.1, 3.1, 2.1), xaxs = "i", yaxs = "i")

  plot(pos.zoom, wt_vals,
       type = "n",
       xlab = NA,
       ylab = "Nucleosome Occupancy",
       main = "TES Profile - Combined (all genes)",
       ylim = c(y_low_3, y_high_3),
       xlim = c(start.plot, end.plot),
       axes = TRUE)

  polygon(c(pos.zoom, rev(pos.zoom)), c(wt_vals, rep(0, length(pos.zoom))),
          col = fill_color, border = NA)

  lines(pos.zoom, wt_vals, col = fill_color, lwd = 2)

  lines(pos.zoom, tko_vals, col = tko_color, lwd = 2.5)

  lines(pos.zoom, ino80_tko_vals, col = ino80_tko_color, lwd = 2.5)

  lines(pos.zoom, ino80_wt_vals, col = ino80_wt_color, lwd = 2.5)

  abline(v = 0, lty = 3, col = "grey70")

  usr <- par("usr")
  y_legend <- usr[3] - (usr[4] - usr[3]) * 0.12
  legend(start.plot, y_legend,
         bty = "n",
         legend = c("WT", "TKO", "TKO + Ino80 depletion", "WT + Ino80 depletion"),
         col = c(fill_color, tko_color, ino80_tko_color, ino80_wt_color),
         lwd = c(2, 2.5, 2.5, 2.5),
         fill = c(fill_color, NA, NA, NA),
         border = c(fill_color, NA, NA, NA),
         xpd = TRUE)

  title(xlab = "distance from TES", line = 3)
}

render_plot(file.path(out.dir, "TES_combined_all_genes"), draw_combined, 7, 7)
