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

mat <- mat[keep_genes, , drop = FALSE]
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

fill_color <- "grey80"
tko_color <- "#A82649"
rpb1_color <- "#E59B41"
ino80_color <- "#2E7DD5"

wt_name    <- "_Wild_type_Rep_1"
tko_name   <- "_TKO_Rep_1"
rpb1_name  <- "_TKO_Rpb1_Rapamycin_120_min"
ino80_name <- "_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep_1"

wt_vals    <- avg_profiles[[wt_name]][profile.index]
tko_vals   <- avg_profiles[[tko_name]][profile.index]
rpb1_vals  <- avg_profiles[[rpb1_name]][profile.index]
ino80_vals <- avg_profiles[[ino80_name]][profile.index]

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

############################################
# Plot: WT (grey filled) vs TKO vs TKO - Pol II vs TKO - Pol II - INO80
############################################

y_low  <- 0
y_high <- max(c(wt_vals, tko_vals, rpb1_vals, ino80_vals))
y_pad  <- (y_high - y_low) * 0.08
y_high <- y_high + y_pad

draw_special <- function() {
  par(mar = c(11, 4.1, 3.1, 2.1), xaxs = "i", yaxs = "i")

  plot(pos.zoom, wt_vals,
       type = "n",
       xlab = NA,
       ylab = "Nucleosome Occupancy",
       main = "TES Profile (> 1000 bp)",
       ylim = c(y_low, y_high),
       xlim = c(start.plot, end.plot),
       axes = TRUE)

  polygon(c(pos.zoom, rev(pos.zoom)), c(wt_vals, rep(0, length(pos.zoom))),
          col = fill_color, border = NA)

  lines(pos.zoom, wt_vals, col = fill_color, lwd = 2)

  lines(pos.zoom, tko_vals, col = tko_color, lwd = 2.5)

  lines(pos.zoom, rpb1_vals, col = rpb1_color, lwd = 2.5)

  lines(pos.zoom, ino80_vals, col = ino80_color, lwd = 2.5)

  abline(v = 0, lty = 3, col = "grey70")

  usr <- par("usr")
  y_legend <- usr[3] - (usr[4] - usr[3]) * 0.12
  legend(start.plot, y_legend,
         bty = "n",
         legend = c("WT", "TKO", "TKO - Pol II", "TKO - Pol II - INO80"),
         col = c(fill_color, tko_color, rpb1_color, ino80_color),
         lwd = c(2, 2.5, 2.5, 2.5),
         fill = c(fill_color, NA, NA, NA),
         border = c(fill_color, NA, NA, NA),
         xpd = TRUE)

  title(xlab = "distance from TES", line = 3)
}

render_plot(file.path(out.dir, "TES_special_comparison_gt1000bp"), draw_special, 7, 7)