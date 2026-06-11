#!/usr/bin/env Rscript

#### PACKAGE LOADING ####
library(RColorBrewer)

######## DATA LOADING ###########
args <- commandArgs(TRUE)

heat.val <- read.delim(file = args[1])
out.prefix <- sub("\\.pdf$", "", args[2])

# Create a subfolder named after the output prefix
out.dir <- file.path(dirname(out.prefix), basename(out.prefix))
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
out.prefix <- file.path(out.dir, basename(out.prefix))

names.pre <- gsub("_monoNucs_profile", "", as.character(heat.val$bin.labels)[-1])
names.pre <- gsub("_1\\.fastq\\.gz", "", names.pre)
row.names(heat.val) <- c("bin", names.pre)
heat.val <- heat.val[order(heat.val[,1]),]

# defining color palette
dark2 <- c(
    RColorBrewer::brewer.pal(8, "Dark2"),
    RColorBrewer::brewer.pal(8, "Set1"),
    RColorBrewer::brewer.pal(8, "Set2")
)

TES.pos <- which(colnames(heat.val) == "tick")

start.plot <- -1250
end.plot <- 1500
pos <- seq(start.plot, end.plot, 10)

profile.index <- (TES.pos + start.plot/10):(TES.pos + end.plot/10)

# Remove pooled
if(sum(names.pre %in% "pooled") > 0){
    names <- names.pre[!(names.pre %in% "pooled")]
} else {
    names <- names.pre
}

############################################
# Individual plots: each sample vs Wild_type_Rep_1
############################################

wt_name <- "_Wild_type_Rep_1"
wt_vals <- as.numeric(heat.val[wt_name, profile.index])

# all samples except WT
mutant_samples <- names[names != wt_name]

for(i in 1:length(mutant_samples)){

    mut      <- mutant_samples[i]
    mut_vals <- as.numeric(heat.val[mut, profile.index])

    plot_names <- c(mut, wt_name)

    low  <- min(c(mut_vals, wt_vals))
    high <- max(c(mut_vals, wt_vals))

    comp_name <- paste0(mut, "_vs_", wt_name)

    pdf(paste0(out.prefix, "_", comp_name, ".pdf"),
        width = 14, height = 8)
    par(mar = c(5.1, 4.1, 4.1, 13.8))

    plot(pos, mut_vals,
         type = "l",
         xlab = "distance from TSS",
         ylab = "MNase fragment density",
         lwd = 2,
         main = paste(mut, "vs", wt_name),
         col = dark2[i],
         ylim = c(low, high))

    lines(pos, wt_vals, col = "grey40", lwd = 2, lty = 2)
    abline(v = 0, lty = 3, col = "grey70")

    legend(1600,
           high + (high - low) / 10,
           bty = "n",
           legend = plot_names,
           col = c(dark2[i], "grey40"),
           lwd = 2,
           lty = c(1, 2),
           xpd = TRUE)

    dev.off()
    cat("Saved:", paste0(out.prefix, "_", comp_name, ".pdf"), "\n")
}

############################################
# Special summary plot: WT + TKO + Rpb1 + Rpb1+Ino80 (rapamycin, rep 1)
############################################

special_samples <- c(
    "_Wild_type_Rep_1",
    "_TKO_Rep_1",
    "_TKO_Rpb1_Rapamycin_120_min",
    "_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep_1"
)

special_colors <- c("#404040", "#A82649", "#E59B41", "#2E7D5B")
special_labels <- c("Wild type", "TKO", "TKO + Rpb1 depletion", "TKO + Rpb1 + Ino80 depletion")

missing <- special_samples[!special_samples %in% names]
if (length(missing) > 0) {
    cat("WARNING: missing special samples:", paste(missing, collapse = ", "), "\n")
    cat("Available names:", paste(names, collapse = ", "), "\n")
} else {

    special_vals <- lapply(special_samples, function(s) as.numeric(heat.val[s, profile.index]))

    all_special <- unlist(special_vals)
    y_pad  <- (max(all_special) - min(all_special)) * 0.1
    y_low  <- min(all_special) - y_pad
    y_high <- max(all_special) + y_pad

    pdf(paste0(out.prefix, "_rapamycin_summary.pdf"), width = 12, height = 6)
    par(mar = c(5.1, 4.1, 4.1, 13.8))

    for (i in seq_along(special_samples)) {
        if (i == 1) {
            plot(pos, special_vals[[i]],
                 type = "l",
                 xlab = "distance from TSS",
                 ylab = "MNase fragment density",
                 lwd  = 2.5,
                 main = "+1 TSS profile - Rapamycin depletion summary",
                 col  = special_colors[i],
                 ylim = c(y_low, y_high))
            abline(v = 0, lty = 2, col = "grey70")
        } else {
            lines(pos, special_vals[[i]], col = special_colors[i], lwd = 2.5)
        }
    }

    legend(1600, y_high,
           bty    = "n",
           legend = special_labels,
           col    = special_colors,
           lwd    = 2.5,
           cex    = 1.2,
           xpd    = TRUE)

    dev.off()
    cat("Saved:", paste0(out.prefix, "_rapamycin_summary.pdf"), "\n")
}
