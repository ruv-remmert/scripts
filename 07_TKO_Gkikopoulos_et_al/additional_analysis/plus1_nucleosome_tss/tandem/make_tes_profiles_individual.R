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
# Individual mutant vs WT_Mnase plots
############################################

wt_name <- "WT_Mnase"
wt_vals <- as.numeric(heat.val[wt_name, profile.index])

# get all mutant samples (everything except WT)
mutant_samples <- names[names != wt_name]

for(i in 1:length(mutant_samples)){

    mut <- mutant_samples[i]
    mut_vals <- as.numeric(heat.val[mut, profile.index])

    plot_names <- c(mut, wt_name)

    # Determine y-axis limits
    low  <- min(c(mut_vals, wt_vals))
    high <- max(c(mut_vals, wt_vals))

    comp_name <- paste0(mut, "_vs_", wt_name)

    pdf(paste0(out.prefix, "_", comp_name, ".pdf"),
        width = 14,
        height = 8)

    par(mar = c(5.1, 4.1, 4.1, 13.8))

    # plot mutant
    plot(pos,
         mut_vals,
         type = "l",
         xlab = "+1 nucleosome",
         ylab = "MNase fragment density",
         lwd = 2,
         main = paste(mut, "vs", wt_name),
         col = dark2[i],
         ylim = c(low, high))

    # overlay WT in grey dashed
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
