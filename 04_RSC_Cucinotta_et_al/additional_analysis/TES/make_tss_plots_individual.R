#!/usr/bin/env Rscript

#### PACKAGE LOADING ####
library(RColorBrewer)

######## DATA LOADING ###########
args <- commandArgs(TRUE)

heat.val <- read.delim(file = args[1])
out.prefix <- sub("\\.pdf$", "", args[2])

# Create a subfolder named after the output prefix and redirect all output into it
out.dir <- file.path(dirname(out.prefix), basename(out.prefix))
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
out.prefix <- file.path(out.dir, basename(out.prefix))

names.pre <- gsub("_monoNucs_profile", "", as.character(heat.val$bin.labels)[-1])
row.names(heat.val) <- c("bin", names.pre)
heat.val <- heat.val[order(heat.val[,1]),]

# defining color palette
dark2 <- c(
    RColorBrewer::brewer.pal(8, "Dark2"),
    RColorBrewer::brewer.pal(8, "Set1"),
    RColorBrewer::brewer.pal(8, "Set2")
)

TSS.pos <- which(colnames(heat.val) == "tick")

start.plot <- -1250
end.plot <- 1500
pos <- seq(start.plot, end.plot, 10)

profile.index <- (TSS.pos + start.plot/10):(TSS.pos + end.plot/10)

# Remove pooled
if(sum(names.pre %in% "pooled") > 0){
    names <- names.pre[!(names.pre %in% "pooled")]
} else {
    names <- names.pre
}

############################################
# Pairwise crl vs IAA comparison plots
############################################

# Find all crl samples
crl_samples <- names[grepl("crl", names)]

for(crl in crl_samples){

    # Create matching IAA sample name
    iaa <- gsub("crl", "IAA", crl)

    # Skip if no matching IAA sample exists
    if(!(iaa %in% names)){
        next
    }

    plot_names <- c(crl, iaa)

    # Determine y-axis limits
    low <- min(sapply(plot_names, function(n)
        min(heat.val[n, profile.index])))

    high <- max(sapply(plot_names, function(n)
        max(heat.val[n, profile.index])))

    # Output name
    comp_name <- paste0(crl, "_vs_", iaa)

    pdf(paste0(out.prefix, "_", comp_name, ".pdf"),
        width = 14,
        height = 8)

    par(mar = c(5.1, 4.1, 4.1, 13.8))

    for(k in 1:length(plot_names)){

        if(k == 1){

            plot(pos,
                 heat.val[plot_names[k], profile.index],
                 type = "l",
                 xlab = "distance from TSS",
                 ylab = "MNase fragment density",
                 lwd = 2,
                 main = paste0(crl, " vs ", iaa),
                 col = dark2[k],
                 ylim = c(low, high))

        } else {

            lines(pos,
                  heat.val[plot_names[k], profile.index],
                  col = dark2[k],
                  lwd = 2)
        }
    }

    legend(1600,
           high + (high-low)/10,
           bty = "n",
           legend = plot_names,
           col = dark2[1:length(plot_names)],
           lwd = 2,
           xpd = TRUE)

    dev.off()
}
