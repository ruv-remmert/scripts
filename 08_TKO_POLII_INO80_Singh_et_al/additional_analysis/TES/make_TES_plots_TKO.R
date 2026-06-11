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
# Combined plot with ALL samples together
############################################

low  <- min(sapply(names, function(n) min(heat.val[n, profile.index])))
high <- max(sapply(names, function(n) max(heat.val[n, profile.index])))

pdf(paste0(out.prefix, "_ALL.pdf"), width = 14, height = 8)
par(mar = c(5.1, 4.1, 4.1, 13.8))

for(k in 1:length(names)){

    if(k == 1){

        plot(pos,
             heat.val[names[k], profile.index],
             type = "l",
             xlab = "distance from TES",
             ylab = "MNase fragment density",
             lwd = 2,
             main = "TES profile: all samples",
             col = dark2[k],
             ylim = c(low, high))

    } else {

        lines(pos,
              heat.val[names[k], profile.index],
              col = dark2[k],
              lwd = 2)
    }
}

legend(1600,
       high + (high-low)/10,
       bty = "n",
       legend = names,
       col = dark2[1:length(names)],
       lwd = 2,
       xpd = TRUE)

dev.off()
cat("Saved:", paste0(out.prefix, "_ALL.pdf"), "\n")
