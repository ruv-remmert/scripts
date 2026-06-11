#!/usr/bin/env Rscript

#### PACKAGE LOADING ####
library(RColorBrewer)

######## DATA LOADING ###########
args   <- commandArgs(TRUE)


heat.val<-read.delim(file = args[1])
out.prefix <- sub("\\.pdf$", "", args[2])

# Create a subfolder named after the output prefix and redirect all output into it
out.dir    <- file.path(dirname(out.prefix), basename(out.prefix))
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
out.prefix <- file.path(out.dir, basename(out.prefix))

names.pre <- gsub("_monoNucs_profile","",as.character(heat.val$bin.labels)[-1])
row.names(heat.val)<-c("bin",names.pre)
heat.val<-heat.val[order(heat.val[,1]),]

#defining color palette
dark2 <- c(RColorBrewer::brewer.pal(8, "Dark2"),RColorBrewer::brewer.pal(8, "Set1"),RColorBrewer::brewer.pal(8, "Set2"))



TSS.pos<-which(colnames(heat.val)=="tick")

start.plot<-c(-1250)
end.plot<-c(1500)
pos<-seq(start.plot,end.plot,10)

profile.index<- (TSS.pos+start.plot/10):(TSS.pos+end.plot/10)


# Remove pooled
if(sum(names.pre %in% "pooled")>0){
    names <- names.pre[!(names.pre %in% "pooled")]
} else {
    names <- names.pre
}

# Overview plot: all samples on one panel
low  <- min(sapply(names, function(n) min(heat.val[n, profile.index])))
high <- max(sapply(names, function(n) max(heat.val[n, profile.index])))

pdf(paste0(out.prefix, ".pdf"), width = 14, height=8)
par(mar=c(5.1, 4.1, 4.1, 13.8))
for(k in 1:length(names)){
    if(k==1){
        plot(pos,
             heat.val[names[k], profile.index],
             type="l",
             xlab="distance from TSS", ylab = "MNase fragment density",
             lwd=2, main = paste0("TSS profile: all samples (", out.prefix, ")"), col=dark2[k], ylim=c(low,high))
    } else {
        lines(pos,
              heat.val[names[k], profile.index],
              col=dark2[k], lwd=2)
    }
}
legend(1600, high+(high-low)/10, bty="n",
       legend=names,
       col=dark2[1:length(names)], lwd=2,
       xpd=TRUE)
dev.off()

# Identify WT_A
wt_name <- names[grepl("WT_A$", names)]

# Identify all mutations (excluding WT)
mutations <- unique(gsub("_[AB]$", "", names[!grepl("WT", names)]))

for(mut in mutations){
    # Find replicates for this mutation
    mut_A <- paste0(mut, "_A")
    mut_B <- paste0(mut, "_B")
    plot_names <- c(wt_name, mut_A, mut_B)
    plot_names <- plot_names[plot_names %in% names] # Only keep those present
    if(length(plot_names) < 2) next # Skip if not enough samples

    # Get min and max for y-axis
    for(j in 1:length(plot_names)){
        curr_min <- min(heat.val[plot_names[j], profile.index ])
        curr_max <- max(heat.val[plot_names[j], profile.index ])
        if(j==1){
            low <- curr_min
            high <- curr_max
        } else {
            low <- min(curr_min,low)
            high <- max(curr_max,high)
        }
    }

    pdf(paste0(out.prefix, "_", mut, ".pdf"), width = 14, height=8)
    par(mar=c(5.1, 4.1, 4.1, 13.8))
    for(k in 1:length(plot_names)){
        if(k==1){
            plot(pos,
                 heat.val[plot_names[k],profile.index ],
                 type="l",
                 xlab="distance from TSS", ylab = "MNase fragment density",
                 lwd=2, main=paste0("TSS profile: ", mut, " vs WT_A"), col=dark2[k], ylim=c(low,high))
        } else {
            lines(pos,
                  heat.val[plot_names[k], profile.index],
                  col=dark2[k],lwd=2)
        }
    }
    legend(1600,high+(high-low)/10, bty="n",
           legend=plot_names,
           col=dark2[1:length(plot_names)], lwd=2,
           xpd=TRUE)
    dev.off()
}
