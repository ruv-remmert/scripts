#!/usr/bin/env Rscript

#### PACKAGE LOADING ####
library(RColorBrewer)

######## DATA LOADING ###########
args   <- commandArgs(TRUE)

heat.val<-read.delim(file = args[1])
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

##remove pooled
if(sum(names.pre %in% "pooled")>0){
    names<-names.pre[!(names.pre %in% "pooled")]
} else {
    names<-names.pre
}

# wild type reference
wt_name <- "WT_A"
wt_vals <- as.numeric(heat.val[wt_name, profile.index])

# create one plot per dataset
for(i in 1:length(names)){

    # skip making a standalone plot for WT itself
    if(names[i] == wt_name) next

    curr_vals <- as.numeric(heat.val[names[i], profile.index])

    # y axis scaled to both lines
    y_pad <- (max(c(curr_vals, wt_vals)) - min(c(curr_vals, wt_vals))) * 0.1
    y_low  <- min(c(curr_vals, wt_vals)) - y_pad
    y_high <- max(c(curr_vals, wt_vals)) + y_pad

    outfile <- paste0("profile_monoNucs_", names[i], ".pdf")
    pdf(outfile, width = 10, height = 6)
    par(mar = c(6, 6, 5, 4),
        cex.axis = 1.3,
        cex.lab  = 1.5,
        cex.main = 1.7)

    # plot mutant line
    plot(pos,
         curr_vals,
         type = "l",
         xlab = "Distance from TSS (bp)",
         ylab = "MNase fragment density",
         lwd  = 2.5,
         main = paste("TSS profile —", names[i], "vs WT_A"),
         col  = dark2[i],
         ylim = c(y_low, y_high))

    # overlay WT in grey
    lines(pos, wt_vals, col = "grey40", lwd = 2.5, lty = 2)

    abline(v = 0, lty = 3, col = "grey70")  # vertical line at TSS

    legend("topright",
           legend = c(names[i], "WT_A"),
           col    = c(dark2[i], "grey40"),
           lwd    = 2.5,
           lty    = c(1, 2),
           bty    = "n",
           cex    = 1.2)

    dev.off()
    cat("Saved:", outfile, "\n")
}
