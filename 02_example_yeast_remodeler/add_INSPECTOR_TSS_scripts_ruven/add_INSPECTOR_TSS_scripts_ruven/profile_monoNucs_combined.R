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

#get min/max values
for(i in names){
    curr_min<-min(heat.val[i,profile.index])
    curr_max<-max(heat.val[i,profile.index])
    if(i==names[1]){
        low<-curr_min
        high<-curr_max
    } else {
        low<-min(curr_min,low)
        high<-max(curr_max,high)
    }
}

# add padding to y axis so lines aren't cramped
y_pad <- (high - low) * 0.1
y_low <- low - y_pad
y_high <- high + y_pad

pdf("profile_monoNucs.pdf", width = 14, height = 7)  # wider, taller figure
par(mar = c(6, 6, 5, 16),                             # more margin space
    cex.axis = 1.3,                                    # bigger axis tick labels
    cex.lab  = 1.5,                                    # bigger axis titles
    cex.main = 1.7)                                    # bigger plot title

for(i in 1:length(names)){
    if(i==1){
        plot(pos,
             heat.val[names[i], profile.index],
             type = "l",
             xlab = "Distance from TSS (bp)",
             ylab = "MNase fragment density",
             lwd = 2.5,                                # thicker lines
             main = "TSS profile",
             col = dark2[i],
             ylim = c(y_low, y_high))
        abline(v = 0, lty = 2, col = "grey50")        # vertical line at TSS
    } else {
        lines(pos,
              heat.val[names[i], profile.index],
              col = dark2[i],
              lwd = 2.5)
    }
}

legend(1600, y_high,
       bty = "n",
       legend = names,
       col = dark2[1:length(names)],
       lwd = 2.5,
       cex = 1.2,                                      # bigger legend text
       xpd = TRUE)

dev.off()
