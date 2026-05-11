#!/usr/bin/env Rscript
 
#### PACKAGE LOADING ####
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(reshape2)
 
#### OUTPUT DIRECTORY ####
outDir <- "."
 
######## DATA LOADING ###########
heatIN <- "/media/linuxmac/Storage2/scripts/02_example_yeast_remodeler/add_INSPECTOR_TSS_scripts_ruven/add_INSPECTOR_TSS_scripts_ruven/heatmap/figures/Fig2/A_TES_nuc_profile/values_Heatmap.txt"
heatmap <- read.delim(heatIN, skip = 1, header = F)
 
# extract bed input
geneList <- heatmap[, 1:6]
colnames(geneList) <- c("chr", "start", "end", "name", "nucPositioning", "strand")
 
# remaining data are heatmap values
heatmapVal <- heatmap[, -(1:6)]
 
# get head info
heatmapHEAD <- read.delim(heatIN, nrows = 1, header = F)
 
# extract info from head
features <- heatmapHEAD %>% str_split("],")
 
# split gene list
gene.groups <- grep("group_labels", features[[1]], value = T) %>%
    str_split_i("\\:\\[", 2) %>% str_split(",") %>% unlist() %>%
    str_split_i("regSorted_", 2) %>% str_split_i("\\.", 1)
 
group.idx <- grep("group_boundaries", features[[1]], value = T) %>%
    str_split_i("\\:\\[", 2) %>% str_split(",") %>% unlist() %>%
    as.numeric()
 
gene.group.list <- list()
for (i in 1:length(gene.groups)) {
    gene.group.list[[gene.groups[i]]] <- (group.idx[i] + 1):(group.idx[i + 1])
}
 
geneList.list <- lapply(gene.group.list, function(x) geneList[x, ])
 
# get positions
binSize <- grep("bin size", features[[1]], value = T) %>% str_split(",") %>%
    unlist() %>% first() %>% str_split_i("\\:\\[", 2) %>% as.numeric()
upstream <- grep("upstream", features[[1]], value = T) %>% str_split(",") %>%
    unlist() %>% first() %>% str_split_i("\\:\\[", 2) %>% as.numeric()
downstream <- grep("downstream", features[[1]], value = T) %>% str_split(",") %>%
    unlist() %>% first() %>% str_split_i("\\:\\[", 2) %>% as.numeric()
 
ncol_sample <- (upstream + downstream) / binSize
pos <- seq(c(-1) * (upstream - binSize / 2), downstream - binSize / 2, binSize)
 
samples <- grep("sample_labels", features[[1]], value = T) %>%
    str_split_i("\\:\\[", 2) %>% str_split(",") %>% unlist()
 
## gene quartiles
geneList.qr <- lapply(geneList.list,
    function(x) {
        qr <- quantile(x$nucPositioning)
        x$regularity <- "regularity_Q1"
        x$regularity[x$nucPositioning > qr["25%"]] <- "regularity_Q2"
        x$regularity[x$nucPositioning > qr["50%"]] <- "regularity_Q3"
        x$regularity[x$nucPositioning > qr["75%"]] <- "regularity_Q4"
        return(x)
    })
 
# merge data sets
df.plot <- rbind(geneList.qr[[1]], geneList.qr[[2]], geneList.qr[[3]])
df.plot$name <- factor(df.plot$name, levels = rev(df.plot$name))
 
######## FUNCTIONS ########
 
getHeatSamp <- function(name, ncol_sample = ncol_sample) {
    idx <- which(samples == name)
    idxHeat <- seq(((idx - 1) * ncol_sample + 1), idx * ncol_sample)
    heat <- heatmapVal[, idxHeat]
    colnames(heat) <- pos
    rownames(heat) <- geneList$name
    return(heat)
}
 
getAverageProfileSample <- function(heat.mx, heat.name) {
    df.heat <- melt(as.matrix(heat.mx), value.name = "nucOccupancy",
        varnames = c("genes", "position"))
    df.complete <- df.heat %>%
        mutate(geneOrientation = case_when(
            genes %in% geneList.list$convergent$name ~ "convergent",
            genes %in% geneList.list$lonely$name     ~ "lonely",
            genes %in% geneList.list$tandem$name     ~ "tandem",
            TRUE ~ NA_character_
        )) %>%
        left_join(df.plot %>% select(name, regularity),
            by = join_by(genes == name))
 
    df.summary <- df.complete %>%
        group_by(position, geneOrientation, regularity) %>%
        summarise(mean_nucOccupancy = mean(nucOccupancy, na.rm = TRUE),
            .groups = "drop") %>%
        mutate(sample = heat.name)
    return(df.summary)
}
 
# palette
Tam <- c("#ffb242", "#de4f33", "#9f2d55", "#341648")
 
get.lm <- function(df) {
    x <- df$position
    y <- df$mean_nucOccupancy_rep
    lm <- loess(y ~ x, span = 0.1)
    loess.line <- predict(lm)
    maxima <- x[which(diff(sign(diff(loess.line))) == -2) + 1]
    n.max <- which.min(abs(maxima))
    df.out <- data.frame(maxima = maxima[1:n.max],
                         peaks = (-n.max):(-1),
                         sample = df$sample[1])
    return(df.out)
}
 
get.slope <- function(df) {
    fit <- lm(maxima ~ peaks, data = df)
    r2 <- summary(fit)$r.squared
    nrl <- fit$coefficients[2]
    return(data.frame(r2 = r2, nrl = nrl, sample = df$sample[1]))
}
 
# Helper: merge two replicates (A and B) and summarise
summarise_condition <- function(rep_A, rep_B, sample_name) {
    df.A <- getAverageProfileSample(heat.mx = getHeatSamp(rep_A, ncol_sample), heat.name = paste0(sample_name, ".rA"))
    df.B <- getAverageProfileSample(heat.mx = getHeatSamp(rep_B, ncol_sample), heat.name = paste0(sample_name, ".rB"))
 
    df.merge <- rbind(df.A, df.B)
 
    df.merge %>%
        group_by(position, geneOrientation, regularity) %>%
        summarise(mean_nucOccupancy_rep = mean(mean_nucOccupancy, na.rm = TRUE),
                  se_nucOccupancy_max = mean_se(mean_nucOccupancy)$ymax,
                  se_nucOccupancy_min = mean_se(mean_nucOccupancy)$ymin,
                  .groups = "drop") %>%
        mutate(sample = sample_name)
}
 
# Helper: run full NRL analysis for one condition and save outputs
run_nrl <- function(df_summary, sample_name) {
    df_group <- df_summary %>%
        group_by(geneOrientation, regularity) %>%
        group_modify(~ get.lm(.x)) %>%
        filter(peaks > (-5))
 
    ggProfile <- ggplot(df_group,
                        aes(x = peaks, y = maxima, color = regularity)) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
        geom_point() +
        scale_color_manual(values = Tam) +
        facet_grid(geneOrientation ~ .) + theme_bw() +
        ylab("distance to T-1 nucleosome") +
        xlab("peaks") +
        ggtitle(sample_name)
 
    pdf(paste0(outDir, "/NRL_", sample_name, ".pdf"), width = 4, height = 3)
    print(ggProfile)
    dev.off()
 
    nrl.table <- df_group %>%
        group_by(geneOrientation, regularity) %>%
        group_modify(~ get.slope(.x))
 
    write_delim(nrl.table, file = paste0(outDir, "/NRL_", sample_name, ".txt"))
 
    return(nrl.table)
}
 
######## RUN PER CONDITION ########
 
df_summary_wt             <- summarise_condition("WT_A",             "WT_B",             "wt")
df_summary_chd1           <- summarise_condition("chd1_A",           "chd1_B",           "chd1")
df_summary_isw1           <- summarise_condition("isw1_A",           "isw1_B",           "isw1")
df_summary_isw2           <- summarise_condition("isw2_A",           "isw2_B",           "isw2")
df_summary_triple         <- summarise_condition("isw1_isw2_chd1_A", "isw1_isw2_chd1_B", "isw1.isw2.chd1")
 
nrl.table.wt              <- run_nrl(df_summary_wt,     "wt")
nrl.table.chd1            <- run_nrl(df_summary_chd1,   "chd1")
nrl.table.isw1            <- run_nrl(df_summary_isw1,   "isw1")
nrl.table.isw2            <- run_nrl(df_summary_isw2,   "isw2")
nrl.table.triple          <- run_nrl(df_summary_triple, "isw1.isw2.chd1")
 
######## SUMMARY PLOT ########
 
nrl.plot.pre <- rbind(nrl.table.wt, nrl.table.chd1,
                      nrl.table.isw1, nrl.table.isw2,
                      nrl.table.triple)
 
nrl.plot <- nrl.plot.pre %>%
    filter(nrl >= 150) %>%
    filter(r2 > 0.99)
 
nrl.plot$regularity <- nrl.plot$regularity %>%
    str_remove("regularity_") %>% as.factor()
 
nrl.plot$sample <- factor(nrl.plot$sample,
                          levels = c("wt", "chd1", "isw1", "isw2", "isw1.isw2.chd1"))
 
gg.nrl <- ggplot(nrl.plot, aes(x = regularity, y = nrl, color = sample, group = sample)) +
    geom_line(aes(linetype = sample)) +
    geom_point() +
    scale_color_manual(values = c("#404040", "#A82649", "#8C5920", "#E59B41", "#2E7D5B")) +
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "twodash", "longdash")) +
    theme_bw() +
    facet_grid(. ~ geneOrientation) +
    ylim(150, 180)
 
pdf(paste0(outDir, "/NRL_summary.pdf"), width = 5, height = 3)
print(gg.nrl)
dev.off()
