#!/usr/bin/env Rscript
 
# Normalized nucleosome occupancy profile plot (mutant minus WT)
# For each mutant, WT mean occupancy is subtracted at every position,
# so the y-axis shows deviation from WT rather than absolute occupancy.
# A flat line at 0 = no change; positive = gain; negative = loss.
 
#### PACKAGE LOADING ####
library(tidyverse)
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
 
# Summarise two replicates into one mean per condition
summarise_condition <- function(rep_A, rep_B, sample_name) {
    df.A <- getAverageProfileSample(heat.mx = getHeatSamp(rep_A, ncol_sample),
                                    heat.name = paste0(sample_name, ".rA"))
    df.B <- getAverageProfileSample(heat.mx = getHeatSamp(rep_B, ncol_sample),
                                    heat.name = paste0(sample_name, ".rB"))
 
    rbind(df.A, df.B) %>%
        group_by(position, geneOrientation, regularity) %>%
        summarise(mean_nucOccupancy = mean(mean_nucOccupancy, na.rm = TRUE),
                  se_max = mean_se(mean_nucOccupancy)$ymax,
                  se_min = mean_se(mean_nucOccupancy)$ymin,
                  .groups = "drop") %>%
        mutate(sample = sample_name,
               position = as.numeric(as.character(position)))
}
 
######## COMPUTE PER-CONDITION SUMMARIES ########
 
df_wt     <- summarise_condition("WT_A",             "WT_B",             "wt")
df_chd1   <- summarise_condition("chd1_A",           "chd1_B",           "chd1")
df_isw1   <- summarise_condition("isw1_A",           "isw1_B",           "isw1")
df_isw2   <- summarise_condition("isw2_A",           "isw2_B",           "isw2")
df_triple <- summarise_condition("isw1_isw2_chd1_A", "isw1_isw2_chd1_B", "isw1.isw2.chd1")
 
######## NORMALIZE: subtract WT at every position ########
# Join WT occupancy onto each mutant by position + geneOrientation + regularity,
# then compute delta = mutant - WT.
 
wt_ref <- df_wt %>%
    select(position, geneOrientation, regularity, wt_occ = mean_nucOccupancy)
 
normalize <- function(df_mutant) {
    df_mutant %>%
        left_join(wt_ref, by = c("position", "geneOrientation", "regularity")) %>%
        mutate(
            delta     = mean_nucOccupancy - wt_occ,
            delta_max = se_max - wt_occ,   # SE band also shifted
            delta_min = se_min - wt_occ
        )
}
 
df_chd1_n   <- normalize(df_chd1)
df_isw1_n   <- normalize(df_isw1)
df_isw2_n   <- normalize(df_isw2)
df_triple_n <- normalize(df_triple)
 
# Combine mutants only (WT is the reference line at 0)
df_norm <- bind_rows(df_chd1_n, df_isw1_n, df_isw2_n, df_triple_n)
 
df_norm$sample <- factor(df_norm$sample,
                         levels = c("chd1", "isw1", "isw2", "isw1.isw2.chd1"))
 
df_norm$regularity <- factor(df_norm$regularity,
                             levels = c("regularity_Q1", "regularity_Q2",
                                        "regularity_Q3", "regularity_Q4"))
 
######## COLORS & THEME ########
 
sample_colors    <- c("chd1" = "#A82649", "isw1" = "#8C5920",
                      "isw2" = "#E59B41", "isw1.isw2.chd1" = "#2E7D5B")
sample_linetypes <- c("chd1" = "dotted",  "isw1" = "dashed",
                      "isw2" = "twodash", "isw1.isw2.chd1" = "solid")
 
######## PLOT 1: all quartiles, faceted by orientation x regularity ########
# Shows the full picture — delta occupancy across all groups
 
gg_norm_full <- ggplot(df_norm,
                       aes(x = position, y = delta,
                           color = sample, linetype = sample)) +
    # shaded SE ribbon per mutant
    geom_ribbon(aes(ymin = delta_min, ymax = delta_max, fill = sample),
                color = NA, alpha = 0.15) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4, linetype = "solid") +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = sample_colors) +
    scale_fill_manual(values = sample_colors) +
    scale_linetype_manual(values = sample_linetypes) +
    facet_grid(geneOrientation ~ regularity, scales = "free_y") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          legend.position = "right") +
    ylab("Δ nucleosome occupancy (mutant − WT)") +
    xlab("distance from T-1") +
    ggtitle("Normalized nucleosome occupancy: mutant − WT")
 
pdf(paste0(outDir, "/normalized_profile_full.pdf"), width = 12, height = 7)
print(gg_norm_full)
dev.off()
 
ggsave(paste0(outDir, "/normalized_profile_full.png"),
       plot = gg_norm_full, width = 12, height = 7, dpi = 300)
 
######## PLOT 2: Q1 vs Q4 only — sharpest contrast ########
# Focuses on the most vs least regular genes to highlight the clearest effect
 
df_norm_q1q4 <- df_norm %>%
    filter(regularity %in% c("regularity_Q1", "regularity_Q4"))
 
gg_norm_q1q4 <- ggplot(df_norm_q1q4,
                        aes(x = position, y = delta,
                            color = sample, linetype = sample)) +
    geom_ribbon(aes(ymin = delta_min, ymax = delta_max, fill = sample),
                color = NA, alpha = 0.15) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4, linetype = "solid") +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = sample_colors) +
    scale_fill_manual(values = sample_colors) +
    scale_linetype_manual(values = sample_linetypes) +
    facet_grid(geneOrientation ~ regularity, scales = "free_y") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          legend.position = "right") +
    ylab("Δ nucleosome occupancy (mutant − WT)") +
    xlab("distance from T-1") +
    ggtitle("Normalized occupancy: Q1 (least regular) vs Q4 (most regular)")
 
pdf(paste0(outDir, "/normalized_profile_Q1vsQ4.pdf"), width = 8, height = 7)
print(gg_norm_q1q4)
dev.off()
 
ggsave(paste0(outDir, "/normalized_profile_Q1vsQ4.png"),
       plot = gg_norm_q1q4, width = 8, height = 7, dpi = 300)
 
######## PLOT 3: per-mutant separate panels, all quartiles ########
# One page per mutant — useful for detailed inspection of each condition
 
for (mut in levels(df_norm$sample)) {
    df_mut <- df_norm %>% filter(sample == mut)
 
    gg_mut <- ggplot(df_mut,
                     aes(x = position, y = delta, color = regularity)) +
        geom_ribbon(aes(ymin = delta_min, ymax = delta_max, fill = regularity),
                    color = NA, alpha = 0.2) +
        geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
        geom_line(linewidth = 1.0) +
        scale_color_manual(values = c("regularity_Q1" = "#ffb242",
                                      "regularity_Q2" = "#de4f33",
                                      "regularity_Q3" = "#9f2d55",
                                      "regularity_Q4" = "#341648")) +
        scale_fill_manual(values  = c("regularity_Q1" = "#ffb242",
                                      "regularity_Q2" = "#de4f33",
                                      "regularity_Q3" = "#9f2d55",
                                      "regularity_Q4" = "#341648")) +
        facet_grid(geneOrientation ~ ., scales = "free_y") +
        theme_bw(base_size = 11) +
        theme(panel.grid.minor = element_blank(),
              strip.background = element_rect(fill = "grey92"),
              legend.position = "right") +
        ylab("Δ nucleosome occupancy (mutant − WT)") +
        xlab("distance from T-1") +
        ggtitle(paste0("Normalized occupancy: ", mut, " − WT"))
 
    pdf(paste0(outDir, "/normalized_profile_", mut, ".pdf"), width = 6, height = 7)
    print(gg_mut)
    dev.off()
 
    ggsave(paste0(outDir, "/normalized_profile_", mut, ".png"),
           plot = gg_mut, width = 6, height = 7, dpi = 300)
}
 
message("Done. Output files written to: ", outDir)
