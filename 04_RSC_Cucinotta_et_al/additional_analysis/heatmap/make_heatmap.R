library(tidyverse)
library(stringr)
library(ggplot2)
library(patchwork)
library(scales)
library(reshape2)

#########################################################
# INPUT
#########################################################

heatIN <- "./figures/Fig2/A_TES_nuc_profile/values_Heatmap.txt"

#########################################################
# LOAD HEATMAP
#########################################################

heatmap <- read.delim(
  heatIN,
  skip = 1,
  header = FALSE
)

geneList <- heatmap[,1:6]

colnames(geneList) <- c(
  "chr",
  "start",
  "end",
  "name",
  "nucPositioning",
  "strand"
)

heatmapVal <- heatmap[, -(1:6)]

#########################################################
# HEADER INFO
#########################################################

heatmapHEAD <- read.delim(
  heatIN,
  nrows = 1,
  header = FALSE
)

features <- str_split(
  as.character(heatmapHEAD[1,1]),
  "],"
)[[1]]

#########################################################
# GENE GROUPS
#########################################################

gene.groups <- grep(
  "group_labels",
  features,
  value = TRUE
) %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1,2] %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  str_replace_all("\"","") %>%
  str_split("regSorted_", simplify = TRUE) %>%
  .[,2] %>%
  str_split("\\.", simplify = TRUE) %>%
  .[,1]

group.idx <- grep(
  "group_boundaries",
  features,
  value = TRUE
) %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1,2] %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  str_replace_all("\\]","") %>%
  as.numeric()

gene.group.list <- list()

for(i in seq_along(gene.groups)){
  gene.group.list[[gene.groups[i]]] <-
    (group.idx[i]+1):(group.idx[i+1])
}

geneList.list <- lapply(
  gene.group.list,
  function(ix) geneList[ix,,drop=FALSE]
)

#########################################################
# FORCE GROUP ORDER
#########################################################

preferred_groups <- c(
  "convergent",
  "lonely",
  "tandem"
)

if(all(preferred_groups %in% names(gene.group.list))){
  gene.group.list <- gene.group.list[preferred_groups]
  geneList.list <- geneList.list[preferred_groups]
}

#########################################################
# WINDOW INFO
#########################################################

binSize <- grep(
  "bin size",
  features,
  value = TRUE
) %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  first() %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1,2] %>%
  as.numeric()

upstream <- grep(
  "upstream",
  features,
  value = TRUE
) %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  first() %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1,2] %>%
  as.numeric()

downstream <- grep(
  "downstream",
  features,
  value = TRUE
) %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  first() %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1,2] %>%
  as.numeric()

ncol_sample <- (upstream + downstream) / binSize

pos <- seq(
  (-1)*(upstream-binSize/2),
  downstream-binSize/2,
  binSize
)

#########################################################
# SAMPLE LABELS
#########################################################

samples <- grep(
  "sample_labels",
  features,
  value = TRUE
) %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1,2] %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  str_replace_all("\\]","") %>%
  str_replace_all("\"","")

#########################################################
# GET SAMPLE MATRIX
#########################################################

getHeatSamp <- function(name){

  idx <- which(samples == name)

  if(length(idx) != 1){
    stop(paste("Sample not found:", name))
  }

  idxHeat <- seq(
    ((idx-1)*ncol_sample+1),
    idx*ncol_sample
  )

  heat <- as.matrix(
    heatmapVal[, idxHeat, drop = FALSE]
  )

  colnames(heat) <- paste0("bin_", seq_len(ncol(heat)))
  rownames(heat) <- geneList$name

  heat
}

#########################################################
# BUILD MATCHED crl vs IAA PAIRS
#########################################################

crl_samples <- samples[grepl("crl", samples)]

sample_pairs <- list()

for(crl in crl_samples){

  iaa <- gsub("crl", "IAA", crl)

  if(!(iaa %in% samples)) next

  pair_name <- gsub("MNase_", "", crl)

  sample_pairs[[pair_name]] <- list(
    crl = crl,
    IAA = iaa
  )
}

#########################################################
# BUILD SAMPLE MATRICES
#########################################################

sample_mats <- list()

for(pair_name in names(sample_pairs)){

  sample_mats[[paste0(pair_name,"_crl")]] <-
    getHeatSamp(sample_pairs[[pair_name]]$crl)

  sample_mats[[paste0(pair_name,"_IAA")]] <-
    getHeatSamp(sample_pairs[[pair_name]]$IAA)
}

#########################################################
# REGULARITY QUARTILES
#########################################################

quartile_levels <- c(
  "regularity_Q1",
  "regularity_Q2",
  "regularity_Q3",
  "regularity_Q4"
)

Tam <- c("#ffb242", "#de4f33", "#9f2d55", "#341648")

df_reg_all <- imap_dfr(
  geneList.list,
  function(df_group, grp){

    qr <- quantile(df_group$nucPositioning, na.rm = TRUE)

    df_group %>%
      mutate(
        group = grp,
        regularity = case_when(
          nucPositioning <= qr["25%"] ~ "regularity_Q1",
          nucPositioning <= qr["50%"] ~ "regularity_Q2",
          nucPositioning <= qr["75%"] ~ "regularity_Q3",
          TRUE ~ "regularity_Q4"
        )
      ) %>%
      select(name, nucPositioning, group, regularity)
  }
)

df_reg_all$regularity <- factor(df_reg_all$regularity, levels = quartile_levels)

#########################################################
# HELPERS
#########################################################

drop_all_na_rows <- function(mat){
  keep <- rowSums(is.na(mat)) < ncol(mat)
  mat[keep,,drop=FALSE]
}

find_tes_col_index <- function(pos, binSize){
  which.min(abs(pos - (-binSize/2)))
}

tes_idx <- find_tes_col_index(pos, binSize)

colWT <- colorRampPalette(c("#FFFFFF", "#381a61", "#381a61"))(100)

#########################################################
# HEATMAP BLOCK
#########################################################

plot_heatmap_block <- function(mat_block,
                               sample_title = NULL,
                               show_title = FALSE,
                               show_x = FALSE,
                               show_legend = FALSE){

  mat_block <- drop_all_na_rows(mat_block)

  row_levels <- rev(rownames(mat_block))
  col_levels <- colnames(mat_block)

  left_col <- col_levels[1]
  right_col <- col_levels[length(col_levels)]
  tes_col <- col_levels[tes_idx]

  df_long <- as.data.frame(mat_block, check.names = FALSE) %>%
    mutate(row_id = factor(rownames(mat_block), levels = row_levels)) %>%
    pivot_longer(-row_id, names_to = "col_id", values_to = "value") %>%
    mutate(col_id = factor(col_id, levels = col_levels))

  p <- ggplot(df_long, aes(x = col_id, y = row_id, fill = value)) +
    geom_raster() +
    scale_fill_gradientn(
      colours = colWT,
      limits = c(0,1000),
      oob = scales::squish,
      na.value = "#FFFFFF"
    ) +
    scale_x_discrete(
      breaks = c(left_col, tes_col, right_col),
      labels = c(as.character(upstream), "TES", as.character(downstream))
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )

  if(show_title) p <- p + ggtitle(sample_title)
  if(!show_x) p <- p + theme(axis.text.x = element_blank())
  if(!show_legend) p <- p + theme(legend.position = "none")

  p
}

#########################################################
# BARPLOT BLOCK
#########################################################

plot_bar_block <- function(df_group, gene_order){

  df_group <- df_group %>%
    filter(name %in% gene_order) %>%
    mutate(name = factor(name, levels = rev(gene_order)))

  ggplot(df_group, aes(x = nucPositioning, y = name, fill = regularity)) +
    geom_col(width = 1) +
    scale_fill_manual(values = Tam) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

#########################################################
# SPLIT BLOCKS
#########################################################

slice_blocks <- function(mat_full){
  imap(gene.group.list, ~drop_all_na_rows(mat_full[.x,,drop=FALSE]))
}

blocks_by_sample <- imap(sample_mats, ~slice_blocks(.x))

#########################################################
# BUILD PLOTS (FIXED LAYOUT)
#########################################################

plots_all <- list()

for(pair_name in names(sample_pairs)){

  crl_name <- paste0(pair_name,"_crl")
  iaa_name <- paste0(pair_name,"_IAA")

  for(g in names(gene.group.list)){

    gene_order <- rownames(blocks_by_sample[[crl_name]][[g]])
    df_group <- df_reg_all %>% filter(group == g)

    plots_all[[length(plots_all)+1]] <-
      plot_bar_block(df_group, gene_order)

    plots_all[[length(plots_all)+1]] <-
      plot_heatmap_block(blocks_by_sample[[crl_name]][[g]],
                         sample_title = crl_name,
                         show_title = TRUE)

    plots_all[[length(plots_all)+1]] <-
      plot_heatmap_block(blocks_by_sample[[iaa_name]][[g]],
                         sample_title = iaa_name,
                         show_title = TRUE,
                         show_legend = TRUE,
                         show_x = TRUE)
  }
}

#########################################################
# FINAL PLOT (FIXED)
#########################################################

final_plot <- wrap_plots(
  plots_all,
  ncol = 3
)

ggsave("combined_TES_heatmap.png",
       plot = final_plot,
       width = 14,
       height = 18,
       dpi = 300)

ggsave("combined_TES_heatmap.pdf",
       plot = final_plot,
       width = 14,
       height = 18)
