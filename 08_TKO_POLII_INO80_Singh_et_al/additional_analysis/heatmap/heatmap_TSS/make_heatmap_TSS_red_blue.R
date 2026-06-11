library(tidyverse)
library(stringr)
library(ggplot2)
library(patchwork)
library(scales)
library(reshape2)
library(grid)

# =========================================================
# INPUT
# =========================================================

heatIN <- "./figures/Fig2/A_TSS_nuc_profile/values_Heatmap.txt"

# =========================================================
# READ HEATMAP
# =========================================================

heatmap <- read.delim(
  heatIN,
  skip = 1,
  header = FALSE
)

geneList <- heatmap[, 1:6]

colnames(geneList) <- c(
  "chr", "start", "end", "name", "nucPositioning", "strand"
)

heatmapVal <- heatmap[, -(1:6)]

# =========================================================
# HEADER METADATA
# =========================================================

heatmapHEAD <- read.delim(heatIN, nrows = 1, header = FALSE)

features <- str_split(as.character(heatmapHEAD[1, 1]), "],")[[1]]

# =========================================================
# GROUP LABELS
# =========================================================

gene.groups <- grep("group_labels", features, value = TRUE) %>%
  str_split("\\:\\[", simplify = TRUE) %>% .[1, 2] %>%
  str_split(",", simplify = FALSE) %>% unlist() %>%
  str_replace_all("\"", "") %>%
  str_split("regSorted_", simplify = TRUE) %>% .[, 2] %>%
  str_split("\\.", simplify = TRUE) %>% .[, 1]

group.idx <- grep("group_boundaries", features, value = TRUE) %>%
  str_split("\\:\\[", simplify = TRUE) %>% .[1, 2] %>%
  str_split(",", simplify = FALSE) %>% unlist() %>%
  str_replace_all("\\]", "") %>%
  as.numeric()

gene.group.list <- list()
for (i in seq_along(gene.groups)) {
  gene.group.list[[gene.groups[i]]] <- (group.idx[i] + 1):(group.idx[i + 1])
}

geneList.list <- lapply(gene.group.list, function(ix) geneList[ix, , drop = FALSE])

# =========================================================
# WINDOW / BINNING
# =========================================================

binSize <- grep("bin size", features, value = TRUE) %>%
  str_split(",", simplify = FALSE) %>% unlist() %>% first() %>%
  str_split("\\:\\[", simplify = TRUE) %>% .[1, 2] %>% as.numeric()

upstream <- grep("upstream", features, value = TRUE) %>%
  str_split(",", simplify = FALSE) %>% unlist() %>% first() %>%
  str_split("\\:\\[", simplify = TRUE) %>% .[1, 2] %>% as.numeric()

downstream <- grep("downstream", features, value = TRUE) %>%
  str_split(",", simplify = FALSE) %>% unlist() %>% first() %>%
  str_split("\\:\\[", simplify = TRUE) %>% .[1, 2] %>% as.numeric()

ncol_sample <- (upstream + downstream) / binSize

pos <- seq((-1) * (upstream - binSize / 2), downstream - binSize / 2, binSize)

# =========================================================
# SAMPLE LABELS
# =========================================================

samples <- grep("sample_labels", features, value = TRUE) %>%
  str_split("\\:\\[", simplify = TRUE) %>% .[1, 2] %>%
  str_split(",", simplify = FALSE) %>% unlist() %>%
  str_replace_all("\\]", "") %>%
  str_replace_all("\"", "")

normalize_sample_name <- function(x) {
  x %>%
    str_remove("^_+") %>%
    str_remove("\\.fastq\\.gz$") %>%
    str_replace_all("Rep_", "Rep") %>%
    str_replace_all("_+", "_")
}

resolve_sample_label <- function(display_name, candidates = character()) {
  wanted      <- normalize_sample_name(c(display_name, unlist(candidates)))
  sample_norm <- normalize_sample_name(samples)
  idx         <- match(wanted, sample_norm, nomatch = 0)
  idx         <- idx[idx > 0]

  if (length(idx) == 0) {
    cat("\nAvailable sample labels:\n")
    print(samples)
    stop(paste("Sample not found for:", display_name))
  }

  samples[idx[1]]
}

sample_label_Wild_type_Rep1 <- resolve_sample_label(
  "Wild_type_Rep1", c("_Wild_type_Rep_1", "Wild_type_Rep_1")
)

sample_label_TKO_Rep1 <- resolve_sample_label(
  "TKO_Rep1", c("_TKO_Rep_1", "TKO_Rep_1")
)

sample_label_TKO_Rpb1_Rapamycin_120_min <- resolve_sample_label(
  "TKO_Rpb1_Rapamycin_120_min", c("_TKO_Rpb1_Rapamycin_120_min")
)

sample_label_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep1 <- resolve_sample_label(
  "TKO_Rpb1_Ino80_Rapamycin_120_min_Rep1",
  c("_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep_1")
)

message("Using samples:")
print(data.frame(
  plot_label = c(
    "Wild_type_Rep1",
    "TKO_Rep1",
    "TKO_Rpb1_Rapamycin_120_min",
    "TKO_Rpb1_Ino80_Rapamycin_120_min_Rep1"
  ),
  header_label = c(
    sample_label_Wild_type_Rep1,
    sample_label_TKO_Rep1,
    sample_label_TKO_Rpb1_Rapamycin_120_min,
    sample_label_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep1
  )
))

# =========================================================
# SAMPLE EXTRACTION
# =========================================================

getHeatSamp <- function(name) {
  idx <- which(samples == name)

  if (length(idx) == 0) {
    cat("\nAvailable sample labels:\n")
    print(samples)
    stop(paste("Sample not found:", name))
  }

  idx     <- idx[1]
  idxHeat <- seq(((idx - 1) * ncol_sample + 1), idx * ncol_sample)
  heat    <- as.matrix(heatmapVal[, idxHeat, drop = FALSE])

  colnames(heat) <- paste0("bin_", seq_len(ncol(heat)))
  rownames(heat) <- geneList$name

  return(heat)
}

# =========================================================
# GROUP ORDER
# =========================================================

preferred_groups <- c("convergent", "lonely", "tandem")

if (all(preferred_groups %in% names(gene.group.list))) {
  gene.group.list <- gene.group.list[preferred_groups]
  geneList.list   <- geneList.list[preferred_groups]
}

# =========================================================
# SAMPLE MATRICES
# =========================================================

sample_mats <- list(
  "Wild_type_Rep1"                        = getHeatSamp(sample_label_Wild_type_Rep1),
  "TKO_Rep1"                              = getHeatSamp(sample_label_TKO_Rep1),
  "TKO_Rpb1_Rapamycin_120_min"            = getHeatSamp(sample_label_TKO_Rpb1_Rapamycin_120_min),
  "TKO_Rpb1_Ino80_Rapamycin_120_min_Rep1" = getHeatSamp(sample_label_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep1)
)

display_titles <- c(
  "Wild_type_Rep1"                        = "WT",
  "TKO_Rep1"                              = "TKO",
  "TKO_Rpb1_Rapamycin_120_min"            = "TKO - Pol II",
  "TKO_Rpb1_Ino80_Rapamycin_120_min_Rep1" = "TKO - Pol II - INO80"
)

# =========================================================
# COLORS
# =========================================================

Tam <- c("#ffb242", "#de4f33", "#9f2d55", "#341648")

quartile_levels <- c("regularity_Q1", "regularity_Q2", "regularity_Q3", "regularity_Q4")

# =========================================================
# QUARTILES
# =========================================================

df_reg_all <- imap_dfr(geneList.list, function(df_group, grp) {
  qr <- quantile(df_group$nucPositioning, na.rm = TRUE)

  df_group %>%
    mutate(
      group = grp,
      regularity = case_when(
        nucPositioning <= qr["25%"] ~ "regularity_Q1",
        nucPositioning <= qr["50%"] ~ "regularity_Q2",
        nucPositioning <= qr["75%"] ~ "regularity_Q3",
        TRUE                        ~ "regularity_Q4"
      )
    ) %>%
    select(name, nucPositioning, group, regularity)
})

df_reg_all$regularity <- factor(df_reg_all$regularity, levels = quartile_levels)

# =========================================================
# HELPERS
# =========================================================

drop_all_na_rows <- function(mat) {
  keep <- rowSums(is.na(mat)) < ncol(mat)
  mat[keep, , drop = FALSE]
}

find_tss_col_index <- function(pos, binSize) {
  target <- -binSize / 2
  which.min(abs(pos - target))
}

heatmap_values      <- unlist(lapply(sample_mats, as.numeric), use.names = FALSE)
heatmap_upper_limit <- as.numeric(quantile(heatmap_values, probs = 0.99, na.rm = TRUE))

if (!is.finite(heatmap_upper_limit) || heatmap_upper_limit <= 0) {
  heatmap_upper_limit <- 1000
}

colWT <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

# =========================================================
# HEATMAP FUNCTION
# =========================================================

plot_heatmap_block <- function(mat_block, sample_title = NULL,
                               show_title = FALSE, show_x = FALSE,
                               show_legend = FALSE, upstream, downstream,
                               tss_label = "TSS", tss_idx) {

  mat_block <- drop_all_na_rows(mat_block)
  row_levels <- rev(rownames(mat_block))
  col_levels <- colnames(mat_block)

  left_col  <- col_levels[1]
  right_col <- col_levels[length(col_levels)]
  tss_col   <- col_levels[tss_idx]

  df_long <- as.data.frame(mat_block, check.names = FALSE) %>%
    mutate(row_id = factor(rownames(mat_block), levels = row_levels)) %>%
    pivot_longer(cols = -row_id, names_to = "col_id", values_to = "value") %>%
    mutate(col_id = factor(col_id, levels = col_levels))

  p <- ggplot(df_long, aes(x = col_id, y = row_id, fill = value)) +
    geom_raster() +
    scale_fill_gradientn(
      colours = colWT,
      limits  = c(0, heatmap_upper_limit),
      oob     = scales::squish,
      na.value = "#FFFFFF",
      name    = "occupancy",
      guide   = guide_colourbar(title.position = "top", barheight = unit(35, "mm"))
    ) +
    scale_x_discrete(
      breaks = c(left_col, tss_col, right_col),
      labels = c(paste0("-", upstream), tss_label, as.character(downstream)),
      drop = FALSE, expand = c(0, 0)
    ) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid    = element_blank(),
      axis.text.y   = element_blank(),
      axis.ticks.y  = element_blank(),
      axis.ticks.x  = element_blank(),
      axis.title    = element_blank(),
      plot.title    = element_text(hjust = 0.5),
      plot.margin   = margin(6, 6, 6, 6)
    )

  if (show_title) p <- p + ggtitle(sample_title)

  if (show_x) {
    p <- p + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0.5))
  } else {
    p <- p + theme(axis.text.x = element_blank())
  }

  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "right")
  }

  return(p)
}

# =========================================================
# BARPLOT FUNCTION
# =========================================================

plot_bar_block <- function(df_group, gene_order, show_x = FALSE, show_legend = FALSE) {

  df_group <- df_group %>%
    filter(name %in% gene_order) %>%
    mutate(name = factor(name, levels = rev(gene_order)))

  p <- ggplot(df_group, aes(x = nucPositioning, y = name, fill = regularity)) +
    geom_col(width = 1) +
    scale_fill_manual(values = Tam, drop = FALSE, name = "regularity") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid   = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin  = margin(6, 6, 6, 6)
    ) +
    ylab(NULL)

  if (!show_x) p <- p + theme(axis.text.x = element_blank())
  if (!show_legend) p <- p + theme(legend.position = "none")

  return(p)
}

# =========================================================
# BUILD BLOCKS
# =========================================================

group_names <- names(gene.group.list)

tss_idx <- find_tss_col_index(pos, binSize)

slice_blocks <- function(mat_full) {
  imap(gene.group.list, function(ix, grp) {
    drop_all_na_rows(mat_full[ix, , drop = FALSE])
  })
}

blocks_by_sample <- imap(sample_mats, function(mat_full, sample_name) slice_blocks(mat_full))

row_counts <- map_int(group_names, function(g) nrow(blocks_by_sample[[1]][[g]]))

# =========================================================
# BARPLOTS
# =========================================================

barplots <- map(group_names, function(g) {
  gene_order  <- rownames(blocks_by_sample[[1]][[g]])
  df_group    <- df_reg_all %>% filter(group == g)
  show_x      <- (g == tail(group_names, 1))
  show_legend <- (g == head(group_names, 1))
  plot_bar_block(df_group, gene_order = gene_order, show_x = show_x, show_legend = show_legend)
})
names(barplots) <- group_names

# =========================================================
# HEATMAPS
# =========================================================

heatmaps <- list()

for (g in group_names) {
  for (s in names(sample_mats)) {

    mat_block   <- blocks_by_sample[[s]][[g]]
    show_title  <- (g == head(group_names, 1))
    show_x      <- (g == tail(group_names, 1))
    show_legend <- (g == tail(group_names, 1) && s == tail(names(sample_mats), 1))

    heatmaps[[paste(g, s, sep = "__")]] <- plot_heatmap_block(
      mat_block,
      sample_title = display_titles[s],
      show_title   = show_title,
      show_x       = show_x,
      show_legend  = show_legend,
      upstream     = upstream,
      downstream   = downstream,
      tss_label    = "TSS",
      tss_idx      = tss_idx
    )
  }
}

# =========================================================
# FINAL HEATMAP LAYOUT
# =========================================================

layout_design <- "
ABCDEZ
FGHIJZ
KLMNOZ
"

plots_in_order <- list(
  # Row 1
  barplots[[group_names[1]]],
  heatmaps[[paste(group_names[1], names(sample_mats)[1], sep = "__")]],
  heatmaps[[paste(group_names[1], names(sample_mats)[2], sep = "__")]],
  heatmaps[[paste(group_names[1], names(sample_mats)[3], sep = "__")]],
  heatmaps[[paste(group_names[1], names(sample_mats)[4], sep = "__")]],

  # Row 2
  barplots[[group_names[2]]],
  heatmaps[[paste(group_names[2], names(sample_mats)[1], sep = "__")]],
  heatmaps[[paste(group_names[2], names(sample_mats)[2], sep = "__")]],
  heatmaps[[paste(group_names[2], names(sample_mats)[3], sep = "__")]],
  heatmaps[[paste(group_names[2], names(sample_mats)[4], sep = "__")]],

  # Row 3
  barplots[[group_names[3]]],
  heatmaps[[paste(group_names[3], names(sample_mats)[1], sep = "__")]],
  heatmaps[[paste(group_names[3], names(sample_mats)[2], sep = "__")]],
  heatmaps[[paste(group_names[3], names(sample_mats)[3], sep = "__")]],
  heatmaps[[paste(group_names[3], names(sample_mats)[4], sep = "__")]],

  guide_area()
)

final_plot <- wrap_plots(plots_in_order, design = layout_design) +
  plot_layout(
    guides  = "collect",
    widths  = c(1.5, 1.5, 1.5, 1.5, 1.5, 0.9),
    heights = as.numeric(row_counts)
  ) &
  theme(legend.position = "right")

ggsave(
  "combined_TSS_ggheatmap_WT_TKO_4samples.png",
  plot   = final_plot,
  width  = 16,
  height = 9,
  dpi    = 300
)

# =========================================================
# QUARTILE AVERAGE PROFILE
# =========================================================

geneList.qr <- lapply(geneList.list, function(x) {
  qr <- quantile(x$nucPositioning)
  x$regularity <- "regularity_Q1"
  x$regularity[x$nucPositioning > qr["25%"]] <- "regularity_Q2"
  x$regularity[x$nucPositioning > qr["50%"]] <- "regularity_Q3"
  x$regularity[x$nucPositioning > qr["75%"]] <- "regularity_Q4"
  return(x)
})

df.plot <- rbind(geneList.qr[[1]], geneList.qr[[2]], geneList.qr[[3]])
df.plot$name <- factor(df.plot$name, levels = rev(df.plot$name))

# =========================================================
# SECOND SAMPLE EXTRACTION
# =========================================================

getHeatSamp2 <- function(name) {
  idx <- which(samples == name)
  if (length(idx) == 0) stop(paste("Sample not found:", name))

  idx     <- idx[1]
  idxHeat <- seq(((idx - 1) * ncol_sample + 1), idx * ncol_sample)
  heat    <- as.matrix(heatmapVal[, idxHeat])

  colnames(heat) <- pos
  rownames(heat) <- make.unique(as.character(geneList$name))

  return(heat)
}

# =========================================================
# AVERAGE PROFILE FUNCTION
# =========================================================

getAverageProfileSample <- function(heat.mx, heat.name) {

  df.heat <- reshape2::melt(
    as.matrix(getHeatSamp2(heat.mx)),
    value.name = "nucOccupancy",
    varnames   = c("genes", "position")
  )

  df.complete <- df.heat %>%
    mutate(
      geneOrientation = case_when(
        genes %in% geneList.list$convergent$name ~ "convergent",
        genes %in% geneList.list$lonely$name     ~ "lonely",
        genes %in% geneList.list$tandem$name     ~ "tandem",
        TRUE ~ NA_character_
      )
    ) %>%
    left_join(df.plot %>% select(name, regularity), by = join_by(genes == name))

  df.summary <- df.complete %>%
    group_by(position, geneOrientation, regularity) %>%
    summarise(mean_nucOccupancy = mean(nucOccupancy, na.rm = TRUE), .groups = "drop") %>%
    mutate(sample = heat.name)

  return(df.summary)
}

# =========================================================
# BUILD PROFILE TABLE
# =========================================================

df_summary <- rbind(
  getAverageProfileSample(sample_label_Wild_type_Rep1,                        "WT"),
  getAverageProfileSample(sample_label_TKO_Rep1,                              "TKO"),
  getAverageProfileSample(sample_label_TKO_Rpb1_Rapamycin_120_min,            "TKO - Pol II"),
  getAverageProfileSample(sample_label_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep1, "TKO - Pol II - INO80")
)

df_summary$position <- as.numeric(df_summary$position)

df_summary$sample <- factor(
  df_summary$sample,
  levels = c(
    "WT",
    "TKO",
    "TKO - Pol II",
    "TKO - Pol II - INO80"
  )
)

# =========================================================
# PROFILE PLOT
# =========================================================

ggProfile <- ggplot(
  df_summary,
  aes(x = position, y = mean_nucOccupancy, color = sample, linetype = sample)
) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("#404040", "#A82649", "#E59B41", "#2E7D5B")) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed", "longdash")) +
  facet_grid(geneOrientation ~ regularity) +
  theme_bw() +
  ylab("nucleosome occupancy") +
  xlab("distance from TSS") +
  ylim(0, 300) +
  theme(
    strip.text      = element_text(size = 10, face = "bold"),
    axis.text       = element_text(size = 9),
    axis.title      = element_text(size = 10),
    legend.text     = element_text(size = 9),
    legend.title    = element_text(size = 10),
    panel.spacing.x = unit(4, "mm"),
    panel.spacing.y = unit(4, "mm")
  )

# =========================================================
# SAVE PROFILE PLOTS
# =========================================================

pdf("reg_score_profile_compare_TSS.pdf", width = 20, height = 12)
print(ggProfile)
dev.off()

ggsave(
  "reg_score_profile_compare_TSS.png",
  plot   = ggProfile,
  width  = 20,
  height = 12,
  dpi    = 300
)