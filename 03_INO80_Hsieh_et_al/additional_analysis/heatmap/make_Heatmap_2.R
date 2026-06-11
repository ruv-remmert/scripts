library(tidyverse)
library(stringr)
library(ggplot2)
library(patchwork)
library(scales)
library(reshape2)

heatIN <- "./figures/Fig2/A_TES_nuc_profile/values_Heatmap.txt"

# Heatmap table: first 6 columns are BED-like gene info; remaining columns are signal
heatmap <- read.delim(heatIN, skip = 1, header = FALSE)

geneList <- heatmap[, 1:6]
colnames(geneList) <- c("chr", "start", "end", "name", "nucPositioning", "strand")

heatmapVal <- heatmap[, -(1:6)]

# Header line with encoded metadata
heatmapHEAD <- read.delim(heatIN, nrows = 1, header = FALSE)
features <- str_split(as.character(heatmapHEAD[1, 1]), "],")[[1]]

# Parse group labels
gene.groups <- grep("group_labels", features, value = TRUE) %>%
  str_split("\\:\\[", simplify = TRUE) %>% .[1, 2] %>%
  str_split(",", simplify = FALSE) %>% unlist() %>%
  str_replace_all("\"", "") %>%
  str_split("regSorted_", simplify = TRUE) %>% .[, 2] %>%
  str_split("\\.", simplify = TRUE) %>% .[, 1]

# Parse group boundaries (row indices)
group.idx <- grep("group_boundaries", features, value = TRUE) %>%
  str_split("\\:\\[", simplify = TRUE) %>% .[1, 2] %>%
  str_split(",", simplify = FALSE) %>% unlist() %>%
  str_replace_all("\\]", "") %>%
  as.numeric()

# Build list of row indices per group (keeps input order exactly)
gene.group.list <- list()
for (i in seq_along(gene.groups)) {
  gene.group.list[[gene.groups[i]]] <- (group.idx[i] + 1):(group.idx[i + 1])
}

# Extract per-group gene meta (for quartiles and barplots)
geneList.list <- lapply(gene.group.list, function(ix) geneList[ix, , drop = FALSE])

# Parse binning and window
binSize <- grep("bin size", features, value = TRUE) %>%
  str_split(",", simplify = FALSE) %>% unlist() %>% first() %>%
  str_split("\\:\\[", simplify = TRUE) %>% .[1, 2] %>%
  as.numeric()

upstream <- grep("upstream", features, value = TRUE) %>%
  str_split(",", simplify = FALSE) %>% unlist() %>% first() %>%
  str_split("\\:\\[", simplify = TRUE) %>% .[1, 2] %>%
  as.numeric()

downstream <- grep("downstream", features, value = TRUE) %>%
  str_split(",", simplify = FALSE) %>% unlist() %>% first() %>%
  str_split("\\:\\[", simplify = TRUE) %>% .[1, 2] %>%
  as.numeric()

ncol_sample <- (upstream + downstream) / binSize
pos <- seq((-1) * (upstream - binSize/2), downstream - binSize/2, binSize)

# Sample labels
samples <- grep("sample_labels", features, value = TRUE) %>%
  str_split("\\:\\[", simplify = TRUE) %>% .[1, 2] %>%
  str_split(",", simplify = FALSE) %>% unlist() %>%
  str_replace_all("\\]", "") %>%
  str_replace_all("\"", "")

getHeatSamp <- function(name) {
  idx <- which(samples == name)
  if (length(idx) != 1) stop("Sample label not uniquely found: ", name)

  idxHeat <- seq(((idx - 1) * ncol_sample + 1), idx * ncol_sample)
  heat <- as.matrix(heatmapVal[, idxHeat, drop = FALSE])

  colnames(heat) <- paste0("bin_", seq_len(ncol(heat)))  # safe / syntactic
  rownames(heat) <- geneList$name

  heat
}

# Prefer a deterministic group order if present
preferred_groups <- c("convergent", "lonely", "tandem")
if (all(preferred_groups %in% names(gene.group.list))) {
  gene.group.list <- gene.group.list[preferred_groups]
  geneList.list <- geneList.list[preferred_groups]
}


pool2 <- function(r1, r2) getHeatSamp(r1) + getHeatSamp(r2)

# --- INO80 dataset: WT and ino80 mutant ---
mat_wt    <- pool2("WT_A",    "WT_B")
mat_ino80 <- pool2("ino80_A", "ino80_B")

sample_mats <- list(
  "wt"    = mat_wt,
  "ino80" = mat_ino80
)

# Same palette you used before (Tam)
Tam <- c("#ffb242", "#de4f33", "#9f2d55", "#341648")
quartile_levels <- c("regularity_Q1", "regularity_Q2", "regularity_Q3", "regularity_Q4")

# Compute quartiles within each group separately (as in your original code)
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


# Drop all-NA rows (only relevant if such rows exist in the real data)
drop_all_na_rows <- function(mat) {
  keep <- rowSums(is.na(mat)) < ncol(mat)
  mat[keep, , drop = FALSE]
}

# Find TES bin index (your original code labels pos == -5 as "T-1"; generalized here)
find_tes_col_index <- function(pos, binSize) {
  target <- -binSize/2
  which.min(abs(pos - target))
}

# Common heatmap colors
colWT <- colorRampPalette(c("#FFFFFF", "#381a61", "#381a61"))(100)

# Create one ggheatmap block that preserves order by using levels_rows/levels_cols
plot_heatmap_block <- function(mat_block, sample_title = NULL, show_title = FALSE,
                               show_x = FALSE, show_legend = FALSE,
                               upstream, downstream, tes_label = "T-1",
                               tes_idx) {

  mat_block <- drop_all_na_rows(mat_block)

  # Preserve exact row & column order
  row_levels <- rev(rownames(mat_block))              # first row on top
  col_levels <- colnames(mat_block)                   # left-to-right order

  left_col  <- col_levels[1]
  right_col <- col_levels[length(col_levels)]
  tes_col   <- col_levels[tes_idx]

  df_long <- as.data.frame(mat_block, check.names = FALSE) %>%
    mutate(row_id = factor(rownames(mat_block), levels = row_levels)) %>%
    pivot_longer(
      cols = -row_id,
      names_to = "col_id",
      values_to = "value"
    ) %>%
    mutate(col_id = factor(col_id, levels = col_levels))

  p <- ggplot(df_long, aes(x = col_id, y = row_id, fill = value)) +
    geom_raster() +
    scale_fill_gradientn(
      colours = colWT,
      limits = c(0, 1000),
      breaks = c(0, 500, 1000, 1500, 2000, 2500),
      oob = scales::squish,
      na.value = "#FFFFFF",
      name = "occupancy",
      guide = guide_colourbar(
        title.position = "top",
        barheight = unit(35, "mm")
      )
    ) +
    scale_x_discrete(
      breaks = c(left_col, tes_col, right_col),
      labels = c(as.character(upstream), tes_label, as.character(downstream)),
      drop = FALSE,
      expand = c(0, 0)
    ) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_blank(),
      legend.text = element_text(colour = "black"),
      legend.title = element_text(colour = "black"),
      plot.title = element_text(colour = "black", hjust = 0.5),
      strip.text = element_blank(),

      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(6, 6, 6, 6)
    )

  if (show_title) {
    p <- p + ggtitle(sample_title)
  } else {
    p <- p + ggtitle(NULL)
  }

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

  p
}

# Barplot per group, aligned to heatmap row order
plot_bar_block <- function(df_group, gene_order, show_x = FALSE, show_legend = FALSE) {

  df_group <- df_group %>%
    filter(name %in% gene_order) %>%
    mutate(name = factor(name, levels = rev(gene_order)))

  p <- ggplot(df_group, aes(x = nucPositioning, y = name, fill = regularity)) +
    geom_col(width = 1) +
    scale_fill_manual(values = Tam, drop = FALSE, name = "regularity") +
    theme_minimal(base_size = 12) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      legend.text = element_text(colour = "black"),
      legend.title = element_text(colour = "black"),

      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(6, 6, 6, 6)
    ) +
    xlab(if (show_x) "regularity score" else NULL) +
    ylab(NULL) +
    scale_x_continuous(
      labels = label_number(scale_cut = cut_short_scale()),
      breaks = breaks_extended(n = 4)
    )

  if (!show_x) {
    p <- p + theme(axis.text.x = element_blank())
  } else {
    p <- p + theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
                   axis.ticks.x = element_line(color = "black"))
  }

  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "right")
  }

  p
}


group_names <- names(gene.group.list)

# Determine TES bin index once (same for all)
tes_idx <- find_tes_col_index(pos, binSize)

# Slice matrices into group blocks
slice_blocks <- function(mat_full) {
  imap(gene.group.list, function(ix, grp) {
    drop_all_na_rows(mat_full[ix, , drop = FALSE])
  })
}

blocks_by_sample <- imap(sample_mats, function(mat_full, sample_name) slice_blocks(mat_full))

# Row counts per group for proportional heights
row_counts <- map_int(group_names, function(g) {
  nrow(blocks_by_sample[[1]][[g]])
})

# Build barplots per group
barplots <- map(group_names, function(g) {
  gene_order <- rownames(blocks_by_sample[[1]][[g]])
  df_group <- df_reg_all %>% filter(group == g)

  show_x      <- (g == tail(group_names, 1))
  show_legend <- (g == head(group_names, 1))

  plot_bar_block(df_group, gene_order = gene_order, show_x = show_x, show_legend = show_legend)
})
names(barplots) <- group_names

# Build heatmaps per group and sample
heatmaps <- list()
for (g in group_names) {
  for (s in names(sample_mats)) {
    mat_block <- blocks_by_sample[[s]][[g]]

    show_title  <- (g == head(group_names, 1))
    show_x      <- (g == tail(group_names, 1))
    show_legend <- (g == tail(group_names, 1) && s == tail(names(sample_mats), 1))

    heatmaps[[paste(g, s, sep = "__")]] <- plot_heatmap_block(
      mat_block,
      sample_title = s,
      show_title = show_title,
      show_x = show_x,
      show_legend = show_legend,
      upstream = upstream,
      downstream = downstream,
      tes_label = "T-1",
      tes_idx = tes_idx
    )
  }
}


# --- Layout: 3 rows (groups) x 4 columns (bar + wt + ino80 + guide) ---
# A=convergent bar, B=convergent wt, C=convergent ino80
# D=lonely bar,     E=lonely wt,     F=lonely ino80
# H=tandem bar,     I=tandem wt,     J=tandem ino80
# Z=guide area spanning all rows
layout_design <- "
ABCZ
DEFZ
HIJZ
"

plots_in_order <- list(
  # Row 1: convergent (A, B, C)
  barplots[[group_names[1]]],
  heatmaps[[paste(group_names[1], "wt",    sep = "__")]],
  heatmaps[[paste(group_names[1], "ino80", sep = "__")]],

  # Row 2: lonely (D, E, F)
  barplots[[group_names[2]]],
  heatmaps[[paste(group_names[2], "wt",    sep = "__")]],
  heatmaps[[paste(group_names[2], "ino80", sep = "__")]],

  # Row 3: tandem (H, I, J)
  barplots[[group_names[3]]],
  heatmaps[[paste(group_names[3], "wt",    sep = "__")]],
  heatmaps[[paste(group_names[3], "ino80", sep = "__")]],

  # Guide area (Z)
  guide_area()
)

final_plot <- wrap_plots(plots_in_order, design = layout_design) +
  plot_layout(
    guides = "collect",
    widths = c(1.5, 1.7, 1.7, 0.9),
    heights = as.numeric(row_counts)
  ) &
  theme(legend.position = "right")

ggsave("combined_TES_ggheatmap_3blocks.png",
       plot = final_plot,
       width = 10, height = 9, dpi = 300)


# ---- quartile average profile ----
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


getHeatSamp2 <- function(name) {
  idx <- which(samples == name)
  idxHeat <- seq(((idx - 1) * ncol_sample + 1), idx * ncol_sample)
  heat <- as.matrix(heatmapVal[, idxHeat])
  colnames(heat) <- pos
  rownames(heat) <- make.unique(as.character(geneList$name))
  return(heat)
}


getAverageProfileSample <- function(heat.mx, heat.name) {
  df.heat <- reshape2::melt(as.matrix(getHeatSamp2(heat.mx)), value.name = "nucOccupancy",
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

summarise_condition <- function(repA, repB, sample_name) {
  rbind(
    getAverageProfileSample(heat.mx = repA, heat.name = paste0(sample_name, ".rA")),
    getAverageProfileSample(heat.mx = repB, heat.name = paste0(sample_name, ".rB"))
  ) %>%
    group_by(position, geneOrientation, regularity) %>%
    summarise(mean_nucOccupancy_rep = mean(mean_nucOccupancy, na.rm = TRUE),
              se_nucOccupancy_max   = mean_se(mean_nucOccupancy)$ymax,
              se_nucOccupancy_min   = mean_se(mean_nucOccupancy)$ymin,
              .groups = "drop") %>%
    mutate(sample = sample_name)
}

# --- INO80 dataset ---
df_summary_plot       <- summarise_condition("WT_A",    "WT_B",    "wt")
df_summary_plot_ino80 <- summarise_condition("ino80_A", "ino80_B", "ino80")

df_summary <- rbind(df_summary_plot, df_summary_plot_ino80)

df_summary$position <- as.numeric(df_summary$position)

df_summary$sample <- factor(df_summary$sample,
                            levels = c("wt", "ino80"))

## plot
ggProfile <- ggplot(df_summary,
                    aes(x = position, y = mean_nucOccupancy_rep, color = sample,
                        linetype = sample)) +
  geom_ribbon(aes(ymin = se_nucOccupancy_min, ymax = se_nucOccupancy_max,
                  fill = regularity),
              color = NA, alpha = 0.5) +
  scale_fill_manual(values = c("grey70", "grey70", "grey70", "grey70")) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = c("#404040", "#2E7D5B")) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  facet_grid(geneOrientation ~ regularity) + theme_bw() +
  ylab("nucleosome occupancy") +
  xlab("distance from T-1") +
  ylim(0, 300)

# Increase panel spacing and text sizes for readability
# 4 regularity columns x 5 inches = 20 wide
# 3 gene orientation rows x 4 inches = 12 tall
ggProfile <- ggProfile +
  theme(
    strip.text      = element_text(size = 10, face = "bold"),
    axis.text       = element_text(size = 9),
    axis.title      = element_text(size = 10),
    legend.text     = element_text(size = 9),
    legend.title    = element_text(size = 10),
    panel.spacing.x = unit(4, "mm"),
    panel.spacing.y = unit(4, "mm")
  )

pdf("reg_score_profile_compare.pdf", width = 20, height = 12)
print(ggProfile)
dev.off()

ggsave("reg_score_profile_compare.png",
       plot = ggProfile,
       width = 20, height = 12, dpi = 300)
