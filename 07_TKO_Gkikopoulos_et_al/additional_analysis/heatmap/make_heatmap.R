library(tidyverse)
library(stringr)
library(ggplot2)
library(patchwork)
library(scales)
library(reshape2)

# =========================================================
# INPUT
# =========================================================

heatIN <- "./figures/Fig2/A_TES_nuc_profile/values_Heatmap.txt"

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
  "chr",
  "start",
  "end",
  "name",
  "nucPositioning",
  "strand"
)

heatmapVal <- heatmap[, -(1:6)]

# =========================================================
# HEADER METADATA
# =========================================================

heatmapHEAD <- read.delim(
  heatIN,
  nrows = 1,
  header = FALSE
)

features <- str_split(
  as.character(heatmapHEAD[1, 1]),
  "],"
)[[1]]

# =========================================================
# GROUP LABELS
# =========================================================

gene.groups <- grep(
  "group_labels",
  features,
  value = TRUE
) %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1, 2] %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  str_replace_all("\"", "") %>%
  str_split("regSorted_", simplify = TRUE) %>%
  .[, 2] %>%
  str_split("\\.", simplify = TRUE) %>%
  .[, 1]

group.idx <- grep(
  "group_boundaries",
  features,
  value = TRUE
) %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1, 2] %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  str_replace_all("\\]", "") %>%
  as.numeric()

gene.group.list <- list()

for (i in seq_along(gene.groups)) {

  gene.group.list[[gene.groups[i]]] <-
    (group.idx[i] + 1):(group.idx[i + 1])
}

geneList.list <- lapply(
  gene.group.list,
  function(ix) geneList[ix, , drop = FALSE]
)

# =========================================================
# WINDOW / BINNING
# =========================================================

binSize <- grep(
  "bin size",
  features,
  value = TRUE
) %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  first() %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1, 2] %>%
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
  .[1, 2] %>%
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
  .[1, 2] %>%
  as.numeric()

ncol_sample <- (upstream + downstream) / binSize

pos <- seq(
  (-1) * (upstream - binSize / 2),
  downstream - binSize / 2,
  binSize
)

# =========================================================
# SAMPLE LABELS
# =========================================================

samples <- grep(
  "sample_labels",
  features,
  value = TRUE
) %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1, 2] %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  str_replace_all("\\]", "") %>%
  str_replace_all("\"", "")

# =========================================================
# SAMPLE EXTRACTION
# =========================================================

getHeatSamp <- function(name) {

  idx <- which(samples == name)

  if (length(idx) == 0) {

    cat("\nAvailable sample labels:\n")
    print(samples)

    stop(
      paste(
        "Sample not found:",
        name
      )
    )
  }

  idx <- idx[1]

  idxHeat <- seq(
    ((idx - 1) * ncol_sample + 1),
    idx * ncol_sample
  )

  heat <- as.matrix(
    heatmapVal[, idxHeat, drop = FALSE]
  )

  colnames(heat) <- paste0(
    "bin_",
    seq_len(ncol(heat))
  )

  rownames(heat) <- geneList$name

  return(heat)
}

# =========================================================
# GROUP ORDER
# =========================================================

preferred_groups <- c(
  "convergent",
  "lonely",
  "tandem"
)

if (all(preferred_groups %in% names(gene.group.list))) {

  gene.group.list <- gene.group.list[
    preferred_groups
  ]

  geneList.list <- geneList.list[
    preferred_groups
  ]
}

# =========================================================
# SAMPLE MATRICES
# =========================================================

sample_mats <- list(

  "WT_MNase" = getHeatSamp(
    "WT_Mnase_1.fastq.gz"
  ),

  "chd1del_isw1del" = getHeatSamp(
    "chd1del_isw1del_MNase_1.fastq.gz"
  ),

  "triple_del" = getHeatSamp(
    "chd1del_isw1del_isw2del_MNase_1.fastq.gz"
  ),

  "chd1del" = getHeatSamp(
    "chd1del_MNase_1.fastq.gz"
  )
)

# =========================================================
# COLORS
# =========================================================

Tam <- c(
  "#ffb242",
  "#de4f33",
  "#9f2d55",
  "#341648"
)

quartile_levels <- c(
  "regularity_Q1",
  "regularity_Q2",
  "regularity_Q3",
  "regularity_Q4"
)

# =========================================================
# QUARTILES
# =========================================================

df_reg_all <- imap_dfr(
  geneList.list,
  function(df_group, grp) {

    qr <- quantile(
      df_group$nucPositioning,
      na.rm = TRUE
    )

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
      select(
        name,
        nucPositioning,
        group,
        regularity
      )
  }
)

df_reg_all$regularity <- factor(
  df_reg_all$regularity,
  levels = quartile_levels
)

# =========================================================
# HELPERS
# =========================================================

drop_all_na_rows <- function(mat) {

  keep <- rowSums(is.na(mat)) < ncol(mat)

  mat[keep, , drop = FALSE]
}

find_tes_col_index <- function(pos, binSize) {

  target <- -binSize / 2

  which.min(abs(pos - target))
}

colWT <- colorRampPalette(
  c(
    "#FFFFFF",
    "#381a61",
    "#381a61"
  )
)(100)

# =========================================================
# HEATMAP FUNCTION
# =========================================================

plot_heatmap_block <- function(
  mat_block,
  sample_title = NULL,
  show_title = FALSE,
  show_x = FALSE,
  show_legend = FALSE,
  upstream,
  downstream,
  tes_label = "T-1",
  tes_idx
) {

  mat_block <- drop_all_na_rows(mat_block)

  row_levels <- rev(rownames(mat_block))

  col_levels <- colnames(mat_block)

  left_col  <- col_levels[1]
  right_col <- col_levels[length(col_levels)]
  tes_col   <- col_levels[tes_idx]

  df_long <- as.data.frame(
    mat_block,
    check.names = FALSE
  ) %>%

    mutate(
      row_id = factor(
        rownames(mat_block),
        levels = row_levels
      )
    ) %>%

    pivot_longer(
      cols = -row_id,
      names_to = "col_id",
      values_to = "value"
    ) %>%

    mutate(
      col_id = factor(
        col_id,
        levels = col_levels
      )
    )

  p <- ggplot(
    df_long,
    aes(
      x = col_id,
      y = row_id,
      fill = value
    )
  ) +

    geom_raster() +

    scale_fill_gradientn(
      colours = colWT,
      limits = c(0, 1000),
      oob = scales::squish,
      na.value = "#FFFFFF",
      name = "occupancy",
      guide = guide_colourbar(
        title.position = "top",
        barheight = unit(35, "mm")
      )
    ) +

    scale_x_discrete(
      breaks = c(
        left_col,
        tes_col,
        right_col
      ),

      labels = c(
        as.character(upstream),
        tes_label,
        as.character(downstream)
      ),

      drop = FALSE,
      expand = c(0, 0)
    ) +

    scale_y_discrete(
      drop = FALSE,
      expand = c(0, 0)
    ) +

    theme_minimal(base_size = 12) +

    theme(
      panel.grid = element_blank(),

      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),

      axis.title = element_blank(),

      plot.title = element_text(
        hjust = 0.5
      ),

      plot.margin = margin(
        6,
        6,
        6,
        6
      )
    )

  if (show_title) {

    p <- p + ggtitle(sample_title)
  }

  if (show_x) {

    p <- p + theme(
      axis.text.x = element_text(
        angle = 270,
        vjust = 0.5,
        hjust = 0.5
      )
    )

  } else {

    p <- p + theme(
      axis.text.x = element_blank()
    )
  }

  if (!show_legend) {

    p <- p + theme(
      legend.position = "none"
    )

  } else {

    p <- p + theme(
      legend.position = "right"
    )
  }

  return(p)
}

# =========================================================
# BARPLOT FUNCTION
# =========================================================

plot_bar_block <- function(
  df_group,
  gene_order,
  show_x = FALSE,
  show_legend = FALSE
) {

  df_group <- df_group %>%
    filter(name %in% gene_order) %>%
    mutate(
      name = factor(
        name,
        levels = rev(gene_order)
      )
    )

  p <- ggplot(
    df_group,
    aes(
      x = nucPositioning,
      y = name,
      fill = regularity
    )
  ) +

    geom_col(width = 1) +

    scale_fill_manual(
      values = Tam,
      drop = FALSE,
      name = "regularity"
    ) +

    theme_minimal(base_size = 12) +

    theme(
      panel.grid = element_blank(),

      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),

      plot.margin = margin(
        6,
        6,
        6,
        6
      )
    ) +

    ylab(NULL)

  if (!show_x) {

    p <- p + theme(
      axis.text.x = element_blank()
    )
  }

  if (!show_legend) {

    p <- p + theme(
      legend.position = "none"
    )
  }

  return(p)
}

# =========================================================
# BUILD BLOCKS
# =========================================================

group_names <- names(gene.group.list)

tes_idx <- find_tes_col_index(
  pos,
  binSize
)

slice_blocks <- function(mat_full) {

  imap(
    gene.group.list,
    function(ix, grp) {

      drop_all_na_rows(
        mat_full[ix, , drop = FALSE]
      )
    }
  )
}

blocks_by_sample <- imap(
  sample_mats,
  function(mat_full, sample_name) {

    slice_blocks(mat_full)
  }
)

row_counts <- map_int(
  group_names,
  function(g) {

    nrow(
      blocks_by_sample[[1]][[g]]
    )
  }
)

# =========================================================
# BARPLOTS
# =========================================================

barplots <- map(
  group_names,
  function(g) {

    gene_order <- rownames(
      blocks_by_sample[[1]][[g]]
    )

    df_group <- df_reg_all %>%
      filter(group == g)

    show_x <- (
      g == tail(group_names, 1)
    )

    show_legend <- (
      g == head(group_names, 1)
    )

    plot_bar_block(
      df_group,
      gene_order = gene_order,
      show_x = show_x,
      show_legend = show_legend
    )
  }
)

names(barplots) <- group_names

# =========================================================
# HEATMAPS
# =========================================================

heatmaps <- list()

for (g in group_names) {

  for (s in names(sample_mats)) {

    mat_block <- blocks_by_sample[[s]][[g]]

    show_title <- (
      g == head(group_names, 1)
    )

    show_x <- (
      g == tail(group_names, 1)
    )

    show_legend <- (
      g == tail(group_names, 1) &&
      s == tail(names(sample_mats), 1)
    )

    heatmaps[[paste(g, s, sep = "__")]] <-
      plot_heatmap_block(
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

# =========================================================
# FINAL LAYOUT
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

final_plot <- wrap_plots(
  plots_in_order,
  design = layout_design
) +

  plot_layout(
    guides = "collect",

    widths = c(
      1.5,
      1.5,
      1.5,
      1.5,
      1.5,
      0.9
    ),

    heights = as.numeric(row_counts)
  ) &

  theme(
    legend.position = "right"
  )

ggsave(
  "combined_TES_ggheatmap_4samples.png",
  plot = final_plot,
  width = 16,
  height = 9,
  dpi = 300
)
