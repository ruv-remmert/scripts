#!/usr/bin/env Rscript

library(tidyverse)
library(reshape2)
library(stringr)

outDir <- "."
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

#########################################################
# LOAD TES HEATMAP
#########################################################

heatIN <- "/media/linuxmac/Storage2/scripts/04_RSC_hypertranscription/additional_analysis/heatmap/figures/Fig2/A_TES_nuc_profile/values_Heatmap.txt"

heatmap <- read.delim(heatIN, skip = 1, header = FALSE)

geneList <- heatmap[, 1:6]
colnames(geneList) <- c("chr","start","end","name","nucPositioning","strand")

heatmapVal <- heatmap[, -(1:6)]

heatmapHEAD <- read.delim(heatIN, nrows = 1, header = FALSE)
features <- str_split(as.character(heatmapHEAD[1,1]), "],")[[1]]

#########################################################
# GENE GROUPS (convergent / lonely / tandem)
#########################################################

gene.groups <- grep("group_labels", features, value = TRUE) %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1,2] %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  str_replace_all("\"","") %>%
  str_split("regSorted_", simplify = TRUE) %>%
  .[,2] %>%
  str_split("\\.", simplify = TRUE) %>%
  .[,1]

group.idx <- grep("group_boundaries", features, value = TRUE) %>%
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

#########################################################
# SAMPLE PARSING
#########################################################

samples <- grep("sample_labels", features, value = TRUE) %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1,2] %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  str_replace_all("\\]","") %>%
  str_replace_all("\"","")

binSize <- grep("bin size", features, value = TRUE) %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  first() %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1,2] %>%
  as.numeric()

upstream <- grep("upstream", features, value = TRUE) %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  first() %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1,2] %>%
  as.numeric()

downstream <- grep("downstream", features, value = TRUE) %>%
  str_split(",", simplify = FALSE) %>%
  unlist() %>%
  first() %>%
  str_split("\\:\\[", simplify = TRUE) %>%
  .[1,2] %>%
  as.numeric()

ncol_sample <- (upstream + downstream) / binSize

pos <- seq((-1)*(upstream-binSize/2),
           downstream-binSize/2,
           binSize)

#########################################################
# HEATMAP EXTRACTOR
#########################################################

getHeatSamp <- function(name){

  idx <- which(samples == name)
  if(length(idx) != 1) stop(paste("Sample not found:", name))

  idxHeat <- seq(((idx-1)*ncol_sample+1), idx*ncol_sample)

  mat <- as.matrix(heatmapVal[, idxHeat])
  colnames(mat) <- pos
  rownames(mat) <- geneList$name

  mat
}

#########################################################
# REGULARITY (Q1–Q4)
#########################################################

df_reg_all <- imap_dfr(
  gene.group.list,
  function(ix, grp){

    df <- geneList[ix, ]

    qr <- quantile(df$nucPositioning, na.rm = TRUE)

    df %>%
      mutate(
        group = grp,
        regularity = case_when(
          nucPositioning <= qr["25%"] ~ "regularity_Q1",
          nucPositioning <= qr["50%"] ~ "regularity_Q2",
          nucPositioning <= qr["75%"] ~ "regularity_Q3",
          TRUE ~ "regularity_Q4"
        )
      )
  }
)

#########################################################
# SAMPLE GROUPING (CRL vs IAA CORRECT)
#########################################################

crl_samples <- samples[grepl("MNase_crl", samples)]
iaa_samples <- samples[grepl("MNase_IAA", samples)]

extract_tag <- function(x){
  str_replace(x, "MNase_(crl|IAA)_", "")
}

crl_tags <- sapply(crl_samples, extract_tag)
iaa_tags <- sapply(iaa_samples, extract_tag)

common_tags <- intersect(crl_tags, iaa_tags)

sample_pairs <- lapply(common_tags, function(tag){
  list(
    tag = tag,
    crl = crl_samples[crl_tags == tag],
    iaa = iaa_samples[iaa_tags == tag]
  )
})
names(sample_pairs) <- common_tags

#########################################################
# PROFILE BUILDING (KEEP FULL STRUCTURE)
#########################################################

getProfile <- function(mat, sample_name){

  df <- melt(as.matrix(mat),
             value.name = "nucOccupancy",
             varnames = c("genes","position"))

  df %>%
    mutate(
      geneOrientation = case_when(
        genes %in% gene.group.list$convergent ~ "convergent",
        genes %in% gene.group.list$lonely     ~ "lonely",
        genes %in% gene.group.list$tandem     ~ "tandem",
        TRUE ~ NA_character_
      )
    ) %>%
    left_join(df_reg_all %>% select(name, regularity),
              by = c("genes" = "name")) %>%
    group_by(position, geneOrientation, regularity) %>%
    summarise(
      mean_nucOccupancy = mean(nucOccupancy, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(sample = sample_name)
}

#########################################################
# BUILD ALL PROFILES
#########################################################

summary_list <- list()

for(tag in names(sample_pairs)){

  crl <- sample_pairs[[tag]]$crl
  iaa <- sample_pairs[[tag]]$iaa

  summary_list[[paste0(tag,"_crl")]] <-
    getProfile(getHeatSamp(crl), paste0(tag,"_crl"))

  summary_list[[paste0(tag,"_IAA")]] <-
    getProfile(getHeatSamp(iaa), paste0(tag,"_IAA"))
}

df_summary <- bind_rows(summary_list)

#########################################################
# NRL CORE FUNCTIONS (RESTORED ORIGINAL LOGIC)
#########################################################

get.lm <- function(df){

  df <- df %>% arrange(position)

  fit <- loess(mean_nucOccupancy ~ position, data = df, span = 0.1)
  smooth <- predict(fit)

  maxima <- df$position[which(diff(sign(diff(smooth))) == -2) + 1]

  if(length(maxima) < 3){
    return(NULL)
  }

  data.frame(
    maxima = maxima,
    peaks = seq_along(maxima),
    sample = df$sample[1]
  )
}

get.slope <- function(df){

  fit <- lm(maxima ~ peaks, data = df)

  data.frame(
    r2 = summary(fit)$r.squared,
    nrl = coef(fit)[2],
    sample = df$sample[1]
  )
}

#########################################################
# RUN NRL (CORRECT FACETING)
#########################################################

nrl_table <- df_summary %>%
  group_by(geneOrientation, regularity, sample) %>%
  group_modify(~ get.lm(.x)) %>%
  group_by(geneOrientation, regularity, sample) %>%
  group_modify(~ get.slope(.x))

#########################################################
# SAVE TABLE
#########################################################

write_delim(nrl_table, file.path(outDir, "NRL_results.txt"))

#########################################################
# PLOT (EXPECTED LINEAR OUTPUT)
#########################################################

gg_nrl <- ggplot(nrl_table,
                 aes(x = peaks, y = nrl,
                     color = regularity,
                     group = regularity)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(geneOrientation ~ sample) +
  theme_bw() +
  ylab("NRL slope") +
  xlab("nucleosome peaks")

ggsave(file.path(outDir, "NRL_plot.pdf"),
       gg_nrl, width = 7, height = 5)
