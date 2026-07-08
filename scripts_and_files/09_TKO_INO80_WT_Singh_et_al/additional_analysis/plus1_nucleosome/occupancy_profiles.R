#!/usr/bin/env Rscript

library(RColorBrewer)

######## DATA LOADING ###########
args <- commandArgs(TRUE)

heat.val <- read.delim(file = args[1])
out.dir <- args[2]
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

names.pre <- gsub("_monoNucs_profile", "", as.character(heat.val$bin.labels)[-1])
names.pre <- gsub("_1\\.fastq\\.gz", "", names.pre)
row.names(heat.val) <- c("bin", names.pre)
heat.val <- heat.val[order(heat.val[,1]),]

TES.pos <- which(colnames(heat.val) == "tick")

start.plot <- -500
end.plot <- 900
pos <- seq(start.plot, end.plot, 10)

profile.index <- (TES.pos + start.plot/10):(TES.pos + end.plot/10)

ino80_color <- "#803DB4"
wt_ino80_color <- "#2CA02C"
fill_color <- "grey80"
polii_ino80_color <- "#E5AB00"

tko_name <- "TKO_Rep_1"
wt_name <- "Wild_type_Rep_2"
tko_ino80_name <- "TKO_Ino80_Rapamycin_90_min_Rep_1"
wt_ino80_name <- "WT_Ino80_Rapamycin_90_min_Rep_1"

tko_vals       <- as.numeric(heat.val[tko_name, profile.index])
wt_vals        <- as.numeric(heat.val[wt_name, profile.index])
tko_ino80_vals <- as.numeric(heat.val[tko_ino80_name, profile.index])
wt_ino80_vals  <- as.numeric(heat.val[wt_ino80_name, profile.index])

######## Load TKO-PolII-INO80 data from second file ###########
if (length(args) >= 3 && args[3] != "") {
    heat.val2 <- read.delim(file = args[3])
    names2 <- gsub("_monoNucs_profile", "", as.character(heat.val2$bin.labels)[-1])
    names2 <- gsub("_1\\.fastq\\.gz", "", names2)
    row.names(heat.val2) <- c("bin", names2)
    heat.val2 <- heat.val2[order(heat.val2[,1]),]

    TES.pos2 <- which(colnames(heat.val2) == "tick")
    profile.index2 <- (TES.pos2 + start.plot/10):(TES.pos2 + end.plot/10)

    tko_polii_ino80_name <- "_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep_1"
    tko_polii_ino80_vals <- as.numeric(heat.val2[tko_polii_ino80_name, profile.index2])
    has_polii_ino80 <- TRUE
} else {
    has_polii_ino80 <- FALSE
}

render_plot <- function(filename_base, draw_fun, width, height) {
    pdf(paste0(filename_base, ".pdf"), width = width, height = height)
    draw_fun()
    dev.off()
    cat("Saved:", paste0(filename_base, ".pdf"), "\n")

    png(paste0(filename_base, ".png"), width = width * 150, height = height * 150, res = 150)
    draw_fun()
    dev.off()
    cat("Saved:", paste0(filename_base, ".png"), "\n")
}

############################################
# Plot 1: TKO (grey filled) vs TKO - INO80 (purple line)
############################################

y_low_1  <- 0
y_high_1 <- max(c(tko_vals, tko_ino80_vals))
y_pad_1  <- (y_high_1 - y_low_1) * 0.08
y_high_1 <- y_high_1 + y_pad_1

draw_plot1 <- function() {
    par(mar = c(8.5, 4.1, 3.1, 2.1), xaxs = "i", yaxs = "i")

    plot(pos, tko_vals,
         type = "n",
         xlab = NA,
         ylab = "Nucleosome Occupancy",
         main = "+1 nucleosome: TKO vs TKO - INO80",
         ylim = c(y_low_1, y_high_1),
         xlim = c(start.plot, end.plot),
         axes = TRUE)

    polygon(c(pos, rev(pos)), c(tko_vals, rep(0, length(pos))),
            col = fill_color, border = NA)

    lines(pos, tko_vals, col = fill_color, lwd = 2)

    lines(pos, tko_ino80_vals, col = ino80_color, lwd = 2.5)

    abline(v = 0, lty = 3, col = "grey70")

    usr <- par("usr")
    y_legend_1 <- usr[3] - (usr[4] - usr[3]) * 0.20
    legend(start.plot, y_legend_1,
           bty = "n",
           legend = c("TKO", "TKO - INO80"),
           col = c(fill_color, ino80_color),
           lwd = c(2, 2.5),
           fill = c(fill_color, NA),
           border = c(fill_color, NA),
           xpd = TRUE)

    title(xlab = "distance from TSS", line = 3)
}

render_plot(file.path(out.dir, "TKO_vs_TKO_Ino80_filled"), draw_plot1, 7, 7)

############################################
# Plot 2: Wild type (grey filled) vs WT - INO80 (green line)
############################################

y_low_2  <- 0
y_high_2 <- max(c(wt_vals, wt_ino80_vals))
y_pad_2  <- (y_high_2 - y_low_2) * 0.08
y_high_2 <- y_high_2 + y_pad_2

draw_plot2 <- function() {
    par(mar = c(8.5, 4.1, 3.1, 2.1), xaxs = "i", yaxs = "i")

    plot(pos, wt_vals,
         type = "n",
         xlab = NA,
         ylab = "Nucleosome Occupancy",
         main = "+1 nucleosome: WT vs WT - INO80",
         ylim = c(y_low_2, y_high_2),
         xlim = c(start.plot, end.plot),
         axes = TRUE)

    polygon(c(pos, rev(pos)), c(wt_vals, rep(0, length(pos))),
            col = fill_color, border = NA)

    lines(pos, wt_vals, col = fill_color, lwd = 2)

    lines(pos, wt_ino80_vals, col = wt_ino80_color, lwd = 2.5)

    abline(v = 0, lty = 3, col = "grey70")

    usr <- par("usr")
    y_legend_2 <- usr[3] - (usr[4] - usr[3]) * 0.20
    legend(start.plot, y_legend_2,
           bty = "n",
           legend = c("WT", "WT - INO80"),
           col = c(fill_color, wt_ino80_color),
           lwd = c(2, 2.5),
           fill = c(fill_color, NA),
           border = c(fill_color, NA),
           xpd = TRUE)

    title(xlab = "distance from TSS", line = 3)
}

render_plot(file.path(out.dir, "WT_vs_WT_Ino80_filled"), draw_plot2, 7, 7)

############################################
# Plot 3: TKO (grey filled) vs TKO - INO80 (purple) vs TKO - Pol II - INO80 (blue)
############################################

if (has_polii_ino80) {

    y_low_3  <- 0
    y_high_3 <- max(c(tko_vals, tko_ino80_vals, tko_polii_ino80_vals))
    y_pad_3  <- (y_high_3 - y_low_3) * 0.08
    y_high_3 <- y_high_3 + y_pad_3

    draw_plot3 <- function() {
        par(mar = c(8.5, 4.1, 3.1, 2.1), xaxs = "i", yaxs = "i")

        plot(pos, tko_vals,
             type = "n",
             xlab = NA,
             ylab = "Nucleosome Occupancy",
             main = "+1 nucleosome: TKO vs TKO - INO80 vs TKO - Pol II - INO80",
             ylim = c(y_low_3, y_high_3),
             xlim = c(start.plot, end.plot),
             axes = TRUE)

        polygon(c(pos, rev(pos)), c(tko_vals, rep(0, length(pos))),
                col = fill_color, border = NA)

        lines(pos, tko_vals, col = fill_color, lwd = 2)

        lines(pos, tko_ino80_vals, col = ino80_color, lwd = 2.5)

        lines(pos, tko_polii_ino80_vals, col = polii_ino80_color, lwd = 2.5, lty = 2)

        abline(v = 0, lty = 3, col = "grey70")

        usr <- par("usr")
        y_legend_3 <- usr[3] - (usr[4] - usr[3]) * 0.20
        legend(start.plot, y_legend_3,
               bty = "n",
               legend = c("TKO", "TKO - INO80", "TKO - Pol II - INO80"),
               col = c(fill_color, ino80_color, polii_ino80_color),
               lwd = c(2, 2.5, 2.5),
               lty = c(1, 1, 2),
               fill = c(fill_color, NA, NA),
               border = c(fill_color, NA, NA),
               xpd = TRUE)

        title(xlab = "distance from TSS", line = 3)
    }

    render_plot(file.path(out.dir, "TKO_vs_TKO_Ino80_vs_TKO_PolII_Ino80_filled"), draw_plot3, 7, 7)
}