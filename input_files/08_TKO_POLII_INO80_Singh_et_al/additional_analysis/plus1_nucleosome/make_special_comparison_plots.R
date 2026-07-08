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

tko_color <- "#A82649"
rpb1_color <- "#E59B41"
ino80_color <- "#2E7DD5"
wt_fill_color <- "grey80"

wt_name <- "_Wild_type_Rep_1"
tko_name <- "_TKO_Rep_1"
rpb1_name <- "_TKO_Rpb1_Rapamycin_120_min"
ino80_name <- "_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep_1"

wt_vals    <- as.numeric(heat.val[wt_name, profile.index])
tko_vals   <- as.numeric(heat.val[tko_name, profile.index])
rpb1_vals  <- as.numeric(heat.val[rpb1_name, profile.index])
ino80_vals <- as.numeric(heat.val[ino80_name, profile.index])

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
# Plot 1: Wild type (filled) vs TKO (line)
############################################

y_low_1  <- 0
y_high_1 <- max(c(wt_vals, tko_vals))
y_pad_1  <- (y_high_1 - y_low_1) * 0.08
y_high_1 <- y_high_1 + y_pad_1

draw_plot1 <- function() {
    par(mar = c(8.5, 4.1, 3.1, 2.1), xaxs = "i", yaxs = "i")

    plot(pos, wt_vals,
         type = "n",
         xlab = NA,
         ylab = "Nucleosome Occupancy",
         main = "+1 nucleosome: Wild type vs TKO",
         ylim = c(y_low_1, y_high_1),
         xlim = c(start.plot, end.plot),
         axes = TRUE)

    polygon(c(pos, rev(pos)), c(wt_vals, rep(0, length(pos))),
            col = wt_fill_color, border = NA)

    lines(pos, wt_vals, col = wt_fill_color, lwd = 2)

    lines(pos, tko_vals, col = tko_color, lwd = 2.5)

    abline(v = 0, lty = 3, col = "grey70")

    usr <- par("usr")
    y_legend_1 <- usr[3] - (usr[4] - usr[3]) * 0.20
    legend(start.plot, y_legend_1,
           bty = "n",
           legend = c("Wild type", "TKO"),
           col = c(wt_fill_color, tko_color),
           lwd = c(2, 2.5),
           fill = c(wt_fill_color, NA),
           border = c(wt_fill_color, NA),
           xpd = TRUE)

    title(xlab = "distance from TSS", line = 3)
}

render_plot(file.path(out.dir, "WT_vs_TKO_filled"), draw_plot1, 7, 7)

############################################
# Plot 2: TKO (line) vs TKO_Rpb1 (line)
############################################

y_low_2  <- 0
y_high_2 <- max(c(tko_vals, rpb1_vals))
y_pad_2  <- (y_high_2 - y_low_2) * 0.08
y_high_2 <- y_high_2 + y_pad_2

draw_plot2 <- function() {
    par(mar = c(8.5, 4.1, 3.1, 2.1), xaxs = "i", yaxs = "i")

    plot(pos, tko_vals,
         type = "l",
         xlab = NA,
         ylab = "Nucleosome Occupancy",
         main = "+1 nucleosome: TKO vs TKO - Pol II",
         ylim = c(y_low_2, y_high_2),
         xlim = c(start.plot, end.plot),
         col = tko_color,
         lwd = 2.5,
         axes = TRUE)

    lines(pos, rpb1_vals, col = rpb1_color, lwd = 2.5)

    abline(v = 0, lty = 3, col = "grey70")

    usr <- par("usr")
    y_legend_2 <- usr[3] - (usr[4] - usr[3]) * 0.20
    legend(start.plot, y_legend_2,
           bty = "n",
           legend = c("TKO", "TKO - Pol II"),
           col = c(tko_color, rpb1_color),
           lwd = 2.5,
           xpd = TRUE)

    title(xlab = "distance from TSS", line = 3)
}

render_plot(file.path(out.dir, "TKO_vs_TKO_Rpb1"), draw_plot2, 7, 7)

############################################
# Plot 3: WT (grey filled) vs TKO (red line) vs TKO_Rpb1 (yellow line)
############################################

y_low_3  <- 0
y_high_3 <- max(c(wt_vals, tko_vals, rpb1_vals))
y_pad_3  <- (y_high_3 - y_low_3) * 0.08
y_high_3 <- y_high_3 + y_pad_3

draw_plot3 <- function() {
    par(mar = c(8.5, 4.1, 3.1, 2.1), xaxs = "i", yaxs = "i")

    plot(pos, wt_vals,
         type = "n",
         xlab = NA,
         ylab = "Nucleosome Occupancy",
         main = "+1 nucleosome: Wild type vs TKO vs TKO - Pol II",
         ylim = c(y_low_3, y_high_3),
         xlim = c(start.plot, end.plot),
         axes = TRUE)

    polygon(c(pos, rev(pos)), c(wt_vals, rep(0, length(pos))),
            col = wt_fill_color, border = NA)

    lines(pos, wt_vals, col = wt_fill_color, lwd = 2)

    lines(pos, tko_vals, col = tko_color, lwd = 2.5)

    lines(pos, rpb1_vals, col = rpb1_color, lwd = 2.5)

    abline(v = 0, lty = 3, col = "grey70")

    usr <- par("usr")
    y_legend_3 <- usr[3] - (usr[4] - usr[3]) * 0.20
    legend(start.plot, y_legend_3,
           bty = "n",
           legend = c("Wild type", "TKO", "TKO - Pol II"),
           col = c(wt_fill_color, tko_color, rpb1_color),
           lwd = c(2, 2.5, 2.5),
           fill = c(wt_fill_color, NA, NA),
           border = c(wt_fill_color, NA, NA),
           xpd = TRUE)

    title(xlab = "distance from TSS", line = 3)
}

render_plot(file.path(out.dir, "WT_vs_TKO_vs_TKO_Rpb1"), draw_plot3, 7, 7)

############################################
# Plot 4: TKO (grey filled) vs TKO_Rpb1 (red line) vs TKO_Rpb1_Ino80 (blue line)
############################################

y_low_4  <- 0
y_high_4 <- max(c(tko_vals, rpb1_vals, ino80_vals))
y_pad_4  <- (y_high_4 - y_low_4) * 0.08
y_high_4 <- y_high_4 + y_pad_4

draw_plot4 <- function() {
    par(mar = c(8.5, 4.1, 3.1, 2.1), xaxs = "i", yaxs = "i")

    plot(pos, tko_vals,
         type = "n",
         xlab = NA,
         ylab = "Nucleosome Occupancy",
         main = "+1 nucleosome: TKO vs TKO - Pol II vs TKO - Pol II - INO80",
         ylim = c(y_low_4, y_high_4),
         xlim = c(start.plot, end.plot),
         axes = TRUE)

    polygon(c(pos, rev(pos)), c(tko_vals, rep(0, length(pos))),
            col = wt_fill_color, border = NA)

    lines(pos, tko_vals, col = wt_fill_color, lwd = 2)

    lines(pos, rpb1_vals, col = tko_color, lwd = 2.5)

    lines(pos, ino80_vals, col = ino80_color, lwd = 2.5)

    abline(v = 0, lty = 3, col = "grey70")

    usr <- par("usr")
    y_legend_4 <- usr[3] - (usr[4] - usr[3]) * 0.20
    legend(start.plot, y_legend_4,
           bty = "n",
           legend = c("TKO", "TKO - Pol II", "TKO - Pol II - INO80"),
           col = c(wt_fill_color, tko_color, ino80_color),
           lwd = c(2, 2.5, 2.5),
           fill = c(wt_fill_color, NA, NA),
           border = c(wt_fill_color, NA, NA),
           xpd = TRUE)

    title(xlab = "distance from TSS", line = 3)
}

render_plot(file.path(out.dir, "TKO_vs_TKO_Rpb1_vs_TKO_Rpb1_Ino80"), draw_plot4, 7, 7)