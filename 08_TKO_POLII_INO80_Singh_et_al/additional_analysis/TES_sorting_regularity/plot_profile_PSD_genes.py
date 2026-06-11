#!/usr/bin/env python3
"""
Plot MNase-seq profile, nucleosome regularity, and gene annotations
for a user-specified genomic region.

Usage:
    python plot_profile_PSD_genes.py --chr XI --start 49399 --end 54078
    python plot_profile_PSD_genes.py --chr XV --start 100000 --end 150000
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
from scipy.ndimage import gaussian_filter1d

BASE_DIR = "/media/linuxmac/Storage2/scripts/08_tko_polII"
PROFILE_BW = os.path.join(BASE_DIR, "INSPECTOR_out/RUN/01_BIGWIG_PROFILES/_Wild_type_Rep_1.bw")
PSD_BW = os.path.join(BASE_DIR, "INSPECTOR_out/RUN/07_REGULARITY/Reference_pos/PSD/PSD_avgNRL_180_PSD_RollingWindow__Wild_type_Rep_1.bw")
GTF_FILE = os.path.join(BASE_DIR, "data/genes_sacCer3.gtf")
PROFILE_COLOR = "#2166AC"
PSD_COLOR = "#D6604D"


def bw_values(bw_path, chrom, start, end):
    bw = pyBigWig.open(bw_path)
    vals = bw.values(chrom, start, end)
    bw.close()
    vals = np.array(vals, dtype=float)
    vals[np.isnan(vals)] = 0
    return vals


def parse_gtf_genes(gtf_path, chrom, start, end):
    genes = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue
            seqname, _, feat_type, g_start, g_end, _, strand, _, attrs = parts
            g_start, g_end = int(g_start), int(g_end)
            if feat_type != "gene" or seqname != chrom:
                continue
            if g_end < start or g_start > end:
                continue
            gene_name = ""
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1]
                    break
            genes.append(dict(start=g_start, end=g_end, strand=strand, name=gene_name))
    return pd.DataFrame(genes) if genes else pd.DataFrame(columns=["start", "end", "strand", "name"])


def main():
    parser = argparse.ArgumentParser(description="Plot MNase profile, nuc regularity and genes for a genomic region.")
    parser.add_argument("--chr", default="XI", help="Chromosome (bigwig convention, e.g. XI)")
    parser.add_argument("--start", type=int, default=49399, help="Region start (0-based)")
    parser.add_argument("--end", type=int, default=54078, help="Region end (0-based, exclusive)")
    parser.add_argument("--outdir", default=None, help="Output directory")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--smooth", type=int, default=20, help="Gaussian smoothing sigma in bp (0 = no smoothing)")
    parser.add_argument("--primary_gene", default="LOS1", help="Gene name to highlight (thick line)")
    args = parser.parse_args()

    chrom = args.chr
    start, end = args.start, args.end
    outdir = args.outdir or os.path.dirname(os.path.abspath(__file__))
    os.makedirs(outdir, exist_ok=True)

    tag = f"{chrom}_{start}_{end}"
    out_path = os.path.join(outdir, f"profile_regularity_{tag}.png")

    positions = np.arange(start, end)

    fig, axes = plt.subplots(3, 1, figsize=(8, 6), gridspec_kw={"height_ratios": [3, 1, 1]}, sharex=True, constrained_layout=True)

    # --- MNase profile -------------------------------------------------------
    ax_prof = axes[0]
    vals = bw_values(PROFILE_BW, chrom, start, end)
    if args.smooth > 0:
        vals = gaussian_filter1d(vals, sigma=args.smooth)
    ax_prof.fill_between(positions, vals, color=PROFILE_COLOR, alpha=0.4)
    ax_prof.plot(positions, vals, color=PROFILE_COLOR, linewidth=0.6)
    ax_prof.set_ylabel("MNase reads", fontsize=9)
    ax_prof.margins(x=0)

    # --- Nucleosome regularity heatmap ----------------------------------------
    ax_nuc = axes[1]
    psd_vals = bw_values(PSD_BW, chrom, start, end)
    im = ax_nuc.imshow(psd_vals[np.newaxis, :], aspect="auto",
                       cmap="RdBu_r", interpolation="nearest",
                       extent=[start, end, 0, 1])
    ax_nuc.set_yticks([])
    ax_nuc.set_ylabel("Regularity", fontsize=9)
    cbar = fig.colorbar(im, ax=ax_nuc, orientation="vertical", fraction=0.02, pad=0.01)
    cbar.ax.tick_params(labelsize=6)

    # --- Gene annotations ----------------------------------------------------
    ax_gene = axes[2]
    gene_df = parse_gtf_genes(GTF_FILE, chrom, start, end)
    if not gene_df.empty:
        for _, row in gene_df.iterrows():
            g_start = max(row["start"], start)
            g_end = min(row["end"], end)
            y = 0.5
            is_primary = row["name"] == args.primary_gene
            lw = 2.5 if is_primary else 1.5
            color = "black" if is_primary else "#888888"
            ax_gene.plot([g_start, g_end], [y, y], color=color, linewidth=lw,
                         solid_capstyle="butt", clip_on=True)
            arrow_dir = 1 if row["strand"] == "+" else -1
            arrow_x = g_end if arrow_dir == 1 else g_start
            dx = 0.015 * (end - start) * arrow_dir
            ax_gene.annotate(
                "", xy=(arrow_x + dx, y), xytext=(arrow_x, y),
                arrowprops=dict(arrowstyle="->,head_width=0.4,head_length=0.3",
                                color=color, lw=lw),
                clip_on=True,
            )
            mid = (g_start + g_end) / 2
            ax_gene.text(mid, y + 0.2, row["name"], ha="center", va="bottom",
                         fontsize=7 if is_primary else 6,
                         fontweight="bold" if is_primary else "normal",
                         color=color)
    ax_gene.set_ylim(0, 1)
    ax_gene.set_xlim(start, end)
    ax_gene.set_yticks([])
    ax_gene.set_ylabel("Genes", fontsize=9)
    ax_gene.set_xlabel(f"Position on chr{chrom} (bp)", fontsize=10)
    ax_gene.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
    fig.suptitle(f"chr{chrom}:{start:,}-{end:,}  (Wild_type_Rep_1)", fontsize=12)
    fig.savefig(out_path, dpi=args.dpi)
    print(f"Plot saved to: {out_path}")


if __name__ == "__main__":
    main()