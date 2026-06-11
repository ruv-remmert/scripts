import pyBigWig
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
from glob import glob

chrom_map = {
    "V": "IX",
    "VI": "V",
    "VII": "VI",
    "VIII": "VII",
    "IX": "VIII",
}

bw_dir_08 = "/media/linuxmac/Storage2/scripts/08_tko_polII/INSPECTOR_out/RUN/07_REGULARITY/Reference_pos/PSD"
bw_dir_09 = "/media/linuxmac/Storage2/scripts/09_tko_and_ino80/INSPECTOR_out/RUN/07_REGULARITY/Reference_pos/PSD"

extra_files_09 = [
    "PSD_avgNRL_180_PSD_RollingWindow_WT_Ino80_Rapamycin_90_min_Rep_1.bw",
    "PSD_avgNRL_180_PSD_RollingWindow_WT_Ino80_Rapamycin_90_min_Rep_2.bw",
    "PSD_avgNRL_180_PSD_RollingWindow_TKO_Ino80_Rapamycin_90_min_Rep_1.bw",
    "PSD_avgNRL_180_PSD_RollingWindow_TKO_Ino80_Rapamycin_90_min_Rep_2.bw",
]

bw_files = sorted(glob(f"{bw_dir_08}/*.bw"))
for f in extra_files_09:
    bw_files.append(f"{bw_dir_09}/{f}")
bw_files = [f for f in bw_files if "Vehicle" not in f]

gtf_path = "/media/linuxmac/Storage2/scripts/08_tko_polII/data/genes_sacCer3.gtf"
gtf = pd.read_csv(gtf_path, sep='\t', comment='#', header=None,
                  names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
gtf = gtf[gtf['feature'] == 'gene'].copy()
gtf['gene_id'] = gtf['attribute'].str.extract(r'gene_id "([^"]+)"')
gtf['gene_name'] = gtf['attribute'].str.extract(r'gene_name "([^"]+)"')
gtf['tag'] = gtf['attribute'].str.extract(r'tag "([^"]+)"')
gtf = gtf[['seqname', 'gene_id', 'gene_name', 'start', 'end', 'tag']]
gtf = gtf[gtf['seqname'] != 'Mito']
gtf['length'] = gtf['end'] - gtf['start']
gtf_sub = gtf[gtf['length'] >= 400]
n_genes = len(gtf_sub)

label_replacements = [
    ("TKO_Rpb1_Ino80", "TKO - Pol II - INO80"),
    ("TKO_Rpb1", "TKO - Pol II"),
    ("TKO_Ino80", "TKO - INO80"),
    ("WT_Ino80", "WT - INO80"),
    ("Wild_type", "WT"),
]

color_map = {
    "TKO - Pol II - INO80": "#6A9FD8",
    "TKO - Pol II": "#DAB648",
    "TKO - INO80": "#9370DB",
    "WT - INO80": "#6DAF6D",
    "WT": "#999999",
    "TKO": "#D87070",
}

def get_color(label):
    for key, color in color_map.items():
        if label.startswith(key):
            return color
    return "steelblue"

file_labels = []
tss_values = []
tes_values = []
bar_colors = []

for n in bw_files:
    bw = pyBigWig.open(n)
    chroms = bw.chroms()

    TSS_vals = np.full((n_genes, 600), np.nan)
    count = 0
    for row in gtf_sub.itertuples(index=False):
        chrom = chrom_map.get(row.seqname, row.seqname)
        start = int(row.start) - 1
        end = start + 600
        if chrom not in chroms or start < 0 or end > chroms[chrom]:
            continue
        vals = bw.values(chrom, start, end)
        TSS_vals[count] = vals
        count += 1
    tss = np.nanmean(TSS_vals)

    TES_vals = np.full((n_genes, 600), np.nan)
    count = 0
    for row in gtf_sub.itertuples(index=False):
        chrom = chrom_map.get(row.seqname, row.seqname)
        end = int(row.end) - 1
        start = end - 600
        if chrom not in chroms or start < 0 or start > chroms[chrom]:
            continue
        vals = bw.values(chrom, start, end)
        TES_vals[count] = vals
        count += 1
    tes = np.nanmean(TES_vals)

    label = n.rsplit('/', 1)[-1]
    label = label.replace('PSD_avgNRL_180_PSD_RollingWindow_', '').replace('.bw', '').lstrip('_')
    label = re.sub(r'_Rapamycin_\d+_min', '', label)
    for old, new in label_replacements:
        label = label.replace(old, new)

    file_labels.append(label)
    tss_values.append(tss)
    tes_values.append(tes)
    bar_colors.append(get_color(label))
    bw.close()

subset_keywords = ["TKO_Rep", "TKO_Rpb1_Ino80", "TKO_Rpb1_Rapamycin", "Wild_type"]
subset_idx = [i for i, f in enumerate(bw_files) if any(kw in f for kw in subset_keywords)]
subset_labels = [file_labels[i] for i in subset_idx]
subset_tss = [tss_values[i] for i in subset_idx]
subset_tes = [tes_values[i] for i in subset_idx]
subset_colors = [bar_colors[i] for i in subset_idx]

plt.figure(figsize=(14, 6))
plt.bar(file_labels, tss_values, color=bar_colors)
plt.xticks(rotation=45, ha='right')
plt.ylabel('regularity')
plt.tight_layout()
plt.savefig("tss_barplot.png")

plt.figure(figsize=(14, 6))
plt.bar(file_labels, tes_values, color=bar_colors)
plt.xticks(rotation=45, ha='right')
plt.ylabel('regularity')
plt.tight_layout()
plt.savefig("tes_barplot.png")

plt.figure(figsize=(14, 6))
plt.bar(subset_labels, subset_tss, color=subset_colors)
plt.xticks(rotation=45, ha='right')
plt.ylabel('regularity')
plt.tight_layout()
plt.savefig("tss_barplot_subset.png")

plt.figure(figsize=(14, 6))
plt.bar(subset_labels, subset_tes, color=subset_colors)
plt.xticks(rotation=45, ha='right')
plt.ylabel('regularity')
plt.tight_layout()
plt.savefig("tes_barplot_subset.png")

legend_items = list(color_map.items())
fig, ax = plt.subplots(figsize=(4, len(legend_items) * 0.6))
ax.axis('off')
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in legend_items]
labels = [name for name, _ in legend_items]
ax.legend(handles, labels, loc='center', frameon=False, fontsize=14)
fig.tight_layout()
fig.savefig("legend.png", bbox_inches='tight')

rep1_idx = [i for i, f in enumerate(bw_files) if "Rep_1" in f or ("Rep" not in f and "Vehicle" not in f)]
rep1_files = [bw_files[i] for i in rep1_idx]
rep1_labels = [re.sub(r'_Rep_\d+', '', file_labels[i]) for i in rep1_idx]
rep1_tss = [tss_values[i] for i in rep1_idx]
rep1_tes = [tes_values[i] for i in rep1_idx]
rep1_colors = [bar_colors[i] for i in rep1_idx]

n_groups = len(rep1_idx)
x = np.arange(n_groups)
width = 0.35

fig, ax = plt.subplots(figsize=(14, 6))
ax.bar(x - width/2, rep1_tss, width, color=rep1_colors, edgecolor='black', linewidth=0.5, label='TSS')
ax.bar(x + width/2, rep1_tes, width, color=rep1_colors, edgecolor='black', linewidth=0.5, hatch='///', label='TES')
ax.set_xticks(x)
ax.set_xticklabels(rep1_labels, rotation=45, ha='right')
ax.set_ylabel('regularity', fontsize=16)
tss_handle = plt.Rectangle((0, 0), 1, 1, facecolor='white', edgecolor='black', linewidth=0.5)
tes_handle = plt.Rectangle((0, 0), 1, 1, facecolor='white', edgecolor='black', linewidth=0.5, hatch='///')
ax.legend(handles=[tss_handle, tes_handle], labels=['TSS', 'TES'], fontsize=16)
for ticklabel, color in zip(ax.get_xticklabels(), rep1_colors):
    ticklabel.set_color(color)
    ticklabel.set_fontsize(18)
fig.tight_layout()
fig.savefig("tss_tes_combined.png")
