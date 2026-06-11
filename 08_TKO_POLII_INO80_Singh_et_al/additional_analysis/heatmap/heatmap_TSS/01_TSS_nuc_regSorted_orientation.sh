#!/usr/bin/env bash

export AnalysisDir=$(pwd)

OutDir=$AnalysisDir"/figures/Fig2/A_TSS_nuc_profile/"
mkdir -p $OutDir

Anno="/media/linuxmac/Storage2/scripts/08_tko_polII/additional_analysis/heatmap/heatmap_TSS"
bwChen="/media/linuxmac/Storage2/scripts/08_tko_polII/INSPECTOR_out/RUN/01_BIGWIG_PROFILES/"

OutData=$OutDir"/data/"
mkdir -p $OutData

computeMatrix reference-point -S $bwChen"/"*".bw" \
 -R $Anno/TSS_regSorted_divergent.bed \
    $Anno/TSS_regSorted_lonely.bed \
    $Anno/TSS_regSorted_tandem.bed \
 -o $OutData"/computeMatrix2plot.txt.gz" \
 --referencePoint TSS \
 -a 750 -b 400 --smartLabels -p 10

plotHeatmap -m  $OutData"/computeMatrix2plot.txt.gz"  \
     -out $OutData'DefaultHeatmap_TSS_regularity.png' \
     --sortRegions no \
     --outFileNameMatrix $OutDir'values_Heatmap.txt'
