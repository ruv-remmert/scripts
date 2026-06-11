#!/usr/bin/env bash

export AnalysisDir=$(pwd)

OutDir=$AnalysisDir"/figures/Fig2/A_TES_nuc_profile/"
mkdir -p $OutDir

Anno=$AnalysisDir

bwChen="/media/linuxmac/Storage2/scripts/03_INO80/INSPECTOR_out/RUN/01_BIGWIG_PROFILES/"
OutData=$OutDir"/data/"

mkdir -p $OutData

computeMatrix reference-point -S $bwChen"/"*".bw" \
 -R $Anno/TES_regSorted_convergent.bed \
    $Anno/TES_regSorted_lonely.bed \
    $Anno/TES_regSorted_tandem.bed \
 -o $OutData"/computeMatrix2plot.txt.gz" \
 --referencePoint TES \
 -a 400 -b 750 --smartLabels -p 10

plotHeatmap -m  $OutData"/computeMatrix2plot.txt.gz"\
     -out $OutData'DefaultHeatmap_TES_regularity.png' \
     --sortRegions no \
     --outFileNameMatrix $OutDir'values_Heatmap.txt'
