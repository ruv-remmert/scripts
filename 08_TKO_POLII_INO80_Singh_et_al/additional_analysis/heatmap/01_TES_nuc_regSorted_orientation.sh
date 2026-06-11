#!/usr/bin/env bash

export AnalysisDir=$(pwd)

## outDir
OutDir=$AnalysisDir"/figures/Fig2/A_TES_nuc_profile/"
mkdir -p $OutDir

Anno="/media/linuxmac/Storage2/scripts/02_example_yeast_remodeler/add_INSPECTOR_TSS_scripts_ruven/add_INSPECTOR_TSS_scripts_ruven/heatmap/"


###############################################################  INO80 dataset   ###########################################################

bwINO80="/media/linuxmac/Storage2/scripts/08_tko_polII/INSPECTOR_out/RUN/01_BIGWIG_PROFILES/"
OutData=$OutDir"/data/"

mkdir -p $OutData

######## TES ########################

computeMatrix reference-point -S $bwINO80"/"*".bw" \
 -R $Anno*'.bed' \
 -o $OutData"/computeMatrix2plot.txt.gz" \
 --referencePoint TES \
 -a 400 -b 750 --smartLabels -p 10


plotHeatmap -m  $OutData"/computeMatrix2plot.txt.gz"  \
     -out $OutData'DefaultHeatmap_TES_regularity.png' \
     --sortRegions no \
     --outFileNameMatrix $OutDir'values_Heatmap.txt'
