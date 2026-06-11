#!/usr/bin/env bash

export AnalysisDir=$(pwd)
set -e

Anno=$AnalysisDir"/data"

## outDir
OutDir=$AnalysisDir"/figures/Fig1/"
mkdir -p $OutDir"/A_featureDistribution/"
cd $OutDir"/A_featureDistribution/"

###############################################################  Chen et al   ###########################################################
mkdir -p $OutDir"/A_featureDistribution/data"

# positions of dynamicNucs changes nucleosomes
ChenDir=$AnalysisDir"/INSPECTOR_out/RUN/"

##regularity
bedtools merge -i "/media/linuxmac/Storage2/scripts/02_example_yeast_remodeler/add_INSPECTOR_TSS_scripts_ruven/add_INSPECTOR_TSS_scripts_ruven/feature_occurrence/INSPECTOR_out/RUN/07_REGULARITY/Reference_pos/Irregular_Regions_bs800_10xlog_reduced_filt.bed
" \
  | bedtools genomecov -bga -i - -g $AnalysisDir"/Profiler_out/RUN/chrom_Sizes.txt" \
      | sort -k1,1 -k2,2n  > data/Irregular_Regions.bedgraph

bedGraphToBigWig data/Irregular_Regions.bedgraph \
      $AnalysisDir"/Profiler_out/RUN/chrom_Sizes.txt" \
      data/Irregular_Regions.bw

## occupancy
sort -k1,1 -k2,2n $ChenDir"/06_DYNAMIC_NUCS/occupancy/positions_varOCC.bed" \
| bedtools merge -i - \
  | bedtools genomecov -bga -i - -g $AnalysisDir"/Profiler_out/RUN/chrom_Sizes.txt" \
      | sort -k1,1 -k2,2n  > data/Occupancy.bedgraph

bedGraphToBigWig data/Occupancy.bedgraph \
      $AnalysisDir"/Profiler_out/RUN/chrom_Sizes.txt" \
      data/Occupancy.bw

## shift
sort -k1,1 -k2,2n $ChenDir"/06_DYNAMIC_NUCS/shift/positions_varSHIFTS.bed" \
| bedtools merge -i - \
  | bedtools genomecov -bga -i - -g $AnalysisDir"/Profiler_out/RUN/chrom_Sizes.txt" \
      | sort -k1,1 -k2,2n  > data/Shifts.bedgraph

bedGraphToBigWig data/Shifts.bedgraph \
      $AnalysisDir"/Profiler_out/RUN/chrom_Sizes.txt" \
      data/Shifts.bw

## fuzziness
sort -k1,1 -k2,2n $ChenDir"/06_DYNAMIC_NUCS/fuzziness/positions_varFUZZ.bed" \
| bedtools merge -i - \
  | bedtools genomecov -bga -i - -g $AnalysisDir"/Profiler_out/RUN/chrom_Sizes.txt" \
      | sort -k1,1 -k2,2n  > data/Fuzziness.bedgraph

bedGraphToBigWig data/Fuzziness.bedgraph \
      $AnalysisDir"/Profiler_out/RUN/chrom_Sizes.txt" \
      data/Fuzziness.bw

########
computeMatrix scale-regions -S data/Irregular_Regions.bw \
     data/Occupancy.bw data/Shifts.bw data/Fuzziness.bw \
      -R $Anno/genes_sacCer3.gtf \
      -o "data/computeMatrix2plot_reg.txt.gz" \
      --outFileNameMatrix "data/computeMatrix2txt_reg.txt.gz" \
      --outFileSortedRegions "data/computeMatrix_geneList_reg.bed" \
      --averageTypeBins sum \
      -b 750 -a 750 --smartLabels -p 10

plotProfile -m "data/computeMatrix2plot_reg.txt.gz" \
     -out 'Default_Reg_Chen.png' \
     --averageType sum \
     --outFileNameData 'values_Profile_Reg.txt' \
     --samplesLabel 'regularity_changes' 'occupancy' 'shifts' 'fuzziness' --perGroup
