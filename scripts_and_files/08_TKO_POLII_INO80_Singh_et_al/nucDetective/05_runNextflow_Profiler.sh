#!/usr/bin/env bash
MY_PWD=$1
cd $MY_PWD


# Run Profiler workflow for MNase-seq samples

nextflow run nucDetective \
--analysis 'profiler' \
--csvInput $MY_PWD'/data/samples_profiler.csv' \
--outDir PROFILER_out \
--genomeIdx $MY_PWD'/data/Bowtie2Index/genome' \
--genomeSize 12157105 \
--publishBamFlt \
--blacklist $MY_PWD'/data//blacklist.bed' \
--peaks_used 3 \
--TSS $MY_PWD'/data/genes_sacCer3.gtf' \
-w ./work_profiler 