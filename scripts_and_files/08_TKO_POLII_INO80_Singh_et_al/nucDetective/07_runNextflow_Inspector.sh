#!/usr/bin/env bash
MY_PWD=$1
cd $MY_PWD


# Run Inspector workflow using output from Profiler workflow

nextflow run nucDetective \
--analysis 'inspector' \
--csvInput $MY_PWD'/data/samples_inspector.csv' \
--outDir INSPECTOR_out \
--genomeSize 12157105 \
--normalize_profiles false \
--chrSizes $MY_PWD'/PROFILER_out/RUN/chrom_Sizes.txt' \
--TSS $MY_PWD'/data/genes_sacCer3.gtf' \
-w ./work_inspector -resume
