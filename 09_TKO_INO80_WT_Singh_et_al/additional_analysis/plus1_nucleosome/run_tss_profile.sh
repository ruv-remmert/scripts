#!/usr/bin/env bash

# Usage: ./run_tss_profile.sh <TSS_annotation.bed> <output_prefix> <bw1> <bw2> ... <bwN>
# Example: ./run_tss_profile.sh TSS_nuc_regSorted.bed TSS_out /path/to/PROFILER_out/RUN/01_NUCLEOSOME_PROFILE/*_monoNucs_profile.bw

TSS_ANNOTATION="$1"
OUTPUT_PREFIX="$2"
shift 2

BW_FILES=("$@")

# Step 1: Compute matrix around TSS
computeMatrix reference-point \
    -S "${BW_FILES[@]}" \
    -R "$TSS_ANNOTATION" \
    --referencePoint TSS \
    -o "${OUTPUT_PREFIX}_matrix.gz" \
    -b 1500 -a 1500 --smartLabels -p 4

# Step 2: Plot profile and emit per-group data table for R scripts
plotProfile \
    -m "${OUTPUT_PREFIX}_matrix.gz" \
    -out "${OUTPUT_PREFIX}_profile.pdf" \
    --outFileNameData "${OUTPUT_PREFIX}_profile.tab" \
    --perGroup \
    --plotHeight 15 \
    --plotWidth 24

# Step 3: Custom R plots
Rscript make_TSS_plots.R "${OUTPUT_PREFIX}_profile.tab" "${OUTPUT_PREFIX}_plot.pdf"
Rscript make_tss_profiles_individual.R "${OUTPUT_PREFIX}_profile.tab" "${OUTPUT_PREFIX}_individual_plot.pdf"

echo "TSS profile and plots generated: ${OUTPUT_PREFIX}_profile.pdf"
