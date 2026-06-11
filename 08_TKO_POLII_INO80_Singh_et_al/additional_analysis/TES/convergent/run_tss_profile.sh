#!/usr/bin/env bash

# Usage: ./run_tes_profile_TKO.sh <TES_annotation.bed> <output_prefix> <bw1> <bw2> ... <bwN>
# Example: ./run_tes_profile_TKO.sh TES_nuc_regSorted.bed TES_out /path/to/INSPECTOR_out/RUN/01_BIGWIG_PROFILES/*.bw

TES_ANNOTATION="$1"
OUTPUT_PREFIX="$2"
shift 2

BW_FILES=("$@")

# Step 1: Compute matrix around TES
computeMatrix reference-point \
    -S "${BW_FILES[@]}" \
    -R "$TES_ANNOTATION" \
    --referencePoint TES \
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
Rscript make_TES_plots_TKO.R "${OUTPUT_PREFIX}_profile.tab" "${OUTPUT_PREFIX}_plot.pdf"
Rscript make_tes_profiles_individual_TKO.R "${OUTPUT_PREFIX}_profile.tab" "${OUTPUT_PREFIX}_individual_plot.pdf"

echo "TES profile and plots generated: ${OUTPUT_PREFIX}_profile.pdf"
