#!/bin/bash

############################################################################
# Standalone Script: Create TSS Profile Plots
# Purpose: Generate TSS nucleosome profile plots from BigWig files and 
#          a custom BED file (regions of interest) without Nextflow
# 
# Usage: ./create_TSS_plots_standalone.sh <BED_file> <output_dir> [options]
############################################################################

set -e

# Default parameters
BIGWIG_DIR="/media/linuxmac/Storage2/scripts/02_example_yeast_remodeler/PROFILER_out/RUN/01_NUCLEOSOME_PROFILE"
R_SCRIPT="/media/linuxmac/Storage2/scripts/02_example_yeast_remodeler/add_INSPECTOR_TSS_scripts_ruven/add_INSPECTOR_TSS_scripts_ruven/profile_monoNucs_combined.R"
NCPUS=4
BEFORE_TSS=1500
AFTER_TSS=1500

# Print usage
usage() {
    echo "Usage: $0 <BED_file> <output_dir> [options]"
    echo ""
    echo "Required arguments:"
    echo "  <BED_file>       Path to BED file with regions of interest (TSS coordinates)"
    echo "  <output_dir>     Directory to save output plots and matrices"
    echo ""
    echo "Optional arguments:"
    echo "  --bigwig-dir DIR     Directory containing BigWig files (default: PROFILER_out/RUN/01_NUCLEOSOME_PROFILE)"
    echo "  --ncpus N            Number of CPUs to use (default: 4)"
    echo "  --before-tss BP      Base pairs before TSS (default: 1500)"
    echo "  --after-tss BP       Base pairs after TSS (default: 1500)"
    echo "  --r-script PATH      Path to profile_monoNucs_combined.R script"
    echo ""
    echo "Example:"
    echo "  $0 my_regions.bed ./TSS_output --bigwig-dir ./bigwig_files --ncpus 8"
    exit 1
}

# Parse arguments
if [ $# -lt 2 ]; then
    usage
fi

BED_FILE="$1"
OUTPUT_DIR="$2"
shift 2

# Parse optional arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --bigwig-dir)
            BIGWIG_DIR="$2"
            shift 2
            ;;
        --ncpus)
            NCPUS="$2"
            shift 2
            ;;
        --before-tss)
            BEFORE_TSS="$2"
            shift 2
            ;;
        --after-tss)
            AFTER_TSS="$2"
            shift 2
            ;;
        --r-script)
            R_SCRIPT="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate inputs
if [ ! -f "$BED_FILE" ]; then
    echo "ERROR: BED file not found: $BED_FILE"
    exit 1
fi

if [ ! -d "$BIGWIG_DIR" ]; then
    echo "ERROR: BigWig directory not found: $BIGWIG_DIR"
    exit 1
fi

if [ ! -f "$R_SCRIPT" ]; then
    echo "ERROR: R script not found: $R_SCRIPT"
    exit 1
fi
BED_FILE="$(realpath "/media/linuxmac/Storage2/scripts/02_example_yeast_remodeler/add_INSPECTOR_TSS_scripts_ruven/add_INSPECTOR_TSS_scripts_ruven/TES/TES_nuc_regSorted.bed")"
# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "=========================================="
echo "TSS Profile Plot Generation (Standalone)"
echo "=========================================="
echo "BED file:       $BED_FILE"
echo "BigWig dir:     $BIGWIG_DIR"
echo "Output dir:     $OUTPUT_DIR"
echo "R script:       $R_SCRIPT"
echo "CPUs:           $NCPUS"
echo "Window:         ${BEFORE_TSS}bp upstream, ${AFTER_TSS}bp downstream"
echo "=========================================="
echo ""

# Step 1: Collect BigWig files
echo "[1/3] Collecting BigWig files..."
BW_FILES=""
for bw in "$BIGWIG_DIR"/*_monoNucs_profile.bw; do
    if [ -f "$bw" ]; then
        BW_FILES="$BW_FILES $bw"
        echo "      Found: $(basename $bw)"
    fi
done

if [ -z "$BW_FILES" ]; then
    echo "ERROR: No BigWig files found in $BIGWIG_DIR"
    exit 1
fi

echo "      Total BigWig files: $(echo $BW_FILES | wc -w)"
echo ""

# Step 2: Run computeMatrix
echo "[2/3] Running computeMatrix..."
echo "      Command: computeMatrix reference-point -S $BW_FILES -R $BED_FILE --referencePoint TSS -b $BEFORE_TSS -a $AFTER_TSS"
computeMatrix reference-point \
    -S $BW_FILES \
    -R "$BED_FILE" \
    --referencePoint TSS \
    -o "computeMatrix.txt.gz" \
    -b $BEFORE_TSS \
    -a $AFTER_TSS \
    --smartLabels \
    -p $NCPUS

if [ ! -f "computeMatrix.txt.gz" ]; then
    echo "ERROR: computeMatrix failed to generate output"
    exit 1
fi
echo "      ✓ Matrix saved: computeMatrix.txt.gz"
echo ""

# Step 3: Run plotProfile
echo "[3/3] Extracting profile data with plotProfile..."
plotProfile \
    -m "computeMatrix.txt.gz" \
    -out "DefaultHeatmap.png" \
    --outFileNameData "values_Profile.txt"

if [ ! -f "values_Profile.txt" ]; then
    echo "ERROR: plotProfile failed to generate output"
    exit 1
fi
echo "      ✓ Heatmap saved: DefaultHeatmap.png"
echo "      ✓ Profile data saved: values_Profile.txt"
echo ""

# Step 4: Run R script to create publication-quality plot
echo "[4/4] Creating publication-quality TSS plot with R..."
Rscript "$R_SCRIPT" "values_Profile.txt"

if [ ! -f "profile_monoNucs.pdf" ]; then
    echo "WARNING: R script may not have generated PDF. Check values_Profile.txt format."
else
    echo "      ✓ TSS plot saved: profile_monoNucs.pdf"
fi

echo ""
echo "=========================================="
echo "✓ TSS profile plots created successfully!"
echo "=========================================="
echo ""
echo "Output files in $OUTPUT_DIR:"
echo "  - computeMatrix.txt.gz         (compressed matrix data)"
echo "  - DefaultHeatmap.png           (deepTools heatmap)"
echo "  - values_Profile.txt           (profile data in text format)"
echo "  - profile_monoNucs.pdf         (publication-quality plot)"
echo ""

