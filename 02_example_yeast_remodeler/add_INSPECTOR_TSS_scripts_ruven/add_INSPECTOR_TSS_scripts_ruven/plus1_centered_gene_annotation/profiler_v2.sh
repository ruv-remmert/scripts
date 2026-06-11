#!/bin/bash

############################################################################
# Standalone Script: Create TSS Profile Plots
# Purpose: Generate TSS nucleosome profile plots from BigWig files and 
#          a custom BED file (regions of interest) without Nextflow
#
# Produces:
#   - Combined plot (all mutants)
#   - Individual pairwise plots: WT vs ISW1, WT vs CHD1, WT vs ISW2,
#     WT vs triple mutant (ISW1-ISW2-CHD1)
#   - Each plot output as both PDF and JPG
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

# Check for ImageMagick (convert) for PDF -> JPG conversion
if ! command -v convert &> /dev/null; then
    echo "ERROR: ImageMagick 'convert' not found. Please install ImageMagick."
    exit 1
fi

BED_FILE="$(realpath "/media/linuxmac/Storage2/scripts/02_example_yeast_remodeler/add_INSPECTOR_TSS_scripts_ruven/add_INSPECTOR_TSS_scripts_ruven/plus1_centered_gene_annotation/TSS_nuc_regSorted.bed")"

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

############################################################################
# Helper: convert a PDF to JPG using ImageMagick
# Usage: pdf_to_jpg <input.pdf> <output.jpg>
############################################################################
pdf_to_jpg() {
    local pdf_in="$1"
    local jpg_out="$2"
    # -density 300: 300 dpi for publication quality
    # -flatten:     merge alpha channel (avoids black backgrounds)
    # [0]:          take only the first page
    convert -density 300 "${pdf_in}[0]" -flatten "$jpg_out"
    echo "      ✓ JPG saved: $(basename $jpg_out)"
}

############################################################################
# Helper: run plotProfile for a specific set of sample labels, then run
#         the R script and convert the resulting PDF to JPG.
#
# Usage: run_comparison <tag> <sample_label1> [<sample_label2> ...]
#
#   <tag>           Short identifier used in output filenames
#                   (e.g. "WT_vs_ISW1")
#   <sample_label>  One or more sample labels exactly as they appear in
#                   the computeMatrix header (--smartLabels uses the
#                   BigWig basename without extension).
############################################################################
run_comparison() {
    local tag="$1"
    shift
    local samples=("$@")   # remaining args are sample labels

    echo "  --> Comparison: $tag"
    echo "      Samples: ${samples[*]}"

    local out_png="${tag}_Heatmap.png"
    local out_data="${tag}_values_Profile.txt"
    local out_pdf="${tag}_profile_monoNucs.pdf"
    local out_jpg="${tag}_profile_monoNucs.jpg"

    # Extract profile data for this subset of samples from the shared matrix
    plotProfile \
        -m "computeMatrix.txt.gz" \
        -out "$out_png" \
        --outFileNameData "$out_data" \
        --samples "${samples[@]}"

    if [ ! -f "$out_data" ]; then
        echo "      WARNING: plotProfile did not produce $out_data — skipping R step."
        return 1
    fi
    echo "      ✓ Profile data saved: $out_data"

    # Run R script to create publication-quality PDF
    Rscript "$R_SCRIPT" "$out_data"

    # The R script always writes to 'profile_monoNucs.pdf' in the working dir;
    # rename it to the comparison-specific name immediately after.
    if [ -f "profile_monoNucs.pdf" ]; then
        mv "profile_monoNucs.pdf" "$out_pdf"
        echo "      ✓ PDF saved: $out_pdf"
        pdf_to_jpg "$out_pdf" "$out_jpg"
    else
        echo "      WARNING: R script did not produce profile_monoNucs.pdf for $tag"
    fi

    echo ""
}

# ---------------------------------------------------------------------------
# Step 1: Collect ALL BigWig files
# ---------------------------------------------------------------------------
echo "[1/3] Collecting BigWig files..."
BW_FILES=()
for bw in "$BIGWIG_DIR"/*_monoNucs_profile.bw; do
    if [ -f "$bw" ]; then
        BW_FILES+=("$bw")
        echo "      Found: $(basename $bw)"
    fi
done

if [ ${#BW_FILES[@]} -eq 0 ]; then
    echo "ERROR: No BigWig files found in $BIGWIG_DIR"
    exit 1
fi

echo "      Total BigWig files: ${#BW_FILES[@]}"
echo ""

# ---------------------------------------------------------------------------
# Step 2: Run computeMatrix ONCE for all samples
# ---------------------------------------------------------------------------
echo "[2/3] Running computeMatrix (combined, all samples)..."
computeMatrix reference-point \
    -S "${BW_FILES[@]}" \
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

# ---------------------------------------------------------------------------
# Step 3: Extract sample labels from the matrix header
#         deepTools stores them as a JSON field "@type" -> "parameters" ->
#         "sample_labels" in the gzipped header line.
# ---------------------------------------------------------------------------
echo "[3/3] Extracting sample labels from matrix..."
SAMPLE_LABELS_JSON=$(zcat "computeMatrix.txt.gz" | head -1 | python3 -c "
import sys, json
header = sys.stdin.read().lstrip('@')
data = json.loads(header)
labels = data.get('sample_labels', [])
print('\n'.join(labels))
")

if [ -z "$SAMPLE_LABELS_JSON" ]; then
    echo "ERROR: Could not extract sample labels from computeMatrix.txt.gz"
    exit 1
fi

echo "      Detected sample labels:"
echo "$SAMPLE_LABELS_JSON" | while read -r lbl; do echo "        - $lbl"; done
echo ""

# Build arrays of matching labels (case-insensitive).
# Triple mutant is identified as any sample whose name contains ALL THREE of
# isw1, isw2, and chd1 — regardless of order or capitalisation.
# This check must come first so such samples are not claimed by the single-
# mutant buckets below.
WT_LABELS=()
ISW1_LABELS=()
CHD1_LABELS=()
ISW2_LABELS=()
TRIPLE_LABELS=()

while IFS= read -r label; do
    lower_label=$(echo "$label" | tr '[:upper:]' '[:lower:]')
    has_isw1=false; has_isw2=false; has_chd1=false
    echo "$lower_label" | grep -q 'isw1' && has_isw1=true
    echo "$lower_label" | grep -q 'isw2' && has_isw2=true
    echo "$lower_label" | grep -q 'chd1' && has_chd1=true

    if $has_isw1 && $has_isw2 && $has_chd1; then
        # Contains all three mutant names → triple mutant
        TRIPLE_LABELS+=("$label")
    elif $has_isw1; then
        ISW1_LABELS+=("$label")
    elif $has_isw2; then
        ISW2_LABELS+=("$label")
    elif $has_chd1; then
        CHD1_LABELS+=("$label")
    elif echo "$lower_label" | grep -qE '(^|[_\-\.])wt([_\-\.]|$)'; then
        WT_LABELS+=("$label")
    fi
done <<< "$SAMPLE_LABELS_JSON"

echo "      WT samples:     ${WT_LABELS[*]:-<none detected>}"
echo "      ISW1 samples:   ${ISW1_LABELS[*]:-<none detected>}"
echo "      CHD1 samples:   ${CHD1_LABELS[*]:-<none detected>}"
echo "      ISW2 samples:   ${ISW2_LABELS[*]:-<none detected>}"
echo "      Triple samples: ${TRIPLE_LABELS[*]:-<none detected>}"
echo ""

# ---------------------------------------------------------------------------
# Step 4: Combined plot (all samples — original behaviour)
# ---------------------------------------------------------------------------
echo "[4/5] Creating COMBINED plot (all samples)..."

plotProfile \
    -m "computeMatrix.txt.gz" \
    -out "DefaultHeatmap.png" \
    --outFileNameData "values_Profile.txt"

if [ ! -f "values_Profile.txt" ]; then
    echo "ERROR: plotProfile failed to generate values_Profile.txt"
    exit 1
fi
echo "      ✓ Heatmap saved: DefaultHeatmap.png"
echo "      ✓ Profile data saved: values_Profile.txt"

Rscript "$R_SCRIPT" "values_Profile.txt"

if [ -f "profile_monoNucs.pdf" ]; then
    mv "profile_monoNucs.pdf" "combined_profile_monoNucs.pdf"
    echo "      ✓ PDF saved: combined_profile_monoNucs.pdf"
    pdf_to_jpg "combined_profile_monoNucs.pdf" "combined_profile_monoNucs.jpg"
else
    echo "      WARNING: R script did not produce profile_monoNucs.pdf for combined plot"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 5: Individual pairwise comparisons (reusing the combined matrix)
# ---------------------------------------------------------------------------
echo "[5/5] Creating individual pairwise comparison plots..."
echo ""

if [ ${#WT_LABELS[@]} -eq 0 ]; then
    echo "WARNING: No WT samples detected — skipping all pairwise comparisons."
    echo "         Check that WT BigWig filenames contain 'wt' (case-insensitive)."
else
    # WT vs ISW1
    if [ ${#ISW1_LABELS[@]} -gt 0 ]; then
        run_comparison "WT_vs_ISW1" "${WT_LABELS[@]}" "${ISW1_LABELS[@]}"
    else
        echo "  --> Skipping WT vs ISW1: no ISW1 samples detected"
    fi

    # WT vs CHD1
    if [ ${#CHD1_LABELS[@]} -gt 0 ]; then
        run_comparison "WT_vs_CHD1" "${WT_LABELS[@]}" "${CHD1_LABELS[@]}"
    else
        echo "  --> Skipping WT vs CHD1: no CHD1 samples detected"
    fi

    # WT vs ISW2
    if [ ${#ISW2_LABELS[@]} -gt 0 ]; then
        run_comparison "WT_vs_ISW2" "${WT_LABELS[@]}" "${ISW2_LABELS[@]}"
    else
        echo "  --> Skipping WT vs ISW2: no ISW2 samples detected"
    fi

    # WT vs Triple mutant
    if [ ${#TRIPLE_LABELS[@]} -gt 0 ]; then
        run_comparison "WT_vs_triple" "${WT_LABELS[@]}" "${TRIPLE_LABELS[@]}"
    else
        echo "  --> Skipping WT vs triple: no triple mutant samples detected"
    fi
fi

echo ""
echo "=========================================="
echo "✓ TSS profile plots created successfully!"
echo "=========================================="
echo ""
echo "Output files in $OUTPUT_DIR:"
echo "  - computeMatrix.txt.gz                   (combined matrix, all samples)"
echo "  - DefaultHeatmap.png                     (deepTools heatmap, all samples)"
echo "  - values_Profile.txt                     (profile data, all samples)"
echo "  - combined_profile_monoNucs.pdf/jpg      (publication plot, all samples)"
echo "  - WT_vs_ISW1_profile_monoNucs.pdf/jpg    (WT vs ISW1)"
echo "  - WT_vs_CHD1_profile_monoNucs.pdf/jpg    (WT vs CHD1)"
echo "  - WT_vs_ISW2_profile_monoNucs.pdf/jpg    (WT vs ISW2)"
echo "  - WT_vs_triple_profile_monoNucs.pdf/jpg  (WT vs triple mutant)"
echo ""
