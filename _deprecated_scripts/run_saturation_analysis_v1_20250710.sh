#!/bin/bash
# Isoform Saturation Curve Analysis Runner
# Version: v1_20250710
# Date: 2025-07-10

# Check if input file is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_file> [output_prefix]"
    echo "Example: $0 your_data.txt results"
    exit 1
fi

INPUT_FILE=$1
OUTPUT_PREFIX=${2:-"saturation_results"}

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    exit 1
fi

echo "=== Running Isoform Saturation Analysis ==="
echo "Input file: $INPUT_FILE"
echo "Output prefix: $OUTPUT_PREFIX"
echo ""

# Run the analysis
python3 saturation_curve_script_v1_20250710.py "$INPUT_FILE" "$OUTPUT_PREFIX"

echo ""
echo "=== Analysis Complete ==="
echo "Check the following output files:"
echo "- ${OUTPUT_PREFIX}_saturation_curve.png (plot)"
echo "- ${OUTPUT_PREFIX}_saturation_curve.pdf (plot)"
echo "- ${OUTPUT_PREFIX}_saturation_results.txt (data)" 