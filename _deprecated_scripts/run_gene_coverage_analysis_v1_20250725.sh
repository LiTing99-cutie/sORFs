#!/bin/bash
# Transcript Coverage Analysis Runner v1_20250725
# 运行转录本覆盖分析的脚本

# 设置输入文件路径
INPUT_FILE="../processed/batch_1/geneCounts/gene.counts.20250725.txt"

# 设置输出目录
OUTPUT_DIR="../processed/batch_1/gene_coverage_analysis"

# 设置参数
PERMUTATIONS=50  # 减少排列次数以加快运行速度
MIN_READS=1      # 最小reads阈值

echo "Starting Transcript Coverage Analysis..."
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Permutations: $PERMUTATIONS"
echo "Min reads threshold: $MIN_READS"
echo ""

# 运行分析脚本
python3 transcript_coverage_analysis_v1_20250725.py \
    --input "$INPUT_FILE" \
    --output "$OUTPUT_DIR" \
    --permutations $PERMUTATIONS \
    --min-reads $MIN_READS

echo ""
echo "Analysis completed!"
echo "Results saved to: $OUTPUT_DIR" 

Rscript plot_transcript_coverage_simple_v1_20250725.R $OUTPUT_DIR/transcript_coverage_results.csv \
    $OUTPUT_DIR/transcript_coverage_curves.pdf