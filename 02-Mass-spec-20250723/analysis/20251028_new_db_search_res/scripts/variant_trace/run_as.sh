#!/bin/bash

# ORF分析运行脚本
# 使用方法：修改路径后运行此脚本

# 输入文件路径
VCF_FILE="/home/user/data3/lit/project/sORFs/07-Genome/results/custom_fa/hom.vcf.gz"
ORF_GTF="/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/results/augment_orf_table/sub.gtf"
ANNOTATION_GTF="/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf"
CLASSIFICATION_FILE="/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/processed/classify/collapsed_classification.filtered_lite_classification.txt"

# 参考基因组文件（用于突变类型判断，需要您提供路径）
# REF_GENOME="/home/user/data/lit/database/public/genome/hg38/hg38.fa"  # 请修改为实际路径

# 输出文件前缀
OUTPUT_PREFIX="/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/processed/variant_trace/"

echo "=========================================="
echo "开始ORF分析"
echo "=========================================="

# Step 1: 运行主分析（SNP重叠和可变剪切来源）
echo "Step 1: 运行主分析..."
python3 orf_as.py \
    --vcf "${VCF_FILE}" \
    --orf-gtf "${ORF_GTF}" \
    --annotation-gtf "${ANNOTATION_GTF}" \
    --classification "${CLASSIFICATION_FILE}" \
    --output-prefix "${OUTPUT_PREFIX}"

# 检查是否成功
if [ $? -ne 0 ]; then
    echo "错误：主分析失败"
    exit 1
fi

echo ""
echo "=========================================="
echo "分析完成！"
echo "输出文件："
echo "  - ${OUTPUT_PREFIX}_snp_overlap.tsv : SNP-ORF重叠分析结果"
echo "  - ${OUTPUT_PREFIX}_as_source.tsv : 可变剪切来源分析结果"
if [ -f "${ANNOTATED_OUTPUT}" ]; then
    echo "  - ${ANNOTATED_OUTPUT} : 带突变类型注释的SNP分析结果"
fi
echo "=========================================="