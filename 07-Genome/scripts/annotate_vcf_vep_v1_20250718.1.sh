#!/bin/bash

# VEP注释脚本 v1_20250710
# 使用VEP对过滤后的VCF文件进行功能注释

set -e

# 激活biotools虚拟环境
source activate biotools

# 设置输入输出路径
INPUT_VCF="/home/user/data3/lit/project/sORFs/07-Genome/processed/vcf_filter/human_brain_21pcw_filtered_pass.vcf"
# INPUT_VCF="/home/user/data3/lit/project/sORFs/07-Genome/processed/vcf_filter/subset_1000.vcf"
OUTPUT_DIR="/home/user/data3/lit/project/sORFs/07-Genome/results/vep"
SAMPLE_NAME="human_brain_21pcw"
# SAMPLE_NAME="subset_1000"
# 创建输出目录
mkdir -p ${OUTPUT_DIR}

echo "[$(date '+%Y-%m-%d %H:%M:%S')] 开始VEP注释流程..."
echo "输入文件: ${INPUT_VCF}"
echo "输出目录: ${OUTPUT_DIR}"

# 检查输入文件
if [ ! -f "${INPUT_VCF}" ]; then
    echo "错误: 输入文件不存在: ${INPUT_VCF}"
    exit 1
fi

# 使用VEP进行注释
vep \
    --input_file ${INPUT_VCF} \
    --output_file ${OUTPUT_DIR}/${SAMPLE_NAME}_vep_annotated.txt \
    --cache \
    --format vcf \
    --offline \
    --dir_cache /home/user/data3/lit/project/sORFs/07-Genome/vep_cache \
    --species homo_sapiens \
    --assembly GRCh38 \
    --force_overwrite \
    --fork 10 \
    --stats_text \
    --stats_html \
    --stats_file ${OUTPUT_DIR}/${SAMPLE_NAME}_vep_stats.html \
    --verbose

echo "[$(date '+%Y-%m-%d %H:%M:%S')] VEP注释完成"