#!/bin/bash

# VCF过滤脚本 v1_20250710
# 基于GATK推荐标准进行SNP和INDEL过滤

set -e

# 设置输入输出路径
INPUT_VCF="/home/user/data3/lit/project/sORFs/07-Genome/results/human_brain_21pcw.vcf.gz"
OUTPUT_DIR="/home/user/data3/lit/project/sORFs/07-Genome/processed/vcf_filter"
mkdir -p ${OUTPUT_DIR}
SAMPLE_NAME="human_brain_21pcw"

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

echo "开始VCF过滤流程..."
echo "输入文件: ${INPUT_VCF}"
echo "输出目录: ${OUTPUT_DIR}"

# 1. 分离SNP和INDEL
echo "步骤1: 分离SNP和INDEL变异..."

gatk SelectVariants \
    -V ${INPUT_VCF} \
    -select-type SNP \
    -O ${OUTPUT_DIR}/${SAMPLE_NAME}_snps.vcf.gz

gatk SelectVariants \
    -V ${INPUT_VCF} \
    -select-type INDEL \
    -O ${OUTPUT_DIR}/${SAMPLE_NAME}_indels.vcf.gz

# 2. 过滤SNP（包含QUAL过滤）
echo "步骤2: 过滤SNP变异..."

gatk VariantFiltration \
    -V ${OUTPUT_DIR}/${SAMPLE_NAME}_snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${OUTPUT_DIR}/${SAMPLE_NAME}_snps_filtered.vcf.gz

# 3. 过滤INDEL（包含QUAL过滤）
echo "步骤3: 过滤INDEL变异..."

gatk VariantFiltration \
    -V ${OUTPUT_DIR}/${SAMPLE_NAME}_indels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O ${OUTPUT_DIR}/${SAMPLE_NAME}_indels_filtered.vcf.gz

# 4. 合并过滤后的SNP和INDEL
echo "步骤4: 合并过滤后的SNP和INDEL..."

gatk MergeVcfs \
    -I ${OUTPUT_DIR}/${SAMPLE_NAME}_snps_filtered.vcf.gz \
    -I ${OUTPUT_DIR}/${SAMPLE_NAME}_indels_filtered.vcf.gz \
    -O ${OUTPUT_DIR}/${SAMPLE_NAME}_filtered.vcf.gz

# 5. 提取PASS变异（去除被过滤的变异）
echo "步骤5: 提取PASS变异..."

gatk SelectVariants \
    -V ${OUTPUT_DIR}/${SAMPLE_NAME}_filtered.vcf.gz \
    --exclude-filtered \
    -O ${OUTPUT_DIR}/${SAMPLE_NAME}_filtered_pass.vcf.gz

# 6. 生成统计报告
echo "步骤6: 生成过滤统计报告..."

echo "=== VCF过滤统计报告 ===" > ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt
echo "生成时间: $(date)" >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt
echo "" >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt

# 原始变异数量
echo "原始变异总数:" >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt
bcftools view -H ${INPUT_VCF} | wc -l >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt

# SNP数量
echo "SNP数量:" >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt
bcftools view -H ${OUTPUT_DIR}/${SAMPLE_NAME}_snps.vcf.gz | wc -l >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt

# INDEL数量
echo "INDEL数量:" >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt
bcftools view -H ${OUTPUT_DIR}/${SAMPLE_NAME}_indels.vcf.gz | wc -l >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt

# 过滤后PASS变异数量
echo "过滤后PASS变异数量:" >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt
bcftools view -H ${OUTPUT_DIR}/${SAMPLE_NAME}_filtered_pass.vcf.gz | wc -l >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt

# 7. 生成详细的过滤统计
echo "步骤7: 生成详细过滤统计..."

# 查看各过滤条件的影响
echo "=== 各过滤条件统计 ===" >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt
bcftools query -f '[%FILTER]\n' ${OUTPUT_DIR}/${SAMPLE_NAME}_filtered.vcf.gz | sort | uniq -c >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt

# 8. 创建索引文件
echo "步骤8: 创建索引文件..."

tabix -p vcf ${OUTPUT_DIR}/${SAMPLE_NAME}_filtered_pass.vcf.gz

echo "VCF过滤完成！"
echo "最终过滤后的文件: ${OUTPUT_DIR}/${SAMPLE_NAME}_filtered_pass.vcf.gz"
echo "过滤报告: ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report.txt"

# 9. 显示过滤结果摘要
echo ""
echo "=== 过滤结果摘要 ==="
echo "原始变异总数: $(bcftools view -H ${INPUT_VCF} | wc -l)"
echo "SNP数量: $(bcftools view -H ${OUTPUT_DIR}/${SAMPLE_NAME}_snps.vcf.gz | wc -l)"
echo "INDEL数量: $(bcftools view -H ${OUTPUT_DIR}/${SAMPLE_NAME}_indels.vcf.gz | wc -l)"
echo "过滤后PASS变异数量: $(bcftools view -H ${OUTPUT_DIR}/${SAMPLE_NAME}_filtered_pass.vcf.gz | wc -l)" 