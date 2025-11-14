#!/bin/bash

# 输入文件
vcf=/home/user/data3/lit/project/sORFs/11-denovo-list/public_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_AllVar_QC_metrics.vcf.gz
CDS_bed=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_ortholog_extraction/noncano.canonical.intergenic_orfs.lnc_orfs.bed
out_dir=$(realpath ../processed/pnps/orf_regions)
mkdir -p $out_dir
# 检查VCF是否有索引
# if [ ! -f "${vcf}.tbi" ]; then
#     echo "ERROR: VCF index not found. Please create index first:"
#     echo "  tabix -p vcf $vcf"
#     exit 1
# fi

# 检查BED文件格式和染色体命名
echo "Checking BED file..."
head -3 $CDS_bed
echo ""

# 检查VCF染色体命名
echo "Checking VCF chromosomes..."
bcftools view -H $vcf | head -1000 | cut -f1 | sort -u | head -5
echo ""

# 一步到位
# -R 只保留BED文件中指定区域的variants
# -z 输出格式: z=bgzip压缩VCF, v=未压缩, b=BCF
bcftools view -R $CDS_bed $vcf -O z -o $out_dir/orf_regions.vcf.gz && \
bcftools index -c $out_dir/orf_regions.vcf.gz && \
echo "Extracted $(bcftools view -H $out_dir/orf_regions.vcf.gz | wc -l) variants"