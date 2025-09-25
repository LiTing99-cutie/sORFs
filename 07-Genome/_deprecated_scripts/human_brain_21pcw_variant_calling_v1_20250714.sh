#!/bin/bash
# v1_20250712
set -e

# 激活biotools环境
source activate biotools

# 路径和文件名定义
SAMPLE=human_brain_21pcw
REF=/home/user/data/lit/database/public/genome/hg38/hg38.fa
PROCESSED=../processed
RESULTS=../results
THREADS=30

mkdir -p $PROCESSED $RESULTS ../tmp

# 1. 基因组重校准（BaseRecalibrator & ApplyBQSR）
# 需要提供已知变异数据库（如dbSNP、1000G等），此处以常见路径为例
DBSNP=/home/user/data3/lit/project/sORFs/07-Genome/public_vcf/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
KNOWN_SITES=$DBSNP

# 1.1 生成初始 recal table
echo "[$(date)] Starting BaseRecalibrator (initial)..."
java -Xmx250G \
     -Djava.io.tmpdir=../tmp \
     -jar /home/user/data2/lit/software/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
     BaseRecalibrator \
     -R $REF \
     -I $PROCESSED/$SAMPLE.markdup.bam \
     --known-sites $KNOWN_SITES \
     -O $PROCESSED/$SAMPLE.recal.table \
     --tmp-dir ../tmp

# 1.2 应用BQSR
echo "[$(date)] Starting ApplyBQSR..."
java -Xmx250G \
     -Djava.io.tmpdir=../tmp \
     -jar /home/user/data2/lit/software/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
     ApplyBQSR \
     -R $REF \
     -I $PROCESSED/$SAMPLE.markdup.bam \
     --bqsr-recal-file $PROCESSED/$SAMPLE.recal.table \
     -O $PROCESSED/$SAMPLE.recal.bam \
     --tmp-dir ../tmp

# 1.3 生成 post-recal table
echo "[$(date)] Starting BaseRecalibrator (post)..."
java -Xmx250G \
     -Djava.io.tmpdir=../tmp \
     -jar /home/user/data2/lit/software/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
     BaseRecalibrator \
     -R $REF \
     -I $PROCESSED/$SAMPLE.recal.bam \
     --known-sites $KNOWN_SITES \
     -O $PROCESSED/$SAMPLE.post_recal.table \
     --tmp-dir ../tmp

# 1.4 分析前后差异
echo "[$(date)] Starting AnalyzeCovariates..."
java -Xmx250G \
     -Djava.io.tmpdir=../tmp \
     -jar /home/user/data2/lit/software/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
     AnalyzeCovariates \
     -before $PROCESSED/$SAMPLE.recal.table \
     -after $PROCESSED/$SAMPLE.post_recal.table \
     -plots $PROCESSED/$SAMPLE.AnalyzeCovariates.pdf \
     --tmp-dir ../tmp

# 2. GATK变异检测
# 2.1 HaplotypeCaller生成gVCF
echo "[$(date)] Starting HaplotypeCaller..."
java -Xmx250G \
     -Djava.io.tmpdir=../tmp \
     -jar /home/user/data2/lit/software/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
     HaplotypeCaller \
     -R $REF \
     -I $PROCESSED/$SAMPLE.recal.bam \
     -O $PROCESSED/$SAMPLE.g.vcf.gz \
     -ERC GVCF \
     --tmp-dir ../tmp

# 2.2 GenotypeGVCFs生成VCF
echo "[$(date)] Starting GenotypeGVCFs..."
java -Xmx250G \
     -Djava.io.tmpdir=../tmp \
     -jar /home/user/data2/lit/software/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
     GenotypeGVCFs \
     -R $REF \
     -V $PROCESSED/$SAMPLE.g.vcf.gz \
     -O $RESULTS/$SAMPLE.vcf.gz \
     --tmp-dir ../tmp

echo "[$(date)] Base quality recalibration and variant calling finished!" 