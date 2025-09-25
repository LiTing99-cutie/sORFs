#!/bin/bash
# v1_20240710
set -e

# 激活biotools环境
source activate biotools

# 路径和文件名定义
SAMPLE=human_brain_21pcw
REF=/home/user/data/lit/database/public/genome/hg38/hg38.fa
R1=../rawdata/L1EJF1602305-p21_Gen2Seq.R1.raw.fastq.gz
R2=../rawdata/L1EJF1602305-p21_Gen2Seq.R2.raw.fastq.gz
PROCESSED=../processed
RESULTS=../results
THREADS=30

mkdir -p $PROCESSED $RESULTS

# 1. 参考基因组索引（如已存在可跳过）
# if [ ! -f "$REF.bwt" ]; then
#     bwa index $REF
# fi
# if [ ! -f "$REF.fai" ]; then
#     samtools faidx $REF
# fi
# if [ ! -f "${REF%.fa}.dict" ]; then
#     gatk CreateSequenceDictionary -R $REF -O ${REF%.fa}.dict
# fi

# 2. 比对并直接转换为排序后的BAM，避免生成中间SAM文件
# bwa mem -M -t $THREADS -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:illumina\tLB:$SAMPLE" $REF $R1 $R2 \
#   | samtools sort -@ $THREADS -o $PROCESSED/$SAMPLE.sorted.bam
# 4. 标记重复
# gatk MarkDuplicates -I $PROCESSED/$SAMPLE.sorted.bam -O $PROCESSED/$SAMPLE.markdup.bam -M $PROCESSED/$SAMPLE.metrics.txt
# samtools index -@ $THREADS $PROCESSED/$SAMPLE.markdup.bam

java -Xmx250G \
     -Djava.io.tmpdir=../tmp \
     -jar /home/user/data2/lit/software/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
     MarkDuplicates \
     -I $PROCESSED/$SAMPLE.sorted.bam \
     -O $PROCESSED/$SAMPLE.markdup.bam \
     -M $PROCESSED/$SAMPLE.metrics.txt \
     --TMP_DIR ../tmp

# java -Xmx250G -Djava.io.tmpdir=../tmp -jar /home/user/data2/lit/software/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar MarkDuplicates -h
