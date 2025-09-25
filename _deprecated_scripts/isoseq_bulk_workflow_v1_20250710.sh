#!/bin/bash
set -e

# Iso-Seq3 bulk workflow v1_20250710
# 参考官方文档：https://isoseq.how/getting-started.html
# 输入文件：经过lima处理的ccs BAM
# 输出目录：../processed

# 激活biotools环境
source activate biotools

# 路径定义
INPUT_BAM=../rawdata/r84130_250703_001_1_A01/p21-IsoSeq.Iso_bc02.bcM0004.ISO.bam
PRIMER_FA=../rawdata/r84130_250703_001_1_A01/primer.fa
GENOME_FA=/home/user/data/lit/database/public/genome/hg38/hg38.fa
GTF_ANN=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
PROCESSED=../processed
THREADS=30

mkdir -p $PROCESSED

# 1. refine：去除polyA尾和concatemers
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start isoseq3 refine"
isoseq3 refine $INPUT_BAM $PRIMER_FA $PROCESSED/p21-IsoSeq.flnc.bam --require-polya --num-threads $THREADS
# 2. cluster2：聚类获得高置信度转录本
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start isoseq3 cluster2"
isoseq3 cluster2 $PROCESSED/p21-IsoSeq.flnc.bam $PROCESSED/p21-IsoSeq.clustered.bam --num-threads $THREADS

# 3. pbmm2：比对到基因组
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start pbmm2 align"
pbmm2 align $GENOME_FA $PROCESSED/p21-IsoSeq.clustered.bam $PROCESSED/p21-IsoSeq.mapped.bam --preset ISOSEQ --sort --log-level INFO --num-threads $THREADS

# 4. isoseq collapse：冗余转录本合并
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start isoseq3 collapse"
isoseq3 collapse --do-not-collapse-extra-5exons $PROCESSED/p21-IsoSeq.mapped.bam $PROCESSED/p21-IsoSeq.flnc.bam $PROCESSED/p21-IsoSeq.collapsed.gff --gtf $GTF_ANN

# 5. pigeon classify：转录本注释分类
pigeon prepare
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start pigeon classify"
pigeon classify $PROCESSED/p21-IsoSeq.collapsed.gff $GTF_ANN $GENOME_FA --fl $PROCESSED/flnc_count.txt --num-threads $THREADS

# 6. pigeon filter：过滤潜在假阳性转录本
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start pigeon filter"
# 修正：pigeon filter 需要classification.txt和--isoforms参数
pigeon filter $PROCESSED/p21-IsoSeq.collapsed.classification.txt --isoforms $PROCESSED/p21-IsoSeq.collapsed.gff --num-threads $THREADS

echo "Iso-Seq3 bulk流程全部完成！" 