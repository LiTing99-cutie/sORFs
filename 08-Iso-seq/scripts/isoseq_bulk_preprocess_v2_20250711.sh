#!/bin/bash
set -e

# Iso-Seq3 bulk preprocess v1_20250710
# 步骤1-3：refine, cluster2, pbmm2 align
# 输入文件：经过lima处理的ccs BAM
# 输出目录：../processed

# 激活biotools环境
source activate biotools

# 路径定义
INPUT_BAM=../rawdata/r84130_250703_001_1_A01/p21-IsoSeq.Iso_bc02.bcM0004.ISO.bam
PRIMER_FA=../rawdata/r84130_250703_001_1_A01/primer.fa
GENOME_FA=/home/user/data/lit/database/public/genome/hg38/hg38.fa
PROCESSED=../processed
THREADS=30

mkdir -p $PROCESSED

# # 1. refine：去除polyA尾和concatemers
# echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start isoseq3 refine"
# isoseq3 refine $INPUT_BAM $PRIMER_FA $PROCESSED/p21-IsoSeq.flnc.bam --require-polya --num-threads $THREADS

# # 2. cluster2：聚类获得高置信度转录本
# echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start isoseq3 cluster2"
# isoseq3 cluster2 $PROCESSED/p21-IsoSeq.flnc.bam $PROCESSED/p21-IsoSeq.clustered.bam --num-threads $THREADS

# 使用conda install -c bioconda pbmm2下载了最新的pbmm2 1.17.0
# 3. pbmm2：比对到基因组
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start pbmm2 align"
pbmm2 align $GENOME_FA $PROCESSED/p21-IsoSeq.clustered.bam $PROCESSED/p21-IsoSeq.mapped.bam --preset ISOSEQ --sort --log-level INFO --num-threads $THREADS

echo "Iso-Seq3 bulk预处理流程（1-3）完成！" 