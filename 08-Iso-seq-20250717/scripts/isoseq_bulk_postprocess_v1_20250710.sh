#!/bin/bash
set -e

# Iso-Seq3 bulk postprocess v1_20250710
# 步骤4-7：collapse, pigeon classify, pigeon filter, pigeon report
# 输出目录：../processed

# 激活biotools环境
source activate biotools

# 路径定义
GENOME_FA=/home/user/data3/lit/project/sORFs/07-Genome/results/custom_fa/custom_ref.fa
GTF_ANN=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
PROCESSED=../processed
THREADS=30

mkdir -p $PROCESSED

# 4. isoseq collapse：冗余转录本合并
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start isoseq3 collapse"
isoseq3 collapse --do-not-collapse-extra-5exons $PROCESSED/p21-IsoSeq.mapped.bam $PROCESSED/p21-IsoSeq.flnc.bam $PROCESSED/collapsed.gff

# 5. pigeon classify：转录本注释分类
source activate pbpigeon
pigeon prepare $GTF_ANN $GENOME_FA
pigeon prepare $PROCESSED/collapsed.gff

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start pigeon classify"
mkdir -p $PROCESSED/classify
classify_dir=$PROCESSED/classify
pigeon classify $PROCESSED/collapsed.sorted.gff ${GTF_ANN%.gtf}.sorted.gtf $GENOME_FA -d $classify_dir --fl $PROCESSED/collapsed.flnc_count.txt --num-threads $THREADS

# 6. pigeon filter：过滤潜在假阳性转录本
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start pigeon filter"
pigeon filter $classify_dir/collapsed_classification.txt --isoforms $PROCESSED/collapsed.sorted.gff --num-threads $THREADS

# 7. pigeon report：报告
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start pigeon report"
pigeon report --exclude-singletons $classify_dir/collapsed_classification.filtered_lite_classification.txt $PROCESSED/saturation.txt

echo "Iso-Seq3 bulk后处理流程（4-7）完成！" 