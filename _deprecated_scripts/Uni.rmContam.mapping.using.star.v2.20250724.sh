#!/usr/bin/sh

################################################
#File Name: Ribo.Run.Uni.test.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 03 Oct 2024 08:33:39 PM CST
################################################

set -eo pipefail

trimmed_fastq=$1
output_path=$2
mkdir -p $output_path
sample=$(basename "$trimmed_fastq" | sed -E 's/_trimmed\.(fastq|fq)\.gz$//')
### 参考基因组注释文件，根据物种替换 ###
#### 修改 ####
star_rRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/ncRNA/STARindex/rRNA/rRNA_STARindex
star_tRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/ncRNA/STARindex/tRNA/tRNA_STARindex
star_snoRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/ncRNA/STARindex/snoRNA/snoRNA_STARindex
STARindex=/home/user/data3/lit/resource/genome/human/hg38/hg38_STARindex_v2.5.2b/
# 0.创立文件夹
cd $output_path
[ -d $sample ] || mkdir -p $sample
cd $sample
[ -d output ] || mkdir output
[ -d log ] || mkdir log

# 激活虚拟环境
source activate biotools
echo -e "***Processing $sample at $(date '+%Y-%m-%d %H:%M:%S')"

# 2.去除contaminant，sequentially
## 2.1 rRNA
echo -e "Remove rRNA comtaminant at $(date '+%Y-%m-%d %H:%M:%S')"
[ -d output/filtered_fq ] || mkdir output/filtered_fq

STAR --runThreadN 50 \
--limitBAMsortRAM 250000000000 \
--genomeDir $star_rRNA_index \
--readFilesIn $trimmed_fastq \
--readFilesCommand zcat \
--outFileNamePrefix output/filtered_fq/$sample.trimmed.rRNA. \
--outSAMtype None \
--outReadsUnmapped Fastx \
--alignEndsType EndToEnd \
--alignIntronMax 1 \
--outFilterMismatchNmax 1 &> log/star.rRNA.log
## 2.2 tRNA
echo -e "Remove tRNA comtaminant at $(date '+%Y-%m-%d %H:%M:%S')"
STAR --runThreadN 50 \
--limitBAMsortRAM 250000000000 \
--genomeDir $star_tRNA_index \
--readFilesIn output/filtered_fq/$sample.trimmed.rRNA.Unmapped.out.mate1 \
--outFileNamePrefix output/filtered_fq/$sample.trimmed.rRNA.tRNA. \
--outSAMtype None \
--outReadsUnmapped Fastx \
--alignEndsType EndToEnd \
--alignIntronMax 1 \
--outFilterMismatchNmax 1 &> log/star.rRNA.tRNA.log
## 2.3 snoRNA
echo -e "Remove snoRNA comtaminant at $(date '+%Y-%m-%d %H:%M:%S')"
STAR --runThreadN 50 \
--limitBAMsortRAM 250000000000 \
--genomeDir $star_snoRNA_index \
--readFilesIn output/filtered_fq/$sample.trimmed.rRNA.tRNA.Unmapped.out.mate1 \
--outFileNamePrefix output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA. \
--outSAMtype None \
--outReadsUnmapped Fastx \
--alignEndsType EndToEnd \
--alignIntronMax 1 \
--outFilterMismatchNmax 1 &> log/star.rRNA.tRNA.snoRNA.log
# 清理fastq文件
echo -e "Statistics of fastq files at $(date '+%Y-%m-%d %H:%M:%S')"
seqkit stats -j 30 $trimmed_fastq output/filtered_fq/$sample.trimmed.rRNA.Unmapped.out.mate1 \
    output/filtered_fq/$sample.trimmed.rRNA.tRNA.Unmapped.out.mate1 \
    output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.Unmapped.out.mate1 &> output/fq.stat.txt
KEEP_INTERMEDIATE=false
if [ "$KEEP_INTERMEDIATE" = false ]; then
    rm -rf output/filtered_fq/$sample.trimmed.rRNA.Unmapped.out.mate1 output/filtered_fq/$sample.trimmed.rRNA.tRNA.Unmapped.out.mate1
fi
cleaned_fastq=output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz
[ -f $cleaned_fastq ] || pigz -p 30 -c output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.Unmapped.out.mate1 > $cleaned_fastq && rm output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.Unmapped.out.mate1

# 3. 比对到参考基因组上
echo -e "Mapping at $(date '+%Y-%m-%d %H:%M:%S')"
[ -d output/alignment ] || mkdir output/alignment
STAR --runThreadN 50 \
--limitBAMsortRAM 250000000000 \
--readFilesIn $cleaned_fastq \
--readFilesCommand zcat \
--seedSearchStartLmax 15 \
--outFilterMismatchNmax 2 \
--genomeDir $STARindex \
--outFileNamePrefix output/alignment/${sample}_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outFilterMultimapNmax 1 \
--outFilterMatchNmin 16 \
--outSAMattributes All \
--alignEndsType EndToEnd &> log/STAR.log
