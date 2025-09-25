#!/usr/bin/sh

################################################
#File Name: Run.mapping.20250227.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 27 Feb 2025 03:50:37 PM CST
################################################

set -eo pipefail

# 输入基本质控后的单个fastq文件，根据文件名输出结构化的输出结果
## 第一层为研究名，第二层为样本名
## 输出output/filtered_fq以及output/alignment

# trimmed fastq file $1
trimmed_fastq=$1
output_path=$2
[ -d $output_path ] || mkdir -p $output_path
### 根据输入的文件名替换 ###
sample=$(basename -s .fq.gz $trimmed_fastq|sed 's/_trimmed//')

### 参考基因组注释文件，根据物种替换 ###
#### 修改 ####
bowtie2_rRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/ncRNA/human/merged/hg38.rRNA
bowtie2_tRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/ncRNA/human/merged/hg38.tRNA
bowtie2_snoRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/ncRNA/human/merged/hg38.snoRNA
STARindex=/home/user/data3/lit/resource/genome/human/hg38/hg38_STARindex_v2.5.2b/

# 0.创立文件夹
cd $output_path
### 根据输入的文件名替换 ###
subdir=$(echo $trimmed_fastq | cut -d "/" -f 12)
[ -d $subdir/$sample ] || mkdir -p $subdir/$sample
cd $subdir/$sample
[ -d output ] || mkdir output
[ -d log ] || mkdir log

# 激活虚拟环境
source activate biotools

# 2.去除contaminant，sequentially
## 2.1 rRNA
echo -e "Remove rRNA comtaminant at $(date '+%Y-%m-%d %H:%M:%S')"
[ -d output/filtered_fq ] || mkdir output/filtered_fq
bowtie2 -q \
--phred33 \
-N 1 \
-x $bowtie2_rRNA_index  \
-U $trimmed_fastq \
-S /dev/null \
-t \
-p 20 \
--no-unal \
--un-gz output/filtered_fq/$sample.trimmed.rRNA.unaligned.fq.gz &> log/bowtie.rRNA.log
## 2.2 tRNA
echo -e "Remove tRNA comtaminant at $(date '+%Y-%m-%d %H:%M:%S')"
bowtie2 -q \
--phred33 \
-N 1 \
-x $bowtie2_tRNA_index  \
-U output/filtered_fq/$sample.trimmed.rRNA.unaligned.fq.gz \
-S /dev/null \
-t \
-p 20 \
--no-unal \
--un-gz output/filtered_fq/$sample.trimmed.rRNA.tRNA.unaligned.fq.gz &> log/bowtie.rRNA.tRNA.log
## 2.3 snoRNA
echo -e "Remove snoRNA comtaminant at $(date '+%Y-%m-%d %H:%M:%S')"
bowtie2 -q \
--phred33 \
-N 1 \
-x $bowtie2_snoRNA_index  \
-U output/filtered_fq/$sample.trimmed.rRNA.tRNA.unaligned.fq.gz \
-S /dev/null \
-t \
-p 20 \
--no-unal \
--un-gz output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz &> log/bowtie.rRNA.tRNA.snoRNA.log

# 清理fastq文件
echo -e "Statistics of fastq files at $(date '+%Y-%m-%d %H:%M:%S')"
seqkit stats -j 20 $trimmed_fastq output/filtered_fq/$sample.trimmed.rRNA.unaligned.fq.gz \
    output/filtered_fq/$sample.trimmed.rRNA.tRNA.unaligned.fq.gz \
    output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz &> output/fq.stat.txt && \
# 检查 seqkit stats 命令的退出状态
if [ $? -eq 0 ]; then
    rm -rf output/filtered_fq/$sample.trimmed.rRNA.unaligned.fq.gz output/filtered_fq/$sample.trimmed.rRNA.tRNA.unaligned.fq.gz
else
    echo "seqkit stats 命令执行失败，不进行文件删除操作。"
fi

# 3. 比对到参考基因组上
echo -e "Mapping at $(date '+%Y-%m-%d %H:%M:%S')"
cleaned_fastq=output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz
[ -d output/alignment ] || mkdir output/alignment
STAR \
--readFilesIn $cleaned_fastq \
--readFilesCommand zcat \
--seedSearchStartLmax 15 \
--runThreadN 40 \
--outFilterMismatchNmax 2 \
--genomeDir $STARindex \
--outFileNamePrefix output/alignment/${sample}_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outFilterMultimapNmax 1 \
--outFilterMatchNmin 16 \
--outSAMattributes All \
--alignEndsType EndToEnd &> log/STAR.log

samtools view output/alignment/${sample}_Aligned.sortedByCoord.out.bam -h -@ 20 > output/alignment/$sample.sam