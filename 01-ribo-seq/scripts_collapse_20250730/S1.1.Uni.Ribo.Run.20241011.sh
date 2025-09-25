#!/usr/bin/sh

################################################
#File Name: Ribo.Run.Uni.test.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 03 Oct 2024 08:33:39 PM CST
################################################

set -eo pipefail

# 20241011
#### 修改输出文件夹 ####
#### 更改了rRNA等污染RNA的index ####

# raw fastq file $1
raw_fastq=$1
sample=$(basename -s .fastq.gz $raw_fastq)
# adapter $2
adapter=$2
# output dir $3
#### 修改 ####
output_path=./mouse_brain_output_20241011

### 参考基因组注释文件，根据物种替换 ###
#### 修改 ####
bowtie2_rRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/ncRNA/mouse/merged/mouse.rRNA
bowtie2_tRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/ncRNA/mouse/merged/mouse.tRNA
bowtie2_snoRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/ncRNA/mouse/merged/mouse.snoRNA
STARindex=/home/user/data3/lit/resource/genome/mouse/mm39/index/mm39_STARindex
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/RiboCode_annot/mm39

# 0.创立文件夹
cd $output_path
subdir=$(echo $raw_fastq | cut -d "/" -f 9)
[ -d $subdir/$sample ] || mkdir -p $subdir/$sample
cd $subdir/$sample
[ -d output ] || mkdir output
[ -d log ] || mkdir log

# 激活虚拟环境
source activate biotools
echo -e "***Processing $subdir $sample at $(date '+%Y-%m-%d %H:%M:%S')"
# 1.原始测序数据质量控制
## fastqc
echo -e "Raw reads fastqc at $(date '+%Y-%m-%d %H:%M:%S')"
[ -d output/fastqc ] || mkdir output/fastqc
fastqc -o output/fastqc -t 10 $raw_fastq &> log/raw_fastqc.log
# 使用trim_galore去接头，过滤低质量reads
echo -e "Trim adapter at $(date '+%Y-%m-%d %H:%M:%S')"
[ -d output/trimmed_fastq ] || mkdir output/trimmed_fastq
if [ -z "$adapter" ]; then
time trim_galore -j 8 -q 20 --length 20 $raw_fastq --gzip -o output/trimmed_fastq &> log/trim_galore.log
else
time trim_galore --adapter $adapter -j 8 -q 20 --length 20 $raw_fastq --gzip -o output/trimmed_fastq &> log/trim_galore.log
fi
## 再次fastqc
echo -e "Trimmed reads fastqc at $(date '+%Y-%m-%d %H:%M:%S')"
trimmed_fastq=$(ls output/trimmed_fastq/*gz)
fastqc -o output/fastqc -t 10 $trimmed_fastq &> log/trimmed_fastqc.log

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
seqkit stats -j 20 $raw_fastq $trimmed_fastq output/filtered_fq/$sample.trimmed.rRNA.unaligned.fq.gz \
output/filtered_fq/$sample.trimmed.rRNA.tRNA.unaligned.fq.gz \
output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz &> output/fq.stat.txt
rm -rf $trimmed_fastq output/filtered_fq/$sample.trimmed.rRNA.unaligned.fq.gz output/filtered_fq/$sample.trimmed.rRNA.tRNA.unaligned.fq.gz

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
--alignEndsType EndToEnd &> log/STAR.log

# 4.使用RiboCode注释ORFs
source activate ribocode
echo -e "Annotate ORFs at $(date '+%Y-%m-%d %H:%M:%S')"
[ -d output/Ribo-ORFs/RiboCode ] || mkdir -p output/Ribo-ORFs/RiboCode
bam=$PWD/output/alignment/${sample}_Aligned.toTranscriptome.out.bam
cd output/Ribo-ORFs/RiboCode && \
metaplots -a $RiboCode_annot -r $bam && \
RiboCode -a $RiboCode_annot -c metaplots_pre_config.txt -l no -g -o ./$sample --output-gtf --output-bed --min-AA-length 6 &> ../../../log/anno_orfs_ribocode.log