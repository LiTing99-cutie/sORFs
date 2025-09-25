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
# 20241211
#### 修改STAR的参数 ####

# raw fastq file $1
raw_fastq_dir=$1
# 提取R1和R2文件路径
raw_fastq_R1=$(ls ${raw_fastq_dir}/*_R1*.fq.gz)
raw_fastq_R2=$(ls ${raw_fastq_dir}/*_R2*.fq.gz)
sample=$(basename -s _R1.fq.gz $raw_fastq_R1)
# output dir $2
#### 修改 ####
output_path=$2
mkdir -p $output_path

### 参考基因组注释文件，根据物种替换 ###
#### 修改 ####
bowtie2_rRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/ncRNA/human/merged/hg38.rRNA
bowtie2_tRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/ncRNA/human/merged/hg38.tRNA
bowtie2_snoRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/ncRNA/human/merged/hg38.snoRNA
STARindex=/home/user/data/lit/database/public/genome/hg38/hg38_STARindex_v2.5.2b
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/RiboCode_annot/hg38

# 0.创立文件夹
cd $output_path
subdir=$(echo $raw_fastq_R1 | cut -d "/" -f 9)
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
fastqc -o output/fastqc -t 10 $raw_fastq_R1 $raw_fastq_R2 &> log/raw_fastqc.log
# 使用trim_galore去接头，过滤低质量reads
echo -e "Trim adapter at $(date '+%Y-%m-%d %H:%M:%S')"
[ -d output/trimmed_fastq ] || mkdir output/trimmed_fastq
cutadapt -j 20 --trim-n --match-read-wildcards -u 3 -n 4 -a AGATCGGAAGAGCACACGTCTG -a AAAAAAAA -a GAACTCCAGTCAC \
-A AGATCGGAAGAGCACACGTCTG -A AAAAAAAA -A GAACTCCAGTCAC \
-e 0.2 --nextseq-trim 20 -m 15 -o output/trimmed_fastq/$sample.trimmed.fastq.1.gz -p output/trimmed_fastq/$sample.trimmed.fastq.2.gz $raw_fastq_R1 $raw_fastq_R2
## 再次fastqc
echo -e "Trimmed reads fastqc at $(date '+%Y-%m-%d %H:%M:%S')"
trimmed_fastq_R1=$(ls output/trimmed_fastq/*trimmed.fastq.1.gz)
trimmed_fastq_R2=$(ls output/trimmed_fastq/*trimmed.fastq.2.gz)
fastqc -o output/fastqc -t 10 $trimmed_fastq_R1 $trimmed_fastq_R2 &> log/trimmed_fastqc.log

# 2.去除contaminant，sequentially
## 2.1 rRNA
echo -e "Remove rRNA comtaminant at $(date '+%Y-%m-%d %H:%M:%S')"
[ -d output/filtered_fq ] || mkdir output/filtered_fq
bowtie2 -q \
--phred33 \
-N 1 \
-x $bowtie2_rRNA_index  \
-1 $trimmed_fastq_R1 -2 $trimmed_fastq_R2 \
-S /dev/null \
-t \
-p 20 \
--no-unal \
--un-conc-gz output/filtered_fq/$sample.trimmed.rRNA.unaligned.fq.gz &> log/bowtie.rRNA.log

## 2.2 tRNA
echo -e "Remove tRNA comtaminant at $(date '+%Y-%m-%d %H:%M:%S')"
bowtie2 -q \
--phred33 \
-N 1 \
-x $bowtie2_tRNA_index  \
-1 output/filtered_fq/$sample.trimmed.rRNA.unaligned.fq.1.gz -2 output/filtered_fq/$sample.trimmed.rRNA.unaligned.fq.2.gz \
-S /dev/null \
-t \
-p 20 \
--no-unal \
--un-conc-gz output/filtered_fq/$sample.trimmed.rRNA.tRNA.unaligned.fq.gz &> log/bowtie.rRNA.tRNA.log
## 2.3 snoRNA
echo -e "Remove snoRNA comtaminant at $(date '+%Y-%m-%d %H:%M:%S')"
bowtie2 -q \
--phred33 \
-N 1 \
-x $bowtie2_snoRNA_index  \
-1 output/filtered_fq/$sample.trimmed.rRNA.tRNA.unaligned.fq.1.gz -2 output/filtered_fq/$sample.trimmed.rRNA.tRNA.unaligned.fq.2.gz \
-S /dev/null \
-t \
-p 20 \
--no-unal \
--un-conc-gz output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz &> log/bowtie.rRNA.tRNA.snoRNA.log

# 清理fastq文件
echo -e "Statistics of fastq files at $(date '+%Y-%m-%d %H:%M:%S')"
seqkit stats -j 20 $raw_fastq_R1 $raw_fastq_R2 $trimmed_fastq_R1 $trimmed_fastq_R2 \
output/filtered_fq/$sample.trimmed.rRNA.unaligned*gz \
output/filtered_fq/$sample.trimmed.rRNA.unaligned*gz \
output/filtered_fq/$sample.trimmed.rRNA.tRNA.unaligned*gz \
output/filtered_fq/$sample.trimmed.rRNA.tRNA.unaligned*gz \
output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.unaligned*gz \
output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.unaligned*gz &> output/fq.stat.txt
rm -rf $trimmed_fastq_R1 $trimmed_fastq_R2 output/filtered_fq/$sample.trimmed.rRNA.unaligned.*.fq.gz output/filtered_fq/$sample.trimmed.rRNA.tRNA.unaligned.*.fq.gz

# 3. 比对到参考基因组上
echo -e "Mapping at $(date '+%Y-%m-%d %H:%M:%S')"
cleaned_fastq_R1=output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.unaligned.fq.1.gz
cleaned_fastq_R2=output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.unaligned.fq.2.gz
[ -d output/alignment ] || mkdir output/alignment
STAR \
--readFilesIn $cleaned_fastq_R1 $cleaned_fastq_R2 \
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

# 4.使用RiboCode注释ORFs
source activate ribocode
echo -e "Annotate ORFs at $(date '+%Y-%m-%d %H:%M:%S')"
[ -d output/Ribo-ORFs/RiboCode ] || mkdir -p output/Ribo-ORFs/RiboCode
bam=$PWD/output/alignment/${sample}_Aligned.toTranscriptome.out.bam
cd output/Ribo-ORFs/RiboCode && \
metaplots -a $RiboCode_annot -r $bam && \
RiboCode -a $RiboCode_annot -c metaplots_pre_config.txt -l no -g -o ./$sample --output-gtf --output-bed --min-AA-length 6 &> ../../../log/anno_orfs_ribocode.log
