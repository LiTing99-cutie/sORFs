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