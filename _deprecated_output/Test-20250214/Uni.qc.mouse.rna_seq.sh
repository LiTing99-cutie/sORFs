#!/usr/bin/sh

################################################
#File Name: Run.qc.mouse.rna_seq.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 13 Feb 2025 04:12:04 PM CST
################################################

set -eo pipefail

# raw fastq file $1
raw_fastq=$1
sample=$(basename -s .fastq.gz $raw_fastq)
# adapter $2
adapter=$2
# output dir $3
output_path=$3
mkdir -p $output_path
fastqc=$4
rm_adpt=$5
trimmed_fastqc=$6

# 0.创立文件夹
cd $output_path
[ -d $sample ] || mkdir -p $sample
cd $sample
[ -d output ] || mkdir output
[ -d log ] || mkdir log
[ -d output/fastqc ] || mkdir output/fastqc
[ -d output/trimmed_fastq ] || mkdir output/trimmed_fastq
# 激活虚拟环境
source activate biotools
echo -e "***Processing $subdir $sample at $(date '+%Y-%m-%d %H:%M:%S')"
# 1.原始测序数据质量控制
## fastqc
if [ $fastqc = "yes" ]; then
	echo -e "Raw reads fastqc at $(date '+%Y-%m-%d %H:%M:%S')"
	[ -d output/fastqc ] || mkdir output/fastqc
	fastqc -o output/fastqc -t 10 $raw_fastq &> log/raw_fastqc.log
fi
# 使用trim_galore去接头，过滤低质量reads
if [ $rm_adpt = "yes" ]; then
	echo -e "Trim adapter at $(date '+%Y-%m-%d %H:%M:%S')"
	[ -d output/trimmed_fastq ] || mkdir output/trimmed_fastq
	if [ $adapter = "none" ]; then
	time trim_galore -j 8 -q 20 --length 20 $raw_fastq --gzip -o output/trimmed_fastq &> log/trim_galore.log
	else
	time trim_galore --adapter "file:$adapter" -j 8 -q 20 --length 20 $raw_fastq --gzip -o output/trimmed_fastq &> log/trim_galore.log
	fi
fi
## 再次fastqc
if [ $trimmed_fastqc = "yes" ]; then
	echo -e "Trimmed reads fastqc at $(date '+%Y-%m-%d %H:%M:%S')"
	trimmed_fastq=$(ls output/trimmed_fastq/*gz)
	fastqc -o output/fastqc -t 10 $trimmed_fastq &> log/trimmed_fastqc.log
fi