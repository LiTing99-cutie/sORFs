#!/usr/bin/sh

################################################
#File Name: Run.mapping.mouse.rna-seq.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 13 Feb 2025 11:28:12 AM CST
################################################

# 这个脚本是对配套的RNA-seq数据进行分析，得到每个样本的bam文件

set -eo pipefail

# output_path=output/mapping_rna_seq
output_path=$1
# STARindex=/home/user/data3/lit/resource/genome/mouse/mm39/index/mm39_STARindex
STARindex=$2
fastq_lst_path=$3
mkdir -p $output_path

cd $output_path
source activate biotools
# 输出所有的SAM attributes
for fastq in $(cat $fastq_lst_path);do
	sample=$(basename $fastq | cut -d "." -f 1)
	echo "***Processing $sample"
	STAR \
	--readFilesIn $fastq \
	--readFilesCommand zcat \
	--runThreadN 20 \
	--genomeDir $STARindex \
	--outFileNamePrefix $sample \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM GeneCounts \
	--outFilterMultimapNmax 1  \
	--outSAMattributes All \
	--limitBAMsortRAM 32000000000
done