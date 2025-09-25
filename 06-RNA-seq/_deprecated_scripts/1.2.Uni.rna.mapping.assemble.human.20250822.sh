#!/usr/bin/sh

################################################
#File Name: 1.1.Uni.rna.mapping.assemble.human.20250227.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 27 Feb 2025 05:05:38 PM CST
################################################

set -eo pipefail

# 输入基本质控后的包含所有文件名的文件
## 对于双端数据自动识别

# trimmed fastq file $1
trimmed_fastq_lst=$1
# human_brain_rna_seq_alignment_assemble
output_path=$2
[ -d $output_path ] || mkdir -p $output_path

### 参考基因组注释文件，根据物种替换 ###
#### 修改 ####
STARindex=/home/user/data3/lit/resource/genome/human/hg38/hg38_custom_STARindex_v2.5.2b/
ref_gtf=/home/user/data3/lit/project/sORFs/09-CustomDb/Test_20250801/processed/mkAnno_for_moPepGen/custom.gtf

# 0.创立文件夹
cd $output_path
[ -d log ] || mkdir log

source activate biotools
# 1.比对
# 输出所有的SAM attributes
for fastq in $(cat $trimmed_fastq_lst);do
	#### 根据文件特点修改，一般为.fastq.gz或者.fq.gz结尾，或者样本名后面带有_trimmed【可选】就可以 ####
	# sample=$(basename $fastq|sed 's/.fastq.gz//;s/.fq.gz//'|sed 's/_trimmed//')
	sample=$(basename $fastq|sed 's/.clean.fastq.gz//')
	echo "***Processing $sample"
	if [[ ${sample#*.} == R1 ]]; then
		echo "Sample $sample ends with R1, assuming paired-end sequencing"
		dir=$(dirname "$fastq");new_base=$(basename "$fastq" | sed 's/R1/R2/');pair_file="$dir/$new_base"
		if [[ -f "$pair_file" ]]; then
            echo "Corresponding R2 file found: $pair_file"
            echo "Running paired-end code for "${sample%.*}""
			STAR \
			--readFilesIn $fastq $pair_file \
			--readFilesCommand zcat \
			--runThreadN 20 \
			--genomeDir $STARindex \
			--outFileNamePrefix ${sample}_ \
			--outSAMtype BAM SortedByCoordinate \
			--outFilterMultimapNmax 1  \
			--outSAMattributes All \
			--limitBAMsortRAM 202080440150
        else
            echo "Warning: Corresponding R2 file not found: $pair_file"
        fi
	elif [[ ${sample#*.} == R2 ]]; then
        echo "Sample $sample ends with R2, skipping processing"
    else
		STAR \
		--readFilesIn $fastq \
		--readFilesCommand zcat \
		--runThreadN 20 \
		--genomeDir $STARindex \
		--outFileNamePrefix ${sample}_ \
		--outSAMtype BAM SortedByCoordinate \
		--outFilterMultimapNmax 1  \
		--outSAMattributes All \
		--limitBAMsortRAM 202080440150
	fi
done &> log/STAR.log



