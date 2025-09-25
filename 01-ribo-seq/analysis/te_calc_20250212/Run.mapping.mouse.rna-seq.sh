#!/usr/bin/sh

################################################
#File Name: Run.mapping.mouse.rna-seq.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 13 Feb 2025 11:28:12 AM CST
################################################

# 这个脚本是对配套的RNA-seq数据进行分析，得到每个样本的bam文件

set -eo pipefail

output_path=output/mapping_rna_seq
STARindex=/home/user/data3/lit/resource/genome/mouse/mm39/index/mm39_STARindex
mouse_brain_data_path=/home/user/data3/licq/peptidomics/PublicData/mouse_brain/
mkdir -p output/mapping_rna_seq

ls $mouse_brain_data_path/2020_NAR_WangH_E15.5_P42/RNA/*fastq.gz \
	$mouse_brain_data_path/2015_Science_ChoJ_hippocampus/RNA/*fastq.gz \
	$mouse_brain_data_path/2019_multispecies_mouse_brain/RNA/*fastq.gz > $output_path/rna_fastq.lst

cd $output_path
source activate biotools
# 输出所有的SAM attributes
for fastq in $(cat rna_fastq.lst);do
	sample=$(basename -s .fastq.gz $fastq)
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