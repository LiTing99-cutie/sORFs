#!/usr/bin/sh

################################################
#File Name: Run.assemble.trans.20250114.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 14 Jan 2025 03:57:36 PM CST
################################################

set -eo pipefail

mkdir -p output/assembled_trans
output_path=output/assembled_trans
mouse_brain_data_path=/home/user/data3/licq/peptidomics/PublicData/mouse_brain/
ls $mouse_brain_data_path/2020_NAR_WangH_E15.5_P42/RNA/*fastq.gz \
	$mouse_brain_data_path/2015_Science_ChoJ_hippocampus/RNA/*fastq.gz \
	$mouse_brain_data_path/2019_multispecies_mouse_brain/RNA/*fastq.gz > output/assembled_trans/rna_fastq.lst
less output/assembled_trans/rna_fastq.lst |xargs zcat |pigz -p 10 > $output_path/merged.fastq.gz

source activate biotools
cd $output_path
mkdir -p alignment/mouse_brain
mkdir log
STARindex=/home/user/data3/lit/resource/genome/mouse/mm39/index/mm39_STARindex
STAR \
--readFilesIn merged.fastq.gz \
--readFilesCommand zcat \
--runThreadN 40 \
--genomeDir $STARindex \
--outFileNamePrefix alignment/mouse_brain \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outFilterMultimapNmax 1  &> log/STAR.log

ref_gtf=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.gtf
mkdir -p stringtie_output
stringtie ./alignment/mouse_brainAligned.sortedByCoord.out.bam \
          -o stringtie_output/mouse_brain.gtf \
          -p 10 \
          -G $ref_gtf

# 全部的转录本条数
less output/assembled_trans/stringtie_output/mouse_brain.gtf |awk '$3=="transcript"'|wc -l
# 新组装的转录本
less output/assembled_trans/stringtie_output/mouse_brain.gtf |grep -v reference_id|awk '$3=="transcript"'|wc -l

# 新输出一个文件，辅助确认是否上面得到的是新组装的转录本
stringtie ./alignment/mouse_brainAligned.sortedByCoord.out.bam \
          -o stringtie_output/mouse_brain.gtf.1 \
          -p 10 \
          -G $ref_gtf \
        -A stringtie_output/gene_abund.tab
# 后续报错，去掉链信息为.的行，大概1%的行链是.
less stringtie_output/mouse_brain.gtf |grep -v reference_id|awk '$7!="."' > stringtie_output/mouse_brain_new_trans.gtf
# run /home/user/data3/lit/project/sORFs/01-ribo-seq/Run.20250126.reformat.stringtie.gtf.R，得到stringtie_output/mouse_brain_new_trans_reformat_ordered.gtf
cat /home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.gtf stringtie_output/mouse_brain_new_trans_reformat_ordered.gtf > \
 stringtie_output/gencode.vM29.add_assemble.gtf
gtfToGenePred -geneNameAsName2 -genePredExt stringtie_output/gencode.vM29.add_assemble.gtf stringtie_output/gencode.vM29.add_assemble.genePred.txt
