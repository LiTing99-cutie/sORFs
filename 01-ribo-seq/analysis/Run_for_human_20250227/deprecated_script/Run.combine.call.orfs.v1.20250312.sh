#!/usr/bin/sh

################################################
#File Name: Run.combine.call.orfs.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 12 Mar 2025 03:31:55 PM CST
################################################

set -eo pipefail

proj_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227
output_path=$proj_path/human_brain_ribo_merge_call_orfs_20250312
mkdir -p $output_path && cd $output_path
# find $proj_path/human_brain_output_20250227 -name *_Aligned.toTranscriptome.out.bam|grep -f $proj_path/human_brain_output_20250227/high_quality_run.txt > bam_lst.txt
# find $proj_path/human_brain_output_20250227 -name *_Aligned.sortedByCoord.out.bam|grep -f $proj_path/human_brain_output_20250227/high_quality_run.txt > bam_1_lst.txt
# find $proj_path/human_brain_output_20250227 -name *sam|grep -f $proj_path/human_brain_output_20250227/high_quality_run.txt|grep -v corrected > sam_lst.txt
# samtools merge -@ 40 -f Merged_Aligned.toTranscriptome.out.bam -b bam_lst.txt
# samtools merge -@ 40 -f Merged_Aligned.sortedByCoord.out.bam -b bam_1_lst.txt
# samtools merge -@ 40 -f Merged.sam -b sam_lst.txt
# samtools index -@ 40 Merged_Aligned.sortedByCoord.out.bam

# 基因组及其注释
fa=/home/user/data3/lit/resource/genome/human/hg38/hg38.fa
gtf=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
# 注意输入的是只有10列的gpe文件
genePred=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gpe.ten_col.txt
STARindex=/home/user/data3/lit/resource/genome/human/hg38/hg38_STARindex_v2.5.2b/
RiboCode_annot=$proj_path/annotation/RiboCode_annot
candidateORF=$proj_path/annotation/RibORF_annot/candidateORF.genepred.txt
pvalue_cutoff=1
PRICE_anno_name=hg38_gencv41
bam=$output_path/Merged_Aligned.toTranscriptome.out.bam
bam_1=$output_path/Merged_Aligned.sortedByCoord.out.bam
bam_1_index=$output_path/Merged_Aligned.sortedByCoord.out.bam.bai
sam=$output_path/Merged.sam
sample=Human_brain

##### RiboCode #####
echo -e "RiboCode start at $(date '+%Y-%m-%d %H:%M:%S')"
source activate ribocode
mkdir -p $output_path/RiboCode && cd $output_path/RiboCode
metaplots -a $RiboCode_annot -r $bam && \
# 增加9种替代的起始密码子
RiboCode --alt_start_codons CTG,ACG,GTG,TTG,ATA,ATC,ATT,AAG,AGG --pval-cutoff $pvalue_cutoff -a $RiboCode_annot -c metaplots_pre_config.txt -l no -g -o ./$sample --output-gtf --output-bed --min-AA-length 6
echo -e "RiboCode done at $(date '+%Y-%m-%d %H:%M:%S')"

##### PRICE #####
echo -e "PRICE start at $(date '+%Y-%m-%d %H:%M:%S')"
mkdir -p $output_path/PRICE && cd $output_path/PRICE
# 先清除文件夹下的文件，再运行
rm -rf ./*
[ -f $bam_1_index ] || samtools index -@ 20 $bam_1
gedi -e Price -reads $bam_1 -genomic $PRICE_anno_name -prefix $sample -plot
gedi -e ViewCIT -m bed -name 'd.getType()'  *.orfs.cit > $sample.orfs.bed
echo -e "PRICE done at $(date '+%Y-%m-%d %H:%M:%S')"

##### RibORF #####
echo -e "RibORF start at $(date '+%Y-%m-%d %H:%M:%S')"
mkdir -p $output_path/RibORF && cd $output_path/RibORF
RibORF_path=/home/user/data2/lit/software/RibORF/RibORF.2.0
# 使用ribocode鉴定得到的read length以及offset
config=../RiboCode/metaplots_pre_config.txt
grep "^#" $config|awk -v OFS='\t' '{print $2,$4+3}'|tail -n +3|head -n -1 > offset.correction.parameters.txt
# 根据offset文件校正
time perl $RibORF_path/offsetCorrect.pl -r $sam -p offset.correction.parameters.txt -o corrected.$sample.sam
# 可视化校正后的结果
time perl $RibORF_path/readDist.pl -f corrected.$sample.sam -g $genePred -o ./ -d 1
# call ORFs
time perl $RibORF_path/ribORF.pl -f corrected.$sample.sam -c $candidateORF -o ./
echo -e "RibORF done at $(date '+%Y-%m-%d %H:%M:%S')"


