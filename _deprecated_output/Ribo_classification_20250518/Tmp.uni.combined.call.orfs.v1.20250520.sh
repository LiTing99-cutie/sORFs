#!/usr/bin/sh

################################################
#File Name: 2.2b.Run.combine.call.orfs.20250331.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 31 Mar 2025 11:48:59 AM CST
################################################

set -eo pipefail

# 基于三种软件来鉴定ORF，要求reads的周期性大于等于50%，且显著性小于0.01
# output_path=$proj_path/human_brain_ribo_merge_call_orfs_20250338
# mkdir -p $output_path && cd $output_path
# bam=$output_path/merged_filtered_reads.toTranscriptome.bam
# bam_1=$output_path/merged_filtered_reads.sortedByCoord.bam
# bam_1_index=$output_path/merged_filtered_reads.sortedByCoord.bam.bai
# sam=$output_path/merged_filtered_reads.sortedByCoord.sam

# 使用方法：
# bash Uni.combine.call.orfs.v1.20250520.sh $PWD/output/orfs $bam_dir/merged.trans.bam $bam_dir/merged.bam $bam_dir/merged.sam 

output_path=$1
mkdir -p $output_path && cd $output_path
# trans bam
bam=$2
bam_1=$3
bam_1_index=${bam_1}.bai
sam=$4

proj_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227

# 基因组及其注释
fa=/home/user/data3/lit/resource/genome/human/hg38/hg38.fa
gtf=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
# 注意输入的是只有10列的gpe文件
genePred=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.10.gpe
STARindex=/home/user/data3/lit/resource/genome/human/hg38/hg38_STARindex_v2.5.2b/
RiboCode_annot=$proj_path/annotation/RiboCode_annot
candidateORF=$proj_path/annotation/RibORF_annot/candidateORF.genepred.txt
pvalue_cutoff=1
PRICE_anno_name=hg38_gencv41
sample=Human_brain

##### RiboCode #####
# echo -e "RiboCode start at $(date '+%Y-%m-%d %H:%M:%S')"
# source activate ribocode
# mkdir -p $output_path/RiboCode && cd $output_path/RiboCode
# # 20250520 修改参数
# metaplots -a $RiboCode_annot -r $bam -f0_percent 0.5 -pv1 0.01 -pv2 0.01 && \
# # 增加9种替代的起始密码子
# RiboCode --alt_start_codons CTG,ACG,GTG,TTG,ATA,ATC,ATT,AAG,AGG --pval-cutoff $pvalue_cutoff -a $RiboCode_annot -c metaplots_pre_config.txt -l no -g -o ./$sample --output-gtf --output-bed --min-AA-length 6
# echo -e "RiboCode done at $(date '+%Y-%m-%d %H:%M:%S')"

##### PRICE #####
echo -e "PRICE start at $(date '+%Y-%m-%d %H:%M:%S')"
mkdir -p $output_path/PRICE && cd $output_path/PRICE
# 先清除文件夹下的文件，再运行
rm -rf ./*
[ -f $bam_1_index ] || samtools index -@ 20 $bam_1
# bam9.4G，java的内存爆了，把这个调整大一点
gedi -mem 128G -e Price -reads $bam_1 -genomic $PRICE_anno_name -prefix $sample -plot
gedi -mem 128G -e ViewCIT -m bed -name 'd.getType()'  *.orfs.cit > $sample.orfs.bed
echo -e "PRICE done at $(date '+%Y-%m-%d %H:%M:%S')"

##### RibORF #####
echo -e "RibORF start at $(date '+%Y-%m-%d %H:%M:%S')"
mkdir -p $output_path/RibORF && cd $output_path/RibORF
RibORF_path=/home/user/data2/lit/software/RibORF/RibORF.2.0
# 使用ribocode鉴定得到的read length以及offset
config=../RiboCode/metaplots_pre_config.txt
grep "^#" $config|awk -v OFS='\t' '{print $2,$4+3}'|tail -n +3|head -n -1 > offset.correction.parameters.txt
# 根据offset文件校正
[ -f $sam ] || samtools view -h $bam_1 > $sam
time perl $RibORF_path/offsetCorrect.pl -r $sam -p offset.correction.parameters.txt -o corrected.$sample.sam
# 可视化校正后的结果
time perl $RibORF_path/readDist.pl -f corrected.$sample.sam -g $genePred -o ./ -d 1
# call ORFs
time perl $RibORF_path/ribORF.pl -f corrected.$sample.sam -c $candidateORF -o ./
echo -e "RibORF done at $(date '+%Y-%m-%d %H:%M:%S')"



