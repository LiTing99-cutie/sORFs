#!/usr/bin/sh

################################################
#File Name: Uni.loose.parameter.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri Dec  6 21:22:50 2024
################################################

set -eo pipefail

source activate biotools

output_path=$1
mkdir -p $output_path && cd $output_path
Ribocode_bam=$2
PRICE_bam=$3

# 基因组及其注释
genePred=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/results/custom.gtf.with_orf.15.gpe
gtf=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/results/custom.gtf.with_orf.gtf
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250813_demo_data_analysis/processed/annotation/RiboCode_annot
candidateORF=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.genepred.txt
PRICE_anno_name=hg38_custom
sample=Human_brain
f0_percent_cutoff=0.6
PVALUE_CUTOFF=0.01
PVAL_CUTOFF=1

##### RiboCode #####
#   -pv1 PVALUE1_CUTOFF, --pvalue1_cutoff PVALUE1_CUTOFF
#                         pvalue cutoff of frame0 > frame2 for automatically
#                         predicting P-site location, default 0.001
#   -pv2 PVALUE2_CUTOFF, --pvalue2_cutoff PVALUE2_CUTOFF
#                         pvalue cutoff of frame0 > frame2 for automatically
#                         predicting P-site location, default 0.001
#   -p PVAL_CUTOFF, --pval-cutoff PVAL_CUTOFF
#                         P-value cutoff for ORF filtering, default 0.05
# echo -e "RiboCode start at $(date '+%Y-%m-%d %H:%M:%S')"
# source activate ribocode
# mkdir -p $output_path/RiboCode && cd $output_path/RiboCode
# metaplots -a $RiboCode_annot -r $Ribocode_bam -f0_percent $f0_percent_cutoff -pv1 $PVALUE_CUTOFF -pv2 $PVALUE_CUTOFF && \
# RiboCode --alt_start_codons CTG,ACG,GTG,TTG --pval-cutoff $PVAL_CUTOFF -a $RiboCode_annot -c metaplots_pre_config.txt -l no -g -o ./$sample --output-gtf --output-bed --min-AA-length 6
# echo -e "RiboCode done at $(date '+%Y-%m-%d %H:%M:%S')"

##### PRICE #####
# 输出的结果已经包括了10种起始密码子以及所有p value值的结果
echo -e "PRICE start at $(date '+%Y-%m-%d %H:%M:%S')"
mkdir -p $output_path/PRICE && cd $output_path/PRICE
# 先清除文件夹下的文件，再运行
# 大bam需要调整内存
# rm -rf ./*
[ -f ${PRICE_bam}.bai ] || samtools index -@ 30 $PRICE_bam
gedi -mem 128G -e Price -reads $PRICE_bam -genomic $PRICE_anno_name -prefix $sample -plot
# gedi -mem 128G -e ViewCIT -m bed -name 'd.getType()'  *.orfs.cit > $sample.orfs.bed
echo -e "PRICE done at $(date '+%Y-%m-%d %H:%M:%S')"

##### RibORF #####
echo -e "RibORF start at $(date '+%Y-%m-%d %H:%M:%S')"
mkdir -p $output_path/RibORF && cd $output_path/RibORF
RibORF_path=/home/user/data2/lit/software/RibORF/RibORF.2.0
# 使用ribocode鉴定得到的read length以及offset
config=../RiboCode/metaplots_pre_config.txt
grep "^#" $config|awk -v OFS='\t' '{print $2,$4+3}'|tail -n +3|head -n -1 > offset.correction.parameters.txt
# 根据offset文件校正
RibORF_sam="${PRICE_bam%.bam}.sam"
[ -f $RibORF_sam ] || samtools view -@ 30 -h $PRICE_bam > $RibORF_sam
time perl $RibORF_path/offsetCorrect.pl -r $RibORF_sam -p offset.correction.parameters.txt -o corrected.$sample.sam
# 可视化校正后的结果
time perl $RibORF_path/readDist.pl -f corrected.$sample.sam -g $genePred -o ./ -d 1
# call ORFs
time perl $RibORF_path/ribORF.pl -f corrected.$sample.sam -c $candidateORF -o ./
echo -e "RibORF done at $(date '+%Y-%m-%d %H:%M:%S')"

