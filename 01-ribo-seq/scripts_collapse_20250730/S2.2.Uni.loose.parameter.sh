#!/usr/bin/sh

################################################
#File Name: Uni.loose.parameter.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri Dec  6 21:22:50 2024
################################################

set -eo pipefail

# 对于每一个样本结果路径进行循环
sample_output_path=$1
Ribo_ORFs_output_path_name=$2
output_path=${sample_output_path}/${Ribo_ORFs_output_path_name}
sample=$(echo $sample_output_path | awk -F'/' '{print $(NF-1)}')

# 基因组及其注释
fa=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
genePred=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.genePred.txt
gtf=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.gtf
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/RiboCode_annot/mm39
STARindex=/home/user/data3/lit/resource/genome/mouse/mm39/index/mm39_STARindex
# 增加9种替代的起始密码子
candidateORF=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/RiboORF/mm39/candidateORF.genepred.txt
bam=$sample_output_path/alignment/*_Aligned.toTranscriptome.out.bam
PRICE_bam=$sample_output_path/Ribo-ORFs/PRICE/*_Aligned.sortedByCoord.out.bam
RibORF_sam=$sample_output_path/Ribo-ORFs/RibORF/$sample.sam

##### RiboCode #####
echo -e "RiboCode start at $(date '+%Y-%m-%d %H:%M:%S')"
source activate ribocode
mkdir -p $output_path/RiboCode && cd $output_path/RiboCode
metaplots -a $RiboCode_annot -r $bam && \
# 增加9种替代的起始密码子，降低p value的cutoff
RiboCode --alt_start_codons CTG,ACG,GTG,TTG,ATA,ATC,ATT,AAG,AGG --pval-cutoff 0.1 -a $RiboCode_annot -c metaplots_pre_config.txt -l no -g -o ./$sample --output-gtf --output-bed --min-AA-length 6
echo -e "RiboCode done at $(date '+%Y-%m-%d %H:%M:%S')"

##### PRICE #####
# 输出的结果已经包括了10种起始密码子以及所有p value值的结果
# echo -e "PRICE start at $(date '+%Y-%m-%d %H:%M:%S')"
# mkdir -p $output_path/PRICE && cd $output_path/PRICE
# gedi -e Price -reads $PRICE_bam -genomic mm39_gencvM29 -prefix $sample -plot
# gedi -e ViewCIT -m bed -name 'd.getType()'  *.orfs.cit > $sample.orfs.bed
# echo -e "PRICE done at $(date '+%Y-%m-%d %H:%M:%S')"

##### RibORF #####
echo -e "RibORF start at $(date '+%Y-%m-%d %H:%M:%S')"
mkdir -p $output_path/RibORF && cd $output_path/RibORF
RibORF_path=/home/user/data2/lit/software/RibORF/RibORF.2.0

# 0. build annotation; see /home/user/data3/lit/project/sORFs/01-ribo-seq/Build_annotation.sh

# 1. bam2sam

# 2. output
# 使用ribocode鉴定得到的read length以及offset
config=../RiboCode/metaplots_pre_config.txt
grep "^#" $config|awk -v OFS='\t' '{print $2,$4+3}'|tail -n +3|head -n -1 > offset.correction.parameters.txt
# 根据offset文件校正
time perl $RibORF_path/offsetCorrect.pl -r $RibORF_sam -p offset.correction.parameters.txt -o corrected.$sample.sam
# 可视化校正后的结果
time perl $RibORF_path/readDist.pl -f corrected.$sample.sam -g $genePred -o ./ -d 1
# call ORFs
time perl $RibORF_path/ribORF.pl -f corrected.$sample.sam -c $candidateORF -o ./
echo -e "RibORF done at $(date '+%Y-%m-%d %H:%M:%S')"

