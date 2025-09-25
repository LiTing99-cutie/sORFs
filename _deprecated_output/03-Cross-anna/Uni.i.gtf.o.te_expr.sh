#!/usr/bin/sh

################################################
#File Name: Uni.tmp.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 28 Feb 2025 04:32:51 PM CST
################################################

set -eo pipefail

# 对于已经注释的ORF，输入gtf
# 计算表达量和翻译效率

# gtf
gtf=$1
# ./output/MS_minus_sampled
output_path=$2
# 是否计算overlapped reads
overlapped=$3
mkdir -p $output_path

project_path=/home/user/data3/lit/project/sORFs
script_2=$project_path/01-ribo-seq/analysis/te_calc_20250212/Uni.rna_reads.20250303.v2.sh
script_3=$project_path/01-ribo-seq/analysis/te_calc_20250212/Uni.cal.p_site.20250303.v2.sh
script_4=$project_path/01-ribo-seq/analysis/te_calc_20250212/Uni.Cal_expr_te.R
cd $output_path
mkdir -p log
bash $script_2 $gtf rna_counts $overlapped &> log/Cal.rna_reads.log
bash $script_3 $gtf ribo_counts $overlapped &> log/Cal.ribo_reads.log
Rscript $script_4 "ribo_counts/p_site.counts.txt" "rna_counts/" \
	"$project_path/01-ribo-seq/analysis/te_calc_20250212/libsize/rna/libsize.txt" \
	"$project_path/01-ribo-seq/analysis/te_calc_20250212/libsize/ribo/libsize.txt" \
	"te"