#!/usr/bin/sh

################################################
#File Name: Uni.tmp.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 28 Feb 2025 04:32:51 PM CST
################################################

set -eo pipefail

# 根据sorf_id，例如ENSMUST00000002379.15+chr17:34066632-34067429
# 1.生成gpe和gtf
# 2.计算表达量和翻译效率

# ./output/MS_minus_sampled.txt
sorfs_id=$1
# ./output/MS_minus_sampled
output_path=$2
# 是否计算overlapped reads
overlapped=$3
mkdir -p $output_path

# 基于给定的ORF_id_trans，生成对应的gpe和gtf
project_path=/home/user/data3/lit/project/sORFs
### 根据参考基因组修改 ###
gpe=$project_path/01-ribo-seq/output/assembled_trans/stringtie_output/gencode.vM29.add_assemble.genePred.txt
fa=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
script_1=$project_path/03-Cross-anna/Uni.gen.genepred.f_sorf_id.py
script_2=$project_path/01-ribo-seq/analysis/te_calc_20250212/Uni.rna_reads.20250303.v2.sh
script_3=$project_path/01-ribo-seq/analysis/te_calc_20250212/Uni.cal.p_site.20250303.v2.sh
script_4=$project_path/01-ribo-seq/analysis/te_calc_20250212/Uni.Cal_expr_te.R
cd $output_path
mkdir -p log
python $script_1 \
	$sorfs_id \
	$gpe \
	gpe
genePredToGtf file gpe gtf
bash $script_2 gtf rna_counts $overlapped &> log/Cal.rna_reads.log
bash $script_3 gtf ribo_counts $overlapped &> log/Cal.ribo_reads.log
Rscript $script_4 "ribo_counts/p_site.counts.txt" "rna_counts/" \
	"$project_path/01-ribo-seq/analysis/te_calc_20250212/libsize/rna/libsize.txt" \
	"$project_path/01-ribo-seq/analysis/te_calc_20250212/libsize/ribo/libsize.txt" \
	"te"