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
sorfs_id=/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/S7/orf.id.txt
# ./output/MS_minus_sampled
output_path=$PWD/output/cnt
# 是否计算overlapped reads
overlapped=yes
mkdir -p $output_path

# 基于给定的ORF_id_trans，生成对应的gpe和gtf
define_annotation_gencode_v41_human
conda activate base
gen_genepred_i_sorf_id=$project_path/03-Cross-anna/Uni.gen.genepred.i_sorf_id.py
project_path=/home/user/data3/lit/project/sORFs
cd $output_path
mkdir -p log
python $gen_genepred_i_sorf_id \
	$sorfs_id \
	$gpe_15 \
	target.gpe
genePredToGtf file target.gpe target.gtf

# 基于生成的gtf以及bam lst，计算表达量（多重比对都计算入内）
script_cal_ribo=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/in_house_phase_I_data_20250520/Uni.cal.p_site.v3.20250522.sh
## 准备bam lst
project_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/in_house_phase_I_data_20250520
ls $project_path/output/alignment/human_brain_Aligned.sortedByCoord.out.bam \
	$project_path/output/bam/bam_1/merged.bam \
	$project_path/output/bam/bam_2/length.24-36.bam \
	$project_path/output/bam/bam_3/p_sites_0.5_1.sam \
	$project_path/output/bam/bam_4/p_sites_0.5_0.01.sam > bam.diff.level.lst
nohup bash $script_cal_ribo target.gtf ribo_counts $overlapped bam.diff.level.lst &> log/Cal.ribo_reads.log &

