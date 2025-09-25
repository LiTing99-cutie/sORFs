#!/usr/bin/sh

################################################
#File Name: 2.2a.Run.filter_and_merge.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 31 Mar 2025 11:47:04 AM CST
################################################

set -eo pipefail

# 使用这个脚本：/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/2.1.1.dataset.three_nt_peri.R，得到每个样本的有周期性的读长以及offset
input_path=$PWD/human_brain_output_20250227
output_path=$PWD/human_brain_ribo_merge_call_orfs_20250338
bash 2.2.1.Uni.filter_and_merge.sh \
    $input_path/Aligned.sortedByCoord.out.bam.lst \
    $input_path/stat/read_length_offset.txt \
    $output_path/merged_filtered_reads.sortedByCoord.bam \
    $output_path/temp_bam_1 &&
bash 2.2.1.Uni.filter_and_merge.sh \
    $input_path/Aligned.toTranscriptome.out.bam.lst \
    $input_path/stat/read_length_offset.txt \
    $output_path/merged_filtered_reads.toTranscriptome.bam \
    $output_path/temp_bam_2
