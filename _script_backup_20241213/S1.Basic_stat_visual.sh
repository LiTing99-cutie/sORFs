#!/usr/bin/sh

################################################
#File Name: S1.Basic_stat_visual.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2024年12月13日 星期五 20时03分10秒
################################################

set -eo pipefail

Rscript S1.Uni.Organize_all_samples.R \
    output/MS/Fragpipe_output/2024_12_13_default_loose_para \
    output/MS/stat/2024_12_13_default_loose_para
Rscript S2.Uni.Visualize_all_samples.R \
    output/MS/stat/2024_12_13_default_loose_para \
    output/MS/visual/2024_12_13_default_loose_para 

Rscript S1.Uni.Organize_all_samples.R \
    output/MS/Fragpipe_output/2024_12_13_default_loose_para \
    output/MS/stat/2024_12_13_default_loose_para
Rscript S2.Uni.Visualize_all_samples.R \
    output/MS/stat/2024_12_13_default_loose_para \
    output/MS/visual/2024_12_13_default_loose_para     