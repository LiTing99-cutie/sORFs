#!/usr/bin/sh

################################################
#File Name: Run.search.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2025年02月24日 星期一 15时31分05秒
################################################

set -eo pipefail

source activate fragpipe
proj_path=/rd1/user/lit/project/sORFs/
# 提取出用trypsin消化的样本
grep ,Trypsin, $proj_path/raw_data/MS_metadata_20250120_new_batch_mouse.csv > output/search/MS_metadata_20250120_new_batch_mouse_Trypsin.csv
# 库路径
database=$PWD/output/2025-02-24-decoys-contam-uniprot.trypsin_all_method_sep_unique.fa.fas
# 样本元数据路径
sample_file=$PWD/output/search/MS_metadata_20250120_new_batch_mouse_Trypsin.csv
# 所有原始数据路径
all_file=$proj_path/raw_data/all_expr_20250114.txt
# 搜库输出路径
output_path=$PWD/output/search
bash $proj_path/Uni.Batch.Fragpipe.Default.v2.sh \
    $all_file \
    $sample_file \
    $database \
    $output_path/default/ \
    nonspecific &> log/default.log
bash $proj_path/Uni.Batch.Fragpipe.Open.v1.sh \
    $all_file \
    $sample_file \
    $database \
    $output_path/open/ \
    nocleavage &> log/open.log