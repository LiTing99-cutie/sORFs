#!/usr/bin/sh

################################################
#File Name: README.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2025年05月04日 星期日 21时20分25秒
################################################

set -eo pipefail

# 运行fragpipe For regular MS file
output_path=$PWD/output/
project_path=/rd1/user/lit/project/sORFs
database_path=$project_path/custom_database/human/trans_based_database_human_20250326/2025-03-26-decoys-uniprot.contam.sorfs.trans.modiHeader.fa.fas
script_path=/rd1/user/lit/project/sORFs/analysis/20250501_publicData_database_searching/20250501_2022-NN_search/Uni.Batch.Fragpipe.v2.20250501.forPublic.sh
raw_path=/rd1/user/lit/project/sORFs/raw_data/public_brain_MS/non-HLA/PXD000561
enzyme=Trypsin
ls $raw_path/*raw > raw_files.path.txt
all_raw_files=$PWD/raw_files.path.txt
[ -f metadata.txt ] && rm -rf metadata.txt
awk -F/ '{gsub(/\.raw$/, "", $NF); print $NF ",'$enzyme'"}' raw_files.path.txt > metadata.txt
metadata_path=$PWD/metadata.txt
workflow_path=/rd1/user/lit/project/sORFs/analysis/20250501_publicData_database_searching/20250501_2022-NN_search/config/

mkdir -p $output_path
mkdir -p log
nohup bash $script_path \
    $all_raw_files \
    $metadata_path \
    $database_path \
    $output_path \
    nocleavage \
    $workflow_path &> log/Uni.Batch.Fragpipe.20250504.log &

