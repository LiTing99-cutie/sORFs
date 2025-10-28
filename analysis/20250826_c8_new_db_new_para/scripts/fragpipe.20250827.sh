#!/bin/bash
# 运行fragpipe
output_path=$(realpath ../processed)
mkdir -p $output_path
project_path=/rd1/user/lit/project/sORFs
## 数据库路径
database_path=$project_path/custom_database/human/custom_db_20250826_v2/2025-08-26-decoys-human_brain_custom_db.fasta.fas
## 质谱文件
all_raw_files=$(realpath ../input/all_expr_files.c8.txt)
metadata_path=$(realpath ../input/metadata.c8.txt)
echo $all_raw_files $metadata_path
## 批量运行脚本
script=/rd1/user/lit/project/sORFs/analysis/20250623_human_ms_rerun/Uni.Batch.Fragpipe.v1.1.20250625.sh
## closed
closed_workflow_path=$project_path/fragpipe_config/LFQ_full
for enzyme in Trypsin_ArgC Trypsin_AspN Trypsin_GluC Trypsin Trypsin_LysN Trypsin_LysC Trypsin_Chymotrypsin;do
bash /rd1/user/lit/project/sORFs/Uni.gene.enzyme.v1.sh $enzyme $closed_workflow_path/${enzyme}_6_50_prot_1.workflow $closed_workflow_path/Basic.workflow
done
bash $script \
    $all_raw_files \
    $metadata_path \
    $database_path \
    $output_path/db_search_c8_closed_full \
    nocleavage \
    $closed_workflow_path

