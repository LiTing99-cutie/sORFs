#!/bin/bash
# 运行fragpipe
output_path=$(realpath ../processed)
mkdir -p $output_path
project_path=/rd1/user/lit/project/sORFs
## 数据库路径
database_path=$project_path/custom_database/human/custom_db_20251024_ribo_filtered/2025-10-24-decoys-contam-candidateORF.filtered.rmdup.renamed.addContam.pep.fa.fas
## 质谱文件
all_raw_files=$project_path/raw_data/In_House_Guomics_SP_MSdata/human_organized_20250625/metadata/all_expr_files.txt
metadata_path=$project_path/raw_data/In_House_Guomics_SP_MSdata/human_organized_20250625/metadata/metadata.formal.v1.txt
echo $all_raw_files $metadata_path
## 批量运行脚本
script=$project_path/analysis/20250623_human_ms_rerun/Uni.Batch.Fragpipe.v1.1.20250625.sh
## closed
closed_workflow_path=$project_path/fragpipe_config/LFQ
for enzyme in Trypsin_ArgC Trypsin_AspN Trypsin_GluC Trypsin Trypsin_LysN Trypsin_LysC Trypsin_Chymotrypsin;do
bash $project_path/Uni.gene.enzyme.v1.sh $enzyme $closed_workflow_path/${enzyme}_6_50_prot_1.workflow $closed_workflow_path/Basic.workflow
done
bash $script \
    $all_raw_files \
    $metadata_path \
    $database_path \
    $output_path/db_search_c8_closed \
    nocleavage \
    $closed_workflow_path
## open
open_workflow_path=$project_path/fragpipe_config/Open
for enzyme in Trypsin_ArgC Trypsin_AspN Trypsin_GluC Trypsin Trypsin_LysN Trypsin_LysC Trypsin_Chymotrypsin;do
bash /rd1/user/lit/project/sORFs/Uni.gene.enzyme.v1.sh $enzyme $open_workflow_path/${enzyme}_6_50_prot_1.workflow $open_workflow_path/Basic.workflow
done
bash $script \
    $all_raw_files \
    $metadata_path \
    $database_path \
    $output_path/db_search_c8_open \
    nocleavage \
    $open_workflow_path
