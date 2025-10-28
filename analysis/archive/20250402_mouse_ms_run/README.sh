# 运行fragpipe
output_path=$PWD/output/db_search_20250402
project_path=/rd1/user/lit/project/sORFs
database_path=$project_path/custom_database/trans_based_database_new_20250325/2025-03-25-decoys-uniprot.contam.sorfs.trans.modiHeader.fa.fas
workflow_path=$project_path/config/Group
all_raw_files=$project_path/raw_data/all_expr_20250114.txt
metadata_path=$project_path/raw_data/MS_metadata_20250120_new_batch_mouse.csv
script=/rd1/user/lit/project/sORFs/analysis/20250326_human_ms_run/Uni.Batch.Fragpipe.v1.20250326.sh
mkdir -p $output_path
mkdir -p log
nohup bash $script \
$all_raw_files \
$metadata_path \
$database_path \
$output_path \
nocleavage \
$workflow_path &> log/Uni.Batch.Fragpipe.20250402.log &