# 运行fragpipe
output_path=$PWD/output/db_search_c8_open_search
project_path=/rd1/user/lit/project/sORFs
database_path=$project_path/custom_database/human/trans_based_database_human_20250326/2025-03-26-decoys-uniprot.contam.sorfs.trans.modiHeader.fa.fas
workflow_path=$project_path/config/Open
all_raw_files=$PWD/output/db_search_c8_open_search/all_expr_files.c8.txt
metadata_path=$PWD/output/db_search_c8_open_search/metadata.c8.txt
script=/rd1/user/lit/project/sORFs/analysis/20250623_human_ms_rerun/Uni.Batch.Fragpipe.v1.1.20250625.sh
echo $all_raw_files $metadata_path
mkdir -p $output_path
bash $script \
    $all_raw_files \
    $metadata_path \
    $database_path \
    $output_path \
    nocleavage \
    $workflow_path