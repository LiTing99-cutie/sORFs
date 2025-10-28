# 运行fragpipe
output_path=$PWD/output/db_search_20250523
project_path=/rd1/user/lit/project/sORFs
database_path=$project_path/custom_database/human/trans_based_database_human_20250326/2025-03-26-decoys-uniprot.contam.sorfs.trans.modiHeader.fa.fas
workflow_path=$project_path/config/Group
all_raw_files=$project_path/raw_data/all_expr_files.20250523.txt
metadata_path=$project_path/raw_data/MS_metadata_20250523_human.txt
script=/rd1/user/lit/project/sORFs/analysis/20250326_human_ms_run/Uni.Batch.Fragpipe.v1.20250326.sh
echo $all_raw_files $metadata_path
mkdir -p $output_path
mkdir -p log
nohup bash $script \
$all_raw_files \
$metadata_path \
$database_path \
$output_path \
nocleavage \
$workflow_path &> log/Uni.Batch.Fragpipe.20250523.log &

# 内存不够，重新运行着两个酶
output_path=$PWD/output/db_search_20250523
egrep 'Trypsin_LysC|Trypsin_LysN' $metadata_path > $project_path/raw_data/MS_metadata_20250523_human.2enzyme.txt
project_path=/rd1/user/lit/project/sORFs
database_path=$project_path/custom_database/human/trans_based_database_human_20250326/2025-03-26-decoys-uniprot.contam.sorfs.trans.modiHeader.fa.fas
workflow_path=$project_path/config/Group
all_raw_files=$project_path/raw_data/all_expr_files.20250523.txt
metadata_path=$project_path/raw_data/MS_metadata_20250523_human.2enzyme.txt
script=/rd1/user/lit/project/sORFs/analysis/20250326_human_ms_run/Uni.Batch.Fragpipe.v1.20250326.sh
nohup bash $script \
$all_raw_files \
$metadata_path \
$database_path \
$output_path \
nocleavage \
$workflow_path &> log/Uni.Batch.Fragpipe.20250525.log &