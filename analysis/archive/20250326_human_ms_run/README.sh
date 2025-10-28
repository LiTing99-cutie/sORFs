
# 修改批量运行fragpipe的脚本
cp /rd1/user/lit/project/sORFs/Uni.Batch.Fragpipe.Default.v2.sh ./Uni.Batch.Fragpipe.v1.20250326.sh

# 生成不同酶切条件下的workflow文件
workflow_path=/rd1/user/lit/project/sORFs/config/Group
for enzyme in Trypsin Trypsin_ArgC Trypsin_AspN Trypsin_GluC Trypsin_LysN Trypsin_LysC Trypsin_Chymotrypsin Nonspecific Nocleavage;do
bash /rd1/user/lit/project/sORFs/Uni.gene.enzyme.v1.sh $enzyme $workflow_path/${enzyme}_6_50_prot_1.workflow $workflow_path/Basic.workflow
done

# 生成人数据的metadata文件
cd /rd1/user/lit/project/sORFs/raw_data/
less MS_metadata_20250114_new_batch.csv |grep ^21pcw > MS_metadata_20250326_new_batch_human.csv

# 运行fragpipe
output_path=$PWD/output/db_search_20250326
project_path=/rd1/user/lit/project/sORFs
database_path=$project_path/custom_database/human/trans_based_database_human_20250326/2025-03-26-decoys-uniprot.contam.sorfs.trans.modiHeader.fa.fas
workflow_path=$project_path/config/Group
all_raw_files=$project_path/raw_data/all_expr_20250114.txt
metadata_path=$project_path/raw_data/MS_metadata_20250326_new_batch_human.csv
mkdir -p $output_path
mkdir -p log
nohup bash Uni.Batch.Fragpipe.v1.20250326.sh \
$all_raw_files \
$metadata_path \
$database_path \
$output_path \
nocleavage \
$workflow_path &> log/Uni.Batch.Fragpipe.20250326.log &

# 服务器重启
mkdir tmp
egrep "Trypsin_Chymotrypsin|Trypsin_GluC|Trypsin_LysN|Trypsin_LysC" $metadata_path > tmp/MS_metadata_20250326_new_batch_human_4_enzymes.csv
nohup bash Uni.Batch.Fragpipe.v1.20250326.sh \
$all_raw_files \
$PWD/tmp/MS_metadata_20250326_new_batch_human_4_enzymes.csv \
$database_path \
$output_path \
nocleavage \
$workflow_path &> log/Uni.Batch.Fragpipe.20250328.log &