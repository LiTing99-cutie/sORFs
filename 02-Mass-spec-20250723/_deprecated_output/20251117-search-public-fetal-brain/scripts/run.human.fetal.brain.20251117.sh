# 运行fragpipe For regular MS file
dir=$(realpath ../processed)
output_path=$dir/fetal_brain/output/
mkdir -p $output_path/closed $output_path/open
## 数据库路径
database_path=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251117-search-public-fetal-brain/custom_db/2025-10-24-decoys-contam-candidateORF.filtered.rmdup.renamed.addContam.pep.fa.fas
script_path=$PWD/Uni.Batch.Fragpipe.v2.20250501.forPublic.sh
all_raw_files_1=/home/user/data3/xieyn/project/peptidome/public_MS_searching/fetal_brain/test/fetal_brain_files.txt
cat $all_raw_files_1 |grep -v Fetal > $output_path/all_raw_files.txt
all_raw_files=$output_path/all_raw_files.txt
[ -f $output_path/metadata.txt ] && rm -rf $output_path/metadata.txt
cat $all_raw_files | xargs -n 1 basename -s ".raw" | awk '{print $0",Trypsin"}' >> $output_path/metadata.txt
metadata_path=$output_path/metadata.txt
workflow_path=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251117-search-public-fetal-brain/fragpipe_config/regular_MS

mkdir -p $output_path
bash $script_path \
    $all_raw_files \
    $metadata_path \
    $database_path \
    $output_path/closed \
    nocleavage \
    $workflow_path/LFQ

bash $script_path \
    $all_raw_files \
    $metadata_path \
    $database_path \
    $output_path/open \
    nocleavage \
    $workflow_path/Open