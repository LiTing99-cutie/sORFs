# 运行fragpipe
output_path=$PWD/output/
project_path=/rd1/user/lit/project/sORFs
database_path=$project_path/custom_database/human/trans_based_database_human_20250326/2025-03-26-decoys-uniprot.contam.sorfs.trans.modiHeader.fa.fas

cp /rd1/user/lit/project/sORFs/analysis/20250326_human_ms_run/Uni.Batch.Fragpipe.v1.20250326.sh Uni.Batch.Fragpipe.v2.20250501.forPublic.sh
script_path=$PWD/Uni.Batch.Fragpipe.v2.20250501.forPublic.sh

ls /rd1/user/lit/project/sORFs/raw_data/public_brain_MS/non-HLA/PXD035950/*raw > raw_files.path.txt
all_raw_files=$PWD/raw_files.path.txt

[ -f metadata.txt ] && rm -rf metadata.txt
awk -F/ '{gsub(/\.raw$/, "", $NF); print $NF ",Trypsin"}' raw_files.path.txt > metadata.txt
metadata_path=$PWD/metadata.txt

mkdir config
cp /rd1/user/lit/project/sORFs/config/Group/Trypsin_6_50_prot_1.workflow ./config
# 更改两行
# workflow.input.data-type.im-ms=false
# workflow.input.data-type.regular-ms=true
workflow_path=$PWD/config/

mkdir -p $output_path
mkdir -p log
source activate fragpipe
nohup bash $script_path \
    $all_raw_files \
    $metadata_path \
    $database_path \
    $output_path \
    nocleavage \
    $workflow_path &> log/Uni.Batch.Fragpipe.20250501.log &

# uniprot database
cd /rd1/user/lit/project/sORFs/analysis/20250501_publicData_database_searching/20250501_2022-NN_search/
database_path=/rd1/user/lit/project/sORFs/custom_database/human/uniprot/2025-05-02-decoys-contam-uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.fasta.fas
output_path=$PWD/output/uniprot_db
all_raw_files=$PWD/raw_files.path.txt
metadata_path=$PWD/metadata.txt
# 将msfragger.group_variable=0更改为0
workflow_path=$PWD/config_1/
script_path=$PWD/Uni.Batch.Fragpipe.v2.20250501.forPublic.sh
mkdir -p $output_path
nohup bash $script_path \
    $all_raw_files \
    $metadata_path \
    $database_path \
    $output_path \
    nocleavage \
    $workflow_path &> log/Uni.Batch.Fragpipe.20250501.uniprot_db.log &

# HLA数据集
cd /rd1/user/lit/project/sORFs/analysis/20250501_publicData_database_searching/20250501_2022-NN_search/
workflow=/rd1/user/lit/project/sORFs/config/HLA/hla_NEW.workflow
project_path=/rd1/user/lit/project/sORFs
database_path=$project_path/custom_database/human/trans_based_database_human_20250326/2025-03-26-decoys-uniprot.contam.sorfs.trans.modiHeader.fa.fas
script=/rd1/user/lit/project/sORFs/Uni.Fragpipe.sh
output_path=$PWD/output/brain_HLA_20250501
mkdir -p $output_path && cd $output_path
mkdir config
sed 's|database.db-path=.*|database.db-path='"$database_path"'|' $workflow > config/fragpipe.workflow
ls /rd1/user/lit/project/sORFs/raw_data/public_brain_MS/HLA/PXD019643/*raw > raw_files.path.brain.HLA.txt
less raw_files.path.brain.HLA.txt| xargs -I {} sh -c 'echo -e "{}\t$(basename {})\t\tDDA"' > config/fragpipe-files.fp-manifest
nohup bash $script config/fragpipe.workflow config/fragpipe-files.fp-manifest ./ 40 &> ../../log/human.hla.log &
# 更改文件输入的个数
less config/fragpipe-files.fp-manifest|head -n1 > config/fragpipe-files.h1.fp-manifest
nohup bash $script config/fragpipe.workflow config/fragpipe-files.h1.fp-manifest ./ 40 &> ../../log/human.hla.h1.log &
# 更改msfragger.misc.slice-db=10
nohup bash $script config/fragpipe.workflow config/fragpipe-files.h1.fp-manifest ./ 40 &> ../../log/human.hla.h1.1.log &
# 更改msfragger.misc.slice-db=100
nohup bash $script config/fragpipe.workflow config/fragpipe-files.h1.fp-manifest ./ 40 &> ../../log/human.hla.h1.1.log &
# 更改msfragger.calibrate_mass=O，msfragger.misc.slice-db=10
mkdir test_10
nohup bash $script config/fragpipe.workflow config/fragpipe-files.h1.fp-manifest ./test_10 40 &> ../../log/human.hla.h1.split10.log &
# 更改msfragger.calibrate_mass=O，msfragger.misc.slice-db=5000
mkdir test
nohup bash $script config/fragpipe.workflow config/fragpipe-files.h1.fp-manifest ./test 40 &> ../../log/human.hla.h1.cali.0.split5000.log &
# 更改msfragger.calibrate_mass=O，msfragger.misc.slice-db=100
mkdir test_100
nohup bash $script config/fragpipe.workflow config/fragpipe-files.h1.fp-manifest ./test 40 &> ../../log/human.hla.h1.cali.0.split100.log &
# 更改database, msfragger.calibrate_mass=O, msfragger.misc.slice-db=10 【运行成功】
database_path=/rd1/user/lit/project/sORFs/custom_database/human/uniprot/2025-05-02-decoys-contam-uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.fasta.fas
sed 's|database.db-path=.*|database.db-path='"$database_path"'|' config/fragpipe.workflow > config/fragpipe.uniprot.db.workflow
nohup bash $script $PWD/config/fragpipe.uniprot.db.workflow $PWD/config/fragpipe-files.h1.fp-manifest $PWD 40 &> ../../log/human.hla.h1.cali.0.split10.uniprot.db.log &
nohup bash $script $PWD/config/fragpipe.uniprot.db.workflow $PWD/config/fragpipe-files.h1.fp-manifest $PWD 40 &> ../../log/human.hla.h1.cali.0.NoSplit.uniprot.db.log &

# 更改msfragger.calibrate_mass=O，msfragger.misc.slice-db=100
mkdir test_100_all_files
nohup bash $script config/fragpipe.workflow config/fragpipe-files.fp-manifest ./test_100_all_files 40 &> ../../log/human.hla.cali.0.split100.all_files.log &

# 更改msfragger.calibrate_mass=O，msfragger.misc.slice-db=10
mkdir test_10_all_files
nohup bash $script config/fragpipe.workflow config/fragpipe-files.fp-manifest ./test_10_all_files 40 &> ../../log/human.hla.cali.0.split10.all_files.log &

# 更改msfragger.calibrate_mass=O，msfragger.misc.slice-db=20
mkdir split_20_all_files
nohup bash $script config/fragpipe.workflow config/fragpipe-files.fp-manifest ./split_20_all_files 40 &> ../../log/human.hla.cali.0.split20.all_files.log &

#删除上述没有运行成功的test*文件夹输出
rm -rf test*
#### 拆分文件进行输入 ####
# 进入 manifest 文件所在目录
cd config/
# 按20行拆分文件，生成 fragpipe-files.fp-manifest.part-aa, part-ab, etc.
split -l 20 fragpipe-files.fp-manifest fragpipe-files.fp-manifest.part-
# 返回原目录
cd ..

# 赋予执行权限
chmod +x run_split_commands.sh

mkdir log
# 运行脚本（先预览命令，确认无误后去掉 --dry-run）
./run_split_commands.sh --dry-run  # 预览
nohup ./run_split_commands.sh &> log/human.hla.cali.0.split20.log &

#### 拆分文件进行输入，修改组特异性FDR控制，20个run需要5个小时，一共需要一天半左右的时间 ####
nohup ./run_split_commands.sh &> log/human.hla.cali.0.split20.log &

uniprot_fasta=/rd1/user/lit/project/sORFs/custom_database/human/uniprot/uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.fasta
seqkit seq -M 150 $uniprot_fasta |seqkit seq -n |awk '{print $1}'> stat_output/Cano.sep.id.txt
