#!/usr/bin/sh

################################################
#File Name: README.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2024年11月01日 星期五 15时42分47秒
################################################

set -eo pipefail

pushd /rd1/user/lit/software/fragpipe/bin

conda activate fragpipe
# 安装java17
# chmod +x /rd1/user/lit/project/sORFs/fragpipe_tools/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8

./fragpipe --headless \
--workflow /rd1/user/lit/project/sORFs/config/2024-11-01-1/fragpipe.workflow \
--manifest /rd1/user/lit/project/sORFs/config/2024-11-01-1/fragpipe-files.fp-manifest \
--workdir /rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/2024-11-01-1 \
--ram 100 \
--threads 30 \
--config-tools-folder /rd1/user/lit/project/sORFs/fragpipe_tools \
--config-diann /rd1/user/lit/project/sORFs/fragpipe_tools/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 \
--config-python ~/rd1/anaconda3/envs/py3/bin/python

popd

mkdir log
nohup bash Run.2024-11-01.sh > log/Run.2024-11-01.log 2>&1 &

bash Uni.Fragpipe.sh config/2024_11_06_1/fragpipe.workflow \
 config/2024_11_06_1/fragpipe-files.fp-manifest \
 output/MS/Fragpipe_output/2024_11_06_1 && bash Uni.Fragpipe.sh config/2024_11_06_2/fragpipe.workflow \
 config/2024_11_06_2/fragpipe-files.fp-manifest \
 output/MS/Fragpipe_output/2024_11_06_2

########## 命令行生成config文件 ##########
# 根据样本生成manifest文件
mkdir output/MS/Fragpipe_output/F2_rerun_new_para/config
cd output/MS/Fragpipe_output/F2_rerun_new_para/config
## 生成manifest文件
ls -d /rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/CAD20241022licq_BSEP_PreExp_DDA_60min/*.d |grep F2 | xargs -I {}  sh -c 'echo -e "{}\t$(basename {})\t\tDDA"' > fragpipe-files.fp-manifest
## 生成workflow文件
workflow=/rd1/user/lit/project/sORFs/config/workflow_LFQ_MBR_trypsin_lysc.workflow
new_database_path=/rd1/user/lit/project/sORFs/custom_database/2024-11-06-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.nonCano.sorf.filtered.fa.fas
sed 's|database.db-path=|database.db-path='$new_database_path'|' $workflow > fragpipe.workflow
## 运行fragpipe
script=/rd1/user/lit/project/sORFs/Uni.Fragpipe.sh
time bash $script fragpipe.workflow fragpipe-files.fp-manifest ../

########## 批量运行 ##########
dos2unix raw_data/MS_metadata_20241113.csv
cat <(ls -d /rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/CAD20241022licq_BSEP_PreExp_DDA_60min/*.d) <(ls -d /rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/CAD20241105licq_BSEP_DDA_60min/*.d) > raw_data/all_expr.txt
bash Uni.Batch.Fragpipe.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_20241113.csv \
$PWD/custom_database/2024-11-06-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.nonCano.sorf.filtered.fa.fas \
$PWD/output/MS/Fragpipe_output/2024_11_14_batch_run &> log/Uni.Batch.Fragpipe.2024_11_14.log

tail -6 $PWD/raw_data/MS_metadata_20241113.csv > $PWD/raw_data/MS_metadata_6_samples.csv
bash Uni.Batch.Fragpipe.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_6_samples.csv \
$PWD/custom_database/2024-11-06-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.nonCano.sorf.filtered.fa.fas \
$PWD/output/MS/Fragpipe_output/2024_11_18_batch_run_6_samples_test_overlap &> log/Uni.Batch.Fragpipe.2024_11_18.log

########## test不加新鉴定的小肽的结果 ##########
nohup bash Run.2024-11-19.sh > log/Run.2024-11-19.log 2>&1 &


########## Test 参数 ##########
nohup bash Run.2024-11-21.test.para.sh > log/Run.2024-11-21.test.para.log 2>&1 &

manifest=/rd1/user/lit/project/sORFs/config/2024-11-01-3-trypsin-lysc-uniprot-custom/fragpipe-files.fp-manifest
condition=trypsin_lysc_6_50_prot_1_no_MBR_no_ionquant
bash Uni.Fragpipe.sh config/workflow_LFQ_MBR_${condition}.workflow \
    $manifest \
    output/MS/Fragpipe_output/${condition} \
    40

########## Test PDV ##########
java -jar FP-PDV-1.2.2.jar your_result_folder_full_paththread_num
java -jar PDV-2.1.0.jar -h
java -jar FP-PDV-1.2.4.jar -h

pepXML=/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/LFQ_MBR_trypsin_lysc_6_50_prot_1_no_MBR_no_ionquant/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562.pepXML
mzML=/rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/CAD20241022licq_BSEP_PreExp_DDA_60min/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562_uncalibrated.mzML
java -jar PDV-1.1.0/PDV-1.1.0.jar -r input_data/SF_200217_U2OS_TiO2_HCD_OT_rep1_myrimatch_mzML.pepXML -rt 2 -s input_data/SF_200217_U2OS_TiO2_HCD_OT_rep1.mzML -st 2 -i input_data/spectrum_scan_number.txt -k s -o output -a 0.05 -c 3 -pw 1 -fw 800 -fh 400 -fu px -ft pdf
java -jar /rd1/user/lit/software/PDV-2.1.0/PDV-2.1.0.jar -r $pepXML -rt 2 -s $mzML -st 2 -k s -o output -a 0.05 -c 3 -pw 1 -fw 800 -fh 400 -fu px -ft pdf
java -jar /rd1/user/lit/software/PDV-2.1.0/PDV-2.1.0.jar -h

########## Test nonspecific ##########
nohup bash Run.2024-11-26.test.para.sh > log/Run.2024-11-26.test.para.log 2>&1 &

########## 统计长度分布，分割fasta文件 ##########
cd /rd1/user/lit/project/sORFs/custom_database
uniprot_fasta=/rd1/user/lit/project/sORFs/custom_database/uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.fasta
# <=30 15 seqs
seqkit seq -M 30 $uniprot_fasta | seqkit stat
# >30 & <=50 35 seqs
seqkit seq -m 31 -M 50 $uniprot_fasta | seqkit stat
# >50 & <=100 620 seqs
seqkit seq -m 51 -M 100 $uniprot_fasta | seqkit stat
# >100 & <=150 1,433 seqs
seqkit seq -m 101 -M 150 $uniprot_fasta | seqkit stat
# <=150 2,103 seqs
seqkit seq -M 150 $uniprot_fasta | seqkit stat

cp uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.nonCano.sorf.filtered.fa segmentation/anno_ribo_sorfs
cp nonCano.sorf.filtered.fa segmentation/ribo_sorfs

cd segmentation/anno_ribo_sorfs
seqkit seq -M 150 uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.nonCano.sorf.filtered.fa > anno_ribo_sorfs_150.fa
seqkit seq -M 100 anno_ribo_sorfs_150.fa > anno_ribo_sorfs_100.fa
seqkit seq -M 50 anno_ribo_sorfs_150.fa > anno_ribo_sorfs_50.fa
seqkit seq -M 30 anno_ribo_sorfs_150.fa > anno_ribo_sorfs_30.fa

cd ../../segmentation/ribo_sorfs
cp nonCano.sorf.filtered.fa ribo_sorfs_150.fa
seqkit seq -M 100 ribo_sorfs_150.fa > ribo_sorfs_100.fa
seqkit seq -M 50 ribo_sorfs_150.fa > ribo_sorfs_50.fa
seqkit seq -M 30 ribo_sorfs_150.fa > ribo_sorfs_30.fa

########## Test input fasta ##########
nohup bash Run.2024-11-26.test.input_fasta.sh > log/Run.2024-11-26.test.input_fasta.log 2>&1 &


########## 运行C8和酸沉淀的结果 ##########
cd output/MS/Fragpipe_output/C8_acid
script=/rd1/user/lit/project/sORFs/Uni.Fragpipe.sh
workflow=/rd1/user/lit/project/sORFs/config/workflow_LFQ_MBR_trypsin.workflow
time bash $script $workflow fragpipe-files.fp-manifest ./ 30

########## 测试two-step ##########
path=/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/F1_6_50_prot_1_non_specific_no_quant_uniprot_step_1
cd $path
time bash $script fragpipe.workflow fragpipe-files.fp-manifest ./ 30

########## 修改参数，重新运行 ##########
# 基于这个workflow修改参数配置，在windows生成，同时也使用Uni.gene.enzyme.sh生成，以便下次新的配置修改
less /rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/Default_F2_prot_1/fragpipe.workflow

find /rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/ -name "*.d" > raw_data/all_expr.txt
bash Uni.Batch.Fragpipe.Default.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_20241125.csv \
$PWD/custom_database/2024-11-06-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.nonCano.sorf.filtered.fa.fas \
$PWD/output/MS/Fragpipe_output/2024_12_03_batch_run &> log/Uni.Batch.Fragpipe.2024_12_03.log

# sample名没对应上文件名，重新运行下trypsin的结果
awk -F',' '$2=="Trypsin"' $PWD/raw_data/MS_metadata_20241125.csv > $PWD/raw_data/MS_metadata_20241125_trypsin.csv
bash Uni.Batch.Fragpipe.Default.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_20241125_trypsin.csv \
$PWD/custom_database/2024-11-06-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.nonCano.sorf.filtered.fa.fas \
$PWD/output/MS/Fragpipe_output/2024_12_05_batch_run &> log/Uni.Batch.Fragpipe.2024_12_05.log

cp -r $PWD/output/MS/Fragpipe_output/2024_12_05_batch_run/Trypsin $PWD/output/MS/Fragpipe_output/2024_12_03_batch_run/

########## 查看uncalibrated.mzML ##########
mzml=/rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/CAD20241120licq_BSEP_C8_PCP_DDA_60min/CAD20241120licq_BSEP_DDA_60min_m14_C8_Slot2-1_1_4073_uncalibrated.mzML
grep -o 'scan=[0-9]*' $mzml \
| awk -F= '{print $2}' \
| sort -n | tail -1

grep -o 'spectrum index="[0-9]*"' $mzml \
| awk -F'"' '{print $2}' \
| sort -n | tail -1

grep "spectrumList count" $mzml

########## 运行open-search ##########
for enzyme in Trypsin Trypsin_LysN Trypsin_LysC Trypsin_Chymotrypsin Nonspecific;do
bash Uni.gene.enzyme.sh $enzyme config/Open/${enzyme}_6_50_prot_1.workflow config/Open/6_50_prot_1.workflow
done

bash Uni.Batch.Fragpipe.Open.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_20241125.csv \
$PWD/custom_database/2024-11-06-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.nonCano.sorf.filtered.fa.fas \
$PWD/output/MS/Fragpipe_output/2024_12_06_batch_run &> log/Uni.Batch.Fragpipe.2024_12_06.log

cat <(head -n1 $PWD/raw_data/MS_metadata_20241125.csv) <(grep F2 $PWD/raw_data/MS_metadata_20241125.csv) > $PWD/raw_data/MS_metadata_F2.csv
bash Uni.Batch.Fragpipe.Open.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_F2.csv \
$PWD/custom_database/2024-11-06-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.nonCano.sorf.filtered.fa.fas \
$PWD/output/MS/Fragpipe_output/2024_12_06_batch_run &> log/Uni.Batch.Fragpipe.2024_12_06.log

########## 更换database，运行default-search ##########
bash Uni.Batch.Fragpipe.Default.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_20241125.csv \
$PWD/custom_database/Ribo_ORFs_loose_para_20241206/2024-12-13-decoys-contam-uniprot.nonCano.sorf.20241213.loose_para.fa.fas \
$PWD/output/MS/Fragpipe_output/2024_12_13_default_loose_para &> log/Uni.Batch.Fragpipe.2024_12_13_default_loose_para.log
# 这个失败了

# 1、test只跑F2，work
bash Uni.Batch.Fragpipe.Default.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_F2.csv \
$PWD/custom_database/Ribo_ORFs_loose_para_20241206/2024-12-13-decoys-contam-uniprot.nonCano.sorf.20241213.loose_para.fa.fas \
$PWD/output/MS/Fragpipe_output/2024_12_13_default_loose_para_F2 &> log/Uni.Batch.Fragpipe.2024_12_13_default_loose_para_F2.log

# 2、test只跑F1，不work
cat <(head -n1 $PWD/raw_data/MS_metadata_20241125.csv) <(grep F1 $PWD/raw_data/MS_metadata_20241125.csv) > $PWD/raw_data/MS_metadata_F1.csv
bash Uni.Batch.Fragpipe.Default.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_F1.csv \
$PWD/custom_database/Ribo_ORFs_loose_para_20241206/2024-12-13-decoys-contam-uniprot.nonCano.sorf.20241213.loose_para.fa.fas \
$PWD/output/MS/Fragpipe_output/2024_12_13_default_loose_para_F1 &> log/Uni.Batch.Fragpipe.2024_12_13_default_loose_para_F1.log

# 3、test只跑F1，修改酶切设置
    workdir=/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/2024_12_13_default_loose_para_F1/Null/
    mkdir -p $workdir/config
    bash Uni.gene.enzyme.sh Trypsin $workdir/config/fragpipe.trypsin.workflow $workdir/fragpipe.workflow
    cp $workdir/fragpipe-files.fp-manifest $workdir/config/fragpipe-files.fp-manifest
    bash Uni.Fragpipe.sh \
        $workdir/config/fragpipe.trypsin.workflow \
        $workdir/config/fragpipe-files.fp-manifest \
        $workdir \
        40
# 跑到一半kill了，发现是原始数据路径下有从msconvert得到mzML，导致fragpipe识别错误，把路径下的*mzML以及*mgf移动到convert文件夹下，重新运行F1发现是ok的

# 4、debug之后重新运行
bash Uni.Batch.Fragpipe.Default.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_20241125.csv \
$PWD/custom_database/Ribo_ORFs_loose_para_20241206/2024-12-13-decoys-contam-uniprot.nonCano.sorf.20241213.loose_para.fa.fas \
$PWD/output/MS/Fragpipe_output/2024_12_13_default_loose_para &> log/Uni.Batch.Fragpipe.2024_12_13_default_loose_para.log