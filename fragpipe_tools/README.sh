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
time bash $script fragpipe.workflow fragpipe-files.fp-manifest ./ 30 &> test.log

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

########## 更换卡松Ribo-Seq参数之后的database，运行default-search ##########
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

##### 在F2上test生成mzML和mgf文件 #####
mkdir output/MS/Fragpipe_output/20241217-test
cd output/MS/Fragpipe_output/20241217-test
mkdir config && cd config
cp /rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/2024_12_03_batch_run/Trypsin_LysC/fragpipe.workflow ./
cp /rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/2024_12_13_default_loose_para_F2/Trypsin_LysC/fragpipe-files.fp-manifest ./
cd ..
bash /rd1/user/lit/project/sORFs/Uni.Fragpipe.sh \
config/fragpipe.workflow \
config/fragpipe-files.fp-manifest \
./ \
40

########## 在F2上test基于转录组的库 ##########
bash Uni.Batch.Fragpipe.Default.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_F2.csv \
$PWD/custom_database/Transcriptome_based/2024-12-18-decoys-contam-uniprot.nonCano.sorf.20241218.trans_based_specific.fa.fas \
$PWD/output/MS/Fragpipe_output/2024_12_18_default_trans_database_F2 &> log/Uni.Batch.Fragpipe.2024_12_18_default_trans_database_F2.log

#20241218，当库过大或者使用基于转录组的库的时候，无法同时使用nonspecific的消化方法进行搜库，因此需要额外制定null_set
for enzyme in Nocleavage;do
bash Uni.gene.enzyme.sh $enzyme config/Default/${enzyme}_6_50_prot_1.workflow config/Default/Trypsin_6_50_prot_1.workflow
done
bash Uni.Batch.Fragpipe.Default.v1.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_20241125.csv \
$PWD/custom_database/Transcriptome_based/2024-12-18-decoys-contam-uniprot.nonCano.sorf.20241218.trans_based.fa.fas \
$PWD/output/MS/Fragpipe_output/2024_12_18_default_trans_database \
nocleavage &> log/Uni.Batch.Fragpipe.2024_12_18_default_trans_database.log

########## F2 uniprot数据库 ##########
bash deprecated_script/Uni.Batch.Fragpipe.Default.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_F2.csv \
$PWD/custom_database/2024-11-21-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.fasta.fas \
$PWD/output/MS/Fragpipe_output/2025_01_07_default_uniprot_database_F2 &> log/Uni.Batch.Fragpipe.2025_01_07_default_uniprot_database_F2.log

bash deprecated_script/Uni.Batch.Fragpipe.Default.sh \
$PWD/raw_data/all_expr.txt \
$PWD/raw_data/MS_metadata_F2.csv \
$PWD/custom_database/2024-11-21-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.fasta.fas \
$PWD/output/MS/Fragpipe_output/2025_01_13_default_uniprot_database_F2_bak_1 &> log/Uni.Batch.Fragpipe.2025_01_13_default_uniprot_database_F2.log

########## 测试输入mzML运行fragpipe ##########
# 在windows上测试

########## default运行新样本,得到MSFragger输出的mzML文件 ##########
# 确认MS_metadata表格中的样本名和文件名可以匹配
find /rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/ -name "*.d" > raw_data/all_expr_20250114.txt
less raw_data/MS_metadata_20250106.csv | sed -n '26,91p' |cut -f 1 -d ','|xargs -i grep {} raw_data/all_expr_20250114.txt
# 查看酶的种类，出现了之前没有使用过的ArgC,AspN以及GluC,需要添加到Uni.gene.enzyme.sh中
less raw_data/MS_metadata_20250106.csv | sed -n '26,91p' | cut -f 2 -d ','|sort|uniq -c
# 修改Trypsin_6_50_prot_1.workflow，将write_calibrated_mzml和write_uncalibrated_mgf都改成true
for enzyme in Trypsin_ArgC Trypsin_AspN Trypsin_GluC Trypsin Trypsin_LysN Trypsin_LysC Trypsin_Chymotrypsin Nonspecific Nocleavage;do
bash Uni.gene.enzyme.v1.sh $enzyme config/Default/${enzyme}_6_50_prot_1.workflow config/Default/Trypsin_6_50_prot_1.workflow
done
# 运行新一批的数据
less raw_data/MS_metadata_20250106.csv | sed -n '26,91p' > raw_data/MS_metadata_20250114_new_batch.csv
bash Uni.Batch.Fragpipe.Default.v2.sh \
$PWD/raw_data/all_expr_20250114.txt \
$PWD/raw_data/MS_metadata_20250114_new_batch.csv \
$PWD/custom_database/2024-11-06-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.nonCano.sorf.filtered.fa.fas \
$PWD/output/MS/Fragpipe_output/2025_01_14_default_ribo_database_formal_data \
nonspecific &> log/Uni.Batch.Fragpipe.2025_01_14_default_ribo_database_formal_data.log

# 上面把人的数据也一起跑了，是不对的，需要重新跑
less raw_data/MS_metadata_20250106.csv | sed -n '26,91p' |grep ^E16 > raw_data/MS_metadata_20250120_new_batch_mouse.csv
bash Uni.Batch.Fragpipe.Default.v2.sh \
$PWD/raw_data/all_expr_20250114.txt \
$PWD/raw_data/MS_metadata_20250120_new_batch_mouse.csv \
$PWD/custom_database/2024-11-06-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.nonCano.sorf.filtered.fa.fas \
$PWD/output/MS/Fragpipe_output/2025_01_20_default_ribo_database_formal_data_mouse \
nonspecific &> log/Uni.Batch.Fragpipe.2025_01_20_default_ribo_database_formal_data_mouse.log

# 20250124 不知道为什么没跑成功，并且Default下面的workflow有些变成了空文件，需要重新生成一下
## 因为把trypsin也拉入循环了
cp output/MS/Fragpipe_output/2025_01_14_default_ribo_database_formal_data/Trypsin/fragpipe.workflow config/Default/Basic.workflow
# 修改Trypsin_6_50_prot_1.workflow，将write_calibrated_mzml和write_uncalibrated_mgf都改成true
for enzyme in Trypsin Trypsin_ArgC Trypsin_AspN Trypsin_GluC Trypsin_LysN Trypsin_LysC Trypsin_Chymotrypsin Nonspecific Nocleavage;do
bash Uni.gene.enzyme.v1.sh $enzyme config/Default/${enzyme}_6_50_prot_1.workflow config/Default/Basic.workflow
done

rm -rf output/MS/Fragpipe_output/2025_01_20_default_ribo_database_formal_data_mouse
less raw_data/MS_metadata_20250106.csv | sed -n '26,91p' |grep ^E16 > raw_data/MS_metadata_20250120_new_batch_mouse.csv
bash Uni.Batch.Fragpipe.Default.v2.sh \
$PWD/raw_data/all_expr_20250114.txt \
$PWD/raw_data/MS_metadata_20250120_new_batch_mouse.csv \
$PWD/custom_database/2024-11-06-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.nonCano.sorf.filtered.fa.fas \
$PWD/output/MS/Fragpipe_output/2025_01_24_default_ribo_database_formal_data_mouse \
nonspecific &> log/Uni.Batch.Fragpipe.2025_01_24_default_ribo_database_formal_data_mouse.log

########## 正式运行新一批数据 ##########

########## Default search ##########

# 20250206 修改config/Default/Basic.workflow中的phi-report.remove-contaminants为true
for enzyme in Trypsin Trypsin_ArgC Trypsin_AspN Trypsin_GluC Trypsin_LysN Trypsin_LysC Trypsin_Chymotrypsin Nonspecific Nocleavage;do
bash Uni.gene.enzyme.v1.sh $enzyme config/Default/${enzyme}_6_50_prot_1.workflow config/Default/Basic.workflow
done

mkdir -p output/MS/Fragpipe_output/Formal/mouse
bash Uni.Batch.Fragpipe.Default.v2.sh \
$PWD/raw_data/all_expr_20250114.txt \
$PWD/raw_data/MS_metadata_20250120_new_batch_mouse.csv \
$PWD/custom_database/Ribo_ORFs_add_assemble_20250125/2025-02-06-decoys-contam-uniprot.nonCano.sorf.20250125.fa.fas \
$PWD/output/MS/Fragpipe_output/Formal/mouse/2025_02_06_default_ribo_database \
nonspecific &> log/Uni.Batch.Fragpipe.2025_02_06_default_ribo_database_mouse.log
# 最后两个酶没有跑成功 [库输入错误，20250209重新跑]
bash Uni.Batch.Fragpipe.Default.v2.sh \
$PWD/raw_data/all_expr_20250114.txt \
$PWD/raw_data/MS_metadata_20250120_new_batch_mouse.csv \
$PWD/custom_database/Transcriptome_based_20250206/2025-02-06-decoys-contam-uniprot.nonCano.sorf.20250206.trans_based.fa.fas \
$PWD/output/MS/Fragpipe_output/Formal/mouse/2025_02_06_default_trans_database/ \
nocleavage &> log/Uni.Batch.Fragpipe.2025_02_06_default_trans_database_mouse.log
# 终端重新跑了一个，脚本重新跑一个
cat <(grep Trypsin_LysN $PWD/raw_data/MS_metadata_20250120_new_batch_mouse.csv) > $PWD/raw_data/MS_metadata_20250120_new_batch_mouse.Trypsin_LysN.csv
bash Uni.Batch.Fragpipe.Default.v2.sh \
$PWD/raw_data/all_expr_20250114.txt \
$PWD/raw_data/MS_metadata_20250120_new_batch_mouse.Trypsin_LysN.csv \
$PWD/custom_database/Transcriptome_based_20250206/2025-02-06-decoys-contam-uniprot.nonCano.sorf.20250206.trans_based.fa.fas \
$PWD/output/MS/Fragpipe_output/Formal/mouse/2025_02_06_default_trans_database/ \
nocleavage &> log/Uni.Batch.Fragpipe.2025_02_06_default_trans_database_mouse.Trypsin_LysN.log

########## Open search ##########
# 同时也修改phi-report.remove-contaminants为true
mv config/Open/6_50_prot_1.workflow config/Open/Basic.workflow
for enzyme in Trypsin Trypsin_ArgC Trypsin_AspN Trypsin_GluC Trypsin_LysN Trypsin_LysC Trypsin_Chymotrypsin Nonspecific Nocleavage;do
bash Uni.gene.enzyme.v1.sh $enzyme config/Open/${enzyme}_6_50_prot_1.workflow config/Open/Basic.workflow
done

# 跑Null的数据会报错，因此需要去掉Null的数据
grep -v Null $PWD/raw_data/MS_metadata_20250120_new_batch_mouse.csv > $PWD/raw_data/MS_metadata_20250120_new_batch_mouse.noNull.csv
bash Uni.Batch.Fragpipe.Open.v1.sh \
$PWD/raw_data/all_expr_20250114.txt \
$PWD/raw_data/MS_metadata_20250120_new_batch_mouse.noNull.csv \
$PWD/custom_database/Ribo_ORFs_add_assemble_20250125/2025-02-06-decoys-contam-uniprot.nonCano.sorf.20250125.fa.fas \
$PWD/output/MS/Fragpipe_output/Formal/mouse/2025_02_06_open_ribo_database \
nocleavage &> log/Uni.Batch.Fragpipe.2025_02_06_open_ribo_database_mouse.log

# [库输入错误，20250209重新跑]
bash Uni.Batch.Fragpipe.Open.v1.sh \
$PWD/raw_data/all_expr_20250114.txt \
$PWD/raw_data/MS_metadata_20250120_new_batch_mouse.csv \
$PWD/custom_database/Transcriptome_based_20250206/2025-02-06-decoys-contam-uniprot.nonCano.sorf.20250206.trans_based.fa.fas \
$PWD/output/MS/Fragpipe_output/Formal/mouse/2025_02_06_open_trans_database \
nocleavage &> log/Uni.Batch.Fragpipe.2025_02_06_open_trans_database_mouse.log

# 重新跑的两个
rm -rf $PWD/output/MS/Fragpipe_output/Formal/mouse/2025_02_06_default_trans_database/
rm -rf $PWD/output/MS/Fragpipe_output/Formal/mouse/2025_02_06_open_trans_database 
bash Uni.Batch.Fragpipe.Default.v2.sh \
$PWD/raw_data/all_expr_20250114.txt \
$PWD/raw_data/MS_metadata_20250120_new_batch_mouse.csv \
$PWD/custom_database/Transcriptome_based_20250209/2025-02-09-decoys-contam-uniprot.nonCano.sorf.20250209.trans_based.fa.fas \
$PWD/output/MS/Fragpipe_output/Formal/mouse/2025_02_09_default_trans_database/ \
nocleavage &> log/Uni.Batch.Fragpipe.2025_02_09_default_trans_database_mouse.log
wait
bash Uni.Batch.Fragpipe.Open.v1.sh \
$PWD/raw_data/all_expr_20250114.txt \
$PWD/raw_data/MS_metadata_20250120_new_batch_mouse.csv \
$PWD/custom_database/Transcriptome_based_20250209/2025-02-09-decoys-contam-uniprot.nonCano.sorf.20250209.trans_based.fa.fas \
$PWD/output/MS/Fragpipe_output/Formal/mouse/2025_02_09_open_trans_database \
nocleavage &> log/Uni.Batch.Fragpipe.2025_02_09_open_trans_database_mouse.log

# 20250216
egrep "Trypsin_GluC|Trypsin_LysN|Trypsin_LysC|Trypsin_Chymotrypsin" $PWD/raw_data/MS_metadata_20250120_new_batch_mouse.csv > $PWD/raw_data/MS_metadata_20250120_new_batch_mouse.4enz.csv
bash Uni.Batch.Fragpipe.Open.v1.sh \
$PWD/raw_data/all_expr_20250114.txt \
$PWD/raw_data/MS_metadata_20250120_new_batch_mouse.4enz.csv \
$PWD/custom_database/Transcriptome_based_20250209/2025-02-09-decoys-contam-uniprot.nonCano.sorf.20250209.trans_based.fa.fas \
$PWD/output/MS/Fragpipe_output/Formal/mouse/2025_02_16_open_trans_database \
nocleavage &> log/Uni.Batch.Fragpipe.2025_02_16_open_trans_database_mouse.log

# 整合下之前跑的open trans的结果
output_path=$PWD/output/MS/Fragpipe_output/Formal/mouse
## 把之前跑了一半中止的结果给去掉
rm -rf $output_path/2025_02_09_open_trans_database/Trypsin_Chymotrypsin
mkdir -p $output_path/2025_02_21_merge_open_trans_database
cp -r $output_path/2025_02_16_open_trans_database/* $output_path/2025_02_09_open_trans_database/* $output_path/2025_02_21_merge_open_trans_database

for enzyme in Trypsin Trypsin_ArgC Trypsin_AspN Trypsin_GluC Trypsin_LysN Trypsin_LysC Trypsin_Chymotrypsin Nonspecific Nocleavage;do
bash Uni.gene.enzyme.v1.sh $enzyme config/Group/${enzyme}_6_50_prot_1.workflow config/Group/Basic.workflow
done