#!/usr/bin/sh

################################################
#File Name: README.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2025年03月21日 星期五 20时37分53秒
################################################

# 分析目的：确认PSM的过滤方案
# 候选的过滤方案：
# Percolator
## Group-specific FDR 不同的阈值：1%，0.1%以及0.01%
### 阈值1%下过滤hyperscore
## Global FDR 卡不同的阈值：1%，3%以及5%
# PeptideProphet

# less /rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/Formal/mouse/2025_02_09_default_trans_database/Trypsin/CAD20241224licq_BSEP_DDA_60min_E16_1_PAGE_T_T_Slot2_49_1_4736_d
cd /rd1/user/lit/project/sORFs/analysis/20250321_PSM_filtering/
cp /rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/Formal/mouse/2025_02_09_default_trans_database/Trypsin/fragpipe.workflow ./fragpipe.workflow
cp /rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/Formal/mouse/2025_02_09_default_trans_database/Trypsin/fragpipe-files.fp-manifest ./fragpipe-files.fp-manifest

# screen
# 设置print_decoys以及print-decoys为true，看是否能在结果中输出decoy，从而可以计算group-specific FDR
mkdir 20250324_print_decoys
time bash /rd1/user/lit/project/sORFs/Uni.Fragpipe.sh fragpipe.workflow fragpipe-files.fp-manifest ./20250324_print_decoys 40

# screen
mkdir 20250324_group_specific && cd 20250324_group_specific
## 替换database
cp /rd1/user/lit/project/sORFs/config/Group/Trypsin_6_50_prot_1.workflow ./
cp ../fragpipe-files.fp-manifest ./
time bash /rd1/user/lit/project/sORFs/Uni.Fragpipe.sh Trypsin_6_50_prot_1.workflow fragpipe-files.fp-manifest ./ 40

cd ..
mkdir 20250324_group_specific_single_sample && cd 20250324_group_specific_single_sample
head -n 2 ../20250324_group_specific/fragpipe-files.fp-manifest|tail -n +2 > fragpipe-files.fp-manifest
## 修改decoys的相关参数
cp ../20250324_group_specific/Trypsin_6_50_prot_1.workflow ./
bash /rd1/user/lit/project/sORFs/Uni.Fragpipe.sh Trypsin_6_50_prot_1.workflow fragpipe-files.fp-manifest ./ 40 &> log.txt

mkdir 20250325_group_specific_single_sample_new_trans_db && cd 20250325_group_specific_single_sample_new_trans_db
cp ../20250324_group_specific_single_sample/fragpipe-files.fp-manifest ./
## 修改db
cp ../20250324_group_specific_single_sample/Trypsin_6_50_prot_1.workflow ./
bash /rd1/user/lit/project/sORFs/Uni.Fragpipe.sh Trypsin_6_50_prot_1.workflow fragpipe-files.fp-manifest ./ 40 &> log.txt