#!/usr/bin/sh

################################################
#File Name: S1.Basic_stat_visual.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2024年12月13日 星期五 20时03分10秒
################################################

set -eo pipefail

Rscript S1a1.Uni.Organize_all_samples.R \
    output/MS/Fragpipe_output/2024_12_13_default_loose_para \
    output/MS/stat/2024_12_13_default_loose_para
Rscript S1b1.Uni.Visualize_all_samples.R \
    output/MS/stat/2024_12_13_default_loose_para \
    output/MS/visual/2024_12_13_default_loose_para 

Rscript S1a1.Uni.Organize_all_samples.R \
    output/MS/Fragpipe_output/2024_12_03_batch_run/ \
    output/MS/stat/2024_12_03_batch_run
Rscript S1b1.Uni.Visualize_all_samples.R \
    output/MS/stat/2024_12_03_batch_run \
    output/MS/visual/2024_12_03_batch_run     

Rscript S1a1.Uni.Organize_all_samples.R \
    output/MS/Fragpipe_output/2024_12_18_default_trans_database_F2 \
    output/MS/stat/2024_12_18_default_trans_database_F2
Rscript S1b1.Uni.Visualize_all_samples.R \
    output/MS/stat/2024_12_18_default_trans_database_F2 \
    output/MS/visual/2024_12_18_default_trans_database_F2    

Rscript S1a1.Uni.Organize_all_samples.R \
    output/MS/Fragpipe_output/2025_01_14_default_ribo_database_formal_data \
    output/MS/stat/2025_01_14_default_ribo_database_formal_data
Rscript S1b1.Uni.Visualize_all_samples.R \
    output/MS/stat/2025_01_14_default_ribo_database_formal_data \
    output/MS/visual/2025_01_14_default_ribo_database_formal_data 

# 运行2025_01_14_default_ribo_database_formal_data的时候需要更改一下代码
## 首先需要提前计算ms2的数量
find raw_data/MS/Guomics_SP_MSdata/CAD20241224licq_BSEP_DDA_60min/ -name "*uncalibrated.mzML" | \
xargs -I {} grep -c "spectrum index=" {} >> output/MS/stat/sample_ms2_count_formal_hm.txt
mv output/MS/stat/sample_ms2_count_formal_hm.txt output/MS/stat/sample_ms2_count_formal_hm.txt.tmp
find raw_data/MS/Guomics_SP_MSdata/CAD20241224licq_BSEP_DDA_60min/ -name "*uncalibrated.mzML" | \
xargs -I {} echo {} | grep -oP '(?<=min_)(.*?)(?=_Slot)' >> output/MS/stat/sample_name_formal_hm.txt
paste output/MS/stat/sample_name_formal_hm.txt output/MS/stat/sample_ms2_count_formal_hm.txt.tmp \
 > output/MS/stat/sample_ms2_count_formal_hm.txt
rm -rf output/MS/stat/sample_name_formal_hm.txt output/MS/stat/sample_ms2_count_formal_hm.txt.tmp
## 其次，由于里面包含了人的样本，所以这部分的结果需要去除

# 20250124
Rscript S1a1.Uni.Organize_all_samples.R \
    output/MS/Fragpipe_output/2025_01_20_default_ribo_database_formal_data_mouse \
    output/MS/stat/2025_01_20_default_ribo_database_formal_data_mouse \
    output/MS/stat/sample_ms2_count_formal_hm.txt
Rscript S1b1.Uni.Visualize_all_samples.R \
    output/MS/stat/2025_01_20_default_ribo_database_formal_data_mouse \
    output/MS/visual/2025_01_20_default_ribo_database_formal_data_mouse 