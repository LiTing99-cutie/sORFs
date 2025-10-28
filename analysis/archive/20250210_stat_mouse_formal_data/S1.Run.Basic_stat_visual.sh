#!/usr/bin/sh

################################################
#File Name: S1.Basic_stat_visual.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2024年12月13日 星期五 20时03分10秒
################################################

set -eo pipefail

proj_path=/rd1/user/lit/project/sORFs
## 计算ms2的数量
pushd $proj_path
[ -f output/MS/stat/sample_ms2_count_formal_hm.txt ] && rm -rf output/MS/stat/sample_ms2_count_formal_hm.txt
find raw_data/MS/Guomics_SP_MSdata/CAD20241224licq_BSEP_DDA_60min/ -name "*uncalibrated.mzML" | \
xargs -I {} grep -c "spectrum index=" {} >> output/MS/stat/sample_ms2_count_formal_hm.txt
mv output/MS/stat/sample_ms2_count_formal_hm.txt output/MS/stat/sample_ms2_count_formal_hm.txt.tmp
find raw_data/MS/Guomics_SP_MSdata/CAD20241224licq_BSEP_DDA_60min/ -name "*uncalibrated.mzML" | \
xargs -I {} echo {} | grep -oP '(?<=min_)(.*?)(?=_Slot)' >> output/MS/stat/sample_name_formal_hm.txt
paste output/MS/stat/sample_name_formal_hm.txt output/MS/stat/sample_ms2_count_formal_hm.txt.tmp \
 > output/MS/stat/sample_ms2_count_formal_hm.txt
rm -rf output/MS/stat/sample_name_formal_hm.txt output/MS/stat/sample_ms2_count_formal_hm.txt.tmp
cat <(head -n1 output/MS/stat/sample_ms2_count.txt) output/MS/stat/sample_ms2_count_formal_hm.txt > output/MS/stat/sample_ms2_count_formal_hm.add_h.txt
popd

mkdir -p output/stat/ output/visual/
# default_ribo 是直接从终端开始跑的
# default_trans
Rscript S1a1.Uni.Organize_all_samples.R \
    $proj_path/output/MS/Fragpipe_output/Formal/mouse/2025_02_09_default_trans_database/ \
    output/stat/default_trans \
    output/sample_metadata.rds
Rscript S1b1.Uni.Visualize_all_samples.R \
    output/stat/default_trans \
    output/visual/default_trans \
    output/sample_order.rds \
    output/sample_metadata.rds
# open_ribo
Rscript S1a1.Uni.Organize_all_samples.R \
    $proj_path/output/MS/Fragpipe_output/Formal/mouse/2025_02_06_open_ribo_database/ \
    output/stat/open_ribo \
    output/sample_metadata.rds
Rscript S1b1.Uni.Visualize_all_samples.R \
    output/stat/open_ribo \
    output/visual/open_ribo \
    output/sample_order.rds \
    output/sample_metadata.rds
# open_trans
cd /rd1/user/lit/project/sORFs/analysis/20250210_stat_mouse_formal_data
for run in open_trans;do
Rscript S1a1.Uni.Organize_all_samples.R \
    $proj_path/output/MS/Fragpipe_output/Formal/mouse/2025_02_21_merge_open_trans_database/ \
    output/stat/$run \
    output/sample_metadata.rds
echo "done"
Rscript S1b1.Uni.Visualize_all_samples.R \
    output/stat/$run \
    output/visual/$run \
    output/sample_order.rds \
    output/sample_metadata.rds
done