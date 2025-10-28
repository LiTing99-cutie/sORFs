mkdir -p stat_output/S0
# 统计所有的二级谱图的数量
calculate_ms2count_script=/rd1/user/lit/project/sORFs/analysis/20250523_human_ms_run/Uni.calculate.ms2count.v1.20250525.sh
mzML_path=/rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/CAD20250514licq_BSEP_DDA_60min_rename
find $mzML_path -name "*_uncalibrated.mzML" > $mzML_path/mzML_list.txt
bash $calculate_ms2count_script $mzML_path/mzML_list.txt stat_output/S0/ms2_count.txt
# 统计
S1-S5
# 合并
mkdir stat_output/S6
cat stat_output/S4/psm_cano_sep_all.txt <(tail -n +2 stat_output/S5/psm_uncano_sep_all.txt) > stat_output/S6/psm_sep_all.txt
