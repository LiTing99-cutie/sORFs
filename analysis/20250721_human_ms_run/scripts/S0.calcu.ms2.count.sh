mkdir -p ../results/S0
# 统计所有的二级谱图的数量
calculate_ms2count_script=Uni.calculate.ms2count.v1.20250525.sh
mzML_path=/rd1/user/lit/project/sORFs/raw_data/In_House_Guomics_SP_MSdata/CAD20250514licq_BSEP_DDA_60min_rename
find $mzML_path -name "*_uncalibrated.mzML" > $mzML_path/mzML_list.txt
bash $calculate_ms2count_script $mzML_path/mzML_list.txt ../results/S0/ms2_count_49.txt
mzML_path=/rd1/user/lit/project/sORFs/raw_data/In_House_Guomics_SP_MSdata/CAD20241224licq_BSEP_DDA_60min
find $mzML_path -name "*_uncalibrated.mzML"|grep 21pcw > $mzML_path/mzML_list_human.txt
bash $calculate_ms2count_script $mzML_path/mzML_list_human.txt ../results/S0/ms2_count_33.txt
cat ../results/S0/ms2_count_49.txt ../results/S0/ms2_count_33.txt > ../results/S0/ms2_count.txt