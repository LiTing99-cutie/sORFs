# 质谱的raw data路径
# 整理所有的all_expr_files（所有的实验文件）以及metadata（元数据）
human_raw_data_path=/rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/human_organized_20250625
ls -d $human_raw_data_path/*.d > $human_raw_data_path/metadata/all_expr_files.txt
cat /rd1/user/lit/project/sORFs/raw_data/MS_metadata_20250326_new_batch_human.csv \
    /rd1/user/lit/project/sORFs/raw_data/MS_metadata_20250523_human.txt > $human_raw_data_path/metadata/metadata.txt
egrep -v "less3K|21pcw_4" $human_raw_data_path/metadata/metadata.txt > $human_raw_data_path/metadata/metadata.formal.txt
less $human_raw_data_path/metadata/metadata.formal.txt|cut -d, -f1-2,3 > $human_raw_data_path/metadata/metadata.formal.v1.txt

# 拷贝脚本
cp /rd1/user/lit/project/sORFs/analysis/20250326_human_ms_run/Uni.Batch.Fragpipe.v1.20250326.sh ./

# 测试c8的open搜库
grep 21pcw_1_C8_T_T $human_raw_data_path/metadata/metadata.formal.v1.txt > output/db_search_c8_open_search/metadata.c8.txt
grep 21pcw_1_C8_T_T $human_raw_data_path/metadata/all_expr_files.txt > output/db_search_c8_open_search/all_expr_files.c8.txt
mkdir log
nohup bash Run.c8.open.20250625.sh &> log/Uni.Batch.Fragpipe.$(date '+%Y%m%d').log &

