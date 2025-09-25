# 整合长度分布文件
cat /home/user/data3/lit/project/sORFs/_deprecated_output/Test-20250509/01-output/stat/reads_length_distribution.txt \
/home/user/data3/lit/project/sORFs/_deprecated_output/Test-20250606/01-output/stat/reads_length_distribution.txt > \
../processed/reads_length_distribution.txt
sed 's/.raw//g' ../processed/reads_length_distribution.txt |grep -E 'p21_0321|p21_40_0422|p21_0523_1|p21_0523_5'|\
awk '$2>24 && $2<35' > ../processed/reads_length_distribution.clean.txt
# 整合周期性分布文件
cat /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/figures/batch_1/qual_assess/stat.all.txt <(tail -n +2 /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250813_demo_data_analysis/figures/batch_1/qual_assess/stat.all.txt) > ../processed/stat.all.txt
cp /home/user/data3/yaoc/translation_model/ribo-seq/01_reads/ribotish.frame_distr_with_peroidicity.tsv ../processed/



