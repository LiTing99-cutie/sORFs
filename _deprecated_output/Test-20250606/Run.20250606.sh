ls /home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/organize_all_test_data_20250515/p21_0523*.fq.gz > raw_fastq.lst
script_dir=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250529
DATE=$(date '+%Y%m%d')
bash $script_dir/Uni.20250529.sh $PWD/raw_fastq.lst &> log/Run.basic.$DATE.log
bash $script_dir/Run.stat.v1.1.20250529.sh &> log/Run.stat.$DATE.log
Rscript $script_dir/Run.organize.v1.20250528.R $PWD &> log/Run.organize.$DATE.log