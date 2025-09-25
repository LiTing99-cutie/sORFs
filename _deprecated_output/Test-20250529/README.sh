mkdir log
echo /home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/organize_all_test_data_20250515/p21_40_0422_add.fq.gz >  raw_fastq.lst
nohup bash Uni.20250529.sh $PWD/raw_fastq.lst &> log/Run.20250529.log &
nohup bash Run.stat.v1.1.20250529.sh &> log/Run.stat.20250529.log &
nohup Rscript Run.organize.v1.20250528.R $PWD &> log/Run.organize.20250529.log &