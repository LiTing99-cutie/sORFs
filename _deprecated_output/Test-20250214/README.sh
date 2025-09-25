#!/usr/bin/sh

################################################
#File Name: README.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 14 Feb 2025 10:04:16 AM CST
################################################

set -eo pipefail

qc_script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/te_calc_20250212/Uni.qc.mouse.rna_seq.sh
fastq_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Mouse_E16_test/cleandata
mkdir output
bash $qc_script $fastq_path/E16_N30_RPF.R2.clean.fastq.gz none output/test yes no no

# 自行去接头
bash $qc_script $fastq_path/E16_N30_RPF.R2.clean.fastq.gz none output/test no yes yes

bash $qc_script $fastq_path/E16_N30_RPF.R2.clean.fastq.gz AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC output/test no yes yes

cd output/test/E16_N30_RPF.R2.clean
trim_galore --adapter "file:/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/E16_mouse_data_test_2025_02_14/multiple_adapters.fa" -j 8 -q 20 --length 20 $fastq_path/E16_N30_RPF.R2.clean.fastq.gz --gzip -o output/trimmed_fastq &> log/trim_galore.log
fastqc -o output/fastqc -t 10 /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/E16_mouse_data_test_2025_02_14/output/test/E16_N30_RPF.R2.clean/output/trimmed_fastq/E16_N30_RPF.R2.clean_trimmed.fq.gz &> log/trimmed_fastqc.log

# 修改adapter为输入文件
cp /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/te_calc_20250212/Uni.qc.mouse.rna_seq.sh ./

bash Uni.qc.mouse.rna_seq.sh $fastq_path/E16_N30_RPF.R1.clean.fastq.gz $PWD/multiple_adapters.fa output/test yes yes yes

bash Uni.qc.mouse.rna_seq.sh $fastq_path/E16_Y60_RPF.R2.clean.fastq.gz $PWD/multiple_adapters.1.fa output/test yes yes yes
bash Uni.qc.mouse.rna_seq.sh $fastq_path/E16_Y60_RPF.R1.clean.fastq.gz $PWD/multiple_adapters.1.fa output/test yes yes yes

bash Uni.qc.mouse.rna_seq.sh $fastq_path/E16_Y60_RPF.R2.clean.fastq.gz $PWD/multiple_adapters.3.fa output/test yes yes yes

# 20250515重新分析
mkdir -p output/test_20250515
mkdir -p output/test_20250515_1
fastq_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Mouse_E16_test/cleandata
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG+GGGGGGGGGGGGGGGG
bash Uni.qc.mouse.rna_seq.sh $fastq_path/E16_N30_RPF.R1.clean.fastq.gz $PWD/adapters.20250515.fa output/test_20250515 no yes yes
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG
bash Uni.qc.mouse.rna_seq.sh $fastq_path/E16_N30_RPF.R1.clean.fastq.gz $PWD/adapters.20250515.1.fa output/test_20250515_1 no yes yes
