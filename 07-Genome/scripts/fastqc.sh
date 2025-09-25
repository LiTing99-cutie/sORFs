# 1.1 对原始数据进行质控
mkdir -p ../results/fastqc_results
source activate biotools
fastqc /home/user/data3/lit/project/sORFs/07-Genome/rawdata/L1EJF1602305-p21_Gen2Seq.R1.raw.fastq.gz \
       /home/user/data3/lit/project/sORFs/07-Genome/rawdata/L1EJF1602305-p21_Gen2Seq.R2.raw.fastq.gz \
       -o ../results/fastqc_results