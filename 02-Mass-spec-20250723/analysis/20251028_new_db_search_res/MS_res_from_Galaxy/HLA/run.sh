# 得到HLA分型
## hlahd
R1=/home/user/data/lit/project/sORFs/07-Genome/rawdata/L1EJF1602305-p21_Gen2Seq.R1.raw.fastq.gz
R2=/home/user/data/lit/project/sORFs/07-Genome/rawdata/L1EJF1602305-p21_Gen2Seq.R2.raw.fastq.gz
# /home/user/data2/lit/software/hlahd.1.7.1/bin/hlahd.sh \
#     -t 8 \              # 线程数
#     -m 100 \            # 最小比对质量
#     -f freq_data \      # 等位基因频率数据
#     $R1 \   # Read 1
#     $R2 \   # Read 2
#     gene_split \        # HLA基因分割数据
#     dictionary \        # HLA等位基因字典（包括HLA-Y等）
#     human_brain\       # 样本名
#     ./ # 输出目录
freq_data=/home/user/data2/lit/software/hlahd.1.7.1/freq_data
dictionary=/home/user/data2/lit/software/hlahd.1.7.1/dictionary
script=/home/user/data2/lit/software/hlahd.1.7.1/bin/hlahd.sh
$script -t 8 -m 100 -f $freq_data $R1 $R2 gene_split $dictionary human_brain ./
## OptiTypePipeline.py
nohup bash run.optitype.sh &
## arcasHLA
nohup bash run.arcasHLA.sh &
# netMHCpan预测
netMHCpan \
  -p MHC-I-peptide.seq.tsv \
  -a HLA-A11:01,HLA-C07:02,HLA-B40:01 \
  -BA \
  -xls -xlsfile netMHCpan-I.xls > netMHCpan-I.log 2>&1
cd JingJie
awk 'length($0) >= 8 && length($0) <= 12' HLA-I.pep.uniq.txt > HLA-I.pep.uniq.8-12.txt
netMHCpan \
  -p HLA-I.pep.uniq.8-12.txt \
  -a HLA-A11:01,HLA-C07:02,HLA-B40:01 \
  -BA \
  -xls -xlsfile netMHCpan-I.xls > netMHCpan-I.log 2>&1
