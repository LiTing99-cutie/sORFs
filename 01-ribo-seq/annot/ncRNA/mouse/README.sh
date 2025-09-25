#!/usr/bin/sh

################################################
#File Name: README.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 11 Oct 2024 05:03:37 PM CST
################################################

set -eo pipefail

# 20241011

# 1. 从数据库中下载序列
mkdir ensembl gtRNAdb merged SILVA
wget https://ftp.ensembl.org/pub/current_fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz -P ensembl
wget --no-check-certificate https://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc39/mm39-tRNAs.tar.gz -P gtRNAdb
tar -xvzf SILVA/LSU_r138.2.RefNR.*.tgz && mv arb-silva* SILVA/LSU_r138.2.RefNR.fa
tar -xvzf SILVA/SSU_r138.2.RefNR.*.tgz && mv arb-silva* SILVA/SSU_r138.2.RefNR.fa
tar -xvzf ./gtRNAdb/*-tRNAs.tar.gz -C gtRNAdb
ensembl_ncRNA=./ensembl/Mus_musculus.GRCm39.ncrna.fa
ensembl_singleLine_ncRNA=./ensembl/Mus_musculus.GRCm39.ncrna.singleLine.fa
ensembl_singleLine_target_ncRNA=./ensembl/Mus_musculus.GRCm39.ncrna.singleLine.target.fa
gunzip -c ./ensembl/*.ncrna.fa.gz > $ensembl_ncRNA

## 2. 从ensembl数据库中提取出rRNA,tRNA,snoRNA
cp ../human/hg38/name.list ./ 
less $ensembl_ncRNA | grep "^>" | awk '{print $5}' |sort |uniq -c
#   23820 gene_biotype:lncRNA
#    2206 gene_biotype:miRNA
#     562 gene_biotype:misc_RNA
#       2 gene_biotype:Mt_rRNA
#      22 gene_biotype:Mt_tRNA
#      22 gene_biotype:ribozyme
#     354 gene_biotype:rRNA
#      51 gene_biotype:scaRNA
#       1 gene_biotype:scRNA
#    1507 gene_biotype:snoRNA
#    1381 gene_biotype:snRNA
#       2 gene_biotype:sRNA
seqkit seq -w 0 $ensembl_ncRNA > $ensembl_singleLine_ncRNA
grep -A1 -f name.list $ensembl_singleLine_ncRNA > $ensembl_singleLine_target_ncRNA
less $ensembl_singleLine_target_ncRNA | grep "^>" | awk '{print $5}' |sort |uniq -c
#       2 gene_biotype:Mt_rRNA
#      22 gene_biotype:Mt_tRNA
#     354 gene_biotype:rRNA
#    1507 gene_biotype:snoRNA

# 3. 合并
LSU=SILVA/LSU_r138.2.RefNR.fa
SSU=SILVA/SSU_r138.2.RefNR.fa
tRNA=gtRNAdb/mm39-tRNAs.fa
cat <(egrep -A1 "gene_biotype:Mt_rRNA|gene_biotype:rRNA" $ensembl_singleLine_ncRNA) $LSU $SSU > merged/mouse.rRNA.fa
cat <(egrep -A1 "gene_biotype:Mt_tRNA" $ensembl_singleLine_ncRNA) $tRNA > merged/mouse.tRNA.fa
egrep -A1 "gene_biotype:snoRNA" $ensembl_singleLine_ncRNA > merged/mouse.snoRNA.fa

# 4. 建立index
cd merged
conda activate biotools
bowtie2-build mouse.rRNA.fa mouse.rRNA
bowtie2-build mouse.tRNA.fa mouse.tRNA
bowtie2-build mouse.snoRNA.fa mouse.snoRNA
