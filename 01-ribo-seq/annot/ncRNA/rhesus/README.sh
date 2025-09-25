#!/usr/bin/sh

################################################
#File Name: README.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 11 Oct 2024 05:19:26 PM CST
################################################

set -eo pipefail

# 20241011
mkdir ensembl gtRNAdb merged SILVA
wget https://ftp.ensembl.org/pub/current_fasta/macaca_mulatta/ncrna/Macaca_mulatta.Mmul_10.ncrna.fa.gz -P ensembl
wget --no-check-certificate https://gtrnadb.ucsc.edu/genomes/eukaryota/Mmula8/rheMac8-tRNAs.tar.gz -P gtRNAdb
