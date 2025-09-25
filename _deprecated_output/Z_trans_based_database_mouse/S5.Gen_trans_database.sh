#!/usr/bin/sh

################################################
#File Name: Tmp.S4.Gen_trans_database.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri Dec 13 11:59:05 2024
################################################

set -eo pipefail

cd /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/mouse_run_20250125/annotation/RibORF_annot
faTrans -stop candidateORF.fa candidateORF.prot.fa
seqkit fx2tab candidateORF.prot.fa > candidateORF.prot.tab
