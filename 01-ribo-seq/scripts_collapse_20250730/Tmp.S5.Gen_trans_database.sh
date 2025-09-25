#!/usr/bin/sh

################################################
#File Name: Tmp.S4.Gen_trans_database.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri Dec 13 11:59:05 2024
################################################

set -eo pipefail

cd annot/RiboORF/mm39
faTrans -stop candidateORF.fa candidateORF.prot.fa
seqkit fx2tab candidateORF.prot.fa > candidateORF.prot.tab
seqkit seq -g -m 6 -M 150 candidateORF.prot.fa >  candidate.sORF.prot.fa
seqkit fx2tab candidate.sORF.prot.fa > candidate.sORF.prot.tab