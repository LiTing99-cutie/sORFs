#!/usr/bin/sh

################################################
#File Name: S1_build_annotation.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Sat 25 Jan 2025 03:44:17 PM CST
################################################

set -eo pipefail

mkdir -p annotation/RiboCode_annot
mkdir -p annotation/RibORF_annot

# ribocode
source activate ribocode
fa=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
gtf=/home/user/data3/lit/project/sORFs/01-ribo-seq/output/assembled_trans/stringtie_output/gencode.vM29.add_assemble.gtf
RiboCode_annot_path=annotation/RiboCode_annot
prepare_transcripts -g $gtf -f $fa -o $RiboCode_annot_path

# riborf
genePred=/home/user/data3/lit/project/sORFs/01-ribo-seq/output/assembled_trans/stringtie_output/gencode.vM29.add_assemble.genePred.txt
RibORF_script_path=/home/user/data2/lit/software/RibORF/RibORF.2.0
RibORF_annot_path=annotation/RibORF_annot
perl $RibORF_script_path/ORFannotate.pl -g $fa -t $genePred -s ATG/CTG/ACG/GTG/TTG/ATA/ATC/ATT/AAG/AGG -l 21 -o $RibORF_annot_path

# price
gedi -e IndexGenome -s $fa -a $gtf -n mm39_gencvM29_add_s