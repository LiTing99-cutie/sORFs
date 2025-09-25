#!/usr/bin/sh

################################################
#File Name: 2.0.Uni.Build_annotation.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 03 Mar 2025 03:44:11 PM CST
################################################

set -eo pipefail


fa=/home/user/data3/lit/resource/genome/human/hg38/hg38.fa
gtf=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
genePred=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gpe
RiboCode_annot_path=annotation/RiboCode_annot
RibORF_script_path=/home/user/data2/lit/software/RibORF/RibORF.2.0
RibORF_annot_path=annotation/RibORF_annot
PRICE_anno_name=hg38_gencv41

mkdir -p $RiboCode_annot_path
mkdir -p $RibORF_annot_path

# ribocode
source activate ribocode
prepare_transcripts -g $gtf -f $fa -o $RiboCode_annot_path

# riborf
perl $RibORF_script_path/ORFannotate.pl -g $fa -t $genePred -s ATG/CTG/ACG/GTG/TTG/ATA/ATC/ATT/AAG/AGG -l 21 -o $RibORF_annot_path

# price
gedi -e IndexGenome -s $fa -a $gtf -n $PRICE_anno_name