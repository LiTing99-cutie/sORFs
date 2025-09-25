#!/usr/bin/sh

################################################
#File Name: README.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 30 May 2025 10:10:02 AM CST
################################################

set -eo pipefail
path=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/
find $path -name "p21*_Aligned.sortedByCoord.out.bam" > in_house_bam.lst

samtools merge -@ 20 -o p21_40_0425_Aligned.sortedByCoord.out.bam $(grep p21_40_0425 in_house_bam.lst)
samtools merge -@ 20 -o p21_40_0422_Aligned.sortedByCoord.out.bam $(grep p21_40_0422 in_house_bam.lst)
ln -s $path/Test-20250509/01-output/call-orfs/p21_40_1_0425.raw/output/alignment/p21_40_1_0425.raw_Aligned.sortedByCoord.out.bam p21_40_1_0425_Aligned.sortedByCoord.out.bam
ln -s $path/Test-20250408/output/p21_0321.raw/output/alignment/p21_0321.raw_Aligned.sortedByCoord.out.bam p21_0321_Aligned.sortedByCoord.out.bam
ls $PWD/*_Aligned.sortedByCoord.out.bam > in_house_bam_org.lst
