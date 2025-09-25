#!/usr/bin/sh

################################################
#File Name: Run.cal.rna_reads.20250213.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 13 Feb 2025 03:17:47 PM CST
################################################

set -eo pipefail

# gtf=test/test.1.gtf
gtf=$1
output_path=$2
if [ "$1" == "-h" ]; then
    echo "Usage: $0 <gtf> <output_path>"
    exit 0
fi
bam_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/te_calc_20250212/output/mapping_rna_seq_filtered
[ -f $bam_path/bam.lst ] || ls $bam_path/*Aligned.sortedByCoord.out.bam > $bam_path/bam.lst
mkdir -p $output_path

# 得到指定ORF的RNA-seq reads数目
bam_lst_1=$(egrep -v  "SRR5262869|SRR5262868|SRR11218257|SRR11218256" $bam_path/bam.lst)
bam_lst_2=$(egrep "SRR5262869|SRR5262868" $bam_path/bam.lst)
bam_lst_3=$(egrep "SRR11218257|SRR11218256" $bam_path/bam.lst)

# 【注意】featureCounts的链特异性参数-s需要根据特定的文库来修改
featureCounts -s 1 -t CDS -g transcript_id -a $gtf -o $output_path/rna_counts_1.txt $bam_lst_1
featureCounts -t CDS -g transcript_id -a $gtf -o $output_path/rna_counts_2.txt $bam_lst_2
featureCounts -s 2 -t CDS -g transcript_id -a $gtf -o $output_path/rna_counts_3.txt $bam_lst_3