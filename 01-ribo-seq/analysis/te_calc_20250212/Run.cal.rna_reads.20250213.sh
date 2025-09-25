#!/usr/bin/sh

################################################
#File Name: Run.cal.rna_reads.20250213.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 13 Feb 2025 03:17:47 PM CST
################################################

set -eo pipefail

# gtf=test/test.1.gtf
gtf=/home/user/data3/lit/project/sORFs/03-Cross-anna/gen_genepred/sorfs.all.gtf
output_path=output/sorfs_all/rna_counts
bam_path=output/mapping_rna_seq/
mkdir -p $output_path

# 得到指定ORF的p-site数目或者RPF reads数目
bam_lst_1=$(find $bam_path -name "*Aligned.sortedByCoord.out.bam"| grep -f <(less $bam_path/*lst |egrep "2015_Science|2019_multispecies"|xargs basename -s ".fastq.gz"))
bam_lst_2=$(find $bam_path -name "*Aligned.sortedByCoord.out.bam"|egrep "SRR5262869|SRR5262868")
bam_lst_3=$(find $bam_path -name "*Aligned.sortedByCoord.out.bam"|egrep "SRR11218257|SRR11218256")

# 【注意】featureCounts的链特异性参数-s需要根据特定的文库来修改
featureCounts -s 1 -t CDS -g transcript_id -a $gtf -o $output_path/rna_counts_1.txt $bam_lst_1
featureCounts -t CDS -g transcript_id -a $gtf -o $output_path/rna_counts_2.txt $bam_lst_2
featureCounts -s 2 -t CDS -g transcript_id -a $gtf -o $output_path/rna_counts_3.txt $bam_lst_3