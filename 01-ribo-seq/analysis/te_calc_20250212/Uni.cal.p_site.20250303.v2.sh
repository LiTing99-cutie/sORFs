#!/usr/bin/sh

################################################
#File Name: Run.cal.p_site.20250213.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 13 Feb 2025 11:02:06 AM CST
################################################

set -eo pipefail

# gtf=test/test.1.gtf
# 需要含有CDS这一feature，否则需要修改代码
gtf=$1
output_path=$2
# 是否把overlapped的reads计算在内,yes or no
overlapped=$3
mkdir -p $output_path

# 得到指定ORF的p-site数目或者RPF reads数目
bam_lst_1=$(find /home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/ -type f -path '*/output/alignment/*_Aligned.sortedByCoord.out.bam')
bam_lst_2=$(find /home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/ -type f -path '*/output/Ribo_ORFs_add_assemble_20250125/RibORF/corrected*sam')

# 【注意】featureCounts的链特异性参数-s需要根据特定的文库来修改
if [[ $overlapped == "yes" ]];then
	featureCounts -O -s 1 -t CDS -g transcript_id -a $gtf -o $output_path/RPF.counts.txt $bam_lst_1
	featureCounts -O -s 1 -t CDS -g transcript_id -a $gtf -o $output_path/p_site.counts.txt $bam_lst_2
else
	featureCounts -s 1 -t CDS -g transcript_id -a $gtf -o $output_path/RPF.counts.txt $bam_lst_1
	featureCounts -s 1 -t CDS -g transcript_id -a $gtf -o $output_path/p_site.counts.txt $bam_lst_2
fi