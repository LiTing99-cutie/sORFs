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
bam_lst=$4
mkdir -p $output_path

# 得到指定ORF的p-site数目或者RPF reads数目

# 【注意】featureCounts的链特异性参数-s需要根据特定的文库来修改
if [[ $overlapped == "yes" ]];then
	featureCounts -O -s 1 -t CDS -g transcript_id -a $gtf -o $output_path/RPF.counts.txt $(cat $bam_lst)
else
	featureCounts -s 1 -t CDS -g transcript_id -a $gtf -o $output_path/RPF.counts.txt $(cat $bam_lst)
fi