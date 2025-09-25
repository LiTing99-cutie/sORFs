#!/usr/bin/sh

################################################
#File Name: Uni.extract.ts.meta.v1.202250325.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 25 Mar 2025 07:52:37 PM CST
################################################

set -eo pipefail

gtf=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
fa=/home/user/data/lit/database/public/genome/hg38/hg38.fa
# gencodev41_human_ts_meta_output
output_path=$3

# 从gtf中提取转录本相关的元数据
## 计算CDS长度
mkdir -p $output_path && cd $output_path
bash /home/user/data3/lit/project/sORFs/01-ribo-seq/S3.0c.Uni.translate_gtf.v2.20250325.sh \
	$gtf \
	$fa
seqkit fx2tab -nl prot.fa > TsId_CDS_length.txt
## 提取其他的元数据
awk -F '\t' '$3 == "transcript" {
    match($9, /transcript_id "([^"]+)"/, tid)
    match($9, /gene_id "([^"]+)"/, gid)
    match($9, /gene_name "([^"]+)"/, gname)
    match($9, /transcript_support_level "([^"]+)"/, tsl)
    match($9, /tag "(appris[_a-z0-9]+)"/, appris)

    gid_val = gid[1] ? gid[1] : ""
    gname_val = gname[1] ? gname[1] : ""
    tsl_val = tsl[1] ? tsl[1] : "NA"
    appris_val = appris[1] ? appris[1] : "None"

    print tid[1]"\t"gid_val"\t"gname_val"\t""tsl"tsl_val"\t"appris_val
}' $gtf |sed 's/appris_//g' > transcript_metadata.txt

awk 'NR==FNR{len[$1]=$2; next} {print $0"\t"(len[$1] ? len[$1] : 0)}' TsId_CDS_length.txt transcript_metadata.txt > ts.meta.txt
