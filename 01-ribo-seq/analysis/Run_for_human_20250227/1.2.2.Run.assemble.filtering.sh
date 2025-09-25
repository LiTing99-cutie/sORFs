#!/usr/bin/sh

################################################
#File Name: 1.2.2.Run.assemble.filtering.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 03 Mar 2025 02:43:45 PM CST
################################################

set -eo pipefail

ref_gtf=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
gffcompare -r $ref_gtf -o gffcmp assembled.gtf

# 24025
less gffcmp.annotated.gtf |grep -v gene_name |awk '$3=="transcript"' | \
awk -F'\t' '{if ($3 == "transcript") {match($9, /class_code "([^"]+)"/, arr); if (arr[1] != "") print arr[1]}}' | sort | uniq -c

less gffcmp.annotated.gtf |awk '$3=="transcript"' | \
awk -F'\t' '{if ($3 == "transcript") {match($9, /class_code "([^"]+)"/, arr); if (arr[1] != "") print arr[1]}}' | sort | uniq -c
# 142373 =
#   10593 c
#    1347 e
#   18486 i
#   30801 j
#    1218 k
#     934 m
#    3296 n
#    2782 o
#    2385 p
#   24025 u
#      39 x
#       2 y

# 142070
less assembled.gtf |grep reference_id|grep -o -P 'transcript_id "[^"]*"' | cut -d'"' -f2 | sort -u | wc -l

# 过滤FPKM大于等于1的类型
awk '$3=="transcript" && match($0, /FPKM "([0-9.]+)"/, a) && a[1] >= 1' assembled.gtf |wc -l

# 提取新转录本
awk '$3=="transcript" && /class_code "u"/' gffcmp.annotated.gtf

# 96211
less assembled.gtf |grep -v reference_id|awk '$3=="transcript"'|wc -l
# 68434
less assembled.gtf |grep -v reference_id|grep -o 'gene_id "[^"]*"' | cut -d'"' -f2 | sort -u | wc -l
# 41110
less assembled.gtf |grep reference_id|grep -o -P 'ref_gene_id "[^"]*"' | cut -d'"' -f2 | sort -u | wc -l
less assembled.gtf |grep reference_id|grep -o -P 'transcript_id "[^"]*"' | cut -d'"' -f2 | sort -u | wc -l