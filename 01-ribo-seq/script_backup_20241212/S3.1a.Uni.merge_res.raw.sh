#!/usr/bin/sh

################################################
#File Name: Uni.merge_res.raw.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed Nov  6 10:21:10 2024
################################################

set -eo pipefail

# 在每一个Ribo-ORFs文件夹下运行
cd merge
sample=$(echo $PWD | awk -F'/' '{print $(NF-3)}')
echo $sample

# PRICE
awk -v OFS='\t' '{print $1,$2,"PRICE"}' PRICE/nonCano.sorf.tab > nonCano.sorf.meta.raw.PRICE.txt

# RiboCode
awk -v OFS='\t' '{print $1,$2,"RiboCode"}' RiboCode/nonCano.sorf.tab > nonCano.sorf.meta.raw.RiboCode.txt

# RibORF
awk -v OFS='\t' '{print $1,$2,"RibORF"}' RibORF/nonCano.sorf.tab > nonCano.sorf.meta.raw.RibORF.txt

# merge
cat nonCano.sorf.meta.raw.PRICE.txt nonCano.sorf.meta.raw.RiboCode.txt nonCano.sorf.meta.raw.RibORF.txt | \
awk -v OFS='\t' '{print $1,$2,$3,"'$sample'"}' > nonCano.sorf.meta.merge.raw.3_ways.txt
