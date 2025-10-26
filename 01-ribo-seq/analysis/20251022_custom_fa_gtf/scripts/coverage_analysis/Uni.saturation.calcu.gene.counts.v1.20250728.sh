#!/bin/bash

source /home/user/data2/lit/bin/lit_utils.sh
gtf=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/results/custom.gtf.with_orf.gtf
bam_dir=$(realpath "$1")
output_dir=$(realpath "$2")
mkdir -p $output_dir/libtype
find $bam_dir -name "*sort*bam" > $output_dir/bam.lst
# bash 2.Uni.libType.v1.20250529.sh $output_dir/bam.lst $output_dir/libtype
# 找到所有bam文件并运行featureCounts 
featureCounts -a "$gtf" -s 1 -o "$output_dir/transcript.counts.txt" -T 30 -O -g transcript_id -t exon $(cat $output_dir/bam.lst)
featureCounts -a "$gtf" -s 1 -o "$output_dir/gene.counts.txt" -T 30 -O -g gene_id -t exon $(cat $output_dir/bam.lst)

