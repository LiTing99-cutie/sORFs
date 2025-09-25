#!/usr/bin/sh

################################################
#File Name: README.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 13 Feb 2025 09:25:16 PM CST
################################################

set -eo pipefail

mkdir test && cd test
# 有些ncORF有可能基因本身是一个protein_coding基因，需要进行intersect来确保不和任何的CDS的外显子区域overlap

## 小肽的bed6
ms_sorfs_gtf=/home/user/data3/lit/project/sORFs/03-Cross-anna/gen_genepred/essen_info_to_gen_genepred.sorfs.gtf
less $ms_sorfs_gtf |awk '$3=="CDS"'|\
awk -v OFS='\t' '{match($0,/gene_id "([^"]+)"/,arr);print $1,$4-1,$5,arr[1],$6,$7}' > ms_sorfs.bed6

## 经典ORF的bed6
anno_gtf=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.gtf
less $anno_gtf |awk '$3=="CDS"'|\
awk -v OFS='\t' '{match($0,/gene_id "([^"]+)"/,arr);print $1,$4-1,$5,arr[1],$6,$7}' > anno.bed6

## intersect
bedtools intersect -a ms_sorfs.bed6 -b anno.bed6 -s | sort -k1,1 -k2,2n > intersect.sorted.txt

mkdir tmp
for orf_id in $(less intersect.sorted.txt |cut -f 4|sort|uniq);do
	grep $orf_id intersect.sorted.txt > tmp/$orf_id.txt
	bedtools merge -i tmp/$orf_id.txt -c 4,5,6 -o distinct > tmp/$orf_id.dedup.txt
done
cat tmp/*.dedup.txt > intersect.dedup.grouped.txt
awk -v OFS='\t' '{sum[$4] += $3 - $2} END {for (key in sum) print key, sum[key]}' intersect.dedup.grouped.txt > ms_sorfs.overlapped.nt.txt.1

## 检查某些case的三碱基周期性
mkdir -p case_tree_nt_period
bam=/home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/2015_Science_ChoJ_hippocampus/SRR2163103/output/Ribo_ORFs_add_assemble_20250125/RibORF/corrected.SRR2163103.bam
grep ENSMUST00000119612.9 /home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.gtf | awk '$3=="CDS"'|\
awk -v OFS='\t' '{match($0,/gene_id "([^"]+)"/,arr);print $1,$4-1,$5,arr[1],$6,$7}' > test.bed6

bedtools coverage -a test.bed6 -b $bam -d > test.bedgraph

sorf_id=ENSMUST00000184255.2-chr17:8785230-8785647
grep $sorf_id ms_sorfs.bed6 > case_tree_nt_period/$sorf_id.bed6
# 大概需要3-4min
bedtools coverage -a case_tree_nt_period/$sorf_id.bed6 -b $bam -d -s > case_tree_nt_period/$sorf_id.bedgraph

# 测试下用bed会不会更加快
bedtools bamtobed -i $bam > output.bed

bedtools intersect -a test.bed6 -b output.bed -wao -s