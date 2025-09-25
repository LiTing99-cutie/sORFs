#!/usr/bin/sh

################################################
#File Name: Run.filter.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue Nov  5 10:12:16 2024
################################################

set -eo pipefail

mkdir -p mouse_brain_output_20241011/filter
output_path=mouse_brain_output_20241011/filter
query=/home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/Ribo_ORFs_merge/nonCano.sorf.fa
db=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/mouse/uniprotkb_Mus_musculus_reviewed_canonical_and_isoform

cd $output_path
# 首先根据uniprot数据库进行过滤，去掉完全一样的蛋白质
blastp -query $query -db $db -out res.out -outfmt '6 qseqid sseqid pident qlen slen length bitscore evalue' -num_threads 20
less res.out | awk '$3==100 && $4==$6' |awk '$5==$6' > exact_same_spep.txt

# 然后根据NCBI refseq蛋白质数据库进行过滤，去掉完全一样的蛋白质
# makeblastdb -in /home/user/data3/lit/project/sORFs/01-ribo-seq/annot/NCBI_refseq/mm39/GCF_000001635.27_GRCm39_protein.rmdup.faa -dbtype prot -out /home/user/data3/lit/project/sORFs/01-ribo-seq/annot/NCBI_refseq/mm39/GCF_000001635.27_GRCm39_protein.rmdup
db_NCBI=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/NCBI_refseq/mm39/GCF_000001635.27_GRCm39_protein.rmdup
blastp -query $query -db $db_NCBI -out res.NCBI.out -outfmt '6 qseqid sseqid pident qlen slen length bitscore evalue' -num_threads 20
less res.NCBI.out | awk '$3==100 && $4==$6' |awk '$5==$6' > exact_same_spep.NCBI.txt

# 生成一个需要去除的list (共462个)
cat exact_same_spep.txt exact_same_spep.NCBI.txt | awk '{print $1}' | sort | uniq > to_be_filtered.txt

# 生成新的tab文件
grep -v -f to_be_filtered.txt ../Ribo_ORFs_merge/nonCano.sorf.tab > nonCano.sorf.filtered.tab

# 生成新的fa文件
seqkit tab2fx nonCano.sorf.filtered.tab > nonCano.sorf.filtered.fa