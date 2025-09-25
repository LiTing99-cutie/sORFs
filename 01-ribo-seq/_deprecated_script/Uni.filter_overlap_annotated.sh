#!/usr/bin/sh

################################################
#File Name: Uni.filter_overlap_annotated.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed Dec 11 16:29:49 2024
################################################

set -eo pipefail

# 弃用，直接在R中通过序列去重即可，不需要使用blastp这么复杂的操作

output_path=$1
# 需要绝对路径
query=$2
db_uniprot=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/mouse/uniprotkb_Mus_musculus_reviewed_canonical_and_isoform
db_NCBI=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/NCBI_refseq/mm39/GCF_000001635.27_GRCm39_protein.rmdup

mkdir -p $output_path/tmp && cd $output_path/tmp
# 1. 根据uniprot数据库进行过滤，去掉完全一样的蛋白质
blastp -query $query -db $db_uniprot -out res.out -outfmt '6 qseqid sseqid pident qlen slen length bitscore evalue' -num_threads 20
less res.out | awk '$3==100 && $4==$6' |awk '$5==$6' > exact_same_spep.txt

# 2. 根据NCBI refseq蛋白质数据库进行过滤，去掉完全一样的蛋白质
blastp -query $query -db $db_NCBI -out res.NCBI.out -outfmt '6 qseqid sseqid pident qlen slen length bitscore evalue' -num_threads 20
less res.NCBI.out | awk '$3==100 && $4==$6' |awk '$5==$6' > exact_same_spep.NCBI.txt

# 3. 过滤
## 生成一个需要去除的list
cat exact_same_spep.txt exact_same_spep.NCBI.txt | awk '{print $1}' | sort | uniq > to_be_filtered.txt
## 生成新的tab文件
grep -v -f to_be_filtered.txt ../Ribo_ORFs_merge/nonCano.sorf.tab > ../nonCano.sorf.filtered.tab
## 生成新的fa文件
seqkit tab2fx nonCano.sorf.filtered.tab > ../nonCano.sorf.filtered.fa