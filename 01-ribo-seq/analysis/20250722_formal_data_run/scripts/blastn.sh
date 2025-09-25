#!/bin/bash

# 创建查询序列文件
mkdir ../tmp/blastn
cd ../tmp/blastn
echo ">query_seq" > query.fa
echo "GGTTTCGTACGTAGCAGAGCA" >> query.fa
echo ">query_seq" > query.1.fa
echo "GTCTAGGCCACACCACCCTGAAGGCGCCTGCTCGCCTCTGATCTGTTGAAGCTAAGCAGGGTCGGTCCTGGTTAGTACTTGGATGGGACTCCGCCTGGTAATAGCCGGTACCGTAGGCT" >> query.1.fa


# 构建BLAST数据库（只需一次，若已构建可跳过）
makeblastdb -in /home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/ncRNA/merged/hg38.rRNA.fa -dbtype nucl -out rRNA_db

# 运行BLAST比对
blastn -query query.fa \
       -db rRNA_db \
       -out blast_result.txt \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -evalue 1
# 打印结果
echo "比对完成，结果如下："
cat blast_result.txt

sed '/^>/!s/U/T/g' /home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/ncRNA/merged/hg38.rRNA.fa > hg38.rRNA.fa
makeblastdb -in hg38.rRNA.fa -dbtype nucl -out rRNA_db_new
blastn -query query.fa \
       -db rRNA_db_new \
       -out blast_result.txt \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -evalue 1

seqkit seq -w 0 hg38.rRNA.fa  > hg38.rRNA.1.fa

cat query.fa | sed 's/T/U/g' > query.T2U.fa

blastn -query query.T2U.fa \
       -subject  /home/user/data3/lit/project/sORFs/01-ribo-seq/ncRNA/LSU_r138.2.RefNR.fa \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
blastn -query query.T2U.fa \
       -subject  /home/user/data3/lit/project/sORFs/01-ribo-seq/ncRNA/SSU_r138.2.RefNR.fa \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

mkdir rRNA_STARindex
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir rRNA_STARindex --genomeFastaFiles ./hg38.rRNA.fa --genomeSAindexNbases 8
mkdir RNAcentral_rRNA_STARindex
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir RNAcentral_rRNA_STARindex --genomeFastaFiles ./RNAcentral.human.rRNA.20250725.fasta --genomeSAindexNbases 8

cat hg38.rRNA.1.fa RNAcentral.human.rRNA.20250725.fasta > rRNA.total.fa

mkdir total_rRNA_STARindex
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir total_rRNA_STARindex --genomeFastaFiles ./rRNA.total.fa --genomeSAindexNbases 8
