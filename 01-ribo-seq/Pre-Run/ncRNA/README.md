# ensembl ncRNA 20240930
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
zless Homo_sapiens.GRCh38.ncrna.fa.gz | grep "^>" | awk '{print $5}' |sort |uniq -c
<!-- 62695 gene_biotype:lncRNA
1945 gene_biotype:miRNA
2419 gene_biotype:misc_RNA
2 gene_biotype:Mt_rRNA
22 gene_biotype:Mt_tRNA
9 gene_biotype:ribozyme
71 gene_biotype:rRNA
51 gene_biotype:scaRNA
1 gene_biotype:scRNA
1020 gene_biotype:snoRNA
2094 gene_biotype:snRNA
6 gene_biotype:sRNA
4 gene_biotype:vault_RNA -->
## 解压，并提取出rRNA,tRNA,snoRNA
seqkit seq -w 0 Homo_sapiens.GRCh38.ncrna.fa > Homo_sapiens.GRCh38.ncrna.singleLine.fa
grep -A1 -f name.list Homo_sapiens.GRCh38.ncrna.singleLine.fa > Homo_sapiens.GRCh38.ncrna.singleLine.target.fa
less Homo_sapiens.GRCh38.ncrna.singleLine.target.fa | grep "^>" | awk '{print $5}' |sort |uniq -c
<!-- 2 gene_biotype:Mt_rRNA
22 gene_biotype:Mt_tRNA
71 gene_biotype:rRNA
1020 gene_biotype:snoRNA -->
# rRNA 
20240930 从SILVA数据库下载人类的sRNA序列
LSU_r138.2.RefNR.zip
SSU_r138.2.RefNR.zip
unzip LSU_r138.2.RefNR.zip 
mv arb-silva.de_2024-09-30_id1350581_tax_silva_trunc.fasta LSU_r138.2.RefNR.fa
unzip SSU_r138.2.RefNR.zip 
mv arb-silva.de_2024-09-30_id1350580_tax_silva_trunc.fasta SSU_r138.2.RefNR.fa
seqkit stat LSU_r138.2.RefNR.fa SSU_r138.2.RefNR.fa
<!-- file                 format  type  num_seqs  sum_len  min_len  avg_len  max_len
LSU_r138.2.RefNR.fa  FASTA   RNA         16   49,137    2,373  3,071.1    3,736
SSU_r138.2.RefNR.fa  FASTA   RNA        364  645,833    1,022  1,774.3    1,909 -->
# tRNA 
20240930 从gtRNAdb下载细胞质tRNA序列
hg38-tRNAs.tar.gz
tar zvxf hg38-tRNAs.tar.gz
seqkit stat hg38-tRNAs.fa
<!-- file           format  type  num_seqs  sum_len  min_len  avg_len  max_len
hg38-tRNAs.fa  FASTA   DNA        432   32,490       70     75.2      108 -->
# 合并所有数据
mkdir merged
cat Homo_sapiens.GRCh38.ncrna.singleLine.target.fa LSU_r138.2.RefNR.fa SSU_r138.2.RefNR.fa hg38-tRNAs.fa > merged/hg38.rRNA.tRNA.snoRNA.fa

cat <(egrep -A1 "gene_biotype:Mt_rRNA|gene_biotype:rRNA" Homo_sapiens.GRCh38.ncrna.singleLine.fa) LSU_r138.2.RefNR.fa SSU_r138.2.RefNR.fa > merged/hg38.rRNA.fa
cat <(egrep -A1 "gene_biotype:Mt_tRNA" Homo_sapiens.GRCh38.ncrna.singleLine.fa) hg38-tRNAs.fa > merged/hg38.tRNA.fa
egrep -A1 "gene_biotype:snoRNA" Homo_sapiens.GRCh38.ncrna.singleLine.fa > merged/hg38.snoRNA.fa
egrep -A1 "gene_biotype:Mt_rRNA|gene_biotype:rRNA" Homo_sapiens.GRCh38.ncrna.singleLine.fa > Homo_sapiens.GRCh38.rRNA.fa