#!/usr/bin/sh

################################################
#File Name: Run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 24 Mar 2025 09:58:24 PM CST
################################################

set -eo pipefail

# 转换candidateORF.prot为输入蛋白质序列模式
candidateORF_fa_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/annotation/RibORF_annot/candidateORF.fa
faTrans -stop $candidateORF_fa_path candidateORF.prot.fa
seqkit fx2tab candidateORF.prot.fa > candidateORF.prot.tab

# 转换uniprot和NCBI为输入模式
uniprot_fasta=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/human/uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.fasta
ncbi_fasta=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/NCBI_refseq/hg38/GCF_000001405.40_GRCh38.p14_protein.faa
seqkit rmdup -s $uniprot_fasta | seqkit seq -s -w0 > /home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/human/uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.rmdup.seq
seqkit rmdup -s $ncbi_fasta | seqkit seq -s -w0 > /home/user/data3/lit/project/sORFs/01-ribo-seq/annot/NCBI_refseq/hg38/GCF_000001405.40_GRCh38.p14_protein.rmdup.seq

# 提取转录本的原始信息
bash Uni.extract.ts.meta.v1.202250325.sh \
    /home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf \
    /home/user/data/lit/database/public/genome/hg38/hg38.fa \
    gencodev41_human_ts_meta_output

# 运行Run.R
Rscript Run.R

# 生成group-specific的数据库
## 修改脚本中的uniprot数据库
cp /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/group_specific_database_20250324/Uni.gen.group_specific_db.20250324.sh ./Uni.gen.group_specific_db.human.20250326.sh
cd output
less trans_based_sorfs.txt | tail -n +2 | awk -v OFS='\t' '{print $1,$3}' |seqkit tab2fx > nonCano.sorf.trans_based.fa
bash ../Uni.gen.group_specific_db.human.20250326.sh \
 $PWD/nonCano.sorf.trans_based.fa ./group_specific_db

# Run.R中有错误，在Galaxy上生成正确的Run.fix.error.20240410.R
## 更新了all_orfs_ORF_id_seq_trans_map.rds，all_orfs_top_trans.rds，以及trans_based_sorfs.txt
cd /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output
path_1=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250404_tmp_from_nanjing/output_20250407/
cp $path_1/all_orfs_top_trans.rds ./
cp $path_1/trans_based_sorfs.txt ./
cp $path_1/correct_incorrect_map.rds /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output