#!/usr/bin/sh

################################################
#File Name: Run.20250201.gen.trans.database.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Sat 01 Feb 2025 02:28:35 PM CST
################################################

set -eo pipefail

mkdir -p output

### Step 0 ###
bash S5.Gen_trans_database.sh

### Step 1 折叠相同终止密码子对应的小肽 ###
Rscript S5.Transcripome_database_collapse_stop_codon.R \
	/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/mouse_run_20250125/annotation/RibORF_annot/candidateORF.genepred.txt \
	output/candidateORF.genepred.callapse.stop_codon.txt

### Step 2 过滤，整理，选择代表性转录本 ###
Rscript S5.Gen_trans_database_data.table.R

uniprot_fasta=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/mouse/uniprotkb_Mus_musculus_reviewed_canonical_and_isoform_2024_10_25.fasta
cd output
less trans_based_sorfs.txt | tail -n +2 | awk -v OFS='\t' '{print $1,$3}' |seqkit tab2fx > nonCano.sorf.20250206.trans_based.fa
cat $uniprot_fasta nonCano.sorf.20250206.trans_based.fa > uniprot.nonCano.sorf.20250206.trans_based.fa

bash $translate_script $prefix.gtf $fa
seqkit fx2tab cds.fa > $prefix.cds.tab
seqkit fx2tab prot.fa > $prefix.prot.tab