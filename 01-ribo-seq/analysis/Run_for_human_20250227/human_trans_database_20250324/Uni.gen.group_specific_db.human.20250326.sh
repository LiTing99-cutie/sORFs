#!/usr/bin/sh

################################################
#File Name: Uni.gen.group_specific_db.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 24 Mar 2025 03:38:20 PM CST
################################################

set -eo pipefail

uniprot_fasta=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/human/uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.fasta
contam_fasta=/home/user/data3/lit/project/sORFs/02-Mass-spec/custom_database/crap.fasta
user_input_fasta=$1
output_path=$2
mkdir -p $output_path
cd $output_path/
cat $uniprot_fasta $contam_fasta > uniprot.contam.fa

seqkit replace -p 'PE=[0-9]+ ' -r '' uniprot.contam.fa > uniprot.contam.rmPEheader.fa
seqkit replace -p $ -r " PE=1" uniprot.contam.rmPEheader.fa > uniprot.contam.modiHeader.fa
seqkit replace -p $ -r " PE=4" $user_input_fasta > user_input_fasta.modiHeader.fa
cat uniprot.contam.modiHeader.fa user_input_fasta.modiHeader.fa > uniprot.contam.sorfs.trans.modiHeader.fa

