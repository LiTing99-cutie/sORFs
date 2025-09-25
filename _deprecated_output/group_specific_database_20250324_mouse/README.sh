#!/usr/bin/sh

################################################
#File Name: README.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 24 Mar 2025 03:41:10 PM CST
################################################

set -eo pipefail
user_input_fasta=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/trans_based_database/output/nonCano.sorf.20250206.trans_based.fa
bash Uni.gen.group_specific_db.20250324.sh $user_input_fasta trans_based_group_specific_20250324
