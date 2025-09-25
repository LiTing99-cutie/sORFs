#!/usr/bin/sh

################################################
#File Name: Run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 24 Mar 2025 09:25:20 PM CST
################################################

set -eo pipefail

cd output
less trans_based_sorfs.txt | tail -n +2 | awk -v OFS='\t' '{print $1,$3}' |seqkit tab2fx > nonCano.sorf.trans_based.20250324.fa
bash /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/group_specific_database_20250324/Uni.gen.group_specific_db.20250324.sh \
 $PWD/nonCano.sorf.trans_based.20250324.fa ./group_specific_db