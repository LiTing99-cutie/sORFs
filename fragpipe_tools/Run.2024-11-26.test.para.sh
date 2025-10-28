#!/usr/bin/sh

################################################
#File Name: Run.2024-11-21.test.para.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2024年11月21日 星期四 19时41分43秒
################################################

set -eo pipefail

manifest=/rd1/user/lit/project/sORFs/config/2024-11-01-2-null-uniprot-custom/fragpipe-files.fp-manifest
for condition in 6_50_prot_1_non_specific_no_quant 6_50_prot_1_non_specific_no_cleavage_no_quant;do
echo -e "***Processing ${condition} at $(date '+%Y-%m-%d %H:%M:%S')"
bash Uni.Fragpipe.sh config/${condition}.workflow \
    $manifest \
    output/MS/Fragpipe_output/${condition} \
    40
done