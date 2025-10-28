#!/usr/bin/sh

################################################
#File Name: Run.2024-11-21.test.para.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2024年11月21日 星期四 19时41分43秒
################################################

set -eo pipefail

manifest=/rd1/user/lit/project/sORFs/config/2024-11-01-3-trypsin-lysc-uniprot-custom/fragpipe-files.fp-manifest
for condition in trypsin_lysc_6_50_prot_1_no_MBR trypsin_lysc_6_50_prot_1_no_MBR_missed_cleavage_3 trypsin_lysc_6_50_prot_1_no_MBR_cysteine_variable;do
echo -e "***Processing ${condition} at $(date '+%Y-%m-%d %H:%M:%S')"
bash Uni.Fragpipe.sh config/workflow_LFQ_MBR_${condition}.workflow \
    $manifest \
    output/MS/Fragpipe_output/${condition} \
    40
done

manifest=/rd1/user/lit/project/sORFs/config/2024-11-01-2-null-uniprot-custom/fragpipe-files.fp-manifest
for condition in 6_50_prot_1_no_cleavage;do
echo -e "***Processing ${condition} at $(date '+%Y-%m-%d %H:%M:%S')"
bash Uni.Fragpipe.sh config/workflow_LFQ_MBR_${condition}.workflow \
    $manifest \
    output/MS/Fragpipe_output/${condition} \
    40
done