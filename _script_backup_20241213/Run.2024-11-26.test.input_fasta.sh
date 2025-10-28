#!/usr/bin/sh

################################################
#File Name: Run.2024-11-21.test.para.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2024年11月21日 星期四 19时41分43秒
################################################

set -eo pipefail

manifest=/rd1/user/lit/project/sORFs/config/2024-11-01-2-null-uniprot-custom/fragpipe-files.fp-manifest
database_path=/rd1/user/lit/project/sORFs/custom_database/segmentation/anno_ribo_sorfs/
workflow=/rd1/user/lit/project/sORFs/config/6_50_prot_1_non_specific_no_quant.workflow
for condition in 0_30 0_50 0_100 0_150;do
    mkdir -p output/MS/Fragpipe_output/test_input_fasta_for_less_3k/${condition} && pushd output/MS/Fragpipe_output/test_input_fasta_for_less_3k/${condition}
    if [ "$condition" == "0_30" ]; then
        database=${database_path}/2024-11-26-decoys-contam-anno_ribo_sorfs_30.fa.fas
    elif [ "$condition" == "0_50" ]; then
        database=${database_path}/2024-11-26-decoys-contam-anno_ribo_sorfs_50.fa.fas
    elif [ "$condition" == "0_100" ]; then
        database=${database_path}/2024-11-26-decoys-contam-anno_ribo_sorfs_100.fa.fas
    elif [ "$condition" == "0_150" ]; then
        database=${database_path}/2024-11-26-decoys-contam-anno_ribo_sorfs_150.fa.fas
    else
        echo "Warning: No specific database found for $condition. Skipping."
        continue
    fi
    echo -e "***Processing ${condition} at $(date '+%Y-%m-%d %H:%M:%S')"
    sed 's|database.db-path=.*|database.db-path='"$database"'|' $workflow > fragpipe.workflow
    bash /rd1/user/lit/project/sORFs/Uni.Fragpipe.sh fragpipe.workflow \
        $manifest \
        ./ \
        40
    popd
done