#!/usr/bin/sh

################################################
#File Name: Uni.qc.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 17 Feb 2025 08:34:04 PM CST
################################################

set -eo pipefail
# run qc in batch
output_path=$1
bam_lst=$2
fastqc=$3
rm_adpt=$4
trimmed_fastqc=$5
export script=Uni.qc.single.sh
export output_path=$output_path
export fastqc=$fastqc
export rm_adpt=$rm_adpt
export trimmed_fastqc=$trimmed_fastqc
mkdir -p $output_path/log
cat $bam_lst| \
parallel -j 5 --joblog $output_path/log/qc.20250217.log 'bash $script {} none $output_path $fastqc $rm_adpt $trimmed_fastqc'
if [ $fastqc = "yes" ]; then
	pushd $output_path && mkdir -p multiqc_before && find ./ -path "*/fastqc/*" |xargs -i cp {} ./multiqc_before && multiqc -n ./multiqc_before.html ./multiqc_before && popd
fi
if [ $trimmed_fastqc = "yes" ]; then
	pushd $output_path && mkdir -p multiqc_after && find ./ -path "*/trimmed_fastqc/*" |xargs -i cp {} ./multiqc_after && multiqc -n ./multiqc_after.html multiqc_after && popd
fi