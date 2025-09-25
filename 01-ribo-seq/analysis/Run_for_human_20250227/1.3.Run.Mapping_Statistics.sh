#!/usr/bin/sh

################################################
#File Name: 1.3.Run.Mapping_Statistics.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 27 Mar 2025 04:13:27 PM CST
################################################

set -eo pipefail

##### 1. 质控步骤过滤掉的reads百分比 #####
Rscript 1.3.Mapping_Statistics.R

##### 2. 唯一比对数量以及唯一比对率 #####
# 定义提取信息的函数
extract_info() {
    local search_path=$1
    local search_term=$2
    local output_file=$3
    local label=$4

    find "$search_path" -name "*Log.final.out" | xargs grep "$search_term" |\
    awk -F'[:\t]' -v label="$label" '{split($1, a, "/"); split(a[length(a)], b, "_Log.final.out"); print b[1] "\t" label "\t" $NF}' > "$output_file"
}

# 处理 ribo-seq 数据
extract_info "$PWD/human_brain_output_20250227/" "Uniquely mapped reads %" \
    "human_brain_output_20250227/stat/Uniquely_mapped_reads_rate.ribo.txt" "Uniquely mapped reads %"

extract_info "$PWD/human_brain_output_20250227/" "Uniquely mapped reads number" \
    "human_brain_output_20250227/stat/Uniquely_mapped_reads_number.ribo.txt" "Uniquely mapped reads number"

# 处理 RNA-seq 数据
extract_info "$PWD/human_brain_rna_seq_alignment_assemble/" "Uniquely mapped reads %" \
    "human_brain_output_20250227/stat/Uniquely_mapped_reads_rate.rna.txt" "Uniquely mapped reads %"

extract_info "$PWD/human_brain_rna_seq_alignment_assemble/" "Uniquely mapped reads number" \
    "human_brain_output_20250227/stat/Uniquely_mapped_reads_number.rna.txt" "Uniquely mapped reads number"

merge_files_by_sample(){
	# 定义输入和输出文件路径
	rate_file=$1
	number_file=$2
	output_file=$3

	# 先提取第一列和第三列，再合并
	awk -F'\t' '{print $1 "\t" $3}' "$rate_file" > temp_rate.txt
	awk -F'\t' '{print $1 "\t" $3}' "$number_file" > temp_number.txt

	# 按第一列（样本名）合并，并添加表头
	echo -e "Sample\tuniquely_mapped_reads_rate\tUniquely_mapped_reads_number" > "$output_file"
	join -t $'\t' -1 1 -2 1 temp_rate.txt temp_number.txt >> "$output_file"

	# 清理临时文件
	rm -rf temp_rate.txt temp_number.txt
}

merge_files_by_sample "human_brain_output_20250227/stat/Uniquely_mapped_reads_rate.ribo.txt" "human_brain_output_20250227/stat/Uniquely_mapped_reads_number.ribo.txt" \
	"human_brain_output_20250227/stat/Merged_Uniquely_mapped_reads.ribo.txt"
merge_files_by_sample "human_brain_output_20250227/stat/Uniquely_mapped_reads_rate.rna.txt" "human_brain_output_20250227/stat/Uniquely_mapped_reads_number.rna.txt" \
"human_brain_output_20250227/stat/Merged_Uniquely_mapped_reads.rna.txt"

###### 3. 读段的长度分布 ######
nohup bash 1.3c.dataset.ribo.Rlenth.sh \
	./human_brain_output_20250227/Aligned.sortedByCoord.out.bam.lst \
	./human_brain_output_20250227/stat/reads_length_distribution.txt &
###### 4. 读段在基因组上的分布 ######
# see 1.3b.dataset.ribo.distri.sh