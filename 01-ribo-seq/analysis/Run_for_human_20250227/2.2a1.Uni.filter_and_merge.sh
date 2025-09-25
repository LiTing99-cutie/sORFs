#!/usr/bin/sh

################################################
#File Name: 2.2.1.Uni.filter_and_merge.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 28 Mar 2025 10:12:56 PM CST
################################################

set -eo pipefail

# 输入文件路径
bam_paths_file=$1       # 存储 BAM 文件路径的文件
length_sample_file=$2  # 长度和样本的文件
output_bam=$3  # 输出的合并 BAM 文件
temp_dir=$4 # 临时存储路径例如 "./temp_bam"
mkdir -p "$temp_dir"

# 读取 BAM 文件路径
declare -A bam_paths
while read -r line; do
    if [[ -n "$line" ]]; then
        sample_name=$(basename $line|sed 's/_Aligned.toTranscriptome.out.bam//;s/_Aligned.sortedByCoord.out.bam//')  # 提取样本名
        bam_paths["$sample_name"]="$line"
    fi
done < "$bam_paths_file"

# 提取指定长度的 reads 并存储临时文件
while IFS=',' read -r length _ sample; do
    if [[ -n "${bam_paths[$sample]}" ]]; then
        bam_file="${bam_paths[$sample]}"
        temp_bam="$temp_dir/${sample}_${length}.bam"

        echo "Processing sample: $sample, length: $length"
        # 使用 samtools 过滤指定长度的 reads
        samtools view -h "$bam_file" | \
        awk -v len="$length" 'BEGIN {OFS="\t"} $1 ~ /^@/ || length($10) == len' | \
        samtools view -b -o "$temp_bam"
    else
        echo "Sample $sample not found in BAM paths file"
    fi
done < "$length_sample_file"

# 合并所有临时 BAM 文件
echo "Merging all filtered BAM files..."
samtools merge -f "$output_bam" "$temp_dir"/*.bam

# 清理临时文件
# rm -rf "$temp_dir"

echo "Filtered and merged BAM file saved to: $output_bam"
