#!/usr/bin/sh

################################################
#File Name: 1.3c.dataset.ribo.Rlenth.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 28 Mar 2025 10:52:02 PM CST
################################################

set -eo pipefail

# 定义输入文件和输出文件
bam_paths_file=$1  # 存储 BAM 文件路径的文件
output_file=$2  # 输出文件，例如"reads_length_distribution.txt"

# 检查输入文件是否存在
if [[ ! -f "$bam_paths_file" ]]; then
    echo "Error: BAM paths file '$bam_paths_file' not found!"
    exit 1
fi

# 清空输出文件
: > "$output_file"

# 遍历 BAM 文件路径
while read -r bam_path; do
    if [[ -n "$bam_path" ]]; then
        # 提取样本名（假设样本名为文件名前缀）
        sample_name=$(basename  -s "_Aligned.toTranscriptome.out.bam" -s "_Aligned.sortedByCoord.out.bam" "$bam_path")
        
        echo "Processing sample: $sample_name"
        
        # 使用 samtools 统计每种长度的 reads 数目
        samtools view -@ 30 "$bam_path" | \
        awk '{print length($10)}' | \
        sort | uniq -c | \
        awk -v sample="$sample_name" '{print sample, $2, $1}' >> "$output_file"
    fi
done < "$bam_paths_file"

echo "Reads length distribution has been saved to: $output_file"
