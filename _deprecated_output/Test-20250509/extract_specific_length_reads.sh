#!/bin/bash

# 功能：从BAM文件中提取指定长度的reads，临时文件保存在当前目录
# 用法：./extract_specific_length_reads.sh <input.bam> <lengths.txt> <output.bam>

# 检查输入参数
if [ $# -ne 3 ]; then
    echo "Usage: $0 <input.bam> <lengths.txt> <output.bam>"
    exit 1
fi

BAM_FILE=$1
LENGTH_FILE=$2
OUTPUT_BAM=$3

# 检查文件是否存在
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file '$BAM_FILE' not found!"
    exit 1
fi

if [ ! -f "$LENGTH_FILE" ]; then
    echo "Error: Lengths file '$LENGTH_FILE' not found!"
    exit 1
fi

# 步骤1：提取目标长度列表
TARGET_LENGTHS=$(awk '{print $1}' "$LENGTH_FILE" | sort -u | tr '\n' ' ')

# 步骤2：在当前目录创建临时文件（而非/tmp）
TEMP_READ_NAMES="$(pwd)/tmp_read_names.$$"  # 使用$$添加进程ID避免冲突
touch "$TEMP_READ_NAMES"

# 步骤3：提取目标reads的名称
samtools view "$BAM_FILE" | \
awk -v lengths="$TARGET_LENGTHS" 'BEGIN {
    split(lengths, arr, " ");
    for (i in arr) target[arr[i]]=1
} 
length($10) in target {print $1}' | \
sort -T . -u > "$TEMP_READ_NAMES"

# 步骤4：生成最终BAM
samtools view -h "$BAM_FILE" | \
awk -v read_names="$TEMP_READ_NAMES" 'BEGIN {
    while ((getline < read_names) > 0) names[$1]=1;
    close(read_names)
}
NR==1 || $0 ~ /^@/ || ($1 in names)' | \
samtools view -b -o "$OUTPUT_BAM"

# 清理临时文件
rm "$TEMP_READ_NAMES"
echo "Success: Extracted reads to $OUTPUT_BAM (temp file: $TEMP_READ_NAMES deleted)"