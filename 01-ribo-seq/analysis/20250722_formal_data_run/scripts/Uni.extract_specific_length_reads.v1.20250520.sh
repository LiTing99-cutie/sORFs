#!/bin/bash

# 功能：高效提取指定长度reads（支持默认25-34nt范围）
# 用法：./extract_reads.sh <input.bam> [lengths.txt] <output.bam>

# 检查输入参数
if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 <input.bam> [lengths.txt] <output.bam>"
    echo "Note: If lengths.txt is omitted, will extract 25-34nt reads by default"
    exit 1
fi

# 参数解析（自动判断参数格式）
if [ $# -eq 2 ]; then
    # 两个参数模式：默认25-34nt
    BAM_FILE=$1
    OUTPUT_BAM=$2
    DEFAULT_RANGE=1
else
    # 三个参数模式：使用自定义长度文件
    BAM_FILE=$1
    LENGTH_FILE=$2
    OUTPUT_BAM=$3
    [ -f "$LENGTH_FILE" ] || { echo "Error: Lengths file '$LENGTH_FILE' not found!"; exit 1; }
fi

[ -f "$BAM_FILE" ] || { echo "Error: BAM file '$BAM_FILE' not found!"; exit 1; }

# 步骤1：设置长度过滤条件
if [ "$DEFAULT_RANGE" = 1 ]; then
    # 默认25-34nt范围
    LENGTH_COND='length($10)>=25 && length($10)<=34'
    echo "Using DEFAULT range 25-34nt"
else
    # 从文件读取目标长度
    declare -A TARGET_LENGTHS
    while read -r len; do
        TARGET_LENGTHS[$len]=1
    done < <(awk '{print $1}' "$LENGTH_FILE" | sort -u)
    
    # 生成AWK条件
    LENGTH_COND="length(\$10) in target"
    AWK_VARS="-v lengths_str=\"${!TARGET_LENGTHS[*]}\""
    echo "Using CUSTOM lengths from $LENGTH_FILE"
fi

# 步骤2：单次处理输出目标reads
samtools view -h "$BAM_FILE" | \
awk $AWK_VARS 'BEGIN {
    if (lengths_str != "") {
        split(lengths_str, arr, " ");
        for (i in arr) target[arr[i]]=1
    }
}
NR==1 || /^@/ {print; next}  # 保留头信息
'"$LENGTH_COND" | \
samtools view -b -o "$OUTPUT_BAM"

echo "Success: Extracted reads to $OUTPUT_BAM"