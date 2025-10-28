#!/bin/bash

# 输入文件
FILE_LIST="$PWD/file_list.txt"       # 文件路径列表
SAMPLE_MAP="$PWD/sample_mapping.txt" # 样本名映射表

# 读取映射表到关联数组
declare -A SAMPLE_MAPPING
while read -r num sample; do
    SAMPLE_MAPPING["$num"]="$sample"
done < "$SAMPLE_MAP"

# 处理每个文件
while read -r filepath; do
    if [[ -d "$filepath" ]]; then
        # 提取原文件名部分
        filename=$(basename "$filepath")

        # 提取数字（如99）
        num=$(echo "$filename" | grep -oP 'min_\K\d+(?=_Slot)')

        # 获取对应的样本名
        sample="${SAMPLE_MAPPING[$num]}"

        if [[ -n "$sample" ]]; then
            # 构造新文件名（替换数字为样本名）
            new_filename=$(echo "$filename" | sed "s/min_${num}_/min_${sample}_/")
            new_path="$PWD/$new_filename"

            # 创建软链接
            ln -sfn "$filepath" "$new_path"
            echo "Created symlink: $new_path -> $filepath"
        else
            echo "No mapping found for number: $num (file: $filename)"
        fi
    fi
done < "$FILE_LIST"