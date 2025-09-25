#!/usr/bin/sh

################################################
#File Name: 2.2.Uni.fc.20250530.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 30 May 2025 11:10:54 AM CST
################################################

set -eo pipefail

#!/bin/bash

# 脚本使用方法说明
usage() {
    echo "Usage: $0 -i input_file -g gtf_file -o output_prefix [-O] [-T threads]"
    echo "  -i 输入文件，包含BamPath、LibraryType和SingleEnd信息"
    echo "  -g GTF注释文件"
    echo "  -o 输出文件前缀"
    echo "  -O 可选参数，计算重叠的reads数量"
    echo "  -T 线程数，默认为10"
    exit 1
}

# 初始化变量
input_file=""
gtf_file=""
output_prefix=""
count_overlaps=0
threads=10

# 解析命令行参数
while getopts ":i:g:o:OT:" opt; do
    case $opt in
        i) input_file="$OPTARG" ;;
        g) gtf_file="$OPTARG" ;;
        o) output_prefix="$OPTARG" ;;
        O) count_overlaps=1 ;;
        T) threads="$OPTARG" ;;
        \?) echo "无效的选项: -$OPTARG" >&2; usage ;;
        :) echo "选项 -$OPTARG 需要一个参数" >&2; usage ;;
    esac
done

# 检查必需参数
if [ -z "$input_file" ] || [ -z "$gtf_file" ] || [ -z "$output_prefix" ]; then
    echo "错误: 缺少必需参数"
    usage
fi

# 检查输入文件是否存在
if [ ! -f "$input_file" ]; then
    echo "错误: 输入文件 $input_file 不存在"
    exit 1
fi

# 检查GTF文件是否存在
if [ ! -f "$gtf_file" ]; then
    echo "错误: GTF文件 $gtf_file 不存在"
    exit 1
fi

# 创建临时目录
temp_dir=$(mktemp -d)
echo "临时目录创建在: $temp_dir"

# 初始化关联数组
declare -A library_map

# 读取输入文件并处理
while IFS=$'\t' read -r bam_path library_type single_end; do
    # 跳过标题行
    if [[ "$bam_path" == "BamPath" ]]; then
        continue
    fi
    
    # 检查BAM文件是否存在
    if [ ! -f "$bam_path" ]; then
        echo "警告: BAM文件 $bam_path 不存在，跳过"
        continue
    fi
    
    # 构造分组键
    key="${library_type}_${single_end}"
    
    # 将BAM文件路径添加到对应的分组
    library_map[$key]+="$bam_path "
done < "$input_file"

# 对每个分组运行featureCounts
for key in "${!library_map[@]}"; do
    # 解析分组信息
    IFS='_' read -r library_type single_end <<< "$key"
    bam_files=${library_map[$key]}
    
    # 构造输出文件名
    output_file="${output_prefix}_${key}.txt"
    
    # 构造featureCounts命令基本参数
    cmd="featureCounts -T $threads -a $gtf_file -o $output_file"
    
    # 添加链特异性参数
    if [ "$library_type" -eq "1" ]; then
        cmd+=" -s 1"
    elif [ "$library_type" -eq "2" ]; then
        cmd+=" -s 2"
    fi
    
    # 添加单端/双端参数
    if [ "$single_end" -eq "0" ]; then
        cmd+=" -p --countReadPairs"
    fi
    
    # 添加重叠reads计数参数
    if [ "$count_overlaps" -eq "1" ]; then
        cmd+=" -O"
    fi
    
    # 添加BAM文件
    cmd+=" $bam_files"
    
    # 执行命令
    echo "运行命令: $cmd"
    eval $cmd
    
    # 检查命令是否成功执行
    if [ $? -ne 0 ]; then
        echo "错误: featureCounts命令执行失败"
        exit 1
    fi
    
    # 处理featureCounts输出文件，去掉路径只保留文件名
    sed -i '1s|.*/||' "${output_file}.summary"
    sed -i '1s|.*/||' "$output_file"
done

# 合并所有结果
echo "合并所有结果到 ${output_prefix}_combined.txt"

# 初始化变量
first_file=1
base_cols=6  # featureCounts结果中前6列是基本信息

# 处理每个分组的结果文件
for key in "${!library_map[@]}"; do
    output_file="${output_prefix}_${key}.txt"
    
    # 跳过summary文件
    if [[ "$output_file" == *"summary"* ]]; then
        continue
    fi
    
    # 处理第一个文件
    if [ "$first_file" -eq 1 ]; then
        # 提取前6列基本信息和所有计数列
        awk -v base_cols="$base_cols" '
        NR==2 {
            # 处理标题行
            printf $1; 
            for(i=2;i<=base_cols;i++) printf "\t"$i; 
            for(i=base_cols+1;i<=NF;i++) printf "\t"$i; 
            print ""
        } 
        NR>2 {
            printf $1; 
            for(i=2;i<=base_cols;i++) printf "\t"$i; 
            for(i=base_cols+1;i<=NF;i++) printf "\t"$i; 
            print ""
        }' "$output_file" > "${output_prefix}_combined.txt"
        
        first_file=0
    else
        # 对于后续文件，只提取计数列（第7列及以后）
        temp_file="${temp_dir}/temp_counts.txt"
        
        awk -v base_cols="$base_cols" '
        NR==2 {
            # 处理标题行
            for(i=base_cols+1;i<=NF;i++) printf "\t"$i; 
            print ""
        } 
        NR>2 {
            for(i=base_cols+1;i<=NF;i++) printf "\t"$i; 
            print ""
        }' "$output_file" > "$temp_file"
        
        # 使用paste命令合并文件
        paste "${output_prefix}_combined.txt" "$temp_file" > "${temp_dir}/temp_combined.txt"
        mv "${temp_dir}/temp_combined.txt" "${output_prefix}_combined.txt"
    fi
done

# 清理临时文件
echo "清理临时目录..."
rm -rf "$temp_dir"

echo "处理完成!"
