#!/usr/bin/env bash
set -euo pipefail

# ==================== 配置部分 ====================
# 脚本路径配置
script_dir="stat_scripts"
Mapping_Statistics_script="$script_dir/Uni.Run.Mapping_Statistics.20250518.R"
# 输出目录配置
# res_dir=$1
# # 参数一：结果路径
# res_dir=$(realpath ../processed/batch_1/)
# bam_dir=$res_dir/filtered_bam
# # 参数二：输出结果路径
# output_dir=$(realpath ../results/batch_1/qual_assess)
# # 参数三：原始fastq路径
# raw_fastq_lst=../processed/raw_fastq_batch_1.20250723.lst

# 参数一：结果路径
res_dir=$1
bam_dir=$res_dir/filtered_bam
# 参数二：输出结果路径
output_dir=$res_dir/qual_assess
# 参数三：原始fastq路径
# raw_fastq_lst=../processed/raw_fastq_batch_1.20250723.lst

mkdir -p "$output_dir"
mkdir -p "$output_dir/ribotish"

# 工具配置
source "/home/user/data2/lit/bin/lit_utils.sh"
define_annotation_gencode_v41_human

# ==================== 功能函数 ====================
# 检查命令是否存在
check_command() {
    if ! command -v "$1" &> /dev/null; then
        echo "错误: 未找到命令 $1"
        exit 1
    fi
}

# 检查文件是否存在
check_file() {
    if [ ! -f "$1" ]; then
        echo "错误: 文件 $1 不存在"
        exit 1
    fi
}

# ==================== 主流程 ====================
echo "====== 开始处理 ======"

# 0. 生成样本列表
echo "生成样本列表..."
ana_dir=/home/user/data/lit/project/sORFs/01-ribo-seq/analysis
find "$ana_dir/20250722_formal_data_run" "$ana_dir/20250813_demo_data_analysis" -name "fq.stat.txt" > "$output_dir/fq.stat.lst"

# 1. 质控统计
echo "处理质控数据..."
# 1.0 原始reads数量
check_command seqkit
raw_reads_1=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/results/batch_1/qual_assess/raw_reads.txt
raw_reads_2=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250813_demo_data_analysis/results/batch_1/qual_assess/raw_reads.txt
cat $raw_reads_1 $raw_reads_2 > "$output_dir/raw_reads.txt"

# 1.1 被其他小RNA污染的比例
check_command Rscript
check_file "$Mapping_Statistics_script"
Rscript "$Mapping_Statistics_script" "$output_dir/fq.stat.lst" "$output_dir"

