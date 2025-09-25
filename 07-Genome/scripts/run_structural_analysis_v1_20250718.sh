#!/bin/bash

# 结构分类分析运行脚本 v1_20250718

set -e

# 激活R环境（如果需要）
# source activate r_env

echo "开始结构分类统计分析..."

# 运行R脚本
Rscript /home/user/data3/lit/project/sORFs/07-Genome/scripts/analyze_structural_categories_v1_20250718.R

echo "分析完成！" 