#!/bin/bash
set -e

# 脚本功能：将所有all_ribo_mapping_reads_20250518目录下的文件移动到data目录下，并删除源文件
# 版本：v1_20250710

SRC_DIR=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Ribo_classification_20250518/output/all_ribo_mapping_reads_20250518
DST_DIR=/home/user/data/lit/project/sORFs/01-ribo-seq/analysis/Ribo_classification_20250518/output/all_ribo_mapping_reads_20250518

mkdir -p $DST_DIR

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start moving files from $SRC_DIR to $DST_DIR"

nohup cp -r $SRC_DIR/* $DST_DIR/ && rm -rf $SRC_DIR/* &
