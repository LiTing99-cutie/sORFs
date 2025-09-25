#!/bin/bash
set -e

# 激活biotools环境
source ~/anaconda3/bin/activate biotools

# 本脚本仅执行VEP注释（第6步）
# 注意：VEP cache数据库是手动从国内镜像下载并上传到服务器的
# 请确保cache已解压到~/.vep或指定目录

SAMPLE=human_brain_21pcw
RESULTS=../results
THREADS=32  # 可根据服务器负载调整

# VEP注释
vep -i $RESULTS/$SAMPLE.vcf.gz \
    -o $RESULTS/$SAMPLE.vep.vcf \
    --vcf \
    --cache \
    --species homo_sapiens \
    --assembly GRCh38 \
    --fork $THREADS

echo "VEP注释完成！" 