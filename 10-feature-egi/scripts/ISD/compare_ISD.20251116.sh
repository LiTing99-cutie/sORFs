#!/bin/bash
################################################
# File Name: run.sh
# Author: rbase    
# Mail: xiaochunfu@stu.pku.edu.cn
# Created Time: Thu Oct 31 15:08:42 2024
# Usage: bash run.sh input1.fa [input2.fa ...]
################################################

# === 基本参数自行修改 ===
workDir=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/ISD
script=/home/user/data3/rbase/opt/aiupred/aiupred.py
calculate_average_script=/home/user/data3/rbase/denovo_tumor/denovo_genes/ISD/calculate_average.sh

outputDir=$workDir/ISD_per_aa
fastaDir=$workDir/ISD_per_aa/fasta
avgDir=$workDir/ISD_avg

mkdir -p "$outputDir" "$fastaDir" "$avgDir"

# 激活环境
source activate PrimeNovo-1

# 只计算 aa_no_C
type="aa_no_C"

if [ $# -lt 1 ]; then
    echo "Usage: $0 input1.fa [input2.fa ...]"
    exit 1
fi

for fa in "$@"; do
    if [ ! -f "$fa" ]; then
        echo "WARNING: $fa not found, skip."
        continue
    fi

    base=$(basename "$fa")
    prefix=${base%.*}            # 去掉扩展名，例如 xxx.fa -> xxx

    echo "==== $fa -> $prefix.$type ===="

    # 去掉终止密码子 *，并去掉 C
    out_fa="$fastaDir/${prefix}.${type}.fa"
    sed 's/*$//g' "$fa" | sed 's/C//g' > "$out_fa"

    # 运行 IUPred
    isd_out="$outputDir/${prefix}.${type}.ISD.out"
    if [ ! -f "$isd_out" ]; then
        python "$script" -i "$out_fa" -o "$isd_out" -v -g 1
    else
        echo "ISD 结果已存在：$isd_out，跳过计算"
    fi

    # 计算平均 ISD
    avg_out="$avgDir/${prefix}.${type}.ISD_avg.txt"
    bash "$calculate_average_script" "$isd_out" "$avg_out"
done

echo "All done!"
