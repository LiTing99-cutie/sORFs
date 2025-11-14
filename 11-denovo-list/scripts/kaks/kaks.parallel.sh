#!/bin/bash
source activate biotools
denovo_dir=$1
output_base=$2

# 创建输出目录
mkdir -p $output_base

# 使用 GNU parallel 并行运行
find $denovo_dir -name "*ortholog.pep.fa" | \
parallel -j 10 '
    pep_file={}
    orf=$(basename {} .ortholog.pep.fa)
    nucl_file=$(dirname {})/$(basename {} .pep.fa).nucl.fa
    output_dir='$output_base'/$orf
    
    echo "Processing: $orf"
    
    bash ./kaks/kaks.bin.sh \
        -p $pep_file \
        -c $nucl_file \
        -o $output_dir \
        -m human \
        -t 8
'

# 合并所有 Ka/Ks 结果
{
    echo -e "ORF\t$(head -1 $(find $output_base -name "kaks.res" | head -1))"
    find $output_base -name "kaks.res" | \
    while read file; do
        orf=$(basename $(dirname $file))
        awk -v orf="$orf" 'NR>1 {print orf"\t"$0}' $file
    done
} > $output_base/all_kaks_results.txt

# 计算中位数和均值
python3 kaks/calc_kaks_stats.py $output_base/all_kaks_results.txt > $output_base/all_kaks_results.summary.txt
