#!/usr/bin/env python3

import sys
import statistics

# 读取文件
input_file = sys.argv[1] if len(sys.argv) > 1 else 'kaks.tsv'

# 存储每个ORF的Ka/Ks值
orf_data = {}

with open(input_file, 'r') as f:
    # 跳过表头
    next(f)
    
    for line in f:
        cols = line.strip().split('\t')
        orf = cols[0]
        ka_ks = cols[4]
        
        # 跳过NA值
        if ka_ks != 'NA' and ka_ks != '-0':
        # if ka_ks != 'NA':
            ka_ks_value = float(ka_ks)
            if orf not in orf_data:
                orf_data[orf] = []
            orf_data[orf].append(ka_ks_value)

# 计算并输出结果
print("ORF\tMedian_Ka_Ks\tMean_Ka_Ks")
for orf, values in orf_data.items():
    if values:
        median = statistics.median(values)
        mean = statistics.mean(values)
        print(f"{orf}\t{median:.6f}\t{mean:.6f}")