#!/usr/bin/env python
#输入小肽ID，以及genePred，得到小肽的genePred，后续可以用于生成protein fasta序列的gtf文件
#Usage: python Uni.generate_sorf_genepred.py sorfs.id.txt gencode.*.annotation.genePred.txt *gpe
#sorfs.id.txt 包含sorf的id，例如ENSMUST00000039333.10+chr8:111821392-111828499

import sys
import pandas as pd
import numpy as np
import re

# 输入文件
inputfile_1 = pd.read_csv(sys.argv[1], sep='\t', header='infer')
# inputfile_1 = pd.read_csv("./gen_genepred/essen_info_to_gen_genepred.sorfs.txt", sep='\t', header='infer')
inputfile_2 = pd.read_csv(sys.argv[2], sep='\t', header=None)
# inputfile_2 = pd.read_csv("/home/user/data3/lit/project/sORFs/01-ribo-seq/output/assembled_trans/stringtie_output/gencode.vM29.add_assemble.genePred.txt", sep='\t', header=None)

# 定义一个函数来提取所需的元素
def extract_elements(orf_id):
    # 使用 re.split() 以多个分隔符分割字符串
    parts = re.split('[:+-]', orf_id)
    # 提取第 1, 3, 4 个元素
    if len(parts) >= 4:
        return parts[0], parts[2], parts[3]
    return None, None, None

# 应用函数并创建新列
inputfile_1[['ENS_id', 'Start', 'End']] = inputfile_1['ORF_id_trans'].apply(extract_elements).apply(pd.Series)

# 合并，构建genepred
merged_dataset=pd.merge(inputfile_1,inputfile_2,how='left',left_on="ENS_id",right_on=0)
final_genepred=merged_dataset.loc[:,["ORF_id_trans",1,2,3,4,"Start","End",7,8,9]]
for col in [3,4,7]:
    final_genepred.iloc[:, col] = final_genepred.iloc[:, col].astype(int) 
# 导出结果
outputfile=sys.argv[3]
# outputfile="./gen_genepred/essen_info_to_gen_genepred.sorfs.gpe"
final_genepred.to_csv(outputfile,header=False,index=False,sep='\t')

