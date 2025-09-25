#!/usr/bin/env python
#输入小肽ID，转录本id，ORF在基因组上的起始位点和终止位点，以及注释的genePred（包括组装的转录本），得到小肽的genePred，后续可以用于生成protein fasta序列的gtf文件
#Usage: python Uni.generate_sorf_genepred.py essen_info_to_gen_genepred.sorfs.txt gencode.*.annotation.genePred.txt *gpe
#essen_info_to_gen_genepred.sorfs.txt 包含"Protein","ENS_id","Start","End"四列

import sys
import pandas as pd
import numpy as np

# 输入文件
inputfile_1 = pd.read_csv(sys.argv[1], sep='\t', header='infer')
# inputfile_1 = pd.read_csv("./gen_genepred/essen_info_to_gen_genepred.sorfs.txt", sep='\t', header='infer')
inputfile_2 = pd.read_csv(sys.argv[2], sep='\t', header=None)
# inputfile_2 = pd.read_csv("/home/user/data3/lit/project/sORFs/01-ribo-seq/output/assembled_trans/stringtie_output/gencode.vM29.add_assemble.genePred.txt", sep='\t', header=None)

# 合并，构建genepred
merged_dataset=pd.merge(inputfile_1,inputfile_2,how='left',left_on="ENS_id",right_on=0)
final_genepred=merged_dataset.loc[:,["ORF_id_trans",1,2,3,4,"Start","End",7,8,9]]
for col in [3,4,7]:
    final_genepred.iloc[:, col] = final_genepred.iloc[:, col].astype(int) 
# 导出结果
outputfile=sys.argv[3]
# outputfile="./gen_genepred/essen_info_to_gen_genepred.sorfs.gpe"
final_genepred.to_csv(outputfile,header=False,index=False,sep='\t')

