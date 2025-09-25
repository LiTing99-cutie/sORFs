#得到一个ID统一化，且后续可以用于生成protein fasta序列的gtf文件
#20241212 增加命令行参数指定FDR阈值以及是否过滤codon为ATG
#Usage: python Generate_genepred_PRICE.py *.orfs.tsv gencode.*.annotation.genePred.txt 0.05 true nonCano.formatted.gpe 

import sys
import pandas as pd
import numpy as np

# 输入文件
inputfile_1 = pd.read_csv(sys.argv[1], sep='\t', header='infer')
inputfile_2 = pd.read_csv(sys.argv[2], sep='\t', header=None)

# 读取用户输入的 FDR 阈值和是否过滤 Codon 的标志
fdr_threshold = float(sys.argv[3])  # 第三个参数作为 FDR 阈值
filter_codon = float(sys.argv[4])  # 第四个参数表示是否过滤 Codon 为 ATG

# 校正p值
import statsmodels.stats.multitest as multi
inputfile_1['FDR']=multi.multipletests(inputfile_1['p value'].values, method="fdr_bh")[1]

# 基于用户输入进行筛选
query_string = f"FDR <= {fdr_threshold}"
if filter_codon:
    query_string += " and Codon == 'ATG'"
query_string += " and Type != 'CDS'"

# 筛选数据
inputfile_1 = inputfile_1.query(query_string)

print(inputfile_1.shape[0])

# 得到ORF在基因组上的起始位点和终止位点
import re
inputfile_1['gstart']=inputfile_1['Location'].apply(lambda x: min(int(num) for num in re.findall(r'[:|-](\d+)',  x)))
inputfile_1['gstop']=inputfile_1['Location'].apply(lambda x: max(int(num) for num in re.findall(r'[:|-](\d+)', x)))

# 得到染色体编号，加上chr字符；得到链的信息以及转录本id
inputfile_1['chr']="chr"+inputfile_1['Location'].apply(lambda x: re.match(r'(\w+)([+-]):', x).group(1))
inputfile_1['chr']=inputfile_1['chr'].replace("chrMT","chrM")
inputfile_1['strand']=inputfile_1['Location'].apply(lambda x: re.match(r'(\w+)([+-]):', x).group(2))
inputfile_1['ens_id']=inputfile_1['Id'].str.split("_", expand=True).iloc[:,0]

# 组建新的ORF_id名称
new_id=inputfile_1['ens_id']+inputfile_1['strand']+inputfile_1['chr']+":"+inputfile_1['gstart'].astype(str)+"-"+inputfile_1['gstop'].astype(str)
inputfile_1['new_id']=new_id

# 合并，构建genepred
merged_dataset=pd.merge(inputfile_1,inputfile_2,how='left',left_on="ens_id",right_on=0)
final_genepred=merged_dataset.loc[:,["new_id","chr","strand",3,4,"gstart","gstop",7,8,9]]

# 导出结果
outputfile=sys.argv[5]
final_genepred.to_csv(outputfile,header=False,index=False,sep='\t')