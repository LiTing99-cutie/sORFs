import sys
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# 输入文件
inputfile1 = sys.argv[1]  # repre.valid.pred.pvalue.parameters.txt
inputfile2 = sys.argv[2]  # repre.valid.ORF.genepred.txt
out1 = sys.argv[3]        # repre.valid.ORF.genepred_formatted.txt

# 读取 RibORF1.0 的预测结果文件
inputfile_inputpred3 = pd.read_csv(inputfile1, sep='\t', header='infer')
inputfile_inputpred3.columns = ['orfID', 'chrom', 'strand', 'codon5', 'codon3', 'length', 'readNum', 'f1', 'f2', 'f3', 'entropy', 'Maxentropy', 'PME', 'codonNum', 'f1Num','f1max', 'pred.pvalue']

# 过滤条件
# inputfile_inputpred1 = inputfile_inputpred[inputfile_inputpred['pred.pvalue']>=0.7]
# inputfile_inputpred2 = inputfile_inputpred1[inputfile_inputpred1['length']<=453]
# inputfile_inputpred3 = inputfile_inputpred2[inputfile_inputpred2['length'] >= 21]

# 调整 codon5 和 codon3 的位置
# inputfile_inputpred3['codon5'] = np.where(inputfile_inputpred3['strand'] == '+', inputfile_inputpred3['codon5'] + 1, inputfile_inputpred3['codon5'])
# inputfile_inputpred3['codon3'] = np.where(inputfile_inputpred3['strand'] == '+', inputfile_inputpred3['codon3'] - 3, inputfile_inputpred3['codon3'])
# inputfile_inputpred3['codon5'] = np.where(inputfile_inputpred3['strand'] == '-', inputfile_inputpred3['codon5'] + 4, inputfile_inputpred3['codon5'])

# 合并 orfID 信息
Merged_orfID = inputfile_inputpred3['orfID'].str.split(":", expand=True)
Merged_formatted = Merged_orfID[0].astype(str) + inputfile_inputpred3['strand'].astype(str) + Merged_orfID[1].astype(str) + ":" + inputfile_inputpred3['codon5'].astype(str) + "-" + inputfile_inputpred3['codon3'].astype(str)

inputfile_inputpred3['ORFid_format'] = Merged_formatted

# 删除不需要的列
finalinput = inputfile_inputpred3.drop(columns=['chrom', 'strand', 'codon5', 'codon3', 'length', 'readNum', 'f1', 'f2', 'f3', 'entropy', 'Maxentropy', 'PME', 'codonNum', 'f1max', 'pred.pvalue'])

# 读取 genepred 文件
genepred_inputfile = pd.read_csv(inputfile2, sep="\t", header=None)
genepred_inputfile.columns = ['orfID', 'chr', 'strand', 'trans_start', 'trans_stop', 'cds_start', 'cds_stop', 'exon', 'exon_start', 'exon_stop']

# 合并数据
genepred_merge = pd.merge(finalinput, genepred_inputfile, how='inner', on='orfID')
genepred_final = genepred_merge['ORFid_format'].astype(str) + "\t" + genepred_merge['chr'].astype(str) + "\t" + genepred_merge['strand'].astype(str) + "\t" + genepred_merge['trans_start'].astype(str) + "\t" + genepred_merge['trans_stop'].astype(str) + "\t" + genepred_merge['cds_start'].astype(str) + "\t" + genepred_merge['cds_stop'].astype(str) + "\t" + genepred_merge['exon'].astype(str) + "\t" + genepred_merge['exon_start'].astype(str) + "\t" + genepred_merge['exon_stop'].astype(str)

# 输出 genepred 文件
output_file1 = open(out1, 'w')
for l in genepred_final:
    output_file1.write(str(l) + "\n")
output_file1.close()