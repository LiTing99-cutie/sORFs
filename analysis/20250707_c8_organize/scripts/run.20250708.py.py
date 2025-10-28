#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import re
pd.set_option('display.max_colwidth', 1000)
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Arial']

def load_msf(path, source):
    df = pd.read_csv(path, sep='\t', low_memory=False)
    df = df[['Spectrum', 'Peptide']]
    df['Peptide_I_L_equal'] = df['Peptide'].str.replace('I', 'L')
    df['Source'] = source
    return df
def load_pfind(path, source):
    df = pd.read_csv(path, sep='\t', header=0, low_memory=False)
    df = df[['File_Name', 'Sequence']]
    df.columns = ['Spectrum', 'Peptide']
    df['Peptide_I_L_equal'] = df['Peptide'].str.replace('I', 'L')
    df['Source'] = source
    return df
def format_msf_spectrum(psm):
    # 拆分 Spectrum 字符串
    split = psm['Spectrum'].str.split('.', n=3, expand=True)
    prefix = split[0]
    scan_number = split[1].str.lstrip('0')  # 去除前导0
    charge = psm['Charge'] if 'Charge' in psm.columns else split[3]  # 兼容Charge列

    # 重新拼接 Spectrum
    psm['Spectrum'] = prefix + '.' + scan_number + '.' + scan_number + '.' + charge.astype(str)
    return psm
def format_pfind_spectrum(df):
    # 去掉_SCANS及后面的内容
    df['Spectrum'] = df['Spectrum'].str.replace(r'_SCANS\d+$', '', regex=True)
    return df
def mean_score(scores):
        if pd.isna(scores):
            return np.nan
        arr = [float(x) for x in scores.split(',')]
        return np.mean(arr)
def load_mztab(path, source):
    # 1. 读取第60行，提取prefix并去掉_uncalibrated
    with open(path) as f:
        for i, line in enumerate(f):
            if i == 59:  # 第60行，0-based
                mgf_path = line.strip().split('\t')[-1]
                mgf_file = os.path.basename(mgf_path.replace("file://", ""))
                prefix = os.path.splitext(mgf_file)[0].replace('_uncalibrated', '')
                break

    # 2. 读取数据部分（header自动识别，跳过前60行）
    df = pd.read_csv(path, sep='\t', header=0, skiprows=60, low_memory=False)
    # 3. 计算肽段水平的置信分数
    df['Score'] = df['opt_ms_run[1]_aa_scores'].apply(mean_score)
    # 4. 选取需要的列
    df = df[['sequence', 'PSM_ID', 'charge','Score']]
    # 5. 生成Spectrum列
    def safe_int_str(x):
        try:
            if float(x).is_integer():
                return str(int(x))
            else:
                return str(x)
        except:
            return str(x)

    df['Spectrum'] = (
        prefix + '.' +
        df['PSM_ID'].apply(safe_int_str) + '.' +
        df['PSM_ID'].apply(safe_int_str) + '.' +
        df['charge'].apply(safe_int_str)
    )
    def remove_mods(seq):
        # 去掉所有+数字.数字（如+57.021、+15.995等）
        return re.sub(r'\+\d+\.\d+', '', seq)

    # 在load_mztab里，df['sequence']处理后再赋值
    df['sequence'] = df['sequence'].apply(remove_mods)       
    # 6. 重命名列，添加Peptide_I_L_equal和Source
    df = df.rename(columns={'sequence': 'Peptide'})
    df['Peptide_I_L_equal'] = df['Peptide'].str.replace('I', 'L')
    df["Source"] = source 

    return df[['Spectrum', 'Peptide', 'Peptide_I_L_equal', 'Score','Source']]
def load_primenovo(path, source):
    # 读取 tsv 文件
    df = pd.read_csv(path, sep='\t', header=0, low_memory=False)
    
    # 处理 prediction 列，移除所有方括号及其内部内容和所有连字符
    df['Peptide'] = df['prediction'].str.replace(r'\[.*?\]|-', '', regex=True)
    
    # 处理 label 列，去掉最后一个下划线及其后内容
    df['Spectrum'] = df['label'].str.replace(r'_[^_]*$', '', regex=True)
    
    # Peptide_I_L_equal 列
    df['Peptide_I_L_equal'] = df['Peptide'].str.replace('I', 'L')
    
    # Source 列
    df['Source'] = source
    df = df.rename(columns={'score': 'Score'})
    # 只保留主要列
    return df[['Spectrum', 'Peptide', 'Peptide_I_L_equal', 'Score','Source']]
def load_punifind(path, source):
    # 读取文件，不指定表头
    df = pd.read_csv(path, sep=',', header=None, low_memory=False)
    
    # 选取需要的列
    df = df[[0, 1, 7]]  # 第一列、第二列、第八列（0-based）
    df.columns = ['Spectrum', 'Score', 'V8']
    
    # 处理肽段列，去掉第一个下划线及其后内容
    df['Peptide'] = df['V8'].str.replace(r'_.*', '', regex=True)
    
    # Peptide_I_L_equal 列
    df['Peptide_I_L_equal'] = df['Peptide'].str.replace('I', 'L')
    
    # Source 列
    df['Source'] = source
    
    # 只保留主要列
    return df[['Spectrum',  'Peptide', 'Peptide_I_L_equal', 'Score','Source']]
def load_novor(path, source):
    # 跳过注释行，找到表头
    with open(path) as f:
        for i, line in enumerate(f):
            if line.startswith('# id'):
                header_line = i
                break

    # 读取数据部分
    df = pd.read_csv(
        path,
        skiprows=header_line+1,
        header=None,
        names=[
            'id', 'scanNum', 'RT', 'mz', 'z', 'pepMass', 'err', 'ppm', 'score', 'peptide', 'aaScore'
        ]
    )

    # 生成 Spectrum 列
    prefix = None
    with open(path) as f:
        for line in f:
            if line.startswith('# input file ='):
                mgf_file = line.split('=')[1].strip()
                import os
                prefix = os.path.splitext(os.path.basename(mgf_file))[0].replace('_uncalibrated', '')
                break
    if prefix is None:
        prefix = 'unknown'

    df['Spectrum'] = (
        prefix + '.' +
        df['scanNum'].astype(str) + '.' +
        df['scanNum'].astype(str) + '.' +
        df['z'].astype(str)
    )

    # 去除Peptide中的括号内容
    df['Peptide'] = df['peptide'].str.replace(r'\([^)]*\)', '', regex=True)

    # Peptide_I_L_equal 列
    df['Peptide_I_L_equal'] = df['Peptide'].str.replace('I', 'L')
    df['Source'] = source
    df = df.rename(columns={'score': 'Score'})
    # 输出包含score
    return df[['Spectrum', 'Peptide', 'Peptide_I_L_equal', 'Score', 'Source']]


# # DB search

# In[2]:


msf_closed = load_msf('../output/db_search/msfragger/closed/psm.tsv', 'msfragger_closed')
msf_closed = format_msf_spectrum(msf_closed)
msf_open = load_msf('../output/db_search/msfragger/open/psm.tsv', 'msfragger_open')
msf_open = format_msf_spectrum(msf_open)
pfind_closed = load_pfind('../output/db_search/pfind/closed/pFind-Filtered.spectra', 'pfind_closed')
pfind_closed = format_pfind_spectrum(pfind_closed)
pfind_open = load_pfind('../output/db_search/pfind/open/pFind-Filtered.spectra', 'pfind_open')
pfind_open = format_pfind_spectrum(pfind_open)


# In[3]:


# msfragger 以 closed 为主
msf_open_extra = msf_open[~msf_open['Spectrum'].isin(msf_closed['Spectrum'])]
msfragger_all = pd.concat([msf_closed, msf_open_extra], ignore_index=True)

# pfind 以 closed 为主
pfind_open_extra = pfind_open[~pfind_open['Spectrum'].isin(pfind_closed['Spectrum'])]
pfind_all = pd.concat([pfind_closed, pfind_open_extra], ignore_index=True)

# 以 msfragger 为主，pfind 补充
pfind_extra = pfind_all[~pfind_all['Spectrum'].isin(msfragger_all['Spectrum'])]
final = pd.concat([msfragger_all, pfind_extra], ignore_index=True)


# In[4]:


all_results = pd.concat([msf_closed, msf_open, pfind_closed, pfind_open], ignore_index=True)
all_results.to_csv('../results/all_dbsearch.csv', index=False)
# 导出到指定文件夹
msfragger_all.to_csv('../results/msfragger_all.csv', index=False)
pfind_all.to_csv('../results/pfind_all.csv', index=False)
final.to_csv('../results/final.csv', index=False)


# In[5]:


final.Source.value_counts()


# In[6]:


final.head()


# In[9]:


gold_df=final


# # De novo search

# In[51]:


casanovo = load_mztab("../output/denovo/casanovo/results.mztab","casanovo")
print(casanovo.head())
primenovo = load_primenovo('../output/denovo/primenovo/denovo.tsv', 'primenovo')
print(primenovo.head())
punifind = load_punifind('../output/denovo/punifind/tims_all_5_merged.csv', 'punifind')
print(punifind.head())
novor = load_novor('../output/denovo/novor/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_uncalibrated.novor.csv', 'novor')
print(novor.head())


# In[52]:


# 统计交集数量
counts = []
method_names = ['casanovo', 'primenovo', 'punifind']
all_methods = [casanovo, primenovo, punifind]

for df in all_methods:
    gold_peptides = set(gold_df['Peptide_I_L_equal'])
    method_peptides = set(df['Peptide_I_L_equal'])
    intersection = gold_peptides & method_peptides
    counts.append(len(intersection))

# Create DataFrame
result_df = pd.DataFrame({
    'Method': method_names,
    'Identified_peptide_sequence': counts
})

# Save to file
result_df.to_csv('../results/denovo_identified_peptide.txt', sep='\t', index=False)


# In[55]:


result_df


# In[56]:


len(gold_peptides)


# In[12]:


# casanovo没有输出谱图ID信息，无法进行对应
def get_recall_coverage_curve(pred_df, gold_df, score_col='Score', peptide_col='Peptide_I_L_equal', spectrum_col='Spectrum'):
    pred_df = pred_df[pred_df[spectrum_col].isin(gold_df[spectrum_col])]
    gold_map = gold_df.set_index(spectrum_col)[peptide_col].to_dict()
    pred_df = pred_df.sort_values(score_col, ascending=False).reset_index(drop=True)
    pred_df['is_correct'] = pred_df.apply(
        lambda row: row[peptide_col] == gold_map.get(row[spectrum_col], None), axis=1
    )
    total = len(gold_df)
    pred_df['coverage'] = (np.arange(1, len(pred_df)+1)) / total
    pred_df['recall'] = pred_df['is_correct'].cumsum() / (np.arange(1, len(pred_df)+1))
    return pred_df[['coverage', 'recall']]

# 计算并导出
primenovo_curve = get_recall_coverage_curve(primenovo, final)
primenovo_curve['Method'] = 'primenovo'
primenovo_curve.to_csv('../results/primenovo_recall_coverage.tsv', sep='\t', index=False)

punifind_curve = get_recall_coverage_curve(punifind, final)
punifind_curve['Method'] = 'punifind'
punifind_curve.to_csv('../results/punifind_recall_coverage.tsv', sep='\t', index=False)


# In[54]:


import numpy as np
import pandas as pd

def get_recall_coverage_curve(pred_df, gold_df, score_col='Score', peptide_col='Peptide_I_L_equal', spectrum_col='Spectrum', n_points=100):
    # 只保留 gold_df 中存在的 spectrum
    pred_df = pred_df[pred_df[spectrum_col].isin(gold_df[spectrum_col])].copy()
    gold_map = gold_df.set_index(spectrum_col)[peptide_col].to_dict()
    
    # 按分数排序
    pred_df = pred_df.sort_values(score_col, ascending=False).reset_index(drop=True)
    
    # 标记是否预测正确
    pred_df['is_correct'] = pred_df.apply(
        lambda row: row[peptide_col] == gold_map.get(row[spectrum_col], None), axis=1
    )
    
    total = len(gold_df)
    pred_df['coverage'] = (np.arange(1, len(pred_df)+1)) / total
    pred_df['recall'] = pred_df['is_correct'].cumsum() / (np.arange(1, len(pred_df)+1))
    
    # 按等间距采样 n_points 个点
    if len(pred_df) > n_points:
        idx = np.linspace(0, len(pred_df) - 1, n_points, dtype=int)
        pred_df = pred_df.iloc[idx]
    
    return pred_df[['coverage', 'recall']]

# 计算并导出
primenovo_curve = get_recall_coverage_curve(primenovo, final, n_points=100)
primenovo_curve['Method'] = 'primenovo'
primenovo_curve.to_csv('../results/primenovo_recall_coverage.tsv', sep='\t', index=False)

punifind_curve = get_recall_coverage_curve(punifind, final, n_points=100)
punifind_curve['Method'] = 'punifind'
punifind_curve.to_csv('../results/punifind_recall_coverage.tsv', sep='\t', index=False)


# ##  合并primenovo和punifind的结果

# In[13]:


merged = pd.merge(
    primenovo,
    punifind,
    on=['Spectrum', 'Peptide_I_L_equal'],
    suffixes=('_primenovo', '_punifind')
)


# In[47]:


print(merged.shape)
print(primenovo.shape)
print(punifind.shape)


# In[14]:


# merged 已有 Spectrum 和 Peptide_I_L_equal
# gold_df 也有 Spectrum 和 Peptide_I_L_equal

# 只保留 merged 中 Spectrum 存在于 gold_df 的
merged_in_gold = merged[merged['Spectrum'].isin(gold_df['Spectrum'])]

# 构建 gold_df 的映射
gold_map = gold_df.set_index('Spectrum')['Peptide_I_L_equal'].to_dict()

# 判断 merged 的肽段是否和 gold_df 一致
merged_in_gold['is_correct'] = merged_in_gold.apply(
    lambda row: row['Peptide_I_L_equal'] == gold_map.get(row['Spectrum'], None), axis=1
)

recall = merged_in_gold['is_correct'].sum() / len(merged_in_gold)
print(f"Recall of merged: {recall:.4f}")


# In[15]:


# 假设 merged 已经有 Spectrum 和 Peptide_I_L_equal
# 选择一列 Peptide（任选 primenovo 或 punifind 的 Peptide）
if 'Peptide_primenovo' in merged.columns:
    merged['Peptide'] = merged['Peptide_primenovo']
elif 'Peptide' in merged.columns:
    # 如果merge时没有加后缀，直接用
    pass
else:
    # 其它情况请补充
    pass

# 添加 Source 列
merged['Source'] = 'primenovo_punifind'

# 只保留需要的四列
merged_for_final = merged[['Spectrum', 'Peptide', 'Peptide_I_L_equal', 'Source']]

# 只保留 merged_for_final 里 Spectrum 不在 final 里的部分
merged_new = merged_for_final[~merged_for_final['Spectrum'].isin(final['Spectrum'])]

# 合并
final_plus_merged = pd.concat([final, merged_new], ignore_index=True)


# In[21]:


# 导出结果
final_plus_merged.to_csv('../results/final_pep_search_res_20250721.tsv', sep='\t', index=False)


# In[17]:


mgf_path="/rd1/user/lit/project/sORFs/analysis/20250707_c8_organize/input/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_uncalibrated.mgf"
mgf_spectra_count = 0
with open(mgf_path) as f:
    for line in f:
        if line.strip() == "BEGIN IONS":
            mgf_spectra_count += 1
print("Total spectra in mgf:", mgf_spectra_count)


# In[50]:


len(final_plus_merged)/mgf_spectra_count


# In[18]:


source_counts = final_plus_merged['Source'].value_counts().reset_index()
source_counts.columns = ['Source', 'Count']
identified_total = source_counts['Count'].sum()
unidentified_count = mgf_spectra_count - identified_total

# 增加一行
source_counts = pd.concat([
    source_counts,
    pd.DataFrame([{'Source': 'Unidentified due to Low quality', 'Count': unidentified_count}])
], ignore_index=True)

rename_dict = {
    "msfragger_closed": "MSFragger closed search",
    "msfragger_open": "MSFragger open search",
    "pfind_closed": "pFind closed search",
    "pfind_open": "pFind open search",
    "primenovo_punifind": "PrimeNovo & pUniFind intersection",
    "Unidentified due to Low quality": "Unidentified (low quality spectra)"
}

source_counts['Source'] = source_counts['Source'].replace(rename_dict)
print(source_counts)

source_counts.to_csv('../results/source_summary.tsv', sep='\t', index=False)


# # 查看修饰类型

# ## pFind修饰

# In[22]:


import pandas as pd
from collections import defaultdict

def contains_other_modification(mod_str):
    """检查修饰字符串是否包含非标准修饰"""
    # 定义允许的标准修饰
    allowed_mods = {"Carbamidomethyl[C]", "Acetyl[ProteinN-term]", "Oxidation[M]"}
    
    # 空值处理
    if pd.isna(mod_str) or mod_str.strip() == "":
        return False
    
    # 分割修饰条目
    mods = [m.strip() for m in mod_str.split(';') if m.strip()]
    
    for mod in mods:
        # 分割位置和修饰类型
        if ',' in mod:
            # 只取修饰类型部分（逗号后面的内容）
            mod_type = mod.split(',', 1)[1].strip()
            # 检查是否是非标准修饰
            if mod_type not in allowed_mods:
                return True
    return False


def get_unique_non_standard_mods(mod_str):
    """获取每行中唯一的非标准修饰集合"""
    # 定义允许的标准修饰
    allowed_mods = {"Carbamidomethyl[C]", "Acetyl[ProteinN-term]", "Oxidation[M]"}
    unique_mods = set()
    
    if pd.isna(mod_str) or mod_str.strip() == "":
        return unique_mods
    
    # 分割修饰条目
    mods = [m.strip() for m in mod_str.split(';') if m.strip()]
    
    for mod in mods:
        if ',' in mod:
            mod_type = mod.split(',', 1)[1].strip()
            if mod_type not in allowed_mods:
                unique_mods.add(mod_type)
    
    return unique_mods


# In[39]:


df = pd.read_csv('../output/db_search/pfind/open/pFind-Filtered.spectra', sep='\t')

# 确保修改列名以匹配实际数据（这里假设列名为"Modification"）
mod_col = 'Modification'  # 根据实际列名修改

# 筛选包含非标准修饰的行
filtered_df = df[df[mod_col].apply(contains_other_modification)]

# 创建一个字典来统计每种修饰出现的行数
mod_counter = defaultdict(int)

# 遍历每一行，收集唯一的非标准修饰
for mod_str in filtered_df[mod_col]:
    unique_mods = get_unique_non_standard_mods(mod_str)
    for mod in unique_mods:
        mod_counter[mod] += 1

# 转换为DataFrame并排序
mod_stats = pd.DataFrame(
    list(mod_counter.items()), 
    columns=['Modification', 'Count']
)
mod_stats = mod_stats.sort_values('Count', ascending=False)

# 添加百分比列
total_rows = len(filtered_df)
mod_stats['Percentage'] = (mod_stats['Count'] / total_rows * 100).round(2)

# 保存统计结果
os.makedirs('../results/modi/', exist_ok=True)
mod_stats.to_csv('../results/modi/pfind_non_standard_modi_row_stats.csv', index=False)


# In[38]:


mod_stats.head()


# ## msFragger修饰

# In[43]:


# 1. 读取新的数据框（请替换为实际文件路径）
mf_modi_file = '/rd1/user/lit/project/sORFs/analysis/20250623_human_ms_rerun/output/db_search_c8_open_search/Trypsin/ptm-shepherd-output/global.modsummary.tsv'  # 修改为实际文件路径
mf_modi_df = pd.read_csv(mf_modi_file, sep='\t')

# 2. 重命名长列名以便操作
mf_modi_df = mf_modi_df.rename(columns={
    'CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2_3_1_7020_d_PSMs': 'Count',
    'CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2_3_1_7020_d_percent_PSMs': 'Percentage'
})

mf_modi_df.head()

mf_modi_df_1 = mf_modi_df.drop(0)

mf_modi_df_1.head()


# In[44]:


mf_modi_df_1.to_csv('../results/modi/msfragger_non_standard_modi_row_stats.csv', index=False)

