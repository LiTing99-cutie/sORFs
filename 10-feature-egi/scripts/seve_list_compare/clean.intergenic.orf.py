#!/usr/bin/env python3

import random
from Bio import SeqIO

# 设置随机种子以便重现
random.seed(42)

# 输入文件
base = "/home/user/data3/rbase/denovo_tumor/denovo_genes/all_ORFs/intergenic_ORFs"
bed_file = f"{base}/intergenic_ORFs.gt30.shuffle.bed"
nucl_file = f"{base}/intergenic_ORFs.gt30.shuffle.fa"
pep_file = f"{base}/intergenic_ORFs.gt30.shuffle.aa.fa"

# 输出文件前缀
out_prefix = "/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/intergenic_orfs.sample.2.5k.lt.150aa"

print("=" * 50)
print("处理intergenic ORFs")
print("=" * 50)

# 步骤1: 读取并筛选蛋白序列
print("\n步骤1: 读取并筛选蛋白序列...")
valid_seqs = []
total_count = 0

for record in SeqIO.parse(pep_file, "fasta"):
    total_count += 1
    # 清理序列名（去掉::之后的部分）
    clean_id = record.id.split("::")[0]
    seq_len = len(record.seq)
    
    if seq_len <= 150:
        record.id = clean_id
        record.description = ""
        valid_seqs.append(record)

print(f"总序列数: {total_count}")
print(f"≤150aa序列数: {len(valid_seqs)}")

# 步骤2: 随机抽取
print("\n步骤2: 随机抽取序列...")
n_sample = min(2500, len(valid_seqs))
selected = random.sample(valid_seqs, n_sample)
selected_ids = {s.id for s in selected}

print(f"抽取数量: {n_sample}")

# 步骤3: 保存蛋白序列
print("\n步骤3: 保存蛋白序列...")
with open(f"{out_prefix}.aa.fa", 'w') as f:
    SeqIO.write(selected, f, "fasta")

# 步骤4: 提取核苷酸序列
print("\n步骤4: 提取核苷酸序列...")
nucl_count = 0
with open(f"{out_prefix}.fa", 'w') as f_out:
    for record in SeqIO.parse(nucl_file, "fasta"):
        clean_id = record.id.split("::")[0]
        if clean_id in selected_ids:
            record.id = clean_id
            record.description = ""
            SeqIO.write(record, f_out, "fasta")
            nucl_count += 1

print(f"提取的核苷酸序列数: {nucl_count}")

# 步骤5: 提取bed记录
print("\n步骤5: 提取bed记录...")
bed_count = 0
with open(bed_file, 'r') as f_in, open(f"{out_prefix}.bed", 'w') as f_out:
    for line in f_in:
        fields = line.strip().split('\t')
        if len(fields) >= 4:
            # bed的第4列通常是名称
            bed_id = fields[3].split("::")[0]
            if bed_id in selected_ids:
                f_out.write(line)
                bed_count += 1

print(f"提取的bed记录数: {bed_count}")

# 步骤6: 统计和验证
print("\n" + "=" * 50)
print("完成！")
print("=" * 50)
print(f"\n输出文件:")
print(f"  蛋白序列: {out_prefix}.aa.fa")
print(f"  核苷酸序列: {out_prefix}.fa")
print(f"  BED文件: {out_prefix}.bed")

print(f"\n序列统计:")
lengths = [len(s.seq) for s in selected]
print(f"  序列数: {len(selected)}")
print(f"  最小长度: {min(lengths)} aa")
print(f"  最大长度: {max(lengths)} aa")
print(f"  平均长度: {sum(lengths)/len(lengths):.1f} aa")
print(f"  中位数长度: {sorted(lengths)[len(lengths)//2]} aa")

# 长度分布
print(f"\n长度分布:")
bins = [0, 50, 100, 150]
for i in range(len(bins)-1):
    count = sum(1 for l in lengths if bins[i] < l <= bins[i+1])
    print(f"  {bins[i]+1}-{bins[i+1]} aa: {count} ({count/len(lengths)*100:.1f}%)")