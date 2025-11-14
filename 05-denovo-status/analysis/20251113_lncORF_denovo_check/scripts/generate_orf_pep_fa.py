#!/usr/bin/env python3

import sys

orf_pep_dir = sys.argv[1]
id_mapping_file = sys.argv[2]
orfs_pep_fa = sys.argv[3]

# 读取映射
id_map = {}
with open(id_mapping_file) as f:
    for line in f:
        orig_id, safe_id = line.strip().split('\t')
        id_map[orig_id] = safe_id

# 处理FASTA
output_files = {}
current_id = None

with open(orfs_pep_fa) as f:
    for line in f:
        if line.startswith('>'):
            seq_id = line[1:].split()[0]
            if seq_id in id_map:
                current_id = id_map[seq_id]
                if current_id not in output_files:
                    output_files[current_id] = open(f"{orf_pep_dir}/{current_id}.ORF_pep.fa", 'w')
                output_files[current_id].write(line)
            else:
                current_id = None
        elif current_id:
            output_files[current_id].write(line)

# 关闭所有文件
for f in output_files.values():
    f.close()