#!/usr/bin/env python3
"""
裁剪核苷酸序列到3的倍数并翻译为蛋白质序列
"""
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def trim_and_translate(nucl_file, pep_file):
    """
    读取核苷酸序列，裁剪到3的倍数，翻译为蛋白质序列
    """
    translated_records = []
    
    for record in SeqIO.parse(nucl_file, 'fasta'):
        # 裁剪序列到3的倍数
        seq_len = len(record.seq)
        trim_length = seq_len % 3
        
        if trim_length > 0:
            trimmed_seq = record.seq[:-trim_length]
        else:
            trimmed_seq = record.seq
        
        # 翻译序列
        try:
            prot_seq = trimmed_seq.translate(table=1,cds=False,to_stop=True)
            # 去除末尾的终止密码子(*)
            prot_seq_str = str(prot_seq).rstrip('*')
            
            # 创建新的序列记录
            record.seq = Seq(prot_seq_str)
            translated_records.append(record)
        except Exception as e:
            print(f"警告: 序列 {record.id} 翻译失败: {e}", file=sys.stderr)
    
    # 写入输出文件
    SeqIO.write(translated_records, pep_file, 'fasta')
    return len(translated_records)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("用法: python3 trim_and_translate.py <input.nucl.fa> <output.pep.fa>")
        sys.exit(1)
    
    nucl_file = sys.argv[1]
    pep_file = sys.argv[2]
    
    try:
        count = trim_and_translate(nucl_file, pep_file)
        print(f"✓ 成功翻译 {count} 条序列: {pep_file}")
    except Exception as e:
        print(f"✗ 翻译失败: {e}", file=sys.stderr)
        sys.exit(1)