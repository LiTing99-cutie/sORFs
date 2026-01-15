#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
检查FASTA文件中的序列长度是否为3的倍数
"""

import sys
import argparse
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict

def check_fasta_file(fasta_path):
    """
    检查单个FASTA文件中的序列
    返回: (总序列数, 不是3倍数的序列列表)
    """
    problems = []
    total = 0
    
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            total += 1
            seq = str(record.seq).replace("-", "").replace("*", "")  # 去除gap和终止密码子
            seq_len = len(seq)
            
            if seq_len % 3 != 0:
                remainder = seq_len % 3
                problems.append({
                    'id': record.id,
                    'length': seq_len,
                    'remainder': remainder
                })
    except Exception as e:
        return None, str(e)
    
    return total, problems

def main():
    parser = argparse.ArgumentParser(
        description='检查FASTA文件中序列长度是否为3的倍数',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 检查单个文件
  python3 check_codon_length.py -i seq.fa
  
  # 检查目录下所有.fa文件
  python3 check_codon_length.py -d /path/to/output_root -p "*.fa"
  
  # 递归检查并输出详细信息
  python3 check_codon_length.py -d /path/to/output_root -p "*.fa" -r -v
  
  # 只显示有问题的文件
  python3 check_codon_length.py -d /path/to/output_root -p "*.fa" -r --only-problems
        """
    )
    parser.add_argument('-i', '--input', type=str, help='单个FASTA文件')
    parser.add_argument('-d', '--directory', type=str, help='要检查的目录')
    parser.add_argument('-p', '--pattern', type=str, default='*.fa', help='文件匹配模式 (默认: *.fa)')
    parser.add_argument('-r', '--recursive', action='store_true', help='递归搜索子目录')
    parser.add_argument('-v', '--verbose', action='store_true', help='显示详细信息')
    parser.add_argument('--only-problems', action='store_true', help='只显示有问题的文件')
    parser.add_argument('-o', '--output', type=str, help='输出报告到文件')
    
    args = parser.parse_args()
    
    # 收集要检查的文件
    files_to_check = []
    
    if args.input:
        files_to_check.append(Path(args.input))
    elif args.directory:
        root_dir = Path(args.directory)
        if not root_dir.exists():
            sys.stderr.write(f"错误: 目录不存在: {root_dir}\n")
            sys.exit(1)
        
        if args.recursive:
            files_to_check = list(root_dir.rglob(args.pattern))
        else:
            files_to_check = list(root_dir.glob(args.pattern))
    else:
        parser.print_help()
        sys.exit(1)
    
    if not files_to_check:
        print(f"未找到匹配的文件 (模式: {args.pattern})")
        sys.exit(0)
    
    print(f"找到 {len(files_to_check)} 个文件待检查...\n")
    
    # 统计信息
    stats = {
        'total_files': 0,
        'total_sequences': 0,
        'problem_files': 0,
        'problem_sequences': 0,
        'error_files': 0
    }
    
    results = []
    
    # 检查每个文件
    for fasta_file in sorted(files_to_check):
        stats['total_files'] += 1
        
        total_seqs, problems = check_fasta_file(fasta_file)
        
        if total_seqs is None:
            # 文件读取错误
            stats['error_files'] += 1
            if not args.only_problems or args.verbose:
                print(f"✗ 错误: {fasta_file}")
                print(f"  {problems}\n")
            results.append(f"ERROR\t{fasta_file}\t{problems}")
            continue
        
        stats['total_sequences'] += total_seqs
        
        if problems:
            stats['problem_files'] += 1
            stats['problem_sequences'] += len(problems)
            
            print(f"⚠ 问题: {fasta_file}")
            print(f"  总序列数: {total_seqs}, 有问题: {len(problems)}")
            
            if args.verbose:
                for p in problems:
                    print(f"    - {p['id']}: 长度={p['length']} (余数={p['remainder']})")
            print()
            
            results.append(f"PROBLEM\t{fasta_file}\t{total_seqs}\t{len(problems)}")
        else:
            if not args.only_problems:
                print(f"✓ 正常: {fasta_file} (序列数: {total_seqs})")
            results.append(f"OK\t{fasta_file}\t{total_seqs}\t0")
    
    # 输出统计
    print("\n" + "="*60)
    print("汇总统计:")
    print(f"  检查文件数: {stats['total_files']}")
    print(f"  总序列数:   {stats['total_sequences']}")
    print(f"  有问题文件: {stats['problem_files']}")
    print(f"  有问题序列: {stats['problem_sequences']}")
    if stats['error_files'] > 0:
        print(f"  读取错误:   {stats['error_files']}")
    print("="*60)
    
    # 保存报告
    if args.output:
        with open(args.output, 'w') as f:
            f.write("状态\t文件路径\t总序列数\t问题序列数\n")
            for line in results:
                f.write(line + "\n")
        print(f"\n报告已保存到: {args.output}")
    
    # 返回状态码
    sys.exit(0 if stats['problem_files'] == 0 else 1)

if __name__ == "__main__":
    main()