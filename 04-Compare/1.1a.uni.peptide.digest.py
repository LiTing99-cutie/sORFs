#!/usr/bin/env python3
import time
import argparse
from pyteomics import fasta, parser

def digest_fasta(input_file, output_file, missed_cleavages=2, min_len=6, max_len=50):
    """使用Trypsin消化FASTA文件"""
    start_time = time.perf_counter()
    
    with open(output_file, "w") as f_out:
        for header, seq in fasta.read(input_file):
            peptides = parser.cleave(
                seq, 
                'trypsin', 
                missed_cleavages=missed_cleavages,
                min_length=min_len,
                max_length=max_len
            )
            for i, pep in enumerate(peptides):
                f_out.write(f">{header}_peptide{i+1}\n{pep}\n")

    elapsed = time.perf_counter() - start_time
    print(f"实际运行时间: {elapsed:.2f}秒")

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(
        description="使用Trypsin消化FASTA文件，生成肽段",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True, 
                        help="输入FASTA文件路径")
    parser.add_argument("-o", "--output", required=True,
                        help="输出FASTA文件路径")
    parser.add_argument("-m", "--missed_cleavages", type=int, default=2,
                        help="允许的漏切位点数")
    parser.add_argument("--min_len", type=int, default=6,
                        help="最小肽段长度")
    parser.add_argument("--max_len", type=int, default=50,
                        help="最大肽段长度")
    
    args = parser.parse_args()
    
    # 运行消化函数
    digest_fasta(
        input_file=args.input,
        output_file=args.output,
        missed_cleavages=args.missed_cleavages,
        min_len=args.min_len,
        max_len=args.max_len
    )

if __name__ == "__main__":
    main()