#!/usr/bin/env python3
"""
筛选候选ORF的脚本
筛选条件:
1. frame0_fraction >= 0.5
2. Psites_codon_coverage >= 0.1
3. inframe_fraction == 0
"""

import pandas as pd
from Bio import SeqIO
import argparse
import sys

def load_data(psite_frame_file, orf_rpf_file, orf_overlap_file):
    """加载三个TSV文件"""
    print("Loading data files...")
    
    # 读取P-site frame统计文件
    df_frame = pd.read_csv(psite_frame_file, sep='\t')
    print(f"Loaded {len(df_frame)} records from {psite_frame_file}")
    
    # 读取RPF P-site文件
    df_rpf = pd.read_csv(orf_rpf_file, sep='\t')
    print(f"Loaded {len(df_rpf)} records from {orf_rpf_file}")
    
    # 读取overlap文件
    df_overlap = pd.read_csv(orf_overlap_file, sep='\t')
    print(f"Loaded {len(df_overlap)} records from {orf_overlap_file}")
    
    return df_frame, df_rpf, df_overlap

def merge_data(df_frame, df_rpf, df_overlap):
    """合并三个数据框"""
    print("\nMerging data...")
    
    # 基于ORF_id合并
    df_merged = df_frame.merge(df_rpf, on='ORF_id', how='inner')
    df_merged = df_merged.merge(df_overlap, on='ORF_id', how='inner')
    
    print(f"Merged data contains {len(df_merged)} records")
    return df_merged

def filter_orfs(df_merged, frame0_threshold=0.5, psite_cov_threshold=0.1, inframe_threshold=0):
    """根据条件筛选ORF"""
    print("\nApplying filters...")
    print(f"Criteria:")
    print(f"  - frame0_fraction >= {frame0_threshold}")
    print(f"  - Psites_codon_coverage >= {psite_cov_threshold}")
    print(f"  - inframe_fraction == {inframe_threshold}")
    
    # 应用筛选条件
    df_filtered = df_merged[
        (df_merged['frame0_fraction'] >= frame0_threshold) &
        (df_merged['Psites_codon_coverage'] >= psite_cov_threshold) &
        (df_merged['inframe_fraction'] == inframe_threshold)
    ]
    
    print(f"\nFiltered results: {len(df_filtered)} ORFs passed the filters")
    print(f"Filtering rate: {len(df_filtered)/len(df_merged)*100:.2f}%")
    
    return df_filtered

def filter_fasta(input_fasta, output_fasta, valid_ids):
    """从FASTA文件中提取符合条件的序列"""
    print(f"\nFiltering FASTA file: {input_fasta}")
    
    valid_ids_set = set(valid_ids)
    filtered_count = 0
    
    with open(output_fasta, 'w') as out_handle:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            if record.id in valid_ids_set:
                SeqIO.write(record, out_handle, 'fasta')
                filtered_count += 1
    
    print(f"Wrote {filtered_count} sequences to {output_fasta}")
    
    # 检查是否所有筛选出的ID都在FASTA文件中找到
    if filtered_count < len(valid_ids):
        print(f"Warning: {len(valid_ids) - filtered_count} ORF IDs were not found in the FASTA file")

def main():
    parser = argparse.ArgumentParser(
        description='Filter candidate ORFs based on frame, P-site coverage, and overlap criteria'
    )
    parser.add_argument(
        '--frame_stats',
        default='/home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/peri/psite_frame_stats.v2.tsv',
        help='Path to P-site frame statistics file'
    )
    parser.add_argument(
        '--rpf_psite',
        default='/home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt',
        help='Path to ORF RPF P-site file'
    )
    parser.add_argument(
        '--overlap',
        default='/home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf_overlap_inframe.txt',
        help='Path to ORF overlap file'
    )
    parser.add_argument(
        '--input_fasta',
        default='/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.pep.fa',
        help='Path to input candidate ORF FASTA file'
    )
    parser.add_argument(
        '--output_fasta',
        default='filtered_candidateORF.fa',
        help='Path to output filtered FASTA file'
    )
    parser.add_argument(
        '--output_table',
        default='filtered_ORF_table.tsv',
        help='Path to output filtered table (optional)'
    )
    parser.add_argument(
        '--frame0_threshold',
        type=float,
        default=0.5,
        help='Minimum frame0_fraction (default: 0.5)'
    )
    parser.add_argument(
        '--psite_cov_threshold',
        type=float,
        default=0.1,
        help='Minimum Psites_codon_coverage (default: 0.1)'
    )
    parser.add_argument(
        '--inframe_threshold',
        type=float,
        default=0.0,
        help='Required inframe_fraction value (default: 0.0)'
    )
    
    args = parser.parse_args()
    
    try:
        # 加载数据
        df_frame, df_rpf, df_overlap = load_data(
            args.frame_stats,
            args.rpf_psite,
            args.overlap
        )
        
        # 合并数据
        df_merged = merge_data(df_frame, df_rpf, df_overlap)
        
        # 筛选ORF
        df_filtered = filter_orfs(
            df_merged,
            args.frame0_threshold,
            args.psite_cov_threshold,
            args.inframe_threshold
        )
        
        # 保存筛选结果表格
        if len(df_filtered) > 0:
            df_filtered.to_csv(args.output_table, sep='\t', index=False)
            print(f"\nSaved filtered table to {args.output_table}")
            
            # 筛选FASTA文件
            filter_fasta(args.input_fasta, args.output_fasta, df_filtered['ORF_id'].tolist())
        else:
            print("\nNo ORFs passed the filters!")
            sys.exit(1)
        
        print("\nDone!")
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()