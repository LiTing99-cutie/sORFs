#!/usr/bin/env python3
"""
分析SNP突变是否被质谱肽段覆盖
"""

import pandas as pd
import sys
from typing import Dict, List, Tuple, Set


def parse_orf_id(orf_id: str) -> Tuple[str, int, int]:
    """
    解析ORF ID获取基本信息
    例如: PB.22.22:chr1:-|23|1332:342:834|noncoding|ATG
    返回: (pb_id, orf_start, orf_end)
    """
    parts = orf_id.split('|')
    pb_id = parts[0].split(':')[0]
    
    # 获取ORF坐标
    coords = parts[2].split(':')
    orf_start = int(coords[1])
    orf_end = int(coords[2])
    
    return pb_id, orf_start, orf_end


def get_peptide_coverage_positions(peptide: str, orf_seq: str) -> List[Tuple[int, int]]:
    """
    在ORF序列中找到肽段的所有匹配位置
    返回: [(start_aa_pos, end_aa_pos), ...]  (1-based氨基酸位置)
    """
    positions = []
    start = 0
    
    while True:
        pos = orf_seq.find(peptide, start)
        if pos == -1:
            break
        # 转换为1-based氨基酸位置
        aa_start = pos + 1
        aa_end = pos + len(peptide)
        positions.append((aa_start, aa_end))
        start = pos + 1
    
    return positions


def check_mutation_coverage(mutation_codon_pos: int, 
                           peptide_positions: List[Tuple[int, int]]) -> bool:
    """
    检查突变位置是否在任何肽段覆盖范围内
    mutation_codon_pos: 1-based氨基酸位置
    peptide_positions: [(start, end), ...] 肽段覆盖范围
    """
    for start, end in peptide_positions:
        if start <= mutation_codon_pos <= end:
            return True
    return False


def analyze_coverage(snp_file: str, peptide_file: str, orf_file: str, output_file: str):
    """
    主分析函数
    """
    # 读取文件
    print("读取文件...")
    snp_df = pd.read_csv(snp_file, sep='\t')
    peptide_df = pd.read_csv(peptide_file, sep='\t')
    orf_df = pd.read_csv(orf_file, sep='\t')
    
    # 创建ORF ID到序列的映射
    print("构建ORF序列索引...")
    orf_seq_dict = dict(zip(orf_df['ORF_id'], orf_df['ORF_seq']))
    
    # 创建ORF ID到肽段的映射
    print("构建肽段索引...")
    orf_peptides = {}
    for _, row in peptide_df.iterrows():
        orf_id = row['Assigned_protein_id']
        peptide = row['Peptide']
        sample = row['Sample']
        
        if orf_id not in orf_peptides:
            orf_peptides[orf_id] = []
        orf_peptides[orf_id].append({
            'peptide': peptide,
            'sample': sample
        })
    
    # 分析每个突变
    print("分析突变覆盖情况...")
    results = []
    
    for idx, snp_row in snp_df.iterrows():
        orf_id = snp_row['orf_id']
        mutation_pos = snp_row['computed_codon_number']
        mutation_type = snp_row['mutation_type']
        aa_change = snp_row['aa_change']
        
        # 获取ORF序列
        if orf_id not in orf_seq_dict:
            results.append({
                'orf_id': orf_id,
                'mutation_position': mutation_pos,
                'mutation_type': mutation_type,
                'aa_change': aa_change,
                'orf_found': False,
                'has_peptides': False,
                'is_covered': False,
                'covering_peptides': '',
                'covering_samples': '',
                'peptide_positions': ''
            })
            continue
        
        orf_seq = orf_seq_dict[orf_id]
        
        # 获取该ORF的肽段
        peptides_info = orf_peptides.get(orf_id, [])
        
        if not peptides_info:
            results.append({
                'orf_id': orf_id,
                'mutation_position': mutation_pos,
                'mutation_type': mutation_type,
                'aa_change': aa_change,
                'orf_found': True,
                'has_peptides': False,
                'is_covered': False,
                'covering_peptides': '',
                'covering_samples': '',
                'peptide_positions': ''
            })
            continue
        
        # 检查每个肽段是否覆盖突变位点
        covering_peptides = []
        covering_samples = []
        peptide_position_strs = []
        is_covered = False
        
        for pep_info in peptides_info:
            peptide = pep_info['peptide']
            sample = pep_info['sample']
            
            # 找到肽段在ORF序列中的位置
            positions = get_peptide_coverage_positions(peptide, orf_seq)
            
            if not positions:
                continue
            
            # 检查是否覆盖突变位点
            if check_mutation_coverage(mutation_pos, positions):
                is_covered = True
                covering_peptides.append(peptide)
                covering_samples.append(sample)
                pos_str = ';'.join([f"{s}-{e}" for s, e in positions])
                peptide_position_strs.append(pos_str)
        
        results.append({
            'orf_id': orf_id,
            'mutation_position': mutation_pos,
            'mutation_type': mutation_type,
            'aa_change': aa_change,
            'ref_aa': snp_row.get('ref_aa', ''),
            'alt_aa': snp_row.get('alt_aa', ''),
            'orf_found': True,
            'has_peptides': True,
            'is_covered': is_covered,
            'covering_peptides': '|'.join(covering_peptides),
            'covering_samples': '|'.join(covering_samples),
            'peptide_positions': '|'.join(peptide_position_strs),
            'total_peptides_for_orf': len(peptides_info)
        })
    
    # 创建结果DataFrame
    result_df = pd.DataFrame(results)
    
    # 保存结果
    print(f"保存结果到 {output_file}...")
    result_df.to_csv(output_file, sep='\t', index=False)
    
    # 打印统计信息
    print("\n" + "="*80)
    print("分析统计:")
    print("="*80)
    print(f"总突变数: {len(result_df)}")
    print(f"找到ORF序列的突变数: {result_df['orf_found'].sum()}")
    print(f"有肽段的突变数: {result_df['has_peptides'].sum()}")
    print(f"被肽段覆盖的突变数: {result_df['is_covered'].sum()}")
    
    if len(result_df) > 0:
        covered_pct = (result_df['is_covered'].sum() / len(result_df)) * 100
        print(f"覆盖率: {covered_pct:.2f}%")
    
    # 按突变类型统计
    if 'mutation_type' in result_df.columns:
        print("\n按突变类型统计:")
        for mut_type in result_df['mutation_type'].unique():
            subset = result_df[result_df['mutation_type'] == mut_type]
            covered = subset['is_covered'].sum()
            total = len(subset)
            print(f"  {mut_type}: {covered}/{total} ({covered/total*100:.1f}%)")
    
    print("\n详细结果已保存到:", output_file)
    
    return result_df


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python analyze_mutation_coverage.py <snp_file> <peptide_file> <orf_file> <output_file>")
        print("\n示例:")
        print("  python analyze_mutation_coverage.py \\")
        print("    snp_annotated.tsv \\")
        print("    all_samples.assignments.unique.post.tsv \\")
        print("    orfs_merged_final.tsv \\")
        print("    mutation_coverage_results.tsv")
        sys.exit(1)
    
    snp_file = sys.argv[1]
    peptide_file = sys.argv[2]
    orf_file = sys.argv[3]
    output_file = sys.argv[4]
    
    try:
        result_df = analyze_coverage(snp_file, peptide_file, orf_file, output_file)
        print("\n分析完成!")
    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)