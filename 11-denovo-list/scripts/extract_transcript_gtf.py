#!/usr/bin/env python3
"""
完整提取denovo ORF的转录本GTF

步骤：
1. 从map文件中提取在denovo_list中的ORF（第二列匹配）
2. Part I: 从map文件第四列提取同基因的同源ORF转录本
3. Part II: 提取不在map第二列的ORF转录本
4. 合并所有部分并输出GTF
"""

import sys
import re
from collections import defaultdict
from pathlib import Path

def extract_pb_id(orf_id):
    """从ORF ID中提取PB转录本ID
    例如: PB.40813.6:chr6:-|127|... -> PB.40813.6
    """
    return orf_id.split(':')[0]

def parse_gtf(gtf_file):
    """
    解析GTF文件，构建：
    1. transcript_id -> gene_id 映射
    2. transcript_id -> GTF行列表
    """
    transcript_to_gene = {}
    transcript_gtf_lines = defaultdict(list)
    
    print(f"  正在解析 {gtf_file}...")
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            attributes = fields[8]
            
            # 提取 gene_id 和 transcript_id
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            transcript_id_match = re.search(r'transcript_id "([^"]+)"', attributes)
            
            if gene_id_match and transcript_id_match:
                gene_id = gene_id_match.group(1)
                transcript_id = transcript_id_match.group(1)
                
                transcript_to_gene[transcript_id] = gene_id
                transcript_gtf_lines[transcript_id].append(line)
    
    return transcript_to_gene, transcript_gtf_lines

def extract_map_matches_with_homologs(denovo_ids, map_file):
    """从map文件中提取第二列在denovo_ids中的行，并解析第四列的同源ORF"""
    matched_lines = []
    orf_data = []
    
    with open(map_file) as f:
        next(f)  # 跳过表头
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                query_orf = fields[1]
                if query_orf in denovo_ids:
                    matched_lines.append(line)
                    
                    # 解析第四列的同源ORF（如果存在）
                    if len(fields) >= 4:
                        homologs_str = fields[3]
                        # 跳过 "No ATG" 等特殊值
                        if not homologs_str.startswith('No '):
                            # 解析同源ORF列表（逗号+空格分隔）
                            homolog_orfs = [h.strip() for h in homologs_str.split(',') if h.strip()]
                            
                            orf_data.append({
                                'query': query_orf,
                                'homologs': homolog_orfs
                            })
    
    return matched_lines, orf_data

def main():
    if len(sys.argv) != 5:
        print("用法: python extract_complete_denovo_gtf.py <denovo_list_file> <map_file> <gtf_file> <output_dir>")
        print("\n参数:")
        print("  denovo_list_file : denovo ORF列表文件（第一列为ORF ID，无表头）")
        print("  map_file         : representative.tsv 文件（第2列=ORF ID，第4列=同源ORF列表，有表头）")
        print("  gtf_file         : Iso-seq GTF文件")
        print("  output_dir       : 输出目录")
        print("\n示例:")
        print("  python extract_complete_denovo_gtf.py \\")
        print("    ../processed/denovo_list.txt \\")
        print("    /path/to/representative.tsv \\")
        print("    /path/to/custom.gtf.with_orf.gtf \\")
        print("    ../processed")
        sys.exit(1)
    
    denovo_list_file = sys.argv[1]
    map_file = sys.argv[2]
    gtf_file = sys.argv[3]
    output_dir = sys.argv[4]
    
    # 创建输出目录
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 输出文件路径
    denovo_id_file = output_dir / "denovo_list.id.txt"
    map_matched_file = output_dir / "denovo_list.in_map.txt"
    output_gtf = output_dir / "denovo_orfs.transcripts.complete.gtf"
    
    print("=" * 70)
    print("完整提取denovo ORF转录本GTF")
    print("=" * 70)
    print(f"\n输入文件:")
    print(f"  Denovo list  : {denovo_list_file}")
    print(f"  Map file     : {map_file}")
    print(f"  GTF file     : {gtf_file}")
    print(f"\n输出目录: {output_dir}")
    print()
    
    # ===== Step 0: 提取denovo ID列表 =====
    print("Step 0: 提取denovo ID列表...")
    denovo_ids = set()
    with open(denovo_list_file) as f:
        # denovo_list.txt 没有表头，直接读取
        for line in f:
            fields = line.strip().split('\t')
            if fields and fields[0]:  # 确保第一列不为空
                denovo_ids.add(fields[0])
    
    # 保存ID列表
    # with open(denovo_id_file, 'w') as f:
    #     for orf_id in sorted(denovo_ids):
    #         f.write(orf_id + '\n')
    
    print(f"  - Denovo ID总数: {len(denovo_ids)}")
    # print(f"  - 保存到: {denovo_id_file}")
    
    # ===== Step 1: 从map中提取匹配的ORF及其同源信息 =====
    print("\nStep 1: 从map文件中提取匹配的ORF及同源信息...")
    matched_lines, orf_data = extract_map_matches_with_homologs(denovo_ids, map_file)
    
    # 保存匹配的行（包含表头）
    with open(map_matched_file, 'w') as f:
        # 先写入表头
        with open(map_file) as mf:
            header = mf.readline()
            f.write(header)
        # 写入匹配的行
        for line in matched_lines:
            f.write(line)
    
    # 提取匹配的ORF ID集合（第二列）
    matched_orf_ids = set()
    for line in matched_lines:
        fields = line.strip().split('\t')
        if len(fields) >= 2:
            matched_orf_ids.add(fields[1])
    
    print(f"  - 在map中找到: {len(matched_orf_ids)} 个ORF")
    print(f"  - 其中有同源信息的: {len(orf_data)} 个ORF")
    print(f"  - 保存到: {map_matched_file}")
    
    # ===== Step 2: 解析GTF文件 =====
    print("\nStep 2: 解析GTF文件...")
    transcript_to_gene, transcript_gtf_lines = parse_gtf(gtf_file)
    print(f"  - 找到 {len(transcript_to_gene)} 个转录本")
    print(f"  - 找到 {len(set(transcript_to_gene.values()))} 个基因")
    
    # ===== Step 3: Part I - 处理同基因同源ORF =====
    print("\nStep 3: Part I - 提取同基因同源ORF转录本...")
    
    part1_transcripts = set()
    stats_part1 = {
        'total_groups': 0,
        'same_gene_groups': 0,
    }
    
    for orf_group in orf_data:
        stats_part1['total_groups'] += 1
        query_orf = orf_group['query']
        homolog_orfs = orf_group['homologs']
        
        # 提取所有相关的转录本ID
        query_pb = extract_pb_id(query_orf)
        homolog_pbs = [extract_pb_id(h) for h in homolog_orfs]
        
        # 获取query转录本的基因ID
        if query_pb not in transcript_to_gene:
            continue
        
        query_gene = transcript_to_gene[query_pb]
        
        # 检查同源ORF的转录本是否属于同一基因
        same_gene_transcripts = []
        for homolog_pb in homolog_pbs:
            if homolog_pb in transcript_to_gene:
                if transcript_to_gene[homolog_pb] == query_gene:
                    same_gene_transcripts.append(homolog_pb)
        
        # 如果有同基因的同源转录本，提取它们
        if same_gene_transcripts:
            stats_part1['same_gene_groups'] += 1
            same_gene_transcripts.append(query_pb)
            part1_transcripts.update(same_gene_transcripts)
    
    print(f"  - 处理的ORF组: {stats_part1['total_groups']}")
    print(f"  - 同基因组数: {stats_part1['same_gene_groups']}")
    print(f"  - Part I转录本数: {len(part1_transcripts)}")
    
    # ===== Step 4: Part II - 提取不在map第二列的ORF =====
    print("\nStep 4: Part II - 提取不在map中的ORF转录本...")
    
    # 读取map文件第二列的所有ID（跳过表头）
    map_ids = set()
    with open(map_file) as f:
        next(f)  # 跳过表头
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                map_ids.add(fields[1])
    
    print(f"  - Map中的ID总数: {len(map_ids)}")
    
    # 找出不在map中的ID
    unmapped_ids = denovo_ids - map_ids
    print(f"  - 不在map中的ID: {len(unmapped_ids)}")
    
    # 提取这些ID对应的转录本
    part2_transcripts = set()
    for orf_id in unmapped_ids:
        pb_id = extract_pb_id(orf_id)
        if pb_id in transcript_gtf_lines:
            part2_transcripts.add(pb_id)
    
    print(f"  - Part II转录本数: {len(part2_transcripts)}")
    
    # ===== Step 5: 合并所有部分 =====
    print("\nStep 5: 合并所有部分...")
    all_transcripts = part1_transcripts | part2_transcripts
    overlap = part1_transcripts & part2_transcripts
    
    print(f"  - Part I (同基因同源): {len(part1_transcripts)} 个转录本")
    print(f"  - Part II (不在map中): {len(part2_transcripts)} 个转录本")
    print(f"  - 重叠: {len(overlap)} 个转录本")
    print(f"  - 合并后总数: {len(all_transcripts)} 个转录本")
    
    # ===== Step 6: 写入输出GTF =====
    print("\nStep 6: 写入输出GTF...")
    
    with open(output_gtf, 'w') as out:
        # 写入注释信息
        out.write(f"# Denovo ORF transcripts GTF\n")
        out.write(f"# Total denovo ORFs: {len(denovo_ids)}\n")
        out.write(f"# In map (representative): {len(matched_orf_ids)}\n")
        out.write(f"# Not in map: {len(unmapped_ids)}\n")
        out.write(f"# Part I (same-gene homologs): {len(part1_transcripts)} transcripts\n")
        out.write(f"# Part II (unmapped ORFs): {len(part2_transcripts)} transcripts\n")
        out.write(f"# Total transcripts: {len(all_transcripts)}\n")
        out.write("#\n")
        
        for transcript_id in sorted(all_transcripts):
            if transcript_id in transcript_gtf_lines:
                for line in transcript_gtf_lines[transcript_id]:
                    out.write(line)
    
    print(f"  - 输出文件: {output_gtf}")
    
    print("\n" + "=" * 70)
    print("✓ 完成!")
    print("=" * 70)
    
    # 输出详细统计
    print("\n最终统计:")
    print(f"  输入:")
    print(f"    - Denovo ORF总数: {len(denovo_ids)}")
    print(f"    - 在map中的ORF: {len(matched_orf_ids)}")
    print(f"    - 不在map中的ORF: {len(unmapped_ids)}")
    print(f"  输出:")
    print(f"    - Part I 转录本: {len(part1_transcripts)}")
    print(f"    - Part II 转录本: {len(part2_transcripts)}")
    print(f"    - 总转录本数: {len(all_transcripts)}")
    print(f"\n  输出文件:")
    print(f"    - {denovo_id_file}")
    print(f"    - {map_matched_file}")
    print(f"    - {output_gtf}")

if __name__ == "__main__":
    main()