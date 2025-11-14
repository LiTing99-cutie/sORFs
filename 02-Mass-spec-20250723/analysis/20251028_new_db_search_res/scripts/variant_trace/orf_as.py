#!/usr/bin/env python3
"""
ORF Analysis Script
分析ORF与纯合SNP的重叠关系，以及ORF的可变剪切来源（UTR延伸基于CDS起止进行判断）
"""

import gzip
import pandas as pd
import numpy as np
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Optional
import sys
import argparse
from dataclasses import dataclass
import re

# 遗传密码表（当前脚本未直接使用，保留以便扩展）
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def complement(base):
    """获取碱基的互补碱基"""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return comp.get(base.upper(), 'N')

def reverse_complement(seq):
    """获取序列的反向互补序列"""
    return ''.join(complement(base) for base in reversed(seq))

@dataclass
class SNP:
    """SNP数据结构"""
    chrom: str
    pos: int  # 1-based
    ref: str
    alt: str

@dataclass
class ORF:
    """ORF数据结构"""
    orf_id: str
    chrom: str
    start: int  # 1-based（ORF整体边界，保留）
    end: int
    strand: str
    frame: int
    pb_id: str
    cds_start: Optional[int] = None  # 新增：ORF的CDS起点（若GTF含CDS条目）
    cds_end: Optional[int] = None    # 新增：ORF的CDS终点

class ORFAnalyzer:
    def __init__(self, vcf_file, orf_gtf, annotation_gtf, classification_file):
        self.vcf_file = vcf_file
        self.orf_gtf = orf_gtf
        self.annotation_gtf = annotation_gtf
        self.classification_file = classification_file

        # 数据容器
        self.snps: List[SNP] = []
        self.orfs: List[ORF] = []
        self.annotation: Dict[str, List[Dict]] = defaultdict(list)  # gene -> list of transcript dicts
        self.classification: Dict[str, Tuple[str, str, str]] = {}   # pb_id -> (associated_gene, structural_category, subcategory)

    def parse_vcf(self):
        """解析VCF文件，提取纯合SNP"""
        print("解析VCF文件...")
        with gzip.open(self.vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')

                # 检查是否为纯合SNP (GT=1/1)
                format_field = parts[8].split(':')
                gt_index = format_field.index('GT')
                sample_data = parts[9].split(':')
                genotype = sample_data[gt_index]

                if genotype == '1/1':  # 纯合
                    ref = parts[3]
                    alt = parts[4]
                    # 只处理单碱基SNP
                    if len(ref) == 1 and len(alt) == 1:
                        self.snps.append(SNP(
                            chrom=parts[0],
                            pos=int(parts[1]),
                            ref=ref,
                            alt=alt
                        ))
        print(f"  找到 {len(self.snps)} 个纯合SNP")

    def parse_orf_gtf(self):
        """解析ORF GTF文件（同时统计ORF整体边界与CDS边界）"""
        print("解析ORF GTF文件...")
        orf_dict: Dict[str, Dict] = {}

        with open(self.orf_gtf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                feature_type = parts[2]
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attr = parts[8]

                # gene_id 作为 ORF 唯一ID
                gene_id_match = re.search(r'gene_id "([^"]+)"', attr)
                if not gene_id_match:
                    continue
                orf_full_id = gene_id_match.group(1)
                pb_id = orf_full_id.split(':')[0]

                # GTF 第8列 phase（0/1/2），没有则置0
                frame = 0
                if parts[7] != '.':
                    try:
                        frame = int(parts[7])
                    except ValueError:
                        frame = 0

                if orf_full_id not in orf_dict:
                    orf_dict[orf_full_id] = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'frame': frame,
                        'pb_id': pb_id,
                        'cds_start': None,
                        'cds_end': None
                    }
                else:
                    # 更新ORF整体边界
                    orf_dict[orf_full_id]['start'] = min(orf_dict[orf_full_id]['start'], start)
                    orf_dict[orf_full_id]['end']   = max(orf_dict[orf_full_id]['end'], end)
                    if frame != 0 and orf_dict[orf_full_id]['frame'] == 0:
                        orf_dict[orf_full_id]['frame'] = frame

                # 汇总CDS边界
                if feature_type == 'CDS':
                    if orf_dict[orf_full_id]['cds_start'] is None:
                        orf_dict[orf_full_id]['cds_start'] = start
                        orf_dict[orf_full_id]['cds_end']   = end
                    else:
                        orf_dict[orf_full_id]['cds_start'] = min(orf_dict[orf_full_id]['cds_start'], start)
                        orf_dict[orf_full_id]['cds_end']   = max(orf_dict[orf_full_id]['cds_end'], end)

        # 转为 ORF 对象
        self.orfs = []
        for orf_id, info in orf_dict.items():
            self.orfs.append(ORF(
                orf_id=orf_id,
                chrom=info['chrom'],
                start=info['start'],
                end=info['end'],
                strand=info['strand'],
                frame=info['frame'],
                pb_id=info['pb_id'],
                cds_start=info['cds_start'],
                cds_end=info['cds_end']
            ))

        print(f"  找到 {len(self.orfs)} 个ORF（去重后）")

    def parse_annotation_gtf(self):
        """解析注释GTF文件，提取基因和转录本信息"""
        print("解析注释GTF文件...")
        gene_to_transcripts: Dict[str, Set[str]] = defaultdict(set)
        transcript_info: Dict[str, Dict] = {}  # transcript_id -> info

        with open(self.annotation_gtf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                feature_type = parts[2]
                if feature_type not in ['exon', 'CDS', 'transcript']:
                    continue

                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attr = parts[8]

                gene_id_match = re.search(r'gene_id "([^"]+)"', attr)
                gene_name_match = re.search(r'gene_name "([^"]+)"', attr)
                transcript_id_match = re.search(r'transcript_id "([^"]+)"', attr)

                if gene_id_match and transcript_id_match:
                    gene_id = gene_id_match.group(1)
                    gene_name = gene_name_match.group(1) if gene_name_match else gene_id
                    transcript_id = transcript_id_match.group(1)

                    gene_to_transcripts[gene_name].add(transcript_id)

                    if transcript_id not in transcript_info:
                        transcript_info[transcript_id] = {
                            'chrom': chrom,
                            'strand': strand,
                            'transcript_start': None,
                            'transcript_end': None,
                            'exons': [],
                            'cds': []
                        }

                    if feature_type == 'transcript':
                        transcript_info[transcript_id]['transcript_start'] = start
                        transcript_info[transcript_id]['transcript_end'] = end
                    elif feature_type == 'exon':
                        transcript_info[transcript_id]['exons'].append((start, end))
                    elif feature_type == 'CDS':
                        transcript_info[transcript_id]['cds'].append((start, end))

        coding_genes = 0
        noncoding_genes = 0

        for gene, transcripts in gene_to_transcripts.items():
            gene_has_coding = False
            for transcript_id in transcripts:
                info = transcript_info[transcript_id]

                # 若缺 transcript 行，用 exon 边界兜底
                if info['transcript_start'] is None and info['exons']:
                    exon_starts = [s for s, e in info['exons']]
                    exon_ends = [e for s, e in info['exons']]
                    info['transcript_start'] = min(exon_starts)
                    info['transcript_end'] = max(exon_ends)

                transcript_data = {
                    'transcript_id': transcript_id,
                    'chrom': info['chrom'],
                    'strand': info['strand'],
                    'transcript_start': info['transcript_start'],
                    'transcript_end': info['transcript_end'],
                    'cds_start': None,
                    'cds_end': None,
                    'is_coding': False
                }

                if info['cds']:
                    cds_starts = [s for s, e in info['cds']]
                    cds_ends = [e for s, e in info['cds']]
                    transcript_data['cds_start'] = min(cds_starts)
                    transcript_data['cds_end'] = max(cds_ends)
                    transcript_data['is_coding'] = True
                    gene_has_coding = True

                self.annotation[gene].append(transcript_data)

            if gene_has_coding:
                coding_genes += 1
            else:
                noncoding_genes += 1

        print(f"  找到 {len(self.annotation)} 个基因的转录本信息")
        print(f"    其中编码基因: {coding_genes}")
        print(f"    非编码基因: {noncoding_genes}")

    def parse_classification(self):
        """解析Iso-seq分类文件"""
        print("解析Iso-seq分类文件...")
        with open(self.classification_file, 'r') as f:
            header = f.readline().strip().split('\t')
            pb_idx = header.index('isoform') if 'isoform' in header else 0
            gene_idx = header.index('associated_gene') if 'associated_gene' in header else -1
            cat_idx = header.index('structural_category') if 'structural_category' in header else -1
            subcat_idx = header.index('subcategory') if 'subcategory' in header else -1

            for line in f:
                parts = line.strip().split('\t')
                if len(parts) <= max(pb_idx, gene_idx, cat_idx, subcat_idx):
                    continue
                pb_id = parts[pb_idx]
                associated_gene = parts[gene_idx] if gene_idx >= 0 else ''
                structural_cat = parts[cat_idx] if cat_idx >= 0 else ''
                subcategory = parts[subcat_idx] if subcat_idx >= 0 else ''

                self.classification[pb_id] = (associated_gene, structural_cat, subcategory)

        print(f"  找到 {len(self.classification)} 个转录本分类信息")

    def analyze_snp_overlap(self):
        """分析ORF与SNP的重叠"""
        print("\n分析ORF与SNP重叠...")
        results: List[Dict] = []

        # 染色体分组
        snp_by_chr: Dict[str, List[SNP]] = defaultdict(list)
        for snp in self.snps:
            snp_by_chr[snp.chrom].append(snp)

        for orf in self.orfs:
            chr_snps = snp_by_chr.get(orf.chrom, [])
            for snp in chr_snps:
                if orf.start <= snp.pos <= orf.end:
                    if orf.strand == '+':
                        rel_pos = snp.pos - orf.start + 1  # 1-based
                        coding_pos = rel_pos - orf.frame
                        if coding_pos > 0:
                            codon_num = (coding_pos - 1) // 3 + 1
                            codon_pos = (coding_pos - 1) % 3 + 1
                            results.append({
                                'orf_id': orf.orf_id,
                                'pb_id': orf.pb_id,
                                'chrom': orf.chrom,
                                'snp_pos': snp.pos,
                                'ref_base': snp.ref,
                                'alt_base': snp.alt,
                                'strand': orf.strand,
                                'orf_rel_pos': rel_pos,
                                'codon_number': codon_num,
                                'codon_position': codon_pos,
                                'frame': orf.frame
                            })
                    else:
                        rel_pos = orf.end - snp.pos + 1  # 1-based
                        coding_pos = rel_pos - orf.frame
                        if coding_pos > 0:
                            codon_num = (coding_pos - 1) // 3 + 1
                            codon_pos = (coding_pos - 1) % 3 + 1
                            results.append({
                                'orf_id': orf.orf_id,
                                'pb_id': orf.pb_id,
                                'chrom': orf.chrom,
                                'snp_pos': snp.pos,
                                'ref_base': snp.ref,
                                'alt_base': snp.alt,
                                'strand': orf.strand,
                                'orf_rel_pos': rel_pos,
                                'codon_number': codon_num,
                                'codon_position': codon_pos,
                                'frame': orf.frame
                            })

        print(f"  找到 {len(results)} 个ORF-SNP重叠")
        return pd.DataFrame(results)

    def analyze_alternative_splicing(self):
        """分析ORF的可变剪切来源（UTR延伸基于CDS起止判断）"""
        print("\n分析ORF可变剪切来源...")
        results: List[Dict] = []

        for orf in self.orfs:
            pb_id = orf.pb_id

            if pb_id in self.classification:
                associated_gene, structural_cat, subcategory = self.classification[pb_id]

                as_source = 'unknown'
                evidence = ''

                # 规则1：非 FSM/ISM 视为 AS
                if structural_cat not in ['incomplete-splice_match', 'full-splice_match']:
                    as_source = 'AS'
                    evidence = f'structural_category={structural_cat}'

                # 规则2：ISM + intron retention
                elif structural_cat == 'incomplete-splice_match' and subcategory and \
                        ('intron' in subcategory.lower() and 'retention' in subcategory.lower()):
                    as_source = 'AS'
                    evidence = f'ISM with subcategory={subcategory}'

                # 规则3：FSM/ISM，基于 CDS 与转录本边界判断 UTR 延伸
                elif structural_cat in ['incomplete-splice_match', 'full-splice_match']:
                    if associated_gene in self.annotation:
                        transcripts = self.annotation[associated_gene]
                        # 同染色体同链
                        relevant_transcripts = [
                            t for t in transcripts
                            if t['chrom'] == orf.chrom and t['strand'] == orf.strand
                        ]

                        if relevant_transcripts:
                            all_transcript_starts = [
                                t['transcript_start'] for t in relevant_transcripts
                                if t['transcript_start'] is not None
                            ]
                            all_transcript_ends = [
                                t['transcript_end'] for t in relevant_transcripts
                                if t['transcript_end'] is not None
                            ]

                            if all_transcript_starts and all_transcript_ends:
                                min_transcript_start = min(all_transcript_starts)
                                max_transcript_end = max(all_transcript_ends)

                                # 使用 CDS 起止；若无 CDS 则回退至 ORF 整体边界
                                orf_cds_start = orf.cds_start if orf.cds_start is not None else orf.start
                                orf_cds_end   = orf.cds_end   if orf.cds_end   is not None else orf.end

                                # 判定延伸类型（注意负链转录方向）
                                is_5utr_extension = False
                                is_3utr_extension = False

                                if orf.strand == '+':
                                    is_5utr_extension = orf_cds_start < min_transcript_start
                                    is_3utr_extension = orf_cds_end   > max_transcript_end
                                else:
                                    is_5utr_extension = orf_cds_end   > max_transcript_end
                                    is_3utr_extension = orf_cds_start < min_transcript_start

                                if is_5utr_extension and is_3utr_extension:
                                    as_source = 'dual_extension'
                                    evidence = "CDS extends beyond both transcript boundaries"
                                elif is_5utr_extension:
                                    as_source = '5UTR_extension'
                                    evidence = "CDS extends beyond 5' transcript boundary"
                                elif is_3utr_extension:
                                    as_source = '3UTR_extension'
                                    evidence = "CDS extends beyond 3' transcript boundary"
                                else:
                                    as_source = 'non_extension'
                                    evidence = "CDS within transcript boundaries"
                            else:
                                as_source = 'unknown'
                                evidence = 'Transcript boundaries not found'
                    else:
                        as_source = 'unknown'
                        evidence = 'Associated gene not found in annotation'

                results.append({
                    'orf_id': orf.orf_id,
                    'pb_id': pb_id,
                    'chrom': orf.chrom,
                    'start': orf.start,
                    'end': orf.end,
                    'strand': orf.strand,
                    'associated_gene': associated_gene,
                    'structural_category': structural_cat,
                    'subcategory': subcategory,
                    'AS_source': as_source,
                    'evidence': evidence
                })
            else:
                results.append({
                    'orf_id': orf.orf_id,
                    'pb_id': pb_id,
                    'chrom': orf.chrom,
                    'start': orf.start,
                    'end': orf.end,
                    'strand': orf.strand,
                    'associated_gene': '',
                    'structural_category': '',
                    'subcategory': '',
                    'AS_source': 'unknown',
                    'evidence': 'PB ID not found in classification file'
                })

        print(f"  分析了 {len(results)} 个ORF的可变剪切来源")
        return pd.DataFrame(results)

    def run(self, output_prefix):
        """运行完整分析流程"""
        # 解析输入文件
        self.parse_vcf()
        self.parse_orf_gtf()
        self.parse_annotation_gtf()
        self.parse_classification()

        # 运行分析
        # snp_results = self.analyze_snp_overlap()
        as_results = self.analyze_alternative_splicing()

        # 保存结果
        # snp_output = f"{output_prefix}_snp_overlap.tsv"
        as_output = f"{output_prefix}_as_source.tsv"

        # snp_results.to_csv(snp_output, sep='\t', index=False)
        as_results.to_csv(as_output, sep='\t', index=False)

        print(f"\n结果已保存：")
        # print(f"  SNP重叠分析：{snp_output}")
        print(f"  可变剪切分析：{as_output}")

        # 输出统计摘要
        print("\n=== 分析摘要 ===")
        # if len(snp_results) > 0:
        #     print(f"SNP-ORF重叠：{len(snp_results)} 个")
        #     print(f"  涉及ORF数：{snp_results['orf_id'].nunique()}")
        #     print(f"  涉及SNP数：{snp_results['snp_pos'].nunique()}")

        if len(as_results) > 0:
            print(f"\n可变剪切来源分布：")
            as_counts = as_results['AS_source'].value_counts()
            for source, count in as_counts.items():
                print(f"  {source}: {count} ({count/len(as_results)*100:.1f}%)")

def main():
    parser = argparse.ArgumentParser(description='ORF分析：SNP重叠和可变剪切来源')
    parser.add_argument('--vcf', required=True, help='纯合SNP VCF文件 (gzipped)')
    parser.add_argument('--orf-gtf', required=True, help='ORF GTF文件')
    parser.add_argument('--annotation-gtf', required=True, help='基因组注释GTF文件')
    parser.add_argument('--classification', required=True, help='Iso-seq分类文件')
    parser.add_argument('--output-prefix', required=True, help='输出文件前缀')

    args = parser.parse_args()

    analyzer = ORFAnalyzer(
        vcf_file=args.vcf,
        orf_gtf=args.orf_gtf,
        annotation_gtf=args.annotation_gtf,
        classification_file=args.classification
    )

    analyzer.run(args.output_prefix)

if __name__ == '__main__':
    main()
