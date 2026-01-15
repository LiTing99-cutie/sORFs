#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
annotate_snp_effect.py (Fixed Version with Debug + refined stop_loss/start_gained)
一体化流程：
1) 从 VCF 读取纯合 SNP（默认首个样本，或 --sample 指定），仅保留单碱基变异；
2) 基于 ORF GTF（优先 CDS，其次 exon，再次 transcript）构建 ORF 的拼接片段（转录方向 5'->3'）；
3) 计算 SNP 与 ORF 片段的重叠（真正在拼接片段内）；
4) 对重叠记录进行跨剪接的密码子定位与功能注释：
   - synonymous / missense / nonsense / stop_loss / stop_retained / start_gained / unknown
   - 本版中：
       * start_gained：仅在 ORF 第一个密码子，非 M → M；
       * stop_loss：仅在 ORF 最后一个密码子“之外”的密码子上发生 * → 非 * 的突变。

修复：
- 负链坐标映射逻辑
- 添加详细调试输出（--debug-negative 选项）
- stop_loss/start_gained 增加首/末密码子相关逻辑
"""

import argparse
import re
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

import pandas as pd
import pysam

# ============== 遗传密码表与工具函数 ==============

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
COMPLEMENT = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

def rc_base(b: str) -> str:
    return COMPLEMENT.get(b.upper(), 'N')

def revcomp(seq: str) -> str:
    return ''.join(COMPLEMENT.get(x.upper(), 'N') for x in reversed(seq))

# ============== 解析 ORF 片段（CDS/exon/transcript） ==============

def parse_orf_segments(orf_gtf: str, debug: bool = False) -> Dict[str, Dict]:
    """
    解析 ORF GTF，返回：
    orf_id -> {
        'chrom': str, 'strand': str,
        'segments': List[Tuple[int,int]],    # 5'->3' 转录方向排序的片段（优先 CDS，其次 exon，再次 transcript）
        'span_start': int, 'span_end': int,  # ORF 片段的总体最小/最大边界（加速粗筛）
        'pb_id': str
    }
    """
    raw: Dict[str, Dict] = {}
    with open(orf_gtf, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9: continue
            chrom, _, ftype, s, e, _, strand, _, attr = p
            start, end = int(s), int(e)
            m = re.search(r'gene_id "([^"]+)"', attr)
            if not m: continue
            oid = m.group(1)
            d = raw.setdefault(oid, {'chrom':chrom,'strand':strand,'cds':[],'exon':[],'tx':[]})
            if ftype == 'CDS':
                d['cds'].append((start,end))
            elif ftype == 'exon':
                d['exon'].append((start,end))
            elif ftype == 'transcript':
                d['tx'].append((start,end))

    out: Dict[str, Dict] = {}
    for oid, info in raw.items():
        chrom, strand = info['chrom'], info['strand']
        segs = info['cds'] if info['cds'] else (info['exon'] if info['exon'] else info['tx'])
        if not segs:
            continue
        # 合并重叠/相邻片段
        segs = sorted(segs, key=lambda x: (x[0], x[1]))
        merged = []
        cur_s, cur_e = segs[0]
        for s,e in segs[1:]:
            if s <= cur_e+1:
                cur_e = max(cur_e, e)
            else:
                merged.append((cur_s,cur_e))
                cur_s,cur_e = s,e
        merged.append((cur_s,cur_e))
        
        # 转录方向排序：正链升序；负链用 5'->3'（高到低）
        if strand == '+':
            ordered = merged
        else:
            ordered = sorted(merged, key=lambda x: (-x[0], -x[1]))  # 按起始坐标降序
        
        pb_id = oid.split(':')[0]
        span_start = min(s for s,_ in ordered)
        span_end   = max(e for _,e in ordered)
        
        if debug and strand == '-' and len(ordered) > 1:
            print(f"\n[DEBUG GTF] 负链ORF: {oid}")
            print(f"  染色体: {chrom}, 链: {strand}")
            print(f"  原始segments(排序前): {sorted(merged, key=lambda x: x[0])}")
            print(f"  5'->3'排序后: {ordered}")
        
        out[oid] = {
            'chrom': chrom,
            'strand': strand,
            'segments': ordered,
            'span_start': span_start,
            'span_end': span_end,
            'pb_id': pb_id
        }
    return out

# ============== VCF 读取（pysam.VariantFile） ==============

def iter_hom_snp_from_vcf(vcf_path: str, sample_name: Optional[str]=None):
    """
    迭代 VCF 中 纯合ALT 的 单碱基 SNP：
    返回字典 {chrom,pos,ref,alt}
    - 若未指定 sample_name，默认取第一个样本；
    - 同时支持 1/1 或 1|1。
    """
    vf = pysam.VariantFile(vcf_path)
    if len(vf.header.samples) == 0:
        raise SystemExit("VCF中没有样本。")
    samp = sample_name if sample_name is not None else list(vf.header.samples)[0]
    if samp not in vf.header.samples:
        raise SystemExit(f"指定样本 {samp} 不在 VCF 中。可选：{list(vf.header.samples)}")

    for rec in vf.fetch():
        if len(rec.alts or []) != 1:  # 仅考虑单ALT（可扩展为遍历ALT）
            continue
        ref = rec.ref
        alt = rec.alts[0]
        if len(ref) != 1 or len(alt) != 1:
            continue
        gt = rec.samples[samp].get('GT', None)  # (1,1) 或 (1,1,…) 多倍体时自行扩展
        if gt is None:
            continue
        # 纯合ALT：所有等位都是 1
        if all(allele == 1 for allele in gt if allele is not None):
            yield {'chrom': rec.chrom, 'pos': rec.pos, 'ref': ref, 'alt': alt}
    vf.close()

# ============== SNP 与 ORF 片段重叠 ==============

def snp_orf_overlap(snp_list: List[Dict], orf_map: Dict[str, Dict]) -> pd.DataFrame:
    """
    输入：
      snp_list: [{'chrom','pos','ref','alt'}, ...]
      orf_map:  parse_orf_segments() 的返回
    输出：
      DataFrame: 每条重叠（SNP 在 ORF 片段内）一行
        orf_id, pb_id, chrom, strand, snp_pos, ref_base, alt_base
    """
    # 按染色体组织 ORF
    orfs_by_chr = defaultdict(list)
    for oid, meta in orf_map.items():
        orfs_by_chr[meta['chrom']].append((oid, meta))

    rows = []
    for s in snp_list:
        chrom, pos, ref, alt = s['chrom'], int(s['pos']), s['ref'], s['alt']
        if chrom not in orfs_by_chr:
            continue
        # 粗筛：先看 span 覆盖
        candidates = [ (oid,meta) for oid,meta in orfs_by_chr[chrom]
                       if meta['span_start'] <= pos <= meta['span_end'] ]
        for oid, meta in candidates:
            # 细筛：检查是否在任一片段内
            in_any = any( (seg_s <= pos <= seg_e) for seg_s,seg_e in meta['segments'] )
            if in_any:
                rows.append({
                    'orf_id': oid,
                    'pb_id': meta['pb_id'],
                    'chrom': chrom,
                    'strand': meta['strand'],
                    'snp_pos': pos,
                    'ref_base': ref,
                    'alt_base': alt
                })
    return pd.DataFrame(rows)

# ============== 拼接坐标转换与取密码子（修复版） ==============

def pos_in_segments_to_spliced_idx(pos: int, segments: List[Tuple[int,int]], strand: str, debug: bool = False) -> Optional[int]:
    """
    给定基因组坐标 pos，求其在拼接后序列中的 0-based index。
    segments 已按转录方向(5'->3')排序。
    正链：片段内部偏移 = pos - s
    负链：片段内部偏移 = e - pos
    """
    acc = 0
    for seg_idx, (s, e) in enumerate(segments):
        if s <= pos <= e:
            if strand == '+':
                idx = acc + (pos - s)
            else:
                idx = acc + (e - pos)
            
            if debug and strand == '-':
                print(f"    [pos_to_idx] segment[{seg_idx}]=({s},{e}), pos={pos}, offset={e-pos}, spliced_idx={idx}")
            
            return idx
        acc += (e - s + 1)
    return None

def spliced_idx_to_genome_triplet(codon_start_idx: int, segments: List[Tuple[int,int]], strand: str, debug: bool = False) -> Optional[List[int]]:
    """
    给定拼接序列中某密码子的起始 0-based 索引，返回三联体对应的 3 个基因组坐标（1-based）。
    可能跨剪接点。
    
    修复：负链坐标映射
    正链：genome = s + offset
    负链：genome = e - offset（offset是在该segment内从5'端的偏移）
    """
    coords = []
    want = [codon_start_idx, codon_start_idx+1, codon_start_idx+2]
    acc = 0
    
    if debug and strand == '-':
        print(f"    [idx_to_triplet] codon_start_idx={codon_start_idx}, want={want}")
    
    for seg_idx, (s, e) in enumerate(segments):
        seg_len = e - s + 1
        for i, wi in enumerate(list(want)):
            if wi is None:
                continue
            if acc <= wi < acc + seg_len:
                offset = wi - acc
                if strand == '+':
                    genome_pos = s + offset
                else:
                    # 负链：segments已按5'->3'排序（从高到低）
                    # offset=0 对应该segment的5'端（最高位置=e）
                    genome_pos = e - offset
                
                coords.append(genome_pos)
                
                if debug and strand == '-':
                    print(f"      segment[{seg_idx}]=({s},{e}), want_idx={wi}, offset={offset}, genome={genome_pos}")
                
                want[i] = None
        acc += seg_len
    
    if len(coords) != 3:
        return None
    return coords

def fetch_coding_codon(ref: pysam.FastaFile, chrom: str, strand: str, gpos_triplet: List[int], debug: bool = False) -> str:
    """基于基因组三个位点，取参考三联体（正链取原序列，负链反向互补为编码方向）。
    
    关键：对于负链，gpos_triplet是按编码顺序（5'->3'）排列的，即从高到低的基因组坐标。
    我们需要按基因组顺序（从低到高）取碱基，然后整体反向互补。
    """
    if strand == '+':
        # 正链：按给定顺序取碱基
        bases = ''.join(ref.fetch(chrom, p-1, p).upper() for p in gpos_triplet)
        codon = bases
    else:
        # 负链：先按基因组顺序（升序）排序，取碱基，再反向互补
        sorted_pos = sorted(gpos_triplet)
        bases = ''.join(ref.fetch(chrom, p-1, p).upper() for p in sorted_pos)
        codon = revcomp(bases)
        
        if debug:
            print(f"    [fetch_codon] 编码顺序positions={gpos_triplet}")
            print(f"    [fetch_codon] 基因组顺序positions={sorted_pos}")
            print(f"    [fetch_codon] 基因组碱基(正链)={bases}, 反向互补(编码)={codon}")
    
    return codon

def annotate_row(row, orf_map, ref: pysam.FastaFile, debug_negative: bool = False) -> Dict[str,str]:
    """
    对单条 SNP-ORF 重叠记录做注释。row 必含：orf_id, chrom, snp_pos, ref_base, alt_base

    本版：
    - start_gained：仅当 alt_aa 为 'M' 且 ref_aa 非 'M'，且该密码子是 ORF 的第一个密码子；
    - stop_loss：仅当 ref_aa 为 '*' 且 alt_aa 非 '*'，且该密码子不是 ORF 的最后一个密码子。
    """
    out = {
        'ref_codon':'NNN','alt_codon':'NNN','ref_aa':'X','alt_aa':'X',
        'mutation_type':'unknown','aa_change':'NA',
        'computed_codon_number': 'NA', 'computed_codon_position':'NA', 'reason':''
    }
    oid = row['orf_id']
    if oid not in orf_map:
        out['reason'] = 'orf_not_in_gtf'
        return out
    chrom = row['chrom']
    snp_pos = int(row['snp_pos'])
    ref_b = str(row['ref_base']).upper()
    alt_b = str(row['alt_base']).upper()
    strand = orf_map[oid]['strand']
    segs = orf_map[oid]['segments']

    # 计算 ORF 总长度（nt）及总密码子数，用于判断首/末密码子
    total_nt = sum(e - s + 1 for s, e in segs)
    total_codons = total_nt // 3

    debug = debug_negative and (strand == '-')
    
    if debug:
        print(f"\n{'='*60}")
        print(f"[DEBUG 负链注释]")
        print(f"  ORF ID: {oid}")
        print(f"  染色体: {chrom}, 链: {strand}")
        print(f"  Segments (5'->3'): {segs}")
        print(f"  SNP位置: {snp_pos}")
        print(f"  VCF REF: {ref_b}, ALT: {alt_b}")

    # SNP 是否仍在片段内（极少数清洗后不一致）
    idx = pos_in_segments_to_spliced_idx(snp_pos, segs, strand, debug=debug)
    if idx is None:
        out['reason'] = 'snp_outside_orf_segments'
        return out

    if debug:
        print(f"  拼接序列索引: {idx}")

    codon_start_idx = (idx // 3) * 3
    codon_pos = (idx % 3) + 1
    codon_number = (codon_start_idx // 3) + 1  # 第几个密码子
    
    if debug:
        print(f"  密码子起始索引: {codon_start_idx}, 密码子内位置: {codon_pos}, 第 {codon_number} 个密码子")
    
    gtriplet = spliced_idx_to_genome_triplet(codon_start_idx, segs, strand, debug=debug)
    if gtriplet is None:
        out['reason'] = 'triplet_cross_failed'
        return out

    if debug:
        print(f"  基因组三联体坐标: {gtriplet}")

    ref_codon = fetch_coding_codon(ref, chrom, strand, gtriplet, debug=debug)
    
    # 找到SNP位置在密码子中的索引
    if strand == '+':
        snp_codon_idx = gtriplet.index(snp_pos)
    else:
        # 负链：SNP在gtriplet中的位置就是编码顺序中的位置
        snp_codon_idx = gtriplet.index(snp_pos)
    
    coding_ref_base = ref_codon[snp_codon_idx]
    expected_ref = ref_b if strand == '+' else rc_base(ref_b)
    
    if debug:
        print(f"  参考密码子(编码方向): {ref_codon}")
        print(f"  SNP在密码子中的索引: {snp_codon_idx} (密码子位置{codon_pos})")
        print(f"  密码子该位置的碱基: {coding_ref_base}")
        print(f"  VCF REF: {ref_b}, 编码方向期望: {expected_ref}")
        print(f"  匹配: {coding_ref_base == expected_ref}")
    
    if coding_ref_base != expected_ref:
        out['reason'] = f'coding_ref_mismatch:{coding_ref_base}!={expected_ref}'
        if debug:
            print(f"  ❌ 不匹配!")
            # 额外调试：检查每个位置
            print(f"  详细检查三个位置(编码顺序):")
            for i, gp in enumerate(gtriplet):
                actual = ref.fetch(chrom, gp-1, gp).upper()
                actual_coding = actual if strand == '+' else rc_base(actual)
                print(f"    编码位置{i+1}: genome_pos={gp}, 正链base={actual}, 编码base={actual_coding}")

    # 构造 ALT 密码子（编码方向）
    alt_codon_list = list(ref_codon)
    alt_codon_list[snp_codon_idx] = alt_b if strand == '+' else rc_base(alt_b)
    alt_codon = ''.join(alt_codon_list)

    ref_aa = CODON_TABLE.get(ref_codon, 'X')
    alt_aa = CODON_TABLE.get(alt_codon, 'X')

    # 判型（加入“首/末密码子”相关逻辑）
    if ref_aa == 'X' or alt_aa == 'X':
        mtype = 'unknown'
    elif ref_aa == alt_aa:
        mtype = 'synonymous'
    elif ref_aa == '*' or alt_aa == '*':
        is_last_codon = (total_codons > 0 and codon_number == total_codons)
        if ref_aa == '*' and alt_aa != '*':
            # 仅“非最后一个密码子”的 * → 非 * 记为 stop_loss，最后一个密码子记为 missense
            if not is_last_codon:
                mtype = 'stop_loss'
            else:
                mtype = 'missense'
        elif ref_aa != '*' and alt_aa == '*':
            # 获得终止密码子，仍然按传统意义标为 nonsense（不要求是最后一个）
            mtype = 'nonsense'
        else:
            # * -> *，保持终止密码子
            mtype = 'stop_retained'
    elif alt_aa == 'M' and ref_aa != 'M':
        # 只有第一个密码子才算 start_gained
        if codon_number == 1:
            mtype = 'start_gained'
        else:
            mtype = 'missense'
    else:
        mtype = 'missense'

    if debug:
        print(f"  ALT密码子: {alt_codon}")
        print(f"  REF氨基酸: {ref_aa}, ALT氨基酸: {alt_aa}")
        print(f"  突变类型: {mtype}")
        print(f"{'='*60}\n")

    out.update({
        'ref_codon': ref_codon,
        'alt_codon': alt_codon,
        'ref_aa': ref_aa,
        'alt_aa': alt_aa,
        'mutation_type': mtype,
        'aa_change': f'{ref_aa}{codon_number}{alt_aa}',
        'computed_codon_number': codon_number,
        'computed_codon_position': codon_pos,
    })
    return out

# ============== 主流程 ==============

def main():
    ap = argparse.ArgumentParser(description='SNP 功能注释（修复负链+调试版）')
    ap.add_argument('--vcf', required=True, help='VCF/BCF（建议索引），将筛选纯合ALT的单碱基SNP')
    ap.add_argument('--orf-gtf', required=True, help='ORF GTF（优先CDS，其次exon/transcript）')
    ap.add_argument('--ref-genome', required=True, help='参考基因组 FASTA（需 .fai）')
    ap.add_argument('--output', required=True, help='最终带注释的TSV')
    ap.add_argument('--sample', default=None, help='VCF样本名（默认取首个样本）')
    ap.add_argument('--write-overlap', default=None, help='可选：输出中间的 snp_overlap TSV')
    ap.add_argument('--debug-negative', action='store_true', help='开启负链详细调试输出')
    ap.add_argument('--debug-limit', type=int, default=5, help='调试模式下只处理前N条负链记录（默认5）')
    args = ap.parse_args()

    # 1) 解析 ORF 片段
    print('[1/4] 解析 ORF 片段 ...')
    orf_map = parse_orf_segments(args.orf_gtf, debug=args.debug_negative)
    print(f'  ORFs: {len(orf_map)}')
    neg_orfs = sum(1 for v in orf_map.values() if v['strand'] == '-')
    print(f'  负链ORFs: {neg_orfs}')

    # 2) 从 VCF 读取纯合 SNP
    print('[2/4] 读取 VCF 纯合ALT 单碱基SNP ...')
    snps = list(iter_hom_snp_from_vcf(args.vcf, args.sample))
    print(f'  SNPs: {len(snps)}')

    # 3) 计算 SNP–ORF 重叠（片段内）
    print('[3/4] 计算 SNP–ORF 重叠 ...')
    ov_df = snp_orf_overlap(snps, orf_map)
    print(f'  overlaps: {len(ov_df)}  (distinct SNP: {ov_df["snp_pos"].nunique() if not ov_df.empty else 0}, ORF: {ov_df["orf_id"].nunique() if not ov_df.empty else 0})')
    
    if not ov_df.empty:
        neg_overlaps = len(ov_df[ov_df['strand'] == '-'])
        print(f'  负链overlaps: {neg_overlaps}')
    
    if args.write_overlap is not None:
        ov_df.to_csv(args.write_overlap, sep='\t', index=False)
        print(f'  已写出中间重叠表: {args.write_overlap}')

    if ov_df.empty:
        print('无重叠记录，写空表并退出。')
        ov_df.to_csv(args.output, sep='\t', index=False)
        return

    # 4) 对重叠进行注释
    print('[4/4] 对重叠进行突变类型注释 ...')
    ref = pysam.FastaFile(args.ref_genome)
    
    # 如果开启调试模式，只处理部分负链记录
    if args.debug_negative:
        neg_mask = ov_df['strand'] == '-'
        neg_df = ov_df[neg_mask].head(args.debug_limit)
        pos_df = ov_df[~neg_mask]
        
        print(f"\n{'='*60}")
        print(f"调试模式：处理前 {len(neg_df)} 条负链记录")
        print(f"{'='*60}")
        
        annos_neg = [annotate_row(r, orf_map, ref, debug_negative=True) for _, r in neg_df.iterrows()]
        annos_pos = [annotate_row(r, orf_map, ref, debug_negative=False) for _, r in pos_df.iterrows()]
        
        anno_df = pd.DataFrame(annos_neg + annos_pos)
        out = pd.concat([pd.concat([neg_df, pos_df]).reset_index(drop=True), anno_df], axis=1)
    else:
        annos = [annotate_row(r, orf_map, ref, debug_negative=False) for _, r in ov_df.iterrows()]
        anno_df = pd.DataFrame(annos)
        out = pd.concat([ov_df.reset_index(drop=True), anno_df], axis=1)
    
    ref.close()

    out.to_csv(args.output, sep='\t', index=False)

    print('\n=== mutation_type counts ===')
    print(out['mutation_type'].value_counts(dropna=False))
    
    # 统计负链的mismatch情况
    if not out.empty:
        neg_out = out[out['strand'] == '-']
        if len(neg_out) > 0:
            print(f'\n=== 负链统计 ===')
            print(f'总负链记录: {len(neg_out)}')
            mismatch = neg_out['reason'].str.contains('coding_ref_mismatch', na=False).sum()
            print(f'coding_ref_mismatch: {mismatch} ({mismatch/len(neg_out)*100:.1f}%)')

    print(f'\n结果已保存到: {args.output}')

if __name__ == '__main__':
    main()
