#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, argparse, time
from collections import defaultdict, namedtuple

import pysam

Segment = namedtuple("Segment", "chrom start end strand cds_offset tid")  # [start,end)

def parse_args():
    ap = argparse.ArgumentParser(
        description="单次遍历BAM统计每个ORF在CDS上的frame0/1/2、frame0_fraction与frame0覆盖codon数（P-site已1bp对齐）"
    )
    ap.add_argument("--bam", required=True, help="P-site已校正BAM（建议1bp对齐）")
    ap.add_argument("--annot", required=True, help="注释文件（GTF或genepred）")
    ap.add_argument("--format", choices=["gtf","genepred"], required=True)
    ap.add_argument("--key-attr", default="gene_id", help="GTF里ORF ID的属性名（默认gene_id）")
    ap.add_argument("--out", required=True, help="输出TSV路径")
    ap.add_argument("--threads", type=int, default=4, help="BGZF解压线程数（仅影响读取速度）")
    ap.add_argument("--log-every", type=int, default=2_000_000, help="每处理N个pileup位点打印一次进度")
    return ap.parse_args()

def parse_gtf_as_cds(annot_path, key_attr="gene_id"):
    cds = defaultdict(list)  # tid -> [(chrom,start,end,strand)]
    with open(annot_path) as f:
        for ln in f:
            if not ln or ln.startswith("#"): continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS": continue
            chrom, start, end, strand, attrs = p[0], int(p[3])-1, int(p[4]), p[6], p[8]
            tid = None
            for kv in attrs.split(";"):
                kv = kv.strip()
                if not kv: continue
                if kv.startswith(key_attr):
                    s = kv.split(None,1)[-1].strip()
                    if s.startswith('"') and s.endswith('"'): s = s[1:-1]
                    tid = s; break
            if tid is None: continue
            cds[tid].append((chrom, start, end, strand))
    return cds

def parse_genepred(annot_path):
    cds = defaultdict(list)
    with open(annot_path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"): continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 11: continue
            name, chrom, strand = p[0], p[1], p[2]
            cdsStart, cdsEnd = int(p[5]), int(p[6])
            exonStarts = [int(x) for x in p[8].rstrip(",").split(",") if x]
            exonEnds   = [int(x) for x in p[9].rstrip(",").split(",") if x]
            if cdsEnd <= cdsStart: continue
            for es, ee in zip(exonStarts, exonEnds):
                s = max(es, cdsStart); e = min(ee, cdsEnd)
                if e > s: cds[name].append((chrom, s, e, strand))
    return cds

def build_segments(cds_by_tid):
    """构建：1) 每个ORF的CDS总长；2) 染色体->起点排序的Segment列表；3) 染色体的[min,max)扫描边界"""
    orf2len = {}
    chrom2segs = defaultdict(list)
    chrom_min = {}
    chrom_max = {}
    for tid, blocks in cds_by_tid.items():
        if not blocks: continue
        strands = {b[3] for b in blocks}
        if len(strands) != 1:
            continue
        strand = blocks[0][3]
        blocks_sorted = sorted(blocks, key=lambda x: x[1], reverse=(strand=="-"))
        offset = 0
        total = 0
        for chrom, s, e, _ in blocks_sorted:
            chrom2segs[chrom].append(Segment(chrom, s, e, strand, offset, tid))
            l = e - s
            offset += l
            total += l
            chrom_min[chrom] = s if chrom not in chrom_min else min(chrom_min[chrom], s)
            chrom_max[chrom] = e if chrom not in chrom_max else max(chrom_max[chrom], e)
        orf2len[tid] = total
    # 每条染色体按start排序
    for chrom in chrom2segs:
        chrom2segs[chrom].sort(key=lambda x: x.start)
    return chrom2segs, orf2len, chrom_min, chrom_max

def main():
    a = parse_args()

    # 读注释
    cds_by_tid = parse_gtf_as_cds(a.annot, a.key_attr) if a.format=="gtf" else parse_genepred(a.annot)
    chrom2segs, orf2len, chrom_min, chrom_max = build_segments(cds_by_tid)

    # 统计容器
    f0 = defaultdict(int); f1 = defaultdict(int); f2 = defaultdict(int)
    cov = {tid: bytearray(max(0, L // 3)) for tid, L in orf2len.items()}  # frame0覆盖的codon位

    # 打开BAM（并行解压仅影响IO）
    bam = pysam.AlignmentFile(a.bam, "rb", threads=a.threads)

    total_cols = 0
    t0 = time.time()

    # 按染色体、限定[min,max)范围做一次pileup扫描（无限深度：max_depth=0）
    for chrom, segs in chrom2segs.items():
        if not segs: continue
        start = chrom_min[chrom]; end = chrom_max[chrom]
        active = []
        s_idx = 0  # segs已按start排序

        for col in bam.pileup(
            chrom, start, end,
            truncate=True,
            stepper="all",
            min_base_quality=0,
            ignore_overlaps=False,
            ignore_orphans=False,
            max_depth=0  # 0 表示不限深度
        ):
            pos = col.reference_pos
            if pos < start or pos >= end:
                continue

            # 扩充/维护 active 段
            while s_idx < len(segs) and segs[s_idx].start <= pos:
                active.append(segs[s_idx]); s_idx += 1
            if active and active[0].end <= pos:
                active = [sg for sg in active if sg.end > pos]
            if not active:
                continue

            # 本位点有效read数（不做MAPQ/FLAG过滤，仅跳过删除/跳跃）
            nreads = 0
            for pr in col.pileups:
                aln = pr.alignment
                if aln.is_unmapped: continue
                if pr.is_del or pr.is_refskip: continue
                nreads += 1
            if nreads == 0:
                continue

            # 把该位点的计数分配到覆盖它的所有段
            for seg in active:
                if not (seg.start <= pos < seg.end):
                    continue
                L = orf2len.get(seg.tid, 0)
                if L <= 0:
                    continue
                offset = seg.cds_offset + ((pos - seg.start) if seg.strand=="+" else (seg.end - 1 - pos))
                if offset < 0 or offset >= L:
                    continue
                frame = offset % 3
                if frame == 0:
                    f0[seg.tid] += nreads
                    ci = offset // 3
                    if ci < len(cov[seg.tid]):
                        cov[seg.tid][ci] = 1
                elif frame == 1:
                    f1[seg.tid] += nreads
                else:
                    f2[seg.tid] += nreads

            total_cols += 1
            if a.log_every and (total_cols % a.log_every == 0):
                elapsed = (time.time() - t0) / 60
                sys.stderr.write(f"[{chrom}] processed {total_cols} pileup positions, elapsed {elapsed:.1f} min\n")

    bam.close()

    # 输出
    with open(a.out, "w") as w:
        w.write("\t".join([
            "ORF_id","CDS_nt_len","n_codons",
            "frame0","frame1","frame2",
            "total_psites","frame0_fraction","frame0_codon_covered"
        ]) + "\n")
        for tid, L in orf2len.items():
            ncod = L // 3
            a0 = f0.get(tid, 0); a1 = f1.get(tid, 0); a2 = f2.get(tid, 0)
            tot = a0 + a1 + a2
            frac = (a0 / tot) if tot > 0 else 0.0
            cov_cnt = int(sum(cov[tid])) if ncod > 0 else 0
            w.write("\t".join(map(str, [
                tid, L, ncod, a0, a1, a2, tot, f"{frac:.6f}", cov_cnt
            ])) + "\n")

if __name__ == "__main__":
    main()
