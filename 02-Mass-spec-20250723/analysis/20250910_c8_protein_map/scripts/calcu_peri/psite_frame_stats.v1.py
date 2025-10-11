#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, argparse, time
from collections import defaultdict, namedtuple

import pysam

Segment = namedtuple("Segment", "chrom start end strand cds_offset")  # [start,end)

def parse_args():
    ap = argparse.ArgumentParser(
        description="逐个ORF统计CDS上P-site的frame0/1/2、frame0_fraction与frame0覆盖codon数"
    )
    ap.add_argument("--bam", required=True, help="P-site已校正BAM（建议1bp对齐）")
    ap.add_argument("--annot", required=True, help="注释文件（GTF或genepred）")
    ap.add_argument("--format", choices=["gtf","genepred"], required=True)
    ap.add_argument("--key-attr", default="gene_id", help="GTF里ORF ID的属性名（默认gene_id）")
    ap.add_argument("--out", required=True, help="输出TSV路径")
    ap.add_argument("--log-every", type=int, default=2000, help="每处理N个ORF打印一次进度")
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

def build_orf_segments(cds_by_tid):
    """为每个ORF构建按转录本方向排序并带cds_offset的段列表"""
    orf2segs = {}         # tid -> [Segment,...]（转录本方向排序）
    tid2len = {}          # tid -> CDS总长nt
    for tid, blocks in cds_by_tid.items():
        if not blocks: continue
        strands = {b[3] for b in blocks}
        if len(strands) != 1:  # 异常跳过
            continue
        strand = list(strands)[0]
        blocks_sorted = sorted(blocks, key=lambda x: x[1], reverse=(strand=="-"))
        offset = 0
        segs = []
        total = 0
        for chrom, s, e, _ in blocks_sorted:
            segs.append(Segment(chrom, s, e, strand, offset))
            l = e - s
            offset += l
            total += l
        orf2segs[tid] = segs
        tid2len[tid] = total
    return orf2segs, tid2len

def main():
    a = parse_args()
    # 读注释
    cds_by_tid = parse_gtf_as_cds(a.annot, a.key_attr) if a.format=="gtf" else parse_genepred(a.annot)
    orf2segs, tid2len = build_orf_segments(cds_by_tid)

    bam = pysam.AlignmentFile(a.bam, "rb")

    f0 = defaultdict(int); f1 = defaultdict(int); f2 = defaultdict(int)
    f0_codons = defaultdict(set)

    # 逐ORF计算
    t0 = time.time()
    for i, (tid, segs) in enumerate(orf2segs.items(), 1):
        L = tid2len.get(tid, 0)
        if L <= 0: continue
        for seg in segs:
            # 仅在该段范围内fetch，避免重复计数（外显子彼此不重叠）
            for rd in bam.fetch(seg.chrom, seg.start, seg.end):
                if rd.is_unmapped: 
                    continue
                pos = rd.reference_start  # P-site位置（已校正为read起点的话）
                if not (seg.start <= pos < seg.end):
                    continue  # fetch可能带出边界，需要再判一次
                # 段内偏移 -> CDS偏移
                if seg.strand == "+":
                    offset = seg.cds_offset + (pos - seg.start)
                else:
                    offset = seg.cds_offset + ((seg.end - 1) - pos)
                if offset < 0 or offset >= L:
                    continue
                frame = offset % 3
                codon_idx = offset // 3
                if frame == 0:
                    f0[tid] += 1
                    f0_codons[tid].add(codon_idx)
                elif frame == 1:
                    f1[tid] += 1
                else:
                    f2[tid] += 1

        if a.log_every and i % a.log_every == 0:
            elapsed = time.time() - t0
            sys.stderr.write(f"[{i}/{len(orf2segs)} ORFs] elapsed {elapsed/60:.1f} min\n")

    bam.close()

    # 输出
    with open(a.out, "w") as w:
        w.write("\t".join([
            "ORF_id","CDS_nt_len","n_codons",
            "frame0","frame1","frame2",
            "total_psites","frame0_fraction","frame0_codon_covered"
        ]) + "\n")
        for tid in sorted(orf2segs.keys()):
            L = tid2len.get(tid, 0)
            ncod = L // 3
            a0 = f0.get(tid, 0); a1 = f1.get(tid, 0); a2 = f2.get(tid, 0)
            tot = a0 + a1 + a2
            frac = (a0 / tot) if tot > 0 else 0.0
            cov = len(f0_codons.get(tid, set()))
            w.write("\t".join(map(str, [tid, L, ncod, a0, a1, a2, tot, f"{frac:.6f}", cov])) + "\n")

if __name__ == "__main__":
    main()
