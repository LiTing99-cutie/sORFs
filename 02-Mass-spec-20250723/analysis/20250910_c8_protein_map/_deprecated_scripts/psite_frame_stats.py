#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, argparse
from collections import defaultdict, namedtuple
from bisect import bisect_right

try:
    import pysam
except ImportError:
    sys.exit("ERROR: 请先安装 pysam：pip install pysam")

Segment = namedtuple("Segment", "start end tid strand cds_offset")  # [start,end)

def parse_args():
    ap = argparse.ArgumentParser(
        description="统计每个ORF的CDS上P-site的frame0/1/2、frame0_fraction与frame0覆盖的codon数"
    )
    ap.add_argument("--bam", required=True, help="P-site已校正的BAM（建议1bp对齐）")
    ap.add_argument("--annot", required=True, help="注释文件（GTF或genepred）")
    ap.add_argument("--format", choices=["gtf","genepred"], required=True, help="注释文件格式")
    ap.add_argument("--key-attr", default="gene_id",
                    help="GTF里ORF/转录本ID的属性名（默认gene_id）")
    ap.add_argument("--out", required=True, help="输出TSV路径")
    return ap.parse_args()

def parse_gtf_as_cds(annot_path, key_attr="gene_id"):
    cds_by_tid = defaultdict(list)  # tid -> list of (chrom, start, end, strand)
    with open(annot_path) as f:
        for line in f:
            if not line or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != "CDS": continue
            # 解析属性
            tid = None
            for kv in attrs.split(";"):
                kv = kv.strip()
                if not kv: continue
                if kv.startswith(key_attr):
                    # 形如 gene_id "XXX"
                    s = kv.split(None, 1)[-1].strip()
                    if s.startswith('"') and s.endswith('"'):
                        s = s[1:-1]
                    tid = s
                    break
            if tid is None: continue
            s = int(start) - 1  # GTF为1-based闭区间，这里转为0-based半开
            e = int(end)
            cds_by_tid[tid].append((chrom, s, e, strand))
    return cds_by_tid

def parse_genepred(annot_path):
    """
    期望标准genepred 12列：name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, ...
    这里把 [cdsStart, cdsEnd) 与外显子求交得到CDS外显子段。
    """
    cds_by_tid = defaultdict(list)
    with open(annot_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11: continue
            name = parts[0]
            chrom = parts[1]
            strand = parts[2]
            txStart = int(parts[3]); txEnd = int(parts[4])
            cdsStart = int(parts[5]); cdsEnd = int(parts[6])
            exonCount = int(parts[7])
            exonStarts = [int(x) for x in parts[8].rstrip(",").split(",")]
            exonEnds   = [int(x) for x in parts[9].rstrip(",").split(",")]

            if cdsEnd <= cdsStart:  # 非编码
                continue

            # 与CDS求交，得到CDS外显子块
            for es, ee in zip(exonStarts, exonEnds):
                s = max(es, cdsStart)
                e = min(ee, cdsEnd)
                if e > s:
                    cds_by_tid[name].append((chrom, s, e, strand))
    return cds_by_tid

def build_segments_with_offsets(cds_by_tid):
    """
    为每个tid的CDS外显子，按转录本方向排序并计算cds_offset。
    返回：by_chrom = {chrom: [Segment... 按start排序]}
    以及：tid2len = {tid: CDS长度（nt）}
    """
    by_chrom = defaultdict(list)
    tid2len = {}
    for tid, blocks in cds_by_tid.items():
        if not blocks: continue
        # 确定链向
        strands = {b[3] for b in blocks}
        if len(strands) != 1:
            # 非典型情况，跳过
            continue
        strand = list(strands)[0]

        # 按转录本方向排序
        if strand == "+":
            blocks_sorted = sorted(blocks, key=lambda x: x[1])
        else:
            blocks_sorted = sorted(blocks, key=lambda x: x[1], reverse=True)

        # 计算cds_offset
        offset = 0
        segs = []
        total_len = 0
        for chrom, s, e, _ in blocks_sorted:
            segs.append((chrom, s, e, strand, offset))
            l = e - s
            offset += l
            total_len += l
        tid2len[tid] = total_len

        # 存入按染色体索引的结构
        for chrom, s, e, strand, off in segs:
            by_chrom[chrom].append(Segment(s, e, tid, strand, off))

    # 对每条染色体的段按start排序，供二分检索
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda seg: seg.start)
    return by_chrom, tid2len

def find_overlaps(by_chrom, chrom, pos):
    """
    给定染色体与单碱基pos（0-based），返回覆盖该pos的所有Segment（可能多个ORF重叠）
    """
    if chrom not in by_chrom:
        return []
    arr = by_chrom[chrom]
    starts = [seg.start for seg in arr]
    i = bisect_right(starts, pos)
    hits = []
    # 向左回退检查可能覆盖的段
    j = i - 1
    while j >= 0 and arr[j].start <= pos:
        if pos < arr[j].end:
            hits.append(arr[j])
        j -= 1
    return hits

def main():
    args = parse_args()

    # 解析注释
    if args.format == "gtf":
        cds_by_tid = parse_gtf_as_cds(args.annot, key_attr=args.key_attr)
    else:
        cds_by_tid = parse_genepred(args.annot)

    by_chrom, tid2len = build_segments_with_offsets(cds_by_tid)

    # 结果容器
    f0 = defaultdict(int); f1 = defaultdict(int); f2 = defaultdict(int)
    f0_codons = defaultdict(set)  # tid -> set of codon_idx

    bam = pysam.AlignmentFile(args.bam, "rb")
    for rd in bam.fetch(until_eof=True):
        if rd.is_unmapped: 
            continue
        # 假设每条read的参考位置已经代表P-site（常见做法：1bp对齐）
        chrom = bam.get_reference_name(rd.reference_id)
        pos = rd.reference_start  # 0-based 左端点

        # 映射到所有覆盖该pos的CDS段（可能多个ORF重叠）
        for seg in find_overlaps(by_chrom, chrom, pos):
            # 计算该pos在CDS中的offset（按转录本方向）
            if seg.strand == "+":
                offset = seg.cds_offset + (pos - seg.start)
            else:
                # 负链：该段方向相反，段内相对偏移要反向计算
                offset = seg.cds_offset + ((seg.end - 1) - pos)

            if offset < 0: 
                continue
            # 防御：若超过该tid的长度则跳过
            total = tid2len.get(seg.tid, 0)
            if not total or offset >= total:
                continue

            frame = offset % 3
            codon_idx = offset // 3
            if frame == 0:
                f0[seg.tid] += 1
                f0_codons[seg.tid].add(codon_idx)
            elif frame == 1:
                f1[seg.tid] += 1
            else:
                f2[seg.tid] += 1

    bam.close()

    # 输出
    with open(args.out, "w") as w:
        w.write("\t".join([
            "ORF_id","CDS_nt_len","n_codons",
            "frame0","frame1","frame2",
            "total_psites","frame0_fraction","frame0_codon_covered"
        ]) + "\n")
        all_tids = set(tid2len.keys()) | set(f0.keys()) | set(f1.keys()) | set(f2.keys())
        for tid in sorted(all_tids):
            L = tid2len.get(tid, 0)
            ncod = L // 3
            a = f0.get(tid, 0); b = f1.get(tid, 0); c = f2.get(tid, 0)
            tot = a + b + c
            frac = (a / tot) if tot > 0 else 0.0
            cov_cod = len(f0_codons.get(tid, set()))
            w.write("\t".join(map(str, [
                tid, L, ncod, a, b, c, tot, f"{frac:.6f}", cov_cod
            ])) + "\n")

if __name__ == "__main__":
    main()
