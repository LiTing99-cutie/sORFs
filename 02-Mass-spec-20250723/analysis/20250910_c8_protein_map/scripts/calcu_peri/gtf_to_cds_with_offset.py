#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, argparse, re, csv
from collections import defaultdict, namedtuple
Seg = namedtuple("Seg", "chrom start end tid strand cds_offset")

def parse_args():
    ap = argparse.ArgumentParser("GTF/GenePred → 带 cds_offset 的 CDS 片段BED")
    ap.add_argument("--annot", required=True)
    ap.add_argument("--format", choices=["gtf","genepred"], required=True)
    ap.add_argument("--out-bed", required=True, help="输出：cds_with_offset.bed")
    ap.add_argument("--out-map", required=True, help="输出：tid_map.tsv（transcript_id gene_id cds_len）")
    return ap.parse_args()

def parse_gtf(p):
    attr_re = re.compile(r'(\S+)\s+"([^"]+)"')
    cds_by_tid = defaultdict(list); tid2gene = {}
    with open(p) as f:
        for ln in f:
            if not ln or ln[0]=='#': continue
            a = ln.rstrip("\n").split("\t")
            if len(a)<9 or a[2]!="CDS": continue
            chrom, start, end, strand, attrs = a[0], int(a[3])-1, int(a[4]), a[6], a[8]
            attr = dict(attr_re.findall(attrs))
            tid = attr.get("transcript_id") or attr.get("transcriptId") or attr.get("transcript")
            gid = attr.get("gene_id") or attr.get("geneId") or attr.get("gene") or ""
            if not tid: continue
            if gid: tid2gene[tid]=gid
            cds_by_tid[tid].append((chrom,start,end,strand))
    segs=[]; tid2len={}
    for tid, xs in cds_by_tid.items():
        if not xs: continue
        strand = xs[0][3]
        xs.sort(key=lambda z:z[1], reverse=(strand=='-'))
        off=0
        for chrom,s,e,_ in xs:
            segs.append(Seg(chrom,s,e,tid,strand,off))
            off += (e-s)
        tid2len[tid]=off
    return segs, tid2gene, tid2len

def parse_genepred(p):
    segs=[]; tid2gene={}; tid2len={}
    with open(p) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"): continue
            a = ln.rstrip("\n").split("\t")
            if len(a)<10: continue
            tid, chrom, strand = a[0], a[1], a[2]
            cdsStart, cdsEnd = int(a[5]), int(a[6])
            if cdsEnd<=cdsStart: continue
            exonStarts = [int(x) for x in a[8].rstrip(",").split(",")]
            exonEnds   = [int(x) for x in a[9].rstrip(",").split(",")]
            chunks=[]
            for s,e in zip(exonStarts,exonEnds):
                ss=max(s,cdsStart); ee=min(e,cdsEnd)
                if ee>ss: chunks.append((chrom,ss,ee,strand))
            if not chunks: continue
            chunks.sort(key=lambda z:z[1], reverse=(strand=='-'))
            off=0
            for chrom,s,e,_ in chunks:
                segs.append(Seg(chrom,s,e,tid,strand,off))
                off += (e-s)
            tid2len[tid]=off; tid2gene[tid]=tid
    return segs, tid2gene, tid2len

def main():
    args=parse_args()
    if args.format=="gtf":
        segs, tid2gene, tid2len = parse_gtf(args.annot)
    else:
        segs, tid2gene, tid2len = parse_genepred(args.annot)

    # 写 BED：chr start end tid 0 strand cds_offset  （7列）
    with open(args.out_bed,"w") as w:
        for s in segs:
            w.write(f"{s.chrom}\t{s.start}\t{s.end}\t{s.tid}\t0\t{s.strand}\t{s.cds_offset}\n")
    # 写映射：transcript_id  gene_id  cds_len
    with open(args.out_map,"w",newline="") as w:
        tsv=csv.writer(w,delimiter="\t")
        tsv.writerow(["transcript_id","gene_id","cds_len"])
        for tid, L in tid2len.items():
            gid = tid2gene.get(tid,"")
            tsv.writerow([tid,gid,L])

if __name__=="__main__":
    main()
