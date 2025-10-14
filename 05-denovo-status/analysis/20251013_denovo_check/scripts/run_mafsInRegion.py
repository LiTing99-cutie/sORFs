#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, argparse, subprocess
from pathlib import Path
from collections import namedtuple, defaultdict
from Bio import AlignIO
from Bio.Seq import Seq

Exon = namedtuple("Exon", ["chrom", "start1b", "end1b", "strand", "bed_name", "order"])

def ensure_dir(p: Path): p.mkdir(parents=True, exist_ok=True)
def sanitize_filename(s: str) -> str: return re.sub(r"[^A-Za-z0-9._-]+", "__", s)

def read_orf_beds_1based(per_chr_bed_dir: Path):
    orf_exons = {}
    chrom_name = per_chr_bed_dir.name  # e.g. chr1
    for bed in sorted(per_chr_bed_dir.glob("*.bed")):
        rows = []
        with bed.open() as f:
            for line in f:
                if not line.strip() or line.startswith("#"): continue
                ps = line.rstrip("\n").split("\t")
                chrom = ps[0] if ps[0].startswith("chr") else "chr"+ps[0]
                s1, e1, name, strand = int(ps[1]), int(ps[2]), ps[3], ps[5]
                rows.append((chrom, s1, e1, strand, name))
        if not rows: continue
        rows.sort(key=lambda x: x[1])
        orf_id = rows[0][4]
        exons = [Exon(r[0], r[1], r[2], r[3], None, i+1) for i, r in enumerate(rows)]
        orf_exons[orf_id] = exons
    if not orf_exons:
        raise SystemExit(f"[ERROR] {per_chr_bed_dir} 未找到 BED")
    return orf_exons, chrom_name

def build_combined_bed(orf_exons, chrom_name, combined_bed: Path, name_map_path: Path):
    ensure_dir(combined_bed.parent); ensure_dir(name_map_path.parent)
    name_map = {}
    with combined_bed.open("w") as bedout, name_map_path.open("w") as mapout:
        for orf_id, exons in orf_exons.items():
            safe_orf = sanitize_filename(orf_id)
            for ex in exons:
                start0, end0 = ex.start1b - 1, ex.end1b
                bed_name = f"{safe_orf}__exon{ex.order}"
                bedout.write(f"{chrom_name}\t{start0}\t{end0}\t{bed_name}\n")
                name_map[bed_name] = (orf_id, ex.order, ex.strand)
                mapout.write(f"{bed_name}\t{orf_id}\t{ex.order}\t{ex.strand}\n")
    return name_map

def run_mafsInRegion(mafsInRegion_bin: str, combined_bed: Path, out_dir: Path, maf_gz: Path):
    ensure_dir(out_dir)
    cmd = [
        mafsInRegion_bin,
        str(combined_bed),     # regions.bed（0-based）
        str(out_dir),          # 输出目录
        str(maf_gz),           # 输入 chrN.maf.gz
        "-outDir"              # 关键开关：将第二个参数解释为目录
        # 若需要保留初始gaps，可再加 "-keepInitialGaps"
    ]
    print(f"[INFO] mafsInRegion: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def exon_maf_path(out_dir: Path, bed_name: str):
    p1 = out_dir / bed_name
    p2 = out_dir / f"{bed_name}.maf"
    if p1.exists(): return p1
    if p2.exists(): return p2
    for p in out_dir.glob(f"{bed_name}*"): return p
    return None

def extract_sequences_from_block(ma, pin_1b, pif_1b, regin0, focal_prefix):
    all_seqs = {}
    for s in ma:
        sp = s.id.split(".")[0]
        chrom = s.id.split(".", 1)[1] if "." in s.id else s.id
        all_seqs[sp] = ["", 0, 0, chrom, int(s.annotations.get("start",0)), int(s.annotations.get("strand",1))]
    p1 = p2 = -1
    for s in ma:
        if not s.id.startswith(focal_prefix + "."): continue
        counter = regin0
        ss = str(s.seq)
        for n, ch in enumerate(ss):
            if (counter >= pin_1b-1) and (counter <= pif_1b-1):
                if p1 == -1: p1 = n
            elif (counter > pif_1b-1) and (p2 == -1):
                p2 = n-1
            if ch != "-": counter += 1
        if p2 == -1: p2 = len(ss)-1
        break
    if p1 == -1: return None
    for s in ma:
        sp = s.id.split(".")[0]
        text = str(s.seq)
        p11 = p22 = 0
        for n, ch in enumerate(text):
            if (n < p1) and (ch != "-"): p11 += 1
            if (n < p2) and (ch != "-"): p22 += 1
        subseq = text[p1:p2+1]
        chrom  = s.id.split(".", 1)[1] if "." in s.id else s.id
        all_seqs[sp] = [subseq, p11, p22, chrom, int(s.annotations.get("start",0)), int(s.annotations.get("strand",1))]
    return all_seqs

def stitch_orf(orf_id, exons, exon_maf_dir: Path, focal: str, cov_thresh: float,
               out_orfs_dir: Path, out_beds_dir: Path):
    orf_seq = defaultdict(lambda: ["", [], [], [], []])  # seq, s_list, e_list, chr_list, strand_list
    L = 0
    coords = []
    for ex in exons:
        L += (ex.end1b - ex.start1b + 1)
        maf_file = exon_maf_path(exon_maf_dir, ex.bed_name)
        if maf_file is None: continue
        with maf_file.open() as fh:
            for ma in AlignIO.parse(fh, "maf"):
                st0 = int(ma[0].annotations.get("start",0))
                usize = int(ma[0].annotations.get("size",0))
                key = f"{st0}/{usize}"
                if key in coords: continue
                all_seqs = extract_sequences_from_block(ma, ex.start1b, ex.end1b, st0, focal)
                if all_seqs is None: continue
                coords.append(key)
                for sp, vals in all_seqs.items():
                    orf_seq[sp][0] += vals[0]
                    orf_seq[sp][1].append(vals[4] + vals[1])
                    orf_seq[sp][2].append(vals[4] + vals[2])
                    orf_seq[sp][3].append(vals[3])
                    orf_seq[sp][4].append(vals[5])

    if focal not in orf_seq or not orf_seq[focal][0]:
        return "MISS"

    safe_orf = sanitize_filename(orf_id)
    fa_out  = out_orfs_dir / f"{safe_orf}.fa"
    maf_out = out_orfs_dir / f"{safe_orf}.maf"
    bed_out = out_beds_dir / f"{safe_orf}.bed"
    ensure_dir(out_orfs_dir); ensure_dir(out_beds_dir)

    m = float(len(orf_seq[focal][0].replace("-","")) / float(L))
    focal_len = len(orf_seq[focal][0])

    with fa_out.open("w") as f1, maf_out.open("w") as f2, bed_out.open("w") as f3:
        for sp, vals in orf_seq.items():
            seq_aln = vals[0]
            if exons[0].strand == "-":
                seq_aln = str(Seq(seq_aln).reverse_complement())
            if focal_len == 0 or (len(seq_aln)/focal_len) < cov_thresh:
                continue
            f1.write(f">{orf_id}_{orf_id}_{sp}_{m:.6f}\n{str(Seq(seq_aln.replace('-','')).translate(cds=False))}\n")
            f2.write(f">{sp}\n{seq_aln}\n")
            if sp != focal:
                t2 = "_splitted" if len(set(vals[3])) != 1 else ""
                for n in range(len(vals[1])):
                    strand_char = "-" if str(vals[4][n]) in ("-1","-") else "+"
                    f3.write(f"{vals[3][n]}\t{vals[1][n]}\t{vals[2][n]}\t{orf_id}\t{sp}{t2}\t{strand_char}\n")
    return "OK"

def main():
    ap = argparse.ArgumentParser(description="单线程：mafsInRegion 一次切片 + 逐 ORF 拼接导出")
    ap.add_argument("--per-chr-bed-dir", required=True, help="该染色体 ORF BED 目录（1-based）")
    ap.add_argument("--maf-root", required=True, help="包含 chrN.maf.gz 的目录")
    ap.add_argument("--output-root", required=True, help="输出根目录")
    ap.add_argument("--work-dir", required=True, help="工作目录（合并BED与mafsInRegion输出）")
    ap.add_argument("--chr", default="", help="染色体名（留空则用目录名）")
    ap.add_argument("--focal", default="hg38", help="焦点物种前缀（默认 hg38）")
    ap.add_argument("--cov-thresh", type=float, default=0.95, help="覆盖阈值（默认0.95）")
    ap.add_argument("--mafsInRegion", default="/home/user/data/lit/Tools/ucsc/bin/mafsInRegion", help="mafsInRegion 路径")
    args = ap.parse_args()

    per_chr_bed_dir = Path(args.per_chr_bed_dir)
    orf_exons, chr_detect = read_orf_beds_1based(per_chr_bed_dir)
    chr_name = args.chr or chr_detect

    # 合并 BED（提交给 mafsInRegion 用 0-based；name 唯一且可回溯）
    tmp_chr_dir = Path(args.work_dir) / chr_name
    ensure_dir(tmp_chr_dir)
    combined_bed = tmp_chr_dir / f"{chr_name}.all_orfs.0based.bed"
    name_map_tsv = tmp_chr_dir / f"{chr_name}.name_map.tsv"
    build_combined_bed(orf_exons, chr_name, combined_bed, name_map_tsv)

    # 运行 mafsInRegion（一次）
    maf_gz = Path(args.maf_root) / f"{chr_name}.maf.gz"
    if not maf_gz.exists():
        raise SystemExit(f"[ERROR] 找不到 {maf_gz}")
    exon_out_dir = tmp_chr_dir / "mafsInRegion_exons"
    run_mafsInRegion(args.mafsInRegion, combined_bed, exon_out_dir, maf_gz)

    # 把 bed_name 填回各 exon
    for orf_id, exons in orf_exons.items():
        safe_orf = sanitize_filename(orf_id)
        for i, ex in enumerate(exons, start=1):
            exons[i-1] = Exon(ex.chrom, ex.start1b, ex.end1b, ex.strand, f"{safe_orf}__exon{i}", ex.order)

    # 顺序处理每个 ORF（无并行）
    out_chr_dir = Path(args.output_root) / chr_name
    ensure_dir(out_chr_dir / "orfs"); ensure_dir(out_chr_dir / "beds")
    total = len(orf_exons); done = 0; ok = miss = fail = 0
    for orf_id, exons in orf_exons.items():
        try:
            status = stitch_orf(orf_id, exons, exon_out_dir, args.focal, args.cov_thresh,
                                out_chr_dir / "orfs", out_chr_dir / "beds")
            if status == "OK": ok += 1
            elif status == "MISS": miss += 1
            else: fail += 1
        except Exception as e:
            fail += 1
            print(f"[FAIL] {orf_id}: {e}")
        done += 1
        if done % 50 == 0 or done == total:
            print(f"[PROG] {chr_name}: {done}/{total}")

    print(f"[DONE] {chr_name}: OK={ok}, MISS={miss}, FAIL={fail}, TOTAL={total}")

if __name__ == "__main__":
    main()
