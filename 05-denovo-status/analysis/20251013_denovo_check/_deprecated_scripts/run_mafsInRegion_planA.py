#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch per-chromosome ORF extraction using UCSC mafsInRegion (one pass),
then per-ORF stitching & export (parallel).

- Read all 1-based ORF BEDs under a chromosome folder
- Merge all exons -> one combined BED (0-based for mafsInRegion), with unique "name" per exon
- mafsInRegion -outDir once: produce per-exon MAF files named by BED name
- Reassemble per-ORF (respect exon order), output:
    <out_root>/<chr>/orfs/<ORF>.fa
    <out_root>/<chr>/orfs/<ORF>.maf
    <out_root>/<chr>/beds/<ORF>.bed
"""

import os
import re
import io
import argparse
import subprocess
from pathlib import Path
from collections import defaultdict, namedtuple
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio import AlignIO
from Bio.Seq import Seq

# ----------------------------- helpers -----------------------------
def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def sanitize_filename(s: str) -> str:
    """safe file name: keep [A-Za-z0-9._-], others -> '__' """
    return re.sub(r"[^A-Za-z0-9._-]+", "__", s)

Exon = namedtuple("Exon", ["chrom", "start1b", "end1b", "strand", "bed_name", "order"])

def read_orf_beds_1based(per_chr_bed_dir: Path):
    """
    读取该染色体目录下所有 ORF 的 1-based BED（每个文件一个 ORF，多行=多外显子）。
    返回：
      orf_exons: dict[orf_id] -> list[Exon]（按 start 升序）
      chrom_name: 'chrN'
    """
    orf_exons = {}
    chrom_name = per_chr_bed_dir.name  # e.g., 'chr1'
    for bed in sorted(per_chr_bed_dir.glob("*.bed")):
        orf_id = bed.stem.replace(".ORF", "")  # 仅用于日志；以 BED 第4列为准
        rows = []
        with bed.open() as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                ps = line.rstrip("\n").split("\t")
                if len(ps) < 6:
                    raise ValueError(f"BED列不足6：{bed}\n{line}")
                chrom = ps[0] if ps[0].startswith("chr") else "chr"+ps[0]
                s1 = int(ps[1]); e1 = int(ps[2])
                name = ps[3]   # 以第4列为真正 ORF ID
                strand = ps[5]
                rows.append((chrom, s1, e1, strand, name))
        if not rows:
            continue
        # 按 start 排序，并为每个外显子分配顺序号
        rows.sort(key=lambda x: x[1])
        exons = []
        for i, (chrom, s1, e1, strand, name) in enumerate(rows, start=1):
            exons.append(Exon(chrom, s1, e1, strand, bed_name=None, order=i))
        # key = ORF ID（取第4列），同一个 ORF 的多文件不应重复
        orf_id = rows[0][4]
        orf_exons[orf_id] = exons
    if not orf_exons:
        raise SystemExit(f"[ERROR] {per_chr_bed_dir} 未找到任何 BED")
    return orf_exons, chrom_name

def build_combined_bed(orf_exons, chrom_name, combined_bed: Path, name_map_path: Path):
    """
    组合所有外显子为一个 0-based 半开 BED，name 字段唯一且安全，用于 mafsInRegion -outDir。
    - name 采用：{SAFE_ORF}__exon{order}
    - 同时写出 name 映射表：name \t ORF_ID \t order \t strand
    返回：name->(orf_id, order, strand) 的 dict
    """
    ensure_dir(combined_bed.parent)
    ensure_dir(name_map_path.parent)

    name_map = {}
    with combined_bed.open("w") as bedout, name_map_path.open("w") as mapout:
        for orf_id, exons in orf_exons.items():
            safe_orf = sanitize_filename(orf_id)
            strand_set = {ex.strand for ex in exons}
            if len(strand_set) != 1:
                # 若出现混合链，仍然按每行的 strand 记录，但给个警告
                print(f"[WARN] ORF {orf_id} 多链：{strand_set}")
            for ex in exons:
                # BED for mafsInRegion: 0-based half-open
                start0 = ex.start1b - 1
                end0   = ex.end1b
                name   = f"{safe_orf}__exon{ex.order}"
                # 写 BED：chrom start end name
                bedout.write(f"{chrom_name}\t{start0}\t{end0}\t{name}\n")
                # 记录映射
                name_map[name] = (orf_id, ex.order, ex.strand)
                mapout.write(f"{name}\t{orf_id}\t{ex.order}\t{ex.strand}\n")
    return name_map

def run_mafsInRegion(mafsInRegion_bin: str, combined_bed: Path, out_dir: Path, maf_gz: Path):
    """
    mafsInRegion regions.bed outDir in.maf.gz
    使用 -outDir：按 BED name 输出独立小 MAF 文件到 out_dir
    """
    ensure_dir(out_dir)
    cmd = [
        mafsInRegion_bin,
        str(combined_bed),
        str(out_dir),          # outDir
        str(maf_gz)            # in.maf.gz
    ]
    # mafsInRegion 会把警告写到 stderr，这里仅检查返回码
    print(f"[INFO] mafsInRegion: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def exon_maf_path(out_dir: Path, bed_name: str):
    """mafsInRegion -outDir 会以 BED 的 name 命名文件。常见两种：name 或 name.maf"""
    p1 = out_dir / bed_name
    p2 = out_dir / f"{bed_name}.maf"
    if p1.exists():
        return p1
    if p2.exists():
        return p2
    # 允许某些版本加前缀或后缀，这里做一次扫描兜底
    for p in out_dir.glob(f"{bed_name}*"):
        return p
    return None

# 与你原脚本一致的列切片逻辑（基于 focal 序列定位列范围），更稳健
def extract_sequences_from_block(ma, pin_1b, pif_1b, regin0, focal_prefix):
    """
    从1个 MAF block中抽取 [pin_1b, pif_1b] 的列范围（闭区间），返回：
      dict[species] = [seq, p11, p22, chr, block_start0, strand(+1/-1)]
    这里 regin0 应为 focal 在该 block 的 0-based 起点（来自 annotations['start']）
    """
    all_seqs = {}
    for s in ma:
        sp = s.id.split(".")[0]
        if sp not in all_seqs:
            chrom = s.id.split(".", 1)[1] if "." in s.id else s.id
            all_seqs[sp] = ["", 0, 0, chrom, int(s.annotations.get("start", 0)), int(s.annotations.get("strand", 1))]

    # 定位 focal 列区间
    p1 = p2 = -1
    for s in ma:
        if not s.id.startswith(focal_prefix + "."):
            continue
        counter = regin0
        for n, ch in enumerate(str(s.seq)):
            if (counter >= pin_1b - 1) and (counter <= pif_1b - 1):
                if p1 == -1:
                    p1 = n
            elif (counter > pif_1b - 1) and (p2 == -1):
                p2 = n - 1
            if ch != "-":
                counter += 1
        if p2 == -1:
            p2 = n
        break
    if p1 == -1:
        return None

    # 切所有物种
    for s in ma:
        sp = s.id.split(".")[0]
        text = str(s.seq)
        p11 = p22 = 0
        for n, ch in enumerate(text):
            if (n < p1) and (ch != "-"): p11 += 1
            if (n < p2) and (ch != "-"): p22 += 1
        subseq = text[p1:p2+1]
        chrom  = s.id.split(".", 1)[1] if "." in s.id else s.id
        all_seqs[sp] = [subseq, p11, p22, chrom, int(s.annotations.get("start", 0)), int(s.annotations.get("strand", 1))]
    return all_seqs

def stitch_orf_from_exon_mafs(orf_id, exons, out_orfs_dir: Path, out_beds_dir: Path,
                              focal: str, cov_thresh: float, exon_maf_dir: Path):
    """
    把该 ORF 的每个外显子对应的小MAF（由 mafsInRegion 生成）按顺序拼接，输出 .fa / .maf / .bed
    """
    orf_seq = defaultdict(lambda: ["", [], [], [], []])  # sp -> [seq, starts[], ends[], chrs[], strands[]]
    L = 0
    coords = []
    focal_prefix = focal

    for ex in exons:
        L += (ex.end1b - ex.start1b + 1)
        maf_file = exon_maf_path(exon_maf_dir, ex.bed_name)
        if maf_file is None:
            continue
        # 该外显子可能切出多个 block
        for ma in AlignIO.parse(maf_file.open(), "maf"):
            st0 = int(ma[0].annotations.get("start", 0))
            usize = int(ma[0].annotations.get("size", 0))
            key = f"{st0}/{usize}"
            if key in coords:
                continue
            all_seqs = extract_sequences_from_block(ma, ex.start1b, ex.end1b, st0, focal_prefix)
            if all_seqs is None:
                continue
            coords.append(key)
            for sp, vals in all_seqs.items():
                orf_seq[sp][0] += vals[0]
                orf_seq[sp][1].append(vals[4] + vals[1])
                orf_seq[sp][2].append(vals[4] + vals[2])
                orf_seq[sp][3].append(vals[3])
                orf_seq[sp][4].append(vals[5])

    # 若缺少焦点物种，视为未命中
    if focal not in orf_seq or not orf_seq[focal][0]:
        return "MISS"

    # 输出
    safe_orf = sanitize_filename(orf_id)
    fa_out  = out_orfs_dir / f"{safe_orf}.fa"
    maf_out = out_orfs_dir / f"{safe_orf}.maf"
    bed_out = out_beds_dir / f"{safe_orf}.bed"

    m = float(len(orf_seq[focal][0].replace("-", "")) / float(L))
    focal_len = len(orf_seq[focal][0])

    with fa_out.open("w") as f1, maf_out.open("w") as f2, bed_out.open("w") as f3:
        for sp, vals in orf_seq.items():
            seq_aln = vals[0]
            # 负链时反向互补
            if exons[0].strand == "-":
                seq_aln = str(Seq(seq_aln).reverse_complement())
            # 覆盖过滤（相对 focal 对齐长度）
            if focal_len == 0 or (len(seq_aln) / focal_len) < cov_thresh:
                continue
            # .fa：去缺口翻译（与原脚本一致 cds=False），标题保留原风格（orf_id重复一次以兼容原逻辑）
            f1.write(f">{orf_id}_{orf_id}_{sp}_{m:.6f}\n{str(Seq(seq_aln.replace('-','')).translate(cds=False))}\n")
            # .maf：每物种一行对齐核酸
            f2.write(f">{sp}\n{seq_aln}\n")
            # .bed：非焦点物种输出坐标
            if sp != focal:
                t2 = "_splitted" if len(set(vals[3])) != 1 else ""
                for n in range(len(vals[1])):
                    strand_char = "-" if str(vals[4][n]) in ("-1", "-") else "+"
                    f3.write(f"{vals[3][n]}\t{vals[1][n]}\t{vals[2][n]}\t{orf_id}\t{sp}{t2}\t{strand_char}\n")

    return "OK"

# ----------------------------- main -----------------------------
def main():
    ap = argparse.ArgumentParser(description="Use mafsInRegion once per chromosome, then per-ORF stitching & export.")
    ap.add_argument("--per-chr-bed-dir", required=True, help="该染色体的 ORF BED 目录（1-based），如 ../processed/get_input/per_chr/ORFs_bed/chr1")
    ap.add_argument("--maf-root", required=True, help="存放 chrN.maf.gz 的目录")
    ap.add_argument("--output-root", required=True, help="输出根目录（每条chr会建子目录）")
    ap.add_argument("--work-dir", required=True, help="工作目录（放合并BED、mafsInRegion输出等临时/可复用文件）")
    ap.add_argument("--chr", default="", help="显式指定染色体名（留空则取 per-chr-bed-dir 的文件夹名）")
    ap.add_argument("--focal", default="hg38", help="焦点物种前缀（默认 hg38）")
    ap.add_argument("--cov-thresh", type=float, default=0.95, help="物种输出覆盖阈值（默认0.95）")
    ap.add_argument("--workers", type=int, default=5, help="染色体内 ORF 并发数（默认5）")
    ap.add_argument("--mafsInRegion", default="/home/user/data/lit/Tools/ucsc/bin/mafsInRegion", help="mafsInRegion 可执行文件路径")
    args = ap.parse_args()

    per_chr_bed_dir = Path(args.per_chr_bed_dir)
    out_root = Path(args.output_root)
    work_dir = Path(args.work_dir)
    ensure_dir(out_root); ensure_dir(work_dir)

    # 读取该 chr 下所有 ORF 的 BED
    orf_exons, chr_name_detected = read_orf_beds_1based(per_chr_bed_dir)
    chr_name = args.chr or chr_name_detected
    print(f"[INFO] 染色体 {chr_name}：ORFs={len(orf_exons)}")

    # 为 mafsInRegion 生成合并 BED（0-based）与 name 映射
    tmp_chr_dir = work_dir / chr_name
    ensure_dir(tmp_chr_dir)
    combined_bed = tmp_chr_dir / f"{chr_name}.all_orfs.0based.bed"
    name_map_path = tmp_chr_dir / f"{chr_name}.name_map.tsv"
    name_map = build_combined_bed(orf_exons, chr_name, combined_bed, name_map_path)

    # 运行 mafsInRegion（一次），输出到 outdir_exons
    maf_gz = Path(args.maf_root) / f"{chr_name}.maf.gz"
    if not maf_gz.exists():
        raise SystemExit(f"[ERROR] 找不到 {maf_gz}")
    outdir_exons = tmp_chr_dir / "mafsInRegion_exons"
    run_mafsInRegion(args.mafsInRegion, combined_bed, outdir_exons, maf_gz)

    # 把 name 映射回 exon 结构，写入 bed_name 以便查找对应的小MAF
    for orf_id, exons in orf_exons.items():
        safe_orf = sanitize_filename(orf_id)
        for i, ex in enumerate(exons, start=1):
            bed_name = f"{safe_orf}__exon{i}"
            if bed_name not in name_map:
                # 兼容性兜底：如果因为同名冲突等问题，尝试继续
                pass
            exons[i-1] = Exon(ex.chrom, ex.start1b, ex.end1b, ex.strand, bed_name, ex.order)

    # 输出目录
    out_chr_dir = out_root / chr_name
    ensure_dir(out_chr_dir / "orfs")
    ensure_dir(out_chr_dir / "beds")

    # 并行逐 ORF 拼接与导出
    results = []
    def _worker(item):
        orf_id, exons = item
        try:
            status = stitch_orf_from_exon_mafs(orf_id, exons, out_chr_dir / "orfs", out_chr_dir / "beds",
                                               args.focal, args.cov_thresh, outdir_exons)
            return f"[{status}] {orf_id}"
        except Exception as e:
            return f"[FAIL] {orf_id}: {e}"

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(_worker, it) for it in orf_exons.items()]
        for i, fu in enumerate(as_completed(futs), 1):
            msg = fu.result()
            results.append(msg)
            if i % 50 == 0 or i == len(orf_exons):
                print(f"[PROG] {chr_name}: {i}/{len(orf_exons)}")

    ok   = sum(1 for x in results if x.startswith("[OK]"))
    miss = sum(1 for x in results if x.startswith("[MISS]"))
    fail = sum(1 for x in results if x.startswith("[FAIL]"))
    print(f"[DONE] {chr_name}: OK={ok}, MISS={miss}, FAIL={fail}, TOTAL={len(orf_exons)})")

if __name__ == "__main__":
    main()
