#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, io, gzip, argparse, subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio import AlignIO
from Bio.Seq import Seq

# ----------------------------- 数据结构 -----------------------------
class TransObject:
    def __init__(self, chrm, gene, strand, starts, ends):
        self.chrm = chrm
        self.gene = gene
        self.strand = strand
        self.start = starts  # list[int], 1-based
        self.end = ends      # list[int], 1-based

# ----------------------------- 工具函数 -----------------------------
def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def parse_orf_bed_1based_single(bed_path: Path) -> TransObject:
    """
    读取“单 ORF 的 1-based BED”（多行=多外显子）。
    格式假定：chrom  start  end  orf_id  score  strand
    """
    chrs, starts, ends, strands, ids = [], [], [], [], []
    with bed_path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                raise SystemExit(f"[ERROR] BED 列不足6列: {bed_path}\n{line}")
            chrs.append(parts[0] if parts[0].startswith("chr") else "chr"+parts[0])
            starts.append(int(parts[1]))
            ends.append(int(parts[2]))
            ids.append(parts[3])
            strands.append(parts[5])
    if not starts:
        raise SystemExit(f"[ERROR] 空BED: {bed_path}")

    # 统一 chr/strand/id
    chrm = chrs[0]
    strand = strands[0]
    gene = ids[0]
    # 成对排序：按 start 升序对 (start,end) 重新排序，避免 start/end 独立排序导致错配
    exons = sorted(zip(starts, ends), key=lambda x: x[0])
    starts = [s for s, e in exons]
    ends = [e for s, e in exons]
    return TransObject(chrm, gene, strand, starts, ends)

def extract_sequences(ma_block, pin, pif, regin, focal):
    """
    从单个 MAF block 中抽取 [pin, pif]（1-based，闭区间）对应的对齐列，返回 dict:
      species -> [seq, p11, p22, chr, block_start, block_strand]
    逻辑基本继承原脚本，但 focal 物种可配置。
    """
    all_seqs = {}
    p1 = p2 = -1  # 对齐列起止（0-based, inclusive）
    # 先定位 focal 在该 block 的列范围
    for s in ma_block:
        sp = s.name.split(".")[0]
        if sp not in all_seqs:
            all_seqs[sp] = ["", 0, 0, s.name.split(".")[1], s.annotations["start"], s.annotations["strand"]]

    # 找 focal
    for s in ma_block:
        sp = s.name.split(".")[0]
        if sp != focal:
            continue
        counter = regin  # block的基因组1-based起点？(依原脚本，regin来自 ma[0].annotations["start"])
        # 注意：PHAST/MAF里 start 多为0-based，这里沿用原脚本相对逻辑
        fseq = ""
        for n, ss in enumerate(str(s.seq)):
            if (counter >= pin - 1) and (counter <= pif - 1):
                if p1 == -1:
                    p1 = n
                fseq += ss
            elif (counter > pif - 1) and (p2 == -1):
                p2 = n - 1
            if ss != "-":
                counter += 1
        if p2 == -1:
            p2 = n
        break

    if p1 == -1:
        # 该 block 不覆盖目标区间
        return None

    # 取所有物种在 [p1, p2] 的列
    for s in ma_block:
        sp = s.name.split(".")[0]
        fseq = ""
        p11 = 0
        p22 = 0
        for n, ss in enumerate(str(s.seq)):
            if (n < p1) and (ss != "-"):
                p11 += 1
            if (n < p2) and (ss != "-"):
                p22 += 1
            if (n >= p1) and (n <= p2):
                fseq += ss
        all_seqs[sp] = [fseq, p11, p22, s.name.split(".")[1], s.annotations["start"], s.annotations["strand"]]
    return all_seqs

def maf_subset_stream(maf_gz: Path, start_1b: int, end_1b: int) -> str:
    """
    通过管道：zcat chrN.maf.gz | maf_parse -s START -e END -o MAF
    返回子对齐文本（MAF 格式字符串）。若无命中，返回空串。
    """
    # 注意：maf_parse 的 -s/-e 均使用 PHAST 坐标约定；这里沿用原代码的调用形式
    # 使用shell管道，便于直接对.gz 流式解压
    cmd = f"zcat -f {maf_gz} | maf_parse -s {start_1b} -e {end_1b} -o MAF"
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        # maf_parse 可能返回非0，但依然打印了空结果，给出告警
        # print(f"[WARN] maf_parse非0：{cmd}\n{p.stderr}", file=sys.stderr)
        pass
    return p.stdout

def get_seq_maf_stream(orf: TransObject, maf_gz: Path, focal: str):
    """
    替代原 get_seq_maf：不落盘子MAF，直接用管道获取子对齐并解析。
    返回：
      orf_seq: {species: [seq, [block_start+p11], [block_start+p22], [chr], [strand]]}
      l: ORF 基因组总长度（累加各exon的 end-start+1）
      coords: ["blockStart/blockSize", ...]（与原脚本保持）
    """
    orf_seq = {}
    coords = []
    total_len = 0

    for idx, (sto, eno) in enumerate(zip(orf.start, orf.end), start=1):
        total_len += (eno - sto + 1)
        # 抓子对齐（MAF文本）
        maf_txt = maf_subset_stream(maf_gz, sto, eno)
        if not maf_txt.strip():
            continue

        # 解析为 block 列表
        for ma in AlignIO.parse(io.StringIO(maf_txt), "maf"):
            # block 元注释（按 Biopython 的 Maf parser）
            st = int(ma[0].annotations.get("start", 0))
            strand = ma[0].annotations.get("strand", "+")
            usize = int(ma[0].annotations.get("size", 0))
            tsize = int(ma[0].annotations.get("srcSize", 0))

            key = f"{st}/{usize}"
            if key in coords:
                continue

            # 与原脚本大致相同的覆盖判断/拼接逻辑
            def append_all(all_seqs):
                coords.append(key)
                for sp, vals in all_seqs.items():
                    if sp not in orf_seq:
                        orf_seq[sp] = ["", [], [], [], []]
                    # seq
                    orf_seq[sp][0] += vals[0]
                    # block-based坐标（起止都以 block.start + 去缺口偏移 ）
                    orf_seq[sp][1].append(vals[4] + vals[1])
                    orf_seq[sp][2].append(vals[4] + vals[2])
                    orf_seq[sp][3].append(vals[3])  # chr
                    orf_seq[sp][4].append(vals[5])  # strand

            if (sto >= st) and (sto <= st + usize):
                all_seqs = extract_sequences(ma, sto, eno, st, focal)
                if all_seqs is None:
                    continue
                append_all(all_seqs)
            elif (eno < st):
                break
            else:
                # 可能是继续拼接的后续 block（保留与原逻辑接近的宽松处理）
                all_seqs = extract_sequences(ma, sto, eno, st, focal)
                if all_seqs is None:
                    continue
                append_all(all_seqs)

    return orf_seq, total_len, coords

# ----------------------------- 单 ORF 任务 -----------------------------
def process_one_orf(bed_path: Path, maf_gz: Path, focal: str, out_root: Path, cov_thresh: float = 0.95):
    """
    输入：单 ORF 的 1-based BED、对应 chr 的 maf.gz、输出根目录、焦点物种、覆盖阈值
    产物：
      out_root/orfs/<orf>.fa
      out_root/orfs/<orf>.maf
      out_root/beds/<orf>.bed
    """
    orf = parse_orf_bed_1based_single(bed_path)
    # 输出路径
    ensure_dir(out_root / "orfs")
    ensure_dir(out_root / "beds")

    # 结果存在则跳过
    fa_out = out_root / "orfs" / f"{orf.gene.replace(':','__')}.fa"
    maf_out = out_root / "orfs" / f"{orf.gene.replace(':','__')}.maf"
    bed_out = out_root / "beds" / f"{orf.gene.replace(':','__')}.bed"
    if fa_out.exists() and maf_out.exists() and bed_out.exists():
        return f"[SKIP] {orf.gene}"

    # 抽取
    orf_seq, L, coords = get_seq_maf_stream(orf, maf_gz, focal)
    if focal not in orf_seq or not orf_seq[focal][0]:
        # 与原脚本行为一致：没有hg38序列就跳过
        return f"[MISS] {orf.gene}\t0"

    # 覆盖度
    m = str(float(len(orf_seq[focal][0].replace("-", "")) / float(L)))

    # 写输出
    with fa_out.open("w") as f1, maf_out.open("w") as f2, bed_out.open("w") as f3:
        focal_len = len(orf_seq[focal][0])
        for sp, vals in orf_seq.items():
            sq = vals[0]
            # 负链取反向互补
            if orf.strand == "-":
                sq = str(Seq(sq).reverse_complement())

            # 覆盖比例过滤（相对于 focal 对齐长度）
            if focal_len == 0:
                continue
            if (len(sq) / focal_len) >= cov_thresh:
                # .fa：去缺口后翻译（与原脚本一致：cds=False）
                f1.write(f">{orf.gene}_{orf.gene}_{sp}_{m}\n{str(Seq(sq.replace('-','')).translate(cds=False))}\n")
                # .maf：对齐核酸序列（带缺口）
                f2.write(f">{sp}\n{sq}\n")
                # .bed：非焦点物种写坐标
                if sp != focal:
                    # 若该物种跨多 block，chr 名不唯一则附加 _splitted（与原逻辑一致）
                    t2 = "_splitted" if len(set(vals[3])) != 1 else ""
                    for n in range(len(vals[1])):
                        strand = str(vals[4][n]).replace("-1", "-").replace("1", "+")
                        f3.write(f"{vals[3][n]}\t{vals[1][n]}\t{vals[2][n]}\t{orf.gene}\t{sp}{t2}\t{strand}\n")

    return f"[OK] {orf.gene}"

# ----------------------------- 主流程：按染色体批次并行 -----------------------------
def main():
    ap = argparse.ArgumentParser(description="按染色体逐个处理，染色体内并行 ORF；支持 chr*.maf.gz 的并行抽取版本")
    ap.add_argument("--per-chr-bed-root", required=True, help="每条染色体下是若干单ORF的1-based BED目录，如 ../processed/per_chr/ORFs_bed/chr1/")
    ap.add_argument("--maf-root", required=True, help="存放 chrN.maf.gz 的目录")
    ap.add_argument("--output-root", required=True, help="输出根目录（每条chr在该目录下建子目录）")
    ap.add_argument("--chr", default="", help="只跑某条染色体（如 chr1）。留空则遍历 per-chr 目录下的所有 chr* 子目录")
    ap.add_argument("--workers", type=int, default=5, help="每条染色体内的并发ORF数（默认5）")
    ap.add_argument("--focal", default="hg38", help="焦点物种名称（默认 hg38）")
    ap.add_argument("--cov-thresh", type=float, default=0.95, help="物种输出的覆盖阈值（默认0.95）")
    args = ap.parse_args()

    per_chr_root = Path(args.per_chr_bed_root)
    maf_root = Path(args.maf_root)
    out_root = Path(args.output_root)

    # 拿到要处理的 chr 列表
    if args.chr:
        chr_dirs = [per_chr_root / args.chr]
    else:
        chr_dirs = sorted([p for p in per_chr_root.glob("chr*") if p.is_dir()])

    if not chr_dirs:
        sys.exit(f"[ERROR] 未找到染色体子目录于: {per_chr_root}")

    for chr_dir in chr_dirs:
        chr_name = chr_dir.name
        maf_gz = maf_root / f"{chr_name}.maf.gz"
        if not maf_gz.exists():
            print(f"[WARN] 缺少 {maf_gz}，跳过 {chr_name}")
            continue

        out_chr = out_root / chr_name
        ensure_dir(out_chr / "orfs")
        ensure_dir(out_chr / "beds")

        bed_files = sorted(chr_dir.glob("*.bed"))
        if not bed_files:
            print(f"[INFO] {chr_name} 无 BED，跳过")
            continue

        print(f"[INFO] 开始处理 {chr_name}：{len(bed_files)} ORFs，workers={args.workers}")
        results = []
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = []
            for bed in bed_files:
                futures.append(
                    ex.submit(process_one_orf, bed, maf_gz, args.focal, out_chr, args.cov_thresh)
                )
            for fu in as_completed(futures):
                try:
                    msg = fu.result()
                except Exception as e:
                    msg = f"[FAIL] {e}"
                results.append(msg)
                if len(results) % 50 == 0:
                    print(f"[PROG] {chr_name}: 已完成 {len(results)}/{len(bed_files)}")

        # 小结
        ok = sum(1 for x in results if x.startswith("[OK]"))
        skip = sum(1 for x in results if x.startswith("[SKIP]"))
        miss = sum(1 for x in results if x.startswith("[MISS]"))
        fail = sum(1 for x in results if x.startswith("[FAIL]"))
        print(f"[DONE] {chr_name}: OK={ok}, SKIP={skip}, MISS={miss}, FAIL={fail}, TOTAL={len(bed_files)}")

if __name__ == "__main__":
    main()
