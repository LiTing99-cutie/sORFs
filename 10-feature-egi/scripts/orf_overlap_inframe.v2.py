#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
新ORF与注释CDS的重叠与in-frame统计 — 并行版（按染色体切分，带日志）
依赖：bedtools, python3, pandas, numpy
"""

import os, sys, argparse, tempfile, subprocess, shutil, time, logging
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import numpy as np

# ---------- 基础工具 ----------

def check_exec(cmd):
    from shutil import which
    if which(cmd) is None:
        sys.exit(f"ERROR: 需要可执行程序 '{cmd}'，请先安装并加入PATH。")

def run_cmd(args, stdout_path=None):
    if stdout_path is None:
        subprocess.check_call(args)
    else:
        with open(stdout_path, "w") as fo:
            subprocess.check_call(args, stdout=fo)

def parse_attrs(attr_str):
    d = {}
    for part in attr_str.strip().strip(';').split(';'):
        if not part.strip():
            continue
        x = part.strip().split(' ', 1)
        if len(x) == 2:
            key = x[0].strip()
            val = x[1].strip().strip('"')
            d[key] = val
    return d

def gtf_to_cds_bed_with_phase(gtf_path, id_attr, tx_attr, gene_attr, out_bed, lengths_out=None):
    groups = defaultdict(list)
    with open(gtf_path) as f:
        for ln in f:
            if not ln or ln.startswith("#"): continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS": continue
            chrom, start, end, strand, attrs = parts[0], int(parts[3])-1, int(parts[4]), parts[6], parts[8]
            a = parse_attrs(attrs)
            the_id = a.get(id_attr) or a.get(tx_attr)
            if not the_id: 
                continue
            tx = a.get(tx_attr, "")
            gene = a.get(gene_attr, "")
            groups[(the_id, strand)].append((chrom, start, end, gene, tx))
    os.makedirs(os.path.dirname(out_bed), exist_ok=True)
    id2len = defaultdict(int)
    with open(out_bed, "w") as out:
        for (the_id, strand), segs in groups.items():
            segs.sort(key=lambda x: x[1], reverse=(strand=="-"))
            phase = 0
            for (chrom, s0, e1, gene, tx) in segs:
                length = e1 - s0
                print(chrom, s0, e1, the_id, 0, strand, phase, gene, sep="\t", file=out)
                phase = (phase + length) % 3
                id2len[the_id] += length
    if lengths_out:
        with open(lengths_out, "w") as fo:
            for _id, L in id2len.items():
                print(_id, L, sep="\t", file=fo)
    return dict(id2len)

def bedtools_sort(in_bed, out_bed, faidx=None):
    args = ["bedtools", "sort"]
    if faidx: args += ["-faidx", faidx]
    args += ["-i", in_bed]
    run_cmd(args, stdout_path=out_bed)

def build_union_bed_nostrand(sorted_ann_bed, out_union_bed):
    threecol = out_union_bed + ".3col"
    merged = out_union_bed + ".merged"
    with open(sorted_ann_bed) as fi, open(threecol, "w") as fo:
        for ln in fi:
            p = ln.rstrip("\n").split("\t")
            print(p[0], p[1], p[2], sep="\t", file=fo)
    run_cmd(["bedtools", "merge", "-i", threecol], stdout_path=merged)
    with open(merged) as fi, open(out_union_bed, "w") as fo:
        for ln in fi:
            c, s, e = ln.rstrip("\n").split("\t")
            print(c, s, e, ".", 0, ".", sep="\t", file=fo)
    os.remove(threecol); os.remove(merged)

# ---------- 按染色体切分 ----------
def split_bed_by_chrom(bed_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    writers = {}
    counts = defaultdict(int)
    with open(bed_path) as f:
        for ln in f:
            if not ln.strip(): continue
            c = ln.split("\t", 1)[0]
            if c not in writers:
                writers[c] = open(os.path.join(out_dir, f"{c}.bed"), "w")
            writers[c].write(ln)
            counts[c] += 1
    for fh in writers.values(): fh.close()
    return {c: os.path.join(out_dir, f"{c}.bed") for c in counts.keys()}

# ---------- 交叠结果汇总 ----------
def agg_union_overlap(wo_path, sums):
    with open(wo_path) as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            oid = p[3]
            sums[oid] += int(p[-1])

def parse_tx_same_and_update(wo_path, sum_tx, sum_gene, inframe_writer):
    with open(wo_path) as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            # A(new bed8)
            chromA = p[0]
            sA = int(p[1]); eA = int(p[2])
            oid = p[3]; strandA = p[5]; phA = int(p[6])
            # B(ann bed8)
            sB = int(p[9]); eB = int(p[10])
            tx = p[11]; strandB = p[13]; phB = int(p[14]); geneB = p[15]
            ov = int(p[-1])

            s0 = max(sA, sB); e0 = min(eA, eB)
            if e0 <= s0: 
                continue
            # 单次判定 in-frame
            if strandA == '+':
                lfA = (phA + (s0 - sA)) % 3
                lfB = (phB + (s0 - sB)) % 3
                is_in = (lfA == lfB)
            else:
                rfA = (phA + (eA - 1 - s0)) % 3
                rfB = (phB + (eB - 1 - s0)) % 3
                is_in = (rfA == rfB)

            sum_tx[oid][tx] += ov
            sum_gene[oid][geneB] += ov

            if is_in:
                # 为了后续每ORF独立merge，这里把“染色体名”改成 "chrom|ORF_id"
                inframe_writer.write(f"{chromA}|{oid}\t{s0}\t{e0}\n")

def merged_len_per_orf_fast(inframe_bed, tmpdir=None):
    if (not os.path.exists(inframe_bed)) or os.path.getsize(inframe_bed) == 0:
        return {}
    sorted_bed = inframe_bed + ".sorted"
    merged_bed = inframe_bed + ".merged"
    bedtools_sort(inframe_bed, sorted_bed)
    run_cmd(["bedtools", "merge", "-i", sorted_bed], stdout_path=merged_bed)
    res = defaultdict(int)
    with open(merged_bed) as f:
        for ln in f:
            c, s, e = ln.rstrip("\n").split("\t")
            s = int(s); e = int(e)
            # c 形如 "chr1|PB.xxx"
            try:
                _, oid = c.split("|", 1)
            except ValueError:
                continue
            res[oid] += (e - s)
    os.remove(sorted_bed); os.remove(merged_bed)
    return dict(res)

# ---------- 并行 worker ----------
def worker_intersects(chrom, new_chr_bed, ann_union_chr_bed, ann_tx_chr_bed, out_dir):
    """每个染色体：跑两次intersect，返回结果文件路径"""
    t0 = time.time()
    unique_wo = os.path.join(out_dir, f"{chrom}.union.wo")
    txsame_wo = os.path.join(out_dir, f"{chrom}.txsame.wo")
    # 无链（唯一覆盖）
    run_cmd(["bedtools", "intersect", "-wo", "-sorted", "-a", new_chr_bed, "-b", ann_union_chr_bed],
            stdout_path=unique_wo)
    # 同链
    run_cmd(["bedtools", "intersect", "-s", "-wo", "-sorted", "-a", new_chr_bed, "-b", ann_tx_chr_bed],
            stdout_path=txsame_wo)
    elapsed = time.time() - t0
    return chrom, unique_wo, txsame_wo, elapsed

# ---------- 主流程 ----------
def main():
    ap = argparse.ArgumentParser(description="新ORF与注释CDS重叠与in-frame统计（并行版）")
    ap.add_argument("--new-gtf", required=True)
    ap.add_argument("--ann-gtf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--new-id-attr", default="ORF_id")
    ap.add_argument("--ann-tx-attr", default="transcript_id")
    ap.add_argument("--ann-gene-attr", default="gene_id")
    ap.add_argument("--tmpdir", default=None)
    ap.add_argument("--workers", type=int, default=8, help="并行进程数（按染色体）")
    ap.add_argument("--faidx", default=None, help="可选：genome.fa.fai，用于更快的bedtools sort")
    args = ap.parse_args()

    # 日志
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S"
    )
    log = logging.getLogger("overlap")

    # 检查依赖
    for cmd in ["bedtools"]:
        check_exec(cmd)

    tmp = tempfile.mkdtemp(prefix="orf_overlaps_parallel_", dir=args.tmpdir)
    try:
        new_bed = os.path.join(tmp, "new.bed8")
        ann_bed = os.path.join(tmp, "ann_tx.bed8")
        new_bed_sorted = new_bed + ".sorted"
        ann_bed_sorted = ann_bed + ".sorted"

        log.info("Step1: GTF -> BED（重建phase）")
        orf_len = gtf_to_cds_bed_with_phase(
            args.new_gtf, args.new_id_attr, "transcript_id", "gene_id",
            new_bed)
        gtf_to_cds_bed_with_phase(
            args.ann_gtf, args.ann_tx_attr, args.ann_tx_attr, args.ann_gene_attr,
            ann_bed)

        log.info("Step2: bedtools sort")
        bedtools_sort(new_bed, new_bed_sorted, faidx=args.faidx)
        bedtools_sort(ann_bed, ann_bed_sorted, faidx=args.faidx)

        log.info("Step3: 构建注释CDS无链并集")
        ann_union_nostrand = os.path.join(tmp, "ann_union_nostrand.bed6")
        build_union_bed_nostrand(ann_bed_sorted, ann_union_nostrand)

        log.info("Step4: 按染色体切分BED")
        split_dir = os.path.join(tmp, "split"); os.makedirs(split_dir, exist_ok=True)
        new_split = split_bed_by_chrom(new_bed_sorted, os.path.join(split_dir, "new"))
        ann_tx_split = split_bed_by_chrom(ann_bed_sorted, os.path.join(split_dir, "ann_tx"))
        ann_union_split = split_bed_by_chrom(ann_union_nostrand, os.path.join(split_dir, "ann_union"))

        chroms = sorted(set(new_split.keys()) & set(ann_union_split.keys()))
        if not chroms:
            sys.exit("ERROR: 没有可并行的染色体分片（new 与 ann_union 没有交集）。")
        log.info(f"并行染色体数：{len(chroms)}；workers={args.workers}")

        # 汇总容器
        unique_sum = defaultdict(int)
        sum_tx = defaultdict(lambda: defaultdict(int))
        sum_gene = defaultdict(lambda: defaultdict(int))
        inframe_intervals_bed = os.path.join(tmp, "inframe_intervals.bed")
        inframe_writer = open(inframe_intervals_bed, "w")

        # 并行执行
        t0 = time.time()
        out_dir = os.path.join(tmp, "intersects"); os.makedirs(out_dir, exist_ok=True)
        futures = []
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            for chrom in chroms:
                # ann_tx 可能没有该染色体，跳过这一条
                if chrom not in ann_tx_split:
                    continue
                f = ex.submit(
                    worker_intersects, chrom,
                    new_split[chrom], ann_union_split[chrom], ann_tx_split[chrom], out_dir
                )
                futures.append(f)

            done_cnt = 0
            for fu in as_completed(futures):
                chrom, unique_wo, txsame_wo, elapsed = fu.result()
                done_cnt += 1
                log.info(f"[{done_cnt}/{len(futures)}] {chrom} finished in {elapsed:.1f}s")

                # 边完成边汇总，避免生成超大中间文件
                agg_union_overlap(unique_wo, unique_sum)
                parse_tx_same_and_update(txsame_wo, sum_tx, sum_gene, inframe_writer)

        inframe_writer.close()
        log.info(f"并行交叠完成，用时 {time.time()-t0:.1f}s；开始in-frame去重合并")

        # in-frame 合并去重（bedtools merge），得到每个ORF的inframe长度
        inframe_uniq = merged_len_per_orf_fast(inframe_intervals_bed)

        log.info("Step5: 汇总统计并输出")
        df = pd.DataFrame({"ORF_id": list(orf_len.keys()), "ORF_len": list(orf_len.values())})
        orf_ids = df["ORF_id"].to_numpy()
        orf_len_arr = df["ORF_len"].to_numpy(dtype=np.int64)

        overlap_bp = np.array([unique_sum.get(oid, 0) for oid in orf_ids], dtype=np.int64)
        inframe_bp = np.array([inframe_uniq.get(oid, 0) for oid in orf_ids], dtype=np.int64)

        # n_tx / n_gene / top_hit
        n_tx_arr = np.fromiter((len(sum_tx.get(oid, {})) for oid in orf_ids), dtype=np.int32)
        n_gene_arr = np.fromiter((len(sum_gene.get(oid, {})) for oid in orf_ids), dtype=np.int32)
        def top_by_sum(d):
            return max(d.items(), key=lambda x: (x[1], x[0]))[0] if d else ""
        top_tx_arr = np.array([top_by_sum(sum_tx.get(oid, {})) for oid in orf_ids], dtype=object)
        top_gene_arr = np.array([top_by_sum(sum_gene.get(oid, {})) for oid in orf_ids], dtype=object)

        with np.errstate(divide='ignore', invalid='ignore'):
            overlap_fraction = np.where(orf_len_arr > 0, overlap_bp / orf_len_arr, 0.0)
            inframe_fraction = np.where(orf_len_arr > 0, inframe_bp / orf_len_arr, 0.0)

        out_df = pd.DataFrame({
            "ORF_id": orf_ids,
            "overlap_bp": overlap_bp,
            "overlap_fraction": overlap_fraction.astype("float32"),
            "inframe_overlap_bp": inframe_bp,
            "inframe_fraction": inframe_fraction.astype("float32"),
            "n_annot_tx_overlapped": n_tx_arr,
            "n_annot_gene_overlapped": n_gene_arr,
            "top_hit_gene": top_gene_arr,
            "top_hit_tx": top_tx_arr
        })
        out_df.to_csv(args.out, sep="\t", index=False)
        log.info(f"[OK] 写出结果：{args.out}")

    finally:
        shutil.rmtree(tmp, ignore_errors=True)

if __name__ == "__main__":
    main()
