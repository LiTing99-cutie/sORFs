#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
新ORF与注释CDS的重叠与in-frame统计（高效+去重版）

输出（每行一个 ORF）：
- ORF_id
- overlap_bp                # 被任意链的注释CDS覆盖的“独特碱基”数（无链并集，去重）
- overlap_fraction          # overlap_bp / ORF_CDS_len   ∈ [0,1]
- inframe_overlap_bp        # 同链且in-frame的“独特碱基”数（按 ORF 合并去重）
- inframe_fraction          # inframe_overlap_bp / ORF_CDS_len  ∈ [0,1]
- n_annot_tx_overlapped     # 同链交叠到的注释转录本数量（去重）
- n_annot_gene_overlapped   # 同链交叠到的注释基因数量（去重）
- top_hit_gene              # 同链按重叠bp汇总最大者
- top_hit_tx                # 同链按重叠bp汇总最大者

依赖：bedtools, python3, pandas, numpy, coreutils(sort)
示例：
python3 orf_overlap_inframe.py --new-gtf new_orf.gtf --ann-gtf gencode.v41.annotation.gtf --out result.tsv
"""

import os, sys, argparse, tempfile, subprocess, shutil
from collections import defaultdict, namedtuple
import pandas as pd
import numpy as np

# ---------- 工具与解析 ----------

def check_exec(cmd):
    from shutil import which
    if which(cmd) is None:
        sys.exit(f"ERROR: 需要可执行程序 '{cmd}'，请先安装并加入PATH。")

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
    """
    读取 GTF，提取 feature=CDS，按 (id_attr, strand) 分组，按转录本方向重建 phase（0/1/2）。
    写出 BED8： chrom start0 end id score(0) strand phase gene
    返回：对于“新ORF”调用时，字典 ORF_id->CDS总长度（bp）
    """
    groups = defaultdict(list)
    with open(gtf_path) as f:
        for ln in f:
            if not ln or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != "CDS":
                continue
            a = parse_attrs(attrs)
            the_id = a.get(id_attr) or a.get(tx_attr)
            if not the_id:
                continue
            tx = a.get(tx_attr, "")
            gene = a.get(gene_attr, "")
            s0 = int(start) - 1
            e1 = int(end)
            groups[(the_id, strand)].append((chrom, s0, e1, gene, tx))

    os.makedirs(os.path.dirname(out_bed), exist_ok=True)
    id2len = defaultdict(int)
    with open(out_bed, "w") as out:
        for (the_id, strand), segs in groups.items():
            if strand == '+':
                segs.sort(key=lambda x: x[1])
            else:
                segs.sort(key=lambda x: x[1], reverse=True)
            phase = 0
            for (chrom, s0, e1, gene, tx) in segs:
                length = e1 - s0
                print(chrom, s0, e1, the_id, 0, strand, phase, gene, sep="\t", file=out)
                phase = (phase + length) % 3
                id2len[the_id] += length

    if lengths_out is not None:
        with open(lengths_out, "w") as fo:
            for _id, L in id2len.items():
                print(_id, L, sep="\t", file=fo)
    return dict(id2len)

def bedtools_sort(in_bed, out_bed):
    with open(out_bed, "w") as fo:
        subprocess.check_call(["bedtools", "sort", "-i", in_bed], stdout=fo)

def build_union_bed_from_ann(sorted_ann_bed, out_union_bed):
    """
    构建“分链并集”：仅为备选用途（本脚本最终 overlap 不使用分链并集）。
    输入：注释逐Tx BED8（chrom start end id 0 strand phase gene）
    输出：BED6（chrom start end name '.' strand）
    """
    tmp_plus = out_union_bed + ".plus.bed"
    tmp_minus = out_union_bed + ".minus.bed"
    tmp_plus_m = out_union_bed + ".plus.merged"
    tmp_minus_m = out_union_bed + ".minus.merged"
    tmp_union = out_union_bed + ".union.tmp"

    with open(sorted_ann_bed) as fi, open(tmp_plus, "w") as fp, open(tmp_minus, "w") as fm:
        for ln in fi:
            p = ln.rstrip("\n").split("\t")
            strand = p[5]
            if strand == '+':
                print(p[0], p[1], p[2], sep="\t", file=fp)
            elif strand == '-':
                print(p[0], p[1], p[2], sep="\t", file=fm)

    with open(tmp_plus_m, "w") as fo:
        subprocess.check_call(["bedtools", "merge", "-i", tmp_plus], stdout=fo)
    with open(tmp_minus_m, "w") as fo:
        subprocess.check_call(["bedtools", "merge", "-i", tmp_minus], stdout=fo)

    with open(tmp_union, "w") as fo:
        with open(tmp_plus_m) as fp:
            for ln in fp:
                c, s, e = ln.rstrip("\n").split("\t")
                print(c, s, e, ".", 0, "+", sep="\t", file=fo)
        with open(tmp_minus_m) as fm:
            for ln in fm:
                c, s, e = ln.rstrip("\n").split("\t")
                print(c, s, e, ".", 0, "-", sep="\t", file=fo)

    bedtools_sort(tmp_union, out_union_bed)
    for p in [tmp_plus, tmp_minus, tmp_plus_m, tmp_minus_m, tmp_union]:
        try: os.remove(p)
        except: pass

def build_union_bed_nostrand(sorted_ann_bed, out_union_bed):
    """
    构建“无链并集”（用于 overlap_bp 的唯一覆盖计算，防止同/反链双计）。
    输入：注释逐Tx BED8
    输出：BED6：chrom start end name '.' strand('.')   —— strand 置 '.'
    """
    threecol = out_union_bed + ".3col"
    merged = out_union_bed + ".merged"
    with open(sorted_ann_bed) as fi, open(threecol, "w") as fo:
        for ln in fi:
            p = ln.rstrip("\n").split("\t")
            print(p[0], p[1], p[2], sep="\t", file=fo)
    with open(merged, "w") as fo:
        subprocess.check_call(["bedtools", "merge", "-i", threecol], stdout=fo)
    with open(merged) as fi, open(out_union_bed, "w") as fo:
        for ln in fi:
            c, s, e = ln.rstrip("\n").split("\t")
            print(c, s, e, ".", 0, ".", sep="\t", file=fo)
    os.remove(threecol); os.remove(merged)

def bedtools_intersect(a_bed, b_bed, out_path, same_strand=True):
    """
    -sorted -wo；same_strand=True -> -s；False -> -S；None -> 不加 -s/-S
    """
    args = ["bedtools", "intersect", "-wo", "-sorted", "-a", a_bed, "-b", b_bed]
    if same_strand is True:
        args.insert(2, "-s")
    elif same_strand is False:
        args.insert(2, "-S")
    # None: 不插入 -s/-S
    with open(out_path, "w") as fo:
        subprocess.check_call(args, stdout=fo)

def agg_union_overlap(union_wo_path):
    """
    读取 “new_bed vs 注释并集 -wo” 的结果，按 ORF_id 汇总 overlap_len。
    A（new_bed）是8列，id在 A.col4
    """
    sums = defaultdict(int)
    with open(union_wo_path) as f:
        for ln in f:
            parts = ln.rstrip("\n").split("\t")
            orf_id = parts[3]
            ov = int(parts[-1])
            sums[orf_id] += ov
    return sums

def merged_len_per_orf(interval_bed):
    """
    对“in-frame区间（chrom start end ORF_id）”按 ORF 逐条流式合并，输出每个 ORF 的不重叠长度。
    依赖：coreutils sort
    """
    if (not os.path.exists(interval_bed)) or os.path.getsize(interval_bed) == 0:
        return {}
    sorted_bed = interval_bed + ".sorted"
    subprocess.check_call(["sort", "-k4,4", "-k1,1", "-k2,2n", interval_bed, "-o", sorted_bed])

    res = {}
    curr_orf = None
    curr_chr = None
    cur_s = -1
    cur_e = -1
    acc = 0

    with open(sorted_bed) as f:
        for ln in f:
            c, s, e, oid = ln.rstrip("\n").split("\t")
            s = int(s); e = int(e)
            if oid != curr_orf:
                if curr_orf is not None:
                    acc += max(0, cur_e - cur_s)
                    res[curr_orf] = res.get(curr_orf, 0) + acc
                curr_orf, curr_chr, cur_s, cur_e, acc = oid, c, s, e, 0
            else:
                if c != curr_chr or s > cur_e:
                    acc += max(0, cur_e - cur_s)
                    curr_chr, cur_s, cur_e = c, s, e
                else:
                    cur_e = max(cur_e, e)

    if curr_orf is not None:
        acc += max(0, cur_e - cur_s)
        res[curr_orf] = res.get(curr_orf, 0) + acc

    os.remove(sorted_bed)
    return res

# ---------- 同链 in-frame 判定 + 计数 ----------

def agg_tx_gene_and_inframe(tx_wo_path, dump_inframe_bed=None):
    """
    读取 “new_bed vs ann_bed（逐Tx） -s -wo”：
      - 计算 in-frame 区间（一次判定整段），可选写出为 BED4：chrom start end ORF_id
      - 统计同链总重叠（按Tx/Gene聚合）用于 top_hit、n_tx、n_gene

    A(new) 8列： 0 chromA,1 startA,2 endA,3 idA,4 score,5 strandA,6 phaseA,7 geneA
    B(ann) 8列： 8 chromB,9 startB,10 endB,11 txB,12 score,13 strandB,14 phaseB,15 geneB
    最末列：overlap_len
    """
    if dump_inframe_bed:
        fo_inf = open(dump_inframe_bed, "w")
    else:
        fo_inf = None

    sum_tx = defaultdict(lambda: defaultdict(int))    # ORF -> {tx: ov_bp}
    sum_gene = defaultdict(lambda: defaultdict(int))  # ORF -> {gene: ov_bp}

    with open(tx_wo_path) as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            # A
            chromA = p[0]
            sA = int(p[1]); eA = int(p[2])
            oid = p[3]
            strandA = p[5]
            phA = int(p[6])
            # B
            sB = int(p[9]); eB = int(p[10])
            tx = p[11]
            strandB = p[13]
            phB = int(p[14])
            geneB = p[15]
            # overlap
            ov = int(p[-1])

            # 交叠起止（0-based half-open）
            s0 = max(sA, sB)
            e0 = min(eA, eB)
            if e0 <= s0:
                continue

            # in-frame 判定（交叠起点一次判定整段）
            if strandA == '+':
                lfA = (phA + (s0 - sA)) % 3
                lfB = (phB + (s0 - sB)) % 3
                is_in = (lfA == lfB)
            else:  # '-'
                rfA = (phA + (eA - 1 - s0)) % 3
                rfB = (phB + (eB - 1 - s0)) % 3
                is_in = (rfA == rfB)

            # 同链重叠用于计数/Top
            sum_tx[oid][tx] += ov
            sum_gene[oid][geneB] += ov

            # in-frame 区间写出，后续按 ORF 合并去重计长度
            if is_in and fo_inf is not None:
                print(chromA, s0, e0, oid, sep="\t", file=fo_inf)

    if fo_inf:
        fo_inf.close()

    n_tx = {oid: len(d) for oid, d in sum_tx.items()}
    n_gene = {oid: len(d) for oid, d in sum_gene.items()}

    top_tx = {}
    for oid, d in sum_tx.items():
        top_tx[oid] = max(d.items(), key=lambda x: (x[1], x[0]))[0] if d else ""

    top_gene = {}
    for oid, d in sum_gene.items():
        top_gene[oid] = max(d.items(), key=lambda x: (x[1], x[0]))[0] if d else ""

    return n_tx, n_gene, top_gene, top_tx

# ---------- 主流程 ----------

def main():
    ap = argparse.ArgumentParser(description="新ORF与注释CDS重叠与in-frame统计（去重版）")
    ap.add_argument("--new-gtf", required=True, help="新ORF的GTF（含CDS，建议有ORF_id属性）")
    ap.add_argument("--ann-gtf", required=True, help="注释GTF（如gencode，取CDS）")
    ap.add_argument("--out", required=True, help="输出TSV路径")
    ap.add_argument("--new-id-attr", default="ORF_id", help="新ORF的ID属性名（默认ORF_id；若无则回退到transcript_id）")
    ap.add_argument("--ann-tx-attr", default="transcript_id", help="注释转录本ID属性名")
    ap.add_argument("--ann-gene-attr", default="gene_id", help="注释基因ID属性名")
    ap.add_argument("--tmpdir", default=None, help="可选的临时目录（默认系统临时目录）")
    args = ap.parse_args()

    check_exec("bedtools")

    tmp = tempfile.mkdtemp(prefix="orf_overlaps_", dir=args.tmpdir)
    try:
        new_bed = os.path.join(tmp, "new.bed8")            # BED8：bed6 + phase + gene
        new_len = os.path.join(tmp, "new.len.tsv")
        ann_bed = os.path.join(tmp, "ann_tx.bed8")         # BED8：bed6 + phase + gene
        new_bed_sorted = new_bed + ".sorted"
        ann_bed_sorted = ann_bed + ".sorted"

        # 1) GTF -> BED（重建phase）
        orf_len = gtf_to_cds_bed_with_phase(
            args.new_gtf, args.new_id_attr, "transcript_id", "gene_id",
            new_bed, lengths_out=new_len)
        gtf_to_cds_bed_with_phase(
            args.ann_gtf, args.ann_tx_attr, args.ann_tx_attr, args.ann_gene_attr,
            ann_bed)

        # 2) sort
        bedtools_sort(new_bed, new_bed_sorted)
        bedtools_sort(ann_bed, ann_bed_sorted)

        # 3) 注释CDS“无链并集”（用于 overlap_bp 的唯一覆盖）
        ann_union_nostrand = os.path.join(tmp, "ann_union_nostrand.bed6")
        build_union_bed_nostrand(ann_bed_sorted, ann_union_nostrand)

        # 4) new vs ann_union_nostrand（不加 -s/-S），得到唯一 overlap_bp
        unique_union_wo = os.path.join(tmp, "unique_union.wo")
        bedtools_intersect(new_bed_sorted, ann_union_nostrand, unique_union_wo, same_strand=None)
        unique_sum = agg_union_overlap(unique_union_wo)

        # 5) 同链逐Tx交叠 -> in-frame 区间 + n_tx/n_gene/top
        tx_same_wo = os.path.join(tmp, "tx_same.wo")
        bedtools_intersect(new_bed_sorted, ann_bed_sorted, tx_same_wo, same_strand=True)

        inframe_intervals_bed = os.path.join(tmp, "inframe_intervals.bed")
        n_tx, n_gene, top_gene, top_tx = agg_tx_gene_and_inframe(tx_same_wo, dump_inframe_bed=inframe_intervals_bed)

        # 6) in-frame 区间按 ORF 合并去重计长度
        inframe_uniq = merged_len_per_orf(inframe_intervals_bed)

        # 7) 汇总输出
        df = pd.DataFrame({"ORF_id": list(orf_len.keys()), "ORF_len": list(orf_len.values())})
        orf_ids = df["ORF_id"].to_numpy()

        overlap_bp = np.array([unique_sum.get(oid, 0) for oid in orf_ids], dtype=np.int64)
        inframe_bp = np.array([inframe_uniq.get(oid, 0) for oid in orf_ids], dtype=np.int64)
        n_tx_arr = np.array([n_tx.get(oid, 0) for oid in orf_ids], dtype=np.int32)
        n_gene_arr = np.array([n_gene.get(oid, 0) for oid in orf_ids], dtype=np.int32)
        top_gene_arr = np.array([top_gene.get(oid, "") for oid in orf_ids], dtype=object)
        top_tx_arr = np.array([top_tx.get(oid, "") for oid in orf_ids], dtype=object)

        orf_len_arr = df["ORF_len"].to_numpy(dtype=np.int64)
        with np.errstate(divide='ignore', invalid='ignore'):
            overlap_fraction = np.where(orf_len_arr > 0, overlap_bp / orf_len_arr, 0.0)
            inframe_fraction = np.where(orf_len_arr > 0, inframe_bp / orf_len_arr, 0.0)

        out_df = pd.DataFrame({
            "ORF_id": orf_ids,
            "overlap_bp": overlap_bp,
            "overlap_fraction": overlap_fraction,
            "inframe_overlap_bp": inframe_bp,
            "inframe_fraction": inframe_fraction,
            "n_annot_tx_overlapped": n_tx_arr,
            "n_annot_gene_overlapped": n_gene_arr,
            "top_hit_gene": top_gene_arr,
            "top_hit_tx": top_tx_arr
        })
        out_df.to_csv(args.out, sep="\t", index=False)
        print(f"[OK] 写出结果：{args.out}")

    finally:
        shutil.rmtree(tmp, ignore_errors=True)

if __name__ == "__main__":
    main()
