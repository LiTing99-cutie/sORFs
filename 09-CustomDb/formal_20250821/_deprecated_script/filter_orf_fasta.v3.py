#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, argparse, gzip, time

def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode, encoding="utf-8", errors="ignore")

def parse_ids_by_threshold(tsv_path, col_keep, thr, cmp_op, id_col="ORF_id", log_every=2_000_000):
    """
    从TSV按阈值筛选ID（单列条件）。
    cmp_op: "ge" -> >=, "eq0" -> 近似等于0（≤1e-12）
    仅保留列: id_col, col_keep
    """
    t0 = time.time()
    passed = set()
    n = 0
    with open_maybe_gzip(tsv_path, "rt") as f:
        header = f.readline()
        if not header:
            return passed
        cols = header.rstrip("\n").split("\t")
        try:
            i_id = cols.index(id_col)
            i_val = cols.index(col_keep)
        except ValueError as e:
            sys.exit(f"[ERROR] {tsv_path}: 缺少列 {e}")
        for line in f:
            n += 1
            if log_every and n % log_every == 0:
                sys.stderr.write(f"[LOG] scanning {os.path.basename(tsv_path)}: {n/1e6:.1f}M rows, kept={len(passed)}\n")
            p = line.rstrip("\n").split("\t")
            if len(p) <= i_val or len(p) <= i_id: 
                continue
            oid = p[i_id]
            try:
                v = float(p[i_val])
            except:
                continue
            if cmp_op == "ge":
                if v >= thr:
                    passed.add(oid)
            elif cmp_op == "eq0":
                if abs(v) <= 1e-12:
                    passed.add(oid)
            else:
                raise ValueError("unsupported cmp_op")
    sys.stderr.write(f"[DONE] {os.path.basename(tsv_path)} scanned in {time.time()-t0:.1f}s, kept={len(passed)}\n")
    return passed

def write_fasta_subset(fa_in, fa_out, keep_ids, log_every=200000):
    """
    从fa_in中导出首字段ID在keep_ids的条目；保留原始描述。
    """
    t0 = time.time()
    n_total = n_kept = 0
    with open_maybe_gzip(fa_in, "rt") as fi, open(fa_out, "w") as fo:
        write = False
        for line in fi:
            if line.startswith(">"):
                n_total += 1
                seq_id = line[1:].strip().split()[0]
                write = (seq_id in keep_ids)
                if write:
                    n_kept += 1
                    fo.write(line)
                # else: skip header and sequence lines until next header
            else:
                if write:
                    fo.write(line)
            if log_every and n_total % log_every == 0 and line.startswith(">"):
                sys.stderr.write(f"[LOG] FASTA parsed {n_total} records, kept={n_kept}\n")
    sys.stderr.write(f"[DONE] FASTA subset written: {fa_out} | total={n_total}, kept={n_kept}, took {time.time()-t0:.1f}s\n")
    missing = len(keep_ids) - n_kept
    if missing > 0:
        sys.stderr.write(f"[WARN] {missing} IDs通过条件但未在FASTA中找到（可能ID不一致或未入库）。\n")

def main():
    ap = argparse.ArgumentParser(description="按RPF阈值筛选ORF并导出FASTA")
    ap.add_argument("--candidate-fa", required=True, help="候选ORF FASTA（pep.fa）")
    ap.add_argument("--psite-stats", required=True, help="psite_frame_stats.v2.tsv")
    ap.add_argument("--rpf-psite", required=True, help="orf.rpf.psite.txt")
    ap.add_argument("--overlap", required=True, help="orf_overlap_inframe.txt")
    ap.add_argument("--out-fa", required=True, help="输出FASTA")
    ap.add_argument("--f0-thr", type=float, default=0.5, help="frame0_fraction 阈值（≥）")
    ap.add_argument("--cov-thr", type=float, default=0.1, help="Psites_codon_coverage 阈值（≥）")
    ap.add_argument("--log-every", type=int, default=2_000_000, help="每处理N行打印一次日志")
    args = ap.parse_args()

    # 条件1：frame0_fraction ≥ f0_thr
    s1 = parse_ids_by_threshold(
        args.psite_stats, col_keep="frame0_fraction", thr=args.f0_thr, cmp_op="ge",
        id_col="ORF_id", log_every=args.log_every
    )
    # 条件2：Psites_codon_coverage ≥ cov_thr
    s2 = parse_ids_by_threshold(
        args.rpf_psite, col_keep="Psites_codon_coverage", thr=args.cov_thr, cmp_op="ge",
        id_col="ORF_id", log_every=args.log_every
    )
    # 条件3：inframe_fraction == 0
    s3 = parse_ids_by_threshold(
        args.overlap, col_keep="inframe_fraction", thr=0.0, cmp_op="eq0",
        id_col="ORF_id", log_every=args.log_every
    )

    # 交集
    inter = s1 & s2 & s3
    sys.stderr.write(f"[SUM] pass_f0={len(s1)}, pass_cov={len(s2)}, pass_inframe0={len(s3)}, FINAL={len(inter)}\n")

    if not inter:
        sys.stderr.write("[WARN] 最终通过集合为空，未生成FASTA。\n")
        # 仍然写一个空文件以示完成
        open(args.out_fa, "w").close()
        return

    # 导出FASTA子集
    write_fasta_subset(args.candidate_fa, args.out_fa, inter)

if __name__ == "__main__":
    main()
