#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, argparse, gzip, time

def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode, encoding="utf-8", errors="ignore")

def parse_ids_by_threshold(tsv_path, col_keep, thr, id_col="ORF_id", log_every=2_000_000):
    """从TSV按“>= thr”筛ID（单列），流式读取。"""
    t0 = time.time()
    passed = set(); n = 0
    with open_maybe_gzip(tsv_path, "rt") as f:
        header = f.readline()
        if not header:
            return passed
        cols = header.rstrip("\n").split("\t")
        try:
            i_id = cols.index(id_col); i_val = cols.index(col_keep)
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
            if v >= thr:
                passed.add(oid)
    sys.stderr.write(f"[DONE] {os.path.basename(tsv_path)} scanned in {time.time()-t0:.1f}s, kept={len(passed)}\n")
    return passed

def parse_inframe_with_override(tsv_path, candidate_ids, canonical_keyword="canonical", id_col="ORF_id",
                                col_inframe="inframe_fraction", tol=1e-12, log_every=2_000_000):
    """
    对候选集合candidate_ids：
      - 若 ORF_id 含 canonical_keyword -> 直接通过
      - 否则要求 inframe_fraction == 0（|v| <= tol）
    返回通过集合
    """
    t0 = time.time()
    passed = {oid for oid in candidate_ids if canonical_keyword in oid}
    need_check = {oid for oid in candidate_ids if canonical_keyword not in oid}
    if not need_check:
        sys.stderr.write("[INFO] 所有候选均为canonical，inframe过滤跳过。\n")
        return passed

    n = 0
    hit = 0
    with open_maybe_gzip(tsv_path, "rt") as f:
        header = f.readline()
        if not header:
            return passed
        cols = header.rstrip("\n").split("\t")
        try:
            i_id = cols.index(id_col); i_val = cols.index(col_inframe)
        except ValueError as e:
            sys.exit(f"[ERROR] {tsv_path}: 缺少列 {e}")
        for line in f:
            n += 1
            if log_every and n % log_every == 0:
                sys.stderr.write(f"[LOG] scanning {os.path.basename(tsv_path)} (inframe): {n/1e6:.1f}M rows\n")
            p = line.rstrip("\n").split("\t")
            if len(p) <= i_val or len(p) <= i_id:
                continue
            oid = p[i_id]
            if oid not in need_check:
                continue
            try:
                v = float(p[i_val])
            except:
                continue
            if abs(v) <= tol:
                passed.add(oid)
                hit += 1
    sys.stderr.write(f"[DONE] {os.path.basename(tsv_path)} inframe checked in {time.time()-t0:.1f}s; "
                     f"canonical_pass={len(passed)-hit}, eq0_pass={hit}\n")
    return passed

def write_fasta_subset(fa_in, fa_out, keep_ids, log_every=200000):
    """从FASTA导出ID在keep_ids的序列；ID取'>'后的首字段。"""
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
    ap = argparse.ArgumentParser(description="按RPF阈值筛选ORF并导出FASTA（canonical跳过inframe==0过滤）")
    ap.add_argument("--candidate-fa", required=True, help="候选ORF FASTA（pep.fa，可.gz）")
    ap.add_argument("--psite-stats", required=True, help="psite_frame_stats.v2.tsv")
    ap.add_argument("--rpf-psite", required=True, help="orf.rpf.psite.txt")
    ap.add_argument("--overlap", required=True, help="orf_overlap_inframe.txt")
    ap.add_argument("--out-fa", required=True, help="输出FASTA")
    ap.add_argument("--f0-thr", type=float, default=0.5, help="frame0_fraction 阈值（≥）")
    ap.add_argument("--cov-thr", type=float, default=0.1, help="Psites_codon_coverage 阈值（≥）")
    ap.add_argument("--canonical-key", default="canonical", help="canonical 关键字（默认 'canonical'）")
    ap.add_argument("--log-every", type=int, default=2_000_000, help="每处理N行打印一次日志")
    args = ap.parse_args()

    # 条件1：frame0_fraction ≥ f0_thr
    s1 = parse_ids_by_threshold(
        args.psite_stats, col_keep="frame0_fraction", thr=args.f0_thr,
        id_col="ORF_id", log_every=args.log_every
    )
    # 条件2：Psites_codon_coverage ≥ cov_thr
    s2 = parse_ids_by_threshold(
        args.rpf_psite, col_keep="Psites_codon_coverage", thr=args.cov_thr,
        id_col="ORF_id", log_every=args.log_every
    )
    cand = s1 & s2
    sys.stderr.write(f"[SUM] pass_f0={len(s1)}, pass_cov={len(s2)}, A∩B candidates={len(cand)}\n")

    if not cand:
        sys.stderr.write("[WARN] A∩B为空；写空FASTA。\n")
        open(args.out_fa, "w").close()
        return

    # 条件3：inframe_fraction==0（非canonical），canonical则直接通过
    s3 = parse_inframe_with_override(
        args.overlap, candidate_ids=cand, canonical_keyword=args.canonical_key,
        id_col="ORF_id", col_inframe="inframe_fraction", tol=1e-12,
        log_every=args.log_every
    )

    final_ids = cand & s3
    sys.stderr.write(f"[SUM] FINAL={len(final_ids)} (after canonical-aware inframe filter)\n")

    # 导出FASTA
    write_fasta_subset(args.candidate_fa, args.out_fa, final_ids)

if __name__ == "__main__":
    main()
