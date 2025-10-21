#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, sys
import pandas as pd

CAND5 = ["codon5","start_codon","codon_5","codon_start","start","start_nt"]
CAND3 = ["codon3","stop_codon","codon_3","codon_stop","end","end_nt"]

def guess_col(cols, candidates):
    lower_map = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    for c in cols:
        cl = c.lower()
        for cand in candidates:
            if cand.lower() in cl:
                return c
    return None

def split_orfid_series(s: pd.Series):
    """按 [|:] 拆分旧 orfID -> (ID, Chr, Strand)"""
    s = s.astype("string")
    parts = s.str.split(r"[|:]", n=8, expand=True)
    for k in range(3):
        if k not in parts.columns:
            parts[k] = pd.Series(pd.NA, index=s.index, dtype="string")
    ID = parts[0].astype("string")
    Chr = parts[1].astype("string")
    Strand = parts[2].astype("string")
    return ID, Chr, Strand

def build_orf_id_custom(df, old_col, codon5_col, codon3_col, aa_col):
    # 旧ID
    old = df[old_col].astype("string").str.strip()
    # 排除以 sp| 开头
    not_sp = ~old.str.startswith(("sp|","tr|"), na=False)

    # 拆分
    ID, Chr, Strand = split_orfid_series(old)

    # 取列
    c5 = df[codon5_col].astype("string")
    c3 = df[codon3_col].astype("string")
    aa = df[aa_col].astype("string")

    # 条件：非 sp 且各字段非缺失
    mask = not_sp & ID.notna() & Chr.notna() & Strand.notna() & c5.notna() & c3.notna() & aa.notna()

    out = pd.Series(pd.NA, index=df.index, dtype="string")
    out.loc[mask] = (
        Strand[mask] + Chr[mask] + ":" + c5[mask] + "-" + c3[mask] + ":" + aa[mask]
    )
    return out, (~not_sp).sum()

def process_once(in_path, out_path, old_col, codon5, codon3, aa_col):
    df = pd.read_csv(in_path, sep="\t", dtype="string", na_filter=False)
    for col in [old_col, aa_col]:
        if col not in df.columns:
            sys.exit(f"[ERROR] 列不存在：{col}")

    if codon5 is None:
        codon5 = guess_col(df.columns, CAND5)
    if codon3 is None:
        codon3 = guess_col(df.columns, CAND3)
    if codon5 is None or codon3 is None:
        sys.exit("[ERROR] 无法识别 codon5/codon3 列名，请用 --codon5/--codon3 指定。")

    df["ORF_id_custom"], sp_excluded = build_orf_id_custom(df, old_col, codon5, codon3, aa_col)
    df.to_csv(out_path, sep="\t", index=False)

    total = len(df)
    ok = df["ORF_id_custom"].notna().sum()
    print(f"[OK] rows={total}, built={ok}, miss={total-ok}, sp_excluded={sp_excluded}; codon5={codon5}, codon3={codon3}, aa={aa_col}")

def process_chunked(in_path, out_path, old_col, codon5, codon3, aa_col, chunksize):
    it = pd.read_csv(in_path, sep="\t", dtype="string", na_filter=False, chunksize=chunksize)
    header_written = False
    total = ok = sp_excluded_total = 0

    for i, df in enumerate(it, 1):
        for col in [old_col, aa_col]:
            if col not in df.columns:
                sys.exit(f"[ERROR] 列不存在：{col}")

        if i == 1:
            if codon5 is None:
                codon5 = guess_col(df.columns, CAND5)
            if codon3 is None:
                codon3 = guess_col(df.columns, CAND3)
            if codon5 is None or codon3 is None:
                sys.exit("[ERROR] 无法识别 codon5/codon3 列名，请用 --codon5/--codon3 指定。")
            print(f"[INFO] Use codon5={codon5}, codon3={codon3}, aa={aa_col}")

        df["ORF_id_custom"], sp_excluded = build_orf_id_custom(df, old_col, codon5, codon3, aa_col)
        df.to_csv(out_path, sep="\t", index=False, header=not header_written, mode=("a" if header_written else "w"))
        header_written = True

        total += len(df)
        ok += df["ORF_id_custom"].notna().sum()
        sp_excluded_total += sp_excluded
        if i % 10 == 0:
            print(f"[INFO] processed {total} rows...")

    print(f"[OK] rows={total}, built={ok}, miss={total-ok}, sp_excluded={sp_excluded_total}")

def main():
    ap = argparse.ArgumentParser(description="从旧 orfID 生成 ORF_id_custom（忽略以 sp| 开头的 ID）")
    ap.add_argument("--in",  dest="in_path",  required=True, help="输入TSV（支持.gz）")
    ap.add_argument("--out", dest="out_path", required=True, help="输出TSV（支持.gz）")
    ap.add_argument("--old-col", required=True, help="旧 orfID 的列名")
    ap.add_argument("--aa-col",  required=True, help="氨基酸序列列名（将拼到新ID末尾）")
    ap.add_argument("--codon5", default=None, help="codon5列名（可不填，自动猜）")
    ap.add_argument("--codon3", default=None, help="codon3列名（可不填，自动猜）")
    ap.add_argument("--chunksize", type=int, default=0, help="分块大小（行数），>0 启用分块")
    args = ap.parse_args()

    if args.chunksize and args.chunksize > 0:
        process_chunked(args.in_path, args.out_path, args.old_col, args.codon5, args.codon3, args.aa_col, args.chunksize)
    else:
        process_once(args.in_path, args.out_path, args.old_col, args.codon5, args.codon3, args.aa_col)

if __name__ == "__main__":
    main()
