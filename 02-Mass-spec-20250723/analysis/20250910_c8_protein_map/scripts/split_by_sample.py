#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import argparse
from pathlib import Path

def detect_sample_col(df):
    for c in df.columns:
        if c.lower() == "sample":
            return c
    raise SystemExit("未找到 sample 列（大小写皆可）。")

def split_one(df, sample_col, out_dir, out_name):
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) 先清洗样本名（string dtype 能保持缺失为 <NA>）
    s = df[sample_col].astype("string").str.strip()

    # 只保留有效样本名：非缺失且非空
    valid_mask = s.notna() & (s != "")
    valid_samples = sorted(s[valid_mask].unique().tolist())

    # 2) 逐样本切分（用清洗后的 s 做筛选，避免 "nan" 被当成样本）
    for samp in valid_samples:
        sub = df[valid_mask & (s == samp)].copy()
        sd = out_dir / samp
        sd.mkdir(parents=True, exist_ok=True)
        sub.to_csv(sd / out_name, sep="\t", index=False)

    # 3) 提示被丢弃的行（sample 缺失/空白）
    dropped = (~valid_mask).sum()
    if dropped:
        print(f"[WARN] {out_name}: 丢弃 sample 缺失/空白的行数 = {dropped}")

    return valid_samples

def main():
    ap = argparse.ArgumentParser(description="按 sample 列拆分两个 merged 输入为每样本文件，并输出样本清单")
    ap.add_argument("--pep_orf_merged", required=True,
                    help=".../MS_res_from_Galaxy/pep.orf.merged.txt")
    ap.add_argument("--intensity_merged", required=True,
                    help=".../MS_res_from_Galaxy/peptide_intensity_IL.merged.tsv")
    ap.add_argument("--out_base", required=True, help="输出基目录，例如 ../processed/by_sample")
    ap.add_argument("--sample_list", required=True, help="输出样本清单 txt")
    args = ap.parse_args()

    out_base = Path(args.out_base)
    out_base.mkdir(parents=True, exist_ok=True)

    pep = pd.read_csv(args.pep_orf_merged, sep="\t", low_memory=False)
    pint = pd.read_csv(args.intensity_merged, sep="\t", low_memory=False)

    col_p = detect_sample_col(pep)
    col_i = detect_sample_col(pint)

    sp1 = split_one(pep,  col_p, out_base, "pep.orf.txt")
    sp2 = split_one(pint, col_i, out_base, "peptide_intensity_IL.tsv")

    # 取交集，提醒差异（通常应一致）
    s1, s2 = set(sp1), set(sp2)
    samples = sorted(s1 | s2)
    (Path(args.sample_list)).write_text("\n".join(samples) + "\n", encoding="utf-8")

    only1 = sorted(s1 - s2)
    only2 = sorted(s2 - s1)
    if only1:
        print(f"[WARN] 仅 pep.orf.merged 中存在的样本: {only1}")
    if only2:
        print(f"[WARN] 仅 intensity.merged 中存在的样本: {only2}")
    print(f"[OK] 拆分完成。样本数={len(samples)}，清单输出：{args.sample_list}")

if __name__ == "__main__":
    main()
