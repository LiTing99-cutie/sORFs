#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def main():
    ap = argparse.ArgumentParser(
        description="按样本归一(relative)后，跨样本对每个 ORF 求均值；度量可选 iBAQ_A 或 iBAQ_B"
    )
    ap.add_argument("--in", dest="in_tsv", required=True,
                    help="输入表（需含列：ORF_id, Sample, 以及 iBAQ_A 或 iBAQ_B）")
    ap.add_argument("--metric", choices=["iBAQ_A", "iBAQ_B"], default="iBAQ_B",
                    help="选择用于计算的度量列（默认 iBAQ_B）")
    ap.add_argument("--out", dest="out_tsv", required=True,
                    help="输出两列：ORF_id, mean_relative_iBAQ")
    args = ap.parse_args()

    usecols = ["ORF_id", "Sample", args.metric]
    try:
        df = pd.read_csv(args.in_tsv, sep="\t", usecols=usecols, engine="pyarrow")
    except Exception:
        df = pd.read_csv(args.in_tsv, sep="\t", usecols=usecols, low_memory=False)

    # 列检查
    for c in ["ORF_id", "Sample", args.metric]:
        if c not in df.columns:
            raise SystemExit(f"缺少必要列：{c}")

    # 数值化度量
    df[args.metric] = pd.to_numeric(df[args.metric], errors="coerce")

    # 仅纳入该度量有效且 >0 的记录
    m = df[args.metric].notna() & (df[args.metric] > 0)
    dfv = df.loc[m].copy()
    if dfv.empty:
        out = pd.DataFrame(columns=["ORF_id", "mean_relative_iBAQ"])
        Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out_tsv, sep="\t", index=False)
        print("[WARN] 无有效数值，输出空表。")
        return

    # 每样本总量（避免 merge，用 transform），并防除零
    sum_per_sample = dfv.groupby("Sample")[args.metric].transform("sum")
    sum_per_sample = sum_per_sample.where(sum_per_sample > 0, np.nan)

    # 样本内相对值
    dfv["relative"] = dfv[args.metric] / sum_per_sample

    # 跨样本对每个 ORF 求均值（仅对该样本里 relative 有值者计入）
    out = (dfv.dropna(subset=["relative"])
              .groupby("ORF_id")["relative"]
              .mean()
              .rename("mean_relative_iBAQ")
              .reset_index())

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[OK] metric={args.metric}  输出：{args.out_tsv}  shape={out.shape}")

if __name__ == "__main__":
    try:
        pd.options.mode.copy_on_write = True
    except Exception:
        pass
    main()
