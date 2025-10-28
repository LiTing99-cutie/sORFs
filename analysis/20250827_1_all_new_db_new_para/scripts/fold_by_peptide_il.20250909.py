#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

def main():
    ap = argparse.ArgumentParser(description="按 Peptide_I_L_equal 聚合并汇总 Source。")
    ap.add_argument("-i", "--input", required=True, help="输入TSV文件路径（含表头）")
    ap.add_argument("-o", "--output", required=True, help="输出TSV文件路径")
    ap.add_argument("--sep", default="\t", help="输入分隔符，默认tab")
    args = ap.parse_args()

    # 读取
    df = pd.read_csv(args.input, sep=',', dtype=str, keep_default_na=False)
    required_cols = ["Peptide_I_L_equal", "Source"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"缺少必要列: {missing}")

    # 清理
    df["Peptide_I_L_equal"] = df["Peptide_I_L_equal"].astype(str).str.strip()
    df["Source"] = df["Source"].astype(str).str.strip()

    # 去重保序拼接
    def join_unique_preserve_order(series):
        seen, out = set(), []
        for x in series:
            x = x.strip()
            if x and x not in seen:
                seen.add(x)
                out.append(x)
        return ",".join(out)

    # 分组聚合
    out = (
        df.groupby("Peptide_I_L_equal", sort=False)
          .agg(
              Source=("Source", join_unique_preserve_order),
              Source_number=("Source", lambda s: pd.Series([x.strip() for x in s if x.strip()]).nunique())
          )
          .reset_index()
    )

    # 保存为TSV
    out.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
