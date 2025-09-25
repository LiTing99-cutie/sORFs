#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

def to_num(s):
    return pd.to_numeric(s, errors="coerce")

def spearman_pair(df, x, y):
    sub = df[[x, y]].apply(to_num).dropna()
    if len(sub) < 3:
        return {"var_x": x, "var_y": y, "spearman_rho": np.nan, "n": len(sub)}
    rho = sub.corr(method="spearman").iloc[0,1]
    return {"var_x": x, "var_y": y, "spearman_rho": float(rho), "n": len(sub)}

def main():
    ap = argparse.ArgumentParser(description="在合并表上计算 Spearman 相关系数")
    ap.add_argument("--merged", required=True, help="合并结果 .tsv（来自 merge_tables.py）")
    ap.add_argument("--out",    required=True, help="输出相关性结果 .tsv")
    args = ap.parse_args()

    df = pd.read_csv(args.merged, sep="\t", low_memory=False)

    pairs = [
        ("FL_TPM", "RPF_RPKM"),
        ("A", "RPF_RPKM"),
        ("A", "iBAQ_B"),
        ("RPF_RPKM", "iBAQ_B"),
    ]
    results = []
    for x, y in pairs:
        for col in (x, y):
            if col not in df.columns:
                print(f"[WARN] 列缺失：{col}；该对将输出 NaN")
        results.append(spearman_pair(df, x, y))

    out_df = pd.DataFrame(results)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(out_df)

if __name__ == "__main__":
    main()
