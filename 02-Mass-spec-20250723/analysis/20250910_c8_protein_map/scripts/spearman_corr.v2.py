#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

def to_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")

def to_bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s
    return s.astype(str).str.strip().str.lower().isin({"true","1","yes","y","t"})

def spearman_pair(df: pd.DataFrame, x: str, y: str, min_n: int = 3):
    sub = df[[x, y]].apply(to_num).dropna()
    n = len(sub)
    if n < min_n:
        return np.nan, n
    rho = sub.corr(method="spearman").iloc[0, 1]
    return float(rho), n

def main():
    ap = argparse.ArgumentParser(description="在合并表上按 canonical 分组计算 Spearman 相关")
    ap.add_argument("--merged", required=True, help="合并结果 .tsv（来自 merge_tables.py）")
    ap.add_argument("--out",    required=True, help="输出相关性结果 .tsv")
    ap.add_argument("--min-n",  type=int, default=3, help="每组最小样本数（默认3）")
    # 可自定义要评估的列对；默认与你原脚本一致
    ap.add_argument("--pairs", nargs="*", default=[
        "FL_TPM:RPF_RPKM",
        "A:RPF_RPKM",
        "A:mean_relative_iBAQ_B",
        "RPF_RPKM:mean_relative_iBAQ_B",
        "Psites_RPKM:mean_relative_iBAQ_B",
    ], help="以 X:Y 的形式给出列对，空格分隔")
    args = ap.parse_args()

    df = pd.read_csv(args.merged, sep="\t", low_memory=False)

    # 识别 is_canonical 列（大小写均可）；若缺失则报错更稳妥
    canon_col = next((c for c in df.columns if c.lower() == "is_canonical"), None)
    if canon_col is None:
        raise SystemExit("输入表缺少 is_canonical/Is_canonical 列，无法分组。")
    is_canon = to_bool_series(df[canon_col])

    # 准备三组：all / canonical / noncanonical
    groups = {
        "all": df,
        "canonical": df.loc[is_canon].copy(),
        "noncanonical": df.loc[~is_canon].copy(),
    }

    # 解析列对
    pairs = []
    for token in args.pairs:
        if ":" not in token:
            raise SystemExit(f"列对格式应为 X:Y，但收到：{token}")
        x, y = token.split(":", 1)
        pairs.append((x, y))

    # 逐组逐对计算
    records = []
    for gname, gdf in groups.items():
        for x, y in pairs:
            # 列存在性检查
            missing = [c for c in (x, y) if c not in gdf.columns]
            if missing:
                records.append({
                    "group": gname, "var_x": x, "var_y": y,
                    "spearman_rho": np.nan, "n": 0,
                    "note": f"missing cols: {','.join(missing)}"
                })
                continue

            rho, n = spearman_pair(gdf, x, y, min_n=args.min_n)
            records.append({
                "group": gname, "var_x": x, "var_y": y,
                "spearman_rho": rho, "n": n, "note": ""
            })

    out_df = pd.DataFrame.from_records(records,
                                       columns=["group","var_x","var_y","spearman_rho","n","note"])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(out_df)

if __name__ == "__main__":
    main()
