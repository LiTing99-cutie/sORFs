#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def to_bool_series(s: pd.Series) -> pd.Series:
    """将多种表示转成布尔：True/1/yes/y/t → True；其余 False"""
    if pd.api.types.is_bool_dtype(s):
        return s
    return s.astype(str).str.strip().str.lower().isin({"true","1","yes","y","t"})

def main():
    ap = argparse.ArgumentParser(
        description="共同锚定集(MoR)归一化 → 直接跨 run 聚合（不再做每个 run 的“除以总和”）"
    )
    ap.add_argument("--in", dest="in_tsv", required=True,
                    help="all_samples.ibaq_b_with_total.tsv（需含 Sample,iBAQ_B,ORF_id；可含 Available,is_canonical）")
    ap.add_argument("--out", dest="out_tsv", required=True,
                    help="输出两列：ORF_id + 聚合后的 iBAQ_norm（见 --agg）")

    # 可选：导出锚定集与缩放因子
    ap.add_argument("--out_anchors", default=None,
                    help="导出锚定 ORF 列表（单列 ORF_id）")
    ap.add_argument("--out_factors", default=None,
                    help="导出每个 Sample 的缩放因子 s_r（两列 Sample, scale_factor）")

    # 锚定与检出判定
    ap.add_argument("--present_by", choices=["ibaq_pos","available_true"],
                    default="ibaq_pos",
                    help="锚定判定：ibaq_pos= iBAQ_B>0（默认）；available_true= Available==TRUE")
    ap.add_argument("--min_anchor_prop", type=float, default=1.0,
                    help="锚定出现比例阈值（默认1.0=所有 run 都有；可设0.8放宽）")
    ap.add_argument("--anchor_only_canonical", action="store_true",
                    help="仅在 is_canonical==True 的 ORF 中选择锚定")

    # ★ 聚合方式（跨 run）
    ap.add_argument("--agg", choices=["geomean","median_log","mean"], default="geomean",
                    help="跨 run 聚合方式：geomean=几何均值(默认)；median_log=对数域中位数；mean=算术均值")

    args = ap.parse_args()

    # 读取列
    must = {"ORF_id","Sample","iBAQ_B"}
    optional = {"Available","is_canonical","Is_canonical"}
    header = pd.read_csv(args.in_tsv, sep="\t", nrows=0)
    missing = must - set(header.columns)
    if missing:
        raise SystemExit(f"缺少必要列：{', '.join(sorted(missing))}")
    usecols = list(must | (optional & set(header.columns)))

    try:
        df = pd.read_csv(args.in_tsv, sep="\t", usecols=usecols, engine="pyarrow")
    except Exception:
        df = pd.read_csv(args.in_tsv, sep="\t", usecols=usecols, low_memory=False)

    # 统一 canonical 列
    canon_colname = next((c for c in df.columns if c.lower()=="is_canonical"), None)
    if canon_colname is not None:
        df["__is_canon__"] = to_bool_series(df[canon_colname])
    else:
        df["__is_canon__"] = False

    # 数值化 iBAQ
    df["iBAQ_B"] = pd.to_numeric(df["iBAQ_B"], errors="coerce")

    # ---------------- 选择“共同锚定集” ----------------
    df_for_anchor = df
    if args.anchor_only_canonical:
        if canon_colname is None:
            raise SystemExit("启用 --anchor_only_canonical 但输入缺少 is_canonical/Is_canonical 列。")
        df_for_anchor = df[df["__is_canon__"]].copy()
        if df_for_anchor.empty:
            raise SystemExit("筛选 canonical 后为空，无法构建锚定集。")

    if args.present_by == "available_true" and "Available" in df_for_anchor.columns:
        present = to_bool_series(df_for_anchor["Available"])
    else:
        present = df_for_anchor["iBAQ_B"].fillna(0) > 0

    total_runs = df["Sample"].nunique()
    det_counts = (df_for_anchor.loc[present, ["ORF_id","Sample"]]
                    .drop_duplicates()
                    .groupby("ORF_id")["Sample"].nunique())
    anchors = det_counts[det_counts >= args.min_anchor_prop * total_runs].index.to_series()

    if anchors.empty:
        raise SystemExit("锚定集为空；请检查 --present_by / --min_anchor_prop / --anchor_only_canonical 设置。")

    if args.out_anchors:
        Path(args.out_anchors).parent.mkdir(parents=True, exist_ok=True)
        anchors.to_frame(name="ORF_id").to_csv(args.out_anchors, sep="\t", index=False)

    # ---------------- 计算 MoR 缩放因子 ----------------
    df_anchor = df[df["ORF_id"].isin(anchors)][["ORF_id","Sample","iBAQ_B"]].copy()
    df_anchor = df_anchor[df_anchor["iBAQ_B"].notna() & (df_anchor["iBAQ_B"] > 0)]
    if df_anchor.empty:
        raise SystemExit("锚定条目均为缺失或零，无法计算缩放因子。")

    # 参考谱：锚定 ORF 的几何平均（log-mean 再 exp）
    df_anchor["logI"] = np.log(df_anchor["iBAQ_B"])
    ref = df_anchor.groupby("ORF_id")["logI"].mean().pipe(np.exp)  # ref_i

    # 每个 run 的 MoR 因子：median_i( I_{i,r}/ref_i )
    df_anchor = df_anchor.merge(ref.rename("refI"), left_on="ORF_id", right_index=True, how="inner")
    df_anchor["ratio"] = df_anchor["iBAQ_B"] / df_anchor["refI"]
    scale = (df_anchor.groupby("Sample")["ratio"]
                     .median()
                     .rename("scale_factor")
                     .reset_index())

    if args.out_factors:
        Path(args.out_factors).parent.mkdir(parents=True, exist_ok=True)
        scale.to_csv(args.out_factors, sep="\t", index=False)

    # ---------------- 应用缩放 + 跨 run 聚合 ----------------
    df = df.merge(scale, on="Sample", how="left")
    df["scale_factor"] = df["scale_factor"].fillna(1.0)
    df["iBAQ_norm"] = df["iBAQ_B"] / df["scale_factor"]

    valid = df["iBAQ_norm"].notna() & (df["iBAQ_norm"] > 0)
    dfv = df.loc[valid, ["ORF_id","iBAQ_norm"]].copy()

    if dfv.empty:
        out = pd.DataFrame(columns=["ORF_id", "agg_iBAQ_norm"])
        Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out_tsv, sep="\t", index=False)
        print("[WARN] 归一化后无有效 iBAQ，输出空表。")
        return

    agg_name = {"geomean":"mean_relative_iBAQ_B",
                "median_log":"mean_relative_iBAQ_B",
                "mean":"mean_relative_iBAQ_B"}[args.agg]

    if args.agg == "geomean":
        out = (dfv.assign(logI=np.log(dfv["iBAQ_norm"]))
                  .groupby("ORF_id")["logI"].mean()
                  .pipe(np.exp)
                  .rename(agg_name)
                  .reset_index())
    elif args.agg == "median_log":
        out = (dfv.assign(logI=np.log(dfv["iBAQ_norm"]))
                  .groupby("ORF_id")["logI"].median()
                  .pipe(np.exp)
                  .rename(agg_name)
                  .reset_index())
    else:  # mean
        out = (dfv.groupby("ORF_id")["iBAQ_norm"]
                  .mean()
                  .rename(agg_name)
                  .reset_index())

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[OK] MoR 归一化 + 跨 run 聚合完成：{args.out_tsv}  shape={out.shape}")
    print(f"    锚定ORF数={len(anchors)}, 运行数={total_runs}, 仅canonical={args.anchor_only_canonical}, 聚合={args.agg}")

if __name__ == "__main__":
    try:
        pd.options.mode.copy_on_write = True
    except Exception:
        pass
    main()
