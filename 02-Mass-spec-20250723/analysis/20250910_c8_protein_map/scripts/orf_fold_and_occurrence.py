#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from pathlib import Path
import pandas as pd

def main():
    ap = argparse.ArgumentParser(description="按 ORF 折叠：先丢列→聚合→计 Run_occurrence（出现即计）")
    ap.add_argument("--in", dest="in_tsv", required=True, help="all_samples.ibaq_b_with_total.tsv")
    ap.add_argument("--out", dest="out_tsv", required=True, help="输出折叠后的 TSV")
    args = ap.parse_args()

    # 只读需要列（保守：读全表头，drop 在内存中完成）
    try:
        df = pd.read_csv(args.in_tsv, sep="\t", engine="pyarrow")
    except Exception:
        df = pd.read_csv(args.in_tsv, sep="\t", low_memory=False)

    # 必要列检查
    for c in ("ORF_id", "Sample"):
        if c not in df.columns:
            raise SystemExit(f"缺少必要列：{c}")

    # 1) 先 drop 指定列
    drop_cols = {"iBAQ_A","iBAQ_B","Available","TotalIntensity","All_peptide_n","Unique_peptide_n","Theo_peptide_n"}
    drop_exist = [c for c in drop_cols if c in df.columns]
    df = df.drop(columns=drop_exist, errors="ignore")

    # 2) Run_occurrence：ORF 出现过的不同样本数（与 iBAQ 无关）
    run_occ = (df[["ORF_id","Sample"]]
               .dropna()
               .drop_duplicates()
               .groupby("ORF_id")["Sample"]
               .nunique()
               .rename("Run_occurrence")
               .reset_index())

    # 3) 折叠：为每个 ORF 选代表行（稳定且高效：按原顺序保留首行）
    #    若你想优先选择某条件（如 is_uniprot==False），可在这里先排序再 drop_duplicates
    df_sorted = df.copy()
    # 不保留 Sample 列到折叠表
    df_sorted = df_sorted.drop(columns=["Sample"], errors="ignore")
    folded = df_sorted.drop_duplicates(subset=["ORF_id"], keep="first")

    # 4) 合并 Run_occurrence
    out = folded.merge(run_occ, on="ORF_id", how="left")
    out["Run_occurrence"] = out["Run_occurrence"].fillna(0).astype("Int64")

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[OK] ORF 折叠完成：{args.out_tsv}  shape={out.shape}")

if __name__ == "__main__":
    # 小优化
    try: pd.options.mode.copy_on_write = True
    except: pass
    main()
