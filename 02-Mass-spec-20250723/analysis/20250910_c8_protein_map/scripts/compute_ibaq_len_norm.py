#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute iBAQ_B (length-normalized total intensity) per ORF and merge file3 fields.

Inputs:
  - file1: TSV with Peptide (already Peptide_IL) and Assigned_protein_id (all peptides are unique by assumption)
  - file2: TSV with Peptide_IL, Intensity (ignore Sample/Source; each peptide appears once globally)
  - file3: TSV with ORF_id, ORF_length, ORF_seq (optional), and any other metadata columns

Outputs:
  - TSV with columns: ORF_id, iBAQ_B, Available, TotalIntensity, [selected columns from file3]
"""

import argparse
import sys
import os
import pandas as pd
import numpy as np
from collections import defaultdict

DEFAULT_FILE3_COLS = [
    "ORF_type", "Start_codon", "ORF_length",
    "is_uniprot", "is_canonical",
    "All_peptide_n", "Unique_peptide_n"
]

def parse_args():
    ap = argparse.ArgumentParser(description="Compute iBAQ_B and merge file3 fields.")
    ap.add_argument("--file1", required=True, help="TSV with columns: Peptide (or Peptide_IL), Assigned_protein_id.")
    ap.add_argument("--file2", required=True, help="TSV with columns: Peptide_IL (or Peptide), Intensity.")
    ap.add_argument("--file3", required=True, help="TSV with columns: ORF_id, ORF_length (optional), ORF_seq (optional), and metadata.")
    ap.add_argument("-o", "--out", required=True, help="Output TSV.")
    ap.add_argument("--chunksize", type=int, default=500000, help="Chunk size for streaming large files.")
    ap.add_argument(
        "--file3-cols",
        default="default",
        help=("Which extra columns from file3 to append. "
              "'default' -> a sensible subset; "
              "'all' -> all columns from file3; "
              "or provide a comma-separated list, e.g. 'ORF_type,Start_codon,is_canonical'."))
    return ap.parse_args()

def load_peptide_to_orf(file1, chunksize):
    head = pd.read_csv(file1, sep="\t", nrows=5, dtype=str)
    cols = [c.strip() for c in head.columns]
    col_pep = "Peptide" if "Peptide" in cols else ("Peptide_IL" if "Peptide_IL" in cols else None)
    col_orf = "Assigned_protein_id" if "Assigned_protein_id" in cols else ("ORF_id" if "ORF_id" in cols else None)
    if col_pep is None or col_orf is None:
        raise ValueError(f"[file1] 必要列缺失，检测到列={cols}")
    pep2orf, all_orfs = {}, set()
    dup_warn = 0
    for chunk in pd.read_csv(file1, sep="\t", dtype=str, usecols=[col_pep, col_orf], chunksize=chunksize):
        chunk[col_pep] = chunk[col_pep].astype(str).str.strip()
        chunk[col_orf] = chunk[col_orf].astype(str).str.strip()
        for p, o in zip(chunk[col_pep].values, chunk[col_orf].values):
            if not p or p == "nan":
                continue
            if p in pep2orf and pep2orf[p] != o:
                dup_warn += 1
                continue
            pep2orf[p] = o
            all_orfs.add(o)
    if dup_warn > 0:
        print(f"[WARN] file1 中发现 {dup_warn} 个肽段映射多个 ORF，已保留首个映射。", file=sys.stderr)
    return pep2orf, all_orfs

def aggregate_intensity_by_orf(file2, pep2orf, chunksize):
    head = pd.read_csv(file2, sep="\t", nrows=5, dtype=str)
    cols = [c.strip() for c in head.columns]
    col_pep = "Peptide_IL" if "Peptide_IL" in cols else ("Peptide" if "Peptide" in cols else None)
    col_int = "Intensity" if "Intensity" in cols else None
    if col_pep is None or col_int is None:
        raise ValueError(f"[file2] 必要列缺失，检测到列={cols}")

    sum_int = defaultdict(float)
    detected_pep_sets = defaultdict(set)

    for chunk in pd.read_csv(file2, sep="\t", dtype=str, usecols=[col_pep, col_int], chunksize=chunksize):
        chunk[col_pep] = chunk[col_pep].astype(str).str.strip()
        chunk[col_int] = pd.to_numeric(chunk[col_int], errors="coerce")
        chunk = chunk[chunk[col_int] > 0]
        if chunk.empty:
            continue
        chunk = chunk[chunk[col_pep].isin(pep2orf.keys())]
        if chunk.empty:
            continue
        for p, inten in zip(chunk[col_pep].values, chunk[col_int].values):
            o = pep2orf.get(p)
            if o is None:
                continue
            sum_int[o] += float(inten)
            detected_pep_sets[o].add(p)

    df = pd.DataFrame({
        "ORF_id": list(sum_int.keys()),
        "TotalIntensity": [sum_int[o] for o in sum_int.keys()],
        "N_detected_unique_peptides": [len(detected_pep_sets[o]) for o in sum_int.keys()]
    })
    return df

def load_file3(file3_path):
    df3 = pd.read_csv(file3_path, sep="\t", dtype=str)
    if "ORF_id" not in df3.columns:
        raise ValueError(f"[file3] 必须包含 ORF_id 列，检测到列={list(df3.columns)}")

    # 处理 ORF_length（数值）与 ORF_seq 长度兜底
    if "ORF_length" in df3.columns:
        df3["ORF_length"] = pd.to_numeric(df3["ORF_length"], errors="coerce")
    else:
        df3["ORF_length"] = np.nan

    if "ORF_seq" in df3.columns:
        mask = df3["ORF_length"].isna() | (df3["ORF_length"] <= 0)
        if mask.any():
            seq_len = df3.loc[mask, "ORF_seq"].astype(str).apply(lambda s: len(s) if (isinstance(s, str) and s.lower() != "nan") else np.nan)
            df3.loc[mask, "ORF_length"] = pd.to_numeric(seq_len, errors="coerce")

    return df3

def select_file3_cols(df3, mode):
    """
    mode:
      - 'default': keep DEFAULT_FILE3_COLS if present (plus ORF_id)
      - 'all': keep all columns
      - comma list: keep specified columns if present
    """
    all_cols = list(df3.columns)
    if mode == "all":
        keep = all_cols
    elif mode == "default":
        keep = ["ORF_id"] + [c for c in DEFAULT_FILE3_COLS if c in all_cols and c != "ORF_id"]
    else:
        req = [c.strip() for c in mode.split(",") if c.strip()]
        keep = ["ORF_id"] + [c for c in req if c in all_cols and c != "ORF_id"]
        missing = [c for c in req if c not in all_cols]
        if missing:
            print(f"[WARN] file3 中缺少以下列，将被忽略: {missing}", file=sys.stderr)
    # 去重保持顺序
    seen, ordered = set(), []
    for c in keep:
        if c not in seen:
            seen.add(c); ordered.append(c)
    return df3[ordered].copy()

def main():
    args = parse_args()

    # 1) file1：肽段->ORF 映射；ORF全集
    pep2orf, all_orfs = load_peptide_to_orf(args.file1, args.chunksize)
    if not pep2orf:
        print("[ERROR] file1 未解析到映射。", file=sys.stderr); sys.exit(1)

    # 2) file2：按 ORF 聚合总强度
    df_sum = aggregate_intensity_by_orf(args.file2, pep2orf, args.chunksize)

    # 3) ORF 骨架（包含未检出的 ORF）
    df_orfs = pd.DataFrame({"ORF_id": sorted(all_orfs)})
    df = df_orfs.merge(df_sum, on="ORF_id", how="left")

    # 4) Available
    df["Available"] = df["N_detected_unique_peptides"].fillna(0).astype(int) > 0

    # 5) 读取 file3 & 选择要附加的列
    df3_full = load_file3(args.file3)
    df3_keep = select_file3_cols(df3_full, args.file3_cols)

    # 6) 并入 file3 字段（含 ORF_length）
    df = df.merge(df3_keep, on="ORF_id", how="left")

    # 7) 计算 iBAQ_B
    df["TotalIntensity"] = pd.to_numeric(df["TotalIntensity"], errors="coerce")
    df.loc[(df["TotalIntensity"].isna()) | (df["TotalIntensity"] <= 0), "TotalIntensity"] = np.nan

    if "ORF_length" in df.columns:
        df["ORF_length"] = pd.to_numeric(df["ORF_length"], errors="coerce")
        valid = (~df["TotalIntensity"].isna()) & (~df["ORF_length"].isna()) & (df["ORF_length"] > 0)
    else:
        valid = pd.Series(False, index=df.index)

    df["iBAQ_B"] = np.nan
    df.loc[valid, "iBAQ_B"] = df.loc[valid, "TotalIntensity"] / df.loc[valid, "ORF_length"]

    # 8) 组织输出列：固定头部 + 选中的 file3 列（去重保序）
    head_cols = ["ORF_id", "iBAQ_B", "Available", "TotalIntensity"]
    # 把 df3_keep 的列（除 ORF_id）接在后面
    tail_cols = [c for c in df3_keep.columns if c != "ORF_id"]
    out_cols = []
    for c in head_cols + tail_cols:
        if c in df.columns and c not in out_cols:
            out_cols.append(c)

    df_out = df[out_cols].copy()

    # 保存
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    df_out.to_csv(args.out, sep="\t", index=False)

    # 汇总
    n_all = df_out.shape[0]
    n_avail = int(df_out["Available"].sum())
    n_ibaq = int(df_out["iBAQ_B"].notna().sum())
    print(f"[DONE] ORFs={n_all}, Available=TRUE={n_avail}, iBAQ_B有值={n_ibaq}")
    print(f"[OUT] {args.out}")

if __name__ == "__main__":
    main()
