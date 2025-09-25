#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import argparse
from pathlib import Path

def main():
    ap = argparse.ArgumentParser(description="以 iBAQ 主表为主，左连接 RPF、Isoform(仅FL列)、RNA，并可追加Gene_type")
    ap.add_argument("--ibaq", required=True, help="../processed/quant/ibaq_b_with_total.tsv")
    ap.add_argument("--rpf",  required=True, help="/.../orf.rpf.psite.txt")
    ap.add_argument("--iso",  required=True, help="/.../isoform.expr.info.txt")
    ap.add_argument("--rna",  required=True, help="/.../rpkm_N_C_A.txt")
    ap.add_argument("--gene_anno", required=False, default=None,
                    help="/home/user/.../Ensembl_...txt（第4列Gene_type, 第5列Geneid）")
    ap.add_argument("--out",  required=True, help="输出合并结果 .tsv")
    args = ap.parse_args()

    # 读入主表
    ibaq = pd.read_csv(args.ibaq, sep="\t", low_memory=False)

    # RPF 全列读取（与 iBAQ 通过 ORF_id 连接）
    rpf  = pd.read_csv(args.rpf,  sep="\t", low_memory=False)

    # iso 仅取 isoform + 需要的两列（若不存在则填0）
    iso_raw = pd.read_csv(args.iso, sep="\t", low_memory=False)
    need_iso_cols = ["isoform", "FL.BioSample_6", "FL_TPM.BioSample_6"]
    for c in need_iso_cols:
        if c not in iso_raw.columns:
            if c == "isoform":
                raise SystemExit("isoform.expr.info.txt 缺少列: isoform")
            # 其余两列缺失则补0
            iso_raw[c] = 0.0
    iso = iso_raw[need_iso_cols].copy()
    iso.rename(columns={"FL.BioSample_6":"FL", "FL_TPM.BioSample_6":"FL_TPM"}, inplace=True)

    # RNA 全列读取（与 iBAQ 通过 Geneid 连接）
    rna  = pd.read_csv(args.rna,  sep="\t", low_memory=False)

    # 以 ibaq 为主，逐步左连接；给潜在重名列加后缀
    merged = ibaq.merge(rpf, how="left", on="ORF_id", suffixes=("", "_rpf"))

    # iBAQ 中若有 Isoform_id，则用 Isoform_id ↔ isoform 连接
    left_iso_key = "Isoform_id" if "Isoform_id" in merged.columns else "isoform"
    if left_iso_key not in merged.columns:
        # 若 iBAQ 没有 Isoform_id/isoform，则不做 iso 合并
        pass
    else:
        # 合并 iso（右表只带 FL/FL_TPM 两列；键为 isoform）
        merged = merged.merge(iso, how="left",
                            left_on=left_iso_key, right_on="isoform",
                            suffixes=("", "_iso"))

        # 清理冗余键列
        # 情况1：pandas 会生成右侧键列 isoform_iso（如果左侧也叫 isoform）
        merged.drop(columns=["isoform_iso"], errors="ignore", inplace=True)
        # 情况2：仍遗留右表的 isoform，而主表已有 Isoform_id，则删掉右表 isoform
        if "Isoform_id" in merged.columns and "isoform" in merged.columns:
            merged.drop(columns=["isoform"], inplace=True)

    merged = merged.merge(rna, how="left", on="Geneid", suffixes=("", "_rna"))

    # 可选追加：Gene_type（从提供的注释文件第4/5列）
    if args.gene_anno:
        # 文件列：第4列=Gene_type, 第5列=Geneid（人类自然数从1计；pandas用0基索引→usecols=[3,4]）
        anno = pd.read_csv(args.gene_anno, sep=" ", header=None, low_memory=False, usecols=[3,4])
        anno.columns = ["Gene_type", "Geneid"]
        # 去重以避免多重匹配
        anno = anno.drop_duplicates(subset=["Geneid"], keep="first")
        merged = merged.merge(anno, how="left", on="Geneid")
        merged["Gene_type"] = merged["Gene_type"].astype("string")
        merged["Gene_type"] = merged["Gene_type"].str.strip()
        merged.loc[merged["Gene_type"].isna() | (merged["Gene_type"] == ""), "Gene_type"] = "novel"

    # 列名规范化：is_uniprot/is_canonical -> 首字母大写
    col_rename = {}
    if "is_uniprot" in merged.columns:
        col_rename["is_uniprot"] = "Is_uniprot"
    if "is_canonical" in merged.columns:
        col_rename["is_canonical"] = "Is_canonical"
    if col_rename:
        merged.rename(columns=col_rename, inplace=True)

    # 输出
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] 合并完成：{args.out}  行数={len(merged)} 列数={merged.shape[1]}")

if __name__ == "__main__":
    main()
