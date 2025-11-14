#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import argparse
from pathlib import Path

def main():
    ap = argparse.ArgumentParser(description="以 iBAQ 主表为主，左连接 RPF、Isoform(仅FL列)、RNA，可追加Gene_type，并合并unique肽数与相对iBAQ均值")
    ap.add_argument("--ibaq", required=True, help="../processed/quant/ibaq_b_with_total.tsv")
    ap.add_argument("--rpf",  required=True, help="/.../orf.rpf.psite.txt")
    ap.add_argument("--iso",  required=True, help="/.../isoform.expr.info.txt")
    ap.add_argument("--rna",  required=True, help="/.../rpkm_N_C_A.txt")
    ap.add_argument("--gene_anno", required=False, default=None,
                    help="/home/user/.../Ensembl_...txt（第4列Gene_type, 第5列Geneid）")

    # 新增：可选合并的两个结果表
    ap.add_argument("--uniq_pep", required=False, default=None,
                    help="../results/unique_peptide_counts.tsv（两列：ORF_id, Unique_peptide_n）")
    ap.add_argument("--rel_ibaq", required=False, default=None,
                    help="../results/ibaq_orf_mean_relative.tsv（两列：ORF_id, mean_relative_iBAQ）")

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
            iso_raw[c] = 0.0
    iso = iso_raw[need_iso_cols].copy()
    iso.rename(columns={"FL.BioSample_6":"FL", "FL_TPM.BioSample_6":"FL_TPM"}, inplace=True)

    # RNA 全列读取（与 iBAQ 通过 Geneid 连接）
    rna  = pd.read_csv(args.rna,  sep="\t", low_memory=False)

    # 以 ibaq 为主，逐步左连接；给潜在重名列加后缀
    merged = ibaq.merge(rpf, how="left", on="ORF_id", suffixes=("", "_rpf"))

    # iBAQ 中若有 Isoform_id，则用 Isoform_id ↔ isoform 连接
    left_iso_key = "Isoform_id" if "Isoform_id" in merged.columns else "isoform"
    if left_iso_key in merged.columns:
        merged = merged.merge(iso, how="left",
                              left_on=left_iso_key, right_on="isoform",
                              suffixes=("", "_iso"))
        merged.drop(columns=["isoform_iso"], errors="ignore", inplace=True)
        if "Isoform_id" in merged.columns and "isoform" in merged.columns:
            merged.drop(columns=["isoform"], inplace=True)

    merged = merged.merge(rna, how="left", on="Geneid", suffixes=("", "_rna"))

    # 可选追加：Gene_type（从注释文件第4/5列）
    if args.gene_anno:
        anno = pd.read_csv(args.gene_anno, delim_whitespace=True, header=None, low_memory=False, usecols=[3,4],
                           dtype={3:"string",4:"string"})
        anno.columns = ["Gene_type", "Geneid"]
        anno["Geneid"] = anno["Geneid"].str.strip()
        anno["Gene_type"] = anno["Gene_type"].str.strip()
        anno = anno.dropna(subset=["Geneid"]).drop_duplicates(subset=["Geneid"], keep="first")
        merged["Geneid"] = merged["Geneid"].astype("string").str.strip()
        merged = merged.merge(anno, how="left", on="Geneid")
        merged["Gene_type"] = merged["Gene_type"].astype("string").str.strip()
        mask_na = merged["Gene_type"].isna() | (merged["Gene_type"].str.strip() == "")
        mask_prefix = merged["ORF_id"].str.startswith(("sp", "tr"), na=False)
        merged.loc[mask_na & ~mask_prefix, "Gene_type"] = "novel"
        merged.loc[mask_na &  mask_prefix, "Gene_type"] = "protein_coding"
        # merged.loc[merged["Gene_type"].isna() | (merged["Gene_type"] == ""), "Gene_type"] = "novel"

    # 列名规范化：is_uniprot/is_canonical -> 首字母大写
    col_rename = {}
    if "is_uniprot" in merged.columns:
        col_rename["is_uniprot"] = "Is_uniprot"
    if "is_canonical" in merged.columns:
        col_rename["is_canonical"] = "Is_canonical"
    if col_rename:
        merged.rename(columns=col_rename, inplace=True)

    # ===== 新增：合并 Unique_peptide_n 与 mean_relative_iBAQ =====
    if args.uniq_pep:
        up = pd.read_csv(args.uniq_pep, sep="\t", usecols=["ORF_id","Unique_peptide_n","Unique_peptide_n_msfragger_closed"])
        merged = merged.merge(up, how="left", on="ORF_id")

    if args.rel_ibaq:
        ri = pd.read_csv(args.rel_ibaq, sep="\t")
        # 兼容列名大小写或不同命名
        if "mean_relative_iBAQ" not in ri.columns:
            # 如果用户改名为 mean_relative_ibaq 或其他，尽量探测
            cand = [c for c in ri.columns if c.lower().startswith("mean_relative")]
            if cand:
                ri = ri.rename(columns={cand[0]: "mean_relative_iBAQ"})
            else:
                raise SystemExit("ibaq_orf_mean_relative.tsv 缺少列：mean_relative_iBAQ")
        ri = ri[["ORF_id","mean_relative_iBAQ"]]
        # 保持为可空浮点，缺失就是 <NA>
        merged = merged.merge(ri, how="left", on="ORF_id")

    # 输出（以 NA 显示缺失）
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out, sep="\t", index=False, na_rep="NA")
    print(f"[OK] 合并完成：{args.out}  行数={len(merged)} 列数={merged.shape[1]}")

if __name__ == "__main__":
    main()
