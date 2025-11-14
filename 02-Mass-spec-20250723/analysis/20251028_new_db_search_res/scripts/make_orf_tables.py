#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import pandas as pd
import numpy as np
from pathlib import Path

def load_candidates(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    need = ["Peptide","Protein_id","orf_type","start_codon","isoform",
            "Geneid","is_uniprot","is_canonical"]
    for c in need:
        if c not in df.columns:
            raise SystemExit(f"[candidates.filtered.tsv] 缺少列：{c}")
    # 标准布尔化
    df["is_uniprot"]   = df["is_uniprot"].astype(str).str.lower().isin(["true","1","yes"])
    df["is_canonical"] = df["is_canonical"].astype(str).str.lower().isin(["true","1","yes"])
    return df

def drop_noncanon_sharing_with_canonical(df):
    # 若某肽存在 canonical 行，则删除其所有 non-canonical 行
    has_canon = df.groupby("Peptide")["is_canonical"].transform("any")
    mask_drop = has_canon & (~df["is_canonical"])
    return df.loc[~mask_drop].copy()

def parse_pb_fields(orfs: pd.Series):
    s = orfs.astype(str)
    is_pb = s.str.startswith(("PB.", "PB|"))
    part = s.where(is_pb, None)

    isoform = part.str.split(":", n=1, expand=True).iloc[:,0]
    after1  = part.str.split(":", n=1, expand=True).iloc[:,1]  # chr..:...
    chr_    = after1.str.split(":", n=1, expand=True).iloc[:,0]
    after2  = after1.str.split(":", n=1, expand=True).iloc[:,1]
    strand  = after2.str.split("|", n=1, expand=True).iloc[:,0]

    isoform = isoform.where(is_pb, "")
    chr_    = chr_.where(is_pb, "")
    strand  = strand.where(is_pb, "")
    return isoform, chr_, strand

def load_orf_seq_len(path):
    # 允许有/无表头；无论如何把第一列按空格截断
    df = pd.read_csv(path, sep="\t", dtype=str, header=0)
    if "ORF_id" not in df.columns:
        df = pd.read_csv(path, sep="\t", dtype=str, header=None,
                         names=["ORF_id","ORF_seq","ORF_length"])
    df["ORF_id"] = df["ORF_id"].astype(str).str.split(" ").str[0]  # 去掉描述
    # ORF_length 转为数值（如后续被填“NA”，写出时会自动转回字符串）
    df["ORF_length"] = pd.to_numeric(df["ORF_length"], errors="coerce")
    return df[["ORF_id","ORF_seq","ORF_length"]].drop_duplicates("ORF_id")

def load_isoform_info(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    for c in ["isoform","structural_category","subcategory"]:
        if c not in df.columns:
            raise SystemExit(f"[isoform.info] 缺少列：{c}")
    return df[["isoform","structural_category","subcategory"]].copy()

def _ensure_cols(df, cols, default="NA"):
    """若整列不存在则新增并填默认值；存在则不动"""
    for c in cols:
        if c not in df.columns:
            df[c] = default
    return df

def make_detail_table(cand_path, orf_seq_len_path, isoform_info_path, out_tsv):
    df = load_candidates(cand_path)
    df = drop_noncanon_sharing_with_canonical(df)

    # 重命名/解析
    df = df.rename(columns={"Protein_id":"ORF_id",
                            "orf_type":"ORF_type",
                            "start_codon":"Start_codon"})
    isoform_id_parsed, chr_, strand = parse_pb_fields(df["ORF_id"])

    # Isoform 优先用已有，其次解析
    df["Isoform_id"] = df["isoform"].where(df["isoform"].notna() & (df["isoform"]!=""), isoform_id_parsed)
    df["Chr"] = chr_
    df["Strand"] = strand

    # 合并序列/长度
    orf_seq = load_orf_seq_len(orf_seq_len_path)
    df = df.merge(orf_seq, how="left", on="ORF_id")

    # 合并 isoform 结构注释
    iso_info = load_isoform_info(isoform_info_path)
    df = df.merge(iso_info, how="left", left_on="Isoform_id", right_on="isoform") \
           .drop(columns=["isoform"], errors="ignore") \
           .rename(columns={
               "structural_category":"Isoform_structural_category",
               "subcategory":"Isoform_subcategory"
           })

    # ====== UniProt 定制：仅改 ORF_type / Start_codon，其它保持现值；空/缺则填 NA ======
    uni = df["is_uniprot"] == True
    df.loc[uni, "ORF_type"]    = "canonical"
    df.loc[uni, "Start_codon"] = "ATG"

    # 这些列若“该行为空值/缺失”，为 UniProt 行填充为 "NA"；已有值不改
    maybe_na_cols = [
        "Isoform_id","Chr","Strand","ORF_seq","ORF_length",
        "Geneid","Isoform_structural_category","Isoform_subcategory"
    ]
    df = _ensure_cols(df, maybe_na_cols, default="NA")
    for c in maybe_na_cols:
        # 仅对 UniProt 行，且该单元格为空/缺失时填 NA
        df.loc[uni & (df[c].isna() | (df[c].astype(str)=="")), c] = "NA"
    # ==========================================================================

    # 整理列与排序
    cols = ["Peptide","ORF_id","ORF_type","Start_codon","Isoform_id",
            "Chr","Strand","ORF_seq","ORF_length","Geneid",
            "Isoform_structural_category","Isoform_subcategory",
            "is_uniprot","is_canonical"]
    for c in cols:
        if c not in df.columns:
            df[c] = "" if c not in ["is_uniprot","is_canonical"] else False

    df = df.reindex(columns=cols)
    df = df.sort_values(by=["is_canonical","ORF_id","Peptide"], ascending=[False, True, True])

    df.to_csv(out_tsv, sep="\t", index=False)
    return df

def fold_to_orf_table(detail_df, out_counts_tsv):
    # 每个 peptide 在明细里命中多少 ORF？
    pep_hits = detail_df.groupby("Peptide")["ORF_id"].nunique()
    unique_pep = set(pep_hits[pep_hits==1].index)

    # 每 ORF 的静态注释（取第一条）
    annot_cols = ["ORF_type","Start_codon","Isoform_id","Chr","Strand",
                  "ORF_seq","ORF_length","Geneid",
                  "Isoform_structural_category","Isoform_subcategory",
                  "is_uniprot","is_canonical"]
    annot = (detail_df
             .sort_values(by=["is_canonical","ORF_id"], ascending=[False,True])
             .drop_duplicates(subset=["ORF_id"])[["ORF_id"]+annot_cols])

    # 统计 All_peptide_n / Unique_peptide_n
    grp = detail_df.groupby("ORF_id")
    counts = grp.agg(All_peptide_n=("Peptide", lambda x: x.nunique()),
                     Unique_peptide_n=("Peptide", lambda x: sum(p in unique_pep for p in set(x))))
    counts = counts.reset_index()

    out = annot.merge(counts, on="ORF_id", how="left")
    out[["All_peptide_n","Unique_peptide_n"]] = out[["All_peptide_n","Unique_peptide_n"]].fillna(0).astype(int)

    # 排序（canonical 优先）
    out = out.sort_values(by=["is_canonical","ORF_id"], ascending=[False, True])
    out.to_csv(out_counts_tsv, sep="\t", index=False)
    return out

def main():
    ap = argparse.ArgumentParser(description="从 candidates.filtered.tsv 生成明细与折叠表（去除与canonical共享的非canonical行）")
    ap.add_argument("--cand_filtered", required=True, help="candidates.filtered.tsv")
    ap.add_argument("--orf_seq_len",   required=True, help="由 FASTA 制作的 orf_seq_len.tsv（三列：ORF_id ORF_seq ORF_length）")
    ap.add_argument("--isoform_info",  required=True, help="isoform.expr.info.txt（含 isoform/structural_category/subcategory）")
    ap.add_argument("--out_detail",    required=True, help="输出：明细表 TSV")
    ap.add_argument("--out_orf",       required=True, help="输出：ORF折叠表 TSV（含 All_peptide_n / Unique_peptide_n）")
    args = ap.parse_args()

    Path(args.out_detail).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_orf).parent.mkdir(parents=True, exist_ok=True)

    detail = make_detail_table(args.cand_filtered, args.orf_seq_len, args.isoform_info, args.out_detail)
    fold_to_orf_table(detail, args.out_orf)

    print("Done.")
    print(f"- {args.out_detail}")
    print(f"- {args.out_orf}")

if __name__ == "__main__":
    main()
