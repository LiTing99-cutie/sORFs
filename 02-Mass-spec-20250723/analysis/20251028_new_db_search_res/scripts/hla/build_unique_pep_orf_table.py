#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HLA 特异性肽段整合（行数一致版）：
1) 基于 HLA 结果筛出 “Mapped Proteins” 为空的行；
2) 以 Protein 作为 ORF_id，按行计数 Unique_peptide_n；
3) 基于 ORF_id→Unique_peptide_n 的映射，左连接 orf_seq_len / RPF / ISO / RNA / GENE_ANNO；
4) 输出两个文件：--out-map 与 --out-extended，行数严格一致。
"""
import argparse
from pathlib import Path
import re
import sys
import pandas as pd
import numpy as np

pd.options.mode.chained_assignment = None

def read_tsv(path, **kwargs):
    return pd.read_csv(path, sep="\t", dtype=str, quoting=3, **kwargs)

def empty_like_nan(x: str) -> bool:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return True
    s = str(x).strip()
    return s == "" or s in {"-", "—", "NA", "na", "None", "null"}

def parse_orf_id(orf_id: str):
    """
    示例：PB.41013.1:chr6:+|1|330:5:44|noncoding|GTG
    返回：Isoform_id, Chr, Strand, ORF_type, Start_codon, Is_uniprot, Is_canonical
    规则：
      - sp| 或 tr| 开头：Is_uniprot=True, Is_canonical=True（按需设定）
      - 其它：从 ORF_id 中解析 ORF_type，ORF_type=='canonical' 则 Is_canonical=True
    """
    out = {
        "Isoform_id": np.nan, "Chr": np.nan, "Strand": np.nan,
        "ORF_type": np.nan, "Start_codon": np.nan,
        "Is_uniprot": False, "Is_canonical": False
    }
    if pd.isna(orf_id):
        return out
    s = str(orf_id).strip()

    if s.startswith(("sp|", "tr|")):
        out["Is_uniprot"] = True
        out["Is_canonical"] = True   # 按你的新要求
        return out

    m = re.match(
        r'^(?P<Isoform_id>[^:]+):(?P<Chr>[^:]+):(?P<Strand>[+-])\|.*\|(?P<ORF_type>[^|]+)\|(?P<Start_codon>[A-Za-z]+)$',
        s
    )
    if m:
        d = m.groupdict()
        d["Start_codon"] = d["Start_codon"].upper()
        d["Is_uniprot"] = False
        d["Is_canonical"] = (d["ORF_type"] == "canonical")
        out.update(d)
        return out

    parts = s.split(":")
    if parts:
        out["Isoform_id"] = parts[0]
    return out

def detect_iso_fl_cols(df_iso):
    fl_col = next((c for c in df_iso.columns if re.match(r'^FL(\.|$)', c)), None)
    fltpm_col = next((c for c in df_iso.columns if re.match(r'^FL_TPM(\.|$)', c)), None)
    return fl_col, fltpm_col

def safe_drop_dups(df, key):
    if key in df.columns:
        return df.drop_duplicates(subset=[key], keep="first")
    return df

def to_numeric(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def main():
    ap = argparse.ArgumentParser(description="HLA 特异性肽段整合并输出双表（行数一致）")
    ap.add_argument("--hla", required=True, help="HLA 搜库结果（TSV，需含：Protein, Peptide, Mapped Proteins）")
    ap.add_argument("--rna", required=True, help="RNA：Geneid, N, C, A")
    ap.add_argument("--iso", required=True, help="ISO：isoform, structural_category, subcategory, associated_gene, FL*, FL_TPM*")
    ap.add_argument("--rpf", required=True, help="RPF：含 ORF_id 及 RPF 指标列")
    ap.add_argument("--orf-seq-len", required=True, help="ORF 序列与长度：ORF_id, ORF_seq, ORF_length")
    ap.add_argument("--gene-anno", required=True, help="基因注释（空格分隔；第4列 Gene_type，第5列 Geneid）")
    ap.add_argument("--out-map", required=True, help="输出：ORF_id 与 Unique_peptide_n 映射表（TSV）")
    ap.add_argument("--out-extended", required=True, help="输出：扩展汇总表（TSV）")
    args = ap.parse_args()

    # 1) HLA → 过滤特异性肽段 → 按 ORF(Protein) 计“行数”作为 Unique_peptide_n
    hla = read_tsv(args.hla)
    hla.columns = [c.strip() for c in hla.columns]
    needed = {"Protein", "Peptide", "Mapped Proteins"}
    if not needed.issubset(hla.columns):
        raise SystemExit(f"HLA 文件需包含列：{needed}")

    mask_unique = hla["Mapped Proteins"].map(empty_like_nan)
    hla_u = hla.loc[mask_unique, ["Protein", "Peptide"]].copy()
    hla_u["Protein"] = hla_u["Protein"].astype(str).str.strip()
    hla_u["Peptide"] = hla_u["Peptide"].astype(str).str.strip()
    hla_u = hla_u[(hla_u["Protein"] != "") & (hla_u["Peptide"] != "")]

    counts = (hla_u.groupby("Protein", as_index=False)
                    .size()
                    .rename(columns={"Protein": "ORF_id", "size": "Unique_peptide_n"}))
    counts["Unique_peptide_n"] = counts["Unique_peptide_n"].astype("Int64")

    # 保存映射表（作为后续扩展表的严格基准）
    Path(args.out_map).parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(args.out_map, sep="\t", index=False)

    # 2) 解析 ORF_id（严格以 counts 为基准，防止行数变化）
    parsed = pd.DataFrame([{"ORF_id": oid, **parse_orf_id(oid)} for oid in counts["ORF_id"]])
    # 保证 ORF_id 唯一
    parsed = safe_drop_dups(parsed, "ORF_id")

    # 3) 合并 orf_seq_len
    orf_seq = read_tsv(args.orf_seq_len)
    orf_seq.columns = [c.strip() for c in orf_seq.columns]
    need_cols_seq = {"ORF_id", "ORF_seq", "ORF_length"}
    if not need_cols_seq.issubset(set(orf_seq.columns)):
        raise SystemExit(f"orf_seq_len 文件应包含列：{need_cols_seq}")
    orf_seq = safe_drop_dups(orf_seq, "ORF_id")

    df = counts.merge(parsed, on="ORF_id", how="left") \
               .merge(orf_seq[list(need_cols_seq)], on="ORF_id", how="left")

    # 4) RPF（去除 ORF_id 重复，优先保留 RPF_reads 大的行）
    rpf = read_tsv(args.rpf)
    rpf.columns = [c.strip() for c in rpf.columns]
    rpf_cols = ["ORF_id", "RPF_reads", "Psites_number", "RPF_RPKM", "Psites_RPKM",
                "RPF_codon_coverage", "Psites_codon_coverage"]
    rpf = rpf[[c for c in rpf_cols if c in rpf.columns]].copy()
    if "ORF_id" in rpf.columns:
        if "RPF_reads" in rpf.columns:
            rpf["_order"] = pd.to_numeric(rpf["RPF_reads"], errors="coerce").fillna(-1)
            rpf = rpf.sort_values("_order", ascending=False).drop(columns="_order")
        rpf = rpf.drop_duplicates(subset=["ORF_id"], keep="first")
    num_cols = [c for c in rpf.columns if c != "ORF_id"]
    rpf = to_numeric(rpf, num_cols)
    df = df.merge(rpf, on="ORF_id", how="left")

    # 5) ISO（通过 Isoform_id；去重以免放大行数）
    iso = read_tsv(args.iso)
    iso.columns = [c.strip() for c in iso.columns]
    need_iso = {"isoform", "structural_category", "subcategory", "associated_gene"}
    if not need_iso.issubset(set(iso.columns)):
        raise SystemExit(f"ISO 文件应至少包含列：{need_iso} 以及 FL*, FL_TPM*")
    fl_col, fltpm_col = detect_iso_fl_cols(iso)
    iso["FL"] = pd.to_numeric(iso[fl_col], errors="coerce") if fl_col else np.nan
    iso["FL_TPM"] = pd.to_numeric(iso[fltpm_col], errors="coerce") if fltpm_col else np.nan
    iso_sub = iso[["isoform", "structural_category", "subcategory", "associated_gene", "FL", "FL_TPM"]].copy()
    iso_sub = iso_sub.rename(columns={
        "isoform": "Isoform_id",
        "structural_category": "Isoform_structural_category",
        "subcategory": "Isoform_subcategory",
        "associated_gene": "Geneid"
    })
    iso_sub = safe_drop_dups(iso_sub, "Isoform_id")
    df = df.merge(iso_sub, on="Isoform_id", how="left")

    # 6) RNA（通过 Geneid；去重）
    rna = read_tsv(args.rna)
    rna.columns = [c.strip() for c in rna.columns]
    if "Geneid" not in rna.columns:
        raise SystemExit("RNA 文件需包含列：Geneid, N, C, A")
    for c in ["N", "C", "A"]:
        if c in rna.columns:
            rna[c] = pd.to_numeric(rna[c], errors="coerce")
        else:
            rna[c] = np.nan
    rna = safe_drop_dups(rna, "Geneid")
    df = df.merge(rna[["Geneid", "N", "C", "A"]], on="Geneid", how="left")

    # 7) GENE_ANNO（空格分隔，第4列 Gene_type，第5列 Geneid；去重）
    gene_anno = pd.read_csv(args.gene_anno, sep=' ', header=None, dtype=str)
    if gene_anno.shape[1] < 5:
        raise SystemExit("GENE_ANNO 至少应有 5 列（第4列 Gene_type，第5列 Geneid）")
    gene_anno = gene_anno.iloc[:, :5].copy()
    gene_anno.columns = ["col1", "col2", "col3", "Gene_type", "Geneid"]
    gene_anno = safe_drop_dups(gene_anno, "Geneid")
    df = df.merge(gene_anno[["Geneid", "Gene_type"]], on="Geneid", how="left")

    # 列顺序
    ordered_cols = [
        "ORF_id", "Unique_peptide_n",
        "Isoform_id", "Chr", "Strand",
        "ORF_type", "Start_codon", "Is_uniprot", "Is_canonical",
        "ORF_length", "ORF_seq",
        "Isoform_structural_category", "Isoform_subcategory",
        "Geneid", "Gene_type",
        "RPF_reads", "Psites_number", "RPF_RPKM", "Psites_RPKM",
        "RPF_codon_coverage", "Psites_codon_coverage",
        "FL", "FL_TPM", "N", "C", "A"
    ]
    final_cols = [c for c in ordered_cols if c in df.columns]
    # 为确保行数与映射表一致，严格以 counts 的 ORF_id 顺序重排
    df_out = df[final_cols]
    # 行数一致性检查
    if len(df_out) != len(counts):
        raise SystemExit(f"[ERR] 扩展表行数({len(df_out)}) 与映射表({len(counts)})不一致，请检查上游表是否存在重复键。")

    Path(args.out_extended).parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(args.out_extended, sep="\t", index=False)

    print(f"[OK] 映射表: {args.out_map}  （{len(counts)} 条 ORF）")
    print(f"[OK] 扩展表: {args.out_extended}  （{len(df_out)} 条 ORF；与映射表一致）")

if __name__ == "__main__":
    main()
