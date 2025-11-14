#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ORF 扩展信息整合（仅需 ORF_id 列）：
- 输入：包含 ORF_id 的文件（tsv/csv/一列文本；可指定列名）
- 过程：解析 ORF_id → 合并 orf_seq_len / RPF / ISO / RNA / GENE_ANNO
- 输出：--out-extended 扩展汇总表（行数与输入唯一 ORF_id 数一致，且顺序与输入一致）

备注：
- 识别 UniProt（sp|/tr|）ID：Is_uniprot=True, Is_canonical=True
- ISO 表自动探测 FL / FL_TPM 列名（以 "FL" 与 "FL_TPM" 开头）
"""

import argparse
from pathlib import Path
import re
import sys
import pandas as pd
import numpy as np

pd.options.mode.chained_assignment = None

def read_table_any(path, sep="\t", dtype=str):
    # 自动按扩展名选择分隔符
    if str(path).endswith(".csv"):
        return pd.read_csv(path, dtype=dtype, quoting=3)
    return pd.read_csv(path, sep=sep, dtype=dtype, quoting=3)

def read_orf_ids(path, colname="ORF_id", sep="\t"):
    """
    读取 ORF 列：
    - 若文件有表头且包含 colname，则取该列；
    - 否则尝试按单列无表头读取；
    - 保留输入顺序，去除空值。
    """
    try:
        df = read_table_any(path, sep=sep)
        cols = [c.strip() for c in df.columns]
        df.columns = cols
        if colname in df.columns:
            s = df[colname].astype(str)
        elif len(df.columns) == 1:
            s = df.iloc[:, 0].astype(str)
        else:
            raise SystemExit(f"[ERR] 输入文件不包含 '{colname}' 列，且列数>1 无法自动识别。请用 --orf-col 指定。现有列：{cols}")
    except pd.errors.EmptyDataError:
        raise SystemExit("[ERR] 输入 ORF 列文件为空。")

    s = s.map(lambda x: x.strip()).replace({"": np.nan, "NA": np.nan, "na": np.nan, "None": np.nan, "null": np.nan})
    s = s.dropna()
    # 保序去重
    orf_ids = pd.unique(s)
    return pd.DataFrame({"ORF_id": orf_ids})

def parse_orf_id(orf_id: str):
    """
    示例：PB.41013.1:chr6:+|1|330:5:44|noncoding|GTG
    返回：Isoform_id, Chr, Strand, ORF_type, Start_codon, Is_uniprot, Is_canonical
    规则：
      - sp| 或 tr| 开头：Is_uniprot=True, Is_canonical=True
      - 其它：从 ORF_id 中解析 ORF_type；ORF_type=='canonical' 则 Is_canonical=True
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
        out["Is_canonical"] = True
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
    ap = argparse.ArgumentParser(description="ORF 扩展信息整合（仅需 ORF_id 列）")
    ap.add_argument("--orf-list", required=True, help="包含 ORF_id 的文件（tsv/csv/一列文本）")
    ap.add_argument("--orf-col", default="ORF_id", help="ORF_id 列名（输入表有表头且不为 ORF_id 时指定）")
    ap.add_argument("--orf-seq-len", required=True, help="ORF 序列与长度：ORF_id, ORF_seq, ORF_length")
    ap.add_argument("--rpf", required=True, help="RPF：含 ORF_id 及 RPF 指标列")
    ap.add_argument("--iso", required=True, help="ISO：isoform, structural_category, subcategory, associated_gene, FL*, FL_TPM*")
    ap.add_argument("--rna", required=True, help="RNA：Geneid, N, C, A")
    ap.add_argument("--gene-anno", required=True, help="基因注释（空格分隔；第4列 Gene_type，第5列 Geneid）")
    ap.add_argument("--out-extended", required=True, help="输出扩展汇总表（TSV）")
    args = ap.parse_args()

    # 0) 读取 ORF 列（保序去重）
    base = read_orf_ids(args.orf_list, colname=args.orf_col)
    if base.empty:
        raise SystemExit("[ERR] 未读取到任何 ORF_id。")
    # 解析 ORF_id
    parsed = pd.DataFrame([{"ORF_id": oid, **parse_orf_id(oid)} for oid in base["ORF_id"]])
    # 1) 合并 orf_seq_len
    orf_seq = read_table_any(args.orf_seq_len)
    orf_seq.columns = [c.strip() for c in orf_seq.columns]
    need_cols_seq = {"ORF_id", "ORF_seq", "ORF_length"}
    if not need_cols_seq.issubset(set(orf_seq.columns)):
        raise SystemExit(f"[ERR] orf_seq_len 需包含列：{need_cols_seq}")
    orf_seq = safe_drop_dups(orf_seq, "ORF_id")

    df = base.merge(parsed, on="ORF_id", how="left") \
             .merge(orf_seq[list(need_cols_seq)], on="ORF_id", how="left")

    # 2) 合并 RPF（同一 ORF_id 取 RPF_reads 最大的记录）
    rpf = read_table_any(args.rpf)
    rpf.columns = [c.strip() for c in rpf.columns]
    rpf_cols = ["ORF_id", "RPF_reads", "Psites_number", "RPF_RPKM", "Psites_RPKM",
                "RPF_codon_coverage", "Psites_codon_coverage"]
    avail_rpf_cols = [c for c in rpf_cols if c in rpf.columns]
    if "ORF_id" not in rpf.columns:
        raise SystemExit("[ERR] RPF 表需包含 ORF_id 列。")
    rpf = rpf[avail_rpf_cols].copy()
    if "RPF_reads" in rpf.columns:
        rpf["_order"] = pd.to_numeric(rpf["RPF_reads"], errors="coerce").fillna(-1)
        rpf = rpf.sort_values("_order", ascending=False).drop(columns="_order")
    rpf = rpf.drop_duplicates(subset=["ORF_id"], keep="first")
    rpf = to_numeric(rpf, [c for c in rpf.columns if c != "ORF_id"])
    df = df.merge(rpf, on="ORF_id", how="left")

    # 3) 合并 ISO（通过 Isoform_id）
    iso = read_table_any(args.iso)
    iso.columns = [c.strip() for c in iso.columns]
    need_iso = {"isoform", "structural_category", "subcategory", "associated_gene"}
    if not need_iso.issubset(set(iso.columns)):
        raise SystemExit(f"[ERR] ISO 至少包含列：{need_iso}，并建议含 FL*/FL_TPM*")
    fl_col, fltpm_col = detect_iso_fl_cols(iso)
    iso["FL"] = pd.to_numeric(iso[fl_col], errors="coerce") if fl_col and fl_col in iso.columns else np.nan
    iso["FL_TPM"] = pd.to_numeric(iso[fltpm_col], errors="coerce") if fltpm_col and fltpm_col in iso.columns else np.nan
    iso_sub = iso[["isoform", "structural_category", "subcategory", "associated_gene", "FL", "FL_TPM"]].copy()
    iso_sub = iso_sub.rename(columns={
        "isoform": "Isoform_id",
        "structural_category": "Isoform_structural_category",
        "subcategory": "Isoform_subcategory",
        "associated_gene": "Geneid"
    })
    iso_sub = safe_drop_dups(iso_sub, "Isoform_id")
    df = df.merge(iso_sub, on="Isoform_id", how="left")

    # 4) 合并 RNA（通过 Geneid）
    rna = read_table_any(args.rna)
    rna.columns = [c.strip() for c in rna.columns]
    if "Geneid" not in rna.columns:
        raise SystemExit("[ERR] RNA 表需包含列：Geneid, N, C, A")
    for c in ["N", "C", "A"]:
        rna[c] = pd.to_numeric(rna[c], errors="coerce") if c in rna.columns else np.nan
    rna = safe_drop_dups(rna, "Geneid")
    df = df.merge(rna[["Geneid", "N", "C", "A"]], on="Geneid", how="left")

    # 5) 合并基因注释（空格分隔，第4列 Gene_type，第5列 Geneid）
    gene_anno = pd.read_csv(args.gene_anno, sep=' ', header=None, dtype=str)
    if gene_anno.shape[1] < 5:
        raise SystemExit("[ERR] GENE_ANNO 至少 5 列（第4列 Gene_type，第5列 Geneid）")
    gene_anno = gene_anno.iloc[:, :5].copy()
    gene_anno.columns = ["col1", "col2", "col3", "Gene_type", "Geneid"]
    gene_anno = safe_drop_dups(gene_anno, "Geneid")
    df = df.merge(gene_anno[["Geneid", "Gene_type"]], on="Geneid", how="left")

    # 输出列顺序
    ordered_cols = [
        "ORF_id",
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
    df_out = df[final_cols].copy()

    # 行数一致性校验
    if len(df_out) != len(base):
        raise SystemExit(f"[ERR] 输出行数({len(df_out)}) ≠ 输入唯一 ORF 数({len(base)})。请检查是否存在重复键导致倍增。")

    Path(args.out_extended).parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(args.out_extended, sep="\t", index=False)
    print(f"[OK] 扩展表: {args.out_extended}  （{len(df_out)} 条 ORF；与输入唯一 ORF 数一致）")

if __name__ == "__main__":
    main()
