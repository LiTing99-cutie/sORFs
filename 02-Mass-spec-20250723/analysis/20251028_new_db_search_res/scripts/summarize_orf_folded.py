#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import pandas as pd
from pathlib import Path

def to_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true","1","yes"])

def stats_series(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce").dropna()
    if x.empty:
        return pd.Series({"count":0,"min":None,"q25":None,"median":None,"mean":None,"q75":None,"max":None})
    return pd.Series({
        "count": int(x.count()),
        "min": float(x.min()),
        "q25": float(x.quantile(0.25)),
        "median": float(x.median()),
        "mean": float(x.mean()),
        "q75": float(x.quantile(0.75)),
        "max": float(x.max()),
    })

def main():
    ap = argparse.ArgumentParser(
        description="在 ORF 折叠表上做统计（默认仅统计 Unique_peptide_n>0 的行）"
    )
    ap.add_argument("-i","--input", default="../processed/pep_assign/orf_folded_counts.tsv",
                    help="输入折叠表（含 ORF_id/ORF_type/Start_codon/ORF_length/is_uniprot/is_canonical/Unique_peptide_n 等），默认 ../processed/pep_assign/orf_folded_counts.tsv")
    ap.add_argument("-o","--outdir", default="../processed/pep_assign/summary",
                    help="输出目录，默认 ../processed/pep_assign/summary")
    ap.add_argument("--min-unique", type=int, default=1,
                    help="筛选阈值：仅统计 Unique_peptide_n >= 此值 的行，默认 1（即>0）")
    ap.add_argument("--sep", default="\t", help="输入分隔符，默认 TAB")
    args = ap.parse_args()

    inp = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(inp, sep=args.sep, dtype=str)

    # 数值化 & 过滤
    df["Unique_peptide_n"] = pd.to_numeric(df.get("Unique_peptide_n"), errors="coerce").fillna(0).astype(int)
    f = df[df["Unique_peptide_n"] >= args.min_unique].copy()

    # 规范布尔列
    f["is_uniprot"] = to_bool(f.get("is_uniprot", pd.Series(False, index=f.index)))
    f["is_canonical"] = to_bool(f.get("is_canonical", pd.Series(False, index=f.index)))

    # 1) ORF_type 计数
    orf_type_counts = f.get("ORF_type", pd.Series(["NA"]*len(f))).fillna("NA").replace({"": "NA"}).value_counts(dropna=False)
    orf_type_counts.to_csv(outdir/"counts_by_ORF_type.tsv", sep="\t", header=["count"])

    # 2) Start_codon 计数
    start_counts = f.get("Start_codon", pd.Series(["NA"]*len(f))).fillna("NA").replace({"": "NA"}).value_counts(dropna=False)
    start_counts.to_csv(outdir/"counts_by_Start_codon.tsv", sep="\t", header=["count"])

    # 3) UniProt 数量
    uniprot_n = int(f["is_uniprot"].sum())

    # 4) Canonical 数量
    canonical_n = int(f["is_canonical"].sum())

    # 5) 按 canonical 分组的 ORF_length 统计
    if "ORF_length" not in f.columns:
        f["ORF_length"] = None
    stats_len = f.groupby("is_canonical")["ORF_length"].apply(stats_series).unstack()
    stats_len.index = stats_len.index.map({True: "canonical", False: "non_canonical"})
    stats_len.to_csv(outdir/"stats_ORF_length_by_canonical.tsv", sep="\t")

    # 6) 按 canonical 分组的 Unique_peptide_n 统计
    stats_up = f.groupby("is_canonical")["Unique_peptide_n"].apply(stats_series).unstack()
    stats_up.index = stats_up.index.map({True: "canonical", False: "non_canonical"})
    stats_up.to_csv(outdir/"stats_Unique_peptide_n_by_canonical.tsv", sep="\t")

    # 7) 总览
    summary_lines = []
    summary_lines.append(f"Input: {inp}")
    summary_lines.append(f"Rows considered (Unique_peptide_n>={args.min_unique}):\t{len(f)}")
    summary_lines.append(f"UniProt rows:\t{uniprot_n}")
    summary_lines.append(f"Canonical rows:\t{canonical_n}")
    summary_lines.append("\n[Counts by ORF_type]\n" + orf_type_counts.to_string())
    summary_lines.append("\n[Counts by Start_codon]\n" + start_counts.to_string())
    summary_lines.append("\n[ORF_length stats by canonical]\n" + stats_len.to_string())
    summary_lines.append("\n[Unique_peptide_n stats by canonical]\n" + stats_up.to_string())

    (outdir/"README.txt").write_text("\n".join(summary_lines), encoding="utf-8")
    print("\n".join(summary_lines))

if __name__ == "__main__":
    main()
