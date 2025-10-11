#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Peptide → Protein 指派（向量化版）
- 合并多组学（isoform 表达、ORF RPF/P-sites、N/C/A）
- 生物学筛选（仅 PB-ORF）：(RPF_reads>=min_rpf OR Psites>=min_psites_alt) AND C>=min_c
  * UniProt(sp|/tr|)不参与筛选，证据列置空
- 唯一指派：逐肽 (候选=1) 或 (canonical 唯一) 视为唯一；不做进一步强行择一
- 产出：候选(筛前/筛后)、唯一(筛前/筛后)、summary
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

# ============== 读取与预处理 ==============

def load_pep2prot(path: str) -> pd.DataFrame:
    """读肽-蛋白映射（TSV，前两列：Peptide, Protein_id；可无表头）"""
    df0 = pd.read_csv(path, sep="\t", header=None, dtype=str, comment="#")
    if df0.shape[1] < 2:
        raise SystemExit("pep2prot 需至少两列")
    h1 = (df0.iloc[0,0] or "").lower()
    h2 = (df0.iloc[0,1] or "").lower()
    looks_header = ("peptide" in h1) or ("protein" in h2)
    if looks_header:
        df = pd.read_csv(path, sep="\t", dtype=str, comment="#").iloc[:, :2]
        df.columns = ["Peptide", "Protein_id"]
    else:
        df = df0.iloc[:, :2].copy()
        df.columns = ["Peptide", "Protein_id"]
    df["Peptide"] = df["Peptide"].astype(str).str.strip()
    df["Protein_id"] = df["Protein_id"].astype(str).str.strip()
    return df

def load_isoform(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    if "isoform" not in df.columns:
        raise SystemExit("isoform 表缺少列 isoform")
    if "associated_gene" not in df.columns:
        df["associated_gene"] = ""
    tpm_col = next((c for c in df.columns if "FL_TPM" in c), None)
    if tpm_col is None:
        tpm_col = "FL_TPM"
        df[tpm_col] = "0"
    out = df[["isoform", "associated_gene", tpm_col]].copy()
    out.rename(columns={tpm_col: "FL_TPM"}, inplace=True)
    out["FL_TPM"] = pd.to_numeric(out["FL_TPM"], errors="coerce").fillna(0.0)
    return out

def _read_csv_fast(path, usecols=None, dtypes=None):
    try:
        return pd.read_csv(path, sep="\t", usecols=usecols, dtype=dtypes,
                           engine="pyarrow", memory_map=True)
    except Exception:
        return pd.read_csv(path, sep="\t", usecols=usecols, dtype=dtypes)

def load_orf_psites(path: str) -> pd.DataFrame:
    """支持 TSV 或 Parquet；仅读三列"""
    cols = ["ORF_id", "RPF_reads", "Psites_number"]
    if path.lower().endswith(".parquet"):
        try:
            df = pd.read_parquet(path, columns=cols)
        except Exception:
            df = pd.read_parquet(path)
    else:
        df = _read_csv_fast(path, usecols=cols, dtypes={"ORF_id":"string"})
    for c in cols:
        if c not in df.columns:
            raise SystemExit(f"ORF表缺少列：{c}")
    df["RPF_reads"] = pd.to_numeric(df["RPF_reads"], errors="coerce").fillna(0).astype("int32")
    df["Psites_number"] = pd.to_numeric(df["Psites_number"], errors="coerce").fillna(0).astype("int32")
    return df[cols].copy()

def load_nca(path: str) -> pd.DataFrame:
    df = _read_csv_fast(path, usecols=["Geneid","N","C","A"])
    for c in ["N","C","A"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    return df

# ============== 唯一指派（向量化） ==============

def unique_assign_vec(df: pd.DataFrame) -> pd.DataFrame:
    """按肽唯一指派：候选=1 或 canonical唯一"""
    if df.empty:
        return pd.DataFrame(columns=["Peptide","Assigned_protein_id","reason"])
    grp_size  = df.groupby("Peptide")["Protein_id"].transform("size")
    canon_cnt = df.assign(_canon=df["is_canonical"].astype(int)) \
                  .groupby("Peptide")["_canon"].transform("sum")
    mask_single = grp_size.eq(1)
    mask_unique_canon = df["is_canonical"] & canon_cnt.eq(1)

    uniq = df.loc[mask_single | mask_unique_canon, ["Peptide","Protein_id"]].copy()
    reason = np.where(mask_single, "single_candidate",
                      np.where(mask_unique_canon, "unique_canonical", ""))
    uniq["reason"] = reason[mask_single | mask_unique_canon]
    uniq.rename(columns={"Protein_id":"Assigned_protein_id"}, inplace=True)
    uniq = uniq.drop_duplicates(subset=["Peptide"], keep="first")
    return uniq

# ============== 主流程 ==============

def main():
    ap = argparse.ArgumentParser(
        description="向量化版：合并→筛选[(PB) (RPF>=min_rpf OR Psites>=min_psites_alt) AND C>=min_c；UniProt跳筛]→唯一指派→统计"
    )
    ap.add_argument("--pep2prot", required=True, help="肽-蛋白TSV（两列，可无表头）")
    ap.add_argument("--isoform", required=True, help="isoform表达TSV（含 isoform, associated_gene；FL_TPM*）")
    ap.add_argument("--orf_psites", required=True, help="ORF证据（TSV/Parquet；ORF_id,RPF_reads,Psites_number）")
    ap.add_argument("--nca", required=True, help="N/C/A 表（TSV；Geneid,N,C,A；C=核外，A=均值）")
    # 新参数：min_rpf 与 min_psites_alt
    ap.add_argument("--min_rpf", type=float, default=6, help="PB 最小 RPF_reads（默认6）")
    ap.add_argument("--min_psites_alt", type=float, default=2, help="PB 最小 Psites（或条件，默认2）")
    ap.add_argument("--min_c", type=float, default=0.2, help="PB 最小核外表达 C（默认0.2）")
    ap.add_argument("--outdir", required=True, help="输出目录")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    # 读取
    pep2prot = load_pep2prot(args.pep2prot)
    iso_tbl  = load_isoform(args.isoform)
    orf_tbl  = load_orf_psites(args.orf_psites)
    nca_tbl  = load_nca(args.nca)

    # 解析 ID（对子集 loc 赋值）
    pep2prot["is_uniprot"] = pep2prot["Protein_id"].str.startswith(("sp|","tr|"))
    mask_pb = ~pep2prot["is_uniprot"]

    pep2prot["isoform"]     = ""
    pep2prot["orf_type"]    = "canonical"  # UniProt 视为 canonical
    pep2prot["start_codon"] = "ATG"

    pb_rsplit = pep2prot.loc[mask_pb, "Protein_id"].str.rsplit("|", n=2, expand=True)
    pb_rsplit = pb_rsplit.reindex(columns=[0,1,2]).fillna("")
    pb_iso      = pb_rsplit[0].str.split(":", n=1, expand=True)[0]
    pb_orf_type = pb_rsplit[1].astype(str)
    pb_start    = pb_rsplit[2].astype(str).str.replace(r",$", "", regex=True).str.upper()

    pep2prot.loc[mask_pb, "isoform"]     = pb_iso.values
    pep2prot.loc[mask_pb, "orf_type"]    = pb_orf_type.values
    pep2prot.loc[mask_pb, "start_codon"] = pb_start.values

    pep2prot["is_canonical"] = pep2prot["is_uniprot"] | (pep2prot["orf_type"].str.lower() == "canonical")

    # 合并
    cand = pep2prot.merge(iso_tbl, how="left", on="isoform") \
                   .merge(orf_tbl, how="left", left_on="Protein_id", right_on="ORF_id") \
                   .drop(columns=["ORF_id"], errors="ignore")
    cand = cand.merge(nca_tbl, how="left", left_on="associated_gene", right_on="Geneid") \
               .drop(columns=["Geneid_y"], errors="ignore")
    if "Geneid_x" in cand.columns:
        cand.rename(columns={"Geneid_x":"Geneid"}, inplace=True)

    # 数值列
    for col, default in [("FL_TPM", 0.0), ("RPF_reads", 0.0), ("Psites_number", 0.0),
                         ("N", 0.0), ("C", 0.0), ("A", 0.0)]:
        if col in cand.columns:
            cand[col] = pd.to_numeric(cand[col], errors="coerce").fillna(default)

    # 保存筛选前候选
    cols_order = ["Peptide","Protein_id","is_uniprot","is_canonical","orf_type","start_codon",
                  "isoform","associated_gene","Geneid","FL_TPM","RPF_reads","Psites_number","N","C","A"]
    cand = cand.reindex(columns=cols_order)
    (outdir/"candidates.pre.tsv").write_text("")
    cand.to_csv(outdir/"candidates.pre.tsv", sep="\t", index=False)

    # ===== 生物学筛选（更新规则） =====
    # keep_pb: (~is_uniprot) AND ( (RPF>=min_rpf) OR (Psites>=min_psites_alt) ) AND (C>=min_c)
    keep_pb  = (~cand["is_uniprot"]) & (
                  (cand["RPF_reads"] >= args.min_rpf) | (cand["Psites_number"] >= args.min_psites_alt)
               ) & (cand["C"] >= args.min_c)
    keep_uni = cand["is_uniprot"]
    kept = cand[keep_pb | keep_uni].copy()

    # UniProt 证据信息置空
    empty_cols = ["FL_TPM","RPF_reads","Psites_number","N","C","A",
                  "associated_gene","isoform","start_codon","orf_type","Geneid"]
    for c in empty_cols:
        kept.loc[kept["is_uniprot"], c] = ""

    kept.to_csv(outdir/"candidates.filtered.tsv", sep="\t", index=False)

    # 唯一指派（筛前/筛后）
    uniq_pre  = unique_assign_vec(cand)
    uniq_post = unique_assign_vec(kept)
    uniq_pre.to_csv(outdir/"assignments.unique.pre.tsv",  sep="\t", index=False)
    uniq_post.to_csv(outdir/"assignments.unique.post.tsv", sep="\t", index=False)

    # 统计
    uniq_prot_pre  = uniq_pre["Assigned_protein_id"].nunique()
    uniq_prot_post = uniq_post["Assigned_protein_id"].nunique()

    pb_only_pre  = uniq_pre[~uniq_pre["Assigned_protein_id"].str.startswith(("sp|","tr|"))]["Assigned_protein_id"].nunique()
    pb_only_post = uniq_post[~uniq_post["Assigned_protein_id"].str.startswith(("sp|","tr|"))]["Assigned_protein_id"].nunique()

    # 新增：Unique non-canonical 的 PB-ORFs
    def count_noncanonical_pb(uniq_df: pd.DataFrame) -> int:
        s = uniq_df["Assigned_protein_id"].astype(str)
        mask_pb = s.str.startswith(("PB.", "PB|"))  # 兼容 PB. 和 PB| 两种前缀
        # orf_type 在 PB ID 倒数第二段；用字符串捕捉 '|canonical|' 来区分
        mask_noncanon = ~s.str.contains(r"\|canonical\|", case=False, regex=True)
        return s[mask_pb & mask_noncanon].nunique()

    nc_pb_pre  = count_noncanonical_pb(uniq_pre)
    nc_pb_post = count_noncanonical_pb(uniq_post)

    with open(outdir/"summary.txt","w") as w:
        w.write(f"Unique peptides (pre) : {len(uniq_pre)}\n")
        w.write(f"Unique peptides (post): {len(uniq_post)}\n")
        w.write(f"Unique proteins (pre, all) : {uniq_prot_pre}\n")
        w.write(f"Unique proteins (post, all): {uniq_prot_post}\n")
        w.write(f"Unique PB-ORFs (pre)  : {pb_only_pre}\n")
        w.write(f"Unique PB-ORFs (post) : {pb_only_post}\n")
        w.write(f"Unique non-canonical PB-ORFs (pre)  : {nc_pb_pre}\n")
        w.write(f"Unique non-canonical PB-ORFs (post) : {nc_pb_post}\n")

    print("Done.")
    print(f"- {outdir/'candidates.pre.tsv'}")
    print(f"- {outdir/'candidates.filtered.tsv'}")
    print(f"- {outdir/'assignments.unique.pre.tsv'}")
    print(f"- {outdir/'assignments.unique.post.tsv'}")
    print(f"- {outdir/'summary.txt'}")

if __name__ == "__main__":
    main()
