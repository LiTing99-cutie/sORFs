#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Peptide → Protein 指派（向量化版，支持自选RPF/PS指标）
- 合并多组学（isoform 表达、ORF RPF/P-sites、N/C/A）
- 生物学筛选（仅 PB-ORF）：
  (RPF_selected>=rpf_min OR Psites_selected>=ps_min) AND C>=min_c
  * UniProt(sp|/tr|)不参与筛选，证据列置空
- 唯一指派：逐肽 (候选=1) 或 (canonical 唯一)
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

# CHG: 支持仅按需读取所选指标列
def load_orf_psites(path: str, rpf_metric: str, ps_metric: str) -> pd.DataFrame:
    """
    支持 TSV 或 Parquet；按需读取 ORF_id + 选择的两列
    rpf_metric ∈ {"RPF_reads","RPF_RPKM","RPF_codon_coverage"}
    ps_metric  ∈ {"Psites_number","Psites_RPKM","Psites_codon_coverage"}
    """
    need_cols = ["ORF_id", rpf_metric, ps_metric]
    if path.lower().endswith(".parquet"):
        try:
            df = pd.read_parquet(path, columns=need_cols)
        except Exception:
            df = pd.read_parquet(path)
    else:
        df = _read_csv_fast(path, usecols=None, dtypes={"ORF_id":"string"})
    for c in need_cols:
        if c not in df.columns:
            raise SystemExit(f"ORF表缺少列：{c}")
    # 数值类型转换：reads/number→整数；RPKM/coverage→浮点
    def _to_num(s, as_int: bool):
        x = pd.to_numeric(s, errors="coerce").fillna(0)
        return x.astype("int64") if as_int else x.astype(float)

    int_like = { "RPF_reads", "Psites_number" }
    df_out = df[need_cols].copy()
    df_out[rpf_metric] = _to_num(df_out[rpf_metric], rpf_metric in int_like)
    df_out[ps_metric]  = _to_num(df_out[ps_metric],  ps_metric  in int_like)
    return df_out

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
        description="向量化：合并→筛选[(PB)(RPF>=rpf_min OR Psites>=ps_min) AND C>=min_c；UniProt跳筛]→唯一指派→统计"
    )
    ap.add_argument("--pep2prot", required=True, help="肽-蛋白TSV（两列，可无表头）")
    ap.add_argument("--isoform", required=True, help="isoform表达TSV（含 isoform, associated_gene；FL_TPM*）")
    ap.add_argument("--orf_psites", required=True, help="ORF证据（TSV/Parquet；含所选列）")
    ap.add_argument("--nca", required=True, help="N/C/A 表（TSV；Geneid,N,C,A；C=核外，A=均值）")

    # NEW: 指标与阈值
    ap.add_argument("--rpf_metric",
                    choices=["RPF_reads","RPF_RPKM","RPF_codon_coverage"],
                    default="RPF_reads",
                    help="选择 RPF 指标用于筛选（默认 RPF_reads）")
    ap.add_argument("--ps_metric",
                    choices=["Psites_number","Psites_RPKM","Psites_codon_coverage"],
                    default="Psites_number",
                    help="选择 P-sites 指标用于筛选（默认 Psites_number）")
    ap.add_argument("--rpf_min", type=float, default=6, help="RPF 指标最小阈值（默认6）")
    ap.add_argument("--ps_min",  type=float, default=2, help="P-sites 指标最小阈值（默认2）")

    ap.add_argument("--min_c", type=float, default=0.2, help="PB 最小核外表达 C（默认0.2）")
    ap.add_argument("--outdir", required=True, help="输出目录")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    # 读取
    pep2prot = load_pep2prot(args.pep2prot)
    iso_tbl  = load_isoform(args.isoform)
    orf_tbl  = load_orf_psites(args.orf_psites, args.rpf_metric, args.ps_metric)  # CHG
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

    # 合并（ORF证据按 Protein_id ↔ ORF_id）
    cand = pep2prot.merge(iso_tbl, how="left", on="isoform") \
                   .merge(orf_tbl, how="left", left_on="Protein_id", right_on="ORF_id") \
                   .drop(columns=["ORF_id"], errors="ignore")
    cand = cand.merge(nca_tbl, how="left", left_on="associated_gene", right_on="Geneid") \
               .drop(columns=["Geneid_y"], errors="ignore")
    if "Geneid_x" in cand.columns:
        cand.rename(columns={"Geneid_x":"Geneid"}, inplace=True)

    # 统一命名选中指标以便筛选与导出
    cand.rename(columns={
        args.rpf_metric: "RPF_selected",
        args.ps_metric:  "Psites_selected"
    }, inplace=True)

    # 数值列
    num_defaults = {
        "FL_TPM": 0.0, "N": 0.0, "C": 0.0, "A": 0.0,
        "RPF_selected": 0.0, "Psites_selected": 0.0
    }
    for col, default in num_defaults.items():
        if col in cand.columns:
            cand[col] = pd.to_numeric(cand[col], errors="coerce").fillna(default)

    # 保存筛选前候选
    cols_order = ["Peptide","Protein_id","is_uniprot","is_canonical","orf_type","start_codon",
                  "isoform","associated_gene","Geneid","FL_TPM",
                  "RPF_selected","Psites_selected","N","C","A"]
    # 尽量保留原始指标名列（若与 selected 不同）
    if args.rpf_metric != "RPF_selected" and args.rpf_metric in pep2prot.columns or args.rpf_metric in cand.columns:
        pass  # 已rename为 RPF_selected
    if args.ps_metric != "Psites_selected" and args.ps_metric in cand.columns:
        pass
    cand = cand.reindex(columns=cols_order)
    cand.to_csv(outdir/"candidates.pre.tsv", sep="\t", index=False)

    # ===== 生物学筛选（更新：使用 RPF_selected / Psites_selected） =====
    keep_pb  = (~cand["is_uniprot"]) & (
                  (cand["RPF_selected"] >= args.rpf_min) | (cand["Psites_selected"] >= args.ps_min)
               ) & (cand["C"] >= args.min_c)
    keep_uni = cand["is_uniprot"]
    kept = cand[keep_pb | keep_uni].copy()

    # UniProt 证据信息置空（尽量覆盖证据相关列）
    evidence_cols = ["FL_TPM","N","C","A","RPF_selected","Psites_selected",
                     "associated_gene","isoform","start_codon","orf_type","Geneid"]
    for c in evidence_cols:
        if c in kept.columns:
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

    # Unique non-canonical 的 PB-ORFs
    def count_noncanonical_pb(uniq_df: pd.DataFrame) -> int:
        s = uniq_df["Assigned_protein_id"].astype(str)
        mask_pb = s.str.startswith(("PB.", "PB|"))
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
