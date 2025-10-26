#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, sys, re
import pandas as pd
from pathlib import Path

def is_uniprot(pid: str) -> bool:
    s = str(pid).strip()
    return s.startswith("sp|") or s.startswith("tr|")

def parse_pb(pid: str):
    """从 PB 样式 ID 解析：isoform_id、orf_type、start_codon"""
    s = str(pid).strip()
    parts = s.split("|")
    if len(parts) >= 5:
        orf_type = parts[-2].strip()
        start_codon = parts[-1].strip().strip(",").upper() or "NA"
        left = parts[0]
        iso = left.split(":", 1)[0]  # PB.x.y
        return iso, orf_type, start_codon
    return None, "unknown", "NA"

def load_pep2prot(path):
    """读肽-蛋白映射：允许无表头（两列），或有表头（前两列用）"""
    df = pd.read_csv(path, sep="\t", header=None, dtype=str, comment="#")
    if df.shape[1] < 2:
        raise SystemExit("pep2prot 输入需至少两列")
    # 若首行像表头，则去表头
    h1 = (df.iloc[0,0] or "").lower()
    h2 = (df.iloc[0,1] or "").lower()
    looks_header = ("peptide" in h1) or ("protein" in h2)
    if looks_header:
        df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
        df = df.iloc[:, :2]
        df.columns = ["Peptide", "Protein_id"]
    else:
        df = df.iloc[:, :2]
        df.columns = ["Peptide", "Protein_id"]
    df["Peptide"] = df["Peptide"].astype(str).str.strip()
    df["Protein_id"] = df["Protein_id"].astype(str).str.strip()
    return df

def load_isoform(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    # 需要 isoform 与 associated_gene；FL_TPM 列可选
    if "isoform" not in df.columns:
        raise SystemExit("isoform 表必须包含列：isoform")
    if "associated_gene" not in df.columns:
        # 允许缺失，但筛选时对 PB 的 C 需要 Geneid，最好有
        df["associated_gene"] = ""
    # 找 FL_TPM 列（如果没有就补0）
    tpm_col = None
    for c in df.columns:
        if "FL_TPM" in c:
            tpm_col = c
            break
    if tpm_col is None:
        tpm_col = "FL_TPM"
        df[tpm_col] = "0"
    df = df[["isoform", "associated_gene", tpm_col]].copy()
    df.rename(columns={tpm_col: "FL_TPM"}, inplace=True)
    df["FL_TPM"] = pd.to_numeric(df["FL_TPM"], errors="coerce").fillna(0.0)
    return df

def load_orf_psites(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    need = ["ORF_id", "RPF_reads", "Psites_number"]
    for c in need:
        if c not in df.columns:
            raise SystemExit(f"ORF表缺少列：{c}")
    df["RPF_reads"] = pd.to_numeric(df["RPF_reads"], errors="coerce").fillna(0).astype(int)
    df["Psites_number"] = pd.to_numeric(df["Psites_number"], errors="coerce").fillna(0).astype(int)
    return df[need].copy()

def load_nca(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    need = ["Geneid", "N", "C", "A"]
    for c in need:
        if c not in df.columns:
            raise SystemExit(f"N/C/A 表缺少列：{c}")
    # 数值列
    for c in ["N","C","A"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    return df[need].copy()

def build_candidates(pep2prot, iso_tbl, orf_tbl, nca_tbl):
    """合并各层信息得到候选表"""
    rows = []
    for _, r in pep2prot.iterrows():
        pep = r["Peptide"]
        pid = r["Protein_id"]
        uni = is_uniprot(pid)
        if uni:
            rows.append({
                "Peptide": pep,
                "Protein_id": pid,
                "is_uniprot": True,
                "isoform": "",
                "orf_type": "canonical",
                "start_codon": "ATG",
                "associated_gene": "",
                "FL_TPM": "",
                "RPF_reads": "",
                "Psites_number": "",
                "Geneid": "",
                "N": "",
                "C": "",
                "A": "",
                "is_canonical": True
            })
        else:
            iso, orf_type, start_codon = parse_pb(pid)
            # isoform 连接
            sub_iso = iso_tbl[iso_tbl["isoform"] == (iso or "")]
            associated_gene = sub_iso["associated_gene"].iloc[0] if len(sub_iso) else ""
            fl_tpm = sub_iso["FL_TPM"].iloc[0] if len(sub_iso) else 0.0

            # ORF连接
            sub_orf = orf_tbl[orf_tbl["ORF_id"] == pid]
            rpf = int(sub_orf["RPF_reads"].iloc[0]) if len(sub_orf) else 0
            psites = int(sub_orf["Psites_number"].iloc[0]) if len(sub_orf) else 0

            # Gene→NCA
            sub_nca = nca_tbl[nca_tbl["Geneid"] == associated_gene] if associated_gene else pd.DataFrame()
            N = float(sub_nca["N"].iloc[0]) if len(sub_nca) else 0.0
            C = float(sub_nca["C"].iloc[0]) if len(sub_nca) else 0.0
            A = float(sub_nca["A"].iloc[0]) if len(sub_nca) else 0.0

            rows.append({
                "Peptide": pep,
                "Protein_id": pid,
                "is_uniprot": False,
                "isoform": iso or "",
                "orf_type": orf_type,
                "start_codon": start_codon,
                "associated_gene": associated_gene,
                "FL_TPM": fl_tpm,
                "RPF_reads": rpf,
                "Psites_number": psites,
                "Geneid": associated_gene or "",
                "N": N, "C": C, "A": A,
                "is_canonical": (orf_type.lower() == "canonical")
            })
    cand = pd.DataFrame(rows)
    # 排序稳定
    cand = cand[[
        "Peptide","Protein_id","is_uniprot","is_canonical","orf_type","start_codon",
        "isoform","associated_gene","Geneid","FL_TPM","RPF_reads","Psites_number","N","C","A"
    ]]
    return cand

def filter_candidates(cand, min_psites, min_c):
    """PB 依据 (Psites>=min_psites & C>=min_c)；UniProt 不筛、证据信息置空"""
    keep_pb = (~cand["is_uniprot"]) & \
              (pd.to_numeric(cand["Psites_number"], errors="coerce").fillna(0) >= min_psites) & \
              (pd.to_numeric(cand["C"], errors="coerce").fillna(0.0) >= min_c)
    keep_uni = cand["is_uniprot"]
    kept = cand[keep_pb | keep_uni].copy()

    # 对 UniProt 置空证据列
    for col in ["FL_TPM","RPF_reads","Psites_number","N","C","A","associated_gene","isoform","start_codon","orf_type","Geneid"]:
        kept.loc[kept["is_uniprot"], col] = ""
    return kept

def unique_assign(df):
    """按肽做唯一指派：仅1条 或 canonical唯一 ⇒ 唯一；返回唯一肽表与唯一蛋白集合大小"""
    uniq_rows = []
    for pep, grp in df.groupby("Peptide", sort=False):
        n = len(grp)
        if n == 1:
            row = grp.iloc[0]
            uniq_rows.append({"Peptide": pep, "Assigned_protein_id": row["Protein_id"], "reason": "single_candidate"})
            continue
        # canonical 唯一？
        canon = grp[grp["is_canonical"] == True]
        if len(canon) == 1:
            row = canon.iloc[0]
            uniq_rows.append({"Peptide": pep, "Assigned_protein_id": row["Protein_id"], "reason": "unique_canonical"})
    uniq_df = pd.DataFrame(uniq_rows)
    unique_protein_count = uniq_df["Assigned_protein_id"].nunique() if not uniq_df.empty else 0
    return uniq_df, unique_protein_count

def main():
    ap = argparse.ArgumentParser(description="Peptide→Protein 指派：合并多组学，筛选（PB用P-sites与C），UniProt跳过筛选并置空证据；筛前/筛后唯一统计。")
    ap.add_argument("--pep2prot", required=True, help="肽-蛋白映射TSV（两列：Peptide, Protein_id；可无表头）")
    ap.add_argument("--isoform", required=True, help="isoform表达TSV（需含 isoform, associated_gene；FL_TPM 列自动探测）")
    ap.add_argument("--orf_psites", required=True, help="ORF证据TSV（列：ORF_id, RPF_reads, Psites_number）")
    ap.add_argument("--nca", required=True, help="核/胞表达TSV（列：Geneid, N, C, A）")
    ap.add_argument("--min_psites", type=int, default=6, help="PB保留的最小Psites_number（默认6）")
    ap.add_argument("--min_c", type=float, default=0.2, help="PB保留的最小核外表达C（默认0.2）")
    ap.add_argument("--outdir", required=True, help="输出目录")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    pep2prot = load_pep2prot(args.pep2prot)
    iso_tbl = load_isoform(args.isoform)
    orf_tbl = load_orf_psites(args.orf_psites)
    nca_tbl = load_nca(args.nca)

    # 合并候选
    cand = build_candidates(pep2prot, iso_tbl, orf_tbl, nca_tbl)
    cand_pre_path = outdir / "candidates.pre.tsv"
    cand.to_csv(cand_pre_path, sep="\t", index=False)

    # 筛选（PB用Psites与C；UniProt不筛并置空证据信息）
    kept = filter_candidates(cand, args.min_psites, args.min_c)
    kept_path = outdir / "candidates.filtered.tsv"
    kept.to_csv(kept_path, sep="\t", index=False)

    # 唯一指派（筛前/筛后）
    uniq_pre, uniq_prot_pre = unique_assign(cand)
    uniq_post, uniq_prot_post = unique_assign(kept)
    uniq_pre_path = outdir / "assignments.unique.pre.tsv"
    uniq_post_path = outdir / "assignments.unique.post.tsv"
    uniq_pre.to_csv(uniq_pre_path, sep="\t", index=False)
    uniq_post.to_csv(uniq_post_path, sep="\t", index=False)

    # 统计（你要求“筛选前和筛选后唯一比对的 ORF id 的数量”）
    # 这里“ORF id 数量”按唯一被指派的 Protein_id 去重计数（包含 PB 与 UniProt 的“全量”口径）
    # 如需“PB-only”可再加一列统计
    pb_only_pre = uniq_pre[~uniq_pre["Assigned_protein_id"].str.startswith(("sp|","tr|"))]["Assigned_protein_id"].nunique()
    pb_only_post = uniq_post[~uniq_post["Assigned_protein_id"].str.startswith(("sp|","tr|"))]["Assigned_protein_id"].nunique()

    with open(outdir / "summary.txt", "w") as w:
        w.write(f"Unique peptides (pre) : {len(uniq_pre)}\n")
        w.write(f"Unique peptides (post): {len(uniq_post)}\n")
        w.write(f"Unique proteins (pre, all) : {uniq_prot_pre}\n")
        w.write(f"Unique proteins (post, all): {uniq_prot_post}\n")
        w.write(f"Unique PB-ORFs (pre)  : {pb_only_pre}\n")
        w.write(f"Unique PB-ORFs (post) : {pb_only_post}\n")

    print("Done.")
    print(f"- {cand_pre_path}")
    print(f"- {kept_path}")
    print(f"- {uniq_pre_path}")
    print(f"- {uniq_post_path}")
    print(f"- {outdir/'summary.txt'}")

if __name__ == "__main__":
    main()
