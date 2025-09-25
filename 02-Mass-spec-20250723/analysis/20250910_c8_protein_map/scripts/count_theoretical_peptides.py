#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from pathlib import Path
import pandas as pd

# ---------- 规则解析 ----------
def detect_rule(sample: str):
    s = sample.strip()
    if s.endswith(("T_T", "LC_T", "T_AC")):
        return {"cterm": set("KR"), "nterm": set(), "proline_block_cterm": True}
    if s.endswith("LN_T"):
        return {"cterm": set("KR"), "nterm": set("K"), "proline_block_cterm": True}
    if s.endswith("T_C"):
        return {"cterm": set("KRFLWY"), "nterm": set(), "proline_block_cterm": True}
    if s.endswith("T_GC"):
        return {"cterm": set("KRDE"), "nterm": set(), "proline_block_cterm": True}
    if s.endswith("T_AN"):
        return {"cterm": set("KR"), "nterm": set("D"), "proline_block_cterm": True}
    raise SystemExit(f"无法根据 Sample='{sample}' 识别酶切规则，请检查后缀。")

# ---------- 生成切点 ----------
def compute_cleavage_sites(seq: str, cterm:set, nterm:set, proline_block_cterm:bool=True):
    s = (seq or "").strip().upper()
    L = len(s)
    cuts = set([0, L])
    if L == 0:
        return sorted(cuts)

    # C 端切：after residue i；若启用阻断且下一个为 P，则不切
    for i in range(L-1):
        if s[i] in cterm:
            if proline_block_cterm and s[i+1] == "P":
                continue
            cuts.add(i+1)

    # N 端切：before residue i
    for i in range(L):
        if s[i] in nterm:
            cuts.add(i)

    return sorted(cuts)

# ---------- 统计允许 ≤2 漏切且长度在区间内的肽段数量 ----------
def count_peptides_with_missed_cleavages(cuts, min_len:int, max_len:int):
    """
    cuts: 递增切点列表（含0与len）
    允许漏切数 ≤2 → 终点索引 j ∈ {i+1,i+2,i+3}
    仅统计长度在 [min_len, max_len] 内的片段
    """
    n = len(cuts)
    if n <= 1:
        return 0
    cnt = 0
    for i in range(n-1):
        for step in (1, 2, 3):  # 0/1/2 漏切
            j = i + step
            if j < n:
                length = cuts[j] - cuts[i]
                if min_len <= length <= max_len:
                    cnt += 1
    return cnt

# ---------- 主流程 ----------
def main():
    ap = argparse.ArgumentParser(description="按 Sample 识别酶切规则，允许2漏切，统计理论肽段数（含长度筛选）")
    ap.add_argument("--in", dest="in_tsv", required=True,
                    help="输入TSV，至少含 ORF_id, ORF_seq；可选含 Sample")
    ap.add_argument("--sample", dest="sample", default=None,
                    help="若输入文件无 Sample 列，用此样本名作为全表规则选择")
    ap.add_argument("--out", dest="out_tsv", required=True,
                    help="输出TSV，新增 Theo_peptide_n 列")
    ap.add_argument("--min-len", type=int, default=6, help="肽段最小长度（默认6）")
    ap.add_argument("--max-len", type=int, default=50, help="肽段最大长度（默认50）")
    args = ap.parse_args()

    df = pd.read_csv(args.in_tsv, sep="\t", dtype=str)
    for c in ("ORF_id", "ORF_seq"):
        if c not in df.columns:
            raise SystemExit(f"输入缺少列：{c}")

    if "Sample" not in df.columns:
        if not args.sample:
            raise SystemExit("输入无 Sample 列，且未提供 --sample。")
        df["Sample"] = args.sample

    theo_counts = []
    for sample, seq in zip(df["Sample"].astype(str), df["ORF_seq"].astype(str)):
        rule = detect_rule(sample)
        cuts = compute_cleavage_sites(seq, rule["cterm"], rule["nterm"], rule["proline_block_cterm"])
        npep = count_peptides_with_missed_cleavages(cuts, args.min_len, args.max_len)
        theo_counts.append(npep)

    df["Theo_peptide_n"] = pd.Series(theo_counts, dtype="Int64")
    # 不输出 Sample 列
    df = df.drop(columns=["Sample"], errors="ignore")
    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[OK] 写出：{args.out_tsv}  行数={len(df)} 列数={df.shape[1]}  长度范围=[{args.min_len},{args.max_len}]")

if __name__ == "__main__":
    main()
