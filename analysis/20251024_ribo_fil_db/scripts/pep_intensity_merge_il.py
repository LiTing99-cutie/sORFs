#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
合并 closed/open 的肽段强度（I↔L 等价），并可按规则择优输出。

默认策略：
  - 同一肽段优先使用 CLOSED 的强度；只有 CLOSED 缺失时，才用 OPEN 补充。
可选开关：
  --treat-zero-as-na        把 0 当作缺失处理（更符合“未定量成功”的语义）
  --closed-zero-fallback    若 CLOSED=0 且 OPEN>0，用 OPEN（需配合你的实际偏好）
聚合方式：
  --agg {max,sum}           同一 search 内，重复肽段强度的聚合（默认 max）

输出：
  out.tsv:  Peptide_I_L_equal  Sample  Intensity  Source
  out.corr.txt:  重叠肽段的相关性报告（log10 强度）

用法示例：
python pep_intensity_merge_il.py \
  --closed /path/to/closed/peptide.tsv \
  --open   /path/to/open/peptide.tsv \
  --sample 21pcw_1_C8_T_T \
  --out    ../processed/quant/peptide_intensity_IL.tsv \
  --agg max --treat-zero-as-na --closed-zero-fallback
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

REQUIRED_COLS = {"Peptide", "Intensity"}

def load_pep_tsv(path: str, agg: str = "max") -> pd.DataFrame:
    """读取 peptide.tsv，取 Peptide/Intensity；做 I→L、聚合到肽段序列级"""
    df = pd.read_csv(path, sep="\t", dtype=str)
    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise SystemExit(f"[ERROR] {path} 缺少列: {', '.join(missing)}")

    pep_il = df["Peptide"].astype(str).str.upper().str.replace("I", "L", regex=False)
    inten = pd.to_numeric(df["Intensity"], errors="coerce").fillna(0.0)
    tmp = pd.DataFrame({"Peptide_IL": pep_il, "Intensity": inten})

    if agg == "max":
        out = tmp.groupby("Peptide_IL", as_index=False)["Intensity"].max()
    elif agg == "sum":
        out = tmp.groupby("Peptide_IL", as_index=False)["Intensity"].sum()
    else:
        raise SystemExit("--agg 仅支持 max/sum")
    return out

def compute_corr(closed_df: pd.DataFrame, open_df: pd.DataFrame):
    """在重叠肽上计算 log10(x+1) 的 Pearson 相关"""
    both = closed_df.merge(open_df, on="Peptide_IL", how="inner",
                           suffixes=("_closed", "_open"))
    if len(both) >= 3:
        x = np.log10(both["Intensity_closed"].to_numpy() + 1.0)
        y = np.log10(both["Intensity_open"].to_numpy() + 1.0)
        r, p = pearsonr(x, y)
    else:
        r, p = (np.nan, np.nan)
    return r, p, len(both)

def main():
    ap = argparse.ArgumentParser(description="合并 closed/open 肽段强度（I↔L 等价）并择优输出")
    ap.add_argument("--closed", required=True, help="closed search 的 peptide.tsv")
    ap.add_argument("--open",   required=True, help="open   search 的 peptide.tsv")
    ap.add_argument("--sample", required=True, help="样本名（写入输出第二列）")
    ap.add_argument("--out",    required=True, help="输出 TSV 文件路径")
    ap.add_argument("--agg",    choices=["max","sum"], default="max",
                    help="同一 search 内重复肽段的聚合方式（默认 max）")
    ap.add_argument("--treat-zero-as-na", action="store_true",
                    help="把 0 当作缺失（NaN）处理（默认否）")
    ap.add_argument("--closed-zero-fallback", action="store_true",
                    help="若 CLOSED 存在但强度=0 且 OPEN>0，则用 OPEN（默认否）")
    args = ap.parse_args()

    # 读取与聚合
    closed = load_pep_tsv(args.closed, agg=args.agg).rename(columns={"Intensity":"Intensity_closed"})
    open_  = load_pep_tsv(args.open,   agg=args.agg).rename(columns={"Intensity":"Intensity_open"})

    # 相关性（仅报告用，不参与决策）
    r, p, n_pair = compute_corr(closed, open_)

    # 合并并按规则选择
    merged = closed.merge(open_, on="Peptide_IL", how="outer")
    ic = pd.to_numeric(merged["Intensity_closed"], errors="coerce")
    io = pd.to_numeric(merged["Intensity_open"],  errors="coerce")

    # 可选：把 0 当作缺失（更贴合“未定量成功”语义）
    if args.treat_zero_as_na:
        ic = ic.mask(ic == 0)
        io = io.mask(io == 0)

    # 默认严格规则：优先 CLOSED，CLOSED 缺失才用 OPEN
    use = ic.copy()
    use = use.fillna(io)

    # 可选回退：CLOSED=0 且 OPEN>0 时用 OPEN（仅当 zero-as-na 未启用或 CLOSED=0 未被置 NaN 时才有意义）
    if args.closed_zero_fallback:
        cond = ic.notna() & (ic == 0) & io.notna() & (io > 0)
        use = use.where(~cond, io)

    merged["Intensity"] = use

    # 源标记：谁提供了最终强度
    src = np.where(ic.notna(), "closed",
                   np.where(io.notna(), "open", "NA"))
    # 若 treat-zero-as-na 导致 closed 变 NaN 而 open 有值，源需要同步调整
    if args.treat_zero_as_na:
        need_open = merged["Intensity"].isna() & io.notna()
        src = np.where(need_open, "open", src)
    merged["Source"] = src

    # 整理输出
    out = merged[["Peptide_IL"]].copy()
    out["Sample"] = args.sample
    out["Intensity"] = merged["Intensity"].astype(float)
    out["Source"] = merged["Source"]
    out = out.sort_values(["Peptide_IL"]).reset_index(drop=True)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)

    # 写相关性报告
    rep = out_path.with_suffix(".corr.txt")
    with open(rep, "w") as w:
        w.write(f"paired_peptides\t{n_pair}\n")
        w.write(f"pearson_r\t{r}\n")
        w.write(f"p_value\t{p}\n")
        w.write(f"agg\t{args.agg}\n")
        w.write(f"treat_zero_as_na\t{args.treat_zero_as_na}\n")
        w.write(f"closed_zero_fallback\t{args.closed_zero_fallback}\n")
        w.write("decision\tprefer_closed_else_open\n")

    # 终端摘要
    n_total = len(out)
    n_from_closed = int((out["Source"] == "closed").sum())
    n_from_open   = int((out["Source"] == "open").sum())
    print(f"Done. Out -> {out_path}")
    print(f"Correlation (paired n={n_pair}): r={r}, p={p}")
    print(f"Rows: total={n_total}, from_closed={n_from_closed}, from_open={n_from_open}")

if __name__ == "__main__":
    main()
