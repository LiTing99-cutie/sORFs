#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

def main():
    ap = argparse.ArgumentParser(description="Compute correlation between iBAQ_A and iBAQ_B.")
    ap.add_argument("--input", required=True, help="TSV from iBAQ computation (must contain iBAQ_A, iBAQ_B).")
    ap.add_argument("--log", action="store_true", help="Compute correlation in log10 space with small pseudocount.")
    ap.add_argument("--plot", default=None, help="Optional: path to save scatter PNG (log10 axes if --log).")
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    if not {"iBAQ_A","iBAQ_B"}.issubset(df.columns):
        raise SystemExit("输入缺少列：iBAQ_A / iBAQ_B")

    x = pd.to_numeric(df["iBAQ_A"], errors="coerce")
    y = pd.to_numeric(df["iBAQ_B"], errors="coerce")
    mask = x.notna() & y.notna() & (x > 0) & (y > 0)
    d = pd.DataFrame({"x": x[mask], "y": y[mask]})
    if len(d) < 3:
        print("[WARN] 有效点不足，无法计算相关。"); return

    if args.log:
        epsx = max(1e-12, np.nanquantile(d["x"], 0.001))
        epsy = max(1e-12, np.nanquantile(d["y"], 0.001))
        d["x"] = np.log10(d["x"] + epsx)
        d["y"] = np.log10(d["y"] + epsy)

    pear = d.corr(method="pearson").iloc[0,1]
    spear = d.corr(method="spearman").iloc[0,1]
    print(f"Pairs={len(d)}")
    print(f"Pearson r = {pear:.4f}")
    print(f"Spearman rho = {spear:.4f}")

    if args.plot:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(5.8,4.6))
        plt.scatter(d["x"], d["y"], s=6, alpha=0.35)
        k, b = np.polyfit(d["x"], d["y"], deg=1)
        xs = np.linspace(d["x"].min(), d["x"].max(), 100)
        plt.plot(xs, k*xs + b, linewidth=1.0)
        plt.xlabel("iBAQ_A (log10)" if args.log else "iBAQ_A")
        plt.ylabel("iBAQ_B (log10)" if args.log else "iBAQ_B")
        plt.title(f"iBAQ_A vs iBAQ_B (n={len(d)}, r={pear:.3f}, ρ={spear:.3f})")
        Path(args.plot).parent.mkdir(parents=True, exist_ok=True)
        plt.tight_layout(); plt.savefig(args.plot, dpi=300)
        print(f"[PLOT] {args.plot}")

if __name__ == "__main__":
    main()
