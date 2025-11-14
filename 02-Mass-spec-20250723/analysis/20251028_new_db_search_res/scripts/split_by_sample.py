#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, sys, re
from pathlib import Path

def find_sample_idx(header_cols):
    for i, c in enumerate(header_cols):
        if c.strip().lower() == "sample":
            return i
    sys.exit("未找到 sample 列（大小写皆可）。")

def split_stream(path, out_dir, out_name):
    path = Path(path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    samples = set()
    bad = 0
    handles = {}  # sample -> file handle

    def open_writer(samp):
        sd = out_dir / samp
        sd.mkdir(parents=True, exist_ok=True)
        fn = sd / out_name
        exists = fn.exists()
        fh = open(fn, "a", encoding="utf-8", newline="")
        if not exists:
            fh.write(header_line + "\n")
        return fh

    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        try:
            header_line = f.readline().rstrip("\n\r")
        except Exception as e:
            sys.exit(f"[FATAL] 读取表头失败: {e}")
        if not header_line:
            sys.exit("[FATAL] 空文件或表头缺失: " + str(path))
        header_cols = header_line.split("\t")
        ncol = len(header_cols)
        sample_idx = find_sample_idx(header_cols)

        for ln, line in enumerate(f, start=2):
            line = line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != ncol:
                # 尝试宽松修复：把多余的列并到最后一列
                if len(parts) > ncol:
                    parts = parts[:ncol-1] + ["\t".join(parts[ncol-1:])]
                else:
                    bad += 1
                    continue

            samp = parts[sample_idx].strip()
            if not samp:
                continue
            # 目录名安全
            samp_dir = re.sub(r"[^A-Za-z0-9._-]+", "_", samp)

            if samp_dir not in handles:
                handles[samp_dir] = open_writer(samp_dir)
            handles[samp_dir].write(line + "\n")
            samples.add(samp)

    for fh in handles.values():
        try:
            fh.close()
        except:
            pass

    if bad:
        print(f"[WARN] {path.name}: 跳过列数不匹配的行数 = {bad}", file=sys.stderr)
    return sorted(samples)

def main():
    ap = argparse.ArgumentParser(description="按 sample 列拆分两个 merged 输入为每样本文件，并输出样本清单（纯Python流式版）")
    ap.add_argument("--pep_orf_merged", required=True, help=".../pep.orf.merged.txt")
    ap.add_argument("--intensity_merged", required=True, help=".../peptide_intensity_IL.merged.tsv")
    ap.add_argument("--out_base", required=True, help="输出基目录，例如 ../processed/by_sample")
    ap.add_argument("--sample_list", required=True, help="输出样本清单 txt")
    args = ap.parse_args()

    out_base = Path(args.out_base)
    sp1 = split_stream(args.pep_orf_merged,  out_base, "pep.orf.txt")
    sp2 = split_stream(args.intensity_merged, out_base, "peptide_intensity_IL.tsv")

    s1, s2 = set(sp1), set(sp2)
    samples = sorted(s1 | s2)
    Path(args.sample_list).write_text("\n".join(samples) + "\n", encoding="utf-8")

    only1 = sorted(s1 - s2)
    only2 = sorted(s2 - s1)
    if only1:
        print(f"[WARN] 仅 pep.orf.merged 中存在的样本: {only1}")
    if only2:
        print(f"[WARN] 仅 intensity.merged 中存在的样本: {only2}")
    print(f"[OK] 拆分完成。样本数={len(samples)}，清单输出：{args.sample_list}")

if __name__ == "__main__":
    main()
