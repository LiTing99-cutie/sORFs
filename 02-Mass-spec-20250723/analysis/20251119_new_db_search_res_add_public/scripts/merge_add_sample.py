#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, sys, re, os, gzip, io
from pathlib import Path

def iter_files(root: Path, name: str):
    for p in root.rglob(name):
        if p.is_file():
            yield p

def open_any(path: Path):
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")

def get_sample_from_path(p: Path, mode: str, depth: int, regex: str):
    if mode == "parent":
        return p.parent.name
    if mode == "grandparent":
        return p.parent.parent.name
    if mode == "dirnameN":
        # depth=1 表示父目录，2 表示祖父目录……
        q = p
        for _ in range(depth):
            q = q.parent
        return q.name
    if mode == "regex":
        m = re.search(regex, str(p))
        if not m or "sample" not in m.groupdict():
            raise SystemExit(f"[ERR] regex 未匹配到命名分组 (?P<sample>...)：{p}")
        return m.group("sample")
    raise SystemExit(f"[ERR] 未知 mode: {mode}")

def main():
    ap = argparse.ArgumentParser(
        description="合并多 TSV，仅保留一个表头，并添加 Sample 列（从文件路径提取样本名）"
    )
    ap.add_argument("--root", required=True, help="搜索根目录")
    ap.add_argument("--name", required=True, help="文件名或模式（如 ibaq_b_with_total.tsv 或 *.tsv）")
    ap.add_argument("--out", required=True, help="输出 TSV 路径")
    ap.add_argument("--mode", choices=["parent","grandparent","dirnameN","regex"],
                    default="grandparent",
                    help="样本名提取方式（默认 grandparent：取文件的祖父目录名）")
    ap.add_argument("--depth", type=int, default=2,
                    help="当 mode=dirnameN 时，第 N 层目录作为样本名（1=父目录，2=祖父目录；默认2）")
    ap.add_argument("--regex", default=r".*/(?P<sample>[^/]+)/.*",
                    help="当 mode=regex 时用于提取样本名的正则，需包含 (?P<sample>...) 命名分组")
    args = ap.parse_args()

    root = Path(args.root)
    files = sorted(iter_files(root, args.name))
    if not files:
        raise SystemExit(f"[ERR] 未找到文件：root={root} name={args.name}")

    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)

    # 读第一个文件头，写入头部 + Sample 列
    with open_any(files[0]) as f0, open(outp, "w", encoding="utf-8", newline="") as w:
        header = f0.readline().rstrip("\n\r")
        if not header:
            raise SystemExit(f"[ERR] 首个文件无表头：{files[0]}")
        w.write(header + "\tSample\n")

        # 写首个文件的内容（跳过表头）
        sample0 = get_sample_from_path(files[0], args.mode, args.depth, args.regex)
        for line in f0:
            line = line.rstrip("\n\r")
            if line:
                w.write(f"{line}\t{sample0}\n")

        # 写其余文件
        for fp in files[1:]:
            sample = get_sample_from_path(fp, args.mode, args.depth, args.regex)
            with open_any(fp) as fr:
                # 丢弃各自的首行表头
                _ = fr.readline()
                for line in fr:
                    line = line.rstrip("\n\r")
                    if line:
                        w.write(f"{line}\t{sample}\n")

    print(f"[OK] 合并完成：{outp}  共 {len(files)} 个文件")

if __name__ == "__main__":
    main()
