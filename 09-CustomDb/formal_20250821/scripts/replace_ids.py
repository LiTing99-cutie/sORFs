#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, argparse, gzip, csv

def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8", errors="ignore")

def load_mapping(tsv_path, key_col="first_id_in_line", val_col="chosen_id"):
    mp = {}
    with open_maybe_gzip(tsv_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if key_col not in reader.fieldnames or val_col not in reader.fieldnames:
            sys.exit(f"[ERROR] {tsv_path} 缺少列：{key_col} 或 {val_col}")
        for row in reader:
            k = (row.get(key_col) or "").strip()
            v = (row.get(val_col) or "").strip()
            if not k or not v:
                continue
            mp[k] = v
    return mp

def rewrite_fasta_ids(fa_in, fa_out, mapping):
    total = 0
    changed = 0
    seen_ids = set()

    with open_maybe_gzip(fa_in, "rt") as fi, open(fa_out, "w") as fo:
        write_seq = False
        for line in fi:
            if line.startswith(">"):
                total += 1
                # 解析原ID与描述
                body = line[1:].rstrip("\n")
                if not body:
                    fo.write(line)  # 空header，原样写出
                    continue
                parts = body.split(None, 1)  # 只拆前两个：ID + 余下描述
                old_id = parts[0]
                desc = (" " + parts[1]) if len(parts) > 1 else ""

                new_id = mapping.get(old_id, old_id)
                if new_id != old_id:
                    changed += 1
                # 可选：检测重名
                if new_id in seen_ids:
                    # 若担心ID碰撞，可在这里报错或添加后缀；当前仅告警
                    pass
                seen_ids.add(new_id)

                fo.write(">" + new_id + desc + "\n")
            else:
                fo.write(line)

    sys.stderr.write(f"[DONE] 处理FASTA记录 {total} 条，替换 {changed} 条。\n")

def main():
    ap = argparse.ArgumentParser(
        description="将FASTA中first_id_in_line替换为chosen_id（只改header首字段）"
    )
    ap.add_argument("--map", required=True, help="分组映射文件（含列 first_id_in_line 与 chosen_id，TSV，可.gz）")
    ap.add_argument("--fasta-in", required=True, help="待替换的FASTA（可.gz）")
    ap.add_argument("--fasta-out", required=True, help="输出FASTA路径")
    ap.add_argument("--key-col", default="first_id_in_line", help="映射表键列名（默认 first_id_in_line）")
    ap.add_argument("--val-col", default="chosen_id", help="映射表值列名（默认 chosen_id）")
    args = ap.parse_args()

    mapping = load_mapping(args.map, key_col=args.key_col, val_col=args.val_col)
    if not mapping:
        sys.stderr.write("[WARN] 映射表为空，原样复制。\n")
    rewrite_fasta_ids(args.fasta_in, args.fasta_out, mapping)

if __name__ == "__main__":
    main()
