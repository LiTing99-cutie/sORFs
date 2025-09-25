#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, re, csv, sys, gzip, os

def open_maybe(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def load_ids(tsv, col, strip_version=False):
    ids=set()
    with open(tsv, newline="") as f:
        reader=csv.DictReader(f, delimiter="\t")
        # 兼容大小写
        key=None
        for k in reader.fieldnames:
            if k.lower()==col.lower():
                key=k; break
        if key is None:
            sys.exit(f"找不到列: {col}")
        for row in reader:
            v=row[key].strip()
            if not v: 
                continue
            if strip_version:
                v=v.split(".")[0]
            ids.add(v)
    return ids

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--ids_tsv", required=True)
    ap.add_argument("--id_col", default="Gene_id")
    ap.add_argument("--strip_version", action="store_true", help="忽略ID版本号（如 .12）")
    args=ap.parse_args()

    keep=load_ids(args.ids_tsv, args.id_col, args.strip_version)
    pat=re.compile(r'gene_id "([^"]+)"')

    with open_maybe(args.gtf) as fin:
        for line in fin:
            if line.startswith("#"):
                sys.stdout.write(line); continue
            m=pat.search(line)
            if not m:
                continue
            gid=m.group(1)
            if args.strip_version:
                gid=gid.split(".")[0]
            if gid in keep:
                sys.stdout.write(line)

if __name__=="__main__":
    main()
