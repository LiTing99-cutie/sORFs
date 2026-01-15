#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
严格版 ORF 同义 / 非同义位点计数

输入：
  1) variants_in_orfs.tsv  : run_pnps_strict.sh 生成的交集结果
     格式（tab）：
       1: chr
       2: start (0-based)
       3: end   (1-based, 等于 VCF POS)
       4: REF
       5: ALT
       6: INFO  （包含 ANN=...）
       7-: ORF bed 列（其中第10列为 ORF_ID）

  2) orf_bed : ORF 区间 bed 文件（第 4 列为 ORF_ID）

输出：
  3) out_tsv : 每行一个 ORF_ID，对应该 ORF 的
               同义位点数 syn_sites、非同义位点数 nonsyn_sites
"""

import sys
from collections import defaultdict


def main():
    if len(sys.argv) != 4:
        sys.stderr.write(
            "Usage: python count_orf_syn_nonsyn_strict.py "
            "<variants_in_orfs.tsv> <orf_bed> <out_tsv>\n"
        )
        sys.exit(1)

    intersect_tsv = sys.argv[1]
    orf_bed = sys.argv[2]
    out_tsv = sys.argv[3]

    # ===== 1. 读取 ORF 列表（bed 第4列）=====
    orf_ids = []
    with open(orf_bed) as bed_f:
        for line in bed_f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            orf_id = cols[3]
            orf_ids.append(orf_id)

    # 去重并保留原始顺序
    orf_ids = list(dict.fromkeys(orf_ids))

    # ===== 2. 遍历交集文件，按 ORF 解析 ANN，统计同义 / 非同义位点 =====
    syn_count = defaultdict(int)
    nonsyn_count = defaultdict(int)

    # 防止同一变异在同一 ORF 上重复计数
    # key = (chr, pos, ref, alt, orf_id)
    seen = set()

    def extract_ann(info: str):
        """从 INFO 字段中提取 ANN= 后面的内容"""
        for part in info.split(";"):
            if part.startswith("ANN="):
                return part[4:]
        return None

    with open(intersect_tsv) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            # 1:chr 2:start 3:end 4:REF 5:ALT 6:INFO 7-: ORF bed 列
            if len(cols) < 10:
                continue

            chrom = cols[0]
            end_pos_1based = cols[2]   # bed end = VCF POS
            ref = cols[3]
            alt_field = cols[4]
            info = cols[5]
            orf_id = cols[9]           # ORF bed 第4列 → 这里第10列

            key = (chrom, end_pos_1based, ref, alt_field, orf_id)
            if key in seen:
                continue
            seen.add(key)

            ann_str = extract_ann(info)
            if not ann_str:
                continue

            alts = alt_field.split(",")
            var_syn = False
            var_nonsyn = False

            # 我们匹配的模式：在 ANN 的一条记录里，包含 "|<ORF_ID>|"
            # 例如：|PB.17.2:chr1:+|64|3647:967:1336|noncoding|ATG|
            pattern = "|" + orf_id + "|"

            for ann_item in ann_str.split(","):
                if not ann_item:
                    continue
                fields = ann_item.split("|")
                if len(fields) < 2:
                    continue

                allele = fields[0]
                effect = fields[1]

                # 1）只看当前 ALT 等位基因
                if allele not in alts:
                    continue

                # 2）只看这一条 ORF 的注释 —— 通过包含 ORF_ID 子串来判断
                if pattern not in ann_item:
                    # 极端情况：ORF_ID 恰好在末尾，没有尾部的 "|"
                    if not ann_item.endswith("|" + orf_id):
                        continue

                # 3）判断效应类型
                if "synonymous_variant" in effect:
                    var_syn = True
                elif "missense_variant" in effect:
                    var_nonsyn = True
                # 如需把 stop_gained / stop_lost 也算非同义，可加：
                # elif "stop_gained" in effect or "stop_lost" in effect:
                #     var_nonsyn = True

            # 一个变异在这个 ORF 上只作一次判定：同义优先，其次非同义
            if var_syn:
                syn_count[orf_id] += 1
            elif var_nonsyn:
                nonsyn_count[orf_id] += 1

    # ===== 3. 输出所有 ORF（没有变异的补 0）=====
    with open(out_tsv, "w") as out_f:
        out_f.write("ORF_ID\tsyn_sites\tnonsyn_sites\n")
        for oid in orf_ids:
            s = syn_count.get(oid, 0)
            n = nonsyn_count.get(oid, 0)
            out_f.write(f"{oid}\t{s}\t{n}\n")

    sys.stderr.write(f"Done. Written to {out_tsv}\n")


if __name__ == "__main__":
    main()
