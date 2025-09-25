#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ORF注释迁移脚本（坐标判定版）
- 规则：
  1) 参考转录本的“CDS剪切位点集合” ⊆ PB转录本剪切位点集合；
  2) 参考转录本每个CDS片段均被PB的某个外显子完全包含（1-based闭区间）。
- FSM/ISM：仅与 associated_transcript 比较；
- 其他类型：与 associated_gene 下所有参考转录本逐一比较（gene_id 或 gene_name 均可）。
- 迁移成功：写 CDS / start_codon / stop_codon（CDS相位统一“.”），并在PB外显子上切分生成 five_prime_UTR / three_prime_UTR。
- 输出：单一 GTF，按 PB 的 gene_id 分组；未迁移者仅保留 exon。
"""

import sys
import os
import argparse
import csv
from collections import defaultdict, namedtuple

# -----------------------------
# 基础数据结构
# -----------------------------
Exon = namedtuple("Exon", "start end")  # 1-based, closed
Feature = namedtuple("Feature", "chrom source ftype start end score strand frame attrs")  # GTF行

def parse_attrs(attr_field: str) -> dict:
    """解析GTF第9列属性为dict"""
    # 兼容形如 key "val"; key "val2"; 的标准GTF
    out = {}
    if not attr_field or attr_field == ".":
        return out
    parts = attr_field.strip().strip(";")
    # 以 ; 分割，再按第一个空格切 key 和 "value"
    for kv in parts.split(";"):
        kv = kv.strip()
        if not kv:
            continue
        # 允许 key "val" 或 key=val（少见）
        if " " in kv:
            k, v = kv.split(" ", 1)
            v = v.strip().strip('"')
            out[k] = v
        elif "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip().strip('"')
    return out

def attrs_str(d: dict) -> str:
    """dict -> GTF属性串（末尾带分号）"""
    # 按稳定顺序写出常见字段
    keys_order = ["gene_id", "transcript_id", "gene_name", "transcript_name",
                  "source_transcript", "source_gene", "orf_status", "fail_reason", "note"]
    seen = set()
    items = []
    for k in keys_order:
        if k in d:
            items.append(f'{k} "{d[k]}"')
            seen.add(k)
    for k, v in d.items():
        if k in seen:
            continue
        items.append(f'{k} "{v}"')
    return "; ".join(items) + ";"

def read_gtf_build_reference(ref_gtf_path):
    """
    读取参考GTF，构建：
      - tx: {tx_id: {chrom,strand,gene_id,gene_name,exons:[Exon],cds:[Exon],start_codons:[Exon],stop_codons:[Exon]}}
      - gene_to_txs_by_id: {gene_id: [tx_id,...]}
      - gene_to_txs_by_name: {gene_name: [tx_id,...]}
      - alt_index: {tx_id_no_version: tx_id_with_version}  # 同时允许无版本查询
    """
    tx = defaultdict(lambda: {"chrom": None, "strand": None, "gene_id": None, "gene_name": None,
                              "exons": [], "cds": [], "start_codons": [], "stop_codons": []})
    gene_to_txs_by_id = defaultdict(list)
    gene_to_txs_by_name = defaultdict(list)
    alt_index = {}  # 无版本 -> 有版本（若存在多个版本，后写覆盖，但我们也会直接建索引到自身）

    with open(ref_gtf_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, frame, attr = cols
            start = int(start); end = int(end)
            a = parse_attrs(attr)
            gid = a.get("gene_id")
            gname = a.get("gene_name")
            tid = a.get("transcript_id")
            if ftype in ("exon", "CDS", "start_codon", "stop_codon") and tid:
                t = tx[tid]
                if t["chrom"] is None:
                    t["chrom"] = chrom
                    t["strand"] = strand
                    t["gene_id"] = gid
                    t["gene_name"] = gname
                # 一致性检查略过
                if ftype == "exon":
                    t["exons"].append(Exon(start, end))
                elif ftype == "CDS":
                    t["cds"].append(Exon(start, end))
                elif ftype == "start_codon":
                    t["start_codons"].append(Exon(start, end))
                elif ftype == "stop_codon":
                    t["stop_codons"].append(Exon(start, end))
            # 索引 gene->tx
            if ftype == "transcript" and tid:
                if gid:
                    gene_to_txs_by_id[gid].append(tid)
                if gname:
                    gene_to_txs_by_name[gname].append(tid)
            # 建立无版本转有版本索引（ENST/ENSG常见）
            if tid:
                tid_nv = tid.split(".")[0]
                alt_index[tid_nv] = tid
            if gid:
                gid_nv = gid.split(".")[0]
                alt_index[gid_nv] = gid

    # 一些 GTF 里没有 transcript 行，确保 gene_to_txs* 完整（从已收集的 exon/CDS 也可推）
    for tid, t in tx.items():
        gid = t["gene_id"]; gname = t["gene_name"]
        if gid and tid not in gene_to_txs_by_id[gid]:
            gene_to_txs_by_id[gid].append(tid)
        if gname and tid not in gene_to_txs_by_name[gname]:
            gene_to_txs_by_name[gname].append(tid)

    # 排序 exon/CDS 片段（按基因组坐标升序）
    for tid, t in tx.items():
        t["exons"].sort(key=lambda e: (e.start, e.end))
        t["cds"].sort(key=lambda e: (e.start, e.end))
        t["start_codons"].sort(key=lambda e: (e.start, e.end))
        t["stop_codons"].sort(key=lambda e: (e.start, e.end))

    return tx, gene_to_txs_by_id, gene_to_txs_by_name, alt_index

def read_pb_exons(pb_gtf_path):
    """
    读取PB的GTF，仅收集 exon：
      pb_tx: {transcript_id: {chrom,strand,gene_id,exons:[Exon]}}
      genes: {gene_id: set(transcript_id)}
    """
    pb_tx = defaultdict(lambda: {"chrom": None, "strand": None, "gene_id": None, "exons": []})
    genes = defaultdict(set)
    with open(pb_gtf_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, frame, attr = cols
            if ftype != "exon":
                continue
            start = int(start); end = int(end)
            a = parse_attrs(attr)
            tid = a.get("transcript_id")
            gid = a.get("gene_id")
            if not tid or not gid:
                # 跳过缺失ID的行
                continue
            d = pb_tx[tid]
            if d["chrom"] is None:
                d["chrom"] = chrom
                d["strand"] = strand
                d["gene_id"] = gid
            d["exons"].append(Exon(start, end))
            genes[gid].add(tid)
    # 排序外显子
    for tid, d in pb_tx.items():
        d["exons"].sort(key=lambda e: (e.start, e.end))
    return pb_tx, genes

def splice_junctions(exons):
    """构造剪切位点集合：{(exon_i.end, exon_{i+1}.start)}，1-based闭区间，链无关"""
    return {(exons[i].end, exons[i+1].start) for i in range(len(exons)-1)}

def overlap(e: Exon, f: Exon) -> bool:
    return not (e.end < f.start or f.end < e.start)

def cds_splice_junctions(ref_exons, ref_cds_blocks):
    """
    参考转录本的“CDS剪切位点集合”：两侧外显子都与某段CDS有交叠的内含子边界
    """
    cds_junc = set()
    if not ref_exons or not ref_cds_blocks:
        return cds_junc
    for i in range(len(ref_exons)-1):
        e1, e2 = ref_exons[i], ref_exons[i+1]
        e1_has = any(overlap(e1, c) for c in ref_cds_blocks)
        e2_has = any(overlap(e2, c) for c in ref_cds_blocks)
        if e1_has and e2_has:
            cds_junc.add((e1.end, e2.start))
    return cds_junc

def cds_contained_by_pb_exons(pb_exons, ref_cds_blocks):
    """
    参考每段CDS需被PB某个外显子完全包含（1-based闭区间）
    """
    if not ref_cds_blocks:
        return False
    j = 0
    for blk in ref_cds_blocks:
        ok = False
        while j < len(pb_exons) and pb_exons[j].end < blk.start:
            j += 1
        k = j
        while k < len(pb_exons) and pb_exons[k].start <= blk.end:
            S, E = pb_exons[k].start, pb_exons[k].end
            if S <= blk.start and blk.end <= E:
                ok = True
                break
            k += 1
        if not ok:
            return False
    return True

def subtract_intervals(interval: Exon, blocks: list):
    """
    从一个区间 interval 中减去 blocks 的并集，返回剩余区间列表（不相交，已排序）
    interval, blocks 为 1-based 闭区间
    """
    if not blocks:
        return [Exon(interval.start, interval.end)]
    # 先裁剪出与 interval 有交叠的CDS
    segs = []
    cur_start = interval.start
    for b in blocks:
        if b.end < interval.start or b.start > interval.end:
            continue
        # 左侧剩余
        if cur_start <= b.start - 1 and cur_start <= interval.end and b.start - 1 >= interval.start:
            segs.append(Exon(cur_start, min(interval.end, b.start - 1)))
        cur_start = max(cur_start, b.end + 1)
        if cur_start > interval.end:
            break
    if cur_start <= interval.end:
        segs.append(Exon(cur_start, interval.end))
    return [s for s in segs if s.start <= s.end]

def merge_intervals(intervals):
    """合并重叠/相邻区间，返回升序列表"""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x.start, x.end))
    out = [intervals[0]]
    for it in intervals[1:]:
        last = out[-1]
        if it.start <= last.end + 1:
            out[-1] = Exon(last.start, max(last.end, it.end))
        else:
            out.append(it)
    return out

def pick_candidates(struct_cat: str, assoc_tx: str, assoc_gene: str,
                    ref_tx, gene_to_txs_by_id, gene_to_txs_by_name, alt_index):
    """
    返回候选参考转录本ID列表（按优先顺序）。
    - FSM/ISM：仅 associated_transcript（允许无版本/有版本匹配）
    - 其他：associated_gene 作为 gene_id 或 gene_name 在参考中展开
    """
    struct_cat = (struct_cat or "").strip()
    assoc_tx = (assoc_tx or "").strip()
    assoc_gene = (assoc_gene or "").strip()
    fsm_ism = {"full-splice_match", "incomplete-splice_match"}

    cands = []

    if struct_cat in fsm_ism:
        if assoc_tx:
            # 支持去版本匹配
            if assoc_tx in ref_tx:
                cands = [assoc_tx]
            else:
                assoc_tx_nv = assoc_tx.split(".")[0]
                if assoc_tx_nv in ref_tx:
                    cands = [assoc_tx_nv]
                elif assoc_tx_nv in alt_index:
                    real = alt_index[assoc_tx_nv]
                    if real in ref_tx:
                        cands = [real]
        return cands

    # 其他类型：用 gene 展开（允许 gene_id 或 gene_name）
    if assoc_gene:
        # 先当 gene_id 找
        if assoc_gene in gene_to_txs_by_id:
            cands = list(gene_to_txs_by_id[assoc_gene])
        else:
            # 尝试去版本
            gn_nv = assoc_gene.split(".")[0]
            if gn_nv in gene_to_txs_by_id:
                cands = list(gene_to_txs_by_id[gn_nv])
        # 再尝试 gene_name
        if not cands:
            if assoc_gene in gene_to_txs_by_name:
                cands = list(gene_to_txs_by_name[assoc_gene])
        # 去重保持顺序
        seen = set(); uniq = []
        for t in cands:
            if t not in seen:
                uniq.append(t); seen.add(t)
        return uniq

    return []

def decide_transfer(pb_d, ref_d):
    """
    判定是否迁移：CDS剪切位点⊆PB剪切位点 且 CDS段被包含
    前置：同染色体、同链
    """
    if ref_d["chrom"] != pb_d["chrom"] or ref_d["strand"] != pb_d["strand"]:
        return False
    pb_junc = splice_junctions(pb_d["exons"])
    ref_cds_junc = cds_splice_junctions(ref_d["exons"], ref_d["cds"])
    if not ref_cds_junc.issubset(pb_junc):
        return False
    if not cds_contained_by_pb_exons(pb_d["exons"], ref_d["cds"]):
        return False
    return True

def make_gene_tx_header(pb_gene_id, pb_tx_id, chrom, strand, exons):
    """构造 gene / transcript 头行（score/phase 均为 .）"""
    tx_start = min(e.start for e in exons)
    tx_end = max(e.end for e in exons)
    gene_feat = Feature(chrom, "PacBio", "gene", tx_start, tx_end, ".", strand, ".",
                        {"gene_id": pb_gene_id})
    tx_feat = Feature(chrom, "PacBio", "transcript", tx_start, tx_end, ".", strand, ".",
                      {"gene_id": pb_gene_id, "transcript_id": pb_tx_id})
    return gene_feat, tx_feat

def derive_codons_from_cds(cds_blocks, strand):
    """若参考没有codon，按CDS边界推导start/stop_codon"""
    if not cds_blocks:
        return [], []
    cds_union = merge_intervals(cds_blocks)
    min_cds = cds_union[0].start
    max_cds = cds_union[-1].end
    if strand == "+":
        start = [Exon(min_cds, min_cds + 2)]
        stop = [Exon(max_cds - 2, max_cds)]
    else:
        start = [Exon(max_cds - 2, max_cds)]
        stop = [Exon(min_cds, min_cds + 2)]
    return start, stop

def make_utr_on_pb_exons(pb_exons, ref_cds_blocks, strand):
    """
    在 PB 的 exon 上生成 five_prime_UTR / three_prime_UTR
    - 简化：以全局 CDS 范围 [min_cds, max_cds] 为界，外侧为 UTR
    """
    if not ref_cds_blocks:
        return [], []
    cds_union = merge_intervals(ref_cds_blocks)
    min_cds = cds_union[0].start
    max_cds = cds_union[-1].end

    five, three = [], []
    for ex in pb_exons:
        # 先从exon中减去CDS联合，得到非编码片段
        noncoding = subtract_intervals(ex, cds_union)
        for seg in noncoding:
            if strand == "+":
                # 位于全局CDS左侧 → 5'；右侧 → 3'；跨界段拆分
                if seg.end < min_cds:
                    five.append(seg)
                elif seg.start > max_cds:
                    three.append(seg)
                else:
                    # 跨界：左半归5'，右半归3'
                    if seg.start < min_cds:
                        five.append(Exon(seg.start, min_cds - 1))
                    if seg.end > max_cds:
                        three.append(Exon(max(min_cds + 1, max_cds + 1), seg.end))
            else:
                # 负链相反：基因组右侧是5'，左侧是3'
                if seg.start > max_cds:
                    five.append(seg)
                elif seg.end < min_cds:
                    three.append(seg)
                else:
                    if seg.end > max_cds:
                        five.append(Exon(max_cds + 1, seg.end))
                    if seg.start < min_cds:
                        three.append(Exon(seg.start, min_cds - 1))
    # 合并可能相邻片段
    return merge_intervals(five), merge_intervals(three)

def write_feature(f: Feature, out):
    """写一条 GTF 记录"""
    attrs = attrs_str(f.attrs)
    row = [f.chrom, f.source, f.ftype, str(f.start), str(f.end), f.score, f.strand, f.frame, attrs]
    out.write("\t".join(row) + "\n")

def main(ref_gtf, pb_gtf, classify_txt):
    # 读取参考与PB
    ref_tx, gene_to_txs_by_id, gene_to_txs_by_name, alt_index = read_gtf_build_reference(ref_gtf)
    pb_tx, pb_genes = read_pb_exons(pb_gtf)

    # 读取分类表
    # 尝试找到列名：isoform / structural_category / associated_transcript / associated_gene
    import pandas as pd
    df = pd.read_csv(classify_txt, sep="\t", dtype=str, quoting=csv.QUOTE_NONE).fillna("")
    # 兼容多种列名
    def pick_col(cands):
        for c in df.columns:
            if c.lower() in cands:
                return c
        return None
    col_iso = pick_col({"isoform", "pb", "pb_id", "transcript_id", "isoform id"})
    col_cat = pick_col({"structural_category"})
    col_atx = pick_col({"associated_transcript"})
    col_ag  = pick_col({"associated_gene", "gene", "gene_name", "gene_id"})
    if not col_iso or not col_cat:
        sys.stderr.write("ERROR: 无法识别 classify_txt 的关键列（isoform/structural_category）。\n")
        sys.exit(1)

    meta = {}
    for _, r in df.iterrows():
        pb_id = (r[col_iso] or "").strip()
        if not pb_id:
            continue
        meta[pb_id] = {
            "structural_category": (r.get(col_cat, "") or "").strip(),
            "associated_transcript": (r.get(col_atx, "") or "").strip(),
            "associated_gene": (r.get(col_ag, "") or "").strip()
        }

    # 处理每个 PB transcript
    # 收集输出：按 gene_id 分组
    gene_bucket = defaultdict(list)  # gene_id -> [Feature,...]

    for pb_id, pb_d in pb_tx.items():
        pb_gene = pb_d["gene_id"]; chrom = pb_d["chrom"]; strand = pb_d["strand"]
        exons = pb_d["exons"]

        # gene / transcript 头行（保持统一输出风格）
        gene_feat, tx_feat = make_gene_tx_header(pb_gene, pb_id, chrom, strand, exons)
        gene_bucket[pb_gene].append(gene_feat)
        gene_bucket[pb_gene].append(tx_feat)

        # 原样输出 exon
        for ex in exons:
            gene_bucket[pb_gene].append(
                Feature(chrom, "PacBio", "exon", ex.start, ex.end, ".", strand, ".",
                        {"gene_id": pb_gene, "transcript_id": pb_id})
            )

        # 找分类信息
        info = meta.get(pb_id, {"structural_category": "", "associated_transcript": "", "associated_gene": ""})
        cands = pick_candidates(info["structural_category"], info["associated_transcript"],
                                info["associated_gene"], ref_tx, gene_to_txs_by_id, gene_to_txs_by_name, alt_index)

        transferred = False
        fail_reason = "no_candidate" if not cands else "no_match"

        for ref_tid in cands:
            # 允许无版本->有版本解析
            ref_tid2 = ref_tid
            if ref_tid2 not in ref_tx:
                nv = ref_tid2.split(".")[0]
                if nv in ref_tx:
                    ref_tid2 = nv
                elif nv in alt_index and alt_index[nv] in ref_tx:
                    ref_tid2 = alt_index[nv]
                else:
                    continue
            ref_d = ref_tx[ref_tid2]
            # 染色体/链必须一致
            if ref_d["chrom"] != chrom or ref_d["strand"] != strand:
                continue
            if decide_transfer(pb_d, ref_d):
                # 迁移 CDS
                for c in ref_d["cds"]:
                    gene_bucket[pb_gene].append(
                        Feature(chrom, "PacBio", "CDS", c.start, c.end, ".", strand, ".",
                                {"gene_id": pb_gene, "transcript_id": pb_id,
                                 "source_transcript": ref_tid2, "source_gene": ref_d["gene_id"], "orf_status": "transferred"})
                    )
                # 迁移 / 推导 codons
                starts = ref_d["start_codons"] if ref_d["start_codons"] else None
                stops  = ref_d["stop_codons"]  if ref_d["stop_codons"]  else None
                if not starts or not stops:
                    ds, dp = derive_codons_from_cds(ref_d["cds"], strand)
                    if not starts: starts = ds
                    if not stops:  stops  = dp
                for s in starts:
                    gene_bucket[pb_gene].append(
                        Feature(chrom, "PacBio", "start_codon", s.start, s.end, ".", strand, ".",
                                {"gene_id": pb_gene, "transcript_id": pb_id})
                    )
                for s in stops:
                    gene_bucket[pb_gene].append(
                        Feature(chrom, "PacBio", "stop_codon", s.start, s.end, ".", strand, ".",
                                {"gene_id": pb_gene, "transcript_id": pb_id})
                    )
                # UTR 生成（基于参考CDS联合在PB外显子上切分）
                five, three = make_utr_on_pb_exons(exons, ref_d["cds"], strand)
                for u in five:
                    gene_bucket[pb_gene].append(
                        Feature(chrom, "PacBio", "five_prime_UTR", u.start, u.end, ".", strand, ".",
                                {"gene_id": pb_gene, "transcript_id": pb_id})
                    )
                for u in three:
                    gene_bucket[pb_gene].append(
                        Feature(chrom, "PacBio", "three_prime_UTR", u.start, u.end, ".", strand, ".",
                                {"gene_id": pb_gene, "transcript_id": pb_id})
                    )
                transferred = True
                break

        if not transferred:
            # 未迁移：可在 transcript 行上标记 orf_status=none（也可不标）
            # 这里不再追加 fail_reason 行，保持输出干净
            pass

    # 写出：按 gene_id 分组，组内按 feature/坐标排序
    out_path = pb_gtf + ".with_orf.gtf"
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("##gff-version 2\n")
        # 对每个gene分组输出
        for gid in sorted(gene_bucket.keys()):
            feats = gene_bucket[gid]
            # 排序：gene(0) < transcript(1) < exon(2) < five_UTR(3) < CDS(4) < three_UTR(5) < start_codon(6) < stop_codon(7)
            order = {"gene":0, "transcript":1, "exon":2, "five_prime_UTR":3, "CDS":4,
                     "three_prime_UTR":5, "start_codon":6, "stop_codon":7}
            feats.sort(key=lambda f: (f.chrom, f.strand, order.get(f.ftype, 99), f.start, f.end, f.attrs.get("transcript_id","")))
            for f in feats:
                write_feature(f, out)

    sys.stderr.write(f"[DONE] Wrote: {out_path}\n")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Transfer ORF annotations (CDS/start/stop + UTR) from reference GTF to PB GTF.")
    ap.add_argument("--ref_gtf", required=True, help="参考注释 GTF（如 gencode.v41.annotation.gtf）")
    ap.add_argument("--pb_gtf", required=True, help="PB 的 GTF（包含 exon 行）")
    ap.add_argument("--classify_txt", required=True, help="Iso-Seq 分类表 TSV（含 isoform/structural_category/associated_* 列）")
    args = ap.parse_args()

    main(args.ref_gtf, args.pb_gtf, args.classify_txt)
