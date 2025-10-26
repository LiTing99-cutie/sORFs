#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
import csv
from collections import defaultdict, namedtuple

Exon = namedtuple("Exon", "start end")  # 1-based, closed
Feature = namedtuple("Feature", "chrom source ftype start end score strand frame attrs")

def parse_attrs(attr_field: str) -> dict:
    out = {}
    if not attr_field or attr_field == ".":
        return out
    parts = attr_field.strip().strip(";")
    for kv in parts.split(";"):
        kv = kv.strip()
        if not kv:
            continue
        if " " in kv:
            k, v = kv.split(" ", 1)
            out[k] = v.strip().strip('"')
        elif "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip().strip('"')
    return out

def attrs_str(d: dict) -> str:
    keys_order = ["gene_id", "transcript_id", "gene_name", "transcript_name",
                  "source_transcript", "source_gene", "orf_status", "fail_reason", "note"]
    seen = set()
    items = []
    for k in keys_order:
        if k in d:
            items.append(f'{k} "{d[k]}"'); seen.add(k)
    for k, v in d.items():
        if k in seen: continue
        items.append(f'{k} "{v}"')
    return "; ".join(items) + ";"

def read_gtf_build_reference(ref_gtf_path):
    tx = defaultdict(lambda: {"chrom": None, "strand": None, "gene_id": None, "gene_name": None,
                              "exons": [], "cds": [], "start_codons": [], "stop_codons": []})
    gene_to_txs_by_id = defaultdict(list)
    gene_to_txs_by_name = defaultdict(list)
    alt_index = {}
    with open(ref_gtf_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9: continue
            chrom, source, ftype, start, end, score, strand, frame, attr = cols
            start = int(start); end = int(end)
            a = parse_attrs(attr)
            gid = a.get("gene_id"); gname = a.get("gene_name"); tid = a.get("transcript_id")
            if ftype in ("exon", "CDS", "start_codon", "stop_codon") and tid:
                t = tx[tid]
                if t["chrom"] is None:
                    t["chrom"] = chrom; t["strand"] = strand
                    t["gene_id"] = gid; t["gene_name"] = gname
                if ftype == "exon": t["exons"].append(Exon(start, end))
                elif ftype == "CDS": t["cds"].append(Exon(start, end))
                elif ftype == "start_codon": t["start_codons"].append(Exon(start, end))
                elif ftype == "stop_codon":  t["stop_codons"].append(Exon(start, end))
            if ftype == "transcript" and tid:
                if gid: gene_to_txs_by_id[gid].append(tid)
                if gname: gene_to_txs_by_name[gname].append(tid)
            if tid:  alt_index[tid.split(".")[0]] = tid
            if gid:  alt_index[gid.split(".")[0]] = gid
    # 弥补无 transcript 行
    for tid, t in tx.items():
        gid = t["gene_id"]; gname = t["gene_name"]
        if gid and tid not in gene_to_txs_by_id[gid]:
            gene_to_txs_by_id[gid].append(tid)
        if gname and tid not in gene_to_txs_by_name[gname]:
            gene_to_txs_by_name[gname].append(tid)
    for t in tx.values():
        t["exons"].sort(key=lambda e: (e.start, e.end))
        t["cds"].sort(key=lambda e: (e.start, e.end))
        t["start_codons"].sort(key=lambda e: (e.start, e.end))
        t["stop_codons"].sort(key=lambda e: (e.start, e.end))
    return tx, gene_to_txs_by_id, gene_to_txs_by_name, alt_index

def read_pb_exons(pb_gtf_path):
    pb_tx = defaultdict(lambda: {"chrom": None, "strand": None, "gene_id": None, "exons": []})
    with open(pb_gtf_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9: continue
            chrom, source, ftype, start, end, score, strand, frame, attr = cols
            if ftype != "exon": continue
            start = int(start); end = int(end)
            a = parse_attrs(attr)
            tid = a.get("transcript_id"); gid = a.get("gene_id")
            if not tid: continue
            d = pb_tx[tid]
            if d["chrom"] is None:
                d["chrom"] = chrom; d["strand"] = strand; d["gene_id"] = gid
            d["exons"].append(Exon(start, end))
    for d in pb_tx.values():
        d["exons"].sort(key=lambda e: (e.start, e.end))
    return pb_tx

def splice_junctions(exons):
    return {(exons[i].end, exons[i+1].start) for i in range(len(exons)-1)}

def overlap(e: Exon, f: Exon) -> bool:
    return not (e.end < f.start or f.end < e.start)

def cds_splice_junctions(ref_exons, ref_cds_blocks):
    cds_junc = set()
    if not ref_exons or not ref_cds_blocks: return cds_junc
    for i in range(len(ref_exons)-1):
        e1, e2 = ref_exons[i], ref_exons[i+1]
        e1_has = any(overlap(e1, c) for c in ref_cds_blocks)
        e2_has = any(overlap(e2, c) for c in ref_cds_blocks)
        if e1_has and e2_has:
            cds_junc.add((e1.end, e2.start))
    return cds_junc

def cds_contained_by_pb_exons(pb_exons, ref_cds_blocks):
    if not ref_cds_blocks: return False
    j = 0
    for blk in ref_cds_blocks:
        ok = False
        while j < len(pb_exons) and pb_exons[j].end < blk.start:
            j += 1
        k = j
        while k < len(pb_exons) and pb_exons[k].start <= blk.end:
            if pb_exons[k].start <= blk.start and blk.end <= pb_exons[k].end:
                ok = True; break
            k += 1
        if not ok: return False
    return True

def merge_intervals(intervals):
    if not intervals: return []
    intervals = sorted(intervals, key=lambda x: (x.start, x.end))
    out = [intervals[0]]
    for it in intervals[1:]:
        last = out[-1]
        if it.start <= last.end + 1:
            out[-1] = Exon(last.start, max(last.end, it.end))
        else:
            out.append(it)
    return out

def subtract_intervals(interval: Exon, blocks: list):
    if not blocks: return [Exon(interval.start, interval.end)]
    segs = []
    cur = interval.start
    for b in blocks:
        if b.end < interval.start or b.start > interval.end:
            continue
        if cur <= b.start - 1:
            segs.append(Exon(cur, min(interval.end, b.start - 1)))
        cur = max(cur, b.end + 1)
        if cur > interval.end: break
    if cur <= interval.end:
        segs.append(Exon(cur, interval.end))
    return [s for s in segs if s.start <= s.end]

def pick_candidates(struct_cat: str, assoc_tx: str, assoc_gene: str,
                    ref_tx, gene_to_txs_by_id, gene_to_txs_by_name, alt_index):
    struct_cat = (struct_cat or "").strip()
    assoc_tx = (assoc_tx or "").strip()
    assoc_gene = (assoc_gene or "").strip()
    fsm_ism = {"full-splice_match", "incomplete-splice_match"}
    cands = []
    if struct_cat in fsm_ism:
        if assoc_tx:
            if assoc_tx in ref_tx:
                cands = [assoc_tx]
            else:
                nv = assoc_tx.split(".")[0]
                if nv in ref_tx: cands = [nv]
                elif nv in alt_index and alt_index[nv] in ref_tx:
                    cands = [alt_index[nv]]
        return cands
    if assoc_gene:
        if assoc_gene in gene_to_txs_by_id:
            cands = list(gene_to_txs_by_id[assoc_gene])
        else:
            nv = assoc_gene.split(".")[0]
            if nv in gene_to_txs_by_id: cands = list(gene_to_txs_by_id[nv])
        if not cands and assoc_gene in gene_to_txs_by_name:
            cands = list(gene_to_txs_by_name[assoc_gene])
        seen=set(); out=[]
        for t in cands:
            if t not in seen: out.append(t); seen.add(t)
        return out
    return []

def decide_transfer(pb_d, ref_d):
    if ref_d["chrom"] != pb_d["chrom"] or ref_d["strand"] != pb_d["strand"]:
        return False
    pb_junc = splice_junctions(pb_d["exons"])
    ref_cds_junc = cds_splice_junctions(ref_d["exons"], ref_d["cds"])
    if not ref_cds_junc.issubset(pb_junc):
        return False
    if not cds_contained_by_pb_exons(pb_d["exons"], ref_d["cds"]):
        return False
    return True

def derive_codons_from_cds(cds_blocks, strand):
    if not cds_blocks: return [], []
    cds_union = merge_intervals(cds_blocks)
    min_cds = cds_union[0].start; max_cds = cds_union[-1].end
    if strand == "+":
        start = [Exon(min_cds, min_cds + 2)]
        stop  = [Exon(max_cds - 2, max_cds)]
    else:
        start = [Exon(max_cds - 2, max_cds)]
        stop  = [Exon(min_cds, min_cds + 2)]
    return start, stop

def make_utr_generic_on_pb_exons(pb_exons, ref_cds_blocks):
    if not ref_cds_blocks: return []
    cds_union = merge_intervals(ref_cds_blocks)
    utr = []
    for ex in pb_exons:
        remains = subtract_intervals(ex, cds_union)
        utr.extend(remains)
    return merge_intervals(utr)

def write_feature(f: Feature, out):
    attrs = attrs_str(f.attrs)
    out.write("\t".join([f.chrom, f.source, f.ftype, str(f.start), str(f.end),
                         f.score, f.strand, f.frame, attrs]) + "\n")

def main(ref_gtf, pb_gtf, classify_txt):
    ref_tx, gene_to_txs_by_id, gene_to_txs_by_name, alt_index = read_gtf_build_reference(ref_gtf)
    pb_tx = read_pb_exons(pb_gtf)

    import pandas as pd
    df = pd.read_csv(classify_txt, sep="\t", dtype=str, quoting=csv.QUOTE_NONE).fillna("")
    def pick_col(cands):
        for c in df.columns:
            if c.lower() in cands: return c
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
        if not pb_id: continue
        meta[pb_id] = {
            "structural_category": (r.get(col_cat, "") or "").strip(),
            "associated_transcript": (r.get(col_atx, "") or "").strip(),
            "associated_gene": (r.get(col_ag, "") or "").strip()
        }

    # 分组容器：按 associated_gene 聚合；顺便记录 gene span/基本信息
    bucket = defaultdict(list)  # gid -> [rows_for_tx, ...]
    gene_span = {}              # gid -> {chrom,strand,start,end}

    out_path = pb_gtf + ".with_orf.gtf"
    # 先构建每个转录本的 rows
    for pb_id, pb_d in pb_tx.items():
        chrom = pb_d["chrom"]; strand = pb_d["strand"]; exons = pb_d["exons"]
        info = meta.get(pb_id, {"structural_category": "", "associated_transcript": "", "associated_gene": ""})
        assoc_gene = info.get("associated_gene") or pb_d["gene_id"] or "PB_UNKNOWN_GENE"

        # 候选参考
        cands = pick_candidates(info["structural_category"], info["associated_transcript"],
                                assoc_gene, ref_tx, gene_to_txs_by_id, gene_to_txs_by_name, alt_index)

        ref_hit = None
        for ref_tid in cands:
            ref_tid2 = ref_tid
            if ref_tid2 not in ref_tx:
                nv = ref_tid2.split(".")[0]
                if nv in ref_tx: ref_tid2 = nv
                elif nv in alt_index and alt_index[nv] in ref_tx: ref_tid2 = alt_index[nv]
                else: continue
            if ref_tx[ref_tid2]["chrom"] != chrom or ref_tx[ref_tid2]["strand"] != strand:
                continue
            if decide_transfer(pb_d, ref_tx[ref_tid2]):
                ref_hit = ref_tid2
                break

        # 该转录本的输出行（顺序：transcript → 每个外显子块 → UTR）
        rows = []
        tx_start = min(e.start for e in exons); tx_end = max(e.end for e in exons)
        rows.append(Feature(chrom, "PacBio", "transcript", tx_start, tx_end, ".", strand, ".",
                            {"gene_id": assoc_gene, "transcript_id": pb_id}))

        cds_blocks = []; starts = []; stops = []
        if ref_hit:
            rd = ref_tx[ref_hit]
            cds_blocks = list(rd["cds"])
            if not rd["start_codons"] or not rd["stop_codons"]:
                s, p = derive_codons_from_cds(cds_blocks, strand)
                starts = rd["start_codons"] if rd["start_codons"] else s
                stops  = rd["stop_codons"]  if rd["stop_codons"]  else p
            else:
                starts = rd["start_codons"]; stops = rd["stop_codons"]

        # —— 修改点 #1：负链外显子按大→小输出 ——
        exon_iter = exons if strand == "+" else list(reversed(exons))

        for ex in exon_iter:
            rows.append(Feature(chrom, "PacBio", "exon", ex.start, ex.end, ".", strand, ".",
                                {"gene_id": assoc_gene, "transcript_id": pb_id}))
            if ref_hit:
                # 该外显子内的 CDS
                cds_here = [c for c in cds_blocks if ex.start <= c.start and c.end <= ex.end]
                # 负链是否需要倒序？示例里 CDS 在 exon 后、start_codon 随后；保持固定顺序即可
                for c in cds_here:
                    rows.append(Feature(chrom, "PacBio", "CDS", c.start, c.end, ".", strand, ".",
                                        {"gene_id": assoc_gene, "transcript_id": pb_id,
                                         "source_transcript": ref_hit,
                                         "source_gene": ref_tx[ref_hit]["gene_id"],
                                         "orf_status": "transferred"}))
                # start / stop codon（同外显子内；“start 在前、stop 在后”）
                for s in starts:
                    if ex.start <= s.start and s.end <= ex.end:
                        rows.append(Feature(chrom, "PacBio", "start_codon", s.start, s.end, ".", strand, ".",
                                            {"gene_id": assoc_gene, "transcript_id": pb_id}))
                for s in stops:
                    if ex.start <= s.start and s.end <= ex.end:
                        rows.append(Feature(chrom, "PacBio", "stop_codon", s.start, s.end, ".", strand, ".",
                                            {"gene_id": assoc_gene, "transcript_id": pb_id}))

        if ref_hit:
            utrs = make_utr_generic_on_pb_exons(exons, cds_blocks)
            for u in utrs:
                rows.append(Feature(chrom, "PacBio", "UTR", u.start, u.end, ".", strand, ".",
                                    {"gene_id": assoc_gene, "transcript_id": pb_id}))

        bucket[assoc_gene].append(rows)

        # 记录 gene 跨度（用于 gene feature 行）
        if assoc_gene not in gene_span:
            gene_span[assoc_gene] = {"chrom": chrom, "strand": strand, "start": tx_start, "end": tx_end}
        else:
            g = gene_span[assoc_gene]
            if chrom == g["chrom"] and strand == g["strand"]:
                g["start"] = min(g["start"], tx_start)
                g["end"]   = max(g["end"],   tx_end)
            else:
                # 若出现同基因跨染色体/链，保守地扩张坐标并保留最先见到的 chrom/strand
                g["start"] = min(g["start"], tx_start)
                g["end"]   = max(g["end"],   tx_end)

    # 写出：每个基因组先写 gene 行，再写该基因内所有转录本（各自内部已按要求排好）
    out_path = pb_gtf + ".with_orf.gtf"
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("##gff-version 2\n")
        for gid in sorted(bucket.keys()):
            # —— 修改点 #2：输出 gene feature 行 ——
            g = gene_span[gid]
            gene_feat = Feature(g["chrom"], "PacBio", "gene", g["start"], g["end"], ".", g["strand"], ".",
                                {"gene_id": gid})
            write_feature(gene_feat, out)

            # 组内转录本按 (chrom,strand,tx_start,tx_end,transcript_id) 排序
            bucket[gid].sort(key=lambda rows: (rows[0].chrom, rows[0].strand, rows[0].start, rows[0].end,
                                               rows[0].attrs.get("transcript_id","")))
            for rows in bucket[gid]:
                for f in rows:
                    write_feature(f, out)

    sys.stderr.write(f"[DONE] Wrote: {out_path}\n")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Transfer ORF annotations to PB GTF with strand-aware ordering and gene feature.")
    ap.add_argument("--ref_gtf", required=True)
    ap.add_argument("--pb_gtf", required=True)
    ap.add_argument("--classify_txt", required=True)
    args = ap.parse_args()
    main(args.ref_gtf, args.pb_gtf, args.classify_txt)
