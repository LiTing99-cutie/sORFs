#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filter ORF FASTA by gene-level cytoplasmic expression (C) and ORF-level Ribo evidence.

保留 ORF 的条件：
  (1) 基因层：C >= c_thresh   （来自 rpkm_N_C_A.txt 的 Geneid, C）
AND
  (2) ORF层：按照所选列（RPF侧、Psites侧）与阈值判断，通过 evidence_mode 组合：
      - any      : (RPF >= rpf_min) OR (Psites >= psite_min)
      - both     : (RPF >= rpf_min) AND (Psites >= psite_min)
      - rpf_only : (RPF >= rpf_min)
      - ps_only  : (Psites >= psite_min)

输入文件：
- rpkm_N_C_A.txt              (列：Geneid, C)
- isoform.expr.info.txt       (列：isoform, associated_gene)
- orf.rpf.psite.txt           (列：ORF_id, 可选的 RPF/Psites 相关列)
- FASTA（头部为完整 ORF_id；isoform = 头部第一个冒号前）

RPF/Psites 可用列名（示例）：
- RPF:    RPF_reads | RPF_RPKM | RPF_codon_coverage
- Psites: Psites_number | Psites_RPKM | Psites_codon_coverage

聚合策略：
- Gene C：同一 Geneid 取最大值（避免误删）
- ORF 证据：同一 ORF_id 对选定列求和（整合多个来源/重复）
"""

import sys, argparse, gzip, os

# ---------- I/O helpers ----------

def smart_open(path, mode='rt'):
    """Auto-detect gzip by suffix."""
    if path.endswith('.gz'):
        return gzip.open(path, mode=mode, encoding=None if 'b' in mode else 'utf-8', newline='')
    return open(path, mode=mode, encoding=None if 'b' in mode else 'utf-8', newline='')

def ensure_file(path, name):
    if not path or not os.path.exists(path):
        sys.exit(f"ERROR: {name} not found: {path}")

# ---------- Readers ----------

def read_gene_C(path):
    """Read Geneid -> max(C)."""
    gene2C = {}
    with smart_open(path, 'rt') as f:
        header = f.readline().rstrip('\n').split('\t')
        try:
            gi = header.index('Geneid')
            ci = header.index('C')
        except ValueError:
            sys.exit("ERROR: rpkm_N_C_A.txt 必须包含列 'Geneid' 与 'C'.")
        for ln, line in enumerate(f, start=2):
            if not line.strip(): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= max(gi, ci): continue
            gene = parts[gi].strip()
            cstr = parts[ci].strip()
            if not gene or cstr == '': continue
            try:
                cval = float(cstr)
            except:
                continue
            if gene not in gene2C or cval > gene2C[gene]:
                gene2C[gene] = cval
    return gene2C

def read_isoform_map(path):
    """Read isoform -> associated_gene (first occurrence)."""
    iso2gene = {}
    with smart_open(path, 'rt') as f:
        header = f.readline().rstrip('\n').split('\t')
        try:
            ii = header.index('isoform')
            gi = header.index('associated_gene')
        except ValueError:
            sys.exit("ERROR: isoform.expr.info.txt 必须包含列 'isoform' 与 'associated_gene'.")
        for line in f:
            if not line.strip(): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= max(ii, gi): continue
            iso = parts[ii].strip()
            gene = parts[gi].strip()
            if iso and gene and iso not in iso2gene:
                iso2gene[iso] = gene
    return iso2gene

def read_orf_ribo_dynamic(path, rpf_col, ps_col):
    """
    Read chosen metrics for ORF: ORF_id -> (sum_rpf, sum_ps)
    Values as float to accommodate RPKM/coverage.
    """
    orf2ribo = {}
    with smart_open(path, 'rt') as f:
        header = f.readline().rstrip('\n').split('\t')
        # normalize header strip
        header = [h.strip() for h in header]
        try:
            oi = header.index('ORF_id')
            ri = header.index(rpf_col)
            pi = header.index(ps_col)
        except ValueError as e:
            cols = ','.join(header)
            sys.exit(f"ERROR: 找不到所需列: {e}. 当前表头: {cols}")
        for ln, line in enumerate(f, start=2):
            if not line.strip(): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= max(oi, ri, pi): continue
            oid = parts[oi].strip()
            if not oid: continue
            def to_float(x):
                try:
                    return float(x.strip())
                except:
                    return 0.0
            r = to_float(parts[ri]) if parts[ri] != '' else 0.0
            p = to_float(parts[pi]) if parts[pi] != '' else 0.0
            if oid in orf2ribo:
                r0, p0 = orf2ribo[oid]
                orf2ribo[oid] = (r0 + r, p0 + p)
            else:
                orf2ribo[oid] = (r, p)
    return orf2ribo

# ---------- FASTA stream ----------

def fasta_stream(path):
    """Yield (header_without_gt, seq) streaming. Header keeps original (without leading '>')."""
    with smart_open(path, 'rt') as f:
        header = None
        seq_chunks = []
        for line in f:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_chunks)
                header = line[1:].strip()  # remove leading '>'
                seq_chunks = []
            else:
                s = line.strip()
                if s:
                    seq_chunks.append(s)
        if header is not None:
            yield header, ''.join(seq_chunks)

def wrap_seq(seq, width=60):
    for i in range(0, len(seq), width):
        yield seq[i:i+width]

# ---------- Logic ----------

def orf_passes(rpf_val, ps_val, rpf_min, ps_min, mode):
    """Return True/False depending on evidence_mode."""
    rp_ok = (rpf_val is not None) and (rpf_val >= rpf_min)
    ps_ok = (ps_val  is not None) and (ps_val  >= ps_min)
    if mode == 'any':
        return rp_ok or ps_ok
    elif mode == 'both':
        return rp_ok and ps_ok
    elif mode == 'rpf_only':
        return rp_ok
    elif mode == 'ps_only':
        return ps_ok
    else:
        raise ValueError(f"Unknown evidence_mode: {mode}")

# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(description="Filter ORFs in FASTA by gene C and ORF Ribo evidence.")
    ap.add_argument('--rpkm', required=True, help="rpkm_N_C_A.txt (cols: Geneid, C)")
    ap.add_argument('--isoform_map', required=True, help="isoform.expr.info.txt (cols: isoform, associated_gene)")
    ap.add_argument('--ribo', required=True, help="orf.rpf.psite.txt (cols include ORF_id and chosen metrics)")
    ap.add_argument('--fasta', required=True, help="input FASTA (.fa/.fa.gz). Header is full ORF_id; isoform = token before first ':'")
    ap.add_argument('--out_fasta', required=True, help="output FASTA after filtering (.fa or .fa.gz)")
    ap.add_argument('--trace_tsv', default=None, help="optional TSV trace for each ORF (keep/drop reason)")

    # thresholds & metrics
    ap.add_argument('--c_thresh', type=float, default=0.2, help="threshold for C (>=)")
    ap.add_argument('--rpf_metric', default='RPF_reads',
                    help="RPF metric column: e.g., RPF_reads|RPF_RPKM|RPF_codon_coverage")
    ap.add_argument('--ps_metric', default='Psites_number',
                    help="P-sites metric column: e.g., Psites_number|Psites_RPKM|Psites_codon_coverage")
    ap.add_argument('--rpf_min', type=float, default=10.0, help="RPF metric threshold (>=)")
    ap.add_argument('--psite_min', type=float, default=2.0, help="P-sites metric threshold (>=)")
    ap.add_argument('--evidence_mode', choices=['any','both','rpf_only','ps_only'], default='any',
                    help="Combine RPF/P-sites: any|both|rpf_only|ps_only (default any)")

    args = ap.parse_args()

    for p, n in [(args.rpkm,'rpkm_N_C_A.txt'), (args.isoform_map,'isoform.expr.info.txt'),
                 (args.ribo,'orf.rpf.psite.txt'), (args.fasta,'FASTA')]:
        ensure_file(p, n)

    # Read tables
    gene2C  = read_gene_C(args.rpkm)
    iso2gene= read_isoform_map(args.isoform_map)
    orf2ribo= read_orf_ribo_dynamic(args.ribo, args.rpf_metric, args.ps_metric)

    n_before = 0
    n_after  = 0

    # trace
    trace_lines = []
    if args.trace_tsv:
        trace_lines.append('\t'.join([
            'ORF_id','isoform','associated_gene','C',
            f"{args.rpf_metric}", f"{args.ps_metric}",
            'keep','drop_reason'
        ]))

    with smart_open(args.out_fasta, 'wt') as fw:
        for orf_id, seq in fasta_stream(args.fasta):
            n_before += 1
            isoform = orf_id.split(':', 1)[0] if ':' in orf_id else orf_id

            drop_reason = []
            gene = iso2gene.get(isoform)
            if not gene:
                drop_reason.append('no_isoform_gene')

            Cval = gene2C.get(gene) if gene else None
            gene_pass = (Cval is not None and Cval >= args.c_thresh)
            if not gene_pass:
                if Cval is None:
                    drop_reason.append('missing_C')
                else:
                    drop_reason.append('low_C')

            rpf_val, ps_val = orf2ribo.get(orf_id, (None, None))
            if rpf_val is None and ps_val is None:
                drop_reason.append('no_ribo_evidence')
                rpf_use, ps_use = 0.0, 0.0
                orf_evidence = False
            else:
                rpf_use = 0.0 if rpf_val is None else rpf_val
                ps_use  = 0.0 if ps_val  is None else ps_val
                orf_evidence = orf_passes(rpf_use, ps_use, args.rpf_min, args.psite_min, args.evidence_mode)
                if not orf_evidence:
                    drop_reason.append('weak_ribo_evidence')

            keep = gene_pass and orf_evidence

            if keep:
                fw.write('>' + orf_id + '\n')
                for chunk in wrap_seq(seq, 60):
                    fw.write(chunk + '\n')
                n_after += 1

            if args.trace_tsv:
                trace_lines.append('\t'.join([
                    orf_id,
                    isoform,
                    gene if gene else '',
                    (f"{Cval:.6g}" if Cval is not None else ''),
                    (f"{rpf_use:.6g}" if rpf_use is not None else ''),
                    (f"{ps_use:.6g}"  if ps_use  is not None else ''),
                    '1' if keep else '0',
                    ';'.join(drop_reason) if drop_reason else ''
                ]))

    kept_ratio = (n_after / n_before) if n_before > 0 else 0.0
    sys.stderr.write(f"[Summary] before={n_before} after={n_after} removed={n_before - n_after} kept_ratio={kept_ratio:.4f}\n")

    if args.trace_tsv:
        with smart_open(args.trace_tsv, 'wt') as tf:
            tf.write('\t'.join([
                '#params',
                f"c_thresh={args.c_thresh}",
                f"rpf_metric={args.rpf_metric}",
                f"ps_metric={args.ps_metric}",
                f"rpf_min={args.rpf_min}",
                f"psite_min={args.psite_min}",
                f"evidence_mode={args.evidence_mode}"
            ]) + '\n')
            tf.write('\n'.join(trace_lines) + '\n')

if __name__ == '__main__':
    main()
