#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parallel MAF extractor using Biopython MafIO.MafIndex.

- Per chromosome: (optionally) unzip chrN.maf.gz -> chrN.maf once, build SQLite index once, then reuse.
- Per ORF BED (1-based, multi-exon): use MafIndex.get_spliced() to export per-ORF alignment/FASTA,
  and MafIndex.search() per exon to derive non-focal species BED coordinates (block-aware).

Requirements:
  - Biopython >= 1.75 (MafIO.MafIndex with search/get_spliced)
  - Large disk space for unzipped .maf (once per chr)
"""

import os
import sys
import io
import re
import gzip
import argparse
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio.AlignIO import MafIO
from Bio import AlignIO
from Bio.Seq import Seq

# ----------------------------- Utils -----------------------------
def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def sanitize_name(s: str) -> str:
    """Safe filename: keep alnum, . _ -, collapse others to '__'."""
    return re.sub(r"[^A-Za-z0-9._-]+", "__", s)

def find_first(iterable, pred):
    for x in iterable:
        if pred(x):
            return x
    return None

def detect_maf_and_prepare(maf_root: Path, chr_name: str, work_dir: Path) -> Path:
    """
    Return path to an UNCOMPRESSED chrN.maf.
    If only .maf.gz exists, unzip once into work_dir/_maf_unzip/chrN.maf (reuse later).
    """
    src_maf = maf_root / f"{chr_name}.maf"
    if src_maf.exists():
        return src_maf

    src_gz = maf_root / f"{chr_name}.maf.gz"
    if not src_gz.exists():
        raise FileNotFoundError(f"Cannot find {src_maf} or {src_gz}")

    dst_dir = work_dir / "_maf_unzip"
    ensure_dir(dst_dir)
    dst_maf = dst_dir / f"{chr_name}.maf"

    # unzip if not exists or gz is newer
    if (not dst_maf.exists()) or (src_gz.stat().st_mtime > dst_maf.stat().st_mtime):
        print(f"[INFO] Unzipping {src_gz} -> {dst_maf}")
        # Use system gzip for speed, fallback to Python if unavailable
        try:
            subprocess.run(["gzip", "-dc", str(src_gz)], check=True, stdout=open(dst_maf, "wb"))
        except Exception:
            with gzip.open(src_gz, "rb") as fin, open(dst_maf, "wb") as fout:
                while True:
                    chunk = fin.read(1024 * 1024)
                    if not chunk:
                        break
                    fout.write(chunk)
    return dst_maf

def build_or_load_maf_index(maf_path: Path, chr_name: str, focal: str, work_dir: Path) -> MafIO.MafIndex:
    """
    Create/reuse SQLite index under work_dir/_maf_index/{chr}.sqlite
    target_seqname must match MAF src like 'hg38.chr1'.
    """
    idx_dir = work_dir / "_maf_index"
    ensure_dir(idx_dir)
    sqlite_path = idx_dir / f"{chr_name}.sqlite"
    target_seqname = f"{focal}.{chr_name}"  # e.g., hg38.chr1

    # MafIndex will build index if sqlite doesn't exist (or is empty)
    idx = MafIO.MafIndex(str(sqlite_path), str(maf_path), target_seqname)
    return idx

def parse_orf_bed_1based_single(bed_path: Path):
    """
    Read a single-ORF 1-based BED (multi-exon rows).
    Return: dict with keys: chrm(str), orf_id(str), strand('+'|'-'), starts(list[int],1-based), ends(list[int],1-based)
    """
    chrs, starts, ends, strands, ids = [], [], [], [], []
    with bed_path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            ps = line.rstrip("\n").split("\t")
            if len(ps) < 6:
                raise ValueError(f"BED columns <6: {bed_path}\n{line}")
            chrs.append(ps[0] if ps[0].startswith("chr") else "chr"+ps[0])
            starts.append(int(ps[1]))
            ends.append(int(ps[2]))
            ids.append(ps[3])
            strands.append(ps[5])

    if not starts:
        raise ValueError(f"Empty BED: {bed_path}")

    chrm = chrs[0]
    strand = strands[0]
    orf_id = ids[0]
    # pairwise sort by start
    exons = sorted(zip(starts, ends), key=lambda x: x[0])
    starts = [s for s, _ in exons]
    ends = [e for _, e in exons]
    return {"chrm": chrm, "orf_id": orf_id, "strand": strand, "starts_1b": starts, "ends_1b": ends}

def one_based_to_zero_half_open(starts_1b, ends_1b):
    """Convert 1-based closed [s,e] to 0-based half-open [s-1, e]."""
    starts0 = [s - 1 for s in starts_1b]
    ends0 = [e for e in ends_1b]
    return starts0, ends0

# ----------------------------- BED generation via block mapping -----------------------------
def species_bed_for_orf(idx: MafIO.MafIndex, chr_name: str, focal: str, orf_meta: dict):
    """
    Reproduce old-script-like BED lines for NON-focal species by scanning blocks per exon.
    Return: dict sp -> list of (chrom, start0, end0, orf_id, sp(+_splitted), strandChar)
    """
    result = {}  # sp -> list of tuples
    starts0, ends0 = one_based_to_zero_half_open(orf_meta["starts_1b"], orf_meta["ends_1b"])
    target_name = f"{focal}.{chr_name}"

    # For _splitted detection
    sp_chr_sets = {}

    for s0, e0 in zip(starts0, ends0):
        # search returns MultipleSeqAlignment per overlapping MAF block
        for aln in idx.search([s0], [e0]):
            # find focal record in this block
            focal_rec = find_first(aln, lambda r: r.id == target_name)
            if focal_rec is None:
                continue
            regin = int(focal_rec.annotations.get("start", 0))  # 0-based start within focal
            # map columns range [p1,p2] that correspond to [s0,e0)
            p1 = -1
            p2 = -1
            counter = regin
            seqF = str(focal_rec.seq)
            for col, ch in enumerate(seqF):
                if (counter >= s0) and (counter <= e0 - 1):
                    if p1 == -1:
                        p1 = col
                elif (counter > e0 - 1) and (p2 == -1):
                    p2 = col - 1
                if ch != "-":
                    counter += 1
            if p1 == -1:
                # no overlap
                continue
            if p2 == -1:
                p2 = len(seqF) - 1

            # now for each species compute genomic spans within this block slice
            for rec in aln:
                sp_full = rec.id  # e.g., 'hg38.chr1'
                sp = sp_full.split(".")[0]
                rec_chr = sp_full.split(".", 1)[1] if "." in sp_full else sp_full
                # count non-gaps up to p1/p2 in THIS species
                text = str(rec.seq)
                p11 = p22 = 0
                for i, ch in enumerate(text):
                    if (i < p1) and (ch != "-"):
                        p11 += 1
                    if (i < p2) and (ch != "-"):
                        p22 += 1
                rec_start0 = int(rec.annotations.get("start", 0))
                strand = rec.annotations.get("strand", 1)  # +1 / -1
                strand_char = "-" if str(strand) in ("-1", "-") else "+"

                # record spans (0-based) using start0+p11 / start0+p22 (consistent with旧脚本)
                tup = (rec_chr, rec_start0 + p11, rec_start0 + p22, orf_meta["orf_id"], sp, strand_char)
                result.setdefault(sp, []).append(tup)
                sp_chr_sets.setdefault(sp, set()).add(rec_chr)

    # append "_splitted" if species spans multiple chromosomes across blocks
    final = {}
    for sp, entries in result.items():
        splitted = "_splitted" if len(sp_chr_sets.get(sp, set())) != 1 else ""
        final[sp] = [(chrom, s, e, orf_id, sp + splitted, strand) for (chrom, s, e, orf_id, sp, strand) in entries]
    return final

# ----------------------------- Per-ORF worker -----------------------------
def process_one_orf(bed_path: Path, idx_sqlite: Path, maf_path: Path, chr_name: str,
                    focal: str, out_chr_dir: Path, cov_thresh: float = 0.95) -> str:
    """
    Process a single ORF bed:
      - parse bed (1-based) -> exons
      - get_spliced alignment for all exons (strand-aware)
      - write per-ORF .fa / .maf
      - write non-focal species .bed using block mapping (search per exon)
    """
    # parse ORF meta
    orf_meta = parse_orf_bed_1based_single(bed_path)
    orf_id = orf_meta["orf_id"]
    orf_id_safe = sanitize_name(orf_id)
    strand = orf_meta["strand"]

    out_orfs = out_chr_dir / "orfs"
    out_beds = out_chr_dir / "beds"
    ensure_dir(out_orfs); ensure_dir(out_beds)

    fa_out  = out_orfs / f"{orf_id_safe}.fa"
    maf_out = out_orfs / f"{orf_id_safe}.maf"
    bed_out = out_beds / f"{orf_id_safe}.bed"
    if fa_out.exists() and maf_out.exists() and bed_out.exists():
        return f"[SKIP] {orf_id}"

    # open index (cheap if sqlite exists)
    target_seqname = f"{focal}.{chr_name}"
    idx = MafIO.MafIndex(str(idx_sqlite), str(maf_path), target_seqname)

    # get spliced alignment for the tumor ORF transcript
    starts0, ends0 = one_based_to_zero_half_open(orf_meta["starts_1b"], orf_meta["ends_1b"])
    strand_flag = 1 if strand == "+" else -1
    try:
        aln = idx.get_spliced(starts0, ends0, strand=strand_flag)
    except Exception as e:
        return f"[FAIL] {orf_id}: get_spliced error: {e}"

    # sequences per species (use species-only key)
    seq_by_sp = {}
    focal_key_prefix = f"{focal}."
    focal_seq = None
    for rec in aln:
        sp_full = rec.id  # e.g., 'hg38.chr1'
        sp = sp_full.split(".")[0]
        seq = str(rec.seq)
        seq_by_sp[sp] = seq
        if sp_full.startswith(focal_key_prefix):
            focal_seq = seq

    if not focal_seq:
        return f"[MISS] {orf_id}\t0"

    # coverage m = |focal ungapped| / genomic length L
    L = sum(e - s + 1 for s, e in zip(orf_meta["starts_1b"], orf_meta["ends_1b"]))
    m = float(len(focal_seq.replace("-", "")) / float(L))

    # write .fa / .maf (reverse_complement all if ORF on '-')
    with fa_out.open("w") as f1, maf_out.open("w") as f2:
        focal_len = len(focal_seq)
        for sp, sq in seq_by_sp.items():
            if focal_len == 0:
                continue
            # coverage filter relative to focal aligned length
            if (len(sq) / focal_len) < cov_thresh:
                continue
            sq_out = str(Seq(sq).reverse_complement()) if strand == "-" else sq
            # .fa: translate ungapped (cds=False, same as old code)
            f1.write(f">{orf_id_safe}_{sp}_{m:.6f}\n{str(Seq(sq_out.replace('-','')).translate(cds=False))}\n")
            # .maf-like fasta: per-species aligned nucleotides
            f2.write(f">{sp}\n{sq_out}\n")

    # write non-focal species BED via block mapping
    sp_bed = species_bed_for_orf(idx, chr_name, focal, orf_meta)
    with bed_out.open("w") as f3:
        for sp, items in sp_bed.items():
            if sp == focal:
                continue
            for chrom, s0, e0, oid, sp_tag, strand_char in items:
                f3.write(f"{chrom}\t{s0}\t{e0}\t{orf_id}\t{sp_tag}\t{strand_char}\n")

    return f"[OK] {orf_id}"

# ----------------------------- Main: per-chromosome loop + inner parallel -----------------------------
def main():
    ap = argparse.ArgumentParser(description="Per-chromosome parallel ORF extractor using Biopython MafIO.MafIndex")
    ap.add_argument("--per-chr-bed-root", required=True,
                    help="Root dir containing per-chr ORF BED subfolders, e.g. ../processed/get_input/per_chr/ORFs_bed")
    ap.add_argument("--maf-root", required=True, help="Folder containing chrN.maf or chrN.maf.gz")
    ap.add_argument("--output-root", required=True, help="Output root; per-chr subfolders will be created")
    ap.add_argument("--work-dir", required=True, help="Working dir for unzipped MAF and SQLite index (persistent)")
    ap.add_argument("--chr", default="", help="Run only this chromosome (e.g., chr1); empty = all chr* under per-chr-bed-root")
    ap.add_argument("--focal", default="hg38", help="Focal species prefix (e.g., hg38)")
    ap.add_argument("--workers", type=int, default=5, help="Parallel workers per chromosome (default 5)")
    ap.add_argument("--cov-thresh", type=float, default=0.95, help="Coverage filter relative to focal alignment length")
    args = ap.parse_args()

    per_chr_root = Path(args.per_chr_bed_root)
    maf_root = Path(args.maf_root)
    out_root = Path(args.output_root)
    work_dir = Path(args.work_dir)

    # determine chr folders
    if args.chr:
        chr_dirs = [per_chr_root / args.chr]
    else:
        # auto-detect subdirs like chr1, chr2, ...
        chr_dirs = sorted([p for p in per_chr_root.glob("chr*") if p.is_dir()])

    if not chr_dirs:
        sys.exit(f"[ERROR] No chromosome subdirs found under {per_chr_root}")

    for chr_dir in chr_dirs:
        chr_name = chr_dir.name
        print(f"[INFO] === Chromosome {chr_name} ===")
        try:
            maf_path = detect_maf_and_prepare(maf_root, chr_name, work_dir)
        except Exception as e:
            print(f"[WARN] Skip {chr_name}: {e}")
            continue

        # Build or load index once; keep sqlite path to pass to workers
        idx = build_or_load_maf_index(maf_path, chr_name, args.focal, work_dir)
        # We re-open index inside workers by giving them sqlite path and maf path
        idx_sqlite = (work_dir / "_maf_index" / f"{chr_name}.sqlite")

        # Gather ORF beds (filenames may contain '|' ':' etc., Pathlib handles it fine)
        bed_files = sorted(chr_dir.glob("*.bed"))
        if not bed_files:
            print(f"[INFO] {chr_name}: no ORF beds found, skip.")
            continue

        out_chr_dir = out_root / chr_name
        ensure_dir(out_chr_dir / "orfs"); ensure_dir(out_chr_dir / "beds")

        print(f"[INFO] {chr_name}: {len(bed_files)} ORFs, workers={args.workers}")
        results = []
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = [
                ex.submit(
                    process_one_orf,
                    bed,
                    idx_sqlite,
                    maf_path,
                    chr_name,
                    args.focal,
                    out_chr_dir,
                    args.cov_thresh,
                )
                for bed in bed_files
            ]
            for i, fu in enumerate(as_completed(futs), 1):
                try:
                    msg = fu.result()
                except Exception as e:
                    msg = f"[FAIL] {e}"
                results.append(msg)
                if i % 50 == 0 or i == len(bed_files):
                    print(f"[PROG] {chr_name}: {i}/{len(bed_files)}")

        ok   = sum(1 for x in results if x.startswith("[OK]"))
        skip = sum(1 for x in results if x.startswith("[SKIP]"))
        miss = sum(1 for x in results if x.startswith("[MISS]"))
        fail = sum(1 for x in results if x.startswith("[FAIL]"))
        print(f"[DONE] {chr_name}: OK={ok}, SKIP={skip}, MISS={miss}, FAIL={fail}, TOTAL={len(bed_files)}")

if __name__ == "__main__":
    main()
