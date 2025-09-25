#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert custom FASTA to TSV with columns:
ID  Seq  Length  Chr  Strand  Transcript_id  Type  Scodon

- ID: the full FASTA name (without '>')
- Transcript_id: substring before the first ':' in the FASTA header
"""
import sys, gzip, argparse
from typing import Iterator, Tuple, Optional

def openmaybe(path: Optional[str]):
    if path is None or path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")

def parse_header(h: str) -> Tuple[str, str, str, str, str, str]:
    """
    Parse headers like:
    >PB.41013.1:chr6:+|1|330:5:44|noncoding|GTG
    >ENST00000832824.1:chr1:+|1|1379:8:230|noncoding|ACG
    Returns: (ID, Chr, Strand, Type, Scodon, Transcript_id)
    """
    h = h.strip()
    if h.startswith(">"):
        h = h[1:]
    _id = h
    transcript_id = h.split(":", 1)[0]
    parts = h.split(":", 2)
    _chr = parts[1] if len(parts) > 1 else ""
    rest = parts[2] if len(parts) > 2 else ""
    rfields = rest.split("|") if rest else []
    strand = rfields[0] if len(rfields) >= 1 else ""
    _type = rfields[-2] if len(rfields) >= 2 else ""
    scodon = rfields[-1] if len(rfields) >= 1 else ""
    return _id, _chr, strand, _type, scodon, transcript_id

def records(handle) -> Iterator[Tuple[str, str]]:
    header = None
    seqlist = []
    for line in handle:
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seqlist).replace(" ", "").replace("\t", "").replace("\r", "").replace("\n", "")
            header = line.strip()
            seqlist = []
        else:
            seqlist.append(line.strip())
    if header is not None:
        yield header, "".join(seqlist).replace(" ", "").replace("\t", "").replace("\r", "").replace("\n", "")

def main():
    ap = argparse.ArgumentParser(description="Custom FASTA -> TSV")
    ap.add_argument("fasta", help="Input FASTA (use - for stdin). .gz supported.")
    ap.add_argument("-o", "--out", default="-", help="Output TSV path (default: stdout)")
    ap.add_argument("--no-seq", action="store_true", help="Do not output the Seq column")
    args = ap.parse_args()
    inp = openmaybe(args.fasta)
    out = sys.stdout if args.out in ("-", None) else open(args.out, "w", encoding="utf-8")
    try:
        if args.no_seq:
            header_line = "\t".join(["ID","Length","Chr","Strand","Transcript_id","Type","Scodon"])
        else:
            header_line = "\t".join(["ID","Seq","Length","Chr","Strand","Transcript_id","Type","Scodon"])
        print(header_line, file=out)
        for h, seq in records(inp):
            _id, _chr, strand, _type, scodon, transcript_id = parse_header(h)
            length = len(seq)
            if args.no_seq:
                row = [_id, str(length), _chr, strand, transcript_id, _type, scodon]
            else:
                row = [_id, seq, str(length), _chr, strand, transcript_id, _type, scodon]
            print("\t".join(row), file=out)
    finally:
        if inp is not sys.stdin:
            inp.close()
        if out is not sys.stdout:
            out.close()

if __name__ == "__main__":
    main()
