#!/usr/bin/env python3
# This script formats RiboCode *.txt to FASTA
# Usage: python RiboCode_fasta_formatting.py Input_RiboCode_output.txt RiboCode_fasta_formatted.fa
# 仅按 ORF_ID 的“倒数 3 个下划线”进行切割（rsplit('_', n=3)），前缀不管

import sys
import pandas as pd

if len(sys.argv) < 3:
    raise SystemExit("Usage: python RiboCode_fasta_formatting.py <input.txt> <output.fa>")

input_txt  = sys.argv[1]
output_fa  = sys.argv[2]

# 读入（严格按制表符；全部按字符串）
df = pd.read_csv(input_txt, sep="\t", header=0, dtype=str, engine="python")

# 必需列检查
need_cols = ["ORF_ID", "transcript_id", "strand", "chrom", "AAseq"]
missing = [c for c in need_cols if c not in df.columns]
if missing:
    raise SystemExit(f"Missing required columns: {missing}")

# 基础清洗
for c in need_cols:
    df[c] = df[c].astype(str).str.strip()

# 正负链拆分
pos_df = df[df["strand"] == "+"].copy()
neg_df = df[df["strand"] == "-"].copy()

# 仅从右侧切 3 次：得到 4 列 = [前缀任意内容, Start, Stop, Extra]
def split_from_right(series: pd.Series) -> pd.DataFrame:
    out = series.str.rsplit("_", n=3, expand=True)
    # 若不足 4 列（极少数异常），补齐
    if out.shape[1] < 4:
        out = out.reindex(columns=range(4))
    out.columns = ["Prefix", "Start", "Stop", "Extra"]
    return out

# 正链解析
if not pos_df.empty:
    pos_split = split_from_right(pos_df["ORF_ID"])
    # 转数值并容错
    pos_split["Start"] = pd.to_numeric(pos_split["Start"], errors="coerce")
    pos_split["Stop"]  = pd.to_numeric(pos_split["Stop"],  errors="coerce")
    mask_pos = pos_split["Start"].notna() & pos_split["Stop"].notna() & pos_df["AAseq"].notna()
    pos_split = pos_split[mask_pos]
    pos_df    = pos_df.loc[mask_pos]

    # 坐标规则：+链 Start = Start - 1，Stop 保持
    pos_start_hdr = (pos_split["Start"] - 1).astype("Int64").astype(str)
    pos_stop_hdr  = pos_split["Stop"].astype("Int64").astype(str)

    pos_headers = (
        ">" + pos_df["transcript_id"] + pos_df["strand"] + pos_df["chrom"] +
        ":" + pos_start_hdr + "-" + pos_stop_hdr
    )
    pos_seqs = pos_df["AAseq"].astype(str)
else:
    pos_headers = pd.Series([], dtype=str)
    pos_seqs    = pd.Series([], dtype=str)

# 负链解析
if not neg_df.empty:
    neg_split = split_from_right(neg_df["ORF_ID"])
    neg_split["Start"] = pd.to_numeric(neg_split["Start"], errors="coerce")
    neg_split["Stop"]  = pd.to_numeric(neg_split["Stop"],  errors="coerce")
    mask_neg = neg_split["Start"].notna() & neg_split["Stop"].notna() & neg_df["AAseq"].notna()
    neg_split = neg_split[mask_neg]
    neg_df    = neg_df.loc[mask_neg]

    # 坐标规则：-链 Start = Stop - 1，Stop = Start（并交换顺序）
    neg_left_hdr  = (neg_split["Stop"] - 1).astype("Int64").astype(str)  # 左端
    neg_right_hdr = neg_split["Start"].astype("Int64").astype(str)       # 右端

    neg_headers = (
        ">" + neg_df["transcript_id"] + neg_df["strand"] + neg_df["chrom"] +
        ":" + neg_left_hdr + "-" + neg_right_hdr
    )
    neg_seqs = neg_df["AAseq"].astype(str)
else:
    neg_headers = pd.Series([], dtype=str)
    neg_seqs    = pd.Series([], dtype=str)

# 写出 FASTA
with open(output_fa, "w") as fh:
    for h, s in zip(pd.concat([pos_headers, neg_headers]).tolist(),
                    pd.concat([pos_seqs,    neg_seqs]).tolist()):
        fh.write(h + "\n")
        fh.write(s + "\n")
