#!/usr/bin/env bash
set -euo pipefail

# ===================== 选项解析 =====================
COL_IDX=""
COL_NAME=""

usage() {
  echo "Usage: $0 [-c N | -n NAME] <input.tsv> <out_dir>"
  echo "  -c, --col N      指定 orf_id 的列号（从1开始）"
  echo "  -n, --name NAME  指定 orf_id 的列名（首行精确匹配）"
  exit 2
}

# 解析可选参数
while [[ $# -gt 0 ]]; do
  case "${1:-}" in
    -c|--col)   COL_IDX="${2:-}"; shift 2 ;;
    -n|--name)  COL_NAME="${2:-}"; shift 2 ;;
    -h|--help)  usage ;;
    --) shift; break ;;
    -*) echo "Unknown option: $1" >&2; usage ;;
    *)  break ;;
  esac
done

# ===================== 位置参数 =====================
IN_TSV="${1:?Usage: $0 [-c N| -n NAME] <input.tsv> <out_dir>}"
OUT_DIR="${2:?Usage: $0 [-c N| -n NAME] <input.tsv> <out_dir>}"

# 你的 genePred 与基因组（按你的路径）
GP="/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.genepred.txt"
GENOME="/home/user/data3/lit/project/sORFs/07-Genome/results/custom_fa/custom_ref.fa"
TRANSLATE_SH="/home/user/data3/lit/project/sORFs/01-ribo-seq/scripts_collapse_20250730/S3.0c.Uni.translate_gtf.v2.20250325.sh"

mkdir -p "$OUT_DIR"

# ===================== 1) 确定 ORF 列号 =====================
echo "[$(date '+%F %T')] Resolve ORF_id column"

if [[ -n "$COL_IDX" && -n "$COL_NAME" ]]; then
  echo "ERROR: 请只使用 -c 或 -n 其中一个。" >&2; exit 1
fi

if [[ -n "$COL_IDX" ]]; then
  # 用户直接给了列号
  ORF_COL="$COL_IDX"
else
  # 列名：优先用 --name，其次回退到默认 'ORF_id'
  WANT_NAME="${COL_NAME:-ORF_id}"
  ORF_COL="$(awk -F'\t' -v key="$WANT_NAME" 'NR==1{for(i=1;i<=NF;i++) if($i==key){print i; exit}}' "$IN_TSV" || true)"
  if [[ -z "${ORF_COL:-}" ]]; then
    echo "ERROR: 表头中未找到列名 '$WANT_NAME'。" >&2
    exit 1
  fi
fi

# 简单合法性检查
if ! [[ "$ORF_COL" =~ ^[0-9]+$ ]] || [[ "$ORF_COL" -lt 1 ]]; then
  echo "ERROR: 非法列号：$ORF_COL" >&2; exit 1
fi
echo "Use ORF_id column index: $ORF_COL"

# ===================== 1) 提取 ORF_id 列 =====================
echo "[$(date '+%F %T')] Extract ORF_id list from table"
# 如你的文件 **没有表头**，把下面 NR>1 改成 NR>=1
awk -v c="$ORF_COL" -F'\t' 'NR>1 && $c!="" && $c!="NA"{print $c}' "$IN_TSV" \
  | sort -u > "$OUT_DIR/ids.txt"

# ===================== 2) 用 ids 子集化 genePred =====================
echo "[$(date '+%F %T')] Filter genePred by ids"
ulimit -v 104857600
awk 'NR==FNR{a[$1]=1; next} ($1 in a)' "$OUT_DIR/ids.txt" "$GP" > "$OUT_DIR/sub.genepred"

# 导出 ORF_id→start/end
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$6,$7}' "$OUT_DIR/sub.genepred" > "$OUT_DIR/starts_ends.tsv"
# 格式：ORF_id  start  end

# ===================== 3) genePred → GTF =====================
echo "[$(date '+%F %T')] Convert genePred to GTF"
genePredToGtf file "$OUT_DIR/sub.genepred" "$OUT_DIR/sub.gtf"

# ===================== 4) 由 GTF+基因组生成 cds.fa/prot.fa =====================
echo "[$(date '+%F %T')] Build cds.fa/prot.fa from GTF and genome"
bash "$TRANSLATE_SH" \
  "$(realpath "$OUT_DIR/sub.gtf")" \
  "$GENOME" \
  "$OUT_DIR/translate_out"

CDS_FA="$OUT_DIR/translate_out/cds.fa"
if [[ ! -s "$CDS_FA" ]]; then
  echo "ERROR: 未找到 cds.fa（检查翻译脚本输出路径）" >&2
  exit 1
fi

# ===================== 5) cds.fa → ORF_id→CDS 映射 =====================
echo "[$(date '+%F %T')] Parse cds.fa to id→CDS map"
awk '
  BEGIN{FS="\t"; OFS="\t"}
  /^>/{
    if(id!=""){print id,seq}
    split(substr($0,2),a," "); id=a[1]; seq=""
    next
  }
  {gsub(/[ \t\r]/,""); seq=seq $0}
  END{ if(id!=""){print id,seq} }
' "$CDS_FA" > "$OUT_DIR/cds_map.tsv"
# 格式：ORF_id  CDS_seq

# ===================== 6) 回填 start/end/CDS =====================
echo "[$(date '+%F %T')] Merge start/end/CDS into original table (fill NA)"
awk -v c="$ORF_COL" -F'\t' -v OFS='\t' '
  BEGIN{
    # 读 starts_ends：ORF_id \t start \t end
    while((getline < ARGV[1])>0){ se[$1]=$2 OFS $3 } close(ARGV[1])
    # 读 cds_map：ORF_id \t CDS_seq
    while((getline < ARGV[2])>0){ cds[$1]=$2 } close(ARGV[2])
    ARGV[1]=ARGV[2]=""  # 丢弃前两个参数，后续处理 IN_TSV
  }
  NR==1{ print $0,"ORF_start","ORF_end","CDS_seq"; next }
  {
    key = $c
    s="NA"; e="NA"; seq="NA"
    if (key in se){
      split(se[key], a, OFS)
      if (a[1] != "" && a[1] != "NA") s=a[1]
      if (a[2] != "" && a[2] != "NA") e=a[2]
    }
    if (key in cds && cds[key] != "" && cds[key] != "NA") seq=cds[key]
    print $0, s, e, seq
  }
' "$OUT_DIR/starts_ends.tsv" "$OUT_DIR/cds_map.tsv" "$IN_TSV" > "$OUT_DIR/augmented.tsv"

echo "Done:"
echo "  IDs file      : $OUT_DIR/ids.txt"
echo "  genePred sub  : $OUT_DIR/sub.genepred"
echo "  GTF           : $OUT_DIR/sub.gtf"
echo "  cds_map       : $OUT_DIR/cds_map.tsv"
echo "  Output table  : $OUT_DIR/augmented.tsv"
