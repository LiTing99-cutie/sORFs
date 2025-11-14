#!/usr/bin/env bash
set -euo pipefail

# ===================== 参数 =====================
IN_TSV="${1:?Usage: $0 <input.tsv> <out_dir>}"
OUT_DIR="${2:?Usage: $0 <input.tsv> <out_dir>}"

# 你的 genePred 与基因组（按你的路径）
GP="/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.genepred.txt"
GENOME="/home/user/data3/lit/project/sORFs/07-Genome/results/custom_fa/custom_ref.fa"
TRANSLATE_SH="/home/user/data3/lit/project/sORFs/01-ribo-seq/scripts_collapse_20250730/S3.0c.Uni.translate_gtf.v2.20250325.sh"

mkdir -p "$OUT_DIR"

# ===================== 1) 提取 ORF_id 列 =====================
echo "[$(date '+%F %T')] Extract ORF_id list from table"
# 自动定位 ORF_id 列号（首行精确匹配）
ORF_COL=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="ORF_id") {print i; exit}}' "$IN_TSV")
if [[ -z "${ORF_COL:-}" ]]; then
  echo "ERROR: 表头中未找到列名 ORF_id" >&2
  exit 1
fi

# 取 ORF_id 非空/非NA 的行
awk -v c="$ORF_COL" -F'\t' 'NR>1 && $c!="" && $c!="NA"{print $c}' "$IN_TSV" \
  | sort -u > "$OUT_DIR/ids.txt"

# ===================== 2) 用 ids 子集化 genePred（仅把 ids 进内存） =====================
echo "[$(date '+%F %T')] Filter genePred by ids"
ulimit -v 104857600
awk 'NR==FNR{a[$1]=1; next} ($1 in a)' "$OUT_DIR/ids.txt" "$GP" > "$OUT_DIR/sub.genepred"

# 说明：按你的文件说明，“第6列和第7列分别是 start 和 end”
# 这里顺便导出一个 ORF_id→start/end 的映射表
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$6,$7}' "$OUT_DIR/sub.genepred" > "$OUT_DIR/starts_ends.tsv"
# 格式：ORF_id  start  end

# ===================== 3) genePred → GTF =====================
echo "[$(date '+%F %T')] Convert genePred to GTF"
genePredToGtf file "$OUT_DIR/sub.genepred" "$OUT_DIR/sub.gtf"

# ===================== 4) 由 GTF+基因组生成 cds.fa/prot.fa =====================
echo "[$(date '+%F %T')] Build cds.fa/prot.fa from GTF and genome"
bash "$TRANSLATE_SH" \
  $(realpath "$OUT_DIR/sub.gtf") \
  "$GENOME" \
  "$OUT_DIR/translate_out"

# 约定 translate 脚本输出：
#   $OUT_DIR/translate_out/cds.fa
#   $OUT_DIR/translate_out/prot.fa
# 若实际路径/名称不同，请相应修改上面调用或下面的路径。

CDS_FA="$OUT_DIR/translate_out/cds.fa"
if [[ ! -s "$CDS_FA" ]]; then
  echo "ERROR: 未找到 cds.fa（检查翻译脚本输出路径）" >&2
  exit 1
fi

# ===================== 5) 将 cds.fa 转为 ORF_id→CDS 序列 映射 =====================
echo "[$(date '+%F %T')] Parse cds.fa to id→CDS map"
# 取fasta id的第一个“空格前字段”为 ORF_id；把序列拼接成一行
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

# ===================== 6) 把 start/end 和 CDS_seq 回填到原表（缺失填 NA） =====================
echo "[$(date '+%F %T')] Merge start/end/CDS into original table (fill NA)"
awk -v c="$ORF_COL" -F'\t' -v OFS='\t' '
  BEGIN{
    # 读 starts_ends：ORF_id \t start \t end
    while((getline < ARGV[1])>0){
      se[$1]=$2 OFS $3
    }
    close(ARGV[1])
    # 读 cds_map：ORF_id \t CDS_seq
    while((getline < ARGV[2])>0){
      cds[$1]=$2
    }
    close(ARGV[2])
    ARGV[1]=ARGV[2]=""   # 丢弃前两个参数，后续处理 IN_TSV
  }
  NR==1{
    print $0,"ORF_start","ORF_end","CDS_seq"; next
  }
  {
    key = $c
    # 默认 NA
    s = "NA"; e = "NA"; seq = "NA"

    # 有 start/end 记录则覆盖；同时排除空字符串
    if (key in se) {
      split(se[key], arr, OFS)
      if (arr[1] != "" && arr[1] != "NA") s = arr[1]
      if (arr[2] != "" && arr[2] != "NA") e = arr[2]
    }
    # 有 CDS 记录则覆盖；同时排除空字符串
    if (key in cds && cds[key] != "" && cds[key] != "NA") {
      seq = cds[key]
    }

    print $0, s, e, seq
  }
' "$OUT_DIR/starts_ends.tsv" "$OUT_DIR/cds_map.tsv" "$IN_TSV" > "$OUT_DIR/augmented.tsv"

