#!/usr/bin/env bash
set -euo pipefail

# 用法：id.convert.sh <file1.tsv> <map.tsv> <out.tsv>
# file1：两列或以上；第2列为要替换的ID
# map  ：两列；第1列=键( first_id_in_line )，第2列=值( chosen_id )
# 输出：
#   1) <out.tsv>        ：file1 的第2列按 map 替换后的结果
#   2) <out.tsv>.stats  ：统计信息

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <file1.tsv> <map.tsv> <out.tsv>" >&2
  exit 2
fi

pep_orf="$1"
id_map="$2"
out="$3"

mkdir -p "$(dirname "$out")"

awk -F'\t' -v OFS='\t' -v STATS="$out.stats" '
  # ------- 读映射表（map.tsv） -------
  FNR==NR{
    # 跳过map的表头（若首行含first_id_in_line/chosen_id/random等关键词）
    if (FNR==1) {
      hdr = tolower($1 OFS $2)
      if (hdr ~ /first_id_in_line/ || hdr ~ /chosen_id/ || hdr ~ /random/) next
    }
    gsub(/\r$/,"",$1); gsub(/\r$/,"",$2)
    if ($1!="") map[$1]=$2
    next
  }

  # ------- 处理 file1（pep_orf） -------
  FNR==1 {
    # 检测并保留file1表头（若有）
    hdr1 = tolower($1 OFS $2)
    has_header = (hdr1 ~ /peptide|first|random|id/)
    if (has_header) { print; next }
  }

  {
    gsub(/\r$/,"",$2)
    total_rows++                              # 统计：file1有效数据行
    if ($2 in map) {
      hit_rows++                              # 命中行数（按行）
      if (!(seen[$2]++)) unique_hits++        # 命中的唯一ID个数（可选）
      $2 = map[$2]
    } else {
      miss_rows++
    }
    print
  }

  END{
    # 写统计
    # 命中率 = 命中行 / 有效总行
    rate = (total_rows ? hit_rows/total_rows : 0)
    printf("Total_rows\t%d\nHit_rows\t%d\nMiss_rows\t%d\nUnique_hit_ids\t%d\nHit_rate\t%.4f\n",
           total_rows, hit_rows, miss_rows, unique_hits, rate) > STATS
  }
' "$id_map" "$pep_orf" > "$out"

echo "Done:"
echo "  Replaced file : $out"
echo "  Stats         : $out.stats"
