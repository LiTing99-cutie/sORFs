#!/usr/bin/env bash
set -euo pipefail

show_help() { cat <<'USAGE'
用法：id.convert.sh [选项] <file1.tsv> <map.tsv> <out.tsv>
说明：将 file1 的某一列（默认第2列）按 map 映射替换，并输出统计信息 <out.tsv>.stats
选项：
  -c, --col N        指定 file1 中待替换 ID 的列号（从1开始，默认2）
      --header MODE  file1 首行是否为表头：auto|yes|no（默认 auto）
  -h, --help         显示帮助
USAGE
}

# 默认参数
COL=2
HEADER_MODE=auto

# 解析参数
args=()
while (($#)); do
  case "${1:-}" in
    -c|--col) COL="${2:?}"; shift 2;;
    --header) HEADER_MODE="${2:?}"; shift 2;;
    -h|--help) show_help; exit 0;;
    --) shift; break;;
    -*)
      echo "未知选项：$1" >&2; show_help; exit 2;;
    *)
      args+=("$1"); shift;;
  esac
done

if (( ${#args[@]} != 3 )); then
  show_help >&2
  exit 2
fi

pep_orf="${args[0]}"
id_map="${args[1]}"
out="${args[2]}"

mkdir -p "$(dirname "$out")"

awk -v FS='\t' -v OFS='\t' \
    -v STATS="$out.stats" \
    -v COL="$COL" \
    -v HEADER_MODE="$HEADER_MODE" '
  BEGIN{
    c = int(COL); if (c<1) c=2; COL=c
    total_rows=hit_rows=miss_rows=unique_hits=0
  }

  # ---------- 读 map ----------
  FNR==NR{
    # 跳过可能的表头
    if (FNR==1){
      line = tolower($1 OFS $2)
      if (line ~ /first_id_in_line|first|chosen_id|target|random|id/) next
    }
    gsub(/\r$/,"",$1); gsub(/\r$/,"",$2)
    if ($1!="") map[$1]=$2
    next
  }

  # ---------- 处理 file1 ----------
  FNR==1 {
    # 表头判定
    is_header = 0
    if (HEADER_MODE=="yes") {
      is_header=1
    } else if (HEADER_MODE=="auto") {
      line = tolower($0)
      # 只基于关键词判断，避免误杀
      if (line ~ /(^|[\t ])[#"]?(peptide|id|orf|protein|first|random)([\t ]|$)/) {
        is_header=1
      }
    }
    if (is_header) { print; next }
  }

  {
    # 清理 CRLF
    for (i=1;i<=NF;i++) sub(/\r$/,"",$i)

    # 列边界检查
    if (COL>NF) {
      total_rows++; miss_rows++; print; next
    }

    i = COL
    orig = $i
    total_rows++
    if (orig in map) {
      hit_rows++
      if (!(seen[orig]++)) unique_hits++
      $i = map[orig]
    } else {
      miss_rows++
    }
    print
  }

  END{
    rate = (total_rows ? hit_rows/total_rows : 0)
    printf("ID_col\t%d\nHeader_mode\t%s\nTotal_rows\t%d\nHit_rows\t%d\nMiss_rows\t%d\nUnique_hit_ids\t%d\nHit_rate\t%.4f\n",
           COL, HEADER_MODE, total_rows, hit_rows, miss_rows, unique_hits, rate) > STATS
  }
' "$id_map" "$pep_orf" > "$out"

echo "Done:"
echo "  Replaced file : $out"
echo "  Stats         : $out.stats"
