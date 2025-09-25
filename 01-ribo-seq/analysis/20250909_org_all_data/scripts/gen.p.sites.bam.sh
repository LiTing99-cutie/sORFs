#!/usr/bin/env bash
set -euo pipefail

# ======== 配置 ========
BAM_LIST="${1:-bam.lst}"                 # 每行一个 BAM 的绝对路径
OFFSET_ALL="${2:-offset.tab.all.txt}"    # 三列：sample  readlen  offset
OUT_DIR="${3:-./riborf_offset_fixed}"    # 输出目录
THREADS="${4:-40}"                        # samtools 线程
RIBORF_PATH="${RibORF_path:-/home/user/data2/lit/software/RibORF/RibORF.2.0}"

mkdir -p "$OUT_DIR" "$OUT_DIR/_tmp"

# 工具检查
command -v samtools >/dev/null 2>&1 || { echo "samtools 未找到"; exit 1; }
[[ -s "$RIBORF_PATH/offsetCorrect.pl" ]] || { echo "找不到 $RIBORF_PATH/offsetCorrect.pl"; exit 1; }

# 逐个样本处理
while IFS= read -r bam; do
  [[ -z "$bam" ]] && continue
  [[ ! -s "$bam" ]] && { echo "WARN: BAM 不存在或为空：$bam"; continue; }

  base="$(basename "$bam")"
  sample="${base%_Aligned.sortedByCoord.out.bam}"
  sample="${sample%.bam}"

  # 为该样本生成两列 offset 表（readlen \t offset）
  off_tmp="$OUT_DIR/_tmp/${sample}.offset.tab"
  awk -v s="$sample" 'BEGIN{OFS="\t"} $1==s {print $2,$3}' "$OFFSET_ALL" > "$off_tmp"
  if [[ ! -s "$off_tmp" ]]; then
    echo "WARN: 未在 $OFFSET_ALL 中找到样本 $sample 的 offset，跳过。"
    rm -f "$off_tmp"
    continue
  fi

  # 中间与输出路径
  sam_in="$OUT_DIR/_tmp/${sample}.in.sam"
  sam_out="$OUT_DIR/${sample}.offsetCorrected.sam"
  bam_out="${sam_out%.sam}.bam"

  echo "[`date +%F' '%T`] Processing: $sample"
  # BAM -> SAM（带头信息）
  samtools view -@ "$THREADS" -h -o "$sam_in" "$bam"

  # 纠偏
  perl "$RIBORF_PATH/offsetCorrect.pl" -r "$sam_in" -p "$off_tmp" -o "$sam_out"

  # SAM -> BAM 并清理 SAM
  samtools view -@ "$THREADS" -bh "$sam_out" > "$bam_out"
  rm -f "$sam_in" "$sam_out" "$off_tmp"

done < "$BAM_LIST"

echo "全部完成。输出 BAM 位于：$OUT_DIR"

# 1) 列出所有待合并 BAM 并统计数量
ls -1 "$OUT_DIR"/*.offsetCorrected.bam > "$OUT_DIR/merge.inputs.txt"
N=$(wc -l < "$OUT_DIR/merge.inputs.txt")
echo "将合并 $N 个 BAM"

# 2) 合并后立刻排序，再建索引（管道更省IO）
samtools merge -@ 40 -b "$OUT_DIR/merge.inputs.txt" -O BAM - \
| samtools sort -@ 40 -o "$OUT_DIR/all.offsetCorrected.merged.sorted.bam" -

samtools index "$OUT_DIR/all.offsetCorrected.merged.sorted.bam"

# 4) 统计合并后总reads数
samtools view -c "$OUT_DIR/all.offsetCorrected.merged.bam"
