#!/usr/bin/env bash
set -euo pipefail

# -------- 路径参数 --------
path_1="/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/batch_1"
path_2="/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250813_demo_data_analysis/processed/batch_1"
STARindex="/home/user/data/lit/database/public/genome/hg38/hg38_custom_STARindex_v2.5.2b"

# -------- 输出目录 --------
outdir="../processed/trans_alignment_change_gtf"
logdir="../log/trans_alignment_change_gtf"
mkdir -p "$outdir" "$logdir"

# -------- 收集 FASTQ 并保存清单 --------
fq_pattern="*trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz"
mapfile -t FQS < <(
  find "$path_1" -type f -name "$fq_pattern" -print
  find "$path_2" -type f -name "$fq_pattern" -print
)
# 去重+排序（可选）
IFS=$'\n' read -r -d '' -a FQS < <(printf "%s\n" "${FQS[@]}" | sort -u && printf '\0')
# 保存清单
printf "%s\n" "${FQS[@]}" | tee "$logdir/fastq_list.txt" >/dev/null

# -------- 逐样本比对（仅输出 TranscriptomeSAM）--------
for cleaned_fastq in "${FQS[@]}"; do
  sample="$(basename "$cleaned_fastq" | sed -E 's/\.trimmed\.rRNA\.tRNA\.snoRNA\.unaligned\.fq\.gz$//')"
  echo "Mapping $sample at $(date '+%Y-%m-%d %H:%M:%S')"

  STAR --runThreadN 50 \
    --limitBAMsortRAM 250000000000 \
    --readFilesIn "$cleaned_fastq" \
    --readFilesCommand zcat \
    --seedSearchStartLmax 15 \
    --outFilterMismatchNmax 2 \
    --genomeDir "$STARindex" \
    --outFileNamePrefix "${outdir}/${sample}_" \
    --outSAMtype None \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --outSAMattributes All \
    --alignEndsType EndToEnd \
    &> "${logdir}/${sample}.STAR.log"
done

echo "Done. FASTQ list saved to ${logdir}/fastq_list.txt"
echo "Transcriptome BAM files are in ${outdir}/, named *_Aligned.toTranscriptome.out.bam"

