# #!/usr/bin/env bash
set -eou pipefail

FA=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.pep.fa
GP=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.genepred.txt
BAM_1=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250827_custom_db_filtering/processed/batch_1_2_merged/batch_1_2_merged.bam
BAM_2=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250909_org_all_data/processed/all.offsetCorrected.merged.sorted.bam
OUT=../processed/ribo
mkdir -p "$OUT"

# 1) 从FASTA提ID（仅取>后第一个字段），去重
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start extracting IDs from FASTA"
seqkit seq -n $FA > "$OUT/ids.txt"

# 2) 以ids作为索引筛选genePred（大文件友好：只把ids进内存）
# 100G
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start filtering genePred"
ulimit -v 104857600
awk 'NR==FNR{a[$1]=1; next} ($1 in a)' "$OUT/ids.txt" "$GP" \
  > "$OUT/sub.genepred"

# 3) genePred→GTF（UCSC kentutils）
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start converting genePred to GTF"
# 若是标准10列genePred：
genePredToGtf file "$OUT/sub.genepred" "$OUT/sub.gtf"
# 若是12列genePredExt（含frame/其它字段），改用：
# genePredToGtf -genePredExt -source=custom "$OUT/sub.genepred" "$OUT/sub.gtf"

# 4) 计数（按CDS，按transcript_id汇总；确认你的文库方向性-s 1是否正确）
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start counting"
source activate biotools
featureCounts -T 30 -a "$OUT/sub.gtf" -t CDS -g transcript_id -O -s 1 -o "$OUT/counts.txt" "$BAM_1" "$BAM_2"