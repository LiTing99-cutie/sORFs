# #!/usr/bin/env bash
set -eou pipefail

FA=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.pep.fa
GP=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.genepred.txt
BAM=../processed/batch_1_2_merged/batch_1_2_merged.bam
OUT=../processed/filtering
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
featureCounts -T 30 -a "$OUT/sub.gtf" -t CDS -g transcript_id -O -s 1 -o "$OUT/counts.txt" "$BAM"

# 5) 取≥10 reads的transcript_id
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start filtering transcripts"
awk 'BEGIN{FS=OFS="\t"} !/^#|^Geneid/ && $7>=10 {print $1}' "$OUT/counts.txt" \
  > "$OUT/keep.ids"

# 6) 回切genePred与FASTA得到最终结果
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start filtering genePred"
awk 'NR==FNR{a[$1]=1; next} ($1 in a)' "$OUT/keep.ids" "$OUT/sub.genepred" \
  > "$OUT/final.genepred"

# 用seqkit速度更快（推荐）；无seqkit可用awk版本见下
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start filtering FASTA"
seqkit grep -n -f "$OUT/keep.ids" "$FA" > "$OUT/final.fa"

# 若没有seqkit，可用awk保留匹配记录：
# awk 'NR==FNR{keep[$1]=1; next}
#      /^>/{id=substr($0,2); printrec=(id in keep)}
#      {if(printrec) print}' "$OUT/keep.ids" "$FA" > "$OUT/final.fa"
