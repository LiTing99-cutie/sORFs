source activate biotools
BAM=/home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-20250621/mapping/merge/Total.bam
# 仅首次需要：更新参考
arcasHLA reference --update

# 变量
# BAM=/path/to/rnaseq.sorted.bam     # 需坐标排序并建立 .bai
OUT=arcas_out_rna                  # 输出目录
T=16                               # 线程

# 若尚未建索引
# samtools sort -@ $T $BAM > /home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-20250621/mapping/merge/Total.sorted.bam
# samtools index -@ $T /home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-20250621/mapping/merge/Total.sorted.bam
BAM=/home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-20250621/mapping/merge/Total.sorted.bam

# 提取 HLA 相关reads → 生成配对FASTQ
mkdir $OUT/tmp
arcasHLA extract "$BAM" -o "$OUT" -t $T --temp $OUT/tmp

# 分型（I+II一起做，II类关键基因已包含）
arcasHLA genotype \
  "$OUT"/*.extracted.1.fq.gz "$OUT"/*.extracted.2.fq.gz \
  -g A,B,C,DRB1,DQB1,DQA1,DPB1,DPA1 \
  -o "$OUT" -t $T --temp /home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/MS_res_from_Galaxy/HLA/arcas_out_rna/tmp

# 查看结果（简表 & JSON）
cat "$OUT"/sample.genotype.txt
# 详细：$OUT/sample.genotype.json
