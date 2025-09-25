# 输入
BAM=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250909_org_all_data/processed/all.offsetCorrected.merged.sorted.bam
GTF=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/results/1/augment_orf_table/sub.gtf
output_path=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/processed/peri
mkdir -p $output_path && cd $output_path
# 2.1 注释 → 带 offset 的 CDS 片段
python /home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/scripts/calcu_peri/gtf_to_cds_with_offset.py \
  --annot "$GTF" --format gtf \
  --out-bed cds_with_offset.bed \
  --out-map tid_map.tsv

# 2.2 BAM → 1bp BED，并排序（两端都要排序）
bedtools bamtobed -i "$BAM" > psites.1bp.bed
sort -k1,1 -k2,2n -S 20% -T /tmp psites.1bp.bed -o psites.1bp.sorted.bed
sort -k1,1 -k2,2n -S 20% -T /tmp cds_with_offset.bed -o cds_with_offset.sorted.bed

# 2.3 交并（C++实现，极快），把每个命中的 p-site 投影为 CDS 相对坐标
bedtools intersect -s -wao -sorted \
  -a psites.1bp.sorted.bed \
  -b cds_with_offset.sorted.bed \
| awk 'BEGIN{OFS="\t"}
# A: psite BED6:  $1..$6
# B: cds  BED7:   $7..$13   (chr start end tid score strand cds_offset)
# 最后一列 $14 是 overlap；psite是1bp，命中时 overlap==1
($14==1){
  pos=$2; strand=$6;
  startB=$8; endB=$9; tid=$10; cds_off=$13;
  if(strand=="+"){ rel = cds_off + (pos - startB); }
  else            { rel = cds_off + (endB - 1 - pos); }
  if(rel<0) next;
  frame = rel % 3;
  codon = int(rel/3);
  print "F",tid,frame,1;        # 帧计数
  if(frame==0) print "C",tid,codon;  # frame0 的 codon 索引
}' > hits.stream.tsv
