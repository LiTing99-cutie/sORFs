output_path=../processed/annotation/RibORF_annot/candidate_ORFs/filtered_v2_20251021
# 添加额外信息
seqkit fx2tab ${output_path}/candidateORF.filtered.rmdup.renamed.pep.fa |cut -f 1 > ${output_path}/candidateORF.filtered.rmdup.renamed.pep.txt
sed -i '1i ORF_id' ${output_path}/candidateORF.filtered.rmdup.renamed.pep.txt
source /home/user/data3/lit/project/sORFs/config.sh
source activate biotools
python extend_id/orf_extend.py \
  --orf-list ${output_path}/candidateORF.filtered.rmdup.renamed.pep.txt \
  --orf-col ORF_id \
  --orf-seq-len $orf_seq_len \
  --rpf $RPF \
  --iso $ISO \
  --rna $RNA \
  --gene-anno $GENE_ANNO \
  --out-extended ${output_path}/candidateORF.filtered.rmdup.renamed.pep.extended.txt
# 增加ORF_start，ORF_end以及CDS_seq列
script=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/scripts/augment_orf_table.sh
bash $script ${output_path}/candidateORF.filtered.rmdup.renamed.pep.extended.txt ${output_path}/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt
# 增加ORF_id_custom列
awk -F'\t' -v OFS='\t' '
NR==1{
  for(i=1;i<=NF;i++) h[$i]=i
  need="Strand Chr ORF_start ORF_end ORF_seq"
  n=split(need,a," ")
  for(i=1;i<=n;i++) if(!(a[i] in h)){print "缺少列: " a[i] > "/dev/stderr"; exit 1}
  print $0, "ORF_id_custom"
  next
}
{
  id = $(h["Strand"]) $(h["Chr"]) ":" $(h["ORF_start"]) "-" $(h["ORF_end"]) ":" $(h["ORF_seq"])
  print $0, id
}' ${output_path}/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/augmented.tsv > ${output_path}/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/augmented.1.tsv
# 增加软件鉴定信息列
orfs_3_ways=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20251022_custom_fa_gtf/processed/orfs/merged/orfs.3_ways.txt
in_with_custom=${output_path}/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/augmented.1.tsv            # 上一步生成，已包含 ORF_id_custom 的表
out_with_flags=${output_path}/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/augmented.2.tsv
python extend_id/merge_orf_flags_and_stats_optimized.py \
  --orfs-3-ways $orfs_3_ways \
  --in-with-custom $in_with_custom \
  --out-with-flags $out_with_flags \
  --out-stats ${output_path}/software.stat.txt \
  -v --tqdm --chunksize 2000000