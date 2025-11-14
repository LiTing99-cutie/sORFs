source activate base
candidateORF=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.pep.fa
# 从中筛选出frame0_fraction大于等于0.5 & Psites_codon_coverage大于等于0.1的ORF；此外，去掉那些和经典CDS重叠的非经典ORF
frame_stats=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/peri/psite_frame_stats.v2.tsv
rpf_psite=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt
overlap=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf_overlap_inframe.txt

output_path=../processed/annotation/RibORF_annot/candidate_ORFs/filtered_v3_20251103
mkdir -p $output_path
python3 filter_orf_fasta.v2.py \
  --candidate-fa /home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.pep.fa \
  --psite-stats /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/peri/psite_frame_stats.v2.tsv \
  --rpf-psite  /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt \
  --overlap    /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf_overlap_inframe.txt \
  --out-fa     $output_path/filtered_candidates.fa \
  --f0-thr 0 --cov-thr 0 &> ../log/filter_orf_fasta.filtered_v3_20251103.log

seqkit rmdup -s -t protein -w 0 -j 20 \
    -D $output_path/candidateORF.filtered.M.dup.txt \
    -o $output_path/candidateORF.filtered.rmdup.pep.fa \
    $output_path/filtered_candidates.fa

python pick_representative.20250910.py \
  -i $output_path/candidateORF.filtered.M.dup.txt \
  -e /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.txt \
  -o $output_path/representative.tsv \
  --seed 42

python3 replace_ids.py \
  --map $output_path/representative.tsv \
  --fasta-in ${output_path}/candidateORF.filtered.rmdup.pep.fa \
  --fasta-out ${output_path}/candidateORF.filtered.rmdup.renamed.pep.fa

echo "--- Count ORF numbers per ORF type --- "
awk '/^>/{ 
    h = substr($0, 2);            # 去掉 >
    n = split(h, a, "\\|");      # 以 | 分割
    t = (n == 5 ? a[n-1] : "Uniprot Reviewed");
    counts[t]++
}
END {
    for (k in counts) print k "\t" counts[k]
}' ${output_path}/candidateORF.filtered.rmdup.renamed.pep.fa | sort -k2,2nr \
    > ${output_path}/candidateORF.filtered.rmdup.renamed.pep.orf_type.txt

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

# 使用shell增加列，但是不统计
# awk -F'\t' -v OFS='\t' '
# # ------- 第1个文件：orfs_3_ways，构建命中表 m[key,software] -------
# FNR==NR{
#   sw = $3; gsub(/^[ \t]+|[ \t]+$/, "", sw)
#   key = $4
#   # 统一软件名
#   tl = tolower(sw)
#   if (tl=="price") m[key,"PRICE"]=1
#   else if (tl=="ribocode") m[key,"RiboCode"]=1
#   else if (tl=="riborf") m[key,"RibORF"]=1
#   next
# }
# # ------- 第2个文件：原结果，加三列 -------
# FNR==1{
#   # 按表头定位 ORF_id_custom 列
#   for(i=1;i<=NF;i++) h[$i]=i
#   if(!("ORF_id_custom" in h)){
#     print "[ERR] 缺少列 ORF_id_custom" > "/dev/stderr"; exit 1
#   }
#   c = h["ORF_id_custom"]
#   print $0, "Is_PRICE", "Is_RiboCode", "Is_RibORF"
#   next
# }
# {
#   key = $(c)
#   print $0, (m[key,"PRICE"]?"TRUE":"FALSE"),
#              (m[key,"RiboCode"]?"TRUE":"FALSE"),
#              (m[key,"RibORF"]?"TRUE":"FALSE")
# }
# ' "$orfs_3_ways" "$in_with_custom" > "$out_with_flags"
