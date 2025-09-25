proj_path=/home/user/data3/lit/project/sORFs/

less ../MS_res_from_Galaxy/out_fasta.len_eq.tsv|grep -v '^#'|cut -f1,2 > ../processed/pick_orf/pep.orf.txt

python pick.orf.20250911.py \
  -i $proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.dup.txt \
  -e $proj_path/10-feature-egi/processed/feature_preprare/isoform.expr.txt \
  -o $proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/representative.tsv \
  --seed 42

awk -v OFS='\t' '{print $3,$2}' $proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/representative.tsv > \
    $proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/id.map.txt

id_map=$proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/id.map.txt
bash id.convert.20250911.sh ../processed/pick_orf/pep.orf.txt \
                   "$id_map" \
                   ../processed/pick_orf/pep.orf.rePicked.txt

python pep_to_protein_assign.20250915.alt.col.py \
  --pep2prot ../processed/pick_orf/pep.orf.rePicked.txt \
  --isoform /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt \
  --orf_psites /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt \
  --nca /home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-isoseq-gtf-20250909/expr/rpkm_N_C_A.txt \
  --rpf_metric RPF_codon_coverage --rpf_min 0.2 \
  --ps_metric Psites_codon_coverage --ps_min 0.1 \
  --min_c 0.2 \
  --outdir ../processed/pep_assign

# 1) 先做 ORF 序列与长度（第一部分命令）
# 生成 ORF_id、ORF_seq、ORF_length（三列，TSV）
FA_1=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.pep.fa
uniprot_fasta=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/human/uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.1.fasta
contam_fasta=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/contaminant_fasta/2022_JPR_contam.fasta
cat $FA_1 $uniprot_fasta $contam_fasta > ../processed/tmp/all.fasta
FA=../processed/tmp/all.fasta
OUT=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/results/custom_db_20250826
(echo -e "ORF_id\tORF_seq\tORF_length"; \
 seqkit fx2tab -l "$FA" \
 | awk -F'\t' 'BEGIN{OFS="\t"}{split($1,a,/ /); print a[1], $2, $4}') \
> "$OUT/orf_seq_len.tsv"
# 2) 生成两个表
python make_orf_tables.py \
  --cand_filtered ../processed/pep_assign/candidates.filtered.tsv \
  --orf_seq_len   "$OUT/orf_seq_len.tsv" \
  --isoform_info  /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt \
  --out_detail    ../processed/pep_assign/orf_detail.noCanonSharing.tsv \
  --out_orf       ../processed/pep_assign/orf_folded_counts.tsv
# 生成iBAQ
cat <(head -n1 ../MS_res_from_Galaxy/peptide_intensity_IL.merged.tsv) \
 <(grep 1_C8_T_T ../MS_res_from_Galaxy/peptide_intensity_IL.merged.tsv) > ../processed/quant/peptide_intensity_IL.tsv

# 计算ibaq
python compute_ibaq_len_norm.py \
  --file1 ../processed/pep_assign/assignments.unique.post.tsv \
  --file2 ../processed/quant/peptide_intensity_IL.tsv \
  --file3 ../processed/pep_assign/orf_folded_counts.tsv \
  --file3-cols all \
  -o ../processed/quant/ibaq_b_with_total.tsv

# 合并定量的信息
python merge_tables.py \
  --ibaq ../processed/quant/ibaq_b_with_total.tsv \
  --rpf  /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt \
  --iso  /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt \
  --rna  /home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-isoseq-gtf-20250909/expr/rpkm_N_C_A.txt \
  --gene_anno /home/user/data2/lit/project/ZNF271/data/annotation/Ensembl_106_Gencode_v41_Human_Transcript_stable_ID_version_Gene_stable_ID_version_Gene_name_Transcript_type_gene_type.txt \
  --out  ../results/ibaq_orf_rpf_iso_rna.tsv

python spearman_corr.py \
  --merged ../results/ibaq_orf_rpf_iso_rna.tsv \
  --out    ../results/spearman_results.tsv
# 如果使用isoform rpkm会不会相关性更高呢
#       var_x     var_y  spearman_rho     n
# 0         A  RPF_RPKM      0.574325  1939
# 1         A    iBAQ_B      0.302757  1559
# 2  RPF_RPKM    iBAQ_B      0.271657  1559  

# 单个样本看看时间
time bash run.all.20250912.sh 21pcw_1_C8_T_T

# 所有样本
python3 split_by_sample.py \
  --pep_orf_merged ../MS_res_from_Galaxy/pep.orf.merged.txt \
  --intensity_merged ../MS_res_from_Galaxy/peptide_intensity_IL.merged.tsv \
  --out_base ../processed/by_sample \
  --sample_list ../processed/by_sample/sample_list.txt

BASE=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/log
nohup cat ../processed/by_sample/sample_list.txt|grep _1_C8_T_T| xargs -P 8 -I{} bash -lc 'd="'"$BASE"'/{}/logs"; mkdir -p "$d";bash ./run.all.20250912.sh "{}">"$d/driver.log" 2>&1' &
