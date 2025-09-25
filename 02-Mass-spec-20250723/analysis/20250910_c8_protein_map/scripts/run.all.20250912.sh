sample=$1
proj_path=/home/user/data3/lit/project/sORFs/
id_map=$proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/id.map.txt
orf_seq_len=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/results/custom_db_20250826/orf_seq_len.tsv
analysis_path=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map
output_dir=$analysis_path/processed/by_sample/$sample
mkdir -p $output_dir/pick_orf $output_dir/pep_assign $output_dir/quant 
pep_orf=$output_dir/pep.orf.txt
peptide_intensity_IL=$output_dir/peptide_intensity_IL.tsv

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start id.convert.20250911.sh"
bash id.convert.20250911.sh $pep_orf \
                   "$id_map" \
                   $output_dir/pick_orf/pep.orf.rePicked.txt

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start pep_to_protein_assign.20250915.alt.col.py"
python pep_to_protein_assign.20250915.alt.col.py \
  --pep2prot $output_dir/pick_orf/pep.orf.rePicked.txt \
  --isoform /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt \
  --orf_psites /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt \
  --nca /home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-isoseq-gtf-20250909/expr/rpkm_N_C_A.txt \
  --rpf_metric RPF_codon_coverage --rpf_min 0.2 \
  --ps_metric Psites_codon_coverage --ps_min 0.1 \
  --min_c 0.2 \
  --outdir $output_dir/pep_assign

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start make_orf_tables.py"
python make_orf_tables.py \
  --cand_filtered $output_dir/pep_assign/candidates.filtered.tsv \
  --orf_seq_len   "/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/results/custom_db_20250826/orf_seq_len.tsv" \
  --isoform_info  /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt \
  --out_detail    $output_dir/pep_assign/orf_detail.noCanonSharing.tsv \
  --out_orf       $output_dir/pep_assign/orf_folded_counts.tsv

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start compute_ibaq_len_norm.py"
python compute_ibaq_len_norm.py \
  --file1 $output_dir/pep_assign/assignments.unique.post.tsv \
  --file2 $peptide_intensity_IL \
  --file3 $output_dir/pep_assign/orf_folded_counts.tsv \
  --file3-cols all \
  -o $output_dir/quant/ibaq_b_with_total.tsv



