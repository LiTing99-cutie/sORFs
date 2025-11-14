python3 hla/build_unique_pep_orf_table.py \
  --hla /home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/MS_res_from_Galaxy/HLA/MHC-I-peptide.tsv \
  --rna /home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-custom-fa-gtf-20251022/expr/rpkm_N_C_A.txt \
  --iso /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt \
  --rpf /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt \
  --orf-seq-len /home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/results/custom_db_20250826/orf_seq_len.tsv \
  --gene-anno /home/user/data2/lit/project/ZNF271/data/annotation/Ensembl_106_Gencode_v41_Human_Transcript_stable_ID_version_Gene_stable_ID_version_Gene_name_Transcript_type_gene_type.txt \
  --out-map ../processed/hla/orf_unique_counts-I.tsv \
  --out-extended ../processed/hla/orf_unique_extended-I.tsv

python3 hla/build_unique_pep_orf_table.py \
  --hla /home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/MS_res_from_Galaxy/HLA/MHC-II-peptide.tsv \
  --rna /home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-custom-fa-gtf-20251022/expr/rpkm_N_C_A.txt \
  --iso /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt \
  --rpf /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt \
  --orf-seq-len /home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/results/custom_db_20250826/orf_seq_len.tsv \
  --gene-anno /home/user/data2/lit/project/ZNF271/data/annotation/Ensembl_106_Gencode_v41_Human_Transcript_stable_ID_version_Gene_stable_ID_version_Gene_name_Transcript_type_gene_type.txt \
  --out-map ../processed/hla/orf_unique_counts-II.tsv \
  --out-extended ../processed/hla/orf_unique_extended-II.tsv