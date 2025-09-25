conda activate base

# 生成随机选择的ID以及代表性ID的对照表格
python pick.orf.20250911.py \
  -i $proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.dup.txt \
  -e $proj_path/10-feature-egi/processed/feature_preprare/isoform.expr.txt \
  -o $proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/representative.tsv \
  --seed 42

awk -v OFS='\t' '{print $3,$2}' $proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/representative.tsv > \
    $proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/id.map.txt

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