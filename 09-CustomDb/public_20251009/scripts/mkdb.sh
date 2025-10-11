# 整合污染
contam_fasta=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/contaminant_fasta/2022_JPR_contam.fasta
candidate_ORFs=/home/user/data3/rbase/translation_pred/models/lib/ORF/candidate_ORFs/candidateORF.6aa.long.M.rmdup.pep.fa
mkdir -p ../results/custom_db_20251009/
cat $candidate_ORFs $contam_fasta > ../results/custom_db_20251009/public_db.fasta
