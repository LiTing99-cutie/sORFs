db_path=/rd1/user/lit/project/sORFs/custom_database/human/trans_based_database_human_20250326/2025-03-26-decoys-uniprot.contam.sorfs.trans.modiHeader.fa.fas
output_path=../results/S5/
seqkit grep -f $output_path/protein.id.txt $db_path > $output_path/protein.id.fa
seqkit fx2tab $output_path/protein.id.fa | awk -v OFS='\t' '{print $1,$NF}' > $output_path/protein.id.tab