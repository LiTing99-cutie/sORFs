db_path=/rd1/user/lit/project/sORFs/custom_database/human/trans_based_database_human_20250326/2025-03-26-decoys-uniprot.contam.sorfs.trans.modiHeader.fa.fas
seqkit grep -f stat_output/S7/protein.id.20250610.txt $db_path > stat_output/S7/protein.id.20250610.fa
seqkit fx2tab stat_output/S7/protein.id.20250610.fa | awk -v OFS='\t' '{print $1,$NF}' > stat_output/S7/protein.id.20250610.tab