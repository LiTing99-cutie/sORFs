cp /rd1/user/lit/project/sORFs/S1* ./

cd sorfs_id
db_path=/rd1/user/lit/project/sORFs/custom_database
seqkit seq -n $db_path/Transcriptome_based_20250209/uniprot.nonCano.sorf.20250209.trans_based.fa | egrep "^ENS|^STRG" > trans_based.sorf.id.txt
seqkit seq -n $db_path/Ribo_ORFs_add_assemble_20250125/uniprot.nonCano.sorf.20250125.fa | egrep "^ENS|^STRG" > ribo_based.sorf.id.txt