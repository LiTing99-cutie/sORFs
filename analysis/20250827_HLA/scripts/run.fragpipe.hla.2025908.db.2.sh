project_path=/rd1/user/lit/project/sORFs
# workflow=$project_path/fragpipe_config/HLA/HLA-DIA.workflow
workflow=$project_path/fragpipe_config/HLA/Nonspecific-HLA-diaPASEF.workflow
script=/rd1/user/lit/project/sORFs/Uni.Fragpipe.sh
# database_path=$project_path/custom_database/human/custom_db_20250826_v2/2025-08-26-decoys-human_brain_custom_db.fasta.fas
# database_path=$project_path/custom_database/human/uniprot/2025-05-02-decoys-contam-uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.fasta.fas
database_path=$project_path/custom_database/human/custom_db_ribo_fil_20250904/2025-09-04-decoys-human_brain_custom_db_ribo_fil.fasta.fas
file=/rd1/user/lit/project/sORFs/raw_data/HLA/MHC-I/XC03541MHC_P21_HLA1_Slot1-1_1_1100.d
output_path=../processed/HLA-I-db-2
mkdir -p $output_path && cd $output_path
mkdir -p config
sed 's|database.db-path=.*|database.db-path='"$database_path"'|' $workflow > config/fragpipe.workflow
echo -e "$file\t$(basename $file)\t\tDIA" > config/fragpipe-files.fp-manifest
bash $script config/fragpipe.workflow config/fragpipe-files.fp-manifest ./ 40