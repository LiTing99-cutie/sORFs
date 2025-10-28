project_path=/rd1/user/lit/project/sORFs
# workflow=$project_path/fragpipe_config/HLA/HLA-DIA.workflow
workflow=$project_path/fragpipe_config/HLA/Nonspecific-HLA-II-diaPASE.workflow
script=/rd1/user/lit/project/sORFs/Uni.Fragpipe.sh
# database_path=$project_path/custom_database/human/custom_db_20250826_v2/2025-08-26-decoys-human_brain_custom_db.fasta.fas
# database_path=$project_path/custom_database/human/uniprot/2025-05-02-decoys-contam-uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.fasta.fas
# database_path=$project_path/custom_database/human/custom_db_ribo_fil_20250904/2025-09-04-decoys-human_brain_custom_db_ribo_fil.fasta.fas
# database_path=$project_path/custom_database/human/custom_db_fil/c_1_rpf_1_psites_0.6/2025-09-17-decoys-candidateORF.filtered.addContam.fa.fas
database_path=$project_path/custom_database/human/custom_db_fil/c_2_rpf_2_psites_1/2025-09-17-decoys-candidateORF.rmdup.pep.addContam.fa.fas
# file=/rd1/user/lit/project/sORFs/raw_data/HLA/MHC-II/Replicate_1/XC03543MHC_P21_HLA2_Slot1-1_1_1102.d
# file=/rd1/user/lit/project/sORFs/raw_data/HLA/MHC-I/XC03541MHC_P21_HLA1_Slot1-1_1_1100.d
file=/rd1/user/lit/project/sORFs/raw_data/HLA/MHC-II/Replicate_2/XC04813MHC_P21_HLA2_Slot1-1_1_18437.d

output_path=../processed/HLA-II-db-2
mkdir -p $output_path && cd $output_path
mkdir -p config
sed 's|database.db-path=.*|database.db-path='"$database_path"'|' $workflow > config/fragpipe.workflow
echo -e "$file\t$(basename $file)\t\tDIA" > config/fragpipe-files.fp-manifest
bash $script config/fragpipe.workflow config/fragpipe-files.fp-manifest ./ 40