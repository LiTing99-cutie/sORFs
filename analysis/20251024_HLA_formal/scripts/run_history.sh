# HLA-I型数据；肽段大小为8-12aa
# HLA-II型数据；肽段大小为8-25aa
project_path=/rd1/user/lit/project/sORFs
database_path=$project_path/custom_database/human/custom_db_20251024_ribo_filtered/2025-10-24-decoys-contam-candidateORF.filtered.rmdup.renamed.addContam.pep.fa.fas
nohup bash run.fragpipe.hla.1.2.20251024.sh $database_path &> ../log/run.fragpipe.hla.1.2.20251024.log &
