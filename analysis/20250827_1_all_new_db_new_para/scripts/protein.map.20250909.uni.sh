sample=$1
out_dir=../processed/protein_map/$sample
mkdir -p $out_dir
cut -f1 ../processed/db_search_merge/$sample.peptide.tsv | tail -n +2 > $out_dir/peptides_il.txt
custom_db=/rd1/user/lit/project/sORFs/custom_database/human/custom_db_20250826_v2/human_brain_custom_db.fasta
pepMatch_jar=/rd1/user/lit/project/sORFs/analysis/archive/20250218_recheck_unique_peptide/PeptideMatchCMD_1.1.jar
[ -d custom_db_index ] || java -jar $pepMatch_jar -a index -d $custom_db -i custom_db_index
java -jar $pepMatch_jar -l -a query -i custom_db_index -Q $out_dir/peptides_il.txt -e -o $out_dir/out_fasta.txt &> $out_dir/out.match_n.txt

awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} NF>=5 && length($1)==($5-$4+1)' \
  "$out_dir/out_fasta.txt" > "$out_dir/out_fasta.len_eq.tsv"

# less "$out_dir/out_fasta.len_eq.tsv"|grep -v '^#'|cut -f1,2 > $out_dir/pep.orf.txt
less "$out_dir/out_fasta.len_eq.tsv"|grep -v '^#'|awk -v OFS='\t' '{print $1,$2,"'$sample'"}' > $out_dir/pep.orf.txt



