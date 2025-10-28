# get pub_non_HLA_non_Duffy.txt from /rd1/user/lit/project/sORFs/analysis/20250331_stat_human_ms_data/3.R
mkdir check_overlap_20250401
cd check_overlap_20250401
cut -f 3 ../output/all_sample_sep_unique_psm.txt|tail -n +2|sort|uniq > in_house_peptide.txt

awk '{print ">seq" NR "\n" $0}' in_house_peptide.txt > queries.fasta
awk '{print ">seq" NR "\n" $0}' pub_non_HLA_non_Duffy.txt > database.fasta

sed 's/I/L/g' queries.fasta > queries_modified.fasta
sed 's/I/L/g' database.fasta > database_modified.fasta

makeblastdb -in database_modified.fasta -dbtype prot -out ./database

blastp -query queries_modified.fasta -db ./database -task blastp-short -evalue 1000 -outfmt "6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore" > blastp.txt

# null
less blastp.txt | awk '$7==1 && $8==0 && $4==$6 && $4==$5' > exact_same_pep_allow_one_mismatch.txt
less blastp.txt | awk '$3==100 && $4==$6 && $4==$5' > exact_same_pep.txt
 
less blastp.txt | awk '$3==100 && $10==$4 && $11==1'
less blastp.txt | awk '$3==100 && $12==$5 && $9==1'

# null
less blastp.txt | awk '$3==100 && $10==$4 && $11==1 && $7==1 && $8==0'
less blastp.txt | awk '$3==100 && $12==$5 && $9==1 && $7==1 && $8==0'