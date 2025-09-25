seqkit tab2fx Public_SEP/iScience_20230.tab > Public_SEP/iScience_20230.fa
less /home/user/data3/lit/project/sORFs/03-Cross-anna/output/MS_sep/merged_final.txt | cut -f1,6|tail -n +2|seqkit tab2fx > in_house.list.fa

query=Public_SEP/iScience_20230.fa
subject=in_house.list.fa
blastp -query $query -subject $subject -out res.out -outfmt '6 qseqid sseqid pident qlen slen length bitscore evalue' -num_threads 20
less res.out | awk '$3==100 && $4==$6' |awk '$5==$6' > exact_same_spep.txt

seqkit tab2fx Public_SEP/iScience_20230.peptide.tab > Public_SEP/iScience_20230.peptide.fa
query=Public_SEP/iScience_20230.peptide.fa
subject=in_house.list.fa
blastp -query $query -subject $subject -out res.out -outfmt '6 qseqid sseqid pident qlen slen length bitscore evalue' -num_threads 20
less res.out | awk '$3==100 && $4==$6' |wc -l