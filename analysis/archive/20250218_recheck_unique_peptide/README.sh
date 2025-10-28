# 检查unique peptide中有哪些是真正的unique peptide
psm_path=/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/Formal/mouse/2025_02_06_default_ribo_database/Trypsin/CAD20241224licq_BSEP_DDA_60min_E16_1_3_30K_T_T_Slot2_4_1_4673_d/psm.tsv
# R脚本中提取小肽对应的unique的psm
# 构建数据库
fa=/rd1/user/lit/project/sORFs/custom_database/Ribo_ORFs_add_assemble_20250125/uniprot.nonCano.sorf.20250125.fa
makeblastdb -in $fa -dbtype prot -out ./database

seqkit tab2fx sep_unique_psm_peptide.txt > sep_unique_psm_peptide.fa
seqkit rmdup -s sep_unique_psm_peptide.fa > sep_unique_psm_peptide.redup.fa
blastp -query sep_unique_psm_peptide.redup.fa -db ./database -outfmt '6 qseqid sseqid pident qlen slen length bitscore evalue' -num_threads 30 -out unique.psm.p.blastp.txt
less unique.psm.p.blastp.txt | awk '$3==100 && $4==$6' > exact_same_pep.txt

grep -A 1 ENSMUST00000199164.2+chr3:58322390-58322504 sep_unique_psm_peptide.redup.fa > tmp.fa
blastp -query tmp.fa -db ./database -task blastp-short -evalue 100 -outfmt '6 qseqid sseqid pident qlen slen length bitscore evalue' > tmp.blastp.txt
blastp -query tmp.fa -subject $fa -task blastp-short -evalue 100 -outfmt '6 qseqid sseqid pident qlen slen length bitscore evalue' > tmp.blastp.txt

blastp -query tmp.fa -db ./database \
  -task blastp-short \
  -word_size 2 \
  -matrix PAM30 \
  -evalue 100 \
  -max_target_seqs 1000 \
  -max_hsps 1 \
  -gapopen 9 -gapextend 1 \
  -outfmt "6 qseqid sseqid pident length mismatch evalue" > tmp.blastp.txt

  java -jar PeptideMatchCMD_1.1.jar -h
  java -jar PeptideMatchCMD_1.1.jar -a index -d uniprot_sprot.fasta -i sprot_index 
  java -jar PeptideMatchCMD_1.1.jar -a query -i sprot_index -q AAFGGSGGR -o out.txt 
  java -jar PeptideMatchCMD_1.1.jar -a query -i sprot_index -q AAFGGSGGR,GVPDIR -o out.txt 
  java -jar PeptideMatchCMD_1.1.jar -a query -i sprot_index -q AAFGGSGGR,GVPDIR -e -o out.txt 
  java -jar PeptideMatchCMD_1.1.jar -a query -i sprot_index -Q query.fasta -e -o out_fasta.txt 
  java -jar PeptideMatchCMD_1.1.jar -a query -i sprot_index -Q query.list -l -e -o out_list.txt 

  java -jar PeptideMatchCMD_1.1.jar -a index -d $fa -i sprot_index 
  java -jar PeptideMatchCMD_1.1.jar -a query -i sprot_index -q RPGPPPPR -e -o out.txt 
  java -jar PeptideMatchCMD_1.1.jar -a query -i sprot_index -Q sep_unique_psm_peptide.redup.fa -e -o out_fasta.txt &> out.match_n.txt