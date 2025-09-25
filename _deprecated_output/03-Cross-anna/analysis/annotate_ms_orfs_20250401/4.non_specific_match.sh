
output_path=/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/S4
ensembl_fa=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/Ensembl-106/Homo_sapiens.GRCh38.pep.all.1.fa
cd $output_path
less /home/user/data3/lit/project/sORFs/02-Mass-spec/human/sep_unique_psm.txt |cut -f3 |tail -n +2|sort|uniq > in_house_peptide.txt
grep -f in_house_peptide.txt $ensembl_fa

jar=/home/user/data2/lit/bin/PeptideMatchCMD_1.1.jar
java -jar $jar -a index -d $ensembl_fa -i ensembl_fa 
java -jar $jar -a query -i ensembl_fa -Q in_house_peptide.txt -l -e -o out_list.txt 
# 237
less out_list.txt |grep "No match"|wc -l
less out_list.txt |grep "No match" > no_match_peptide.txt
java -jar $jar -a query -i ensembl_fa -Q in_house_peptide.txt -l -o out_list.txt 
# 242
less out_list.txt |grep "No match"|wc -l
