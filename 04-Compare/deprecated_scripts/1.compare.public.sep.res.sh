output_path=$PWD/output/S1
database_fa=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/group_specific_db/uniprot.contam.sorfs.trans.modiHeader.fa
jar=/home/user/data2/lit/bin/PeptideMatchCMD_1.1.jar
# pub_non_HLA_peptide.txt 由1.compare.public.sep.res.Rmd生成
mkdir -p $output_path && cd $output_path
java -jar $jar -a index -d $database_fa -i database_fa 
java -jar $jar -a query -i database_fa -Q pub_non_HLA_peptide.txt -l -e -o out_list.txt  &> match_n.txt
less out_list.txt |grep "No match"|wc -l
less out_list.txt |grep "No match" > no_match_peptide.txt

