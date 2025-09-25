seqkit seq -w 0 uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.fasta > uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.1.fasta
seqkit grep -f output/anno_sep.id.txt uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.1.fasta > output/anno_sep.id.fasta

pep_match(){
    jar=/rd1/user/lit/project/sORFs/analysis/20250218_recheck_unique_peptide/PeptideMatchCMD_1.1.jar
    query=$1
    db=$2
    [ -d db ] || java -jar $jar -a index -d $db -i db
    java -jar $jar -a query -i db -Q $query -o out_fasta.txt 
}
pep_match output/anno_sep.id.fasta ./Homo_sapiens.GRCh38.pep.all.1.fa 
less out_fasta.txt |grep "no match"|wc -l
less out_fasta.txt |awk '$4==1 && $3==$5' > perfect_match.txt
# 有106个只有一个匹配
less perfect_match.txt  |cut -f 1 |sort |uniq -c|awk '$1==1'|wc -l

join -1 2 -2 6 <(sort -k2,2 perfect_match.txt) <(sort -k6,6 Ensembl_106_id_map.txt) > perfect_match.add.trans.gene.id.txt


less /rd1/user/lit/project/sORFs/analysis/20250331_stat_human_ms_data/output/all_sample_sep_unique_psm.txt|cut -f 3|tail -n +2|sort |uniq  > in_house_pep.txt
pep_match(){
    # 输入，每行一个肽段
    # I/L equivalent
    jar=/rd1/user/lit/project/sORFs/analysis/20250218_recheck_unique_peptide/PeptideMatchCMD_1.1.jar
    query=$1
    db=$2
    db_name=$3
    output=$4
    [ -d $db_name ] || java -jar $jar -a index -d $db -i $db_name
    java -jar $jar -a query -i $db_name -Q $query -l -e -o $output
}

db=/rd1/user/lit/project/sORFs/custom_database/human/trans_based_database_human_20250326/uniprot.contam.sorfs.trans.modiHeader.fa
pep_match in_house_pep.txt $db uniprot_human_add_sorfs tmp_match/out_list.txt &> tmp_match/match_n.txt