mkdir -p S3
define_annotation_gencode_v41_human
ls $uniprot_prot_fasta
seqkit seq -n -M 150 $uniprot_prot_fasta |cut -f 1 -d " " > S3/uniprot.human.sep.id.txt

Rscript sORFs/02-Mass-spec/human/3.inspect.anno_sep.R

###### 实现ID的对应 ######
# 提取detected undetected annotated protein
seqkit grep -f <(cut -f 1 /home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/anno_sep_ms_info.txt|tail -n +2) \
 $uniprot_prot_fasta > S3/anno.detected.fa
grep -v -f <(cut -f 1 /home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/anno_sep_ms_info.txt|tail -n +2) \
S3/uniprot.human.sep.id.txt > S3/uniprot.undetected.sep.id.txt
seqkit grep -f S3/uniprot.undetected.sep.id.txt $uniprot_prot_fasta > S3/anno.undetected.fa
seqkit seq -s S3/anno.undetected.fa > S3/anno.undetected.seq
seqkit seq -s S3/anno.detected.fa > S3/anno.detected.seq

# uniprot中所有小于150aa的蛋白质
cd S3
define_annotation_gencode_v41_human
seqkit grep -f uniprot.human.sep.id.txt $uniprot_prot_fasta|seqkit seq -w 0 > uniprot.human.sep.fa
seqkit fx2tab uniprot.human.sep.fa|awk -v OFS='\t' '{print $1, $NF}' > uniprot.human.sep.tab
pep_match(){
    jar=/home/user/data2/lit/bin/PeptideMatchCMD_1.1.jar
    query=$1
    db=$2
    [ -d db ] || java -jar $jar -a index -d $db -i db
    java -jar $jar -a query -i db -Q $query -o out_fasta.txt 
}
pep_match uniprot.human.sep.fa $ensembl_prot_fasta 
# 913/4376没有匹配
less out_fasta.txt |grep "No match"|wc -l
less out_fasta.txt |awk '$4==1 && $3==$5' > perfect_match.txt
# 3300/4376有完美的匹配
sort -u -k1,1 perfect_match.txt > perfect_match_uniq.txt
# # 2230/4376有完美且唯一的匹配
# less perfect_match.txt  |cut -f 1 |sort |uniq -c|awk '$1==1'|wc -l

# 看是否没有isoform的小于150aa的匹配上的比例更加高
uniprot_fa_no_iso=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/human/uniprotkb_taxonomy_id_9606_AND_reviewed_2025_04_26.fasta.gz
seqkit seq -M 150 $uniprot_fa_no_iso > uniprot.human.sep.no_iso.fa
pep_match uniprot.human.sep.no_iso.fa $ensembl_prot_fasta 
# 1827
less out_fasta.txt |awk '$4==1 && $3==$5'|sort -u -k1,1 |wc -l

grep -f <(cut -f2 perfect_match_uniq.txt) $gtf > uniprot.human.sep.gtf
gtfToGenePred -genePredExt -geneNameAsName2 uniprot.human.sep.gtf uniprot.human.sep.15.gpe 

awk -F '\t' '$3 == "transcript" {
    match($9, /transcript_id "([^"]+)"/, tid)
    match($9, /protein_id "([^"]+)"/, pid)
    print tid[1]"\t"pid[1]
}' uniprot.human.sep.gtf > uniprot.human.sep.trans_id.pro_id.txt

###### 整合annotated SEP以及 novel SEP的信息######
# /home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/anno_sep_ms_info.txt
# /home/user/data3/lit/project/sORFs/02-Mass-spec/human/ms_info_updated.txt


