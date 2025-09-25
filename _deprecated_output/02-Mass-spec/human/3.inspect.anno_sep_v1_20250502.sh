mkdir -p S3
define_annotation_gencode_v41_human
ls $uniprot_prot_fasta
seqkit seq -n -M 150 $uniprot_prot_fasta |cut -f 1 -d " " > S3/uniprot.human.sep.id.txt

# uniprot中所有小于150aa的蛋白质
cd S3
define_annotation_gencode_v41_human
seqkit grep -f uniprot.human.sep.id.txt $uniprot_prot_fasta|seqkit seq -w 0 > uniprot.human.sep.fa
seqkit rmdup -s uniprot.human.sep.fa > uniprot.human.sep.rmdup.fa
seqkit fx2tab uniprot.human.sep.rmdup.fa|awk -v OFS='\t' '{print $1, $NF}' > uniprot.human.sep.rmdup.tab

# 将uniprot蛋白质id和ensembl transcript id对应起来
bash /home/user/data3/lit/project/sORFs/01-ribo-seq/S3.0c.Uni.translate_gtf.v2.20250325.sh $gtf $fa
pep_match(){
    jar=/home/user/data2/lit/bin/PeptideMatchCMD_1.1.jar
    query=$1
    db=$2
    db_name=$3
    [ -d $db_name ] || java -jar $jar -a index -d $db -i $db_name
    java -jar $jar -a query -i $db_name -Q $query -o out_fasta.txt 
}
pep_match uniprot.human.sep.rmdup.fa prot.fa ensembl.prot
# 966/4376没有匹配
less out_fasta.txt |grep "No match"|wc -l
# 3268 有匹配
less out_fasta.txt |awk '$4==1 && $3==$5' | sort -u -k1,1 > perfect_match_to_ensId.uniq.txt
wc -l perfect_match_to_ensId.uniq.txt

# 看是否没有isoform的小于150aa的匹配上的比例更加高【并没有】
uniprot_fa_no_iso=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/human/uniprotkb_taxonomy_id_9606_AND_reviewed_2025_04_26.fasta.gz
seqkit seq -M 150 $uniprot_fa_no_iso > uniprot.human.sep.no_iso.fa
seqkit rmdup -s uniprot.human.sep.no_iso.fa > uniprot.human.sep.no_iso.rmdup.fa
pep_match uniprot.human.sep.no_iso.rmdup.fa $ensembl_prot_fasta ensembl.prot
# 1797
less out_fasta.txt |awk '$4==1 && $3==$5'|sort -u -k1,1 |wc -l

grep -f <(cut -f2 perfect_match_to_ensId.uniq.txt) $gpe_15 > uniprot.human.sep.15.gpe 



