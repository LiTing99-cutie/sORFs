mkdir -p S3
define_annotation_gencode_v41_human
ls $uniprot_prot_fasta
seqkit seq -n -M 150 $uniprot_prot_fasta |cut -f 1 -d " " > S3/uniprot.human.sep.id.txt

# uniprot中所有小于150aa的蛋白质
cd S3
define_annotation_gencode_v41_human
seqkit grep -f uniprot.human.sep.id.txt $uniprot_prot_fasta|seqkit seq -w 0 > uniprot.human.sep.fa
seqkit rmdup -s uniprot.human.sep.fa > uniprot.human.sep.rmdup.fa
# 4345
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
# 966/4345没有匹配
less out_fasta.txt |grep "No match"|wc -l
# 3243 有匹配
less out_fasta.txt |awk '$4==1 && $3==$5' | sort -u -k1,1 > perfect_match_to_ensId.uniq.txt
wc -l perfect_match_to_ensId.uniq.txt
# 3243
grep -f <(cut -f2 perfect_match_to_ensId.uniq.txt) $gpe_15 > uniprot.human.sep.15.gpe 

# 2024-06-05 得到uniprot蛋白的id和类型，如果是非isoform的蛋白，那么就得到PE蛋白质证据类型；如果是isoform蛋白，那么返回isoform
grep -f S3/uniprot.human.sep.id.txt $uniprot_prot_fasta >  S3/uniprot.human.sep.full.id.txt
awk '
{
    if ($0 ~ /Isoform/) {
        print "Isoform"
    } else {
        match($0, /PE=([0-9]+)/, arr)
        pe_value = (arr[1] != "" ? arr[1] : 1)
        print "PE_" pe_value
    }
}' S3/uniprot.human.sep.full.id.txt > S3/uniprot.human.sep.type.txt
paste -d '\t' S3/uniprot.human.sep.id.txt S3/uniprot.human.sep.type.txt > S3/uniprot.human.sep.id.type.txt