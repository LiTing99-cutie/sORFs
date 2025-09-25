#Convert gtf to cds.fa and prot.fa
#Usage translate_gtf.sh *.gtf *.fa [genome fasta]; output cds.fa prot.fa in the current directory 

gtf=$1
genome_fasta=$2
output_dir=$3
mkdir -p $output_dir
cd $output_dir
# cds gtf
grep -P "\tCDS" $gtf > cds.gtf
# cds bed (1-based)
awk 'BEGIN{OFS=FS="\t"}{split($9,A,"\"");print $1,$4,$5,A[4],".",$7}' cds.gtf > cds.1-based.bed6
# cds bed (0-based)
awk 'BEGIN{OFS=FS="\t"}{print $1,$2-1,$3,$4,$5,$6}' cds.1-based.bed6 > cds.bed6
# note that for reverse strand; if cds regions are in diffrent exons; the larger coordinate should be put before
awk -v OFS=FS="\t" '$6=="+"' cds.bed6 > cds.f.bed6
awk -v OFS=FS="\t" '$6=="-"' cds.bed6 |sort -k4,4 -k2,2nr > cds.r.bed6
# cds fasta
bedtools getfasta -s -name -fi $genome_fasta -bed ./cds.f.bed6 > ./cds.chopped.f.fa
bedtools getfasta -s -name -fi $genome_fasta -bed ./cds.r.bed6 > ./cds.chopped.r.fa
cat cds.chopped.f.fa cds.chopped.r.fa > cds.chopped.fa
# connect chopped cds fasta sequence
awk '/^>/ { split($0,tA,"::"); if (current_id != tA[1]) { current_id = tA[1]; print seq; seq=""; print $0;}} 
	!/^>/ { seq = seq $0 } END { print seq }' ./cds.chopped.fa > ./cds.fa
sed -i '1d; s/::chr.*//g' ./cds.fa
# Translate DNA .fa file to peptide .fa
# 20250325修改
faTrans -stop ./cds.fa ./prot.multiLine.fa
awk '/^>/ {if (NR!=1) {printf("\n")}; printf("%s\n",$0); next} {printf("%s",$0)} END {printf("\n")}' ./prot.multiLine.fa > ./prot.fa
# clean
mkdir -p translate_tmp
mv cds.gtf cds.1-based.bed6 cds.bed6 cds.f.bed6 cds.r.bed6 cds.chopped.f.fa cds.chopped.r.fa cds.chopped.fa prot.multiLine.fa translate_tmp
rm -rf translate_tmp