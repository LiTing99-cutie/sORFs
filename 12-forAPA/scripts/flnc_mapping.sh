source activate biotools
flnc_bam=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/processed/p21-IsoSeq.flnc.bam
threads=40
fa=/home/user/data/lit/database/public/genome/hg38/hg38.fa
outdir=../processed
script_dir=/home/user/data3/lit/project/sORFs/12-forAPA/scripts/pa_iden_script
# 将flnc.bam比对到基因组上
echo "*** bam2fastq 3.5.0"
time bam2fastq -o ${outdir}/flnc -u $flnc_bam
echo "*** minimap2 2.15-r905 map to genome"
time minimap2 -t $threads -ax splice -uf --MD $fa ${outdir}/flnc.fastq | \
     samtools sort -@ $threads -o ${outdir}/Aligned.out.sorted.bam
[ -f ${outdir}/flnc.fastq ] && rm -rf ${outdir}/flnc.fastq
echo "*** Index BAM"
samtools index -@ $threads ${outdir}/Aligned.out.sorted.bam
# PA位点鉴定
# time samtools view -h ${outdir}/Aligned.out.sorted.bam| \
#     $script_dir/samAddTag.pl --coverage --identity --unmapped unmapped.sam 2>lengthInconsistent.sam | \
#     samtools view -buS - |samtools sort -m 5G -@ 40 -o mapped.sorted.bam - && $script_dir/sam2bed.pl -s -t CV,ID mapped.sorted.bam >mapped.bed12+
# time bash $script_dir/PAratio.OnGenes.modi.1.sh -g  $genePredExt -r mapped.bed12+ -o ./ -f $fa
# awk -v OFS="\t" '{split($4,a,",");split($5,readN,",");split($6,usage,",");for(i=1;i<=length(a);i++){print $1,a[i]-1,a[i],$3,readN[i],$2,usage[i]}}' PA.onGene.tsv > PAusage.bed6+