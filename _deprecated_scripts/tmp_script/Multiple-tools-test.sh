#!/usr/bin/sh

################################################
#File Name: Multiple-tools-test.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 15 Oct 2024 06:29:34 PM CST
################################################

set -eo pipefail

##### RibORF #####

RibORF_path=/home/user/data2/lit/software/RibORF/RibORF.2.0

# project path
project_path=/home/user/data3/lit/project/sORFs/01-ribo-seq

# 基因组及其注释
fa=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
genePred=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.genePred.txt
sample_output_path=$project_path/mouse_brain_output_20241011/2015_Science_ChoJ_hippocampus/SRR2163083/output
log_path=$project_path/mouse_brain_output_20241011/2015_Science_ChoJ_hippocampus/SRR2163083/log
[ -d $sample_output_path/Ribo-ORFs/RibORF ] || mkdir $sample_output_path/Ribo-ORFs/RibORF

# 1. Annotation (每个基因组运行一次)
mkdir -p annot/RiboORF/mm39
perl $RibORF_path/ORFannotate.pl -g $fa -t $genePred -o annot/RiboORF/mm39
candidateORF=$project_path/annot/RiboORF/mm39/candidateORF.genepred.txt

# 2. bam2sam
cd $sample_output_path
bam=alignment/SRR2163083_Aligned.sortedByCoord.out.bam
samtools view $bam -h -@ 40 > alignment/SRR2163083.sam
sam=alignment/SRR2163083.sam

# 3. output
# time 59min
cd Ribo-ORFs/RibORF
time perl $RibORF_path/readDist.pl -f ../../alignment/SRR2163083.sam -g $genePred -o ./ -d 25,26,27,28,29,30,31,32,33,34, -l 40 -r 70 &> log/RibORF.readDist.1.log
perl $RibORF_path/parameterOffset.pl ./sta.read.dist.25,26,27,28,29,30,31,32,33,34,.txt
less offset.corretion.parameters.txt | awk -v OFS='\t' '{print $1,$6}' > offset.corretion.parameters.2.txt
# 6min
time perl $RibORF_path/offsetCorrect.pl -r ../../alignment/SRR2163083.sam -p offset.corretion.parameters.2.txt -o corrected.SRR2163083.mapping.sam &> $log_path/offsetCorrect.log
# 6min
time perl $RibORF_path/readDist.pl -f *.sam -g $genePred -o ./ -d 1
# several hours
time perl $RibORF_path/ribORF.pl -f *.sam -c $candidateORF -o ./

# test do not use custom offset files
mkdir test && cd test
perl $RibORF_path/parameterOffset.pl ../sta.read.dist.25,26,27,28,29,30,31,32,33,34,.txt
less offset.corretion.parameters.txt | awk -v OFS='\t' '{print $1,$6}' > offset.corretion.parameters.2.txt
# 6min
time perl $RibORF_path/offsetCorrect.pl -r ../../../alignment/SRR2163083.sam -p offset.corretion.parameters.2.txt -o corrected.SRR2163083.mapping.sam &> $log_path/offsetCorrect.log
# 6min
time perl $RibORF_path/readDist.pl -f *.sam -g $genePred -o ./ -d 1
# several hours
time perl $RibORF_path/ribORF.pl -f *.sam -c $candidateORF -o ./

compare_script=/home/user/data3/lit/project/sORFs/01-ribo-seq/ref/Ribo-seq-Tool-Comparison-Scripts-v2.0/
# 对repre.valid.pred.pvalue.parameters.txt进行过滤
python $compare_script/Scripts\ for\ RibORF1.0\ Analysis/ReFormat_genepred.py repre.valid.pred.pvalue.parameters.nonCano.sorf.txt repre.valid.ORF.genepred.txt nonCano.sorf.genepred.formatted.txt
genePredToGtf file nonCano.sorf.genepred.formatted.txt nonCano.sorf.formatted.gtf
genome_fa=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
gffread nonCano.sorf.formatted.gtf -g $genome_fa -y nonCano.sorf.fa
seqkit fx2tab nonCano.sorf.fa > nonCano.sorf.tab
##### PRICE #####
# 1. annotation
mkdir -p annot/PIRCE/mm39 && pushd annot/PIRCE/mm39
fa=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
gtf=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.gtf
# 24min
time gedi -e IndexGenome -s $fa -a $gtf -n mm39_gencvM29
popd
mkdir mouse_brain_output_20241011/2015_Science_ChoJ_hippocampus/SRR2163083/output/Ribo-ORFs/PRICE
cd mouse_brain_output_20241011/2015_Science_ChoJ_hippocampus/SRR2163083/output/Ribo-ORFs/PRICE
fq=../../filtered_fq/SRR2163083.trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz
sample=SRR2163083
STARindex=/home/user/data3/lit/resource/genome/mouse/mm39/index/mm39_STARindex
STAR \
--readFilesIn $fq \
--readFilesCommand zcat \
--seedSearchStartLmax 15 \
--runThreadN 40 \
--outFilterMismatchNmax 2 \
--genomeDir $STARindex \
--outFileNamePrefix ${sample}_ \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--outSAMattributes All

samtools index -@ 30 SRR2163083_Aligned.sortedByCoord.out.bam

# 3min
gedi -e Price -reads SRR2163083_Aligned.sortedByCoord.out.bam -genomic mm39_gencvM29 -prefix SRR2163083 -plot

gedi -e ViewCIT -m bed -name 'd.getType()'  SRR2163083.orfs.cit > SRR2163083.orfs.bed

genePredToGtf file nonCano.formatted.gpe nonCano.formatted.gtf
genome_fa=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
gffread nonCano.formatted.gtf -g $genome_fa -y nonCano.fa
seqkit fx2tab -l nonCano.fa | awk '{print $1,$3}' > nonCano.pro.l.txt
awk '$2>=6 && $2<=150' nonCano.pro.l.txt > nonCano.sorf.pro.l.txt
seqkit grep -n -f <(cut -f 1 nonCano.sorf.pro.l.txt) nonCano.fa > nonCano.sorf.fa
seqkit fx2tab nonCano.sorf.fa > nonCano.sorf.tab
##### RiboCode #####
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/RiboCode_annot/mm39
plot_orf_density -a $RiboCode_annot -c metaplots_pre_config.txt -t ENSMUST00000070533.5 -s 3741571 -e 3286245 --start-codon ATG

# stranded-specific
/home/user/BGM/lit/anaconda3/envs/py2/bin/infer_experiment.py -i ../../alignment/SRR2163083_Aligned.sortedByCoord.out.bam \
-q 255 -r /home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.bed -s 200000

# 得到非经典的小肽
compare_script=/home/user/data3/lit/project/sORFs/01-ribo-seq/ref/Ribo-seq-Tool-Comparison-Scripts-v2.0/
cat <(head -n1 SRR2163083_collapsed.txt) <(awk '$2 != "annotated" && $11 <= 450 ' SRR2163083_collapsed.txt) > SRR2163083_collapsed.nonCano.sorf.txt
tail -n +2 SRR2163083_collapsed.nonCano.sorf.txt | cut -f 1 > SRR2163083_collapsed.nonCano.sorf.id.txt
time grep -F -f SRR2163083_collapsed.nonCano.sorf.id.txt SRR2163083_collapsed.gtf > SRR2163083_collapsed.nonCano.sorf.gtf
python $compare_script/Scripts\ for\ RiboCode\ Analysis/Formatting_RiboCode_gtf.py SRR2163083_collapsed.nonCano.sorf.gtf SRR2163083_collapsed.nonCano.sorf.formatted.gtf
python $compare_script/Scripts\ for\ RiboCode\ Analysis/Generate_Fasta_RiboCode.py SRR2163083_collapsed.nonCano.sorf.txt SRR2163083_collapsed.nonCano.sorf.fa
gtfToGenePred SRR2163083_collapsed.nonCano.sorf.formatted.gtf SRR2163083_collapsed.nonCano.sorf.formatted.gpe

cat <(head -n1 SRR2163083.txt) <(awk '$2 != "annotated" && $10 <= 450 ' SRR2163083.txt) > SRR2163083.nonCano.sorf.txt
tail -n +2 SRR2163083.nonCano.sorf.txt | cut -f 1 > SRR2163083.nonCano.sorf.id.txt
time grep -F -f SRR2163083.nonCano.sorf.id.txt SRR2163083.gtf > SRR2163083.nonCano.sorf.gtf
python $compare_script/Scripts\ for\ RiboCode\ Analysis/Formatting_RiboCode_gtf.py SRR2163083.nonCano.sorf.gtf SRR2163083.nonCano.sorf.formatted.gtf
python $compare_script/Scripts\ for\ RiboCode\ Analysis/Generate_Fasta_RiboCode.py SRR2163083.nonCano.sorf.txt SRR2163083.nonCano.sorf.fa
gtfToGenePred SRR2163083.nonCano.sorf.formatted.gtf SRR2163083.nonCano.sorf.formatted.gpe
seqkit fx2tab SRR2163083.nonCano.sorf.fa > SRR2163083.nonCano.sorf.tab
##### Compare #####
cd merge
# 7087
wc -l ../RiboCode/SRR2163083_collapsed.nonCano.sorf.formatted.gpe
# 35413
wc -l ../RiboCode/SRR2163083.nonCano.sorf.formatted.gpe
# 35263
wc -l ../RibORF/nonCano.sorf.formatted.genepred.txt
# 3826
wc -l ../PRICE/nonCano.sorf.pro.l.txt 
# 3477 RiboCode vs RibORF
comm -12 <(cut -f1 ../RiboCode/SRR2163083_collapsed.nonCano.sorf.formatted.gpe| sort -k1,1) <(cut -f1 ../RibORF/nonCano.sorf.genepred.formatted.txt| sort -k1,1) | wc -l
# 1514 RiboCode vs PRICE
comm -12 <(cut -f1 ../RiboCode/SRR2163083_collapsed.nonCano.sorf.formatted.gpe| sort -k1,1) <(cut -f1 ../PRICE/nonCano.sorf.pro.l.txt| sort -k1,1) | wc -l
# 1922 RibORF vs PRICE
comm -12 <(cut -f1 ../RibORF/nonCano.sorf.formatted.genepred.txt| sort -k1,1) <(cut -f1 ../PRICE/nonCano.sorf.pro.l.txt| sort -k1,1) | wc -l

# merge以及去重
cat ../RibORF/nonCano.sorf.fa ../PRICE/nonCano.sorf.fa ../RiboCode/SRR2163083_collapsed.nonCano.sorf.fa > merge.fa
seqkit rmdup -n merge.fa -o merge.rmdup.n.fa
seqkit rmdup -s merge.fa -o merge.rmdup.fa -d dup.fa
seqkit rmdup -s ../RiboCode/SRR2163083_collapsed.nonCano.sorf.fa -o /dev/null -d RiboCode.dup.fa
seqkit seq -s RiboCode.dup.fa -o RiboCode.dup.seq.txt
grep -B 1 -f RiboCode.dup.seq.txt ../RiboCode/SRR2163083_collapsed.nonCano.sorf.fa > RiboCode.dup.all.fa

seqkit rmdup -s ../PRICE/nonCano.sorf.fa -o /dev/null -d PRICE.dup.fa
seqkit seq -s PRICE.dup.fa -o PRICE.dup.seq.txt
grep -B 1 -f PRICE.dup.seq.txt ../PRICE/nonCano.sorf.fa > PRICE.dup.all.fa

seqkit rmdup -s ../RibORF/nonCano.sorf.fa -o /dev/null -d RibORF.dup.fa
seqkit seq -s RibORF.dup.fa -o RibORF.dup.seq.txt
grep -B 1 -f RibORF.dup.seq.txt ../RibORF/nonCano.sorf.fa > RibORF.dup.all.fa

seqkit tab2fx RibORF.nonCano.sorf.rmDup.tab > RibORF.nonCano.sorf.rmDup.fa
seqkit tab2fx PRICE.nonCano.sorf.rmDup.tab > PRICE.nonCano.sorf.rmDup.fa
seqkit tab2fx RiboCode.nonCano.sorf.rmDup.tab > RiboCode.nonCano.sorf.rmDup.fa

cat *nonCano.sorf.rmDup.fa > merge.1.fa
seqkit rmdup -s  merge.1.fa -o  merge.1.rmDup.fa
seqkit rmdup -n merge.1.fa -o  merge.1.n.rmDup.fa