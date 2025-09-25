#!/bin/bash
sample=E16_0321.raw
sample_output_path=$PWD/output/${sample}
bam=$sample_output_path/output/alignment/${sample}_Aligned.toTranscriptome.out.bam
bam_1=$sample_output_path/output/alignment/${sample}_Aligned.sortedByCoord.out.bam
script=/home/user/BGM/lit/anaconda3/envs/py2/bin/read_distribution.py
cleaned_fastq=output/filtered_fq/$sample.trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz
source /home/user/data2/lit/bin/lit_utils.sh
define_annotation_gencode_v41_mouse

# source activate ribocode
# cd $PWD/output/$sample/output/Ribo-ORFs/RiboCode/
# metaplots -a $RiboCode_annot -r $bam -f0_percent 0.5 -pv1 0.01 -pv2 0.01 -o meta_1
# metaplots -a $RiboCode_annot -r $bam -f0_percent 0.5 -pv1 1 -pv2 1 -o meta_2
# metaplots -a $RiboCode_annot -r $bam -f0_percent 0 -pv1 1 -pv2 1 -o meta_3

cd $sample_output_path
$script -i $bam_1 -r $bed > $sample.r_distri.txt

source activate biotools
[ -d output/alignment ] || mkdir output/alignment
STAR \
--readFilesIn $cleaned_fastq \
--readFilesCommand zcat \
--seedSearchStartLmax 15 \
--runThreadN 40 \
--outFilterMismatchNmax 2 \
--genomeDir $STARindex \
--outFileNamePrefix output/alignment/${sample}_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outFilterMultimapNmax 1 \
--outFilterMatchNmin 16 \
--outSAMattributes All \
--alignEndsType EndToEnd &> log/STAR.log

source activate ribocode
samtools index $bam_1
ribotish quality -b $bam_1 -g $gtf