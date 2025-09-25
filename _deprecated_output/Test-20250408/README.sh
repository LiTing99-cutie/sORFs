script_1=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/Uni.qc.single.sh
script_2=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250306/Uni.call.orfs.ribocode.20250306.sh
cp /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250306/Uni.call.orfs.ribocode.20250306.sh ./
script_3=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/Uni.call.orfs.ribocode.human.20250408.sh

##### mouse
bash $script_1 $PWD/E16_0321.raw.fastq.gz AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC $PWD/output yes yes yes
trimmed_fq=$PWD/output/E16_0321.raw/output/trimmed_fastq/E16_0321.raw_trimmed.fq.gz
bash $script_2 $trimmed_fq $PWD/output
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/RiboCode_annot/mm39
bam=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/output/E16_0321.raw/output/alignment/E16_0321.raw_Aligned.toTranscriptome.out.bam
metaplots -a $RiboCode_annot -r $bam -f0_percent 0 -pv1 1 -pv2 1
metaplots -a $RiboCode_annot -r $bam -f0_percent 0.5 -pv1 0.01 -pv2 0.01 -o meta_1
##### human
sample=p21_0321
bash $script_1 $PWD/$sample.raw.fastq.gz AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC $PWD/output yes yes yes
trimmed_fq=$PWD/output/$sample.raw/output/trimmed_fastq/$sample.raw_trimmed.fq.gz
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/annotation/RiboCode_annot
bam=$PWD/output/$sample.raw/output/alignment/$sample.raw_Aligned.toTranscriptome.out.bam
bash $script_3 $trimmed_fq $PWD/output
cd $PWD/output/$sample.raw/output/Ribo-ORFs/RiboCode/
metaplots -a $RiboCode_annot -r $bam -f0_percent 0.5 -pv1 0.01 -pv2 0.01 -o meta_1
metaplots -a $RiboCode_annot -r $bam -f0_percent 0.5 -pv1 1 -pv2 1 -o meta_2
metaplots -a $RiboCode_annot -r $bam -f0_percent 0 -pv1 1 -pv2 1 -o meta_3

bam_1=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/output/p21_0321.raw/output/alignment/p21_0321.raw_Aligned.sortedByCoord.out.bam
script=/home/user/BGM/lit/anaconda3/envs/py2/bin/read_distribution.py
define_annotation_gencode_v41_human
$script -i $bam_1 -r /home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.bed > $sample.r_distri.txt

samtools view "$bam_1" | \
        awk '{print length($10)}' | \
        sort | uniq -c | \
        awk -v sample="$sample" '{print sample, $2, $1}' > $sample.r_len_distri.txt
define_annotation_gencode_v41_human
samtools index $bam_1

conda activate biotools
sample_1=p21_0321.raw
STARindex=/home/user/data3/lit/resource/genome/human/hg38/hg38_STARindex_v2.5.2b/
cd /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/output/p21_0321.raw
cleaned_fastq=output/filtered_fq/$sample_1.trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz
[ -d output/alignment ] || mkdir output/alignment
STAR \
--readFilesIn $cleaned_fastq \
--readFilesCommand zcat \
--seedSearchStartLmax 15 \
--runThreadN 40 \
--outFilterMismatchNmax 2 \
--genomeDir $STARindex \
--outFileNamePrefix output/alignment/${sample_1}_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outFilterMultimapNmax 1 \
--outFilterMatchNmin 16 \
--outSAMattributes All \
--alignEndsType EndToEnd &> log/STAR.log
conda activate ribocode
samtools index $bam_1
ribotish quality -b $bam_1 -g $gtf

$script -i $bam_1 -q 255 -r $anno_bed > $sample.stranded.txt

# 具体check下和tRNA的比对情况
conda activate biotools
bowtie2_tRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/ncRNA/human/merged/hg38.tRNA
bowtie2 -q \
--phred33 \
-N 1 \
-x $bowtie2_tRNA_index  \
-U output/trimmed_fastq/p21_0321.raw_trimmed.fq.gz \
-S output/filtered_fq/$sample_1.trimmed.tRNA.aligned.sam \
-t \
-p 20 \
--no-unal \
--un-gz /dev/null

bowtie2 -q \
--phred33 \
-N 1 \
-x $bowtie2_tRNA_index  \
-U output/trimmed_fastq/p21_0321.raw_trimmed.fq.gz \
-S output/filtered_fq/$sample_1.trimmed.tRNA.aligned.sam \
-t \
-p 20 \
--no-unal \
--un-gz /dev/null \
--al-gz output/filtered_fq/$sample_1.trimmed.tRNA.aligned.fq.gz

STAR \
--readFilesIn /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/output/p21_0321.raw/output/filtered_fq/p21_0321.raw.trimmed.tRNA.aligned.fq.gz \
--readFilesCommand zcat \
--seedSearchStartLmax 15 \
--runThreadN 40 \
--outFilterMismatchNmax 2 \
--genomeDir $STARindex \
--outFileNamePrefix output/alignment/${sample_1}_tRNA \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outFilterMultimapNmax 1 \
--outFilterMatchNmin 16 \
--outSAMattributes All \
--alignEndsType EndToEnd &> log/STAR.tRNA.log

# 根据meta_1_pre_config.txt来预测ORF
conda activate ribocode
define_annotation_gencode_v41_human
sample=p21_0321.raw
cd $PWD/output/$sample/output/Ribo-ORFs/RiboCode/
bam=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/output/p21_0321.raw/output/alignment/p21_0321.raw_Aligned.toTranscriptome.out.bam
[ -f p21_0321.raw_Aligned.toTranscriptome.out_psites.hd5 ] && rm -rf p21_0321.raw_Aligned.toTranscriptome.out_psites.hd5
metaplots -a $RiboCode_annot -r $bam -f0_percent 0.5 -pv1 0.01 -pv2 0.01 -o meta_1 && \
RiboCode -a $RiboCode_annot -c meta_1_pre_config.txt -l no -g -o ./${sample}_config_1 --output-gtf --output-bed --min-AA-length 6 

##### 在mouse上也补上这些分析
bash Run.mouse.add.sh