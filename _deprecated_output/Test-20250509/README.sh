bin_dir=/home/user/data3/lit/project/sORFs/06-RNA-seq/bin
script_1=$bin_dir/Uni.qc.single.sh
script_2=$bin_dir/Uni.qc.batch.v1.sh
script_3=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/Uni.call.orfs.ribocode.human.20250408.sh
cd /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250509
mkdir -p 01-output/fastqc
# ln -s 
mkdir 02-rawdata
raw_fastq_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/data-20250509/raw_data/
find $raw_fastq_path -name "*fastq.gz" >  01-output/fastqc/raw_fastq.lst
# qc
bash $script_2 01-output/fastqc \
<(cat 01-output/fastqc/raw_fastq.lst) \
yes yes yes $script_1 "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
# call orfs
for trimmed_fq in $(find $PWD -name "*raw_trimmed.fq.gz");do
bash $script_3 $trimmed_fq $PWD/01-output/call-orfs
done
# check on read length distribution and quality
define_annotation_gencode_v41_human
conda activate ribocode
## each 4min
for bam in $(find $PWD -name "*Aligned.sortedByCoord.out.bam");do
[ -f $bam.bai ] || samtools index $bam
ribotish quality -b $bam -g $gtf -p 30
done

# 卡松阈值
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/annotation/RiboCode_annot
for bam in $(find $PWD -name "*_Aligned.toTranscriptome.out.bam");do
sample=$(basename -s _Aligned.toTranscriptome.out.bam $bam)
pushd $PWD/01-output/call-orfs/$sample
metaplots -a $RiboCode_annot -r $bam -f0_percent 0 -pv1 1 -pv2 1 -o meta_3
popd
done
