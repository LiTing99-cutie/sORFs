script=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/Uni.qc.single.sh
bash $script /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test_20250306/E16_Y40_RPF.fq.gz none $PWD/output yes no no
bash $script /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test_20250306/E16_Y40_RPF.fq.gz AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC $PWD/output no yes yes

bash $script /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test_20250306/E16_control.fq.gz AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC $PWD/output yes yes yes

cp /home/user/data3/lit/project/sORFs/01-ribo-seq/S1.1.Uni.Ribo.Run.20241011.sh ./Uni.call.orfs.ribocode.20250306.sh
fq=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test_20250306/output/E16_Y40_RPF/output/trimmed_fastq/E16_Y40_RPF_trimmed.fq.gz
bash Uni.call.orfs.ribocode.20250306.sh $fq $PWD/output

fq=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test_20250306/output/E16_Y40_RPF/output/filtered_fq/E16_Y40_RPF.trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz
fastqc -o output/E16_Y40_RPF/output/rm_contam_fastqc -t 10 $fq


# [lit@rhesusbase sORFs]$ less /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test_20250306/E16_control.fq.gz |grep GCGTACACGGAGTCAAGACCCGCAACGCAC
# GCGTACACGGAGTCAAGACCCGCAACGCACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTG
# GCGTACACGGAGTCAAGACCCGCAACGCACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTA
# GCGTACACGGAGTCAAGACCCGCAACGCACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTG
# GCGTACACGGAGTCAAGACCCGCAACGCACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTA
# [lit@rhesusbase sORFs]$ seqkit stat /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test_20250306/E16_control.fq.gz
# file                                                                                     format  type    num_seqs        sum_len  min_len  avg_len  max_len
# /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test_20250306/E16_control.fq.gz  FASTQ   DNA   23,953,175  1,796,488,125       75       75       75
fq=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test_20250306/output/E16_Y40_RPF/output/trimmed_fastq/E16_Y40_RPF_trimmed.fq.gz
seqkit seq -s $fq |sort|uniq -c > E16_Y40.uniq.seq.txt

sort -k1,1nr E16_Y40.uniq.seq.txt > E16_Y40.uniq.sorted.seq.txt

grep -E 'ATGTACACGGAGTCGACCCGCAACGCTG|GCGTACACGGAGTCAAGACCCGCAACGCAC' E16_Y40.uniq.sorted.seq.txt|sort -k1,1nr
grep -E ' ATGTACACGGAGTCGACCCGCAACGCTG$| GCGTACACGGAGTCAAGACCCGCAACGCAC$' E16_Y40.uniq.sorted.seq.txt

fq=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test_20250306/output/E16_control.fq.gz/output/trimmed_fastq/E16_control_trimmed.fq.gz
seqkit seq -s $fq |sort|uniq -c > E16_control.uniq.seq.txt
sort -k1,1nr E16_control.uniq.seq.txt > E16_control.uniq.sorted.seq.txt
grep -E 'ATGTACACGGAGTCGACCCGCAACGCTG|GCGTACACGGAGTCAAGACCCGCAACGCAC' E16_control.uniq.sorted.seq.txt|sort -k1,1nr
grep -E ' ATGTACACGGAGTCGACCCGCAACGCTG$| GCGTACACGGAGTCAAGACCCGCAACGCAC$' E16_control.uniq.sorted.seq.txt

# 放松阈值 20250409
sample=E16_Y40_RPF
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/RiboCode_annot/mm39
bam=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250306/output/E16_Y40_RPF/output/alignment/E16_Y40_RPF_Aligned.toTranscriptome.out.bam
cd $PWD/output/$sample/output/Ribo-ORFs/RiboCode/
conda activate ribocode
metaplots -a $RiboCode_annot -r $bam -f0_percent 0.5 -pv1 0.01 -pv2 0.01 -o meta_1
metaplots -a $RiboCode_annot -r $bam -f0_percent 0.5 -pv1 1 -pv2 1 -o meta_2
metaplots -a $RiboCode_annot -r $bam -f0_percent 0 -pv1 1 -pv2 1 -o meta_3