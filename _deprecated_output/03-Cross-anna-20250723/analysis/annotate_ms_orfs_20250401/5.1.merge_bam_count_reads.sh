# 根据链特异性来计算表达量
## 首先合并相同类型的libary的文件，其次计算表达量
bam_lst=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_rna_seq_alignment_assemble/bam.lst
source /home/user/data2/lit/bin/lit_utils.sh 
define_annotation_gencode_v41_human
mkdir -p output/S5/merged_bam && cd output/S5/merged_bam
samtools merge -@ 20 -o merged.1.bam $(egrep human_brain_rna_[0-9] $bam_lst)
samtools merge -@ 20 -o merged.2.bam $(egrep SRR15906 $bam_lst)
samtools merge -@ 20 -o merged.3.bam $(egrep SRR1551322 $bam_lst)
samtools merge -@ 20 -o merged.4.bam $(egrep SRR156254 $bam_lst)
## featureCounts的链特异性参数-s需要根据特定的文库来修改
featureCounts -s 1 -T 10 -a $gtf -o rna_counts_1.txt merged.1.bam
featureCounts -s 2 -T 10 -a $gtf -o rna_counts_2.txt merged.2.bam
featureCounts -s 2 -p --countReadPairs -T 10 -a $gtf -o rna_counts_3.txt merged.3.bam
featureCounts -T 10 -a $gtf -o rna_counts_4.txt merged.4.bam