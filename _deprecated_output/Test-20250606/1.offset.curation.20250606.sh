bam_1=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250606/01-output/call-orfs/p21_0523_1/output/alignment/p21_0523_1_Aligned.toTranscriptome.out.bam
bam_2=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250606/01-output/call-orfs/p21_0523_1/output/alignment/p21_0523_1_Aligned.sortedByCoord.out.bam
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/annotation/RiboCode_annot
script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250606/Test.offset.curation.20250606.R

mkdir -p test_output_20250606 && cd test_output_20250606

# 方法一
source activate ribocode
# Warning:The predicted P-site location(12) for length 31 is not the highest peak(13)
metaplots -a $RiboCode_annot -r $bam_1 -f0_percent 0.5 -pv1 0.01 -pv2 0.01 -o ribocode

# 方法二
source /home/user/data2/lit/bin/lit_utils.sh
define_annotation_gencode_v41_human
[ -f $bam_2.bai ] || samtools index $bam_2
ribotish quality -b $bam_2 -g $gtf --th 0.5 -p 30 -l 24,36 -f ribotish_qual.pdf -r ribotish_para.py -o ribotish_qual.txt

# 合并
grep "^#" ribocode_pre_config.txt | awk -v OFS='\t' '{print $2,$4}' | tail -n +3 |head -n -1|sort -k1,1n > ribocode_offset_tab.txt
cat ribotish_para.py | grep -oE '[0-9]+: [0-9]+' | awk -F': ' '!seen[$1]++ {print $1"\t"$2}' > ribotish_offset_tab.txt

Rscript $script ribocode_offset_tab.txt ribotish_offset_tab.txt $bam_2 all_offset_tab.txt

