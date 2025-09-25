cd /home/user/data3/lit/project/sORFs/08-Iso-seq/rawdata/merge
source activate biotools
bam_1=../r84130_250703_001_1_A01/p21-IsoSeq.Iso_bc02.bcM0004.ISO.bam
bam_2=../sequencing_add_20250716/p21-IsoSeq.Iso_bc02.bcM0004.ISO.bam
merged_bam=p21-IsoSeq.Iso_bc02.bcM0004.ISO.merged.20250716.bam

echo "Merging BAM files..."
samtools merge -@ 30 $merged_bam $bam_1 $bam_2

echo "Generating PBI index..."
pbindex -j 30 $merged_bam

echo "Done!"

# 检查reads数量
echo "BAM1 reads: $(samtools view -@ 30 -c $bam_1)"
echo "BAM2 reads: $(samtools view -@ 30 -c $bam_2)"
echo "Merged reads: $(samtools view -@ 30 -c $merged_bam)"

# 检查是否有重复
echo "Expected total: $(($(samtools view -@ 30 -c $bam_1) + $(samtools view -@ 30 -c $bam_2)))"

# 合并前后的文件大小并不是简单的加和，但是reads数目经过相加是正确的
# =====================
# 合并结果记录（20250716）
# BAM1 reads: 18270512
# BAM2 reads: 2761568
# Merged reads: 21032080
# Expected total: 21032080
# 结论: 合并成功，reads数量完全匹配，无重复reads
# =====================