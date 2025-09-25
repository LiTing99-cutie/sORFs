cd /home/user/data3/lit/project/sORFs/08-Iso-seq/rawdata/merge
source activate biotools
bam_1=../r84130_250703_001_1_A01/p21-IsoSeq.Iso_bc02.bcM0004.ISO.bam
bam_2=../sequencing_add_20250716/p21-IsoSeq.Iso_bc02.bcM0004.ISO.bam

# 生成时间戳函数
timestamp() {
    echo "$(date +"%Y-%m-%d %H:%M:%S")"
}

echo "=== 脚本开始运行 ==="
echo "开始时间: $(timestamp)"

echo ""
echo "=== 步骤1: 合并BAM文件 ==="
echo "合并开始时间: $(timestamp)"
samtools merge -@ 30 -h $bam_1 p21-IsoSeq.Iso_bc02.bcM0004.ISO.merged.tagged.20250717.bam $bam_1 $bam_2
echo "合并完成时间: $(timestamp)"

echo ""
echo "=== 步骤2: 生成PBI索引 ==="
echo "索引生成开始时间: $(timestamp)"
pbindex -j 30 p21-IsoSeq.Iso_bc02.bcM0004.ISO.merged.tagged.20250717.bam
echo "索引生成完成时间: $(timestamp)"

echo ""
echo "=== 步骤3: 验证结果 ==="
echo "验证开始时间: $(timestamp)"
samtools view -H p21-IsoSeq.Iso_bc02.bcM0004.ISO.merged.tagged.20250717.bam | grep -i sample > grep.sample.txt
echo "验证完成时间: $(timestamp)"

echo ""
echo "=== 脚本执行完成 ==="
echo "结束时间: $(timestamp)"