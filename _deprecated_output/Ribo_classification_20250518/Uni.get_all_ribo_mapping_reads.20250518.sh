# 输入基本质控后的单个fastq文件，根据文件名输出结构化的输出结果
## 第一层为研究名，第二层为样本名
## 输出output/filtered_fq以及output/alignment

# cleaned fastq file $1
cleaned_fastq=$1
output_path=$2
[ -d $output_path ] || mkdir -p $output_path
### 根据输入的文件名替换 ###
# 20250518修改 #
sample=$(basename -s .fq.gz $cleaned_fastq|sed 's/.trimmed.rRNA.tRNA.snoRNA.unaligned//')

### 参考基因组注释文件，根据物种替换 ###
#### 修改 ####
bowtie2_rRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/ncRNA/human/merged/hg38.rRNA
bowtie2_tRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/ncRNA/human/merged/hg38.tRNA
bowtie2_snoRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/ncRNA/human/merged/hg38.snoRNA
STARindex=/home/user/data3/lit/resource/genome/human/hg38/hg38_STARindex_v2.5.2b/

# 0.创立文件夹
cd $output_path
### 根据输入的文件名替换 ###
subdir=$(echo $cleaned_fastq | cut -d "/" -f 12)
[ -d $subdir/$sample ] || mkdir -p $subdir/$sample
cd $subdir/$sample
[ -d output ] || mkdir output
[ -d log ] || mkdir log

# 激活虚拟环境
source activate biotools

# 3. 比对到参考基因组上
echo -e "Mapping at $(date '+%Y-%m-%d %H:%M:%S')"
[ -d output/alignment ] || mkdir output/alignment
# 20250521 内存不足增加--limitBAMsortRAM 35833093247选项，后续不需要可以去除
STAR \
--readFilesIn $cleaned_fastq \
--limitBAMsortRAM 35833093247 \
--readFilesCommand zcat \
--seedSearchStartLmax 15 \
--runThreadN 40 \
--outFilterMismatchNmax 2 \
--genomeDir $STARindex \
--outFileNamePrefix output/alignment/${sample}_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outFilterMultimapNmax 9999 \
--outFilterMatchNmin 16 \
--outSAMattributes All \
--alignEndsType EndToEnd &> log/STAR.log

# samtools view output/alignment/${sample}_Aligned.sortedByCoord.out.bam -h -@ 20 > output/alignment/$sample.sam