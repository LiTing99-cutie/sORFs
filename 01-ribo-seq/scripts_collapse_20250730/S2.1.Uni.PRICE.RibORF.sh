set -eo pipefail

# 结果输出路径
output_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/

# 基因组及其注释
fa=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
genePred=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.genePred.txt
gtf=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.gtf
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/RiboCode_annot/mm39
STARindex=/home/user/data3/lit/resource/genome/mouse/mm39/index/mm39_STARindex
candidateORF=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/RiboORF/mm39/candidateORF.genepred.txt

# 对于每一个样本结果路径进行循环
sample_output_path=$1
sample=$(echo $sample_output_path | awk -F'/' '{print $(NF-1)}')

##### PRICE #####
echo -e "Price start at $(date '+%Y-%m-%d %H:%M:%S')"
mkdir -p $sample_output_path/Ribo-ORFs/PRICE
cd $sample_output_path/Ribo-ORFs/PRICE
fq=../../filtered_fq/*.trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz
source activate biotools
STAR \
--readFilesIn $fq \
--readFilesCommand zcat \
--seedSearchStartLmax 15 \
--runThreadN 5 \
--outFilterMismatchNmax 2 \
--genomeDir $STARindex \
--outFileNamePrefix ${sample}_ \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 1 \
--alignEndsType EndToEnd \
--outSAMattributes All
samtools index -@ 5 *_Aligned.sortedByCoord.out.bam
gedi -e Price -reads *_Aligned.sortedByCoord.out.bam -genomic mm39_gencvM29 -prefix $sample -plot
gedi -e ViewCIT -m bed -name 'd.getType()'  *.orfs.cit > $sample.orfs.bed
echo -e "Price done at $(date '+%Y-%m-%d %H:%M:%S')"

##### RibORF #####
echo -e "RibORF start at $(date '+%Y-%m-%d %H:%M:%S')"
mkdir -p $sample_output_path/Ribo-ORFs/RibORF
RibORF_path=/home/user/data2/lit/software/RibORF/RibORF.2.0

# 0. build annotation; see /home/user/data3/lit/project/sORFs/01-ribo-seq/Build_annotation.sh
# 1. bam2sam
cd $sample_output_path/Ribo-ORFs/RibORF
samtools view ../../alignment/*_Aligned.sortedByCoord.out.bam -h -@ 5 > ./$sample.sam

# 2. output
# 使用ribocode鉴定read length以及offset
source activate ribocode
metaplots -a $RiboCode_annot -r ../../alignment/*_Aligned.toTranscriptome.out.bam
grep "^#" metaplots_pre_config.txt|awk -v OFS='\t' '{print $2,$4+3}'|tail -n +3|head -n -1 > offset.correction.parameters.txt
pwd
cat -A offset.correction.parameters.txt
# 根据offset文件校正
echo "perl $RibORF_path/offsetCorrect.pl -r $sample.sam -p offset.correction.parameters.txt -o corrected.$sample.sam"
time perl $RibORF_path/offsetCorrect.pl -r $sample.sam -p offset.correction.parameters.txt -o corrected.$sample.sam
# 可视化校正后的结果
time perl $RibORF_path/readDist.pl -f corrected.$sample.sam -g $genePred -o ./ -d 1
# call ORFs
time perl $RibORF_path/ribORF.pl -f corrected.$sample.sam -c $candidateORF -o ./
echo -e "RibORF done at $(date '+%Y-%m-%d %H:%M:%S')"

