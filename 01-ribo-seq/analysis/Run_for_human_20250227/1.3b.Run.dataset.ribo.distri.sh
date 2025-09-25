# 查看是什么原因无法call orfs
bam=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227/2022_NN/SRR15906424/output/alignment/SRR15906424_Aligned.sortedByCoord.out.bam
/home/user/BGM/lit/anaconda3/envs/py2/bin/read_distribution.py -i $bam -r /home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.bed

### 测试homer
mkdir tmp
bedtools bamtobed -i $bam > tmp/your_bam_file.bed
annotatePeaks.pl tmp/your_bam_file.bed hg38 > tmp/output.txt
annotatePeaks.pl tmp/your_bam_file.bed hg38 -gtf /home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf > tmp/output.1.txt

### 测试rseqc 【选择了这个】
cd human_brain_output_20250227
mkdir reads_distri_20250328
find ./ -name *_Aligned.sortedByCoord.out.bam > reads_distri_20250328/bam.lst
script=/home/user/BGM/lit/anaconda3/envs/py2/bin/read_distribution.py
for bam in $(cat reads_distri_20250328/bam.lst);do
sample=$(echo $bam|awk -F'alignment/|_Aligned' '{print $2}')
$script -i $bam -r /home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.bed > reads_distri_20250328/$sample.r_distri.txt
done

cd reads_distri_20250328
out_file="summary_tag_count.tsv"
[ -f $out_file ] && rm -rf $out_file
echo -e "Sample\tType\tTag_count" > $out_file
for file in *r_distri.txt; do
    sample_name=$(basename "$file" .r_distri.txt)
    awk -v sample="$sample_name" 'NR>=6 && NR<=15 {print sample"\t"$1"\t"$3}' "$file" >> $out_file
done

### 测试alfa
cd tmp
conda activate base
GTF_FILE=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
CHR_LENGTHS_FILE=/home/user/data/lit/database/public/genome/hg38/hg38_comChr.chrom.sizes
sed 's/gene_type/gene_biotype/' $GTF_FILE > hg38.gtf
python ALFA-master/alfa.py -a hg38.gtf -g hg38 --chr_len $CHR_LENGTHS_FILE -p 10

bam=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227/2022_NN/SRR15906424/output/alignment/SRR15906424_Aligned.sortedByCoord.out.bam
script=/home/user/BGM/lit/anaconda3/envs/py2/bin/infer_experiment.py
anno_bed=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.bed
$script -i $bam -q 255 -r $anno_bed
cd alfa_output
python ../ALFA-master/alfa.py -g ../hg38 --bam $bam bam_1 -s unstranded -d 3 --pdf output.pdf -p 20