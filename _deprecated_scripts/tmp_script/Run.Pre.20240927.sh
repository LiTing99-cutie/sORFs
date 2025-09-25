mkdir Pre-Run
cd Pre-Run
# 0.创立文件夹
[ -d output ] || mkdir output
[ -d log ] || mkdir log
# 1.原始测序数据质量控制
## fastqc
[ -d output/fastqc ] || mkdir output/fastqc
mouse_example_fastq=/home/user/data/lit/project/ZNF271/data/ribo-seq/2015-Science-GSE72064/mouse/hippocampus/SRR2163103_GSM1854003_RPF_control_rep3_Mus_musculus_RNA-Seq.fastq.gz
fastqc -o output/fastqc -t 10 $mouse_example_fastq
//这里可以加上log输出
## 去接头，过滤低质量reads
### trim_galore
conda activate biotools
nohup bash -c "time trim_galore -j 8 -q 20 --length 20 $mouse_example_fastq --gzip -o /home/user/data/lit/project/ZNF271/data/ribo-seq/2015-Science-GSE72064/mouse/ --fastqc_args "-o output/fastqc -t 10"" &> log/trim_galore.log &

### fastp
nohup bash -c "time fastp \
  -i $mouse_example_fastq \
  -o /home/user/data/lit/project/ZNF271/data/ribo-seq/2015-Science-GSE72064/mouse/trimmed_reads.fq.gz \
  --cut_mean_quality 20 \
  --cut_window_size 4 \
  -5 \
  -3 \
  --length_required 20 \
  --thread 8 \
  -z 4 \
  --html /home/user/data/lit/project/ZNF271/data/ribo-seq/2015-Science-GSE72064/mouse/fastp_report.html \
  --json /home/user/data/lit/project/ZNF271/data/ribo-seq/2015-Science-GSE72064/mouse/fastp_report.json " &> log/fastp.log &

#### 可能--fastqc_args中指定的-o输出文件需要和trim后的文件在同一个相对路径下，需要改成绝对路径
fastqc -o output/fastqc -t 10 /home/user/data/lit/project/ZNF271/data/ribo-seq/2015-Science-GSE72064/mouse/SRR2163103_GSM1854003_RPF_control_rep3_Mus_musculus_RNA-Seq_trimmed.fq.gz

# 2.去除contaminant
trimmed_fastq=/home/user/data/lit/project/ZNF271/data/ribo-seq/2015-Science-GSE72064/mouse/SRR2163103_GSM1854003_RPF_control_rep3_Mus_musculus_RNA-Seq_trimmed.fq.gz
cd /home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/ncRNA/merged 
bowtie2-build hg38.rRNA.tRNA.snoRNA.fa hg38.rRNA.tRNA.snoRNA
nohup bowtie2 -q \
--phred33 \
-N 1 \
-x /home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/ncRNA/merged/hg38.rRNA.tRNA.snoRNA \
-U $trimmed_fastq \
-S output/filtered_fq/SRR2163103.trimmed.aligned.sam \
-t \
-p 10 \
--no-unal \
--un-gz output/filtered_fq/SRR2163103.trimmed.unaligned.fq.gz &> log/bowtie.log &

## 统计去除前和去除后，计算比例
n1=`seqkit stats $trimmed_fastq | tail -n1 |awk '{print $4}'|sed -e 's/,//g'`
n2=`seqkit stats output/filtered_fq/SRR2163103.trimmed.unaligned.fq.gz | tail -n1 |awk '{print $4}'|sed -e 's/,//g'`
awk "BEGIN{print $n2/$n1}"

## 使用STAR比对应该更快，逐步比对去除污染
### 并没有更快
trimmed_fastq=/home/user/data/lit/project/ZNF271/data/ribo-seq/2015-Science-GSE72064/mouse/SRR2163103_GSM1854003_RPF_control_rep3_Mus_musculus_RNA-Seq_trimmed.fq.gz
rRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/ncRNA/STARindex/rRNA/rRNA_STARindex
mkdir filtered_fq_STAR
nohup STAR --runThreadN 10 \
--genomeDir $rRNA_index \
--readFilesIn $trimmed_fastq \
--readFilesCommand zcat \
--outFileNamePrefix output/filtered_fq_STAR/SRR2163103.trimmed.filtered.rRNA \
--outSAMtype BAM Unsorted \
--outReadsUnmapped Fastx \
--outFilterMismatchNmax 1 &> log/STAR.rRNA.log &

nohup STAR --runThreadN 10 \
--genomeDir $rRNA_index \
--readFilesIn $trimmed_fastq \
--readFilesCommand zcat \
--outFileNamePrefix output/filtered_fq/${sample}.rRNA. \
--outSAMtype None \
--outReadsUnmapped Fastx \
--outFilterMismatchNmax 1 &> log/STAR.rRNA.20250724.log &


## 只去除ensembl rRNA
trimmed_fastq=/home/user/data/lit/project/ZNF271/data/ribo-seq/2015-Science-GSE72064/mouse/SRR2163103_GSM1854003_RPF_control_rep3_Mus_musculus_RNA-Seq_trimmed.fq.gz
cd /home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/ncRNA/ 
bowtie2-build Homo_sapiens.GRCh38.rRNA.fa Homo_sapiens.GRCh38.rRNA
nohup bowtie2 -q \
--phred33 \
-N 1 \
-x /home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/ncRNA/Homo_sapiens.GRCh38.rRNA  \
-U $trimmed_fastq \
-S output/filtered_fq/SRR2163103.trimmed.aligned.ensembl.rRNA.sam \
-t \
-p 10 \
--no-unal \
--un-gz output/filtered_fq/SRR2163103.trimmed.unaligned.ensembl.rRNA.fq.gz &> log/bowtie.ensembl.rRNA.log &

n3=`seqkit stats output/filtered_fq/SRR2163103.trimmed.unaligned.ensembl.rRNA.fq.gz | tail -n1 |awk '{print $4}'|sed -e 's/,//g'`
awk "BEGIN{print $n3/$n1}"

# 只去除tRNA,rRNA,snoRNA并进行比较
pushd /home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/ncRNA/merged
bowtie2-build hg38.tRNA.fa hg38.tRNA
bowtie2-build hg38.rRNA.fa hg38.rRNA
bowtie2-build hg38.snoRNA.fa hg38.snoRNA
popd

mkdir -p output/filtered_fq/tRNA output/filtered_fq/rRNA output/filtered_fq/snoRNA

trimmed_fastq=/home/user/data/lit/project/ZNF271/data/ribo-seq/2015-Science-GSE72064/mouse/SRR2163103_GSM1854003_RPF_control_rep3_Mus_musculus_RNA-Seq_trimmed.fq.gz
n1=151866345
[ -f output/filtered_fq/types.percentage.txt ] && rm -rf output/filtered_fq/types.percentage.txt
for type in tRNA rRNA snoRNA;do
echo $type
bowtie2 -q \
--phred33 \
-N 1 \
-x /home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/ncRNA/merged/hg38.$type  \
-U $trimmed_fastq \
-S /dev/null \
-t \
-p 20 \
--no-unal \
--un-gz output/filtered_fq/${type}/SRR2163103.trimmed.unaligned.fq.gz &> log/bowtie.hg38.$type.log
n=`seqkit stats output/filtered_fq/$type/SRR2163103.trimmed.unaligned.fq.gz | tail -n1 |awk '{print $4}'|sed -e 's/,//g'`
awk "BEGIN{print $n/$n1}" >> output/filtered_fq/types.percentage.txt
echo "Done for $type"
done

samtools flagstat output/filtered_fq/SRR2163103.trimmed.aligned.ensembl.rRNA.sam > output/filtered_fq/flagstat.ensembl.rRNA.txt
samtools flagstat output/filtered_fq/SRR2163103.trimmed.aligned.sam > output/filtered_fq/flagstat.txt

# 3. 比对到参考基因组上
mkdir output/alignment
sample=SRR2163103
index=/home/user/data3/lit/resource/genome/mouse/mm39/index/mm39_STARindex
nohup bash -c "time STAR \
--readFilesIn output/filtered_fq/SRR2163103.trimmed.unaligned.fq.gz \
--readFilesCommand zcat \
--seedSearchStartLmax 15 \
--runThreadN 40 \
--outFilterMismatchNmax 2 \
--genomeDir $index \
--outFileNamePrefix output/alignment/${sample}_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outFilterMultimapNmax 1 \
--outFilterMatchNmin 16 \
--alignEndsType EndToEnd" &> log/STAR.log &

FA=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
GTF=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.gtf
mkdir output/Ribo-ORFs
mkdir output/Ribo-ORFs/SRR2163103

RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/RiboCode_annot/mm39
bam=/home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/output/alignment/SRR2163103_Aligned.toTranscriptome.out.bam
# 正确
nohup bash -c "time (cd output/Ribo-ORFs/SRR2163103 && metaplots -a $RiboCode_annot -r $bam && \
RiboCode -a $RiboCode_annot -c metaplots_pre_config.txt -l no -g -o ./SRR2163103 --output-gtf --output-bed --min-AA-length 6)" &> log/anno_orfs.log &

wc -l output/Ribo-ORFs/SRR2163103/SRR2163103_collapsed.txt
less output/Ribo-ORFs/SRR2163103/SRR2163103_collapsed.txt | awk '$11<100'|wc -l