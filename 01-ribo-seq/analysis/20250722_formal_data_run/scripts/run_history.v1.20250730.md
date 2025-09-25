# 1一、qc，去除污染，比对
## 生成原始数据list文件
raw_data_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/in_house_human_organized_20250625
ls $raw_data_path/*.fq.gz > ../processed/raw_fastq_batch_1.20250723.lst
## qc去接头，去除污染，比对
nohup bash Uni.qc.mapping.using.star.v3.20250725.sh /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/raw_fastq_batch_1.20250723.lst /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/batch_1 &> ../log/Run.qc.mapping.using.star.v3.20250725.log &
- 使用STAR替换之前的bowtie比对，但是更换了rRNA的index（添加了RNA central的rRNA数据以及将以前的rRNA中的U替换为T，虽然后者不太影响结果）
- STAR参数为：--alignEndsType EndToEnd，--alignIntronMax 1，--outFilterMultimapNmax 10000，--outFilterMismatchNmax 1 

# 二、提取特定长度的reads
nohup bash run.extract.25-34.length.reads.20250725.sh &
- 提取特异性比对的reads中符合特定长度的reads

# 三、测序饱和分析
## 计算reads个数
nohup bash Uni.saturation.calcu.gene.counts.v1.20250728.sh ../processed/batch_1/filtered_bam ../processed/batch_1/geneCounts/ &> ../log/saturation.calcu.gene.counts.20250728.log &
- gene counts以及transcript counts
## 饱和分析
bash Uni.transcript_coverage_analysis_v1_20250725.sh ../processed/batch_1/geneCounts/transcript.counts.txt  ../results/batch_1/coverage_analysis/transcript
bash Uni.transcript_coverage_analysis_v1_20250725.sh ../processed/batch_1/geneCounts/gene.counts.txt  ../results/batch_1/coverage_analysis/gene
- 进行饱和度分析

# 四、评估ribo-seq文库结果
bash ribo.qual.assess.20250730.sh
Rscript Uni.organize.v1.20250730.R /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/results/batch_1/qual_assess \
    /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/figures/batch_1/qual_assess


